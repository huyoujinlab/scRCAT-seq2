mod fmi;
use bio::pattern_matching;
use gzp::par::compress::Compression;
use gzp::deflate::Gzip;
use gzp::ZBuilder;
use std::fs::File;
use std::io::{self, BufReader};
use flate2::read::MultiGzDecoder;
use bio::io::fastq::{self, FastqRead};
use clap::Parser;
use bio::alphabets::dna::{self, alphabet};
use threadpool::ThreadPool;
// use std::sync::mpsc::channel;
use crossbeam::channel::{bounded, Receiver};
use std::sync::{Arc, Mutex};


const BATCH_SIZE: usize = 100_000;

/// Searches for sequence patterns in paired-end FASTQ files.
#[derive(Parser, Debug)]
#[clap(author = "Your Name", version = "1.0", about)]
struct Args {
    /// Input FASTQ file for R1, optionally gzipped. R1 contains the cDNA sequence.
    #[clap(short = 'i', long, required = true)]
    input1: String,
    /// Input FASTQ file for R2, optionally gzipped. R2 contains the Cell Barcode and UMI sequence.
    #[clap(short = 'I', long, required = true)]
    input2: String,
    /// Output FASTQ file for R1, containing the cDNA sequence from the R1 file.
    #[clap(short = 'o', long, required = true)]
    output1: String,
    /// Output FASTQ file for R2, containing the cDNA sequence from the R2 file.
    #[clap(short = 'O', long, required = true)]
    output2: String,
    /// DNA sequence to search for in R2, used to locate the UMI sequence.
    #[clap(short = 'p', long, required = true)]
    pattern1: String,
    /// If false, search for the CB and UMI pattern in the forward direction; if true, search in the reverse direction.
    #[clap(short, long, default_value = "false", parse(try_from_str), possible_values = &["true", "false"])]
    reverse: bool,
    /// If true, remove the polyA sequence from the cDNA sequence.
    #[clap(short = 'a', long, default_value = "false", parse(try_from_str), possible_values = &["true", "false"])]
    remove_polya: bool,
     /// If true, only trim the sequence ahead of adapter.
     #[clap(short = 'd', long, default_value = "false", parse(try_from_str), possible_values = &["true", "false"])]
     only_cut_adpater: bool,
    /// DNA sequence to search for in R2, used as a helper sequence to locate the CB and UMI.
    #[clap(short = 'q', long, default_value = "")]
    pattern2: String,
    /// DNA sequence to search for in R1, used for trimming the TSO and other sequences in R1 sequence and find the cDNA sequence in R2.
    #[clap(short = 'P', long, required_if_eq("remove_polya", "false"))]
    pattern3: Option<String>,
    /// Length of UMI & Cell Barcode sequence between pattern1 and pattern2 (if exists).
    #[clap(short, long, required = true)]
    umi: usize,
    /// Maximum number of mismatches allowed when searching for patterns.
    #[clap(short, long, default_value = "1")]
    max_mismatches: u8,
    /// Number of threads to use for parallel processing.
    #[clap(short, long, default_value = "8")]
    threads: u8,
    #[clap(long, default_value = "false", parse(try_from_str), possible_values = &["true", "false"])]
    verbose: bool,
}

fn main() -> io::Result<()> {
    let args = Args::parse();
    println!("{:?}", args);

    let input1_path = args.input1;
    let input2_path = args.input2;
    let output1_path = args.output1;
    let output2_path = args.output2;
    // let output3_path = args.output3;
    let pattern1 = args.pattern1.as_bytes().to_vec();
    let pattern2 = args.pattern2.as_bytes().to_vec();
    let pattern3 = if let Some(p3) = args.pattern3 {
        p3.as_bytes().to_vec()
    } else {
        Vec::new()
    };
    let reverse = args.reverse;
    let umi_len = args.umi;
    let max_mismatches = args.max_mismatches as usize;
    let num_threads = args.threads as usize;
    let remove_polya = args.remove_polya;
    let only_cut_adpater = args.only_cut_adpater;
    let verbose = args.verbose;

    let pattern3_rev = if !pattern3.is_empty() {
        dna::revcomp(&pattern3)
    } else {
        Vec::new()
    };
    let alphabet = dna::iupac_alphabet();

    // print!("pattern1: {:?}, pattern2: {:?}, pattern3: {:?}, pattern3_rev: {:?}, umi_len: {}, max_mismatches: {}, reverse: {}, alphabet: {:?}", pattern1, pattern2, pattern3, pattern3_rev, umi_len, max_mismatches, reverse, alphabet);

    let file1 = File::open(input1_path.clone()).expect("Failed to open R1 file");
    let file2 = File::open(input2_path.clone() ).expect("Failed to open R2 file");
    let reader1: Box<dyn io::Read> = if input1_path.clone().ends_with(".gz") {
        Box::new(MultiGzDecoder::new(file1))
    } else {
        Box::new(file1)
    };
    let reader2: Box<dyn io::Read> = if input2_path.clone().ends_with(".gz") {
        Box::new(MultiGzDecoder::new(file2))
    } else {
        Box::new(file2)
    };
    let mut fastq_reader1 = fastq::Reader::new(BufReader::new(reader1));
    let mut fastq_reader2 = fastq::Reader::new(BufReader::new(reader2));

    let output_threads = num_threads / 2;
    let par_com1 = ZBuilder::<Gzip, _>::new()
        .num_threads(output_threads)
        .compression_level(Compression::default())
        .from_writer(File::create(output1_path)?);
    let mut fastq_writer1 = fastq::Writer::new(par_com1);
    let par_com2 = ZBuilder::<Gzip, _>::new()
        .num_threads(output_threads)
        .compression_level(Compression::default())
        .from_writer(File::create(output2_path)?);
    let mut fastq_writer2 = fastq::Writer::new(par_com2);

    // output file for CB and UMI
    // let par_com3 = ZBuilder::<Gzip, _>::new()
    //     .num_threads(output_threads)
    //     .compression_level(Compression::default())
    //     .from_writer(File::create(output3_path)?);
    // let mut fastq_writer3 = fastq::Writer::new(par_com3);

    let pool = ThreadPool::new(num_threads - 1);
    let (sender_a, receiver_a) = bounded(1024);
    let (sender_b, receiver_b) = bounded(1024);

    // 创建一个线程读取数据，然后传递给通道 a
    for (record1, record2) in fastq_reader1.records().zip(fastq_reader2.records()) {
        // let (Ok(record1), Ok(record2)) = (record1, record2) {
            let record1 = record1.unwrap();
            let record2 = record2.unwrap();
            sender_a.send((record1, record2));
    }

    // 创建多个线程处理通道 a 中的数据，然后传递给通道 b
    for _ in 0..num_threads - 1 {
        let receiver_a = receiver_a.clone();
        let sender_b = sender_b.clone();
        let pattern1_cloned = pattern1.clone();
        let pattern2_cloned = pattern2.clone();
        let pattern3_cloned = pattern3.clone();
        let alphabet_cloned = alphabet.clone();

        pool.execute(move || {
            for (record1, record2) in receiver_a.iter() {
                let (new_record1, new_record2) = fmi::process_record(
                    &record1,
                    &record2,
                    &pattern1_cloned,
                    &pattern2_cloned,
                    &pattern3_cloned,
                    &alphabet_cloned,
                    umi_len,
                    max_mismatches,
                    remove_polya,
                    reverse,
                    only_cut_adpater,
                    verbose,
                );
                sender_b.send((new_record1, new_record2)).unwrap();
            }
        });
    }

    drop(sender_b);

    // 在主线程中读取通道 b 中的数据，写入文件
    for received in receiver_b.iter() {
        if let (Some(new_record1), Some(new_record2)) = (&received.0, &received.1) {
            fastq_writer1.write_record(new_record1)?;
            fastq_writer2.write_record(new_record2)?;
        }
    }
    pool.join();
    Ok(())
}