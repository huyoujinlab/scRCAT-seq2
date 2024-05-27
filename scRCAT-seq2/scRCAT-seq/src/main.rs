use rust_htslib::{bam::{self, Read}, bcf::header};
use std::collections::HashMap;
use clap::Parser;
use crossbeam::channel::{bounded, Receiver};
use std::sync::{Arc, Mutex};
use threadpool::ThreadPool;

// file: struct.rs
include!("struct.rs");
include!("gtf_parse.rs");


#[derive(Parser, Debug)]
#[clap(author = "Your Name", version = "1.0", about)]
struct Args {
    /// Input bam files.
    #[clap(short = 'b', long, required = true)]
    bam: String,
    /// Input gtf files.
    #[clap(short = 'g', long, required = true)]
    gtf: String,
    /// Output file path, output is csv format.
    #[clap(short = 'o', long, required = true)]
    out: String,
    /// reverse complementary  CB and UMI
    #[clap(short = 'r', long, default_value = "false", parse(try_from_str), possible_values = &["true", "false"])]
    revcmp: bool,
    // /// reverse complementary  CB and UMI
    // #[clap(short = 'm', long, required = false)]
    // mrevcmp: String,
}

// fn main() {
//     let args = Args::parse();
//     let revcmp = args.revcmp;

//     println!("{:?}", args);

//     let gtf_path = args.gtf.as_str();
//     let (mut gene_tree, mut gene_map) = parse_gtf(gtf_path).unwrap();
//     println!("Gene tree: {:?}, Gene map: {:?}", gene_tree.len(), gene_map.len());
//     let mut bam = bam::Reader::from_path(args.bam).unwrap();
//     // let mut bam = bam::Reader::from_path("ACACCCAACTAC.bam").unwrap();
//     let header: bam::HeaderView = bam.header().clone();
//     let tid_to_name: HashMap<_, _> = header.target_names().iter().enumerate().map(|(tid, &name)| (tid as u32, String::from_utf8_lossy(name).into_owned())).collect();

//     let mut recMap: RecordMap = HashMap::new();

//     for r in bam.records() {
//         let record = r.unwrap();
//         rec2info(record, &mut recMap, &tid_to_name, &gene_tree, revcmp);
//     }
//     println!("Parse BAM done!");
//     let res_recMap = add_transcript(&mut recMap, &mut gene_map);
//     println!("Add transcript done!");

//     let mut df = recMap2df(&res_recMap).unwrap();
//     println!("{:?}", df);
//     let mut csv = File::create(args.out).unwrap();
//     CsvWriter::new(&mut csv).finish(&mut df).unwrap();
// }

fn main() {
    let args = Args::parse();
    let revcmp = args.revcmp;

    println!("{:?}", args);

    let gtf_path = args.gtf.as_str();
    let (mut gene_tree, mut gene_map) = parse_gtf(gtf_path).unwrap();
    println!("Gene tree: {:?}, Gene map: {:?}", gene_tree.len(), gene_map.len());


    let mut recMap: RecordMap = HashMap::new();
    let rec_map: Arc<Mutex<HashMap<String, RecordInfo>>> = Arc::new(Mutex::new(HashMap::new()));

    if args.bam.contains(",") {
        let all_bam = args.bam.split(",").collect::<Vec<&str>>();
        // let all_rev = args.mrevcmp.split(",").collect::<Vec<&str>>();
        let all_rev = vec!["false"; all_bam.len()];
        for (bam_path, rev) in all_bam.iter().zip(all_rev.iter()) {
            println!("Parsing {:?}!", bam_path);
            let revss = if *rev == "true" {
                true
            } else {
                false
            };
            let mut bam = bam::Reader::from_path(bam_path).unwrap();
            let header = bam.header().clone();
            let tid_to_name: HashMap<_, _> = header.target_names().iter().enumerate().map(|(tid, &name)| (tid as u32, String::from_utf8_lossy(name).into_owned())).collect();
            for r in bam.records() {
                let record = r.unwrap();
                rec2info(record, &mut recMap, &tid_to_name, &gene_tree, revss);
            }
            println!("Parsing {:?} finished!", bam_path);
        }
    // } else {
    //     let mut bam = bam::Reader::from_path(args.bam).unwrap();
    //     // let mut bam = bam::Reader::from_path("ACACCCAACTAC.bam").unwrap();
    //     let header: bam::HeaderView = bam.header().clone();
    //     let tid_to_name: HashMap<_, _> = header.target_names().iter().enumerate().map(|(tid, &name)| (tid as u32, String::from_utf8_lossy(name).into_owned())).collect();


    //     for r in bam.records() {
    //         let record = r.unwrap();
    //         rec2info(record, &mut recMap, &tid_to_name, &gene_tree, revcmp);
    //     }
    //     println!("Parse BAM done!");
    } else {
        let mut bam = bam::Reader::from_path(args.bam).unwrap();
        // let mut bam = bam::Reader::from_path("ACACCCAACTAC.bam").unwrap();
        let header: bam::HeaderView = bam.header().clone();

        let tid_to_name: Arc<HashMap<_, _>> = Arc::new(header.target_names().iter().enumerate().map(|(tid, &name)| (tid as u32, String::from_utf8_lossy(name).into_owned())).collect());
        let gene_tree: Arc<GeneTree> = Arc::new(gene_tree);

        let num_threads = 8;
        let pool = ThreadPool::new(num_threads);
        let (sender, receiver): (_, Receiver<_>) = bounded(1024);

        pool.execute(move || {
            for record in bam.records() {
                let record = record.unwrap();
                sender.send(record).unwrap();
            }
            drop(sender);
        });

        for _ in 0..num_threads {
            let receiver = receiver.clone();
            let tid_to_name = tid_to_name.clone();
            let gene_tree = gene_tree.clone();
            let mut rec_map = rec_map.clone();

            pool.execute(move || {
                for record in receiver.iter() {
                    if let Some(res) = record2info(&record, &tid_to_name, &gene_tree) {
                        // rec_maentry(res.key).or_insert(res.value);
                        add_rec2map(res, &mut rec_map);
                    }
                }
            });
        }
        pool.join();
        println!("Parse BAM done! len: {}", rec_map.lock().unwrap().len());
        
    }
    let mut map = match Arc::try_unwrap(rec_map) {
        Ok(mutex) => mutex.into_inner().unwrap(),
        Err(rec_map) => {
            // 如果还有其他的 `Arc` 引用,`try_unwrap` 会失败并返回原始的 `Arc`
            // 在这种情况下,你可能需要重新考虑你的设计
            panic!("Failed to unwrap Arc!");
        }
    };
    let mut res_recMap = add_transcript(map, &mut gene_map);
    println!("Add transcript done!");
    // let mut res_recMaps = match Arc::try_unwrap(res_recMap) {
    //     Ok(mutex) => mutex.into_inner().unwrap(),
    //     Err(res_recMap) => {
    //         // 如果还有其他的 `Arc` 引用,`try_unwrap` 会失败并返回原始的 `Arc`
    //         // 在这种情况下,你可能需要重新考虑你的设计
    //         panic!("Failed to unwrap Arc!");
    //     }
    // };
    let mut df = recMap2df(&res_recMap).unwrap();
    println!("{:?}", df);
    let mut csv = File::create(args.out).unwrap();
    CsvWriter::new(&mut csv).finish(&mut df).unwrap();
}