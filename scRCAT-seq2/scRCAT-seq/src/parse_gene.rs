use bio::io::gff;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use bio_types::strand::Strand;
use flate2::read::GzDecoder;

#[derive(Debug, Clone)]
struct Gene {
    id: String,
    chr: String,
    start: u64,
    end: u64,
    strand: Strand,
    transcripts: HashMap<String, Transcript>,
}

#[derive(Debug, Clone)]
struct Transcript {
    id: String,
    gene_id: String,
    chr: String,
    start: u64,
    end: u64,
    strand: Strand,
    exons: HashMap<String, Exon>,
    exon_intervals: IntervalNode,
}

#[derive(Debug, Clone)]
struct Exon {
    id: String,
    transcript_id: String,
    chr: String,
    start: u64,
    end: u64,
    strand: Strand,
}

#[derive(Debug, Clone)]
struct IntervalNode {
    root: Option<Box<IntervalTreeNode>>,
}

impl IntervalNode {
    fn new() -> Self {
        IntervalNode { root: None }
    }

    fn insert(&mut self, id: String, start: u64, end: u64) {
        if let Some(ref mut root) = self.root {
            root.insert(id, start, end);
        } else {
            self.root = Some(Box::new(IntervalTreeNode::new(id, start, end)));
        }
    }
}

#[derive(Debug, Clone)]
struct IntervalTreeNode {
    id: String,
    start: u64,
    end: u64,
    left: Option<Box<IntervalTreeNode>>,
    right: Option<Box<IntervalTreeNode>>,
}

impl IntervalTreeNode {
    fn new(id: String, start: u64, end: u64) -> Self {
        IntervalTreeNode {
            id,
            start,
            end,
            left: None,
            right: None,
        }
    }

    fn insert(&mut self, id: String, start: u64, end: u64) {
        if end < self.start {
            if let Some(ref mut left) = self.left {
                left.insert(id, start, end);
            } else {
                self.left = Some(Box::new(IntervalTreeNode::new(id, start, end)));
            }
        } else if start > self.end {
            if let Some(ref mut right) = self.right {
                right.insert(id, start, end);
            } else {
                self.right = Some(Box::new(IntervalTreeNode::new(id, start, end)));
            }
        } else {
            if start < self.start {
                self.start = start;
            }
            if end > self.end {
                self.end = end;
            }
            if let Some(ref mut left) = self.left {
                if left.end < end {
                    left.insert(id.clone(), start, end);
                }
            }
            if let Some(ref mut right) = self.right {
                if right.start > start {
                    right.insert(id, start, end);
                }
            }
        }
    }
}

impl Gene {
    fn new(id: String, chr: String, start: u64, end: u64, strand: Strand) -> Self {
        Gene {
            id,
            chr,
            start,
            end,
            strand,
            transcripts: HashMap::new(),
        }
    }

    fn add_transcript(&mut self, transcript: Transcript) {
        self.transcripts.insert(transcript.id.clone(), transcript);
    }
}

impl Transcript {
    fn new(id: String, gene_id: String, chr: String, start: u64, end: u64, strand: Strand) -> Self {
        Transcript {
            id,
            gene_id,
            chr,
            start,
            end,
            strand,
            exons: HashMap::new(),
            exon_intervals: IntervalNode::new(),
        }
    }

    fn add_exon(&mut self, exon: Exon) {
        self.exons.insert(exon.id.clone(), exon.clone());
        self.exon_intervals.insert(exon.id.clone(), exon.start, exon.end);
    }
}

impl Exon {
    fn new(id: String, transcript_id: String, chr: String, start: u64, end: u64, strand: Strand) -> Self {
        Exon {
            id,
            transcript_id,
            chr,
            start,
            end,
            strand,
        }
    }

    fn is_overlap(&self, other: &Exon) -> bool {
        self.start <= other.end && self.end >= other.start
    }
}

// fn add_exon_to_transcript(transcripts: &mut HashMap<String, Vec<Transcript>>, exons: &HashMap<String, Vec<Exon>>) {
//     for (gene_id, ts_vec) in transcripts.iter_mut() {
//         for ts in ts_vec.iter_mut() {
//             if let Some(exon_vec) = exons.get(&ts.id) {
//                 for ex in exon_vec {
//                     ts.add_exon(ex.clone());
//                 }
//             }
//         }
//     }
// }

// fn add_transcript_to_gene(genes: &mut HashMap<String, Gene>, transcripts: &HashMap<String, Vec<Transcript>>) {
//     for (gene_id, ts_vec) in transcripts.iter() {
//         if let Some(gene) = genes.get_mut(gene_id) {
//             for ts in ts_vec {
//                 gene.add_transcript(ts.clone());
//             }
//         }
//     }
// }

fn build_gene_structure(
    genes: &mut HashMap<String, Gene>,
    transcripts: &mut HashMap<String, Vec<Transcript>>,
    exons: &HashMap<String, Vec<Exon>>,
) {
    for (gene_id, gene) in genes.iter_mut() {
        if let Some(transcript_vec) = transcripts.get(gene_id) {
            for transcript in transcript_vec {
                if let Some(exon_vec) = exons.get(&transcript.id) {
                    for ex in exon_vec {
                        transcript.add_exon(ex.clone());
                    }
                }
                gene.add_transcript(transcript.clone());
            }
        }
    }
}

fn main() {
    let file = File::open("gencode.v45.annotation.gtf.gz").unwrap();
    let bufread = BufReader::new(GzDecoder::new(file));
    let mut gtf = gff::Reader::new(bufread, gff::GffType::GTF2);
    let mut genes: HashMap<String, Gene> = HashMap::new();
    let mut transcripts: HashMap<String, Vec<Transcript>> = HashMap::new();
    let mut exons: HashMap<String, Vec<Exon>> = HashMap::new();
    for record in gtf.records() {
        let record = record.unwrap();
        match record.feature_type() {
            "gene" => {
                let gene_id = record.attributes().get("gene_id").unwrap().to_string();
                let chr = record.seqname().to_string();
                let start = record.start();
                let end = record.end();
                let strand = record.strand().unwrap();
                let gene = Gene::new(gene_id.clone(), chr, *start, *end, strand);
                genes.insert(gene_id.clone(), gene);
            }
            "transcript" => {
                let transcript_id = record.attributes().get("transcript_id").unwrap().to_string();
                let gene_id = record.attributes().get("gene_id").unwrap().to_string();
                let chr = record.seqname().to_string();
                let start = record.start();
                let end = record.end();
                let strand = record.strand().unwrap();
                let transcript = Transcript::new(transcript_id.clone(), gene_id.clone(), chr, *start, *end, strand);
                transcripts.entry(gene_id.clone()).or_insert(Vec::new()).push(transcript);
            }
            "exon" => {
                let exon_id = record.attributes().get("exon_id").unwrap().to_string();
                let transcript_id = record.attributes().get("transcript_id").unwrap().to_string();
                let chr = record.seqname().to_string();
                let start = record.start();
                let end = record.end();
                let strand = record.strand().unwrap();
                let exon = Exon::new(exon_id.clone(), transcript_id.clone(), chr, *start, *end, strand);
                exons.entry(transcript_id.clone()).or_insert(Vec::new()).push(exon);
            }
            _ => {}
        }
    }
    // add_exon_to_transcript(&mut transcripts, &exons);
    // add_transcript_to_gene(&mut genes, &transcripts);
    build_gene_structure(&mut genes, &mut transcripts, &exons);
    println!("{:?}", genes.len());
    println!("{:?}", transcripts.len());
    println!("{:?}", exons.len());
    if let Some(ts) = transcripts.get("ENSG00000243485.5") {
        println!("{:?}", ts);
    }
    if let Some(gene) = genes.get("ENSG00000243485.5") {
        println!("{:?}", gene);
    }
}