use bio::io::gff;
use std::fs::File;
use std::io::BufReader;
use bio_types::strand::Strand;
use flate2::read::GzDecoder;
use std::env;
use std::process;
use std::io::Error;
// use std::collections::HashMap;

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
    exon_intervals: IntervalTree,
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
            exon_intervals: IntervalTree::new(),
        }
    }

    fn add_exon(&mut self, exon: Exon) {
        self.exons.insert(exon.id.clone(), exon.clone());
        self.exon_intervals.insert(exon.start, exon.end, exon.id.clone());
    }

    fn to_ranges(&mut self) {
        self.exon_intervals.to_ranges();
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

fn build_gene_structure(genes: &mut HashMap<String, Gene>, transcripts: &mut HashMap<String, Vec<Transcript>>, exons: &HashMap<String, Vec<Exon>>) {
    for (gene_id, gene) in genes.iter_mut() {
        if let Some(transcript_vec) = transcripts.get_mut(gene_id) {
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

#[derive(Debug, Clone)]
struct IntervalTree {
    root: Option<Box<IntervalNode>>,
    range: Vec<(u64, u64)>
}

#[derive(Debug, Clone)]
struct IntervalNode {
    start: u64,
    end: u64,
    max: u64,
    height: i32,
    gene_name: String,
    left: Option<Box<IntervalNode>>,
    right: Option<Box<IntervalNode>>,
}

impl IntervalNode {
    fn new() -> Self {
        IntervalNode {
            start: 0,
            end: 0,
            max: 0,
            height: 0,
            gene_name: "".to_string(),
            left: None,
            right: None,
        }
    }

    fn attrs(start: u64, end: u64, gene_name: String) -> Self {
        IntervalNode {
            start,
            end,
            max: end,
            height: 1,
            gene_name,
            left: None,
            right: None,
        }
    }

    fn update_max(&mut self) {
        self.max = self.end;
        if let Some(left) = &self.left {
            self.max = self.max.max(left.max);
        }
        if let Some(right) = &self.right {
            self.max = self.max.max(right.max);
        }
    }

    fn update_height(&mut self) {
        let left_height = self.left.as_ref().map_or(0, |node| node.height);
        let right_height = self.right.as_ref().map_or(0, |node| node.height);
        self.height = 1 + left_height.max(right_height);
    }

    fn balance_factor(&self) -> i32 {
        let left_height = self.left.as_ref().map_or(0, |node| node.height);
        let right_height = self.right.as_ref().map_or(0, |node| node.height);
        left_height - right_height
    }
}

impl IntervalTree {
    fn new() -> Self {
        IntervalTree { root: None, range: vec![] }
    }

    fn insert(&mut self, start: u64, end: u64, gene_name: String) {
        let root = self.root.take();
        self.root = self.insert_node(root, start, end, gene_name);
    }
    fn insert_node(
        &mut self,
        node: Option<Box<IntervalNode>>,
        start: u64,
        end: u64,
        gene_name: String,
    ) -> Option<Box<IntervalNode>> {
        match node {
            None => Some(Box::new(IntervalNode::attrs(start, end, gene_name))),
            Some(mut node) => {
                if start < node.start {
                    node.left = self.insert_node(node.left.take(), start, end, gene_name);
                } else {
                    node.right = self.insert_node(node.right.take(), start, end, gene_name);
                }

                node.update_max();
                node.update_height();

                self.balance(node)
            }
        }
    }

    fn balance(&mut self, mut node: Box<IntervalNode>) -> Option<Box<IntervalNode>> {
        let balance_factor = node.balance_factor();

        if balance_factor > 1 {
            if node.left.as_ref().unwrap().balance_factor() < 0 {
                node.left = self.rotate_left(node.left.take().unwrap());
            }
            return self.rotate_right(node);
        }

        if balance_factor < -1 {
            if node.right.as_ref().unwrap().balance_factor() > 0 {
                node.right = self.rotate_right(node.right.take().unwrap());
            }
            return self.rotate_left(node);
        }

        Some(node)
    }

    fn rotate_left(&mut self, mut node: Box<IntervalNode>) -> Option<Box<IntervalNode>> {
        let mut right = node.right.take().unwrap();
        node.right = right.left.take();
        right.left = Some(node);

        right.left.as_mut().unwrap().update_height();
        right.update_height();
        right.left.as_mut().unwrap().update_max();
        right.update_max();

        Some(right)
    }

    fn rotate_right(&mut self, mut node: Box<IntervalNode>) -> Option<Box<IntervalNode>> {
        let mut left = node.left.take().unwrap();
        node.left = left.right.take();
        left.right = Some(node);

        left.right.as_mut().unwrap().update_height();
        left.update_height();
        left.right.as_mut().unwrap().update_max();
        left.update_max();

        Some(left)
    }

    fn search(&self, start: u64, end: u64) -> Vec<String> {
        let mut result = Vec::new();
        self.search_node(self.root.as_ref(), start, end, &mut result);
        result
    }

    fn search_node(
        &self,
        node: Option<&Box<IntervalNode>>,
        start: u64,
        end: u64,
        result: &mut Vec<String>,
    ) {
        if let Some(node) = node {
            if start <= node.max {
                if node.start <= end && start <= node.end {
                    result.push(node.gene_name.clone());
                }

                self.search_node(node.left.as_ref(), start, end, result);
                self.search_node(node.right.as_ref(), start, end, result);
            }
        }
    }

    fn to_ranges(&mut self) {
        let mut ranges = Vec::<(u64, u64)>::new();
        self.to_ranges_node(self.root.as_ref(), &mut ranges);
        ranges.sort_by_key(|&(start, _)| start);
        self.range = ranges;
    }
    
    fn to_ranges_node(
        &self,
        node: Option<&Box<IntervalNode>>,
        ranges: &mut Vec<(u64, u64)>,
    ) {
        if let Some(node) = node {
            ranges.push((node.start, node.end));
            self.to_ranges_node(node.left.as_ref(), ranges);
            self.to_ranges_node(node.right.as_ref(), ranges);
        }
    }

    fn is_sub(&self, readcoord: &ReadCoord, start: usize, strand: char) -> bool {
        let mut a_ptr = start;
        let mut b_ptr = 0;
        let ranges: Vec<(u64, u64)> = self.range.clone();
        let interval = readcoord.to_value();
        // println!("{:?}", ranges);
         // 如果 intervals_b 的第一个区间的开始位置大于 intervals_a 的第一个区间的结束位置, intervals_b 跳过这个区间
        if interval[b_ptr].0 > ranges[a_ptr].1 {
            if a_ptr == ranges.len() - 1 {
                return false;
            }
            return self.is_sub(readcoord, a_ptr + 1, strand);
        }
        while a_ptr < ranges.len() && b_ptr < interval.len() {
            if interval.len() == 1 {
                // 如果read的info中含有TES或者TES，那么reads的起始位点不能超过转录本的起始位点10bp或者终止位点不能超过转录本的终止位点10bp
                if readcoord.info == "TSS" && a_ptr == 0 && strand == '+' {
                    if ranges[a_ptr].0 - 10 <= interval[b_ptr].0 && ranges[a_ptr].1 + 1 >= interval[b_ptr].1 {
                        return true;
                    } else {
                        // println!("interval is one, {:?}, {:?}", interval, ranges[a_ptr]);
                        return false;
                    }
                } else if readcoord.info == "TSS" && a_ptr == ranges.len() - 1 && strand == '-' {
                    if ranges[a_ptr].0 - 1 <= interval[b_ptr].0 && ranges[a_ptr].1 + 10 >= interval[b_ptr].1 {
                        return true;
                    } else {
                        // println!("interval is one, {:?}, {:?}", interval, ranges[a_ptr]);
                        return false;
                    }
                } else if readcoord.info == "TES" && a_ptr == ranges.len() - 1 && strand == '+' {
                    if ranges[a_ptr].0 - 1 <= interval[b_ptr].0 && ranges[a_ptr].1 + 10 >= interval[b_ptr].1 {
                        return true;
                    } else {
                        // println!("interval is one, {:?}, {:?}", interval, ranges[a_ptr]);
                        return false;
                    }
                } else if readcoord.info == "TES" && a_ptr == 0 && strand == '-' {
                    if ranges[a_ptr].0 - 10 <= interval[b_ptr].0 && ranges[a_ptr].1 + 1 >= interval[b_ptr].1 {
                        return true;
                    } else {
                        // println!("interval is one, {:?}, {:?}", interval, ranges[a_ptr]);
                        return false;
                    }
                } else {
                    // 如果是单区间，允许出现1个碱基的错配，因为在比对的时候，可能会出现错配
                    if ranges[a_ptr].0 - 1 <= interval[b_ptr].0 && ranges[a_ptr].1 + 1 >= interval[b_ptr].1 {
                        return true;
                    } else {
                        // println!("interval is one, {:?}, {:?}", interval, ranges[a_ptr]);
                        return false;
                    }
                }
            } else if interval.len() < 1 {
                println!("error, interval length less than 1");
            } else {
                // 和单区间的判断方法一致，但是这里需要注意的是转录起始和终止位点就只出现在最后一个区间和第一个区间，所以分开判断
                if b_ptr == 0 {
                    // 在第一个区间中，只可能出现正链的转录起始位点和负链的转录终止位点，因为所有的reads都是按照染色体的起始位点进行排序的
                    if readcoord.info == "TSS" && strand == '+' && a_ptr == 0 {
                        if !(ranges[a_ptr].0 - 10 <= interval[b_ptr].0 && ranges[a_ptr].1 == interval[b_ptr].1) {
                            return false;
                        }
                    } else if readcoord.info == "TES" && strand == '-' && a_ptr == 0 {
                        if !(ranges[a_ptr].0 - 10 <= interval[b_ptr].0 && ranges[a_ptr].1 == interval[b_ptr].1) {
                            return false;
                        }
                    } else {
                        if !(ranges[a_ptr].0 <= interval[b_ptr].0 && ranges[a_ptr].1 == interval[b_ptr].1) {
                            // println!("type name of range: {:?}, interval: {:?}", std::any::type_name(ranges[a_ptr].0), std::any::type_name(interval[b_ptr].0));
                            // println!("error, interval not match at 0, {:?}, {:?}", interval[b_ptr], ranges[a_ptr]);
                            return false;
                        }
                    }
                } else if b_ptr == interval.len() - 1 {
                    // 在最后一个区间内，只可能出现正链的转录终止位点和负链的转录起始位点
                    if readcoord.info == "TES" && strand == '+' && a_ptr == ranges.len() - 1 {
                        if !(ranges[a_ptr].0 == interval[b_ptr].0 && ranges[a_ptr].1 + 10 >= interval[b_ptr].1) {
                            return false;
                        }
                    } else if readcoord.info == "TSS" && strand == '-' && a_ptr == ranges.len() - 1 {
                        if !(ranges[a_ptr].0 == interval[b_ptr].0 && ranges[a_ptr].1 + 10 <= interval[b_ptr].1) {
                            return false;
                        }
                    }  else {
                        if !(ranges[a_ptr].0 == interval[b_ptr].0 && ranges[a_ptr].1 >= interval[b_ptr].1) {
                            // println!("error, interval not match at last, {:?}, {:?}", interval[b_ptr], ranges[a_ptr]);
                            return false;
                        }
                    }
                } else {
                    if !(ranges[a_ptr].0 == interval[b_ptr].0 && ranges[a_ptr].1 == interval[b_ptr].1) {
                        // println!("error, interval not match at {:?}, {:?}, {:?}", b_ptr, interval[b_ptr], ranges[a_ptr]);
                        return false;
                    }
                }
            }
            a_ptr += 1;
            b_ptr += 1;
        }
        // check if the last interval is the same
        if b_ptr != interval.len() {
            // println!("error, interval have more exon at {:?}, {:?}, {:?}", b_ptr, interval[b_ptr], ranges[a_ptr]);
            return false;
        } else {
            return true;
        }
    }
}

type GeneTree = HashMap<String, IntervalTree>;

// fn main() -> Result<(), Error> {
//     let args: Vec<String> = env::args().collect();
//     if args.len() != 2 {
//         eprintln!("Usage: {} <gff_file>", args[0]);
//         process::exit(1);
//     }
//     let gff_file = &args[1];
//     println!("Reading GFF file: {}", gff_file);

//     // define variables to store the data parsed from the GFF file
//     let (gene_tree, genes) = parse_gtf(gff_file)?;
//     if let Some(gene) = genes.get("ENSG00000243485.5") {
//         println!("{:?}", gene);
//     }
//     Ok(())
// }

fn parse_gtf(path: &str) -> Result<(GeneTree, HashMap<String, Gene>), Error> {
    // let reader = BufReader::new(GzDecoder::new(File::open(path)?));
    let file = File::open(path)?;
    let reader: Box<dyn std::io::Read> = if path.ends_with(".gz") {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let mut gtf_reader = gff::Reader::new(BufReader::new(reader), gff::GffType::GTF2);

    let mut gene_tree = GeneTree::new();
    let mut genes: HashMap<String, Gene> = HashMap::new();
    let mut transcripts: HashMap<String, Vec<Transcript>> = HashMap::new();
    let mut exons: HashMap<String, Vec<Exon>> = HashMap::new();

    let mut exon_no = 0;
    for result in gtf_reader.records() {
        let record = result?;
        match record.feature_type() {
            "gene" => {
                let gene_id = record.attributes().get("gene_id").unwrap().to_string();
                let chr = record.seqname().to_string();
                let start = record.start();
                let end = record.end();
                let strand = record.strand().unwrap();
                let gene = Gene::new(gene_id.clone(), chr.clone(), *start, *end, strand);
                genes.insert(gene_id.clone(), gene);
                gene_tree.entry(chr)
                    .and_modify(|tree| {
                        tree.insert(*start, *end, gene_id.clone());
                    })
                    .or_insert_with(|| {
                        let mut tree = IntervalTree::new();
                        tree.insert(*start, *end, gene_id.clone());
                        tree
                    });
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
                let transcript_id = record.attributes().get("transcript_id").unwrap().to_string();
                let exon_id = record.attributes().get("exon_id").map(|v| v.to_string())
                    .unwrap_or_else(|| {
                        exon_no += 1;
                        format!("{}_{}", transcript_id, exon_no)
                    });
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
    build_gene_structure(&mut genes, &mut transcripts, &mut exons);
    for (_, g) in genes.iter_mut() {
        for (_, t) in g.transcripts.iter_mut() {
            t.exon_intervals.to_ranges();
        }
    }
    Ok((gene_tree, genes))
}


fn strand_to_symbol(strand: &Strand) -> char {
    match strand {
        Strand::Forward => '+',
        Strand::Reverse => '-',
        _ => '.',
    }
}