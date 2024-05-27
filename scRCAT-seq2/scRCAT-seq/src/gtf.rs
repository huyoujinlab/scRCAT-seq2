use flate2::read::GzDecoder;
use std::{collections::{HashSet}, fs::File, io::BufReader};
use bio::io::gff;

// fn main() -> Result<(), Box<dyn std::error::Error>> {
//     let gtf_path = "gencode.v45.annotation.gtf.gz";
//     let gtf_file = std::fs::File::open(gtf_path)?;
//     let gtf_reader = BufReader::new(GzDecoder::new(gtf_file));
//     let mut gtf_reader = gff::Reader::new(gtf_reader, gff::GffType::GTF2);

//     let mut tree = GeneTree::new();
//     let mut processed_chr = HashSet::new();
    
//     for result in gtf_reader.records() {
//         let record = result?;
//         if record.feature_type() != "gene" {
//             continue;
//         }
//         let gene = &record.attributes().get("gene_id").unwrap();
//         let chr = record.seqname().to_owned();
//         let start = record.start().clone();
//         let end = record.end().clone();

//         if !processed_chr.contains(&chr) {
//             println!("Processing new chromosome: {}", chr);
//             processed_chr.insert(chr.clone());
//         }

//         tree.entry(chr.to_string())
//             .and_modify(|interval_tree| {
//                 interval_tree.insert(start, end, gene.to_string());
//             }).or_insert_with(|| {
//                 let mut interval_tree = IntervalTree::new();
//                 interval_tree.insert(start, end, gene.to_string());
//                 interval_tree
//             });
//     }
//     println!("{:?}", tree.len());
//     Ok(())
// }

#[derive(Debug, Clone)]
struct IntervalTree {
    root: Option<Box<IntervalNode>>,
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
    fn new(start: u64, end: u64, gene_name: String) -> Self {
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
        IntervalTree { root: None }
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
            None => Some(Box::new(IntervalNode::new(start, end, gene_name))),
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
}

type GeneTree = HashMap<String, IntervalTree>;


fn parse_gtf(gtf_path: &str) -> Result<GeneTree, Box<dyn std::error::Error>> {
    let gtf_file = File::open(gtf_path)?;
    // if the file is gzipped, we need to use GzDecoder
    let mut gtf_readers = BufReader::new(&gtf_file);
    if gtf_path.ends_with(".gz") {
        let mut gtf_readers = BufReader::new(GzDecoder::new(&gtf_file));
    }
    let mut gtf_reader = gff::Reader::new(gtf_readers, gff::GffType::GTF2);

    let mut tree = GeneTree::new();
    let mut processed_chr = HashSet::new();
    
    for result in gtf_reader.records() {
        let record = result?;
        if record.feature_type() != "gene" {
            continue;
        }
        let gene = &record.attributes().get("gene_id").unwrap();
        let chr = record.seqname().to_owned();
        let start = record.start().clone();
        let end = record.end().clone();
        
        if !processed_chr.contains(&chr) {
            println!("Processing new chromosome: {}", chr);
            processed_chr.insert(chr.clone());
        }

        tree.entry(chr.to_string())
            .and_modify(|interval_tree| {
                interval_tree.insert(start, end, gene.to_string());
            }).or_insert_with(|| {
                let mut interval_tree = IntervalTree::new();
                interval_tree.insert(start, end, gene.to_string());
                interval_tree
            });
    }
    Ok(tree)
}