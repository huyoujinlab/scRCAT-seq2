use polars::prelude::*;
use bio::alphabets::dna;
use rayon::prelude::*;

#[derive(Debug, Clone)]
struct RecordInfo {
    umi: String,
    cell_barcode: String,
    gene: String,
    chromosome: String,
    coordinates: Vec<ReadCoord>,
    potential_gene: HashMap<String, i32>,
    all_trans: HashMap<String, usize>,
    uniq_trans: String,
    select_trans: String,
}
    
#[derive(Debug, Clone)]
struct ReadCoord {
    coord: Vec<Vec<usize>>,
    info: String,
}

impl ReadCoord {
    fn new() -> ReadCoord {
        ReadCoord {
            coord: Vec::new(),
            info: String::new(),
        }
    }

    fn to_string(&self) -> String {
        let mut result = Vec::new();
        let mut tmp = Vec::new();
        for coord in &self.coord {
            tmp = Vec::new();
            for i in 0..coord.len() {
                tmp.push(coord[i].to_string());
            }
            result.push(tmp.join("-"));
        }
        result.join(",")
    }

    fn to_value(&self) -> Vec<(u64, u64)> {
        let mut result = Vec::new();
        for coord in &self.coord {
            if coord.len() != 2 {
                println!("Error: coord length is not 2");
                process::exit(1);
            }
            result.push((coord[0] as u64, coord[1] as u64));
        }
        result
    }
}

type RecordMap = HashMap<String, RecordInfo>;


impl RecordInfo {
    fn new() -> RecordInfo {
        RecordInfo {
            umi: String::new(),
            cell_barcode: String::new(),
            gene: String::new(),
            chromosome: String::new(),
            coordinates: Vec::new(),
            potential_gene: HashMap::<String, i32>::new(),
            all_trans: HashMap::<String, usize>::new(),
            uniq_trans: String::new(),
            select_trans: String::new(),
        }
    }

    fn with_attr(chrom: String, gene: String, umi: String, cell_barcode: String) -> RecordInfo {
        RecordInfo {
            umi: umi,
            cell_barcode: cell_barcode,
            gene: gene,
            chromosome: chrom,
            coordinates: Vec::new(),
            potential_gene: HashMap::<String, i32>::new(),
            all_trans: HashMap::<String, usize>::new(),
            uniq_trans: String::new(),
            select_trans: String::new(),
        }
    }

    fn to_string(&self) -> String {
        let mut result = String::new();
        result.push_str(&self.umi);
        result.push_str("\t");
        result.push_str(&self.cell_barcode);
        result.push_str("\t");
        result.push_str(&self.gene);
        result.push_str("\t");
        result.push_str(&self.chromosome);
        result.push_str("\t");
        for coord in &self.coordinates {
            result.push_str(&coord.info);
            result.push_str("\t");
        }
        result
    }

    fn add_coord(&mut self, coord: ReadCoord) {
        self.coordinates.push(coord);
    }

    fn add_potential_gene(&mut self, gene: Vec<String>) {
        for g in gene {
            self.potential_gene.entry(g).and_modify(|count| {
                *count += 1;
            }).or_insert(1);
        }
    }

    fn pg_to_string(&self) -> String {
        let mut result = Vec::new();
        for (gene, count) in &self.potential_gene {
            result.push(format!("{}:{}", gene, count));
        }
        result.join(";")
    }

    fn set_umi(&mut self, umi: String) {
        self.umi = umi;
    }

    fn set_cb(&mut self, cell_barcode: String) {
        self.cell_barcode = cell_barcode;
    }

    fn set_gene(&mut self, gene: String) {
        if self.gene == "-" {
            self.gene = gene;
        } else if self.gene.is_empty() {
            self.gene = gene;
        } else {
            return;
        }
    }

    fn set_chromosome(&mut self, chromosome: String) {
        self.chromosome = chromosome;
    }

    fn coord_to_string(&self) -> String {
        let mut result = String::new();
        for coord in &self.coordinates {
            result.push_str(&coord.to_string());
            result.push_str(";");
        }
        result
    }

    fn all_trans_to_string(&self) -> String {
        let mut result = vec![];
        if self.all_trans.len() > 0 {
            for (key, value) in self.all_trans.iter() {
                result.push(format!("{}:{}", key, value));
            }
            result.join(";")
        } else {
            result.push("NA".to_string());
            result.join(";")
        }
    }

    fn select_gene(&self) -> String {
        if self.gene != "-" {
            if self.potential_gene.contains_key(&self.gene) || !self.potential_gene.contains_key(&self.gene) {
                return self.gene.clone();
            }
        }
        if self.gene == "-" {
            let mut max_value = 0;
            let mut max_key = String::new();
            
            for (key, value) in &self.potential_gene {
                if *value > max_value {
                    max_value = *value;
                    max_key = key.clone();
                }
            }
            
            return max_key;
        }
        String::new()
    }
}

fn get_aux_value(record: &bam::Record, tag: &[u8]) -> String {
    match record.aux(tag) {
        Ok(aux) => {
            match aux {
                bam::record::Aux::String(s) => String::from(s),
                // bam::record::Aux::Int8(i) => Some(i.to_string()),
                // bam::record::Aux::Int16(i) => Some(i.to_string()),
                // bam::record::Aux::Int32(i) => Some(i.to_string()),
                // bam::record::Aux::Uint8(u) => Some(u.to_string()),
                // bam::record::Aux::Uint16(u) => Some(u.to_string()),
                // bam::record::Aux::Uint32(u) => Some(u.to_string()),
                bam::record::Aux::Float(f) => f.to_string(),
                _ => "".to_string(),
            }
        }
        Err(_) => "".to_string(),
    }
}

fn rec2info(record: bam::Record, rec_map: &mut RecordMap, header: &HashMap<u32, String>, gene_tree: &HashMap<String, IntervalTree>, revcmp: bool) {
    if record.mapq() != 255 {
        return;
    }
    let gene: String = get_aux_value(&record, b"GX");
    let mut umi = get_aux_value(&record, b"UB");
    let mut cell_barcode = get_aux_value(&record, b"CB");
    if revcmp {
        umi = String::from_utf8_lossy(&dna::revcomp(umi.into_bytes())).to_string();
        cell_barcode = String::from_utf8_lossy(&dna::revcomp(&cell_barcode.into_bytes())).to_string();
    };
    // let chromosome = String::from_utf8_lossy(header.tid2name(record.tid() as u32)).to_string();
    let chromosome = header.get(&(record.tid() as u32)).unwrap().to_string();
    let coord = get_coord(&record);
    let start = record.pos() as u64;
    let end = record.cigar().end_pos() as u64;
    // let poten_gene = gene_tree.get(&chromosome).unwrap().search(start, end);
    let poten_gene = if let Some(gene) = gene_tree.get(&chromosome) {
        gene.search(start, end)
    } else {
        vec!["".to_string()]
    };
    // println!("Potentail Gene: {:?}", poten_gene);

    // println!("{:?}, {:?}", coord, coord.to_string());

    let key = format!("{}-{}-{}", umi, cell_barcode, chromosome);
    // println!("Key: {}, gene: {:?}, poten_gene: {:?}", key, gene, poten_gene);
    rec_map.entry(key.clone()).and_modify(|full_info| {
        full_info.add_coord(coord.clone());
        full_info.add_potential_gene(poten_gene);
        full_info.set_gene(gene.clone());
    }).or_insert_with(|| {
        let mut rec_info = RecordInfo::new();
        rec_info.set_umi(umi);
        rec_info.set_cb(cell_barcode);
        rec_info.set_gene(gene);
        rec_info.set_chromosome(chromosome);
        rec_info.add_coord(coord);
        rec_info
    });
    // println!("{:?}, {:?}", rec_map.len(), rec_map)
}

fn record2info(record: &bam::Record, chrom_map: &HashMap<u32, String>,  gene_tree: &HashMap<String, IntervalTree>) -> Option<RecordInfo> {
    if record.mapq() != 255 {
        return None;
    }
    let chrom_name = chrom_map.get(&(record.tid() as u32)).unwrap().to_string();
    let gene = get_aux_value(record, b"GX");
    let umi = get_aux_value(record, b"UB");
    let cell_barcode = get_aux_value(record, b"CB");
    let coord = get_coord(record);
    let start = record.pos() as u64;
    let end = record.cigar().end_pos() as u64;
    let pgene = if let Some(gene) = gene_tree.get(&chrom_name) {
        gene.search(start, end)
    } else {
        vec![]
    };
    let mut rec_info = RecordInfo::with_attr(chrom_name.clone(), gene.clone(), umi.clone(), cell_barcode.clone());
    rec_info.add_potential_gene(pgene);
    rec_info.add_coord(coord);
    return Some(rec_info);
}


fn add_rec2map(rec: RecordInfo, rec_map: &mut Arc<Mutex<HashMap<String, RecordInfo>>>) {
    let key = format!("{}-{}-{}", rec.umi, rec.cell_barcode, rec.chromosome);
    let mut map = rec_map.lock().unwrap();
    map.entry(key.clone()).and_modify(|full_info| {
        full_info.add_coord(rec.coordinates[0].clone());
        full_info.add_potential_gene(rec.potential_gene.keys().cloned().collect());
        full_info.set_gene(rec.gene.clone());
    }).or_insert(rec);
}

fn get_coord(record: &bam::Record) -> ReadCoord {
    let mut ranges: Vec<Vec<usize>> = Vec::new();
    let mut current_range: Vec<usize> = Vec::new();

    let cigar = record.cigar();
    let mut pos = record.pos() as usize;
    // bam is 0-based, but gtf is 1-based, so we need to add 1 to make it 1-based
    pos += 1;
    for op in &cigar {
        match op {
            bam::record::Cigar::Match(len) | bam::record::Cigar::Equal(len) | bam::record::Cigar::Diff(len) => {
                if current_range.is_empty() {
                    current_range.push(pos);
                }
                pos += *len as usize;
            }
            bam::record::Cigar::Del(len) => {
                // if !current_range.is_empty() {
                //     current_range.push(pos);
                // }
                pos += *len as usize;
            }
            bam::record::Cigar::RefSkip(len) => {
                if !current_range.is_empty() {
                    current_range.push(pos-1); // -1 to macth gtf
                    ranges.push(current_range);
                    current_range = Vec::new();
                }
                pos += *len as usize;
            }
            bam::record::Cigar::Ins(_len) | bam::record::Cigar::SoftClip(_len) | bam::record::Cigar::HardClip(_len) | bam::record::Cigar::Pad(_len) => {
                // 这些操作不消耗参考序列,不影响范围计算
            }
        }
    }

    if !current_range.is_empty() {
        current_range.push(pos);
        ranges.push(current_range);
    }
    // let qname = String::from_utf8_lossy(record.qname());
    // // println!("{:?}", qname.split("_").last().unwrap());
    // let qname = qname.split("_").last().unwrap().to_string();
    // 如果read是双端测序文件，并且read是在R1链上，根据是否含有TSS或者TES字样判断是否TSS还是TES
    let qname = if record.is_paired() && record.is_first_in_template() {
        let name = get_aux_value(record, b"RG");
        if name.contains("TES") {
            "TES".to_string()
        } else if name.contains("TSS") {
            "TSS".to_string()
        } else {
            "other".to_string()
        }
    } else {
        "other".to_string()
    };

    
    ReadCoord {
        coord: ranges.clone(),
        info: qname,
    }
}

fn recMap2df(rec_map: &RecordMap) -> Result<DataFrame, PolarsError> {
    let mut chr_vec: Vec<String> = Vec::new();
    let mut gene_vec: Vec<String> = Vec::new();
    let mut pgene_vec: Vec<String> = Vec::new();
    let mut cb_vec: Vec<String> = Vec::new();
    let mut umi_vec: Vec<String> = Vec::new();
    let mut coord_vec: Vec<String> = Vec::new();
    let mut all_trans: Vec<String> = Vec::new();
    let mut uniq_trans: Vec<String> = Vec::new();
    let mut select_trans: Vec<String> = Vec::new();

    for (_key, value) in rec_map.iter() {
        chr_vec.push(value.chromosome.clone());
        gene_vec.push(value.gene.clone());
        pgene_vec.push(value.pg_to_string());
        cb_vec.push(value.cell_barcode.clone());
        umi_vec.push(value.umi.clone());
        coord_vec.push(value.coord_to_string());
        all_trans.push(value.all_trans_to_string());
        uniq_trans.push(value.uniq_trans.clone());
        select_trans.push(value.select_trans.clone());
    }

    let df = DataFrame::new(vec![
        Series::new("chr", chr_vec),
        Series::new("gene", gene_vec),
        Series::new("potential_gene", pgene_vec),
        Series::new("cell_barcode", cb_vec),
        Series::new("umi", umi_vec),
        Series::new("select_trans", select_trans),
        Series::new("uniq_trans", uniq_trans),
        Series::new("freq", all_trans),
        Series::new("coord", coord_vec),
    ])?;

    Ok(df)
}

// fn add_transcript(recMap : &mut RecordMap, gene_map: &mut HashMap<String, Gene>) -> RecordMap {
//     let mut res_record_map = RecordMap::new();
//     for (key, umi) in recMap.into_iter() {
//         let gene_id = umi.select_gene();
//         if gene_id == "" {
//             res_record_map.insert(key.clone(), umi.clone());
//             continue;
//         }
//         // let gene = gene_map.get_mut(&gene_id).unwrap();
//         let gene = if let Some(gene) = gene_map.get_mut(&gene_id) {
//             gene
//         } else {
//             println!("gene_id not found: {}", gene_id);
//             continue;
//         };
//         let mut all_trans_id = vec![];
//         for read in umi.coordinates.iter() {
//             let mut trans_id_vec: Vec<String> = vec![];
//             for (_, trans) in gene.transcripts.iter() {
//                 if trans.exon_intervals.is_sub(read, 0, strand_to_symbol(&trans.strand)) {
//                     // println!("read: {:?}, trans: {:?}", read, trans.id);
//                     trans_id_vec.push(trans.id.clone());
//                 }
//             }
//             if trans_id_vec.len() > 0 {
//                 all_trans_id.push(trans_id_vec);
//             }
//         }
//         let mut freq_map: HashMap<String, usize> = HashMap::new();
//         for trans_vec in &all_trans_id {
//             for trans in trans_vec.clone() {
//                 let count = freq_map.entry(trans).or_insert(0);
//                 *count += 1;
//             }
//         }
//         umi.all_trans = freq_map.clone();
//         if all_trans_id.len() > 0 {
//             // find unique transcripts
//             let uniq_trans: Vec<String> = all_trans_id.clone()
//                 .into_iter()
//                 .fold(all_trans_id[0].clone(), |mut acc, vec| {
//                     acc.retain(|x| vec.contains(x));
//                     acc
//                 });
//             if uniq_trans.len() > 0 {
//                 umi.uniq_trans = uniq_trans.join("|");
//                 if uniq_trans.len() == 1 {
//                     umi.select_trans = uniq_trans[0].clone();
//                 } else {
//                     let max_trans = uniq_trans.into_iter()
//                         .max_by_key(|key| freq_map.get(key).unwrap_or(&0)).unwrap().clone();
//                     umi.select_trans = max_trans;
//                 }
//             } else {
//                 umi.uniq_trans = "NA".to_string();
//             }
//         } else {
//             umi.uniq_trans = "NA".to_string();
//         }
//         res_record_map.insert(key.clone(), umi.to_owned());
//     }
//     res_record_map
// }


fn add_transcript(recMap: RecordMap, gene_map: &mut HashMap<String, Gene>) -> RecordMap {
    let (sender, receiver): (_, Receiver<(_,_)>) = bounded(1024);
    let num_threads = 8;
    let pool = ThreadPool::new(num_threads);
    let mut res_recMap = Arc::new(Mutex::new(RecordMap::new()));
    pool.execute(move ||{
        for (key, umi) in recMap.into_iter() {
            sender.send((key, umi)).unwrap();
        }
        drop(sender);
    });
    for _ in 0..num_threads {
        let receiver = receiver.clone();
        let mut gene_map = gene_map.clone();
        let mut recMaps = res_recMap.clone();
        pool.execute(move || {
            for (key, value) in receiver.iter() {
                let mut value = value.clone();
                update(&key, &mut value, &mut recMaps, &mut gene_map);
            }
        });
    }
    pool.join();
    println!("{}", res_recMap.lock().unwrap().len());
    let mut res_recMaps = match Arc::try_unwrap(res_recMap) {
        Ok(mutex) => mutex.into_inner().unwrap(),
        Err(res_recMap) => {
            // 如果还有其他的 `Arc` 引用,`try_unwrap` 会失败并返回原始的 `Arc`
            // 在这种情况下,你可能需要重新考虑你的设计
            panic!("Failed to unwrap Arc!");
        }
    };
    return res_recMaps;
}

fn update(key: &String, umi: &mut RecordInfo, recMap: &mut Arc<Mutex<RecordMap>>, gene_map: &mut HashMap<String, Gene>) {
    let gene_id = umi.select_gene();
    if gene_id == "" {
        let mut map = recMap.lock().unwrap();
        map.insert(key.clone(), umi.clone());
        return;
    }
    let gene = if let Some(gene) = gene_map.get_mut(&gene_id) {
        gene
    } else {
        println!("gene_id not found: {}", gene_id);
        return;
    };
    let mut all_trans_id = vec![];
    for read in umi.coordinates.iter() {
        let mut trans_id_vec: Vec<String> = vec![];
        for (_, trans) in gene.transcripts.iter() {
            if trans.exon_intervals.is_sub(read, 0, strand_to_symbol(&trans.strand)) {
                // println!("read: {:?}, trans: {:?}", read, trans.id);
                trans_id_vec.push(trans.id.clone());
            }
        }
        if trans_id_vec.len() > 0 {
            all_trans_id.push(trans_id_vec);
        }
    }
    let mut freq_map: HashMap<String, usize> = HashMap::new();
    for trans_vec in &all_trans_id {
        for trans in trans_vec.clone() {
            let count = freq_map.entry(trans).or_insert(0);
            *count += 1;
        }
    }
    umi.all_trans = freq_map.clone();
    if all_trans_id.len() > 0 {
        // find unique transcripts
        let uniq_trans: Vec<String> = all_trans_id.clone()
            .into_iter()
            .fold(all_trans_id[0].clone(), |mut acc, vec| {
                acc.retain(|x| vec.contains(x));
                acc
            });
        if uniq_trans.len() > 0 {
            umi.uniq_trans = uniq_trans.join("|");
            if uniq_trans.len() == 1 {
                umi.select_trans = uniq_trans[0].clone();
            } else {
                let max_trans = uniq_trans.into_iter()
                    .max_by_key(|key| freq_map.get(key).unwrap_or(&0)).unwrap().clone();
                umi.select_trans = max_trans;
            }
        } else {
            umi.uniq_trans = "NA".to_string();
        }
    } else {
        umi.uniq_trans = "NA".to_string();
    }
    let mut map = recMap.lock().unwrap();
    map.insert(key.clone(), umi.to_owned());
    return;
}