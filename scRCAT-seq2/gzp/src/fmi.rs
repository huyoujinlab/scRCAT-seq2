use bio::alphabets::dna::{self, alphabet};
use bio::data_structures::bwt::{bwt, less, Occ};
use bio::data_structures::fmindex::{BackwardSearchResult, FMIndex, FMIndexable};
use bio::data_structures::suffix_array::suffix_array;
use bio::io::fastq::{self, Record};
use bio::pattern_matching::ukkonen::{self, unit_cost, Ukkonen};
use regex::Regex;
use bio::alphabets::Alphabet;
use strsim::hamming;

pub fn searchTag(record: &Record, pattern1: &[u8], alphabet: &Alphabet) -> Vec<usize> {
    let mut sequence = record.seq().to_vec();
    sequence.push(b'$');
    let sa = suffix_array(&sequence);
    let bwt = bwt(&sequence, &sa);
    let less = less(&bwt, alphabet);
    let occ = Occ::new(&bwt, 3, alphabet); // Using a small sampling rate for Occ table
    let fmindex = FMIndex::new(&bwt, &less, &occ);

    // let positions1 = fmindex_search(&fmindex, &sa, pattern1);
    // let positions2 = fmindex_search(&fmindex, &sa, pattern2);
    let search1 = fmindex.backward_search(pattern1.iter());
    let positions1 = match search1 {
        BackwardSearchResult::Complete(saint) => {
            saint.occ(&sa)
        },
        BackwardSearchResult::Partial(saint, l) => {
            // println!("Partial {:?} {:?}", saint.occ(&sa), l);
            for pos in saint.occ(&sa) {
                if pos + l >= sequence.len() {
                    let mut sequence = record.seq().to_vec();
                    let mut sequence2 = sequence[pos..].to_vec();
                    print!("Partial: raw: {}, sub: {}, length: {}", String::from_utf8_lossy(&sequence), String::from_utf8_lossy(&sequence2), l);
                }
            }
            Vec::new()
        },
        BackwardSearchResult::Absent => {
            Vec::new()
        },
    };
    positions1
}

pub fn search_sequence(record: &Record, pattern1: &[u8], pattern2: &[u8], gap: usize, _max_mismatch: u8, alphabet: &Alphabet) -> Vec<usize> {
    let mut sequence = record.seq().to_vec();
    sequence.push(b'$');
    let sa = suffix_array(&sequence);
    let bwt = bwt(&sequence, &sa);
    let less = less(&bwt, alphabet);
    let occ = Occ::new(&bwt, 3, alphabet); // Using a small sampling rate for Occ table
    let fmindex = FMIndex::new(&bwt, &less, &occ);

    // let positions1 = fmindex_search(&fmindex, &sa, pattern1);
    // let positions2 = fmindex_search(&fmindex, &sa, pattern2);
    let search1 = fmindex.backward_search(pattern1.iter());
    let search2 = fmindex.backward_search(pattern2.iter());
    let positions1 = match search1 {
        BackwardSearchResult::Complete(saint) => {
            saint.occ(&sa)
        },
        BackwardSearchResult::Partial(saint, _) => {
            Vec::new()
        },
        BackwardSearchResult::Absent => {
            Vec::new()
        },
    };
    let positions2 = match search2 {
        BackwardSearchResult::Complete(saint) => {
            saint.occ(&sa)
        },
        BackwardSearchResult::Partial(saint, _) => {
            Vec::new()
        },
        BackwardSearchResult::Absent => {
            Vec::new()
        },
    };

    let mut locations = Vec::new();
    if positions1.is_empty() || positions2.is_empty() {
        return locations;
    }
    // println!("id: {}, positions1: {:?}, positions2: {:?}", record.id(), positions1, positions2);
    for &pos1 in &positions1 {
        for &pos2 in &positions2 {
            if pos2 > pos1 && pos2 - pos1 - pattern1.len() == gap {
                locations.push(pos1);
                break; // Only add the first match for each location in pattern1
            }
        }
    }

    locations
}


pub fn search_ukkonen(record: &Record, pattern1: &[u8], max_mismatch: usize) -> Vec<usize> {
    // revised it to return the position of the first character of the pattern
    let mut ukkonen = Ukkonen::with_capacity(150, unit_cost);
    let occ: Vec<(usize, usize)> = ukkonen.find_all_end(pattern1, record.seq(), max_mismatch).collect();
    let seqs = record.seq();
    let mut res = Vec::new();
    for (mut pos, _) in &occ {
        // println!("pos: {}", pos);
        pos = pos + 1;
        if (pos < pattern1.len()) {
            continue;
        }
        let matched = &seqs[pos - pattern1.len()..pos];
        // let Ok(ss) = compare_strings(&String::from_utf8_lossy(matched), &String::from_utf8_lossy(pattern1), max_mismatch);
        // if  ss {
        //     println!("pos: {}, seq: {}, matched: {}, pattern: {}", pos, String::from_utf8_lossy(&seqs), String::from_utf8_lossy(matched), String::from_utf8_lossy(pattern1));
        //     res.push(pos - pattern1.len());
        // }
        let ss = hamming(&String::from_utf8_lossy(matched), &String::from_utf8_lossy(pattern1));
        if let Ok(mismatch_count) = ss {
            if mismatch_count <= max_mismatch {
                // println!("pos: {}, seq: {}, matched: {}, pattern: {}", pos, String::from_utf8_lossy(&seqs), String::from_utf8_lossy(matched), String::from_utf8_lossy(pattern1));
                res.push(pos - pattern1.len());
            }
        }
    }
    res
}


pub fn search_ukkonen_two(record: &Record, pattern1: &[u8], pattern2: &[u8], gap: usize, max_mismatch: usize, verbose: bool) -> Vec<usize> {
    let occ1 = search_ukkonen(record, pattern1, max_mismatch);
    let seqs = record.seq();
    let (mut pos1, pos2) = (0, 0);
    let mut res = Vec::new();
    if pattern2.len() > 0 {
        let occ2 = search_ukkonen(record, pattern2, max_mismatch);
        for pos1 in &occ1 {
            for pos2 in &occ2 {
                // if verbose {
                //     println!("pos1: {}, pos2: {}, seq: {}", pos1, pos2, String::from_utf8_lossy(&seqs));
                // }
                if pos2 > pos1 && pos2 - pos1 - pattern1.len() == gap {
                    // res.push(*pos1 - pattern1.len());
                    res.push(*pos1);
                    break; // Only add the first match for each location in pattern1
                }
            }
        }
        return res;
    } else {
        for pos1 in &occ1 {
            res.push(*pos1);
        }
        return res;
    }
}


pub fn find_t_positions(sequence: &str) -> Vec<(usize, usize)> {
    let re = Regex::new(r"T{5,}(?:(?:[ACGN]T{4,}){1,2})?").unwrap();
    let mut positions = Vec::new();

    for mat in re.find_iter(sequence) {
        let matched_str = &sequence[mat.start()..mat.end()];
        let t_count = matched_str.chars().filter(|&c| c == 'T').count();
        if t_count >= 15 {
            positions.push((mat.start(), mat.end()));
        }
    }

    positions
}

pub fn find_a_positions(sequence: &str) -> Vec<(usize, usize)> {
    let re = Regex::new(r"A{5,}(?:(?:[CGNT]A{4,}){1,2})?").unwrap();
    let mut positions = Vec::new();

    for mat in re.find_iter(sequence) {
        let matched_str = &sequence[mat.start()..mat.end()];
        let a_count = matched_str.chars().filter(|&c| c == 'A').count();
        if a_count >= 15 {
            positions.push((mat.start(), mat.end()));
        }
    }

    positions
}


fn trim_read1(record: &Record, pattern3: &[u8], alphabet: &Alphabet, max_mismatches: usize, remove_polya: bool, verbose: bool) -> fastq::Record {
    if pattern3.len() > 0  && remove_polya == false {
        let rev_pattern3 = dna::revcomp(pattern3);
        let loc_tso_r1 = searchTag(record, &rev_pattern3, alphabet);
        if loc_tso_r1.len() > 0 && loc_tso_r1[0] > 20 {
            return fastq::Record::with_attrs(record.id(),
                record.desc(), 
                &record.seq()[..loc_tso_r1[0]].to_owned(),
                &record.qual()[..loc_tso_r1[0]].to_owned());
        } else {
            return record.clone();
        }
    } else if pattern3.len() == 0 && remove_polya {
        let loc_tso_r1 = find_a_positions(std::str::from_utf8(record.seq()).unwrap());
        if loc_tso_r1.len() > 0 && loc_tso_r1[0].0 > 20 {
            return fastq::Record::with_attrs(record.id(),
                record.desc(), 
                &record.seq()[..loc_tso_r1[0].0].to_owned(),
                &record.qual()[..loc_tso_r1[0].0].to_owned());
        } else {
            return record.clone();
        }
    } else {
        println!("pattern 3 and remove_polya are mutually exclusive, please provide only one");
        return fastq::Record::new();
    }
}

fn trim_read2(record: &Record, pattern1: &[u8], pattern2: &[u8], pattern3: &[u8], umi_len: usize, max_mismatches: usize, remove_polya: bool, only_cut_adpter: bool, verbose: bool) -> fastq::Record {
    // if verbose {
    //     println!("pattern1: {:?}, pattern2: {:?}, pattern3: {:?}", );
    // }
    if only_cut_adpter {
        let loc_umi = search_ukkonen_two(record, pattern1, pattern2, umi_len, max_mismatches, verbose);
        if verbose {
            println!("pattern1: {:?}, pattern2: {:?}, pattern3: {:?}\nloc_umi: {:?}, rs_id: {:?}", String::from_utf8_lossy(pattern1), String::from_utf8_lossy(pattern2), String::from_utf8_lossy(pattern3), loc_umi, record.id());
        }
        if loc_umi.len() > 0 {
            return fastq::Record::with_attrs(record.id(),
                record.desc(), 
                &record.seq()[loc_umi[0]..].to_owned(),
                &record.qual()[loc_umi[0]..].to_owned());
        } else {
            return fastq::Record::new();
        }
    } else {
        let loc_umi = search_ukkonen_two(record, pattern1, pattern2, umi_len, max_mismatches, verbose);
        if remove_polya {
            let loc_tso_r2 = find_t_positions(std::str::from_utf8(record.seq()).unwrap());
            let mut cb_records = fastq::Record::new();
            let mut cDNA_record = fastq::Record::new();
            if verbose {
                println!("remove_polya, Read2: loc_umi: {:?}, loc_tso_r2: {:?}, id: {:?}, seq: {:?}", loc_umi, loc_tso_r2, record.id(), String::from_utf8_lossy(record.seq()));
            }
            if loc_umi.len() > 0 && loc_tso_r2.len() > 0 {
                cb_records = fastq::Record::with_attrs(record.id(),
                    record.desc(), 
                    &record.seq()[loc_umi[0]..loc_umi[0] + pattern1.len() + umi_len + pattern2.len()].to_owned(),
                    &record.qual()[loc_umi[0]..loc_umi[0] + pattern1.len() + umi_len + pattern2.len()].to_owned());
                cDNA_record = fastq::Record::with_attrs(record.id(),
                    record.desc(), 
                    &record.seq()[loc_tso_r2[0].1..],
                    &record.qual()[loc_tso_r2[0].1..]);
                return fastq::Record::with_attrs(record.id(),
                    record.desc(),
                    &[cb_records.seq(), cDNA_record.seq()].concat(),
                    &[cb_records.qual(), cDNA_record.qual()].concat());
            } else if loc_umi.len() > 0 && loc_tso_r2.len() == 0 {
                return fastq::Record::with_attrs(record.id(),
                    record.desc(), 
                    &record.seq()[loc_umi[0]..].to_owned(),
                    &record.qual()[loc_umi[0]..].to_owned());
            } else {
                return fastq::Record::new();
            }
        } else {
            let loc_tso_r2 = search_ukkonen(record, pattern3, max_mismatches);
            let mut cb_records = fastq::Record::new();
            let mut cDNA_record = fastq::Record::new();
            if verbose {
                println!("Read2: loc_umi: {:?}, loc_tso_r2: {:?}, id: {:?}, seq: {:?}", loc_umi, loc_tso_r2, record.id(), String::from_utf8_lossy(record.seq()));
            }
            if loc_umi.len() > 0 && loc_tso_r2.len() > 0 {
                cb_records = fastq::Record::with_attrs(record.id(),
                    record.desc(), 
                    &record.seq()[loc_umi[0]..loc_umi[0] + pattern1.len() + umi_len + pattern2.len()].to_owned(),
                    &record.qual()[loc_umi[0]..loc_umi[0] + pattern1.len() + umi_len + pattern2.len()].to_owned());
                cDNA_record = fastq::Record::with_attrs(record.id(),
                    record.desc(), 
                    &record.seq()[loc_tso_r2[0]..],
                    &record.qual()[loc_tso_r2[0]..]);
                return fastq::Record::with_attrs(record.id(),
                    record.desc(),
                    &[cb_records.seq(), cDNA_record.seq()].concat(),
                    &[cb_records.qual(), cDNA_record.qual()].concat());
            } else if loc_umi.len() > 0 && loc_tso_r2.len() == 0 {
                return fastq::Record::with_attrs(record.id(),
                    record.desc(), 
                    &record.seq()[loc_umi[0]..].to_owned(),
                    &record.qual()[loc_umi[0]..].to_owned());
            } else {
                return fastq::Record::new();
            }
        }
    }
}

fn trim_read2_reverse(record: &Record, rev_record: &Record, pattern1: &[u8], pattern2: &[u8], pattern3: &[u8], umi_len: usize, max_mismatches: usize, remove_polya: bool, only_cut_adpter: bool, verbose: bool) -> fastq::Record {
    if only_cut_adpter {
        let loc_umi = search_ukkonen_two(rev_record, pattern1, pattern2, umi_len, max_mismatches, verbose);
        let pattern1_rev = dna::revcomp(pattern1);
        let pattern2_rev = dna::revcomp(pattern2);
        let loc_umi_revcmp = search_ukkonen_two(record, &pattern1_rev, &pattern2_rev, umi_len, max_mismatches, verbose);
        if verbose {
            println!("loc_umi: {:?}, loc_umi_revcmp: {:?}", loc_umi, loc_umi_revcmp);
        }
        if loc_umi.len() > 0 && loc_umi_revcmp.len() > 0 {
            let cb_records = fastq::Record::with_attrs(rev_record.id(),
                rev_record.desc(), 
                &rev_record.seq()[loc_umi[0]..loc_umi[0] + pattern1.len() + umi_len + pattern2.len()].to_owned(),
                &rev_record.qual()[loc_umi[0]..loc_umi[0] + pattern1.len() + umi_len + pattern2.len()].to_owned());
            let cDNA_record = fastq::Record::with_attrs(record.id(),
                record.desc(), 
                &record.seq()[loc_umi_revcmp[0] + pattern1.len()..],
                &record.qual()[loc_umi_revcmp[0] + pattern1.len()..]);
            if verbose {
                println!("cb_records: {:?}, cDNA_record: {:?}", cb_records, cDNA_record);
            }
            return fastq::Record::with_attrs(rev_record.id(),
                rev_record.desc(),
                &[cb_records.seq(), cDNA_record.seq()].concat(),
                &[cb_records.qual(), cDNA_record.qual()].concat());
        } else {
            return fastq::Record::new();
        }
    } else {
        let loc_umi = search_ukkonen_two(rev_record, pattern1, pattern2, umi_len, max_mismatches, verbose);
        let loc_tso_r2 = search_ukkonen(record, pattern3, max_mismatches);
        let pattern1_rev = dna::revcomp(pattern1);
        let loc_umi_rev = search_ukkonen(record, &pattern1_rev, max_mismatches);
        // if record.id() == "V350233155L2C002R02100781353:0:0:0:0" {
        //     println!("Read2: loc_umi: {:?}, loc_umi_rev: {:?} loc_tso_r2: {:?}, id: {:?}, seq: {:?}", loc_umi, loc_umi_rev, loc_tso_r2, record.id(), String::from_utf8_lossy(record.seq()));
        // }
        let loc_tso_r2 = if loc_umi_rev.len() > 0 && loc_tso_r2.len() > 1 {
            let mut res = vec!{};
            for x in &loc_tso_r2 {
                if *x > loc_umi_rev[0] {
                    res.push(*x);
                }
            }
            // println!("{:?}, {:?}, {:?}", loc_tso_r2, res, record);
            res
        } else {
            search_ukkonen(record, pattern3, max_mismatches)
        };

        if verbose {
            println!("Read2: loc_umi: {:?}, loc_tso_r2: {:?}, id: {:?}, seq: {:?}", loc_umi, loc_tso_r2, record.id(), String::from_utf8_lossy(record.seq()));
        }
        if loc_umi.len() > 0 && loc_tso_r2.len() > 0 {
            // if the TSO is found before the UMI, return empty record
            if loc_tso_r2[0] < loc_umi_rev[0] {
                return fastq::Record::new();
            }
            let umi_end = loc_umi[0] + pattern1.len() + umi_len + pattern2.len();
            if umi_end > rev_record.seq().len() {
                return fastq::Record::new();
            }
            
            let cb_records = fastq::Record::with_attrs(
                rev_record.id(),
                rev_record.desc(),
                &rev_record.seq()[loc_umi[0]..umi_end].to_owned(),
                &rev_record.qual()[loc_umi[0]..umi_end].to_owned(),
            );
            let cDNA_record = fastq::Record::with_attrs(
                rev_record.id(),
                rev_record.desc(),
                &record.seq()[loc_tso_r2[0]..],
                &record.qual()[loc_tso_r2[0]..],
            );
            return fastq::Record::with_attrs(
                rev_record.id(),
                rev_record.desc(),
                &[cb_records.seq(), cDNA_record.seq()].concat(),
                &[cb_records.qual(), cDNA_record.qual()].concat(),
            );
        } else {
            return fastq::Record::new();
        }
    }
}


pub fn process_record(record1: &fastq::Record, record2: &fastq::Record, pattern1: &[u8], pattern2: &[u8], pattern3: &[u8], alphabet: &Alphabet, umi_len: usize, max_mismatches: usize, remove_polya: bool, reverse: bool, only_cut_adpter: bool, verbose: bool) -> (Option<fastq::Record>, Option<fastq::Record>) {
    // pattern 1 is the UMI prefix pattern, pattern 2 is the UMI suffix pattern, pattern 3 is the TSO pattern
    if reverse {
        let rev_record2 = fastq::Record::with_attrs(record2.id(),
            record2.desc(), 
            &dna::revcomp(record2.seq()),
            &record2.qual().iter().rev().cloned().collect::<Vec<u8>>());
        // return process_record(record1, &rev_record2, pattern1, pattern2, pattern3, alphabet, umi_len, max_mismatches, remove_polya, false, only_cut_adpter, verbose);
        let trimmed_record1 = trim_read1(record1, pattern3, alphabet, max_mismatches, remove_polya, verbose);
        let trimmed_record2 = trim_read2_reverse(record2, &rev_record2, pattern1, pattern2, pattern3, umi_len, max_mismatches, remove_polya, only_cut_adpter, verbose);
        if trimmed_record1.is_empty() || trimmed_record2.is_empty() {
            return (None, None);
        }
        return (Some(trimmed_record1), Some(trimmed_record2));
    }
    let trimmed_record1 = trim_read1(record1, pattern3, alphabet, max_mismatches, remove_polya, verbose);
    let trimmed_record2 = trim_read2(record2, pattern1, pattern2, pattern3, umi_len, max_mismatches, remove_polya, only_cut_adpter, verbose);
    if trimmed_record1.is_empty() || trimmed_record2.is_empty() {
        return (None, None);
    }
    return (Some(trimmed_record1), Some(trimmed_record2));
}
 