use noodles::bam::{self, record, Record};
use std::collections::HashMap;


fn main() {
    let mut reader = bam::io::reader::Builder::default().build_from_path("/home/zjw/zjw/20210331circle/script_and_log/20221026_circle_0.2-0.6-2-splint-gradient/c0.2/ES-1/hsa/isorestructure/preprocess/ES.filtered.Aligned.GeneTagged.UBfix.coordinateSorted_unique.bam").unwrap();
    let header = reader.read_header().unwrap();

    let mut i = 0;
    for result in reader.records() {
        let record = result.unwrap();
        println!("{:?}", record.alignment_start());
        i += 1;
        if i > 10 {
            break;
        }
        
    }
}
