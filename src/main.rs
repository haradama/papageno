use bio::bio_types::sequence::SequenceRead;
use bio::io::fastq;
use std::collections::HashMap;
use std::env;
use std::path::Path;

fn phred_to_p_error(phred: f64) -> f64 {
    10f64.powf(-phred / 10.0)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let path = Path::new(&args[1]);

    let offset = HashMap::from([
        ("Sanger", 33),
        ("Solexa", 64),
        ("Illumina_1.3+", 64),
        ("Illumina_1.5+", 64),
        ("Illumina_1.8+", 33),
    ]);

    let reader: fastq::Reader<std::io::BufReader<std::fs::File>> = fastq::Reader::from_file(&path).unwrap();

    for result in reader.records() {
        let record: fastq::Record = result.unwrap();
        let quality: &[u8] = record.qual();

        println!("{} {}", record.id(), record.len());
        let mut cumulative_expected_error = 0.0;
        let mut v: Vec<f64> = Vec::new();
        for phred_score in quality {
            let p_error = phred_to_p_error(f64::from(*phred_score - offset["Sanger"]));
            cumulative_expected_error += p_error;
            v.push(cumulative_expected_error);
        }
        println!("{}", cumulative_expected_error);

        // println!(":? = {:?}", quality);
    }
}
