use clap::{Arg, App};
use rust_htslib::{bam, bam::Read};
use std::sync::Arc;
use std::thread;

fn main() {
    let matches = App::new("NoisExtractor")
        .version("0.1.0")
        .author("Nacho Garcia <iggl@fhi.no>")
        .about("Noise Extraction of BAM Files")
        .arg(Arg::with_name("file")
                 .short("f")
                 .long("file")
                 .takes_value(true)
                 .help("BAM file"))
        .get_matches();

    let myfile = matches.value_of("file");
    println!("Hello world my friend!");

    

}