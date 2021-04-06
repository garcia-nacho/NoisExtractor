use rust_htslib::bam::{IndexedReader, Read, Record};
use clap::{Arg, App};
use std::sync::Arc;

//use std::time::Instant;


fn main() {
    //let before = Instant::now();
    let matches = App::new("NoisExtractor")
        .version("0.1.0")
        .author("Nacho Garcia <iggl@fhi.no>")
        .about("Noise Extraction of BAM Files")
        .arg(Arg::with_name("file")
                 .short("f")
                 .long("file")
                 .takes_value(true)
                 .help("BAM file"))
        .arg(Arg::with_name("cutoff")
                .short("c")
                .long("cutoff")
                .takes_value(true)
                .help("Cut-off for calling secondary sequence"))
        .arg(Arg::with_name("output")
                .short("o")
                .long("out")
                .takes_value(true)
                .help("Output files to save"))
        .get_matches();

    let myfile = matches.value_of("file").unwrap();
//cutoff for calling the secondary sequence
    let c_o = matches.value_of("cutoff");
    let mut cutoff= 0.0;

    match c_o {
        None => cutoff = 0.0,
        Some(s) => {
            match s.parse::<f64>() {
                Ok(c) => cutoff = c,
                Err(_) => cutoff = 0.0,
            }
        }
    }
//println!("cutoff{}", cutoff);
//Output files
    let output_file = matches.value_of("output").unwrap_or("none");

    let mut bam = IndexedReader::from_path(&myfile).unwrap();
    //bam.fetch((0, 10500, 10510)).unwrap(); 
    //MN908947.3

    for p in bam.pileup() {
        let pileup = p.unwrap();

        let mut total_base = Vec::new();

        for alignment in pileup.alignments() {
            if !alignment.is_del() && !alignment.is_refskip() {

            
                total_base.push(alignment.record().seq()[alignment.qpos().unwrap()]);

            }
            else{
                total_base.push(0);
            }
                    
        }

        let mut read_vector = Vec::new();
        let mut majority =0;
        let mut majority2 =0;
        let mut majority_base="N";
        let mut majority_base2="N";
        let bases = vec!["A", "T", "C","G","D"];

        // Base A
        let base_A = total_base.iter().filter(|&n| *n == 65).count();
        read_vector.push(base_A);

        if base_A > majority {
            majority=base_A;
            majority_base = "A";
        }

        // Base T
        let base_T = total_base.iter().filter(|&n| *n == 84).count();
        read_vector.push(base_T);

        if base_T > majority {
            majority=base_T;
            majority_base = "T";    
        }

        //Base C
        let base_C = total_base.iter().filter(|&n| *n == 67).count();
        read_vector.push(base_C);

        if base_C > majority {
            majority=base_C;
            majority_base = "C";       
        }


        //Base G
        let base_G = total_base.iter().filter(|&n| *n ==71).count();
        read_vector.push(base_G);

        if base_G > majority {
            majority=base_G;
            majority_base = "G";
        }

        // Deletions
        let base_D= total_base.iter().filter(|&n| *n == 0).count();
        read_vector.push(base_D);

        if base_D > majority {
            majority=base_D;
            majority_base = "D";
        }

        //Sequence 2 calculation
        let permutation = permutation::sort(&read_vector[..]);
        let ordered_bases = permutation.apply_slice(&bases[..]);
        let ordered_reads = permutation.apply_slice(&read_vector[..]);
        let mut freq_mr2: f64 = 0.0;

        if ordered_reads[3]==0 {
            majority_base2 = ordered_bases[4];
        }else{
            majority_base2 = ordered_bases[3];
            majority2= ordered_reads[3];
            freq_mr2 = majority2 as f64 / pileup.depth() as f64;
            
            if freq_mr2 < cutoff {
                majority_base2= ordered_bases[4];
            }
        }
        
        //Noise calculation
        let noise: f64 =  (pileup.depth() as u16 - majority as u16).into();
        let noise_t= noise / pileup.depth() as f64;

        println!("{}\t{}\t{}\t{}\t{}\t{}", pileup.pos(),noise_t,pileup.depth() as u32,majority_base, majority_base2, freq_mr2);
        //println!("Total {}", pileup.depth() as u16);

        //println!("{:?}", ordered_bases);
        

    }

    //println!("Elapsed time: {:.2?}", before.elapsed());

}