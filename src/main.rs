/* 
FHINEx 
FHI Noise Extractor

Licenced under GPL v3 licence

Ignacio Garcia Llorente 2021 
iggl@fhi.no
Folkehelseinstitute 
Oslo-Norway

*/

//use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::{bam, bam::Read, bam::IndexedReader};
use clap::{Arg, App};
use std::path::Path;
use std::fs::File;
use std::io::prelude::*;
use std::iter;
use regex::Regex;
use std::fs;

//For speed testing purposes uncomment down
//use std::time::Instant; 

fn main() {
    //For speed testing purposes uncomment down
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
        .arg(Arg::with_name("del_mode")
                .short("d")
                .long("del_mode")
                .takes_value(false)
                .help("Mask deletions"))
        .get_matches();

    //matches.is_present("dyn_mode")

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
    let output_file = matches.value_of("output").unwrap_or("temp");
    let name = &myfile;
    let name = name.replace(".bam", "");
    let re = Regex::new(r".*/").unwrap();

    let name = re.replace(&name, "");

    let mut header_mr: String = ">".to_owned();
    let mut header_ct: String = ">".to_owned();
    
    header_mr.push_str(&name);
    header_ct.push_str(&name);

    header_mr.push_str("_MajorityRule\r");
    header_ct.push_str("_Contaminant\r");

    //Consensus
        let mut output_fa = String::from(output_file);
        output_fa.push_str("_MR.fa");
        //println!("{}", output_fa);

        let path_fa = Path::new(&output_fa);
        let display_fa = path_fa.display();
        let mut file_1 = match File::create(&path_fa) {
            Err(why) => panic!("couldn't create {}: {}", display_fa, why),
            Ok(file_1) => file_1,
        };
        match file_1.write_all(&header_mr.as_bytes()) {
            Err(why) => panic!("couldn't write to {}: {}", display_fa, why),
            Ok (_) => (),
        }

    //Contaminant
        let mut output_ct = String::from(output_file);
        output_ct.push_str("_contaminant.fa");
        //println!("{}", output_ct);

        let path_ct = Path::new(&output_ct);
        let display_ct = path_ct.display();
        let mut file_2 = match File::create(&path_ct) {
            Err(why) => panic!("couldn't create {}: {}", display_ct, why),
            Ok(file_2) => file_2,
        };
        match file_2.write_all(&header_ct.as_bytes()) {
            Err(why) => panic!("couldn't write to {}: {}", display_ct, why),
            Ok(_) => (),
        }

    let mut bam = IndexedReader::from_path(&myfile).unwrap();
    //bam.fetch((0, 10500, 10510)).unwrap(); 

    let mut index=1;

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
              
            //Insertion
            // * fetch(FetchDefinition::Unmapped) or fetch("*") -> Fetch unmapped (as signified by the 'unmapped' flag in the BAM - might be unreliable with some aligners.
            
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

        println!("{}\t{}\t{}\t{}\t{}\t{}", pileup.pos()+1,noise_t,pileup.depth() as u32,majority_base, majority_base2, freq_mr2);
        
        let pad_n: u32 =&pileup.pos()+1-&index;
        //println!("{}",pad_n);
        if output_file != "temp"{

            if pad_n >1{
            
                let pad: String = iter::repeat("N").take(pad_n as usize).collect();

                match file_2.write_all(pad.as_bytes()) {
                    Err(why) => panic!("couldn't write to {}",  why),
                    Ok(_) => (),
                }
    
                match file_1.write_all(pad.as_bytes()) {
                    Err(why) => panic!("couldn't write to {}",  why),
                    Ok(_) => (),
                }

        }
            //Needed a loop to fill gaps here plus if

            
            if matches.is_present("del_mode") && majority_base2=="D" {
                majority_base2="";
            }
            if matches.is_present("del_mode") && majority_base=="D" {
                majority_base="";
            }
            

            match file_2.write_all(majority_base2.as_bytes()) {
                Err(why) => panic!("couldn't write to {}",  why),
                Ok(_) => (),
            }

            match file_1.write_all(majority_base.as_bytes()) {
                Err(why) => panic!("couldn't write to {}",  why),
                Ok(_) => (),
            }
        }
        
        index =  pileup.pos()+1;
    }
    //Delete temp files
    if output_file == "temp"{

        fs::remove_file(&output_fa).unwrap();
        fs::remove_file(&output_ct).unwrap();

    }

    //println!("Elapsed time: {:.2?}", before.elapsed());

}