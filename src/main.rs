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
        .get_matches();

    let myfile = matches.value_of("file").unwrap();


    let mut bam = IndexedReader::from_path(&myfile).unwrap();
    //bam.fetch((0, 10500, 10510)).unwrap(); 
    //MN908947.3

    for p in bam.pileup() {
        let pileup = p.unwrap();
        //println!("{}:{} depth {}", pileup.tid(), pileup.pos(), pileup.depth());
        let mut seq_string = String::from("Seq:");
        let mut foo = Vec::new();

        for alignment in pileup.alignments() {
            if !alignment.is_del() && !alignment.is_refskip() {
                //println!("{}", alignment.record().seq()[alignment.qpos().unwrap()] as char);
                //seq_string.push(alignment.record().seq()[aligSnment.qpos().unwrap()] as char);
            
                foo.push(alignment.record().seq()[alignment.qpos().unwrap()]);
                //println!("{}", alignment.record().seq());
            }
            else{
                foo.push(0);


            }

                        

        }
        let mut read_vector = Vec::new();
        let mut majority =0;
        let mut majority_base="N";

        let base_A = foo.iter().filter(|&n| *n == 65).count();
        read_vector.push(base_A);

        if base_A > majority {
            majority=base_A;
            majority_base = "A";
        }

        let base_T = foo.iter().filter(|&n| *n == 84).count();
        read_vector.push(base_T);

        if base_T > majority {
            majority=base_T;
            majority_base = "T";
        }
        
        let base_C = foo.iter().filter(|&n| *n == 67).count();
        read_vector.push(base_C);

        if base_C > majority {
            majority=base_C;
            majority_base = "C";
        }

        let base_G = foo.iter().filter(|&n| *n ==71).count();
        read_vector.push(base_G);

        if base_G > majority {
            majority=base_G;
            majority_base = "G";
        }
 
        let base_D= foo.iter().filter(|&n| *n == 0).count();
        read_vector.push(base_D);

        if base_D > majority {
            majority=base_D;
            majority_base = "D";
        }


        //let noise: f32 =  (pileup.depth() as u16 - majority as u16) / pileup.depth() as u16;
        let noise: f64 =  (pileup.depth() as u16 - majority as u16).into();
        let noise_t= noise / pileup.depth() as f64;



        println!("{}\t{}\t{}\t{}", pileup.pos(),noise_t,pileup.depth() as u32,majority_base);
        //println!("Total {}", pileup.depth() as u16);
        //println!("Noise {:?}", noise);
        //println!("Noise_t {:?}", noise_t);
        //println!("MR {:?}", majority_base);
        

    }

    //println!("Elapsed time: {:.2?}", before.elapsed());

}