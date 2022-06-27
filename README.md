# NoiseExtractor (aka. FINex)   
[![Linux](https://svgshare.com/i/Zhy.svg)](https://svgshare.com/i/Zhy.svg)   [![made-with-rust](https://img.shields.io/badge/Made%20with-Rust-1f425f.svg)](https://www.rust-lang.org/)   
      
Retrieve the noise for all the positions of aligened and sorted reads

# Installing NoisExtractor:   
Compile the binary from the source and copy the binary (noise_extractor) to your favourite location.    
See [here](https://rustc-dev-guide.rust-lang.org/building/how-to-build-and-run.html) about how to compile Rust source code.   

# Running NoiseExtractor:   
<code>noise_extractor -f Bam.file.bam -c 0.15 -o output > output.tsv</code>   
Where *-f* is the indexed bam file you want to explore. *-c* is the cutoff under which a position is considered "noisy" (Here we call noise to the sum of the frequencies of bases that are not the one consider the major variant). *-o* is the prefix of the fasta files generated in the output (considering the cutoff defined on *-c*). *-o* and *-c* are not required and if not provided NoiseExtractor will not generate fasta output files (Actually it will generate them anyway, but they will be deleted afterwards)
   
# Output:   
You will get a tsv with 6 columns: Position, Noise, Depth, Major Variant, Minor Variant, Frequency of minor variant. Note that the deletions and insertions are exported as "D" and "I" respectively in the output file.   
Additionally. NoisExtractor will generate two fasta files for the major and minor sequence (Only if the flags *-c* and *-o* are used) 

# Limitations: 
Noise_extractor can find insertions but they are not inserted in the fasta files. 
