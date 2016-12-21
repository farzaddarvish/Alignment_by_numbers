# Alignment_by_numbers

## Description
Alignment_by_numbers is a prototype reference-based nucleotide sequence aligner that created to evaluate the
proposed methodology in the manuscript 'Alignment by numbers: sequence assembly using compressed numerical representations'.

## Installation (Linux, Mac OS X):
- Download repository as a zip or with `git clone https://github.com/Avramis/Alignment_by_numbers.git`  
- `cd Code`  
- `clang++ -Wall -O3 -o alignment_by_numbers *.cpp *.c *.cc` (compile)  
- `cd ../Exec`  
- `./Alignment_by_numbers`  

## Examples
```
Usage:
    Parameters:
        -fa [~/reference_genome.fasta]    -- (Required parameter) fasta file directory and file name
        -fq [~/reads_file.fastq]          -- (Required parameter) fastq file directory and file name
        -o [~/outputdirectory]            -- (Optional parameter) output directory for sam file. If not provided the sam file is save in the same directory as the fastq file.
        -kmer [x]                         -- (Optional parameter) kmer size for the analysis [default: 100]
        -rep [Representation method]      -- (Optional parameter) Representation method for converting nucleotide sequences to numerical sequences [default: Voss_indicators]
        Options:
            01) Atomic_Numbers
            02) Complex_Numbers
            03) Dna_Walk
            04) EIIP_numbers
            05) Integer_numbers
            06) Pair_numbers
            07) Real_numbers
            08) Tetrahedron
            09) Voss_indicators   --   default option
            10) Z_curve
        -tra [Transformation method]      -- (Optional parameter) Transformation method for compressing data [default: DFT]
        Options:
            01) DFT   --   default option
            02) DWT
            03) PAA
        -clvl [y]                         -- (Optional parameter) Compression level to apply [default: 2]
        -knn [kn]                         -- (Optional parameter)the number of kn neighbor data to use for the KNN search [default: 100]
	-s [true/false]                   -- (Optional parameter) run a sensitive search [default: fasle]
		-bs [x]                           -- (Optional parameter) set the size of the tree [default: 100000]
	-rd [~/report output directory]   -- (Optional parameter) output directory to save the timing report file. If not provided the report text file is save in the same directory as the fastq file.
       -h                                -- Print Help Options


    Example:
        01) Alignment_by_numbers -fa ~/genome.fasta -fq ~/reads.fastq
        02) Alignment_by_numbers -fa ~/genome.fasta -fq ~/reads.fastq -kmer 300 -rep Tetrahedron -tra DWT -clvl 4 -knn 30
        03) Alignment_by_numbers -fa ~/genome.fasta -fq ~/reads.fastq -kmer 300 -rep Tetrahedron -tra DWT -o ~/samfiledirectory -clvl 4 -knn 30 -s true -bs 10000 -rd ~/reportsfile
```
