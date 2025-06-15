# blast2counts
This repository contains scripts to analyze RNA-seq data with BLASTn.

## analyze_susp_fastqs.py
This script allows the extraction of valid reads from a given FASTQ file. Reads with long N stretches are ignored.


```
Usage:
  python3 analyze_susp_fastqs.py --in <FASTQ_GZ_FILE> --out <FASTA_FILE>
  --in      STR    Input FASTQ file (.fastq.gz)
  --out     STR    Output FASTA file (.fasta)
```

`--in` full path to the input FASTQ file (gzip-compressed).

`--out` full path to the output FASTA file.


## blast_PE_reads.py

```
Usage:
  python3 blast_PE_reads.py --ref <FILE> --gff <FILE> --reads1 <FILE> --reads2 <FILE> --out <FOLDER>

Mandatory:
  --ref        STR    Reference genome sequence file (FASTA)
  --gff        STR    Reference annotation file (GFF)
  --reads1     STR    Reads input file 1 (FASTA)
  --reads2     STR    Reads input file 2 (FASTA)
  --out        STR    Output folder

Optional:
  --cpus       INT    Number of threads for BLAST [10]
  --wordsize   INT    Word size in BLASTn [10]
  --minsim     INT    Minimal BLAST hit similarity [80]%
  --minlen     INT    Minimal BLAST hit length [25]
  --maxeval    FLOAT  Maximal BLAST hit evalue [0.001]
  --minscore   INT    Minimal BLAST score [30]
```

`--ref` full path to reference genome sequence file (FASTA).

`--gff` full path to reference annotation file (GFF).

`--reads1` full path to FASTA file containing R1 reads.

`--reads2` full path to FASTA file containing R2 reads.

`--out` full path to output folder. If the folder does not exist, it will be created.

`--cpus` specifies the number of threads to use in the BLASTn analysis. Increasing this value can substantially reduce the run time. Default: 10.

`--wordsize` specifies the word_size used by BLASTn. Increasing this value substantially reduces the computational costs, but also the sensitivity. Default: 10.

`--minsim` specifies the minimal similarity of BLAST hits (in percent) to be considered in the filtering step. Default: 80.

`--minlen` specifies the minimal BLAST hit length to be considered in the filtering step. Default: 25.

`--maxeval` specifies the maximal e-value acceptablein the BLAST hit filtering. Default: 0.001.

`--minscore` specifies the minimal BLAST hit score to be considered during the filtering. Default: 30.


## References


