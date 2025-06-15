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
  --ref      STR    Reference genome sequence file (FASTA)
  --gff      STR    Reference annotation file (GFF)
  --reads1   STR    Reads input file 1 (FASTA)
  --reads2   STR    Reads input file 2 (FASTA)
  --out      STR    Output folder
```

`--ref` full path to reference genome sequence file (FASTA).

`--gff` full path to reference annotation file (GFF).

`--reads1` full path to FASTA file containing R1 reads.

`--reads2` full path to FASTA file containing R2 reads.

`--out` full path to output folder. If the folder does not exist, it will be created.


## References


