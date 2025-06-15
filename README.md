# blast2counts
analyze RNA-seq data with BLASTn

## analyze_susp_fastqs.py
This script allows the extraction of valid reads from a given FASTQ file. Reads with long N stretches are ignored.


```
Usage:
  python3 analyze_susp_fastqs.py --in <FASTQ_GZ_FILE> --out <FASTA_FILE>
--in <FASTQ.GZ_FILE>
					--out <FASTA_FILE>
  --in      STR    Input FASTQ file (.fastq.gz)
  --out     STR    Output FASTA file (.fasta)

```

`--in` full path to the input FASTQ file (gzip-compressed).

`--out` full path to the output FASTA file.


## 

## References


