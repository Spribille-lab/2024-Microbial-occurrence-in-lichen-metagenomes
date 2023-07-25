#!/bin/bash

# creates the initial sourmash 'signatures' used for genome comparisons. 
# input file is the list of downloaded SRA files
# input file is found at `output/sra_accession_for_sra_download.txt`

mkdir ./data/sourmash_signatures # make directory, if one doesn't already exist

while IFS= read -r LINE; do
  FILE_1=./data/fastq/"$LINE"_1.fastq.gz
  FILE_2=./data/fastq/"$LINE"_2.fastq.gz
  
  sourmash sketch dna -p k=21,k=31,k=51 "$FILE_1" "$FILE_2" --name "$LINE" -o ./data/sourmash_signatures/"$LINE".zip
  
done < "$1"

# run sourmash compare to check for duplicate metagenomes

sourmash compare ./data/sourmash_signatures/*.zip --ksize 21 --csv ./results/sourmash_compare_matrix.csv

Rscript ./scripts/sourmash_compare_converter.R