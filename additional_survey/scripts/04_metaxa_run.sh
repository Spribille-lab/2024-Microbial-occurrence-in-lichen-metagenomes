#!/bin/bash

# run metaxa over all files that match our input SRA numbers
# input file is found at `output/sra_accession_for_sra_download.txt`

mkdir ./results/metaxa_runs # make output directory, if one doesn't already exist
mkdir ./results/metaxa_bacteria # make output directory, if one doesn't already exist

cd ./results/metaxa_runs/

while IFS= read -r LINE; do
  FILE_1=../../data/fastq/"$LINE"_1.fastq.gz
  FILE_2=../../data/fastq/"$LINE"_2.fastq.gz
  
  mkdir "$LINE"
  cd ./"$LINE"
  
  metaxa2 -1 ../"$FILE_1" -2 ../"$FILE_2" -f p --mode m --plus --graphical F --split_pairs F -o "$LINE" -cpu 32
  
  cp "$LINE".bacteria.fast ../../metaxa_bacteria # copy the bacterial .fasta file
  
  cd ..
  
done < "$1"

cd ../.. # get back to project home directory


