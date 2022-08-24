#!/bin/bash
#SBATCH --account=def-tspribi
#SBATCH --time=1-0:0
#SBATCH --cpus-per-task=64
#SBATCH --job-name=Script
#SBATCH --output=Spript.logs.out
#SBATCH --mem=249G

iqtree -s gtdbtk.bac120.user_msa.fasta --seqtype AA -bb 50000 -pre rhizobiales -seed 12345 -m TEST -T 64

