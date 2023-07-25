#!/bin/bash

# run pre-fetch on your file with a list of sra_ids
# this list can be found in 'output/sra_accession_for_sra_download.txt'

prefetch --option-file "$1" --max-size 21g # added file size

# now run the fasterq-dump on these cached prefetched files

mkdir ./data/fastq # make an output directory, if one doesn't already exist

# adjust path to be reflective of where your cache is located
for file in /data/andrewc/scripts/sra_tools/sra_cache/sra/*sra

do
fasterq-dump --outdir ./data/fastq/ --split-files "$file"
done

# zip files for sourmash script
pigz ./data/fastq/*.fastq

# NOTE: CHECK THAT ALL FILES HAVE DOWNLOADED. SOME MAY BE >20GB
# FILES LARGER THAN 20GB NEED TO HAVE `--max-size` FLAG ADDED TO PREFETCH
# LOOK AT THE `prefetch --help` FOR MORE DETAILS.