# Coverage

## Subsampling
### G1 (Bryoria tortuosa)
Filtered reads and trimmed so all are 125 bp
```
metawrap read_qc -1 G1_1.fastq -2 G1_2.fastq -o G1_filtered  -t 14
trimmomatic PE G1_filtered_1.fastq.gz G1_filtered_2.fastq.gz G1_trimmomatic_1_PE.fastq G1_trimmomatic_1_SE.fastq G1_trimmomatic_2_PE.fastq G1_trimmomatic_2_SE.fastq HEADCROP:18 CROP:125 MINLEN:125
```
Subsampled, assembles, and anazyed: see Snakefile

### X12 (Alectoria sarmentosa)

Only filtered (all reads are 125 bp already)
```
nice -n 10 metawrap read_qc -1 14313X12_170811_D00294_0331_ACBBVEANXX_7_1.fastq -2 14313X12_170811_D00294_0331_ACBBVEANXX_7_2.fastq -o filtered_reads  -t 14
```

Subsampled, assembles, and anazyed: see Snakefile

### GT57 (Alectoria slurry)
Filtered reads and trimmed so all are 125 bp
```
metawrap read_qc -1 GT57_1.fastq -2 GT57_2.fastq -o GT57_filtered  -t 14
trimmomatic PE GT57_filtered/final_pure_reads_1.fastq  GT57_filtered/final_pure_reads_2.fastq GT57_trimmomatic_1_PE.fastq GT57_trimmomatic_1_SE.fastq GT57_trimmomatic_2_PE.fastq GT57_trimmomatic_2_SE.fastq HEADCROP:10 CROP:125 MINLEN:125
```
Subsampled, assembles, and anazyed: see Snakefile

## Visualizing G1
See `code/G1_visualization.R`

## Reanalysis of lichen metagenomes
### 1. Setting up the pipeline
* Donwload and install [blast_getLCA](https://github.com/frederikseersholm/blast_getLCA)
* Download and install the NCBI nt database:
```
module load StdEnv/2020
module load  blast+/2.11.0
update_blastdb.pl nt --decompress <- won't be using, but will keep for now jic

/cvmfs/ref.mugqic/genomes/blast_db/nt <- will use


```
* Config sra-toolkit
```
module load StdEnv/2020  gcc/9.3.0
module load sra-toolkit/2.10.8
echo '/repository/user/main/public/root = "/scratch/gultagr/sratoolkit-cache"' > $HOME/.ncbi/user-settings.mkfg
```
    
* setting up a virtual env for snakemake
```
 module load python/3.7
 virtualenv --no-download ~/env_snakemake
 source ~/env_snakemake/bin/activate
 pip install --no-index snakemake
  pip install --no-index requests
  pip install --no-index jsonschema
  pip install --no-index nbformat
  pip install --no-index pyparsing
```
* config snakemake to work with slurm (followed https://www.sichong.site/2020/02/25/snakemake-and-slurm-how-to-manage-workflow-with-resource-constraint-on-hpc/)
```
mkdir -p ~/.config/snakemake/slurm
mkdir logs_slurm
nano  ~/.config/snakemake/slurm/config.yaml

jobs: 30
cluster: "sbatch --account=def-tspribi -t {resources.time_min} --mem={resources.mem_mb} -c {resources.cpus} -o logs_slurm/{rule}_{wildcards} -e logs_slurm/{rule}_{wildcards}"

```

* megahit
```
module load megahit/1.2.9
```
* bbtools
```
module load bbmap/38.86
```
* seqtk
```
module load StdEnv/2020
module load seqtk/1.3
```

* kraken
```
module load StdEnv/2020
module load gcc/9.3.0
module load kraken2/2.1.1
```
build the db:
```
diamond --outfmt '6 qseqid sacc sseqid pident qlen length mismatch gapopen gaps evalue bitscore nident' --max-target-seqs 100  --db /cvmfs/ref.mugqic/genomes/blast_db/nt --threads 48  --query ../bbduk/cypho_SRR14722311.fasta --out tmp2
```
* CCMetgen
set up
```
git clone https://github.com/vrmarcelino/CCMetagen
module load  intel/2020.1.217
module load  kma/1.3.0
```
build the db:
```
 cp -r /cvmfs/ref.mugqic/genomes/blast_db/2020-04-22 ncbi_nt_2020-04-22
 wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
 gunzip nucl_gb.accession2taxid.gz
 cut -f 2-3 nucl_gb.accession2taxid > accession_taxid_nucl.map
 
cd /scratch/gultagr/ncbi_nt_2021_07_21
 wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz
 gzip -d nt.gz
 mv nt nt.fa
 wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
 gunzip nucl_gb.accession2taxid.gz
 cut -f 2-3 nucl_gb.accession2taxid > accession_taxid_nucl.map
 awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' < nt.fa > nt_sequential.fa
 #edit the rename_nt.py: change the input file from "nt.fa" to "nt_sequential.fa"
 python ../CCMetagen/benchmarking/rename_nt/rename_nt.py
 #had t remove the first line because it was empty
  sed '1d' nt_w_taxid.fas > nt_w_taxid2.fas
  kma index -i nt_w_taxid2.fas -o kma_db
 
```

* metaxa2

Installation
```
module load StdEnv/2020
module load hmmer/3.2.1
module load gcc/9.3.0
module load  blast+/2.11.0
module load mafft/7.471
module load vsearch/2.15.2

wget https://microbiology.se/sw/Metaxa2_2.2.3.tar.gz
tar -xvf Metaxa2_2.2.3.tar.gz
cd Metaxa2_2.2.3
```
Testing
```
cp /scratch/gultagr/coverage/data/fastq/SRR13685127.sra_*  .
metaxa2 -1 SRR13685127.sra_1.fastq.gz -2 SRR13685127.sra_2.fastq.gz -f p --mode m --plus -o SRR13685127
metaxa2 -1 SRR5808930.sra_1.fastq.gz -2 SRR5808930.sra_2.fastq.gz -f p --mode m --plus -o SRR5808930
```


### 1a. Designing part of pipeline for measuring lecanoro coverage
#### Plan A: de-novo binning plus taxonomy assignment

* EukCC
```
mkdir eukccdb
cd eukccdb
wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.tar.gz
tar -xzvf eukcc2_db_ver_1.tar.gz
export EUKCC2_DB=~/bin/eukccdb/eukcc2_db_ver_1

salloc --time=4:0:0 --ntasks=10 --mem-per-cpu 50G --account=def-tspribi
singularity pull docker://openpaul/eukcc2
singularity exec ~/bin/eukcc2_latest.sif eukcc single -h
singularity exec ~/bin/eukcc2_latest.sif eukcc single test.fa
singularity exec  ~/bin/eukcc2_latest.sif eukcc folder fasta_eukcc/ -o eukcc_out



source ~/env_snakemake/bin/activate
pip install eukcc
pip install guppy3
module load hmmer/3.2.1
module load gcc/9.3.0
module load openmpi/4.0.3
module load  metaeuk/4-a0f584d


mkdir eukcc_out
eukcc folder fasta_eukcc/ -o eukcc_out
```
INstalling epa-ng
```

git clone https://github.com/Pbdas/epa-ng
make
export PATH=$PATH:~/bin/epa-ng/bin/
```

* concoct
```
source ~/env_snakemake/bin/activate
pip install concoct
pip install pandas==0.25.1 #Successfully uninstalled pandas-1.3.0+computecanada
module load bwa/0.7.17
module load samtools/1.12

```
* testing the approach
```
cd testing
bwa index -a bwtsw ../assemblies/SRR14722310/final.contigs.fa
bwa mem ../assemblies/SRR14722310/final.contigs.fa ../../../data/fastq/SRR14722310.sra_1.fastq.gz ../../../data/fastq/SRR14722310.sra_2.fastq.gz > SRR14722310.sam
samtools sort SRR14722310.sam -o SRR14722310.bam
samtools index SRR14722310.bam

cut_up_fasta.py ../assemblies/SRR14722310/final.contigs.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa
concoct_coverage_table.py contigs_10K.bed  SRR14722310.bam > coverage_table.tsv
concoct --composition_file contigs_10K.fa --coverage_file coverage_table.tsv -b concoct_output/
merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv
mkdir concoct_output/fasta_bins
extract_fasta_bins.py ../assemblies/SRR14722310/final.contigs.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins


```


#### Plan B: Quast read alignment

Testing the approach. Will the reads aling to a MAG of a different lecanoro?
```
cd testing
module load StdEnv/2020
module load gcc/9.3.0
 module load  quast/5.0.2
 metaquast.py ../assemblies/SRR14722310/final.contigs.fa -o quast -r lecanoro_mag.fa --fungus -1 ../../../data/fastq/SRR14722310.sra_1.fastq.gz -2 ../../../data/fastq/SRR14722310.sra_2.fastq.gz --no-gc --no-plots --no-html --no-icarus --no-snps --no-read-stats -t 3
```
Only 7% aligned with avg. coverage of 1. Mostl likely this doesn't work

#### Plan C: supervised binning
```
cd ~/bin
git clone https://github.com/DRL/blobtools.git
cd blobtools
python setup.py install --user
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
./blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp
module load diamond/2.0.9
cd ~/scratch/diamonddb/
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
diamond makedb --in nr.gz -d nr
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
gunzip prot.accession2taxid.gz
#realized that need to build diamond db with taxonomy
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xfz taxdump.tar.gz names.dmp nodes.dmp
diamond makedb --in nr.gz -d nr_tax --taxonmap prot.accession2taxid --taxonnodes nodes.dmp --threads 20

cd ~/scratch/coverage/analysis/03_metagenome_reanalysis/testing/
diamond blastx -d ~/scratch/diamonddb/nr.dmnd -q ../assemblies/SRR14722310/final.contigs.fa -o diamond.tsv --outfmt 6  --sensitive --max-target-seqs 1 --evalue 1e-25 -b4 -c2 --memory-limit 187 --threads 48
~/bin/blobtools/blobtools taxify -f diamond.tsv  -m ~/scratch/diamonddb/prot.accession2taxid  -s 0 -t 2
~/bin/blobtools/blobtools create -i ../assemblies/SRR14722310/final.contigs.fa -b SRR14722310.bam -t diamond.tsv  -o SRR14722310 <- fails
~/bin/blobtools/blobtools create -i ../assemblies/SRR14722310/final.contigs.fa -b SRR14722310.bam -t SRR14722310.diamond.tsv.taxified.out -o SRR14722310
~/bin/blobtools/blobtools view -i SRR14722310.blobDB.json -o SRR14722310 --cov


diamond blastx -d ~/scratch/diamonddb/nr_tax.dmnd -q ../assemblies/SRR14722310/final.contigs.fa -o diamond_tax.tsv --outfmt 6 qseqid staxids bitscore  --sensitive --max-target-seqs 1 --evalue 1e-25 -b4 -c2 --memory-limit 187 --threads 48
~/bin/blobtools/blobtools create -i ../assemblies/SRR14722310/final.contigs.fa -b SRR14722310.bam -t diamond_tax.tsv  -o SRR14722310
~/bin/blobtools/blobtools view -i SRR14722310.blobDB.json -o SRR14722310 --cov  && \
~/bin/blobtools/blobtools plot -i SRR14722310.blobDB.json -o SRR14722310

diamond blastx -d ~/scratch/diamonddb/nr_tax.dmnd -q ../assemblies/SRR14722310/final.contigs.fa -o diamond_tax.tsv --outfmt  6  --sensitive --max-target-seqs 1 --evalue 1e-25 -b4 -c2 --memory-limit 187 --threads 48

diamond makedb --in nr.gz -d nr_tax2 --taxonmap names.dmp --taxonnodes nodes.dmp --threads 20
```

Lost track, will start from a clean sleet here:
```
 mkdir new
 cd new/
 ~/bin/blobtools/blobtools create -i ../../assemblies/SRR14722310/final.contigs.fa -b ../SRR14722310.bam -t ../diamond_tax.tsv  -o SRR14722310
 Traceback (most recent call last):
 File "/home/gultagr/bin/blobtools/blobtools", line 7, in <module>
 main()
 File "/home/gultagr/bin/blobtools/lib/interface.py", line 60, in main
 create.main()
 File "/home/gultagr/bin/blobtools/lib/create.py", line 108, in main
 blobDb.parseHits(hit_libs)
 File "/home/gultagr/bin/blobtools/lib/BtCore.py", line 420, in parseHits
 for hitDict in BtIO.readTax(hitLib.f, set(self.dict_of_blobs)):
 File "/home/gultagr/bin/blobtools/lib/BtIO.py", line 529, in readTax
 'score' : float(col[2])
 IndexError: list index out of range
```
 Maybe the problem with missing staxids? WIll try again to produce full outfmt 6 diamond output:
```
diamond blastx -d ~/scratch/diamonddb/nr_tax.dmnd -q ../../assemblies/SRR14722310/final.contigs.fa -o diamond_tax2.tsv --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore  --sensitive --max-target-seqs 1 --evalue 1e-25 -b2 -c1 --memory-limit 187 --threads 48
head diamond_tax2.tsv
>k141_33432              111     k141_33432      MBV8116294.1    76.7    86      11      2       1       252     13      91      2.97e-29        111
>k141_43183              118     k141_43183      MBI3884212.1    56.8    111     47      1       330     1       316     426     2.25e-28        118
>k141_25074      2316528 176     k141_25074      WP_129227759.1  97.8    92      2       0       1       276     49      140     2.33e-54        176
>k141_27860      2483403 128     k141_27860      WP_172349139.1  91.9    74      6       0       361     140     772     845     1.90e-31        128
```
same problem, not all contig have taxid even when they have sseqid. maybe something is wrong with db configuration.

will test if the prolem indeed in missing values in diamond_tax.tsv
```
grep -Ev $'^\t|\t\t|\t$' ../diamond_tax.tsv > diamond_tax_remove_missing.tsv
 ~/bin/blobtools/blobtools create -i ../../assemblies/SRR14722310/final.contigs.fa -b ../SRR14722310.bam -t diamond_tax_remove_missing.tsv  -o SRR14722310
 
 ~/bin/blobtools/blobtools view -i SRR14722310.blobDB.json -o SRR14722310 --cov  && \
 ~/bin/blobtools/blobtools plot -i SRR14722310.blobDB.json -o SRR14722310
```
Worked! Exploring output:
```
$ head  SRR14722310.SRR14722310.blobDB.json.bestsum.phylum.p8.span.100.blobplot.stats.txt
## 1.1.1
## bam0=/scratch/gultagr/coverage/analysis/03_metagenome_reanalysis/testing/SRR14722310.bam
# name  colour  count_visible   count_visible_perc      span_visible    span_visible_perc       n50     gc_mean gc_std  bam0_mean       bam0_std        bam0_read_map   bam0_read_map_p
all     None    65,934  100.0%  72,226,316      100.0%  2,481   0.61    0.093   6.1     12.1    5,678,526       51.1%
Ascomycota      #1f77b4 4,613   100.0%  33,473,493      100.0%  13,054  0.49    0.045   14.8    18.5    4,006,758       36.1%
no-hit  #d3d3d3 34,275  100.0%  22,058,335      100.0%  653     0.61    0.099   5.5     14.1    986,537 8.9%
Proteobacteria  #ff7f0e 13,886  100.0%  8,246,684       100.0%  604     0.67    0.045   5.8     5.4     347,411 3.1%
Acidobacteria   #d62728 7,024   100.0%  3,911,368       100.0%  566     0.62    0.038   4.1     4.4     116,253 1.0%
Verrucomicrobia #9467bd 2,077   100.0%  1,911,581       100.0%  1,184   0.65    0.038   5.5     3.1     85,025  0.8%
Actinobacteria  #e377c2 1,818   100.0%  864,771 100.0%  482     0.69    0.036   5.2     7.0     32,256  0.3%
```
Ascomycota mean coverage makes sense, but how doest it look relative to all coverages of ascomycota contigs?
```
grep "Ascomycota"  SRR14722310.SRR14722310.blobDB.table.txt
grep "Ascomycota"  SRR14722310.SRR14722310.blobDB.table.txt | awk '{print $5}' | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'
```
Median coverage according to this oneliner is 15.5.

Designing report-creating system

```
echo "SRR14722310" > report1
grep "Ascomycota"  SRR14722310.SRR14722310.blobDB.table.txt | awk '{print $5}' | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' > report2
paste report1 report2 > report3
```

Double check with correct paths:
```
 ~/bin/blobtools/blobtools create -i assemblies/SRR14722310/final.contigs.fa -b blobtools/SRR14722310.bam -t blobtools/SRR14722310_diamond_final.tsv -o blobtools/SRR14722310
  ~/bin/blobtools/blobtools view -i blobtools/SRR14722310.blobDB.json -o blobtools/ --cov
  ~/bin/blobtools/blobtools plot -i blobtools/SRR14722310.blobDB.json -o blobtools/
```

### 2. Prepare dataset:
* Got SRA ids from Lendemer:
    * Downloaded SRA run info for PRJNA700635 and PRJNA731936 > `analysis/03_metagenome_reanalysis/SraRunInfo_PRJNA731936.csv analysis/03_metagenome_reanalysis/SraRunInfo_PRJNA700635.csv`
    * Downloaded [Appendix S2](https://bsapubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fajb2.1339&file=ajb21339-sup-0002-AppendixS2.xlsx) from Lendemer et al. 2019 with the voucehr metadata for their metagenomes
    * Made finel list of SRA IDs using `code/getting_SRA_id_for_lendemer_data.R`. Only kept IDs that matched between the metadata from the NCBI and from the Appendix > `analysis/03_metagenome_reanalysis/sra_ids_lendemer.txt`

* Got SRA ids from other sources
    * See `analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.csv`. Copied all non-lendemer SRA ids into `analysis/03_metagenome_reanalysis/sra_ids_other.txt`
    * Combined the two files into the final list of SRA IDs to use:
```
cat analysis/03_metagenome_reanalysis/sra_ids_lendemer.txt analysis/03_metagenome_reanalysis/sra_ids_other.txt >analysis/03_metagenome_reanalysis/sra_ids_all.txt
```
* Got a list of libraries from our lab that are not likely to be used on anything else, the table is [here](https://docs.google.com/spreadsheets/d/10WtB3Vm5ZzAPAiXhnSwq6kx6AbhfV6w4smsEvDLOfos/edit#gid=0)

Downloaded the data using sratoolkit
```
cd data/fastq
prefetch --option-file analysis/03_metagenome_reanalysis/sra_ids_all.txt

for file in /scratch/gultagr/sratoolkit-cache/sra/*sra

do
fasterq-dump --split-files "$file"
done

gzip *.fastq


```



### 3. RUn the snakefile
```
snakemake --profile slurm -f --max-jobs-per-second 3 --max-status-checks-per-second 3


snakemake --cores 10 -f
```

Tried to run some diamond on debary
```
diamond blastx -d /data/databases/diamond/nr.dmnd -q ../assemblies/T1867/final.contigs.fa -o T1867_diamond_tax.tsv --outfmt 6 qseqid staxids bitscore  --sensitive --max-target-seqs 1 --evalue 1e-25 -b4 -c2 --memory-limit 150 --threads 24
```




### 4. Add rule to count bp for the metagenomes (issue_4)

Tried in `bbuk/tmp`:
```
reformat.sh in=../../../../data/lendemer_data/SRR13685123_1.fastq in2=../../../../data/lendemer_data/SRR13685123_2.fastq >out_bp 2>&1
```
Modified   `bbduk/Snakefile` to include rules producing table `bp_report.txt`



### 5. Added more sequences as search queries, to include a wider diversity of fungi and loci (ITS, LSU, SSU) (issue_8).\
Found more hits

### 6. Added more metagenomes to reanalyze (issue_12)\
First, added the rest of Lendemer data as it became public. Updated script `getting_SRA_id_for_lendemer_data.R` and generated with it:
* table `lendemer_table_sra.txt` contains all lendemer data from the first and second batch`getting_SRA_id_for_lendemer_data
* list `sra_ids_lendemer_first_batch.txt` with all lendemer data reanalyzed in the first batch
* list `sra_ids_lendemer_second_batch.txt` with second batch of Lendemer data`sra_ids_first_batch

Also, added table `all_metagenome_reanalysis.csv` that includes all data for reanalysis. This includes not only Lendemer data (both batches), but other metagenomes: from Philipp's paper and already published from other papers (ATP, Lutzoni paper, Leavitt Symbiosis paper, skimming, etc)

Moved libraries produced by our lab from debary manually, used a bach script to rename them to match the format of SRA-derived files:
```
bash code/rename_our_fastq.sh
```

### 7. Added more metagenomes to reanalyze yet
Added third batch including: 
* SRA files from Leavitt and Lutzoni labs, which I accidentally left out the last time. They are in the `ra_ids_addition_2021.10.15.txt`
```
cd data/fastq
prefetch --option-file ../../analysis/03_metagenome_reanalysis/sra_ids_addition_2021.10.15.txt

for file in /scratch/gultagr/sratoolkit-cache/sra/*sra

do
fasterq-dump --split-files "$file"
done

gzip *.fastq
```
* Files from debary. The list is [here](https://docs.google.com/spreadsheets/d/10WtB3Vm5ZzAPAiXhnSwq6kx6AbhfV6w4smsEvDLOfos/edit#gid=0)

All files are added into the `all_metagenome_reanalysis.txt`


### Preparing dat ato send to Paul
1. Made temp folder `data/fatsq/our_data` and copied there all fastq files from our lab (~350 Gb)
2. Made temp folder `analysis/03_metagenome_reanalysis/assemblies_send_to_paul`. Deleted all extra files (.bwt files etc) and only left assemblies and log files (~50 Gb)

## Visualizing X12 and GT57

compile reports
```
mkdir analysis/01_subsampling/X12/reports
cat analysis/01_subsampling/X12/blast/tmp_report_* >  analysis/01_subsampling/X12/reports/blast_report.txt
cat analysis/01_subsampling/X12/quast/tmp_report_* > analysis/01_subsampling/X12/reports/quast_report.txt
 cat analysis/01_subsampling/X12/bbduk/tmp_report_* > analysis/01_subsampling/X12/reports/bbduk_report.txt
mkdir analysis/01_subsampling/GT57/reports
cat analysis/01_subsampling/GT57/quast/tmp_report_* > analysis/01_subsampling/GT57/reports/quast_report.txt
cat analysis/01_subsampling/GT57/bbduk/tmp_report_* > analysis/01_subsampling/GT57/reports/bbduk_report.txt
```

run `code/alectoria_visualization.R`

## Nonpareil
1. Sampled 500 Mbp of each metagenome
```
seqtk sample -s 123 analysis/01_subsampling/X12/X12_1.fastq.gz 4000000 > analysis/04_nonpareil/X12/X12_1_500Mbp.fastq
seqtk sample -s 123 analysis/01_subsampling/G1/G1_trimmomatic_1.fastq.gz 4000000 > analysis/04_nonpareil/G1/G1_1_500Mbp.fastq
seqtk sample -s 123 analysis/01_subsampling/GT57/GT57_trimmomatic_1_PE.fastq.gz 4000000 > analysis/04_nonpareil/GT57/GT57_1_500Mbp.fastq
```
2. Checked bp in each sample
```
reformat.sh in=analysis/04_nonpareil/X12/X12_1_500Mbp.fastq > analysis/04_nonpareil/X12/out_bp 2>&1
reformat.sh in=analysis/04_nonpareil/G1/G1_1_500Mbp.fastq > analysis/04_nonpareil/G1/out_bp 2>&1
reformat.sh in=analysis/04_nonpareil/GT57/GT57_1_500Mbp.fastq > analysis/04_nonpareil/GT57/out_bp 2>&1
```
All good.

3. Ran Nonpareil
```
nonpareil -s analysis/04_nonpareil/X12/X12_1_500Mbp.fastq -T kmer -f fastq -b analysis/04_nonpareil/X12/X12_1_500Mbp -t 10 -R 50000
nonpareil -s analysis/04_nonpareil/G1/G1_1_500Mbp.fastq -T kmer -f fastq -b analysis/04_nonpareil/G1/G1_1_500Mbp -t 10 -R 50000
nonpareil -s analysis/04_nonpareil/GT57/GT57_1_500Mbp.fastq -T kmer -f fastq -b analysis/04_nonpareil/GT57/GT57_1_500Mbp.fastq -t 10 -R 50000
```
4. Visualize curves with `code/nonpareil.R`


