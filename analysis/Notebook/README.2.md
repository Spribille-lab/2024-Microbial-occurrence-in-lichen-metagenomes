# Coverage 2.0: working with Paul

## Data

### 1. Prepared dataset:
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

### 2. Downloaded the data using sratoolkit
```
cd data/fastq
prefetch --option-file analysis/03_metagenome_reanalysis/sra_ids_all.txt

for file in /scratch/gultagr/sratoolkit-cache/sra/*sra

do
fasterq-dump --split-files "$file"
done

gzip *.fastq


```
### 3. Filtered dataset

* Using sourmash detected that some SRA libraries are identical despite having different IDs (Paul). All duplicated are shared between two Lendemer studies, PRJNA731936 and PRJNA700635. The list is in `analysis/03_metagenome_reanalysis/similar_datasets.csv`, in total there are 42 identical pairs.
* Removed the duplicated libraries from the smaller PRJNA700635. The other copies remain in the analysis
```
for ID_SAMPLE in `cut -d ',' -f1 analysis/03_metagenome_reanalysis/similar_datasets.csv`; do rm data/fastq/"$ID_SAMPLE"*; done
```

## Analysis
### 1. Incorporated Paul's assemblies and reran parts of analysis dependent on them (yeast detection in the assemblies + assembly size)
1. Moved them to `analysis/03_metagenome_reanalysis/assemblies_paul`

2. Modified Snakefile:
* Changed input for the blast rule
* Changed how assembly lengths are calculated: now used stats.sh from bbtools
* Removed part of script calculating lecanoro coverage 
* Didn't touch the other parts of the pipeline

3. Reran snakemake starting at the rule blast and every rule dependent on it:
```
snakemake --cores 1 -f -R blast
```
4. Added analysis of presence/absence of Trebouxia
* Modified Snakefile to replace balst rule with metaxa on assembly
* Modified the former `code/yeast_presence_metagenome_reanalysis.R` -> `code/symbiont_presence_metagenome_reanalysis.R`. Now it makes the dothist plot too.
* Renamed `code/metagenome_reanalysis_viz.R` -> `code/exploring_yeast_pres_abs.R`. This script contains my exploration of different slices of the dataset and yeast presence/absence in them. Might use it for a diff paper (in lichenologist?), where I will focus on the yeasts.

5. Calculating the most common lineages in lichens
* Used metaxa functionality to summarize metaxa output from all reads and all assemblies. Added rules metaxa_tax_profile_assembly and metaxa_tax_profile_reads to the Snakefile that run metaxa2_ttt on each `*.taxonomy.txt`
* Used metaxa2_dc to assemble occurrence matrix:
```
metaxa2_dc *.level_3.txt -o metaxa_level_3_combined.txt
metaxa2_dc *.level_4.txt -o metaxa_level_4_combined.txt
metaxa2_dc *.level_5.txt -o metaxa_level_5_combined.txt
metaxa2_dc *.level_6.txt -o metaxa_level_6_combined.txt
```
6. Added `code/metaxa.R` to analyze the metaxa results
    * Saved `analysis/03_metagenome_reanalysis/occurrence_assembly_metaxa_level5.tsv` and `analysis/03_metagenome_reanalysis/occurrence_reads_metaxa_level5.tsv` showing the lineage with the highest occurrence counts
    * Made figure `results/figures/metaxa.png` showing number of metaxa-detected lineages depending on depth and number of MAGs
7. Prepped data for David to make the trees
    * Used `code/make_list_fungal_algal_mags.R` to make files with lists of fungal and algal mags. Saved them as `analysis/05_MAGs/tables/list_fungal_MAGs.txt` and `analysis/05_MAGs/tables/list_algal_MAGs.txt`
    * Copy them into a folder on cedar
```
xargs -a ../../tables/list_algal_MAGs.txt cp -t ~/projects/def-tspribi/mags_for_phylogenimcs_jan22/algal
xargs -a ../../tables/list_fungal_MAGs.txt cp -t ~/projects/def-tspribi/mags_for_phylogenimcs_jan22/fungal
```
    * Manually copied reference genomes from `data/share/ref_genomes/`
    
    
 8. Phylogenomic trees for the eukariotes
 * setting up
 ```
 module load singularity/3.8
source ~/env_snakemake/bin/activate
git clone --recursive https://github.com/reslp/phylociraptor.git
```

9. Metaxa-style analysis with curated database of 16S and 18S
* The idea is to keep metaxa analysis for euks, and for bacteria use metaxa-extracted rDNA fastas but reclassify them according to GTDB taxonomy using IdTaxa
* Compiled database
    * downloaded [SBDI Sativa curated 16S GTDB database](https://scilifelab.figshare.com/articles/dataset/SBDI_Sativa_curated_16S_GTDB_database/14869077), Version 3, posted on 10.11.2021 accessed on 2022-02-01
    * decompressed
    ```
    gzip -d data/db/gtdb-sbdi-sativa.r06rs202.assignTaxonomy.fna.gz
    gzip -d data/db/gtdb-sbdi-sativa.r06rs202.fna.gz
    ```
    * compiled with  `code/idtaxa_compiling_db.R`.
    * Used `analysis/03_metagenome_reanalysis/assembly_SRR14722232.extraction.fasta` as a test dataset. There is a problem with acetos: can only get genus-level annotations if the confidence threshold is at 40%!!
    * As of now, this is nt the full db (will need to add 16S of the undescribed lineages).

10. testing ways to obtain 16S of MAGs from the data set using MarkerMAG. Used SRR13125477  for now, since it has relatively few MAGs + a MAG from the undescribed Aceto clade sister to the  BOG908 and CAJCIS01 clade + another Aceto MAG - to check how that affects the outcome
    * Made a new analysis folder `analysis/06_16S_mag_matching`
    * Testing out on debary before trying to move to computecanada.
```
pip3 install MarkerMAG
PATH=$PATH:/data/tagirdzh/miniconda3/envs/metawrap-env/bin/
PATH=$PATH:/data/tagirdzh/bin/
cd analysis/06_16S_mag_matching
mkdir selected_mags/SRR13125477
cp ../05_MAGs/MAGs/bacs/public_SRR13125477_metawrap_bin.2.fa.gz selected_mags/SRR13125477/
cp ../05_MAGs/MAGs/bacs/public_SRR13125477_metawrap_bin.5.fa.gz selected_mags/SRR13125477/
gzip -d selected_mags/SRR13125477/*.fa.gz
seqtk seq -a ../../data/fastq/SRR13125477.sra_1.fastq.gz > SRR13125477.sra_1.fa
seqtk seq -a ../../data/fastq/SRR13125477.sra_2.fastq.gz > SRR13125477.sra_2.fa
MarkerMAG rename_reads -r1 SRR13125477.sra_1.fa -r2 SRR13125477.sra_2.fa -p SRR13125477 -t 5
```
* First, tried running with metaxa-obtained 16S sequences. Failed.
```
MarkerMAG link -p tmp_SRR13125477 -marker ../03_metagenome_reanalysis/assembly_SRR13125477.bacteria.fasta -mag selected_mags/SRR13125477/ -x fa -r1 SRR13125477_R1.fasta -r2 SRR13125477_R2.fasta -t 12 -o tmp_SRR13125477 -no_polish
[2022-02-02 06:17:12] parameters for linking
+ mismatch:    2%
+ min_M_len:   45bp
+ min_M_pct:   35%
+ min_link_num_gnm:    9
+ min_link_num_ctg:    3
+ rd2_end_seq_len:     1000bp
+ max_short_cigar_pct: 75,85
[2022-02-02 06:17:12] parameters for estimating copy number
+ MAG_cov_subsample_pct:       25%
+ min_insert_size_16s: -1000bp
+ ignore_ends_len_16s: 150bp
+ ignore_lowest_pct:   25%
+ ignore_highest_pct:  25%
+ both_pair_mapped:    False
[2022-02-02 06:17:15] Rd1: quality control provided 16S rRNA gene sequences to:
[2022-02-02 06:17:15] Rd1: remove sequences shorter than 1200 bp
[2022-02-02 06:17:15] Rd1: cluster at 99% identity and keep only the longest one in each cluster
[2022-02-02 06:17:15] Rd1: qualified 16S rRNA gene sequences exported to:
[2022-02-02 06:17:15] assembly_SRR13125477.bacteria_unpolished_min1200bp_c99.fasta
[2022-02-02 06:17:16] Rd1: mapping input reads to marker genes
[2022-02-02 06:17:16] Rd1: sorting mappping results
Traceback (most recent call last):
File "/data/tagirdzh/miniconda3/bin/MarkerMAG", line 133, in <module>
link_16s.link_16s(args, config_dict)
File "/data/tagirdzh/miniconda3/lib/python3.8/site-packages/MarkerMAG/link_16s.py", line 2765, in link_16s
os.remove(input_reads_to_16s_sam_bowtie)
FileNotFoundError: [Errno 2] No such file or directory: 'tmp_SRR13125477/tmp_SRR13125477_rd1_wd/tmp_SRR13125477_input_reads_to_16S.sam'
[E::hts_open_format] Failed to open file tmp_SRR13125477/tmp_SRR13125477_rd1_wd/tmp_SRR13125477_input_reads_to_16S.sam
samtools sort: can't open "tmp_SRR13125477/tmp_SRR13125477_rd1_wd/tmp_SRR13125477_input_reads_to_16S.sam": No such file or directory
```
* Second, tried identifying 16S in the MAGs - didn't find anything niether with MarkerMAG nor with blast
```
MarkerMAG barrnap_16s -p Test -g selected_mags/SRR13125477/ -x fa -t 6 -force

 blastn -query aceto_16s_sample.fa -subject ../../selected_mags/SRR13125477/public_SRR13125477_metawrap_bin.2.fa -outfmt 6 -evalue 1e-5
```


Setting up:
```
conda install -c bioconda matam
index_default_ssu_rrna_db.py
```
Prepping data:
```
cp ../../data/fastq/SRR13125477.sra_*.fastq.gz .
gzip -d SRR13125477.sra_*.fastq.gz

MarkerMAG matam_16s -p soil -r1 SRR13125477_R1.fasta -r2 SRR13125477_R2.fasta -pct 1,5,10,25,50,75,100 -i 0.999 -d /data/tagirdzh/miniconda3/opt/matam-1.6.0/db/ # failed because thier read renaming module doesn't work

#renamed manually
awk '{ if ('NR%4==1' || 'NR%4==3'){ $1=$1".1" } print }' SRR13125477.sra_1.fastq  > SRR13125477_renamed.sra_1.fastq
awk '{ if ('NR%4==1' || 'NR%4==3'){ $1=$1".2" } print }' SRR13125477.sra_2.fastq  > SRR13125477_renamed.sra_2.fastq
#make fastas with matching names
seqtk seq -a SRR13125477_renamed.sra_1.fastq > SRR13125477_renamed.sra_1.fasta
seqtk seq -a SRR13125477_renamed.sra_2.fastq > SRR13125477_renamed.sra_2.fasta

#tried again
MarkerMAG matam_16s -p soil -r1 SRR13125477_renamed.sra_1.fastq -r2 SRR13125477_renamed.sra_2.fastq -pct 1,5,10,25,50,75,100 -i 0.999 -d /data/tagirdzh/miniconda3/opt/matam-1.6.0/db/SILVA_128_SSURef_NR95

```

* Tried using the output for the link module
```
MarkerMAG link -p tmp_SRR13125477 -marker soil_Matam16S_wd/soil_assembled_16S_uclust_0.999.fasta -mag selected_mags/SRR13125477/ -x fa -r1 SRR13125477_R1.fasta -r2 SRR13125477_R2.fasta -t 12 -o tmp_SRR13125477 -no_polish
```

* Third, tried to assemble 16S sequences de novo using matam module of MarkerMAG. Failed with the same samtools error, which I raies as an issue on [github](https://github.com/songweizhi/MarkerMAG/issues/2); no response from the authors

* Trying throwing spagetti at the wall

```
pip3 install MarkerMAG==1.1.3

#make sam file manually
bowtie2-build soil_Matam16S_wd/soil_assembled_16S_uclust_0.999.fasta  16S.fai
bowtie2 -x 16S.fai -1 SRR13125477_renamed.sra_1.fastq -2 SRR13125477_renamed.sra_2.fastq -S manual_linking/SRR13125477_input_reads_to_16S.sam -p 8
samtools sort -n -O sam --threads 8 -o manual_linking/SRR13125477_input_reads_to_16S.sorted.sam manual_linking/SRR13125477_input_reads_to_16S.sam

MarkerMAG link -p tmp_SRR13125477 -marker soil_Matam16S_wd/soil_assembled_16S_uclust_0.999.fasta -mag selected_mags/SRR13125477/ -x fa -r1 SRR13125477_R1.fasta -r2 SRR13125477_R2.fasta -t 12 -o  manual_linking/tmp_SRR13125477 -no_polish -sorted_sam16s manual_linking/SRR13125477_input_reads_to_16S.sorted.sam
>Traceback (most recent call last):
File "/data/tagirdzh/miniconda3/bin/MarkerMAG", line 133, in <module>
link_16s.link_16s(args, config_dict)
File "/data/tagirdzh/miniconda3/lib/python3.8/site-packages/MarkerMAG/link_16s.py", line 2811, in link_16s
pool_parse_sam16s.map(parse_sam16s_worker, list_for_parse_sam16s_worker)
File "/data/tagirdzh/miniconda3/lib/python3.8/multiprocessing/pool.py", line 364, in map
return self._map_async(func, iterable, mapstar, chunksize).get()
File "/data/tagirdzh/miniconda3/lib/python3.8/multiprocessing/pool.py", line 771, in get
raise self._value
KeyError: 'soil_subsample_10_90'

pip3 install MarkerMAG==1.1.5 -> same
pip3 install MarkerMAG==1.0.41 -> same

conda create -n markermag python==3.7
 conda activate markermag
 pip3 install MarkerMAG
 PATH=$PATH:/data/tagirdzh/bin/
 conda install seqtk
 conda install samtools
 conda install bowtie2
conda install -c bioconda blast
conda install -c bioconda spades
conda install -c bioconda barrnap

MarkerMAG link -p tmp_SRR13125477 -marker soil_Matam16S_wd/soil_assembled_16S_uclust_0.999.fasta -mag selected_mags/SRR13125477/ -x fa -r1 SRR13125477_renamed.sra_1.fasta -r2 SRR13125477_renamed.sra_2.fasta -t 12 -o  manual_linking/tmp_SRR13125477 -sorted_sam16s manual_linking/SRR13125477_input_reads_to_16S.sorted.sam  -> same




 

```


12. Identifying 16S genes found in the MAGs by PROKKA
```
grep "rRNAs" analysis/07_annotate_MAGs/*/*.log
```

Found:
* private_T1888_metawrap_bin.6 (CAHJWL01): CAHJWL01T1888_02966  5S ribosomal RNA
* private_T1904_metawrap_bin.1 (CAIMSN01): CAIMSN01T1904_02346  16S ribosomal RNA
* private_T1916_metawrap_bin.2 (CAHJWL01): CAHJWL01T1916_00822  5S ribosomal RNA
* private_TS1974_metawrap_bin.6 (LMUY01): LMUY01TS1974_01766  5S ribosomal RNA
* private_VT16_metawrap_bin.4 (CAHJWL01): CAHJWL01VT16_00620    16S ribosomal RNA
CAHJWL01VT16_00623  23S ribosomal RNA
CAHJWL01VT16_00624  5S ribosomal RNA
* private_VT1_metawrap_bin.2 (CAHJWO01): CAHJWO01VT1_04001 16S ribosomal RNA (partial)
* private_VT22_metawrap_bin.8 (EB88): EB88VT22_03259   5S ribosomal RNA
EB88VT22_03260  23S ribosomal RNA (partial)
EB88VT22_03348  16S ribosomal RNA (partial)
* private_VT34_metawrap_bin.3 (CAHJXG01): CAHJXG01VT34_01425 5S ribosomal RNA
CAHJXG01VT34_01426 23S ribosomal RNA
* private_X16_metawrap_bin.17 (CAHJWL01): CAHJWL01X16_01058   5S ribosomal RNA
CAHJWL01X16_01331 16S ribosomal RNA (partial)
* public_ERR4179390_metawrap_bin.2: CAHJWL01ERR4179390_02634 16S ribosomal RNA (partial)
CAHJWL01ERR4179390_02996 5S ribosomal RNA
* public_SRR11456915_metawrap_bin.3: NOSTOCSRR11456915_07647  16S ribosomal RNA (partial)
* public_SRR11456915_metawrap_bin.9: CAHJWO01SRR11456915_01547  5S ribosomal RNA
* public_SRR11456918_metawrap_bin.7: SPHINGOMSRR11456918_03063  5S ribosomal RNA
* public_SRR14722095_metawrap_bin.2: TERRIGLOSRR14722095_03551   5S ribosomal RNA
* public_SRR14722108_metawrap_bin.1 (CAHJWO01): CHTONIOSRR14722108_02680 16S ribosomal RNA
* public_SRR14722117_metawrap_bin.2: CAHJWO01SRR14722117_00405 23S ribosomal RNA
CAHJWO01SRR14722117_00408  16S ribosomal RNA
* public_SRR14722125_metawrap_bin.3: TERRIGLOBUSSRR14722125_01218  5S ribosomal RNA
TERRIGLOBUSSRR14722125_01219  23S ribosomal RNA
TERRIGLOBUSSRR14722125_01222  16S ribosomal RNA
* public_SRR14722154_metawrap_bin.3: EB88SRR14722154_01884  16S ribosomal RNA (partial)
EB88SRR14722154_02563 23S ribosomal RNA
EB88SRR14722154_02564  5S ribosomal RNA
* public_SRR14722157_metawrap_bin.1: EB88SRR14722157_02068   5S ribosomal RNA
* public_SRR14722229_metawrap_bin.6: TERRIGLOSRR14722229_03588 5S ribosomal RNA
* public_SRR2387885_metawrap_bin.6: CAHJWL01SRR2387885_03273  23S ribosomal RNA (partial)
CAHJWL01SRR2387885_03274  5S ribosomal RNA
* public_SRR7232211_metawrap_bin.1: CAIMSN01SRR7232211_02935  16S ribosomal RNA
CAIMSN01SRR7232211_02939  23S ribosomal RNA
CAIMSN01SRR7232211_02940  5S ribosomal RNA

undescribed lineages:
* Acetobacteraceae gen. sp.: none
* CAHJWL01: complete
    * private_VT16_metawrap_bin.4 (CAHJWL01): CAHJWL01VT16_00620 16S ribosomal RNA
    * public_ERR4179390_metawrap_bin.2: CAHJWL01ERR4179390_02634 16S ribosomal RNA (partial)
* CAHJWO01: complete
    * public_SRR14722108_metawrap_bin.1 (CAHJWO01): CHTONIOSRR14722108_02680 16S ribosomal RNA
    * public_SRR14722117_metawrap_bin.2: CAHJWO01SRR14722117_00408  16S ribosomal RNA
    * private_VT1_metawrap_bin.2 (CAHJWO01): CAHJWO01VT1_04001 16S ribosomal RNA (partial)
* CAHJXG01: none

Added code to `code/find_dominant_bacteria.R` to make a mag table for the remaining unannotated Acetobacteraceae gen. sp. and CAHJXG01. Saved the table as `analysis/07_annotate_MAGs/extra_mag_table.txt`

Made a second Snakefile2 to make annotations for the extra MAGs (put them in an extra_ann folder)

```
snakemake --cores 10 -f -s Snakefile2
```

New rRNAs - not a single 16S
* private_T1867_metawrap_bin.2.log: CAHJXG01T1867_03424 5S ribosomal RNA
* private_X16_metawrap_bin.6.log: CAHJXG01X16_00904 5S ribosomal RNA
* public_SRR14722049_metawrap_bin.3.log: ACETOBACSRR14722049_02114 5S ribosomal RNA
* public_SRR14722049_metawrap_bin.4.log: CAHJXG01SRR14722049_00866  5S ribosomal RNA
CAHJXG01SRR14722049_00867 23S ribosomal RNA (partial)
* public_SRR14722071_metawrap_bin.2.log: ACETOBACSRR14722071n1_00857  5S ribosomal RNA
* public_SRR14722082_metawrap_bin.2.log: ACETOBACSRR14722082_01468 5S ribosomal RNA
* public_SRR14722143_metawrap_bin.2.log: ACETOBACSRR14722143_01499 23S ribosomal RNA (partial)
ACETOBACSRR14722143_01500 5S ribosomal RNA
* public_SRR14722165_metawrap_bin.2.log: CAHJXG01SRR14722165_00786  5S ribosomal RNA


13. Going manual with extracting 16S from Acto gen.sp and CAHJXG01, using approach described [here](https://github.com/strowig-lab/iMGMC/blob/master/linking/README.md)
Will continue with SRR13125477
```
mkdir analysis/06_16S_mag_matching/all_manual
cd analysis/06_16S_mag_matching/all_manual
conda install bbtools

mkdir SRR13125477_mags
cp ../../05_MAGs/MAGs/bacs/public_SRR13125477_metawrap_bin.* SRR13125477_mags/
cp ../../05_MAGs/MAGs/euks/public_SRR13125477* SRR13125477_mags/
gzip -d SRR13125477_mags/*
mkdir SRR13125477_split

#aling all reads to all MAGs from this metagenome and extract reads unique to each mag
bbsplit.sh in1=../SRR13125477_renamed.sra_1.fastq in2=../SRR13125477_renamed.sra_2.fastq ref=SRR13125477_mags basename=out_%.fq outu1=clean1.fq outu2=clean2.fq ambig2=split path=SRR13125477_split t=15

wc -l  out_public_SRR13125477_metawrap_bin.2.fq
> 1762032 out_public_SRR13125477_metawrap_bin.2.fq

# mapping the reads unique to out_public_SRR13125477_metawrap_bin.2.fq to 16S from this metagenome
bbmap.sh in1=out_public_SRR13125477_metawrap_bin.2.fq ref=../soil_Matam16S_wd/soil_assembled_16S_uclust_0.999.fasta t=15 minid=0.90 sortscafs=f nzo=f ambiguous=all local=t statsfile=SRR13125477_split/public_SRR13125477_metawrap_bin.2.statsfile \
scafstats=SRR13125477_split/public_SRR13125477_metawrap_bin.2.scafstats \
covstats=SRR13125477_split/public_SRR13125477_metawrap_bin.2.covstat \
rpkm=SRR13125477_split/public_SRR13125477_metawrap_bin.2.rpkm
#a handful of ref sequences have 1-2 reads mapped to them. The only 16S with 2 reads mapped to it blasts as lecanoro mitochondrial dna

#tried with another Aceto MAG
bbmap.sh in1=out_public_SRR13125477_metawrap_bin.5.fq ref=../soil_Matam16S_wd/soil_assembled_16S_uclust_0.999.fasta t=15 minid=0.90 sortscafs=f nzo=f ambiguous=all local=t statsfile=SRR13125477_split/public_SRR13125477_metawrap_bin.5.statsfile \
scafstats=SRR13125477_split/public_SRR13125477_metawrap_bin.5.scafstats \
covstats=SRR13125477_split/public_SRR13125477_metawrap_bin.5.covstat \
rpkm=SRR13125477_split/public_SRR13125477_metawrap_bin.5.rpkm
#same. The 16S with highest # of reads blasts as lecanoro 18S

#tried with metaxa 16s, for kicks
bbmap.sh in1=out_public_SRR13125477_metawrap_bin.2.fq ref=../../03_metagenome_reanalysis/assembly_SRR13125477.bacteria.fasta t=15 minid=0.90 sortscafs=f nzo=f ambiguous=all local=t statsfile=SRR13125477_split/public_SRR13125477_metawrap_bin.2_metaxa.statsfile \
scafstats=SRR13125477_split/public_SRR13125477_metawrap_bin.2_metaxa.scafstats \
covstats=SRR13125477_split/public_SRR13125477_metawrap_bin.2_metaxa.covstat \
rpkm=SRR13125477_split/public_SRR13125477_metawrap_bin.2_metaxa.rpkm
#nothing mapped2
```
Desisively doesn't work

14. Trying [Jorg](https://github.com/lmlui/Jorg)

* Downloaded MIRA mira_V5rc1_linux-gnu_x86_64_static.tar.bz2 from [here](https://github.com/bachev/mira/releases/tag/V5rc1)
```
cd /data/tagirdzh/bin/
tar -xf mira_V5rc1_linux-gnu_x86_64_static.tar.bz2
PATH=$PATH:/data/tagirdzh/bin/mira_V5rc1_linux-gnu_x86_64_static/bin
cd dbdata
./mira-install-sls-rrna.sh rfam_rrna-21-12.sls.gz
```
* settin up
```
conda create -n jorg
conda activate jorg
conda install -c bioconda seqtk
conda install bwa
conda install -c bioconda -c conda-forge barrnap
cd /data/tagirdzh/bin/
git clone https://github.com/lmlui/Jorg
PATH=$PATH:/data/tagirdzh/bin/Jorg
```
* trying out
* NB: here I switched to a different MAG from Acet gen.sp (private_X1_metawrap_bin.1) because it's has highest completeness (99%) and coverage (52X). Luckily X1 also has only one bacterial MAG!
* For CAHJXG01 the mag with hichest completenes (98%) and coverage (68X) is private_VT34_metawrap_bin.3. Unfortunately, VT34 has 2 CAHJXG01 MAGs, but the second has coverage ~ 6X. The other metagenome with high quality CAHJXG01 MAGs have more Aceto MAGs, so this still is the better option
```
cd /data/tagirdzh/coverage/analysis/06_16S_mag_matching/
mkdir jorg
cd jorg
cp ../../05_MAGs/MAGs/bacs/private_X1_metawrap_bin.1.fa.gz .
gzip -d private_X1_metawrap_bin.1.fa.gz
cp /data/tagirdzh/bin/Jorg/Example/manifest_template.conf .
jorg -b private_X1_metawrap_bin.1.fa  --forward ../../../data/fastq/X1.sra_1.fastq.gz  --reverse ../../../data/fastq/X1.sra_2.fastq.gz -k 33 -c 40 -i 5
```

error: private_X1_metawrap_bin.1.out.fasta is empty after the first iteration. will try with lower coverage threshold
```
jorg -b private_X1_metawrap_bin.1.fa  --forward ../../../data/fastq/X1.sra_1.fastq.gz  --reverse ../../../data/fastq/X1.sra_2.fastq.gz -k 33 -c 20 -i 5
```

Same result. The authors suggest the bin is too fragemented. Will parse PROKKA log files to see which MAG from the groups of interest has lowest # of contigs
```
grep "contigs totalling" */*.log
```
Aceto gen sp: private_X1_metawrap_bin.1 = 227 contigs, private_X11_metawrap_bin.3 = 153 contigs, private_VT1_metawrap_bin.3 = 245 contigs, public_SRR14722071_metawrap_bin.1 = 221 contigs

CAHJXG01: private_VT34_metawrap_bin.3 = 129 contigs, private_T1887_metawrap_bin.15 = 451 contigs, private_X15_metawrap_bin.2 = 74 contigs, public_SRR11456914_metawrap_bin.12 = 89 contigs

WIll try less fragmented bins
* SRR11456914 has 3 CAHJXG01 MAGs: public_SRR11456914_metawrap_bin.12 in 99.98% complete and has coverage 55X, while the others have 13X and 3X. Might work, though there is a LOT of bacterial MAGs here
* X15 has only one CAHJXG01 MAG, but it's only 15X might try too
* X11 has one Aceto gen sp MAG at 28X other MAGs from Aceto family are 6X and 1X. It's not dramatically better than private_X1_metawrap_bin.1 (153 vs 227 contigs) still will try - what else can I do?

```
mkdir public_SRR11456914_metawrap_bin.12
cd public_SRR11456914_metawrap_bin.12
cp ../../../05_MAGs/MAGs/bacs/public_SRR11456914_metawrap_bin.12.fa.gz .
gzip -d public_SRR11456914_metawrap_bin.12.fa.gz
jorg -b public_SRR11456914_metawrap_bin.12.fa  --forward ../../../../data/fastq/SRR11456914.sra_1.fastq.gz   --reverse ../../../../data/fastq/SRR11456914.sra_2.fastq.gz  -k 33 -c 20 -i 5
```
MIRA complains about long read names, renamed and tried again
```
cp ../../../../data/fastq/SRR11456914.sra_* .
 gzip -d SRR11456914.sra_*.fastq.gz
sed '/^@S/ s/ .*//' SRR11456914.sra_1.fastq > SRR11456914.sra_1_tmp.fastq
sed '/^@S/ s/ .*//' SRR11456914.sra_2.fastq > SRR11456914.sra_2_tmp.fastq
sed '/^+S/ s/ .*//' SRR11456914.sra_1_tmp.fastq > SRR11456914.sra_1_renamed.fastq
sed '/^+S/ s/ .*//' SRR11456914.sra_2_tmp.fastq > SRR11456914.sra_2_renamed.fastq
rm SRR11456914.sra_1.fastq SRR11456914.sra_2.fastq SRR11456914.sra_1_tmp.fastq SRR11456914.sra_2_tmp.fastq
jorg -b public_SRR11456914_metawrap_bin.12.fa  --forward SRR11456914.sra_1_renamed.fastq   --reverse SRR11456914.sra_2_renamed.fastq  -k 33 -c 20 -i 5
```
same problem. added COMMON_SETTINGS -NW:cmrnl=no to the manifest template -> this problem is solved

Still getting  Zero length fasta file 'public_SRR11456914_metawrap_bin.12.out.fasta'. mira_iteration_1.log says "Contig does not meet requirement of minimum reads per contig." Will try to map reads manually and see how many are mapped
```
conda deactivate
 conda activate bbtools
bbmap.sh in1=SRR11456914.sra_1_renamed.fastq in2=SRR11456914.sra_2_renamed.fastq ref=public_SRR11456914_metawrap_bin.12.fa t=15 minid=0.90 out=bbmap_public_SRR11456914_metawrap_bin.12.sam
awk '{if ($5 > 1)cnt[$3]++}END{for (x in cnt){print x,cnt[x]}}' bbmap_public_SRR11456914_metawrap_bin.12.sam
#lowest:
public_SRR11456914_metawrap_bin.12_88 904
#highest12
public_SRR11456914_metawrap_bin.12_1 261266
```



14. dNdS for bacterial genes to assess if they are pseudogenized
* Used `code/pre_dnds.R` to get fastas of certan KOs from ceratn taxa
    * Plan to run on all genes from Calvin cycle, bacteriochlorophyll pathway, and anoxyphotosystem II + pentose-phosphate pathway as a control
    * Prepared a table for all Acetos and Beijerinckos and all genes from these pathways
    * For starters picked 3 KOs: K00033 (PPP), K11532 (Calvin, is missing in RH-AL1), K08928 (photosystem)
    * For starters made three sets: all Acetos, all Bejerinkos, just Lichenihabitans, just CAHJXG01
    
* Followed [this protocol](https://www.protocols.io/view/introduction-to-calculating-dn-ds-ratios-with-code-qhwdt7e?step=4)
    * Produce alignments
    ```
   cd analysis/08_dNdS
   mafft --genafpair --maxiterate 10000 K08928.Lichenihabitans.ffn > K08928.Lichenihabitans.aligned.fna
    mafft --genafpair --maxiterate 10000 K08928.Lichenihabitans.faa > K08928.Lichenihabitans.aligned.faa
    ```
    * Get dN/dS ratios
    ```
    ~/bin/pal2nal.v14/pal2nal.pl K08928.Lichenihabitans.aligned.faa K08928.Lichenihabitans.aligned.fna -output paml -nogap > K08928.Lichenihabitans.pal2nal
    ``


15. Used `../../code/multi_level_screening_viz.R` to vizualize presence of key symbionts by three level of screening. Save the table showing presence/absence of rDNA of groups of interest as `03_metagenome_reanalysis/occurrence_rDNA_groups_of_interest.tsv`

16. Rhizobiales phylogeny
* Used dataset from [Volpiano et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8026895/#FS1)
	* `analysis/09_rhizobiales_phylogeny/Table_1.xlsx` is the supplementary from this paper. Table S1 gives the list of genomes they used
	* Saved the list of ftp addresses into a separate file `analysis/09_rhizobiales_phylogeny/ftp_list.txt`
	* Removed \r from the file
	```
	sed -i.bak 's/\r$//g' ftp_list.txt
	```
	* Downloaded protein and assembly fastas
	```
	mkdir genomes
	cd genomes
	while read p; do wget "$p"/*_genomic.fna.gz; done < ../ftp_list.txt
	rm *_cds_from_genomic.fna.gz
	rm *_rna_from_genomic.fna.gz
	
	mkdir ../annotations
	cd ../annotations
	while read p; do wget "$p"/*_protein.faa.gz; done < ../ftp_list.txt
	```
	
* Installed GTDB-Tk
	```
	conda create -n gtdb-tk
	conda install -c conda-forge -c bioconda gtdbtk
	download-db.sh
	```
	
* Copied Rhizobiales MAGs into the same folders
```
cd ../genomes
while read p; do cp ../../../"$p" . ; done < ../list_rhizobiales_mags.txt

cd ../annotations
while read p; do cp ../../../"$p" . ; done < ../list_rhizobiales_mags.txt
```
* Added outgroup from rhodobacter
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Rhodobacter_amnigenus/latest_assembly_versions/GCF_009908265.2_ASM990826v2/GCF_009908265.2_ASM990826v2_genomic.fna.gz
```

* Ran GTDB-Tk to identify and align marker genes
```
gtdbtk identify --genome_dir genomes --out_dir gtdbtk_identify --extension gz --cpus 10
gtdbtk align --identify_dir gtdbtk_identify --out_dir gtdbtk_align --cpus 20
gtdbtk classify --genome_dir genomes --align_dir gtdbtk_align --out_dir gtdbtk_classify -x gz --cpus 20

```

* Run IQTree
```
conda activate iqtree

#in narval
module load StdEnv/2020
module load gcc/9.3.0
module load iq-tree/2.1.2
iqtree -s ../gtdbtk_align/gtdbtk.bac120.user_msa.fasta --seqtype AA -bb 50000 -pre rhizobiales -seed 12345 -m TEST -T 12
```

* Screen protein fasta files to check if they have Nif genes
```
 grep "nitrogenase iron protein" ncbi_annotations/*faa > grep_nitrogenase_iron_protein.txt
```


* Do blast search on all genomes
	* used  `analysis/09_rhizobiales_phylogeny/Snakefile_tblastn` for searching nucleotide fastas
	```
	 snakemake --cores 10 -n -s Snakefile_tblastn 
	```
	
	* used  `analysis/09_rhizobiales_phylogeny/Snakefile_blastp` for searching protein fastas of *only ncbi genomes*
	* concatenated all into one file 
	```
	cat  blast_nifh/tmp_*_report* > blast_nifh/blast_nifh.txt
	```

* used `code/rhizobiales_nihf.R` to analyze results of NifH search:
	* methods of search (grep on the ncbi annotations vs tblastn vs blastp) are almost entirely consistent. 
	* The only inconsistency: in GCF_002879535.1 grep doesn't show NifH but tblastn&blastp do. The protein (WP_143973967.1) is labelled as nitrogenase reductase, partial in the ncbi annotations.
	* Will use tblastn search results for annotating the tree
	
	
	
17. Prepped data for metaGEM
* Used `code/metadata_for_metagem.R` to prepare metadata for all MAGs from VT1 and X11. Manually filled in info for taxonomy and QC scores for the eukaryotic MAGs
* COpied files
```
 mkdir MAG_fastas
 cd MAG_fastas
 perl -pe 's/\r\n|\n|\r/\n/g' ../lichen_metadata_metaGEM.txt > ../lichen_metadata_metaGEM_fixed_for_unix.txt
 awk 'NR!=1 {print $1} ' ../lichen_metadata_metaGEM_fixed_for_unix.txt > ../t
 while read f;  do cp ../../05_MAGs/MAGs/bacs/${f}* .; done < ../t
```
* Selected 15 chloro and 15 cyano lichens, for metaGEM analysis. Selected them by ranging metagenomes based on sequencing depth
	* saved list of metagenomes `analysis/10_metaGEM/selected_mtg.txt`, list of mags in `analysis/10_metaGEM/mags_in_selected_mtg.txt`
* Copied files
```
 mkdir MAG_fastas_july2022
 cd MAG_fastas_july2022
 perl -pe 's/\r\n|\n|\r/\n/g' ../mags_in_selected_mtg.txt > ../mags_in_selected_mtg_fixed_for_unix.txt
 awk 'NR!=1 {print $1} ' ../mags_in_selected_mtg_fixed_for_unix.txt > ../t
 while read f;  do cp ../../05_MAGs/MAGs/bacs/${f}* .; done < ../t
  while read f;  do cp ../../05_MAGs/MAGs/euks/${f}* .; done < ../t
```

18. TO REPLACE METAXA used `code/idtaxa_compiling_db.R` to assign taxonomic positions to bacterial rDNA (assemblies and reads) based on GTDB
*  modified  `../../code/multi_level_screening_viz.R` . used metaxa for eukaryotes and idtaxa for bacteria
* added `../../03_metagenome_reanalysis/occurrence_rDNA_stats.tsv` showing % of metagenomes given lineage was detected with given method

19. Three-level screening of selected lineages
* Used `code/multilevel_screening.R` to make a heatmap and a table for the presence of key lineages (4 bacterial and 3 eukaryotic)
* Used `code/multilevel_screening_less_freq.R` to check how less frequent bacteria (the next 6 families) compare in their frequencies

20. Geopgraphy and metagenomes
* Used `code/get_coordinates.sh` to get coordinates from ENA. Saved temp files as `PRJ*_locations.txt`, concatenated file as `locations.txt`
* Cleaned up and put the new table in `locations_compiled_manually.txt`. added there coordinates from other sources:
	* from papers that didn't have coordinates in ENA, but had it elsewhere
	* from de novo metagenomes
	* picked coordinates based on the text description of the location, in the cases where it's not available
	* one metagenome doesn't have coordinates (protected species)
	* coordinate formats aren't stnadardized!
* Used `code/map.R` to proccess the coordinates and plot the map
	* saved final coordinates in `analysis/03_metagenome_reanalysis/locations_final.txt`
	* saved `analysis/03_metagenome_reanalysis/map.pdf`

21. Fishing expedition into the algal MAGs
* started a new folder `analysis/11_algal_MAGs`
* check first on one MAG `private_T1889_concoct_bin.16` (Coccomyxa according to BAT; 98% complete) using a MetH from NCBI (BAU71143.1)
```
cp ../05_MAGs/MAGs/euks/private_T1889_concoct_bin.16.fa.gz .
gzip -d private_T1889_concoct_bin.16.fa.gz 
tblastn -query meth_genbank.fa -subject private_T1889_concoct_bin.16.fa -outfmt 6 -evalue 1e-5 <-nothing
```
* check first on one MAG `public_SRR7232211_concoct_bin.8` (Trebouxia according to BAT; 95% complete) using a MetH from NCBI (BAU71143.1)
```
cp ../05_MAGs/MAGs/euks/public_SRR7232211_concoct_bin.8.fa.gz .
gzip -d public_SRR7232211_concoct_bin.8.fa.gz 
tblastn -query meth_genbank.fa -subject public_SRR7232211_concoct_bin.8.fa -outfmt 6 -evalue 1e-5 -out tblastn_meth_public_SRR7232211_concoct_bin.8.txt
```
result: found many hits, at close inspection all but one are split between 3 contigs, probably just hits to the same gene? Isolated one "area", which encompassed many hits
```
samtools faidx public_SRR7232211_concoct_bin.8.fa  public_SRR7232211_concoct_bin.8_2872:376-2004
```
Blasted against NCBI, got hits to MetH. 

* Now will try with MetE (BAU71146.1), to check maybe the hits are going to be the same?
```
tblastn -query mete_genbank.fa -subject private_T1889_concoct_bin.16.fa -outfmt 6 -evalue 1e-5 
BAU71146.1      private_T1889_concoct_bin.16_73 60.000  80      32      0       668     747     7618    7857    2.59e-24        110
BAU71146.1      private_T1889_concoct_bin.16_73 79.688  64      13      0       498     561     6368    6559    1.55e-23        107
BAU71146.1      private_T1889_concoct_bin.16_73 62.319  69      23      2       235     303     4479    4676    1.64e-17        88.2
BAU71146.1      private_T1889_concoct_bin.16_73 61.039  77      29      1       85      161     3449    3676    3.99e-17        87.0
BAU71146.1      private_T1889_concoct_bin.16_73 75.510  49      12      0       606     654     7214    7360    1.48e-15        82.0
BAU71146.1      private_T1889_concoct_bin.16_73 51.190  84      34      2       161     237     3949    4200    5.35e-15        80.1
BAU71146.1      private_T1889_concoct_bin.16_73 69.388  49      15      0       557     605     6823    6969    3.63e-12        70.9
BAU71146.1      private_T1889_concoct_bin.16_73 44.262  61      24      1       437     497     5958    6110    1.79e-06        52.0
samtools faidx private_T1889_concoct_bin.16.fa private_T1889_concoct_bin.16_73:3949-7857 

tblastn -query mete_genbank.fa -subject public_SRR7232211_concoct_bin.8.fa -outfmt 6 -evalue 1e-5 
BAU71146.1      public_SRR7232211_concoct_bin.8_350     56.842  95      37      3       160     251     19779   20060   5.37e-22        103
BAU71146.1      public_SRR7232211_concoct_bin.8_350     69.118  68      21      0       498     565     22301   22504   1.22e-20        99.0
BAU71146.1      public_SRR7232211_concoct_bin.8_350     74.510  51      13      0       605     655     23100   23252   8.33e-18        89.7
BAU71146.1      public_SRR7232211_concoct_bin.8_350     68.182  66      21      0       57      122     18420   18617   5.78e-14        77.4
BAU71146.1      public_SRR7232211_concoct_bin.8_350     58.929  56      23      0       707     762     23943   24110   2.50e-13        75.1
BAU71146.1      public_SRR7232211_concoct_bin.8_350     82.051  39      7       0       668     706     23570   23686   2.20e-12        72.0
BAU71146.1      public_SRR7232211_concoct_bin.8_350     58.491  53      20      1       251     303     20463   20615   5.68e-11        67.4
BAU71146.1      public_SRR7232211_concoct_bin.8_350     64.444  45      16      0       561     605     22709   22843   2.25e-10        65.5
BAU71146.1      public_SRR7232211_concoct_bin.8_350     43.590  78      33      2       420     497     21780   21980   4.56e-08        57.8
samtools faidx public_SRR7232211_concoct_bin.8.fa public_SRR7232211_concoct_bin.8_350:18420-24110
```
results: hits are not the same (as to be expected actually, as metH and metE have no sequence similarity, see [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC539065/))! both MAGs have hits, in both all hits come from one contig, evalue is lower than for the hit from MetH-Trebouxia

When blasted against NCBI, hits come from MetE (=5-methyltetrahydropteroyltriglutamate--homocysteine methyltransferase)

* prelim conclusion: the one Coccomyxa mag I tried has only metE, the one Trebouxia has both metE and metH (Kroft et al. [show](https://www.nature.com/articles/nature04056) that some algae can have both, and use metE in the absence of B12)
	* are the MAGs that only have metH?
*  screened all algal mags
	* used `get_algal_mag_list.R` to get algal mags with >90% completeness
	* copied all these mags into the folder
	```
	cat good_algal_mags_list.txt | xargs -I {} cp ../05_MAGs/MAGs/euks/{}.fa.gz .
	gzip -d *.gz
	```
	* made Snakemake pipeline ` analysis/11_algal_MAGs/Snakefile_tblastn`
	* results: of 19 MAGs, all had MetE (B12-independent), 13 had MetH (B12-dep)