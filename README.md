# A global survey of lichen symbionts from metagenomes 
This repository contains scripts and intermediate results for the manuscript (Tagirdzhanova et al. 2022, biorXiv)

## Abstract


## Overview

```
project
├── README.md							# this doc; description of the repo and the project log
├── code 								# all scripts generated for the analysis, with the exception of snakemake pipelines (those can be found in exploratory/)
├── analysis 							# exploratory analysis, trees and tables generated for the analysis. Only folders relevant for this publications are included. Some files are designated as Supplementary data (see below)
│   ├── 03_metagenome_reanalysis			# information related to the used metagenomes (metadata, SRA IDs, location information) and analysis on the metagenome-level (i.e. rDNA screenening) 
│   ├── 05_MAGs 							# analysis on the level of MAGs: phylogenomic trees, tables related to MAG occurrences and coverage, and exploratory figures
│   ├── 07_annotate_MAGs					# snakemake pipelines for annotating selected bacterial MAGs and summarized outputs of annotations: KO annotations summarized by bacterial genus, CAZy annotations
│   ├── 09_rhizobiales_phylogeny			# analysis of Rhizobiales MAGs and reference genomes: phylogenmoci trees, and blast-based screening for several genes associated with C1 metabolism and nitrogen fixation
│   ├── 11_algal_MAGs						# screening of algal MAGs for methionin synthases MetE and MetH
│   └── Notebook 							# temp log files for various parts of the analysis; the cleaned-up version of the same logs is in this file
└── results 							# results included in the manuscript
    ├── figures 							# figures
    └── tables 								# tables, included as Supplementary tables and key data tables used to produce figures

```

## 1. Dataset construction
Software used:
* [sratoolkit](https://github.com/ncbi/sra-tools/wiki)
* sourmash v4.2.2 (Pierce et al., 2019)
* R libraries: tidyverse, stringr, ggmap, rnaturalearth, sf, patchwork, scales

### 1.1. Prepared dataset:
* Got SRA ids from Lendemer et al. 2019:
    * Downloaded SRA run info for PRJNA700635 and PRJNA731936 > `analysis/03_metagenome_reanalysis/SraRunInfo_PRJNA731936.csv analysis/03_metagenome_reanalysis/SraRunInfo_PRJNA700635.csv`
    * Downloaded [Appendix S2](https://bsapubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fajb2.1339&file=ajb21339-sup-0002-AppendixS2.xlsx) from Lendemer et al. 2019 with the voucher metadata for their metagenomes
    * Made finel list of SRA IDs using `code/getting_SRA_id_for_lendemer_data.R`. Only kept IDs that matched between the metadata from the NCBI and from the Appendix > `analysis/03_metagenome_reanalysis/sra_ids_lendemer.txt`

* Got SRA ids from other sources
    * See `analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.csv`. Copied all SRA ids into `analysis/03_metagenome_reanalysis/sra_ids_other.txt`
    * Combined the two files into the final list of SRA IDs to use:
```
cat analysis/03_metagenome_reanalysis/sra_ids_lendemer.txt analysis/03_metagenome_reanalysis/sra_ids_other.txt >analysis/03_metagenome_reanalysis/sra_ids_all.txt
```
* Added de novo generated metagenomic data (see Table XXX)

### 1.2. Downloaded data using sratoolkit
```
cd data/fastq
prefetch --option-file analysis/03_metagenome_reanalysis/sra_ids_all.txt

for file in sratoolkit-cache/sra/*sra

do
fasterq-dump --split-files "$file"
done

gzip *.fastq

```
### 1.3. Filtered dataset to remove duplicate metagenomes

* Using sourmash detected that some SRA libraries are identical despite having different IDs. All duplicated are shared between two studies, PRJNA731936 and PRJNA700635. The list is in `analysis/03_metagenome_reanalysis/similar_datasets.csv`, in total there are 42 identical pairs.
* Removed the duplicated libraries from the smaller PRJNA700635 study. The other copies remain in the analysis
```
for ID_SAMPLE in `cut -d ',' -f1 analysis/03_metagenome_reanalysis/similar_datasets.csv`; do rm data/fastq/"$ID_SAMPLE"*; done
```

### 1.4. Compled information on the metagenomes
* Saved the table as `results/tables/all_metagenome_reanalysis.txt`
* This table is compiled manually, using information from the NCBI metadata, and literature on lichens
* In the text it's reffered to as **Table SXXX**

### 1.5. Produced **Fig. SXXX**: the map of sampling locations for all metagenomes used in the study
* Used `code/get_coordinates.sh` to obtain location information from ENA
* To the samples that didn't have coordinates, I assigned coordinates based on the description of the sampling location in the original paper
* Cleaned-up version of the output is saves as `analysis/03_metagenome_reanalysis/locations_compiled_manually.txt`
* Used `code/map.R` to draw map of all sampling locations
* The map is saved as `results/figures/map.pdf`, in the text it's reffered to as **Fig. SXXX**

### 1.6. Produced **Fig. SXXX**: dot hitogram of metagenomes arranged by sequencing depth
* Used `code/metagenome_by_source_fig.R`
* Saved the figure as `results/figures/metagenome_dothist_by_source.png`


## 2. Metagenomic assembly and binning
Software used:
* fastp (Chen et al., 2018)
* metaWRAP pipeline v.1.2 (Uritskiy et al. 2018)
* metaSPAdes (Nurk et al. 2017
* CONCOCT (Alneberg et al. 2014)
* metaBAT2 (Kang et al. 2015)
* CheckM v1.1.3 (Parks et al. 2015)
* dRep v3 (Olm et al. 2017)

```
XXX to be added after asking Paul for details
```

## 3. Taxonomic assignments of MAGs
Software used:
* GTDB-Tk v1.5.0 (Chaumeil et al. 2020)
* IQ-TREE (Nguyen et al. 2015)
* BAT (CAT v5.2.3, database version: 20210107; von Meijenfeldt et al. 2019)
* iTOL (Letunic & Bork 2019)

### 3.1. Prokaryotic MAGs
```
XXX to be added after asking Paul for details
GTDB-Tk
IQTREE
```
* The tree is available as Supplementary data: XXX

### 3.2. Eukaryotic MAGs
```
XXX BAT analysis: to be added after asking Paul for details
```
Computed two phylogenomic trees: one for fungi, one for algae. 
	* The list of reference genomes is in `results/tables/reference_genomes_full_table.csv`
	* This table is compiled manually, in the text it's referred as **Table SXXX**
```
XXX euk phylogenomic: to be added after asking David for details
```

### 3.3. Produced figure Fig. XXX: eukaryiotic phylogenomic trees


## 4. Occurrence analysis
Software used:
* BWA (Li & Durbin 2009)
* Metaxa2 (Bengtsson‐Palme et al. 2015)
* [BBTools](https://sourceforge.net/projects/bbmap/)
* IDTAXA (Murali et al. 2018)
* iTOL (Letunic & Bork 2019)
* Snakemake (Mölder et al. 2021)
* R libraries: tidyverse, ape, phytools, stringr, taxize, myTAI, igraph, qgraph, plotly, DECIPHER, R.utils, treeio, seriation, ComplexHeatmap, DECIPHER, circlize, conflicted

### 4.1. Aligned all metagenomic datasets to all MAGs
```
XXX to be added after asking Paul for details

```
### 4.2.-4.3 Assigned MAGs putative roles and removed potentially misidentified samples
* Combined taxonomy for prokaryotic MAGs (i.e. outputs of GTDB-Tk) and eukaryotic MAGs (i.e. outputs of BAT) in one table
	* Used `code/combine_mag_taxonomy_annotation.R`
	* Saved the table as `analysis/05_MAGs/tables/MAG_taxonomy_combined.tsv`
* Assigned putative roles in the symbiosis based on the BAT taxonomy and coverage
	* Used `code/assign_putative_mag_roles.R`
	* As input used the metagenome metadata `results/all_metagenome_reanalysis.txt` and BWA alignments `analysis/05_MAGs/tables/read_mapping/`
	* Saved the table with putative assignments as `analysis/05_MAGs/tables/MAG_putative_roles_bwa.tsv`
	
* Manually corrected assignments of the LFS (a.k.a. 'mycobiont') MAGs based on the phylogenomic tree
	* Used `rename_euk_tree.R` rename the tips in the phylogenomic trees to improve readability by adding putative assigments and the names of lichen symbioses.
		* saved the renamed fungal tree as `analysis/05_MAGs/trees/eukaryotes/Fungi/phylogeny/Concatenated_IQTREE/concat_putative_renamed.contree`
		* saved the renamed algal tree as `analysis/05_MAGs/trees/eukaryotes/Algae/phylogeny/Concatenated_IQTREE/concat_renamed.contree`
	* Analyzed the phylogenomic trees to check is there are inconsistencies between the lichen names and the LFS MAG placement on the tree 
	* manually corrected assignments in several metagenomes:
		* Changed LFS in SRR14722059 (tuckermanopsis) public_SRR14722038_metabat2_bin.1 -> public_SRR14722026_metabat2_bin.4, because it's grouped with other tuckermannopsis'es
		* Changed SRR14722135 (Mycocalicium) LFS assignemnet from public_SRR14722135_concoct_merged.1 to public_SRR14722135_metabat2_bin.8 (lower coverage but groups with the right group)
		* Changed SRR14722098 (Acarospora sinopica) LFS from public_SRR14722098_metabat2_bin.9 to public_SRR14722090_metabat2_bin.14 (lower coverage but groups with the right group)
	* Checked against NCBI to confirm that the organism as listed in the NCBI metadata, FEN number, and placement in the tree are consistent. Made changes:
		* SRR14722289 changed Diploschistes (acc. to NCBI metadata)-> Parmotrema (acc. to the NCBI name, FEN number and the tree)
		* SRR14721950 changed Myalospora (acc. to NCBI name and metadata) -> Pseudosagedia (acc. to FEN number and the tree)
	* Identified problematic cases, that couldn't be resolved. All of them have "mycobiont_missassigned" in the `results/tables/MAG_confirmed_roles_bwa.tsv` table
		* SRR14722032: Platismatia tuckermanii, the only fungal MAG groups with Ochrolehcia/Pertusaria/Lepra
		* SRR14722092: Chrysotrix, grouped with trapeliopsis/gomphillus
		* SRR14722327: Catillaria, grouped with Leprocaulon
		* SRR14722303: Rinodina, groups with Lecanora
		* SRR14722033: Lecanora, groups with Parmeliaceae
		* SRR14722324: Lepraria,  grouped with Leprocaulon
		* SRR14722131: Parmotrema, groups with Thelotrema
		* SRR14722208: Ricasolia, groups with Parmeliaceae
		* SRR14722229: Lecidea, groups with Lecanora
		* SRR14722034: Lepraria,  grouped with Parmeliaceae
		* SRR14722085: Psedosagedia, grouped with theloschistales
		* SRR14722185: Mycobilimbia, groups with byssoloma
		* SRR14722222: Arthonia, groupswith eurotios
		* SRR14722160: Herteliana, groups with Lepraria
	* used `code/make_list_excluded_mtg.R` to save this list as a table in `results/tables/excluded_metagenomes.txt`
		* This table is referred to in the text as **Table SXXX**
	* saved the updated assignments as `results/tables/MAG_confirmed_roles_bwa.txt`

**NB:** This table is one of the key tables used for producing figures and tables downstream. Each line represents a MAG occurrence, i.e. each instant a MAG is present in a metagenome. It provides info on: 
* the MAG: size, taxonomy 
* the metagenome: what lichen symbiosis it's made from)
* the occurrence: breadth and depth of coverage of this MAG in this metagenome
	
### 4.5. Identifying most frequent bacterial groups
* Used `code/find_dominant_bacteria.R` to identify bacterial lineages of interest
* Ranked bacterial groups by their frequency, defined as the total number of occurrences
	* Summarized on four levels: individual MAGs, bacterial genera, families, and orders
	* Saved the resulting tables as `analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_*_frequency.tsv`. Here I counted total number of occurrences per bacterial group 
* Ranked bacterial groups by their diversity, defined as the number of unique MAGs
	* Summarized on three levels: bacterial genera, families, and orders
	* Save the resulting tables as `analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_*_diversity.tsv`. Here is just a number of MAGs from a given group
* Only included the metagenomes that a) yielded an LFS MAG b) wasn't removed as potential misidentification 
* Made lists of top frequency genera and families by groups of lichens: photobionts, mycobionts, and compinations.  
	* saved them as `results/tables/bacterial_*_by_lichen_group.txt`
	* These two tables combined are referred to in the text as **Table SXXX**
* Visualized bacterial occurrences
	* Used `code/rename_bac_tree.R`
	* Renamed tip labels in the bacterial phylogenomic tree produced (see 3.1) to add taxonomic assignments
	* Produced tables with annotations: phylum-level taxonomy and the number of occurrences per MAG. The tables are desinged to be compatible with iTOL
	* Produced the figure using iTOL, the figure is reference in the table as **Fig. XXX**
	
### 4.6. Co-occurrence analysis
* Plot cooccurrence networks for selected groups:
	* LFS + eukaryotic groups (green algae, Cystobasidiomycetes and Tremellomycetes fungi) + Cyanobacteria
	* LFS + top genus in Beijerinckiaceae (Lichenihabitans)
	* LFS + 2 top genera in Acetobacteraceae (CAHJXG01 and LMUY01)
	* LFS + 3 top genera in Acidobacteriaceae (Terriglobus, EB88, and CAHJWL01)
* Only included the metagenomes that a) yielded an LFS MAG b) wasn't removed as potential misidentification 
* Cooccurrence is defined as two MAGs occurring together in one metagenome
* The graphs only show LFS MAGs and MAGs from the selected groups. Each node is a MAG, edges represent cooccurrences. The thicker the edge, the more frequently the two MAGs cooccur
* Used `code/cooc_graph_*.R`
* Saved the figures as `results/figures/coocc_graph_.pdf`
* The figures are referenced in the text as **Fig. XXX**

### 4.7. rDNA-based screening
* Screened metagenomic assemblies and raw unassembled reads for the rDNA sequences
* Made a Snakemake pipeline `analysis/03_metagenome_reanalysis/Snakefile`. This pipeline:
	* Ran Metaxa2 on all metagenomic assemblies and raw read sets. Metaxa uses a HMM-based algorithm to detect all rDNA sequences, and then compares them to the DIAMOND database for the taxonomic classification. Output is saved in the form `reads_{sample}.level_5.txt`
	* Ran reformat.sh (BBTools) to calculate sequencing depth for each sample
	* Summarized all sequencing depths in one table (`analysis/03_metagenome_reanalysis/bp_report.txt`). This table is used in several scripts, including `code/metagenome_by_source_fig.R`
	* Ran stats.sh (BBTools) to calculate the assembly length for each sample
	* Summarized all assembly lengths in one table (`analysis/03_metagenome_reanalysis/assembly_report.txt`)
```
cd analysis/03_metagenome_reanalysis/
snakemake --cores 10 
cd ../../
```
* Compiled all Metaxa2 reports into one big rDNA-based occurrence matrix. Used the metaxa_dc module
```
metaxa2_dc *.level_5.txt -o metaxa_level_5_combined.txt
```
* Used IDTAXA to reclassify bacterial sequences
	* While, I kept Metaxa classifications for eukaryotes, I re-classified bacterial rDNA sequences with IDTAXA
	* Rationale: IDTAXA uses GTDB system, and therefore IDTAXA assignments would be consistent with the system we used for MAG taxonomic assigments
	* Used `code/idtaxa_compiling_db.R` to set up the IDTAXA database and to reclassify bacterial sequences. As input, used fasta files produced by Metaxa `analysis/03_metagenome_reanalysis/{sample}.bacteria.fasta.txt`
	* Saved IDTAXA output as `analysis/03_metagenome_reanalysis/idtaxa_reads.txt` and `analysis/03_metagenome_reanalysis/idtaxa_assemblies.txt`
* Visualized presence/absence of selected lineages according to the different types of screening
	* Used `code/multi_level_screening_viz.R` to make a heatmap and a table for the presence of key lineages (4 bacterial and 3 eukaryotic) in metagenomes
	* The heatmap shows detection of the lineages on three "levels": as MAGs, as rDNA in the assemblies, and as rDNA in the reads
	* Only included the metagenomes that a) yielded an LFS MAG b) wasn't removed as potential misidentification 
	* The dendrogram  shows the fungal phylogenomic tree (see 3.2). Tips that are not LFS are dropped from the tree
	* Saved the figure as `results/figures/multilevel_screening.pdf`. The figure is referenced in the text as **Fig.XXX**
	* Made a table with the prevalence (i.e. percentage of metagenomes that contained a group of organisms) of the selected groups, according to each "level" of detection. Save the table as `analysis/03_metagenome_reanalysis/occurrence_rDNA_stats.tsv`
	* Calculated prevalence of Lichenihabitans using the same approach as above
* Calculated prevalence for all bacterial families, based on the rDNA screening of assemblies and reads
	* Only included the metagenomes that a) yielded an LFS MAG b) wasn't removed as potential misidentification 
		* If all metagenomes are included, the ranking stays the same
	* Saved the table as `results/tables/idtaxa_top_bac_families.tsv` 
	* This table is referenced in the text as **Table SXXX**

## 5. Estimating relative abundance of symbiont groups
Software used:
* R libraries: tidyverse, scales

**NB:** some of results in this subsection are relevant to the next section (6)

### 5.1. Analyzed the number and total coverage per MAG category in each metagenome
* MAG 'categories' correspond to the putative roles (see 4.2)
* Used `code/summarize_mag_cov_counts.R`  
* As a source, used `results/tables/MAG_confirmed_roles_bwa.tsv`
* Calculated the number of MAGs per category, saved the results as `results/tables/MAG_counts_summary.tsv`. This table is referenced in the text as **Table SXXX**
* Calculated the total coverage depth per  MAGs per category, saved results into `results/tables/MAG_coverage_summary.tsv` 
* Saved metagenomes that had an LFS MAG, but its coverage was lower than combined coverage of bacteria in `analysis/05_MAGs/tables/MAG_coverage_high_bacteria.tsv`
* Produced a plot showing the relationship between the number of MAGs and sequencing depth, saved as `results/figures/mags_vs_depth.png`. This figure is referenced in the text as **Fig. SXXX**

### 5.2. Analyzed relative abundance for the key bacterial genera
* Abundance inferred from the coverage depth. Relative abundance is relative to the LFS abundance
* Relative abundance of a lineage in a sample = (coverage depth of its MAG)/(coverage depth of the LFS MAG) 
* Only the metagenomes that yielded an LFS MAG were used
* Used `code/relative_cov_bac.R` to plot relative abundance of the key symbiont groups: 13 top bacterial genera and the tree eukaryotic symbiont groups
* Saved as "results/figures/relative_cov_boxplot_genus.png". This figure is referenced in the text as **Fig. SXXX**
* Saved the table that contains median relative abundanceper bacterial genus as `results/tables/median_relative_coverage_by_genus.tsv`. This table is referenced in the text as **Table SXXX**

## 6. Analyzing how sequencing depth affects MAG recovery
Software used:
* R libraries: tidyverse, stringr, patchwork

Analyzed the recovery of the LFS and photosynthetic partner MAGs as a function of sequencing depth
* Used `code/recovery_myco_photo_MAGs.R`
* Used the functional assignments from `results/tables/MAG_confirmed_roles_bwa.tsv`
* Only included MAGs extracted by binning (including those discarded during dereplication). 
* Didn't count a MAG as present if it was only detected via BWA alignments. Rationale: a genome can be present in a metagenome but isn't recovered as MAG due to e.g. low coverage. In this context, count this genome as present does not agree with the goal of this analysis (i.e. what sequencing depth is efficient to produce an LFS/photosynthetic partner MAG?)
* Saved the figure `results/figures/myco_photobiont_vs_depth.png`. This figure is referenced in the text as **Fig. SXXX**

## 7. Functional analysis

