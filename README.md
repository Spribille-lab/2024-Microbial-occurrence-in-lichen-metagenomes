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
    * Made finel list of SRA IDs using `code/getting_SRA_id.R`. Only kept IDs that matched between the metadata from the NCBI and from the Appendix > `analysis/03_metagenome_reanalysis/sra_ids_lendemer.txt`

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
* Preliminary taxonomic assignments based on the database search
```
XXX BAT analysis: to be added after asking Paul for details
```
* To refine taxonomic assignments, computed two phylogenomic trees: one for fungi, one for algae. 
	* The list of reference genomes is in `results/tables/reference_genomes_full_table.csv`
	* This table is compiled manually, in the text it's referred as **Table SXXX**
```
XXX euk phylogenomic: to be added after asking David for details
```

* Produced figure Fig. XXX: eukaryiotic phylogenomic trees
	* Used iTOL
	* Manually annotated trees taxonomically, based on the reference genomes


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
	* the metagenome: what lichen symbiosis it's made from
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
Software used:
* PROKKA v1.13 (Seemann 2014)
* Snakemake (Mölder et al. 2021)
* KEGG Orthology Database (Kanehisa et al. 2002)
* KofamScan (Aramaki et al. 2020)
* BLAST (Altschul et al. 1990)
* [dbcan2 v3.0.2](https://github.com/linnabrown/run_dbcan)
* FeGenie (Garber et al. 2020)
* Emerald
* [getLCA](https://github.com/frederikseersholm/getLCA)
* R libraries: tidyverse, stringr, seriation, ComplexHeatmap, DECIPHER, circlize, RColorBrewer, patchwork, scales, waffle, extrafont, hrbrthemes, simplifyEnrichment, 

### 7.1. General functional annotations and Protein space analysis
* Annotated all bacterial MAGs with PROKKA

```
prokka --compliant --centre UoA --outdir {MAG_ID} --prefix {MAG_ID} {MAG_ID}.fa
```

* Compared all predicted proteins against the MGnify database
```
Get more details from Ellen
```

### 7.2. Functional clustering
* Annotated predicted proteins against the KEGG Orthology Database using KofamScan
```
/bin/kofam_scan/exec_annotation -o {MAG_ID}/{MAG_ID}.kegg.mapper.txt  {MAG_ID}/{MAG_ID}.faa -f mapper --tmp-dir tmp_{MAG_ID}
```

* Modified script from Zoccarato et al. (2022) to reconstruct KEGG modules
	* The script is saved as `code/kegg_module_reconstruct.R`
* The main script is `code/bacteria_functional_clustering.R`. Using this script did the following:
	* Combined all KEGG annotations in one table `analysis/07_annotate_MAGs/summarized_outputs/all_mags_kegg_combined.txt`
	* For the MAGs with >90% completeness, reconstructed KEGG modules, saved the table as `analysis/07_annotate_MAGs/summarized_outputs/kegg_module_matrix.txt`
	* Compared different clustering methods (kmeans, dynamicTreeCut, mclust, apcluster, hdbscan, fast_greedy, leading_eigen, louvain, walktrap, MCL, binary_cut). Saved the output figures as `analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_comparison.pdf` and `analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_comparison2.pdf` (referred in the text as **Fig. SXXX**)
	* For the selected methods of clustering (kmeans, apcluster, hdbscan), saved results (i.e. cluster assignments for all used MAGs), as a table (`analysis/07_annotate_MAGs/summarized_outputs/kegg_modules_clustering.txt`) and as a figure (`analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_results.pdf`)
	* Calculated taxonomic coherence for the selected three methods (following the definitions from Zoccarato et al. 2022). Saved the resulting figure as `analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_modules_clustering_taxonomic_coherence.pdf` (referred in the text as **Fig. SXXX**)

### 7.3. Selected bacterial MAGs for in depth annotation
* Picked the 13 bacterial genera with highest numbers of occurrences, that together accounted for 53% of all bacterial occurrences
* From them, selected all highest quality MAGs (completenes >95%, contamination <10%, according to CheckM)
* Used `code/select_mags_for_annotaiont.R`
* Saved the list of selected MAGs as `analysis/07_annotate_MAGs/mag_table.txt` and `analysis/07_annotate_MAGs/mag_list.txt`. In total, 63 MAGs were selected
* Prepared the table with information on the selected MAGs, saved as `results/tables/selected_mag_info.txt`. This table is referenced in the text as Table SXXX
* Made figures to illustrate the process of MAG selection
	* Used `code/draw_mag_selection_fig.R`
	* Saved the figure as `results/figures/mag_selection.svg`, it is referred in the text as **Fig. XXX**

### 7.4. In-depth functional annotations
Analyzed the 63 selected MAGs

* Created a KEGG table for the selected MAGs
	* Used `code/combine_kegg.R`
	* Saved KEGG annotations for all selected MAGs as `analysis/07_annotate_MAGs/summarized_outputs/selected_mags_kegg_combined.txt`
	* Saved combined KEGG annotations for each genus as `analysis/07_annotate_MAGs/summarized_outputs/{genus}.kegg.combined.txt`. Those were used for exploratory analysis only
* Visualized key metabolic traits
	* Screened KEGG annotations for the annotation related to the key metabolic traits
		* Anoxygenic Photosystem II (K08928, K08929, K13991, K13992, K08926, K08927)
		* Calvin cycle (K00855, K01601, K01602, K00927, K00134, (K01623 or K01624), K00615, K03841, (K01807 or K01808))
		* C1 metabolism: methanol dehydrogenase (K23995) and methane monooxygenase (K10946 and K16157)
		* Bacteriochlorophyll biosynthesis pathway (K04035, K04037, K04038, K04039, K11333, K11334, K11335, K11336, K11337, K04040, K10960)
		* Carotenoid biosynthesis pathway (K02291, K10027, K09844, K09844, K09845, K09846)
		* Nitrogenase (K02588)
		* Cofactor biosynthesis pathway
			* Cobalamin ((K00768 or K02226) or (K02232, K02231, (K02227 or K02225), (K00798 or K19221))
			* Biotin (K00652, K00833, K01935, K01012)
			* Riboflavin ((K00794, K00793, K11753) or((K01497 or K14652), K11752, K21064))
			* Thiamine (K00878, K00941, K00788)
		* Urease ((K14048 or (K01430, K01429)), K01428)
		* Transport systems
			* sorbitol/mannitol transporter (K10227, K10228, K10229, K10111)
			* urea transporter (K11959, K11960, K11961, K11962, K11963)
			* erythritol transporter (K17202, K17203, K17204)
			* xylitol transporter (K17205, K17206, K17207)
			* inositol  transporter (K17208, K17209, K17210)
  			* glycerol  transporter (K17321 K17322, K17323, K17324, 17325)
			* fucose transporter (K02429) 
			* glycerol aquaporin transporter (K02440)
			* glycerol/sorbitol transporter (K02781, K02782, K02783)
			* ammonium transporter (K03320) 
			* ribose transporter (K10439, K10440, K10441)
			* xylose transporter (K10543, K10544, K10545)
			* multiple sugar transporter (K10546, K10547, K10548)
			* fructose transporter (K10552, K10553, K10554)
			* arabinose transporter (K10537, K10538, K10539)
			* branched-chain amino acid transporter (K01999, K01997, K01998, K01995, K01996)
			* L-amino acid transporter (K09969, K09970, K09971, K09972)
			* glutamate transporter (K10001, K10002, K10003, K10004)
			* capsular transporter (K10107, K09688, K09689)
	* Used `code/pathway_fig.R`
	* Saved the figure as `results/figures/pathway_fig.svg`, it is referenced in the text as **Fig. XXX**
	 
* Screened the metagenomic assemblies for NifH, a gene involved in nitrogen fixation
	* Made a Snakemake pipeline: `analysis/07_annotate_MAGs/Snakefile_tblastn_metagenome`
	* Searched the metagenomic assemblies using tblastn, as a query used a sequence from NCBI (ABZ89802.1)
	* Extracted the hits as fasta files
	* Obtained taxonomic assignments for the hits by searching them against the NCBI_nt database and using getLCA
* Annotated CAZymes using run_dbcan  
 
```
run_dbcan {MAG_ID}/{MAG_ID}.faa protein --out_dir {MAG_ID}_dbcan --db_dir /bin/run_dbcan/db/
```
* Processed the run_dbcan outputs using `code/combine_dbcan.R` to combine dbcan annotations. 
	* Followed [Krüger et al. 2019](https://www.nature.com/articles/s41396-019-0476-y#Sec21) to filter annotations:
		* HMMER: E.Value<1e-20,Coverage>0.3
		* DIAMOND: E.Value<1e-20,X..Identical>30
	* Only kept the annotations that were consistent between the two tools
	* Saved outputs
		* Genes annotated as CAZymes `analysis/07_annotate_MAGs/summarized_outputs/cazymes_gene_assignments.txt`
		* Number of CAZymes assigned to different families, pee MAG
			* As a table: `results/tables/cazymes_summarized.txt`. This table is referenced in the text as **Table SXXX**
			* As a heatmap: `analysis/07_annotate_MAGs/summarized_outputs/cazy_heatmap.pdf`
		* Table showing median # of genes from diff. CAZy classes summarized by genus + median total # of CAZymes: `results/tables/median_cazy_by_genus.txt`. This table is referenced in the text as **Table SXXX**
   		* Same by family: `analysis/07_annotate_MAGs/summarized_outputs/cazy_class_percentage_by_bac_family.txt`
    	* Figure showing # of genes from diff. CAZy classes grouped by family: `results/figures/cazy_bac_genus.svg`
* Annotated iron metabolism genes using FeGenie
```
will add details after asking Arkadiy
```
* Annotated biosynthetic gene clusters using Emerald
```
will add details after asking Ellen
```
## 8. Loss of function in Rhizobiales
Software used:
* GTDB-Tk v1.5.0 (Chaumeil et al. 2020)
* IQ-TREE v2.1.2 (Nguyen et al. 2015)
* BLAST (Altschul et al. 1990)
* iTOL (Letunic & Bork 2019)
* R libraries: tidyverse, ape, phytools, stringr, RColorBrewer

### 8.1. Rhizobiales phylogenomic analysis
* Assembled dataset based on the selection from [Volpiano et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8026895/#FS1)
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
	
	* Added outgroup from rhodobacter
	```
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/Rhodobacter_amnigenus/latest_assembly_versions/GCF_009908265.2_ASM990826v2/GCF_009908265.2_ASM990826v2_genomic.fna.gz
	```
* Added Rhizobiales MAGs from our data
	* Used `list_rhizobiales_mags.R`
	* Saved the list as `analysis/09_rhizobiales_phylogeny/list_rhizobiales_mags.txt`
	* Copied Rhizobiales MAGs into the same folder
```
cd ../genomes
while read p; do cp ../../../"$p" . ; done < ../list_rhizobiales_mags.txt

cd ../annotations
while read p; do cp ../../../"$p" . ; done < ../list_rhizobiales_mags.txt
```
* Ran GTDB-Tk to identify and align marker genes
```
gtdbtk identify --genome_dir genomes --out_dir gtdbtk_identify --extension gz --cpus 10
gtdbtk align --identify_dir gtdbtk_identify --out_dir gtdbtk_align --cpus 20
gtdbtk classify --genome_dir genomes --align_dir gtdbtk_align --out_dir gtdbtk_classify -x gz --cpus 20

```

* Run IQTree
```
iqtree -s ../gtdbtk_align/gtdbtk.bac120.user_msa.fasta --seqtype AA -bb 50000 -pre rhizobiales -seed 12345 -m TEST -T 12
```

### 8.2. Screening for genes related to C1 metabolism and nitrogen fixation
* Screened the fasta files from NCBI to check if they have nitrogenase listed 
```
 grep "nitrogenase iron protein" ncbi_annotations/*faa > grep_nitrogenase_iron_protein.txt
```


* Searched all genomes using blast, used a NifH (ABZ89802.1) as a blast query
	* Used  `analysis/09_rhizobiales_phylogeny/Snakefile_tblastn` for searching nucleotide fastas
```
cd analysis/09_rhizobiales_phylogeny/
snakemake --cores 10 -n -s Snakefile_tblastn 

```
	* Used  `analysis/09_rhizobiales_phylogeny/Snakefile_blastp` for searching predicted protein fastas of *only ncbi genomes* 
	* Concatenated all into one file 
```
cat  blast_nifh/tmp_*_report* > blast_nifh/blast_nifh.txt
```

* Used `code/rhizobiales_nihf.R` to analyze the results of nitrogenase search:
	* methods of search (grep on the ncbi annotations vs tblastn vs blastp) are almost entirely consistent. 
	* The only inconsistency: in GCF_002879535.1 grep doesn't show NifH but tblastn&blastp do. The protein (WP_143973967.1) is labelled as nitrogenase reductase, partial in the ncbi annotations.
	* Decided to use tblastn search results for annotating the tree

* Searched for genes related to methane and methanol metabolism 
	* Used tblastn
	* Searched for 4 genes: PmoC (WP_016921575.1), MmoX (ABD13903.1), XxoF (VVC56072.1), MxaF (CAD91828.2)
```
snakemake --cores 10 -n -s Snakefile_tblastn_metahne 
cd ../../
```
	
* Summarized all searches
	* Used `code/rhizobiales_nihf.R`	
	* Used tblastn results for the five genes: NifH, PmoC, MmoX, XxoF, MxaF
	* Prepared annotation files for iTOL, saved them as:
		* bacterial family: `analysis/09_rhizobiales_phylogeny/iqtree/itol_fam.txt`
		* presence/absence of the genes of interest: `analysis/09_rhizobiales_phylogeny/iqtree/itol_{gene_name}.txt`
		* source (NCBI or our data): `analysis/09_rhizobiales_phylogeny/iqtree/itol_source.txt`
		* number of occurrences (for the MAGs only): `analysis/09_rhizobiales_phylogeny/iqtree/itol_occurrences.txt`
* Visualized the tree using iTOL. This figure is referenced in the text as **Fig. SXXX**

## 9. Searching for cobalamin-dependent genes in algal MAGs
* Made a list of all algal MAGs with >90% completeness
	* Used `get_algal_mag_list.R`
	* Saved the list as `analysis/11_algal_MAGs/good_algal_mags_list.txt`
* Copied them into the folder

```
cat good_algal_mags_list.txt | xargs -I {} cp ../05_MAGs/MAGs/euks/{}.fa.gz .
gzip -d *.gz
```
* ran Snakemake pipeline, which executes the tblastn search
```
cd analysis/11_algal_MAGs/
snakemake --cores 10 -n -s Snakefile_tblastn 
cd ../../
```

