# Annotating MAGs

## Aims
* Select MAGs in top 12 genera (by # of occurrences)
* Annotate them
* Make functional predictions

## Progress
1. Made a new folder `analysis/07_annotate_MAGs/`
2. Used `code/find_dominant_bacteria.R` to make a list of MAGs to annotate
* For 12 general (top genera by occurrence, together account for >50% of occurrences, see `bacterial_genera_frequency.tsv` for the list) selected all MAGs that satisfy strigent quality requrements (>95% complete)
* Added to the list three more MAGs that didn't meet the criteria but are still interesting:
    * private_T1916_metawrap_bin.6 and private_T1889_metawrap_bin.7 are two Lichenihabitans MAGs that often co-occurr with the dominant MAG from this genus, public_SRR14722130_metawrap_bin.2. They both have lower completenes: 93% and 69%
    * private_T1894_metawrap_bin.4 is a MAG from VCDI01 (=Lichenicoccus). It meets the quality standard but wasn't selected since it's below the occurrence threshold.
* Saved this list in `analysis/07_annotate_MAGs/mag_list.txt`
* saved table with file paths and locustags to use for annotation as `mag_table.txt`
* saved table with inof about mags (taxonomy and QC) as `mag_table_info.txt`

3. Annotated the MAGs with prokka

Setting up on CC
```
module load  bioperl/1.7.7
git clone https://github.com/tseemann/prokka.git $HOME/prokka
prokka/bin/prokka --setupdb
```
Setting up on debary - had to fix "solving environment" problem
```
conda update conda
conda config --set channel_priority strict

conda create -n prokka2  -c conda-forge -c bioconda prokka=1.13=3
conda activate prokka2
conda update perl-bioperl
conda install gxx_linux-64
conda install gcc_impl_linux-64
conda install perl-app-cpanminus
cpanm install Bio::SearchIO::hmmer3 -f
conda install -c conda-forge mamba
mamba install -c conda-forge -c bioconda snakemake
```

4. Trying out with one MAG
```
cd analysis/07_annotate_MAGs/
cp ../05_MAGs/MAGs/bacs/public_SRR14722130_metawrap_bin.2.fa.gz .
~/prokka/bin/prokka --compliant --centre UoN --outdir public_SRR14722130_metawrap_bin.2 --locustag LICHENIHABSRR14722130 --prefix public_SRR14722130_metawrap_bin.2 public_SRR14722130_metawrap_bin.2.fa
```
5. Annotated MAG by [KAAS](https://www.genome.jp/kaas-bin/kaas_main). Used .faa files and GENES data set for prokaryotes (hsa, dme, ath, sce, pfa, eco, sty, hin, pae, nme, hpy, rpr, mlo, bsu, sau, lla, spn, cac, mge, mtu, ctr, bbu, syn, aae, mja, afu, pho, ape), BBH assignment method

6. Did KEGG annotation with [kofamscan](https://github.com/takaram/kofam_scan)
```
cd /data/databases/KOfam
wget ftp://ftp.genome.jp/pub/db/kofam/*
gzip -d /data/databases/KOfam/ko_list.gz
tar -xvzf profiles.tar.gz
cd /data/tagirdzh/bin/
git clone "https://github.com/takaram/kofam_scan"

#put profile: /data/databases/KOfam/profiles/prokaryote in the config.yml

/data/tagirdzh/bin/kofam_scan/exec_annotation -o public_SRR14722130_metawrap_bin.2/public_SRR14722130_metawrap_bin.2.ko.txt public_SRR14722130_metawrap_bin.2/public_SRR14722130_metawrap_bin.2.faa -f  detail-tsv
/data/tagirdzh/bin/kofam_scan/exec_annotation -o public_SRR14722130_metawrap_bin.2/public_SRR14722130_metawrap_bin.2.kegg.mapper.txt public_SRR14722130_metawrap_bin.2/public_SRR14722130_metawrap_bin.2.faa -f mapper
```

7. Set up dbcan analysis
```
conda create -n run_dbcan python=3.8 diamond hmmer prodigal -c conda-forge -c bioconda
 conda activate run_dbcan
 
cd /data/tagirdzh/bin/
git clone https://github.com/linnabrown/run_dbcan.git
cd run_dbcan
test -d db || mkdir db
cd db \
&& wget http://bcb.unl.edu/dbCAN2/download/CAZyDB.09242021.fa && diamond makedb --in CAZyDB.09242021.fa -d CAZy \
&& wget https://bcb.unl.edu/dbCAN2/download/Databases/V10/dbCAN-HMMdb-V10.txt && mv dbCAN-HMMdb-V10.txt dbCAN.txt && hmmpress dbCAN.txt \
&& wget http://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
&& wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm && hmmpress tf-1.hmm \
&& wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm && hmmpress tf-2.hmm \
&& wget http://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm && hmmpress stp.hmm \
&& cd ../ && wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.fna \
&& wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.faa \
&& wget http://bcb.unl.edu/dbCAN2/download/Samples/EscheriaColiK12MG1655.gff
pip install dbcan==3.0.1
conda install snakemake

```
made a snakefile to annotate other mags

8. Used `code/combine_kegg.R` to combine kegg annotations for multiple MAGs from the same genus
    * saved files for KEGG mapper in `analysis/07_annotate_MAGs/summarized_outputs/*.kegg.combined.txt`
    * saved table combining all kegg annotations as `analysis/07_annotate_MAGs/summarized_outputs/kegg_combined.txt`
    * Made an NMDS for all KEGG annotations, saved as `results/figures/kegg_nmds.png`
9. Used `code/combine_dbcan.R` to combine dbcan annotations. Follwed [this paper](https://www.nature.com/articles/s41396-019-0476-y#Sec21) to filter annotations:
* HMMER: E.Value<1e-20,Coverage>0.3
* DIAMOND: E.Value<1e-20,X..Identical>30
* Only kept annotations cinsistent between the two tools
* saved output tables: `analysis/07_annotate_MAGs/summarized_outputs/cazymes_gene_assignments.txt` and `analysis/07_annotate_MAGs/summarized_outputs/cazymes_summarized.txt`
* saved heatmap as `analysis/07_annotate_MAGs/summarized_outputs/cazy_heatmap.pdf`
* summarized #of CAZymes per taxonomic unit:
    *  Table showing average # of genes from diff. CAZy classes summarized by geneus + average of total # of CAZymes: `analysis/07_annotate_MAGs/summarized_outputs/avg_cazy_by_genus.txt`
    * Same by family: `analysis/07_annotate_MAGs/summarized_outputs/cazy_class_percentage_by_bac_family.txt`
    * Figure showing # of genes from diff. CAZy classes grouped by family: `results/figures/cazy_bac_genus.png`

10. Testing: let's check if genes marked as missing by kofam_scan can be found by other means
* Took private_X1_metawrap_bin.1 (CAHJXG01), which, according to kofam_scan has all Calvyn-Bensen genes except rbc (4.1.1.39) and 3.1.3.37
* Annotated with KAAS the protein fasta produced by PROKKA, saved as `private_X1_metawrap_bin.1.kaas.txt`. Actually, it identified fewer genes
* Is rbc present in the unfiltered KO table? yes
```
grep "K01601" analysis/07_annotate_MAGs/private_X1_metawrap_bin.1/private_X1_metawrap_bin.1.ko.txt
ACETOBACX1_02165    K01601    376.07    18.1    0.00016    "ribulose-bisphosphate carboxylase large chain [EC:4.1.1.39]
```
In the filtered `private_X1_metawrap_bin.1.kegg.mapper.txt` it is listed as K00788, thiamine-phosphate pyrophosphorylase, e value of 1.1e-60. Most probably not what we are looking for


11. Looking for pseudogenes
Setting up pseudofinder
```
cd /data/tagirdzh/bin/
git clone https://github.com/filip-husnik/pseudofinder.git
cd pseudofinder
bash setup.sh
conda activate pseudofinder
```
Test run
```
mkdir pseudofinder
cd pseudofinder/
python3 /data/tagirdzh/bin/pseudofinder/pseudofinder.py annotate --genome ../private_T1916_metawrap_bin.6/private_T1916_metawrap_bin.6.gbk --outprefix private_T1916_metawrap_bin.6_ref --database /data/databases/NCBI_nr/nr --threads 16 --diamond --reference reference_genomes/GCF_020616535.1_ASM2061653v1_genomic.gbff
```
To run on all mags, added a new rule to snakefile

Turns out, some prokka-produced gbk files are broken (in LOCUS lines, seq ID and length are collided). Have to produce new fgbk files from gff. Used a script from github ([here](https://github.com/chapmanb/bcbb/blob/master/gff/Scripts/gff/gff_to_genbank.py))
```
conda deactivate
conda create -n  bcbio-gff -c conda-forge -c bioconda  bcbio-gff
conda activate bcbio-gff

conda create -n  emboss -c conda-forge -c bioconda  emboss
conda activate emboss
conda deactivate
conda activate pseudofinder
python3 /data/tagirdzh/bin/pseudofinder/pseudofinder.py annotate --genome public_SRR11456914_metawrap_bin.12.genbank --outprefix public_SRR11456914_metawrap_bin.12_pseudofinder --database /data/databases/NCBI_nr/nr --threads 16 --diamond
```

Visualize how parameter space affect % of identified pseudogenes:
```
cp pseudofinder/private_T1916_metawrap_bin.6_pseudofinder_* .
pseudofinder.py visualize -g private_T1916_metawrap_bin.6/private_T1916_metawrap_bin.6.gbk -log pseudofinder/private_T1916_metawrap_bin.6_pseudofinder_log.txt -op private_T1916_metawrap_bin.6_visualize_parameters
```

* Used `../../code/get_pseudogene_positions.R` to analyze:
    * Results mostly depend on the lenght cutoff. The default is 0.65 of "normal" length. Are the pseudogenes mostly on the ends of contigs? Seems not. checked the % of pseudogenes detected less than 10 bp away from a contig start/end; the highest was 8%, most were < 2%. In the same script, made a table to count pseudogenes and their % plus total # of ORFs
    * Made the list of KOs that are pseudogenized. Only selected those KOs that have ALL genes in a given MAG pseudogenized. Didn't see much of a pattern, but TonB is listed as a candidate pseudogene very often. The reason is the length: the ORF are either longer or shorter than its supposed orthologues. Not sure o which extent I can trust this.

* Used refgenome for Lichenihabitans <- this command didn't seem to work: hang up on codeml step for 24h
```
python3 /data/tagirdzh/bin/pseudofinder/pseudofinder.py annotate --genome ../private_T1916_metawrap_bin.6/private_T1916_metawrap_bin.6.gbk
```
Tried sleuth module. To get gff3 and ORF nucleotide fastas for the reference genome, used a [webtool](http://genome2d.molgenrug.nl/g2d_tools_conversions.html)
```
python3 /data/tagirdzh/bin/pseudofinder/pseudofinder.py sleuth -ctg ../private_T1916_metawrap_bin.6/private_T1916_metawrap_bin.6.fna -rg reference_genomes/GCF_020616535.1_ASM2061653v1_manual.gff3 -rn reference_genomes/GCF_020616535.1_ASM2061653v1_cds.fna -t 10 -out private_T1916_metawrap_bin.6_sleuth2
```
produced empty output

```
python3 /data/tagirdzh/bin/pseudofinder/pseudofinder.py sleuth -ctg ../private_T1916_metawrap_bin.6/private_T1916_metawrap_bin.6.fna -rg reference_genomes/GCF_004137085.1_ASM413708v1_genomic.gff -rn reference_genomes/GCF_004137085.1_ASM413708v1_cds_from_genomic.fna   -t 10 -out private_T1916_metawrap_bin.6_sleuth3
```



12. annotate iron genes
Setting up FeGenie
```
conda create -n fegenie -c conda-forge -c bioconda -c defaults fegenie=1.0 --yes
conda activate fegenie
```

Test run <- realized that the output only gives genes assigned to clusteres.
```
cd /data/tagirdzh/coverage/analysis/07_annotate_MAGs/
mkdir fegenie
cd fegenie
cp ../private_T1916_metawrap_bin.6/private_T1916_metawrap_bin.6.faa .
FeGenie.py -bin_dir . -bin_ext faa -t 16 --orfs --meta  --heme  --makeplots  -out private_T1916_metawrap_bin.6
```
Include all results:
```
FeGenie.py -bin_dir . -bin_ext faa -t 16 --orfs --meta  --heme  --makeplots  --all_results -out private_T1916_metawrap_bin.6_all
```

Run for all mags.  -all_results are questionable since it's prone to false positives (see arkadiy's comment on [github](https://github.com/Arkadiy-Garber/FeGenie/issues/27)). still included it
```
while read mag; do cp ../"$mag"/*.faa . ; done < ../mag_list.txt
cp /data/databases/NCBI_nr/nr .
FeGenie.py -bin_dir . -bin_ext faa -t 16 --orfs --meta  --heme  --makeplots -out fegenie_all -ref nr  --all_results

```

To solve the problem with blast
```
cd /data/tagirdzh/bin
git clone https://github.com/Arkadiy-Garber/FeGenie.git
cd FeGenie
bash setup.sh

/data/tagirdzh/bin/FeGenie/FeGenie.py -bin_dir . -bin_ext faa -t 16 --orfs --meta  --heme  --makeplots -out fegenie_all -ref nr  --all_results
/data/tagirdzh/bin/FeGenie/FeGenie.py -bin_dir . -bin_ext faa -t 16 --orfs --meta  --heme  --makeplots -out fegenie_strict -ref nr

```

Run one Nostoc MAG in strict mode, as faa (same, no genes for iron acquisition)
```
mkdir nostoc_test
 cd nostoc_test/
cp ../private_X5_metawrap_bin.1.faa .
/data/tagirdzh/bin/FeGenie/FeGenie.py -bin_dir . -bin_ext faa -t 16 --orfs --meta  --heme  --makeplots -out fegenie_private_X5_metawrap_bin.1

 cp ../../private_X5_metawrap_bin.1/private_X5_metawrap_bin.1.fna .
/data/tagirdzh/bin/FeGenie/FeGenie.py -bin_dir . -bin_ext fna -t 16 --meta  --heme  --makeplots -out fegenie_private_X5_metawrap_bin.1.fna

```

Asked Arkadiy to run the analysis. Visualized his results with `../../code/fegenie_viz.R` Saved the figure as `../../results/figures/fegenie_bubbleplot.png`

13. Manually searched for rhodopsin and associated genes
```
for file in ../*/*.faa; do blastp -query rhodopsin.fasta -subject "$file" -outfmt 6 -evalue 1e-5 ; done >> blastp_rhodopsin.txt
```
Manually extracted and searched against NCBI. All blast as rphodpsins


```
for file in ../*/*.faa; do blastp -query carotene-dioxygenase.fasta -subject "$file" -outfmt 6 -evalue 1e-5 ; done >> blastp_blh.txt
```

Manually extracted and searched against NCBI. All blast as Brp/Blh Beta-carotene 15,15'-dioxygenase.

Saved the results as `analysis/07_annotate_MAGs/summarized_outputs/manual_blast_search_results.txt`

14. antismash
```
 conda create -n antismash -c conda-forge -c bioconda antismash
  conda activate antismash
download-antismash-databases
antismash ../private_T1916_metawrap_bin.6/private_T1916_metawrap_bin.6.gbk --taxon bacteria --clusterhmmer  --cb-general --cb-subclusters   --cb-knownclusters --rre  --asf  --output-dir private_T1916_metawrap_bin.6

```
Added a snakemake rule. For the broken gbk files, annotated instead nucleotide fastas using the web server (see folder `analysis/07_annotate_MAGs/antismash/reran_fails`)

15. Used `code/kegg_heatmap.R` to make a heatmap of keggs of interest + results of blast search and fegenie results from Arkadiy

saved it as `analysis/07_annotate_MAGs/summarized_outputs/kegg_heatmap.pdf`

16. Arkadiy used SprayNPray on Lichenihabitans MAGs. In his table, the last column shows top blast hits for each ORF on this contig. While the majority of hits are to Lichinihabitans and other Rhizobiales, there are some hits to other groups common in lichens (Sphingomonas, Acetobacteraceae and Acidobacteria). Judging from the fact that these assignments are often present in the contigs where the majority of hits are to Rhizobiales, I don't think it's a binning error. However, this isn't evidence of HGT either. Arkadiy says: "BLAST hits from different taxa are caused by the fact that this MAG does not have a closely sequenced representative in NBCI's NR database (which I used as a reference for this SprayNPray run). The amino acid identity (AAI) values seem to average around 70-75%, which is pretty low. This makes the MAG's ORFs more likely to randomly recruit top hits from multiple taxa, instead of a single closely related taxa."

Asked him to run it on eukaryotic MAGs from VT1 (Platismatia) metagenome: private_X3_concoct_bin.10.fa (mycobiont), public_SRR7232211_concoct_bin.8.fa (alga), and private_T1904_concoct_bin.38.fa (another alga), to see if there are HGT from bacteria.

17. Ellen and Paul will compare the protein space of our bacterial mags to their collection. Sent them the PROKKA command I used, so they will annotate all bacterial MAGs, not just the ones that I selected for annotation

18. Ellen used emerald (Rob's lab secondary metabolism tool). I used `code/emeralg_heatmap.R` to analyze the results and make heatmap and a table

19. Manually searched NifH genes in metagenomic assemblies
* [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8026895/) says in Rhizobiales nif genes can be transmitted via plasmids
* used blast to search metagenomic assemblies
* Query: NifH of Rhodoblastus (Beijerinckiaceae): SNB82122.1
```
tblastn -subject ../../05_MAGs/assemblies_paul/VT1_scaffolds.fasta -query nifh.faa -evalue 1e-10 -outfmt 6
SNB82122.1      NODE_833_length_12008_cov_16.419058     45.106  235     113     6       2       230     8022    8696    8.89e-19        89.7 <- blasts to protochlorophyllide reductase
SNB82122.1      NODE_294_length_35245_cov_234.966382    34.553  246     145     7       28      265     17629   16916   1.62e-16        82.8 <- blasts to protochlorophyllide reductase
SNB82122.1      NODE_733_length_14060_cov_16.796787     47.977  173     86      3       8       177     4417    4932    7.14e-12        68.9 <- blasts to protochlorophyllide reductase

tblastn -subject ../../05_MAGs/assemblies_paul/TS1974_scaffolds.fasta  -query nifh.faa -evalue 1e-10 -outfmt 6
SNB82122.1      NODE_59_length_279649_cov_99.727458     41.328  271     144     7       4       267     82781   81993   9.13e-22        98.2 <- blasts to protochlorophyllide reductase
SNB82122.1      NODE_1369_length_24129_cov_5.960912     41.176  272     138     7       4       267     22967   22194   2.24e-18        88.2 <- blasts to protochlorophyllide reductase
SNB82122.1      NODE_3254_length_6203_cov_6.059369      40.221  271     147     7       4       267     790     1578    3.30e-18        87.4
SNB82122.1      NODE_14367_length_1455_cov_2.450000     44.206  233     114     6       4       230     564     1232    1.34e-17        85.1
SNB82122.1      NODE_442_length_84115_cov_695.967880    34.553  246     145     7       28      265     74999   75712   1.21e-16        82.8
SNB82122.1      NODE_67642_length_483_cov_0.887850      44.375  160     80      5       28      183     465     1       4.32e-14        70.9
SNB82122.1      NODE_596_length_64693_cov_98.914060     38.224  259     147     7       8       261     11000   11752   6.80e-13        71.6
SNB82122.1      NODE_982_length_39011_cov_6.125552      36.996  273     158     7       1       265     25428   26228   1.07e-11        68.2
```

20. Made kegg annotation and clustering for all MAGs
* Ellen annotated genomes with prokka. I used files she send along with emerald results
* made a different folder for the new annotations `mkdir analysis/07_annotate_MAGs/annotation_ALL_mags`
* made a new snakefile to run kegg annotations
* used `code/bacteria_functional_clustering.R` for analysis + `code/kegg_module_reconstruct.R` as a helper script to reconstruct kegg modules
* saved tables with kegg families and modules occurrecnes as `analysis/07_annotate_MAGs/summarized_outputs/all_mags_kegg_combined.txt` and `kegg_module_matrix.txt`. Module complteness is only for mags that are >90% complete
* nMDS:
	* nmds for all MAGs that occur > 2 times based on kegg families: `analysis/07_annotate_MAGs/summarized_outputs/nmds_all_mags_kegg.png`
	* tried on full set of MAGs, failed to converge
	* nmds based on kegg modules (>2 occurrence, >90% completeness): `analysis/07_annotate_MAGs/summarized_outputs/nmds_all_mags_kegg_module.png`
* hierarchical clusting
	* compared different ways of clustering, for kegg families and kegg modules. For both, the three methods that made most sense in terms of # of clusters were kmeans, hdbscan and apcluster: `analysis/07_annotate_MAGs/genomes_by_kegg_*_clustering_comparison*` shows 
	* tables with results for 3 methods of clustering, based on families and based on modules: `analysis/07_annotate_MAGs/summarized_outputs/kegg_families_clustering.txt` and `analysis/07_annotate_MAGs/summarized_outputs/kegg_modules_clustering.txt`
	* heatmaps showing the 3 clustering methods Trees on the side are phylogenomic trees: analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_*_clustering_results.pdf
	* taxonomic coherence (how well the clustering maps to taxonomy): `analysis/07_annotate_MAGs/summarized_outputs/genomes_by_kegg_*_clustering_taxonomic_coherence.pdf`

21. Made a table for frequency of several key modules in key groups
	* `code/modules_of_interest_freq.R`
	* saved it as `analysis/07_annotate_MAGs/summarized_outputs/summarized_outputs/modules_of_interest_freq.txt`
	* will use for making a figure with presence/absence of pathways
	
22. Made a figure showing the presence of key pathways
	* `code/pathway_fig.R`
	* used `analysis/07_annotate_MAGs/kegg_of_interest_fig.txt`
	* save result as `results/figures/pathway_fig.svg`
	
## Results
### Carboon fixation
* Most analyzed had some parts of Calvin-Bansen pathway, usually one or two missing genes from different parts of the pathway:
    * often, it's phosphoribulokinase - the same gene reposrted as missing by Pankratov
    * ribulose-bisphosphate carboxylase missing from all MAGs of the new clade of Acetobacteriaceae.
    * Lichenihabitans, CAIMSN01 (=Lichenicola), VCDI01 (=Lichenicoccus), Sphingomonas, most LMUY01 (putatively Acidisphaera, Acetobacteriaceae)  lack both  ribulose-bisphosphate carboxylase  and phosphoribulokinase
    * RH-AL1 lacks sedoheptulose 1,7-bisphosphatase and ribulose-bisphosphate carboxylase, and in all except one MAG, phosphoribulokinase
    * CAHJWO01 (Chthoniobacterales) lack phosphoribulokinase, ribulose-bisphosphate carboxylase, and sedoheptulose 1,7-bisphosphatase
    
* There were some Aceto MAGs with complete Calvin-Bansen pathway:  private_VT34_metawrap_bin.3  (CAHJXG01), public_SRR13125477_metawrap_bin.5 (LMUY01)
* Acidobacteriaceae (CAHJWL01, EB88 and Terriglobus) had 5 missing genes
* Nostoc is also shown to lack sedoheptulose 1,7-bisphosphatase, but it's function is performed by fructose-1,6-bisphosphatase II, which is present according to prokka

* Tried to search for rubisco manually
    * Used BAJ80171.1 and GAN73793.1 (large and small chain rubisco from Acidiphilium) and searched it against LMUY01, VCDI01, and Acetobacteraceae gen.sp.
    * In most MAGs didn't find anything. The only exception is a candidate for large subunit in VCDI01 (VCDI01T1894_02700, Ribulose bisphosphate carboxylase-like protein 2 according to PROKKA, K24268=3-oxoisoapionate-4-phosphate decarboxylase). Nothing similar to the small subunit
    * Used QQM06217.1 (rubisco from Rhodopseudomonas) to search agains Lichinhabitans and RH-AL1 MAGs. Nothing
    * tblastn didn't yield anythin new

* NiFe Hydrogenases:
    * present in nostoc
    * CAHJXG01 have K00441 frhB; coenzyme F420 hydrogenase subunit beta [EC:1.12.98.1]
    * absent in others
    * Aceto and Beijerinckos have NADP-reducing hydrogenases according to PROKKA - but these blast as NADH-quinone oxidoreductases or formate dehydrogenases, so likely are not hydrogenases
    
* Alternative pathways that are missing (list from [here](https://www.biorxiv.org/content/10.1101/2021.04.29.441244v1.full.pdf):
    * reductive tricarboxylic acid (rTCA) cycle,
    * 3-hydroxypropionate bi-cycle (HBC),
    * Wood-Ljungdahl (WL) cycle
    * dicarboxylate/4-hydroxybutyrate (DH) cycle
    * 4-hydroxypropionate cycle

### bacteriochlorophyll
* Many Aceto and Beijerincko had (almost) complete pathways for bacteriochloropphyl.
    * All Lichenihabitans (Beijerincko) and almost all CAHJXG01 (Aceto) and RH-AL1 (Beijerincko) lack divinyl chlorophyllide a 8-vinyl-reductase (BchJ/1.3.1.75/1.3.7.13 = K04036/K19073/K21231). It is present though in private_T1887_metawrap_bin.15 (CAHJXG01) and private_T1894_metawrap_bin.9 (RH-AL1)
    * The new clade of Acetobacteriaceae mostly had all genes requeired. One MAG (private_X1_metawrap_bin.1) didn't have any, despite being 99% complete and having an almost complete Calvin-Bansen pathway!
    * VCDI01 (=Lichenicoccus) has complete pathway
    * CAIMSN01 (=Lichenicola) lacks the pathway completely
    * all LMUY01 lack it completely, except public_SRR13125477_metawrap_bin.5  which has all except divinyl chlorophyllide a 8-vinyl-reductase. The same MAG have a compete Calvin-Bansen pathway
* Acidobacteriaceae, Sphingomonas  and CAHJWO01 (Chthoniobacterales) lack it

* Tried to search for DVR manually
    * Used VVC55152.1 (protein sequence of divinyl chlorophyllide a 8-vinyl-reductase from RH-AL1) as a blastp query for searching Beijerincko MAGs. In RH-AL1 didn't get any hits with evalue<1e-5, except in the RH-AL1 MAG that had DVR (there, the hit has an evalue about 1e-170). In Lichenihabitans, got hits with evalue about 1e-7 (LICHENIHSRR11456919n1_01627, LICHENIHSRR11456919_03125, LICHENIHABT1916_04068). The genes retrned as hits didn't have a KEGG annotation, and were classified by PROKKA as either 2-alkyl-3-oxoalkanoate reductase or 3 beta-hydroxysteroid dehydrogenase/Delta 5-->4-isomerase
    * Used SIR59060.1 (DVR from an Acetobacteracaea Acidiphilium) to search CAHJXG01. The MAG that had DVR according to KEGG yielded one hit evalue 1e-12. The other MAGs returned either nothing, or a hit to a protein identified as Aurachin B dehydrogenase by PROKKA (CAHJXG01VT34_02907 and CAHJXG01X15_03481, evalues 1e-7 and 1e-10). Unlike the search on Beijerinckos, the difference between the evalues of the "true" and the "false" hits wasn't that big. If searched against NCBI the "true" hit has a hit to other DVR genes, but the others (aka Aurachin B dehydrogenase) blasts as a NADH dehydrogenase.
    *  Same result, regardless of whether I searched protein fastas produced by PROKKA or the whole MAGs.

### Nitrogen fixation
* According to KEGG, Only in Nostoc
* [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8026895/) says in Rhizobiales nif genes can be transmitted via plasmids
    * used blast to search metagenomic assemblies
    


### Anoxygenic photosynthesis in purple bacteria (module M00612)
* This module consists of two sub-modules: Anoxygenic photosystem II module (M00597) and Calvin cycle (module M00165).
* For Calvin cycle, see above. It’s actually easier to view it here, under the modules tab though.
* All Acetos except CAIMSN01 (=Lichenicola) plus both Beijerincko (Lichenihabitans and RH-AL1) have both genes from Anoxygenic photosystem II. Nobody else does
* Anoxygenic photosystem II RC-LH1 is repsresented well, and Anoxygenic photosystem II RC-LH2 is missing. Not surprising - according to Yurkov LH2 is present sporadiclly in AAPs

### Quorum sensing
* I expected more to be present, actually
* Some Acetos have parts of DSF pathway: all CAIMSN01 MAGs and 2 MAGs from CAHJXG01 have RpfF (DSF synthase) and RpfC (sensor histidine kinase), the first two steps of the pathway. The intracellular part is hit-and-miss (mostly miss)

### Transporters and secretion systems
* The secretion systems are not very consistent within a genus
* Type 2 is usually mostly present (though gspS and gspO often are missing)
* Type 4 and Type 6 are present sporadically: e.g. among 6 CAIMSN01 MAGs, two had Type 4 and a third had an almost-complete Type 6
* ABC Transporters are also fairly variable within a genus
* Polyol transporters are usually present
* Also usually present: ribose, heme, phospholipids, phosphates, lipopolysacharides, alkanesulfonate
* Nickel transporters NikABCDE are missing from all
* Iron-siderophore transporter FepBCDG are missing
* Ferric hydroxamate transporter fhuBCD are missing
* Vitamin B12 transporter btuCDF are missing
* For efeUOB iron uptake system see https://onlinelibrary.wiley.com/doi/10.1111/j.1365-2958.2007.05802.x
* ethylmalonyl-CoA (EMC) pathway isn't complete in any mag
* Acetos, Sphingos, and Nostoc had generic fur regulator of iron uptake, and Beijerinkos had Rhizobiales-specific irr. See [this paper](https://journals.asm.org/doi/full/10.1128/mbio.02900-21) for info. Acidos seem to have none?

### Two-component systems: list here more-or less complete ones
* In  VCDI01 (=Lichenicoccus): K+ limitation (Kdp proteins), Acidic conditions (ChvG and ChvI), low Nitrogen availability (proteins Gln and Ntr), chemotaxis
* In CAHJXG01, Aceto gen.sp, LMUY01: same plus phosphate limitation (Pho proteins)
* Lichenihabitans:  Acidic conditions (ChvG and ChvI), low Nitrogen availability (proteins Gln and Ntr), chemotaxis, Cell cycle family, Redox signal
*  RH-AL1:  K+ limitation, Acidic conditions,  low Nitrogen availability, chemotaxis, Cell cycle family, Redox signal
* CAHJWO01 (Chthoniobacterales):  weirdly, no complete system,  Acidic conditions is missing completely
* Acidos: K+ limitation
* Sphingos: K+ limitation, Acidic conditions,  low Nitrogen availability, Redox signal

* I think it’s interesting that both Bereinkos have low Nitrogen availability system mostly present, but both lack NifA (Nif-specific regulatory protein)

### Antibiotic resistance
* Acidos, Sphingomonas, and 1 MAG from Beijerinckos (RH-AL1) have a complete beta-Lactam resistance module (M00627)
* CAHJXG01 have Multidrug resistance, efflux pump BpeEF-OprC

### CAZymes
* Acetos and Beijerinckos have fewer CAZymes compared to other families
    * Acetos: min=47, max=84, avg=62
    * Beijerinckos: min=26, max=84, avg=61
* Acidos and Nostoc have a lot, but while Nostoc is very rich in GTs, Acidos are rich in GHs
    * Acidos: min=87, max=149, avg=115,
    * Nostoc: min=117, max=140, avg=129
* Acidos are the only group that have GH as their dominant class, others tend to have more GTs
* There's variation between diff. MAGs, but mostly MAGs from one family have simlar profiles
* Across all groups, most GTs come from GT2 and GT4
* Acidos have some mannose-acting GHs that are mostly lacking in other groups:
    * GH92 (mannosidase): present in almost all Acidos (16 out of 17 MAGs), 3.5 genes/MAG on average. Outside of Acidos present in 3 Sphingomonas MAGs and 3 CAIMSN01 (Acetos)
    * GH125 (exo-alpha-1,6-mannosidase): present in all Acidos, always one gene, outside of Acidos is present only in 1 CAHJWO01 (Chtoniobacteriales)
    * GH38 (alpha-mannosidase): present in all Acidos, outside of Acidos is present in Nostocs and CAHJWO01 (Chtoniobacteriales)
    * GH76 (alpha-1,6-mannosidase/alpha-glucosidase): present in almost all Acidos (16 out of 17 MAGs), outside of Acidos is present in one CAHJWO01 (Chtoniobacteriales)
* Also overrepsented in Acidos: GH51 (endoglucanase, endoxylanase, cellobiohydrolase), GH1, GH2 (multifunctional), GH27, GH144 (endo-beta-1,2-glucanase), GH146 (beta-L-arabinofuranosidase),  GH18 (chitinase, also present in Nostoc and Sphingomonas), GH28 (polygalaturonases), GH55 (beta-1,3-glucanase), GH43_26 (L-arabinofuranosidase),
* One family present in all Beijerinckos and Sphingos, and all but one Acetos - and missing from Acidos and Nostocs is GH103 (peptidoglycan lytic transglycosylase family 3)
* One family that is present in almost all MAGs and isn't present at all in Acidos is GH13_9 (a-1,4-glucan branching enzyme)
* In general, GH5s are not well represented in any group

### KEGG overview
* On a large scale, in line with taxonnomy: Alphaproteobacteria are all together; Acidobcateria, Verrucomicrobia, and Cyanobacteria form well-separated clouds
* Still, interesting that within Alphaproteobacteria Acetos and Beijerinckos are mixed, and Sphingos are a bit separate, which is not in accordance with taxonomy (Beirinckos and Sphingos are more closely related to each other than to Acetos). Evidence that Acetos and Beijerinckos form a single functional "guild"?


### Steroids
* No MAG has any genes from steroid metabolism according to KEGG. Searched for osc (K01852, mentioned in [this paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4919349/)) - didn't find

### Secondary metabolism
* antiSMASH
* Found some homoserine lactone clusters: these are known to be involved in quorum sensing (see [this](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3425365/) paper)
    * searched for all Homoserine lactone receptors as listed by kegg found:
    * K18098 LuxR family transcriptional regulator, quorum-sensing system regulator BjaR1
    * K18096 acyl-homoserine lactone synthase - produces same lactone. both are present in one RH-AL1 and one Lichenihab MAG + LuxR alone is in one CAIMSN01 MAG
    * K20268 rhizosphere induced protein - is part of the same system in Rhizobium. Present in the same Lichenihab MAG
* Carotenoid GC in  CAHJXG01 (private_VT34_metawrap_bin.3, private_X15_metawrap_bin.2, private_T1887_metawrap_bin.15, all 100% similarity), in LMUY01 (public_SRR13125477_metawrap_bin.5 100% similarity), in Sphingo (public_ERR4179389_metawrap_bin.2, 25% similarity, but this includes key genes), in RH-AL1 (private_T1894_metawrap_bin.9, private_TS1974_metawrap_bin.3, private_VT26_metawrap_bin.1)
* zeaxanthin GC in Sphingo (private_VT34_metawrap_bin.4, private_X3_metawrap_bin.3, public_SRR11456918_metawrap_bin.7, public_SRR11456919_metawrap_bin.1), Lichenihabitans (public_SRR11456919_metawrap_bin.3, private_T1916_metawrap_bin.6, public_SRR14722130_metawrap_bin.2,private_T1889_metawrap_bin.7),

* Emerald
* Nostoc has most of them, especially in terpens, RiPP and NRPs. Won't look too deeply into that and focus on other bacteria instead 
* Carotenoids (several BGCs among terpens)
* Antibiotics and toxins (some terpens, polyketides. most common is BGC0000228: granaticin)
* Hormons: 
	* BGC0001605 (gibberellin), 
	* BGC0002005 (RaxX, in plant pathogens mimcs host hormone, see [here](https://academic.oup.com/jxb/article/70/16/4267/5522262)), 
	* BGC0001463 animal hormon (chemical cue in marine animals, see [here](https://royalsocietypublishing.org/doi/full/10.1098/rspb.2013.3086))
	* BGC0001363 (induced chlamidospores in fungi, see [here](https://pubmed.ncbi.nlm.nih.gov/26943626/))
* Biosurfactants (NRP) [This paper](https://academic.oup.com/jac/article/57/4/609/669417) says: "some may play essential roles for the survival of biosurfactant-producing microorganisms through facilitating nutrient transport or microbe–host interactions or by acting as biocide agents. Biosurfactant roles include increasing the surface area and bioavailability of hydrophobic water-insoluble substrates, heavy metal binding, bacterial pathogenesis, quorum sensing and biofilm formation"
	* enabling bacterial endosymbionts to colonize fungi: (holrhizin). see doi:10.1039/c8ob01515e 
* Exopolysaccharides (in most Acetos and Acidos and Beijerinckos), and teichuronic acid (cell wall polisacchiride, in Acidos and Nostocs + few in others). Acidos have also capsular polysaccharide BGCs (not on heatmap)

### Pseudogenezation
* Nostoc has higher # of genes (6939 on avg) and higher 5 of pseudogenes (18%) compared to others (3000-5000 ORFs and 6-12% pseudogenes).
* Pseudogenes include non-ORFs (NOT a part of PRKKA annotation) and some PRFs (those that were either fragmented, or too short compared to known genes from these families). No genes were found with dNdS > 0.3
* TonB is listed as a candidate pseudogene a lot

### Iron metabolism
* siderophore gene clusters are absent according to antiSMASH
* fegenie finds potential siderophore clusters in all Nostocs and few others, but these might not be real
* siderophpre transporters are present in all, but they are similar to other transporters, can say for sure
* iron transporters are present in Acetos (except 3 MAGs, one CAIMS and two LMU), Beijerinckos, Sphingos, all but one Nostocs, and some others


### Cofactors
1. Riboflavin
* no MAG has 5-amino-6-(5-phospho-D-ribitylamino)uracil phosphatase. However [this paper](https://www.tandfonline.com/doi/full/10.1080/1040841X.2016.1192578) says that this step is performed by broad spectrum hydrolases

2. thiamine
* Nobody has thiamine biosynthesis pathway, but some have Thiamine salvage pathway (thiMDE). See [this paper](https://journals.asm.org/doi/10.1128/JB.00641-06#:~:text=Overall%2C%20the%20results%20show%20that,survival%20when%20thiamine%20is%20limiting.)


### Functional clustering
* nMDS didn't work well. Will focus on hierarchical clustering
* Analysis based on kegg modules make more sense that kegg families. In families, the correlation matrix was too "saturated" (most values > 0.95), and the resulting grouping was weird (e.g. most of non-proteobacteria ended up together in one cluster, including cyanos and acido!). Will focus on modules in the next points
* Used three methods of clustering (hdbscan, kmeans, apcluster). Each was roughly consistent with phylogeny
	* hdbscan was mostly like "one phylum - one cluster", with some exceptions (one Actonomyceta ended up together with Verrucomicrobia, and Acidobacteriota was split between two clusters). NB: cluster "0" contains unassigned MAGs!!!
	* apcluster and kmeans were more granular, but they were inconsistent among each other. E.g., in apcluster, the majority of Rhizobiales grouped together with a portion of Acetobacterales; in kmeans, many of the same MAGs were together with Caulobacterales
	* probably the smartest choice is to stick with taxonomy






### Literature
* good papers about marine AAPs: https://journals.asm.org/doi/full/10.1128/AEM.03678-15,
* a bunch of MAGs have a phage terminase K06909 (https://journals.asm.org/doi/10.1128/JVI.01328-19) can it be for horizontal gene transfere?


