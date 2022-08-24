# MAGs

## Goal: Identify functional roles of each MAG in each metagenome

### Process:
* Two sources of information: euk tree and coverage
1. Used  to collect all mag taxonomy annotation in one table called `analysis/05_MAGs/tables/MAG_taxonomy_combined.tsv`

2. Used `assigne_putaitve_mag_roles.R` to assign putative roles to all MAG/metagenome occurrences
	* Used the bwa tables: this way I can detect ALL MAG occurrences, including those MAGs that were removed during dereplication
	* The idea here that one MAG can play different roles in diff metagenomes, e.g. the same Nostoc can be a phtobiont in one lichen and live in cephalodia in another or a lecanoromycete MAG is the mycobiont in one lichen and a low-coverage contaminant in another. Therefore each MAG-metagenome pairing should be reviewed separately
	* To take into account that, I first make a table here with putative assignments as follows:
		* All chlorophyta are assumed "photobiont_chloro"
		* All Cyanobactera are assumed "photobiont_cyano", except
		* if in a metagenome that has a green alga, in that case Caynobacteria are assumed to be in cephalodia
		* All non-cyanobacteria bacteria are "bacteria_other"
		* If a fungal MAG is the only fungal MAG in a metagenome, it's assumed to be "mycobiont"
		* If there >1 fungal MAGs, the one with deepest coverage is a "mycobiont", and others are "fungi_other"
	* The table is saved as `analysis/05_MAGs/tables/MAG_putative_roles_bwa.tsv`
	
2. Added metadata to the trees:
	* Used `code/rename_bac_tree.R` to:
    	* rename the bacterial tree so each MAG has its taxonomic annotation in the name. new tree is in `analysis/05_MAGs/trees/renamed/gtdbtk.bac120.user_msa.fasta_renamed.treefile`
    	* Same script also makes upadted iTOL annotation file (`analysis/05_MAGs/trees/renamed/itol_gtdb-layer_renamed.txt`).
    	* Manually added first 4 lines from the original iTOL file (`analysis/05_MAGs/trees/bac_itol/itol_gtdb-layer.txt`) to the new

	* Used `code/rename_euk_tree.R` to add info on the source lichen to the tips. 
		* Looked at all both trees sent by Paul: old, without ref. genomes, and new.
		* To the tips added: the BAT taxonomy assignment from `code/combine_mag_taxonomy_annotation.R` and the name of the lichen the MAG is from
		* Wrote renamed tree as `analysis/05_MAGs/trees/euk_tree_with_ref/lichen_mags_iqtree_renamed.tree` and `analysis/05_MAGs/trees/eukaryotes_renamed.treefile`
		* The old tree (eukaryotes_renamed.treefile) is more consistent with published lecanoromycete phylogeny, will use it from now on
		
3. Manually checked the functional assignments and wrote them into `analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.tsv`	
	* Checked against the biology of the lichens
		* Renamed some Nostoc as from "cyano_cephalodia" -> "bacteria_other", in the cases where the lichen isn't known to form cephalodia and the cyanobacteria are likely a contaminant (public_SRR14722307_metawrap_bin.1 and public_SRR14722176_metawrap_bin.1 in Xylographa SRR14722188, public_SRR14722086_metawrap_bin.4 in Dimelaena SRR14722167)
		* Renamed some Nostoc "phtotbiont_cyano" -> "bacteria_other", in cases where a chlorolichen had a cyano MAG but its actual photobiont wasn't recovered (public_SRR14722307_metawrap_bin.1	SRR14722176, public_SRR14722230_metawrap_bin.2	GTX0158, public_SRR14722086_metawrap_bin.7	SRR14722098, private_X5_metawrap_bin.1	SRR14722043, private_X5_metawrap_bin.1	SRR14722044, private_X5_metawrap_bin.1	SRR14722054, private_X5_metawrap_bin.1	SRR14722141, private_X5_metawrap_bin.1	SRR14722144, public_SRR14722086_metawrap_bin.3	SRR14722098, public_SRR14722176_metawrap_bin.1	SRR14722176, public_SRR14722145_metawrap_bin.2	SRR14722054, 
		* Chnaged photobiont in ERR4179391 (P. polydactylon). Reclassified Coccomyxa (public_SRR14722135_metabat2_bin.5, cov 0.7) as algae_other, reclassified Nostoc (public_SRR11456915_metawrap_bin.3 cov 12, public_SRR11456921_metawrap_bin.3 cov 35) as a photobion_cyano
	* Checked against the renamed tree	
		* Changed mycobiont in SRR14722059 (tuckermanopsis) public_SRR14722038_metabat2_bin.1 -> public_SRR14722026_metabat2_bin.4, because it's grouped with other tuckermannopsis'es
		* Changed SRR14722135 (Mycocalicium) mycobioont assignemnet from public_SRR14722135_concoct_merged.1 to public_SRR14722135_metabat2_bin.8 (lower coverage but gorups with the right gorup)
		* Changed SRR14722098 (Acarospora sinopica) mycobiont from public_SRR14722098_metabat2_bin.9 to public_SRR14722090_metabat2_bin.14 (lower coverage but groups with the right group)
	* Checked against NCBI to confirm that the organism as listed in the NCBI metadata, FEN number, and placement in the tree are consistent. Made changes:
			* SRR14722289 changed Diploschistes (acc. to NCBI metadata)-> Parmotrema (acc. to the NCBI name, FEN number and the tree)
			* SRR14721950 changed Myalospora (acc. to NCBI name and metadata) -> Pseudosagedia (acc. to FEN number and the tree)
	* Identified problematic cases, that couldn't be resolved. All of them have "mycobiont_missassigned" in the `analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.tsv` table
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
	* used `code/make_list_excluded_mtg.R` to save this list as a table in `analysis/05_MAGs/tables/excluded_metagenomes.txt`

4. Renamed the euk tree again to reflect the latest changes
	* Used the second half of `code/rename_euk_tree.R`
	* Wrote the new tree into `analysis/05_MAGs/trees/eukaryotes_iterative_renaming.treefile`
	* In FigTree flagged the problematic cases and saved the image as `analysis/05_MAGs/trees/eukaryotes_iterative_renaming.treefile.pdf`

5. Save iTol bacterial tree 
	* Used Paul's original files from `analysis/05_MAGs/trees/bac_itol`
	* saved as `results/figures/bacterial_tree_itol.svg`
	* used `code/rename_bac_tree.R` to add barchart annotation of # of occurrences per MAG
    

### Prelim results
* Tree mostly makes sense on a local level, but the deep nodes are seemingly unresolved. Might want to run a more comprehensive phylogenomic tree 
* 14 metagenomes will be removed from co-occurrence analysis, because I suspect they were misidentified/mislabelled

## Goal: MAG recovery based on depth for diff classes of lichens
### Aims
* What coverage enough for recovery the mycobiont?
* What coverage enough for recovery the mycobiont + photobiont?
* Number of MAGs as a function of coverage
* Bacterial MAGs in macro vs crust
* Coverage of different types of MAGs. Is it common to have higher coverage bacterial cov than mycobiont cov?


### Progress:

1. Used `code/recovery_myco_photo_MAGs.R` to analyze recovery of mycobiont and photobiont MAG as a function of sequencing depth
	* Used the functional assignments from `analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.tsv`
	* But only included MAGs extracted by binning (including those pre dereplication). Didn't count a MAG as present if it didn't show up in binning but only in bwa alignments to avoid situations where a genome is present in a metagenome but could not be recovered as MAG due to e.g. low coverage
	* Saved figures as `analysis/05_MAGs/exploratory_fig/myco_photobiont_vs_depth_denseplots.png` and `analysis/05_MAGs/exploratory_fig/myco_photobiont_vs_depth.png`


2. Manually looked into all MAGs discarded during dereplication (`analysis/05_MAGs/tables/MAG_taxonomy_and_role_add_drep_corrected.txt`, see `code/recovery_myco_photo_MAGs.R`):
	* Makes sense: 
		* Tremella public_SRR7232212_concoct_bin.5.fa  grouped with another Tremella public_SRR7232214_concoct_bin.7.fa  
		* Four Letharia MAGs all grouped into public_SRR7232211_concoct_bin.3.fa
		* 8 photobiont MAGs grouped together with other algae
		* Protousnea MAG private_X7_concoct_bin.1.fa grouped with the other protousnea
		* Two Nephromas together (SRR14722113 and SRR14722175)
		* Two Dermatocarpons (SRR14721937 and SRR14721961)
		* Two Phaeophysceas together (public_SRR14721931_concoct_merged.0.fa and public_SRR14722054_metabat2_bin.5.fa)
		* Three Heterodermias together (private_GTX0163_concoct_bin.14.fa, public_SRR14721925_concoct_bin.8.fa, public_SRR14722138_metabat2_bin.7.fa)
		* Two other Heterodermias (public_SRR14722081_metabat2_bin.3.fa, public_SRR14722144_metabat2_bin.3.fa)
		* Two other Heterodermias (public_SRR13685161_metabat2_bin.1.fa, public_SRR14722017_metabat2_bin.2.fa)
		* Two Gomphillus (private_GTX0158_concoct_bin.31.fa, public_SRR13685154_metabat2_bin.1.fa)
		* Two Phyllopsoras (public_SRR13685157_metabat2_bin.1.fa, public_SRR14721960_metabat2_bin.2.fa)
		* Two Trapelias (public_SRR14721945_metabat2_bin.3.fa, public_SRR14722251_metabat2_bin.1.fa)
		* Two Umbilicarias (public_SRR14722016_metabat2_bin.3.fa, public_SRR14722058_metabat2_bin.4.fa)
		* Two Cladonias (public_SRR14722296_metabat2_bin.7.fa, public_SRR14722301_metabat2_bin.5.fa), (public_SRR14722057_concoct_merged.0.fa, public_SRR14722084_metabat2_bin.3.fa),(public_SRR14722170_metabat2_bin.4.fa, public_SRR14722292_metabat2_bin.5.fa), (public_SRR14722088_concoct_bin.14.fa, public_SRR14722171_concoct_bin.7.fa, public_SRR14722298_metabat2_bin.8.fa),(public_ERR4179389_metabat2_bin.32.fa, public_SRR14722277_metabat2_bin.3.fa, public_SRR14722319_metabat2_bin.1.fa), (public_SRR14721965_metabat2_bin.1.fa, public_SRR14722009_concoct_bin.5.fa, public_SRR14722087_metabat2_bin.6.fa), (public_SRR14722082_metabat2_bin.6.fa, public_SRR14722117_metabat2_bin.4.fa)
		* Two Micareposis (public_SRR14722275_concoct_merged.0.fa, public_SRR14722304_metabat2_bin.6.fa)
		* Peltigeras (public_SRR11456922_metabat2_merged.1.fa,public_SRR11456923_metabat2_merged.0.fa,public_SRR11456913_concoct_bin.13.fa,public_SRR11456915_metabat2_merged.0.fa,public_SRR11456921_concoct_merged.0.fa,public_SRR11456921_metabat2_merged.0.fa,public_SRR14722320_concoct_bin.17.fa)
		* Two Evernias (private_VT26_concoct_bin.9.fa, private_X9_concoct_bin.20.fa)
		* Two Tuckermannopsis (public_SRR14722026_metabat2_bin.4.fa, public_SRR14722291_concoct_bin.3.fa)
		* Two Imshaugias (private_TS1974_concoct_bin.10.fa,public_SRR14721985_metabat2_bin.4.fa)
		* Three Platismatias (private_VT1_metabat2_bin.10.fa, private_X3_concoct_bin.10.fa, public_SRR14722007_metabat2_bin.10.fa)
		* Two Hypotrachynas (public_SRR14722014_metabat2_bin.1.fa, public_SRR14722219_metabat2_bin.6.fa)
		* Parmotremas (public_SRR14722020_concoct_merged.0.fa, public_SRR14722021_metabat2_bin.1.fa), (public_SRR14721997_metabat2_bin.1.fa, public_SRR14722028_metabat2_bin.1.fa), (public_SRR14722077_concoct_merged.0.fa, public_SRR14722139_metabat2_bin.5.fa, public_SRR14722329_metabat2_bin.1.fa), (public_SRR14721943_metabat2_bin.1.fa, public_SRR14722068_metabat2_bin.3.fa, public_SRR14722108_metabat2_bin.1.fa)
		* Two Punctelias (private_VT16_concoct_bin.11.fa, public_SRR14722300_metabat2_bin.4.fa)
		* Xanthoparmelias (public_SRR12240174_metabat2_bin.2.fa,public_SRR12240175_concoct_merged.0.fa, public_SRR12240177_metabat2_bin.4.fa, public_SRR12240178_metabat2_bin.6.fa, public_SRR12240179_concoct_bin.1.fa, public_SRR12240180_concoct_bin.4.fa, public_SRR12240187_metabat2_bin.8.fa, public_SRR12240188_concoct_bin.1.fa, public_SRR13126796_metabat2_bin.1.fa, public_SRR13167197_metabat2_bin.1.fa, public_SRR14722127_concoct_bin.0.fa), (public_SRR13125985_concoct_bin.8.fa, public_SRR13126647_concoct_merged.0.fa, public_SRR13126828_metabat2_bin.1.fa)
		* Two Flavoparmalias (public_SRR14721992_concoct_merged.0.fa, public_SRR14722106_metabat2_bin.2.fa)
		* Lepatogium (public_SRR14722015_metabat2_bin.2.fa, public_SRR14722197_metabat2_bin.3.fa)
		* Phlyctis (public_SRR13685158_concoct_merged.0.fa, public_SRR14722311_concoct_merged.0.fa)
		* Pertusaria (public_SRR14721984_concoct_merged.0.fa, public_SRR14722056_metabat2_bin.4.fa)
		* Buellia (public_SRR14721949_concoct_bin.4.fa, public_SRR14722083_concoct_bin.22.fa, public_SRR14722105_concoct_merged.0.fa)
	* Weird ones: the same ones that were flagged in goal 1
		* The only fungal MAG from Parmotrema public_SRR14722131_metabat2_bin.1.fa grouped with the only Acanthothecis MAG public_SRR14721923_metabat2_bin.2.fa. Identity 0.9966565
		* The only MAGs from Platismatia and Lepra (public_SRR14722032_metabat2_bin.2.fa public_SRR14722186_metabat2_bin.3.fa), identity 0.986
		* Reassigned public_SRR14722188_concoct_bin.0.fa as "fungi other" since it already has a mycobiont assigned (public_SRR14722188_metabat2_bin.25.fa) <- need to check placement!
		* Ochrolechia and Harteliana (SRR14722159 SRR14722160). Reassigned public_SRR14722159_metabat2_bin.8.fa as fungi_other. Most probably a contaminant from the nighbouring well, since SRR14722159 already have a mycobiont MAG assigned
		* Two only fungal MAGs from Rinodina and Lecanora (public_SRR14722303_metabat2_bin.3.fa and public_SRR14722309_metabat2_bin.1.fa), 0.99 identical
		* Two only fungal MAGs from Leprocaulon and Lepraria (public_SRR14722274_metabat2_bin.2.fa, public_SRR14722324_concoct_bin.1.fa)
		* Two only fungal MAGs from Cladonia and Melanelia (public_SRR14722121_metabat2_bin.4.fa, public_SRR14722156_concoct_bin.48.fa)
		* Two only fungal MAGs from Lecanora and Pseudevernia (public_SRR14722033_concoct_bin.2.fa, public_SRR14722129_metabat2_bin.3.fa)
		* Two only fungal MAGs from Diploschistes and Parmotrema (public_SRR14722289_concoct_bin.0.fa, public_SRR14722326_metabat2_bin.1.fa)
		* Two only fungal MAGs from Pseudosagedia  and Megalospora (public_SRR14722085_metabat2_bin.6.fa, public_SRR14722221_metabat2_bin.2.fa)
		
			
3. Used `code/summarize_mag_cov_counts.R` to analyze number and total coverage per MAG category in each metagenome. 
	* Used `analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.tsv`, which is based on the bwa tables, to include all genomes present in a metagenome
	* Saved results into `analysis/05_MAGs/tables/MAG_coverage_summary.tsv` and `analysis/05_MAGs/tables/MAG_counts_summary.tsv`
	* Saved metagenomes that had a mycobiont MAG, but its coverage was lower than combined coverage of bacteria in `analysis/05_MAGs/tables/MAG_coverage_high_bacteria.tsv`
	* Identified metagenomes that had more than 1 mycobiont genome (not counting those already identifed as fungi_other in step 1), in `analysis/05_MAGs/tables/extra_myco_mtg.tsv`
    4. Used `code/completeness_breadth_vs_cov.R` to replicate Paul's result and check how coverage depth affects:
        * EukCC completeness for de novo produced MAGs
        * breadth of coverage for the bwa alignments
        * saved the result as `results/figures/completenes_and_breadth_vs_cov.png`
        * Used pauls' spreadsheet for the coverage and completeness stats, saved it as `analysis/05_MAGs/tables/euks_cov_QC.txt`

### Prelim results
* **Required depth:** Min depth for mycobiont MAG to be detected ~ 550 Mbp. Min depth for both myco and photo ~ 2 Gbp.
    * Weirdly, had one low-depth metagenome (SRR13685159, ~ 600 Mbp) that had algal MAG, but not myco. Pretty sure, this one geniunely didn't have myco MAG
* **Number of Mags as a function of coverage:** Doesn't seem to plato. Don't see any patterns based on archtecture
* Combined bacterial coverage was higher than the mycobiont coverage in 81 metagenome! These metagenomes included both high and low coverage, both crusts and macros
* A bunch of metagenomes included low-coverage genomes of othe mycobiont MAGs, e.g. private_X14_concoct_bin.2.fa (H. physodes mycobiont) present at >90% in X3 (Platismatia). Is it because Hypogymnia genome is too close to Platismatia and that's why it maps to it? Most probably not: in another Platismatia metagenome (VT1) it's only present at 0.4%. Most likely, the X3 sample had a little bit of Hypogymnia, and was sequenced deep enough to catch it.
* MAGs can evade detecton during binning because:
    * for a mag to reach 90% completeness, it need to have 10X coverage. For bwa aligments, you need only 4X to reach 90% breadth
    * presumably, mixing the two algal MAGs results in high contamination rates, similar to what we had in Alectoria. None of the 8 metagenomes where the two alga co-occurred yeilded a usable algal mag


## Goal: cooccurrence of everything
### Aims
* What bacteria occur together?
* Which MAGs occur most commonly?
* Any bacteria-alga stable cooccurrence
* Two-alga complex: what lichens they occur in?

### To do:
* Ask Paul how he made a cut-off for cooccurrence
* Apply the same procedure to all MAGs, include all metagenomes with more than one MAG
* Make an occurrence table for all MAGs vs all metagenomes

### Progress
1. Used `code/find_dominant_bacteria.R` to identify bacterial lineages of interest:
	* Put numbers of occurrence in `analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_*_frequency.tsv`. Here I used the bwa table to identify all metagenomes each MAG occurred in
	* Put numbers of MAGs in `analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_*_diversity.tsv`. Here is just a number of MAGs from a given group
        * Put fraction of of bacterial reads coming from a MAG/taxon. Only counted bacterial reads that are aligned to bacterial MAGs, i.e. summarized read counts for bacterial MAGs in all metagenomes. Saved results as  `analysis/05_MAGs/tables/bacteria_dominant_groups/bacterial_*_reads_percentage.tsv`.
        * Identified top-MAGs in each genus (by occurrence). Saved this list as  `analysis/05_MAGs/tables/bacteria_dominant_groups/top_mag_from_top_genera.tsv`
		* Made lists of top frequence genera and families by groups of lichens: photobionts, mycobionts, and compinations. Saved as 
			* saved them as `../05_MAGs/tables/bacteria_dominant_groups/bacterial_*_by_lichen_group.txt`
2. Used `code/cooc_graph_euk.R` to make a coocurrence map for mycobionts, photobionts and yeasts. Saved the image in `results/figures/coocc_graph_euk.png`
	* Used `analysis/05_MAGs/tables/MAG_confirmed_roles_bwa.tsv`
	* Filtered:
		* Removed metagenomes that were identified as missasigned
		* Removed metagenomes that didn't have a mycobiont mag
		* Only showed mycobionts, photobionts, cephalodia, and yeasts. Removed "fungi_other" to exclude cases where a lecanoromycete MAG is present a low coverage probably as a contaminant (those would result in erroneous co-occurrences)
3. Used `code/checking_coverage_coocurring.R` to compare coverage of the two coocurring algal MAGs. 
	* Saved the results into `analysis/05_MAGs/tables/two_alga_coverage_depth.tsv`
	* Checked that the metagenomes showing only one MAG of the two, genuinely don't have the other. Checked coverage breadths of the missing MAGs, and they are well below the threshold
	* Made a figure showing the relationship between breadhs of depths of coverage for the two mags. Left graph shows all metagenomes as dots, arranged based on the breadth of coverage of the two mags. The right graph shows depths as axes and only displays metagenomes that contained both mags, line is for the context, it shows y=x
4. Used `code/cooc_graph_licheninhab.R` to make a coocurrence map for mycobionts and licheninhbintans. 
	* Saved the image in `results/figures/results/figures/coocc_graph_licheninhab.png`
	* Explored coverages of the three most common Licheninhabitans MAGs: public_SRR14722130_metawrap_bin.2, private_T1916_metawrap_bin.6, private_T1889_metawrap_bin.7
5.  Used `code/cooc_graph_batceria_other.R` to make a coocurrence map for mycobionts and two sets of other bacteria: two genera of Acetobacteraceae (CAHJXG01 and LMUY01) and two genera of Acidobacteriaceae (EB88 and Terriglobus)
6. used `code/heatmaps_exploratory.R` to make first heatmaps:
	* Selected three bacterial families: Aceto, Acido, and Beijerinckiaceae (contians Lciehinhabitans). Made a heatmap with all metagenomes that yelded a mycobiont mag + all MAGs from these families + all cyano and chlorophotobiont. Saved as `analysis/05_MAGs/exploratory_fig/heatmap_full_families_of_interest`. Looks very granular, hard to pick patterns 
	* Selected five top bacterial genera:"Lichenihabitans","CAHJXG01","LMUY01","EB88","RH-AL1". Made a heatmap with all metagenomes that yelded a mycobiont mag + all MAGs from these genera + all cyano and chlorophotobiont. Saved as `analysis/05_MAGs/exploratory_fig/heatmap_selected_genera`. This is much better, will cross-reference with the mycobiont tree and search for patterns
	* Made a hetmap not on mag level but on genus level for bacteria.Selected 4 families: Aceto, Acido, Beijerinckiacea, and Sphoingomonadaceae. Saved as `analysis/05_MAGs/exploratory_fig/heatmap_bac_genus_level`. 
	* for the existing mag-level heatmaps will try to add phylogeny for the columns, or at least re-order the columns to manually add these phylogenyes in illustrator
7. used `code/heatmaps.R` for better heatmaps:
	* Added trees on the top to show relationship between the mags
	* Did it for the same families as above, added genus assignments as a bottom annotation
	* Streamlined the code and DRYed it
	* Checked the table `analysis/03_metagenome_reanalysis/all_metagenome_reanalysis.txt` to make sure that the order and class assignments are correct.
	* Added right annotation with the mycobiont orders. Double checked the assignments, the tree looks iffy: **asked david to run a proper tree**
	* Added another right annotation for the sequencing depth of the metagenome
	
8. Added analysis of bacterial abundance inferred from the depth of coverage. 
	* First tried to add this infp to the heatmap, failed (`code/heatmaps_abundance_failed.R`). 
	* Then used `code/relative_cov_bac.R` to plot it separately
	* Saved as `analysis/05_MAGs/exploratory_fig/relativa_bac_cov.png`
	* Made a box plot, by genus (only including the top 13 genera) + eukaryotes, saved as "results/figures/relative_cov_boxplot_genus.png"

9. Used `code/nmds_occur.R` to make nMDS for bacterial occurence (bactrial MAG/genus as 'species', lichen metagenomes as 'sites')
    * Most analyses I tried didn't converge or produced ordinations with 1-2 outliers skewing the whole plot
    * The only thing that worked more or less was to:
        * Removie low-depth metagenomes (<2Gbp)
        * Select 5 key families: Cyanos, Acetos, Acidos, Beijerinckos, Sphingos and UBA10450 (Chthoniobacterales)
        * Group occurences by bacterial genus
    * Still, lowish quality ordination (stress = 0.16, converged after 1435 tries)
    * Saved figures as:
        * ordintation showing sites (colored by lichen order) and bacterial genera:  `analysis/05_MAGs/exploratory_fig/nmds_sites_species.png`
        * ordintation showing sites colored by seq. depth: `analysis/05_MAGs/exploratory_fig/nmds_sites_depth.png`
10. Used `code/heatmap_parmeliaceae.R` to make heatmaps specifically for parmeliaceae metagenomes 
	* used metagenomes with > 2Gbp of data
	* made heatmaps for the 4 key families
11. Used `code/multivar_lichens_by_mag_occ.R` to try to group lichens and bacterial MAGs by their occurrence using multivar methods other than nMDS: hclust, cca
	* saved 2 cca plots. in each used occurrenc matrixes + predictors (photobion, mycobiont order, and depth)
		* `05_MAGs/exploratory_fig/cca_metagenomes_vs_bac_families_removed_outliers.pdf` is a CCA on all metagenomes with >2Gbp data, based on bacterial families, with removed outliers (Verrucariales and Lichinales)
		* `05_MAGs/exploratory_fig/cca_peltigerales_vs_bac_genera.pdf` is a CCA on Peltigerales metagenomes, based on bacterial genera
	* used clustering and heatmaps 
		* didn't make sense on the whole dataset (even patterns of co-occurrence I knew weren't showing up)
		* worked better on a subset: only used metagenomes with >2Gbp of data and bacterial genera that occurred in >4 metagenomes
		* saved exploratory clustering plots as `analysis/05_MAGs/exploratory_fig/metagenomes_by_mag_occurrences_clustering_comparison*.pdf` and ` analysis/05_MAGs/exploratory_fig/bacteria_by_mag_occurrences_clustering_comparison*.pdf`
		* used kmeans and apcluster, both for metagenomes and bacteria
		* saved heatmap as `05_MAGs/exploratory_fig/clustering_metagenomes_vs_bac_genera.pdf`
		* this heatmap shows clustering by kmeans and apcluster as a side annotation. taxonomy and photobiont are also side annotations. dendrograms are hclust!
		* tried to make a version with phylogenomic tree, but it wasn't interpretable

### Results
1. Identified lineages of interest. By occurrence, they are:
	* d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Beijerinckiaceae;g__Lichenihabitans;s__ (139 times)
	* f__Acetobacteraceae;g__CAHJXG01 (103 times)
	* f__Acetobacteraceae;g__LMUY01;s__ (79 times)
	* f__Acidobacteriaceae;g__EB88 (59 times)
	* d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Beijerinckiaceae;g__RH-AL1;s_ (59 times)
	* f__Acetobacteraceae;g__CAIMSN01. (56)
	* d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiae;o__Chthoniobacterales;f__UBA10450;g__CAHJWO01; (54)
	* f__Acidobacteriaceae;g__Terriglobus (45)
	* Nostoc (43)
	
	By number of MAGs: EB88 (33), CAHJXG01 (27), Sphingomonas (27), Tous-C9LFEB (21).
	
	**This mostly matches the bacterial lineages that had expansions on the tree (see below)**
2. Recreated Paul's result with two alga. Reasonably sure that this isn't an artifact. Both MAGs appear on their own occasionally (public_SRR7232211_concoct_bin.8 once, private_T1904_concoct_bin.38 8 times), but when they occur together (8 times), the coverage depth is very similar. All occurrences appear to be withi Lecanorales
3. Letharia spp. (probably) occurs with 4 diff algae, but cypho and tremella are the same
4. Licheninhabitans is everywhere, few MAGs but consistent occurence. 
	* Explored coverages of the three most common Licheninhabitans MAGs: public_SRR14722130_metawrap_bin.2, private_T1916_metawrap_bin.6, private_T1889_metawrap_bin.7
	* They often co-occur in diferent combinations, private_T1916_metawrap_bin.6 rarily occurs without either of the two
	* Relative coverage is all other the place
5. Acetobacteriacea also form a dense network, and Acidobacteriaceae make a much looser graph (see `results/figures/coocc_graph_acetobacteraceae.png` and `results/figures/coocc_graph_acidobacteriaceae.png`)
6. RH-AL1 (genus sister to Lichenihabitans) emerged as another genus of interest. Mostly occur in Lecanorales, sometimes together with Lichenihabitans, sometimes by itself. There are two closely related MAGs that ofet co-occur	
7. Enterovirga: another Beijerinckiaceae genus, found only in Peltigerales
8. Most MAGs in the dominant groups have lower abundance than the mycobiont, but there's a lot of variation. Only a few MAGs have hight abundance. One exception: Sphingomonadales genus UBA1936 had two MAGs and three ocurrences, each time had relative coverage around 1.
9. Comparing lists of dominant taxa by different metric (diversity = # of MAGs, # of occurrences, percentage of bacterial reads): more or less consistent of the order and family level (i.e. same taxa but in a different order). Genus, and especialy MAG levels are different. Concluded that reads percentage isn't a good metric as it's biased towards more deeply sequenced metagenomes
10. 12 bacterial genera account for 50% of bacterial reads. On the MAG level, 90 bacterial MAGs account for 50% bacterial reads. Top genera come from the same 4 families + one from Chthoniobacterales (Verrucomicrobia). Most top MAGs come from the same groups, but there are some others (e.g. Caulobacterales, Pyrinomonadales, Polyangiales, and more)
11. NB: our Lichenihabitans are selected according to GTDB. They can be theoretically split into:
    * Lichenihabitans sensu NCBI ("closest_placement_taxonomy" in GTDB_Tk output is Lichenihabitans psoromatis). The clade that contains all frequent MAGs
    * Lichenibacterium sensu NCBI ("closest_placement_taxonomy" in GTDB_Tk output is Lichenihabitans minor = Lichenibacterium minor in NCBI). The clade that contains all frequent MAGs
12. Looked up which top-genera that in the annotation have GTDB IDs have proper NCBI names:
    * CAIMSN01 = Lichenicola Taxonomy check:Inconclusive
    * VCDI0 = Lichenicoccus Taxonomy check:Inconclusive
    * LMUY01 = Acidisphaera Taxonomy check:Inconclusive
    * EB88 = Acidipila Taxonomy check:Inconclusive
    * HMF7410 = Mucilaginibacter Taxonomy check:IOK
    * UBA5172 = Granulicella Taxonomy check:Inconclusive
    * LB1R16 = Glacieibacterium Taxonomy check:Inconclusive

13. Looked up sources of MAGs, included MAGs from top families
    * CAHJ* MAGs from lichen metagenome (Lobaria)
    * CAJCIS01 MAG from Graz, couldnt find the description of sample, project is called "Microbial community analysis from the vegetation of Alpine bogs"
    * BOG-908 from soil
    * RH-AL1 from soil
    * Tous-C9LFEB from bog and freshwater
    * DTNP01 from smud sediment
    * KBS-83 couldn't find
    * CAIQPK01 lake water
    * UBA1936 wastewater
    * CADCVW01 soil biocrust
    
14. nMDS ordination isn't conclusive: some grouping is present but not much. Peltigerales are concntrated in one part of the cloud, and some bacteria are with them that looked associted with Peltigerals on heatmaps: Nostoc, Entervigra, Methybocaterium, Sphingomonas
    
    
15. Most frequent bacteria by lichen group (photobiont/mycobiont/combination) have clear differences:
	* Sphingomonadaceae are higher in cyanolichens (in Peltigerales to be exact). Lichenihabitans is much lower, but there are othe Rhizobiales (Methylobacterium and Enterovigra)
	* The only two Lichinomycetes have Ktedonobacteraceae (Chloroflexota) as the most common group! It's also present in a bunch of others, incl. parmeliaeceae and cladoniaceae
	* In arthoniomycetes second after Acetobacteraceae is Jatrophihabitantaceae (Actinobacteriota, Mycobacteriales). it's also present in peltigerales and lecanorales
	* In Umbilicariales, beijerinckos are rare, acido and aceto are at the top
	* Lecideales similar to general picture, but Lichinhabitans is rare (only one occurrence)
	* Lecanorales, Caliciales, Pertusariales, Baeomycetales are similar to the general picture. Teloschistales more or less too
	* Comparing chloro and cyano peltigerales shows same patterns: Lichenihabitans is most common in chloro, and far down in cyano. Sphingo and CAHJXG01 are in the top in both.
	* other lichen groups with mixed photobiont types (gyalectales, arthoniomycetes, eurotiomycetes) don't have minimal sample size to compare
	
16. Multivar methods of grouping metagenomes based on MAG occurrence: not sure
	* Hierarchical clustering worked only on a subset (removing shallow metagenome and rare bacteria). Tried it on 3 levels (MAGs, genera, families). maybe the problem is with uneven depth? may be individual metagenomes with singleton bacteria skew the picture?
	* In clustering of bacteria: Nostoc, Sphingomonas, Enterovigra and Methylobacterium are in one cluster according to kmeans
	* Clustering of metagenomes - hard to see a pattern. no clear link to photobiont
	* CCA with predictors (photobiont, mycobiont order, depth) made more sense. Tried it on the same three levels. 
	* Verrucariales and Lichinales were outliers on CCA, skewing the plot
	* In general, there are some differences between diff. mycobionts and photobionts. BUT not sure how to account for uneven sampling and outliers
	* In Peltigerales, can see depth as one main axis, and photobiont type as the second
	* in general, don't feel good about this analysis: the problem with uneven depth, biased sampling (lecanorales way overrepresented), rare bacteria potentially skewing the picture
