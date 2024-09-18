# Package names
packages <- c("tidyverse",
              "rentrez",
              "readxl")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# load in libraries
library(tidyverse)
library(rentrez)
# also used 'readxl', but call it directly (only used once)



# supplementary material --------------------------------------------------

# pull data from full list of metagenomes in supplementary materials
supp_mater <- readxl::read_xlsx('data/media-1.xlsx', 
                                sheet = 'S1_metagenomeList',
                                skip =1) %>% 
  dplyr::select(-10) %>% # single SRR code in here
  rename_all(~tolower(.x) %>% 
               str_replace_all(., ' ', '_')) 

# Rentrez search ----------------------------------------------------------


# check databses on NCBI to confirm SRA is on there (it is!)
entrez_dbs()

# check searchable terms of SRA
entrez_db_searchable(db = 'sra')

# PDAT/MDAT are going to be what we want
# format searchable fields can be found here: https://www.ncbi.nlm.nih.gov/books/NBK49540/
# PDAT: 
# - The date that records were made public in Entrez. 
# - The date format is YYYY/MM/DD. 
# - The colon ( : ) separates the beginning and end of a date range.
# MDAT:
# - The date of most recent modification of a sequence record. 
# - The date format is YYYY/MM/DD.
# - Only the year is required. The Modification Date is often used as a range of dates. 
# - The colon ( : ) separates the beginning and end of a date range.

# Lichen associated search keys:
# - Lichen
# - Lichen metagenome
# - Lichens
# - Pezizomycotina

# Method associated search keys:
# - WGS (whole genome sequencing)
# - Illumina/illumina

# changed date from last search to overlap with previous search.
# a few records must have been processing, and not caught


sra_search <- entrez_search(db = 'sra',
                                           term = '(Lichen[ALL] OR lichen metagenome[ORGN] OR Lichens[ALL] OR Pezizomycotina[ORGN]) AND WGS[STRA] AND 2021/10/01:2023[PDAT] AND (Illumina[ALL] OR illumina[ALL])', 
                                           use_history = TRUE,
                                           retmax = 99999)


# check number of records
sra_search$count


# because of large number of records, need to run this over a loop
# note: loop is faster than map/lapply because you can do batches of records
# each batch is one 'request'

# set up objects to fill
new_summaries <- vector()
sra_records <- tibble()
new_output <- tibble()

# run for length of lichen records in batches of 50
for(seq_start in seq(1, sra_search$count, 50)){
  Sys.sleep(3) # NCBI limits to 1 request per 3 seconds, this is just in case
  new_summaries <- entrez_summary(db = 'sra',
                                  web_history = sra_search$web_history,
                                  retmax = 50, retstart = seq_start)
  new_output <- tibble(sra_uid = extract_from_esummary(new_summaries, elements = 'uid'),
                       expxml_text = extract_from_esummary(new_summaries, elements = 'expxml'),
                       lichen = str_extract(expxml_text, '(?<=ScientificName=\\").+(?=\\"/><Sample)'),
                       produced_by = str_extract(expxml_text, '(?<=center_name=\\").+(?=\\" contact_name)') %>% 
                         if_else(is.na(.), 
                                 str_extract(expxml_text, '(?<=contact_name=\\").+(?=\\" lab_name)'),
                                 .),
                       #date = extract_from_esummary(., elements = 'createdate'),
                       runs_text = extract_from_esummary(new_summaries, elements = 'runs'),
                       sra_id =str_extract(runs_text, pattern = '(?<=acc=\\").+(?=\\" total_spots)'),
                       bioproject = str_extract(expxml_text, '(?<=Bioproject>).+(?=</Bioproject)')) %>% 
    dplyr::select(-expxml_text, -runs_text)
  sra_records <- rbind(sra_records, new_output)
  cat(seq_start+49, "summaries downloaded\r")
}

# check it (should display nothing)
sra_records %>% filter(is.na(bioproject))


# Check SRA records against linked published articles  -----------------------

# Check for linked articles on pubmed
# NOTE: THIS RETURNS VARIABLE RESULTS. THIS IS A KNOWN API ISSUE WITH NCBI
# PLEASE USE THE .CSV FILE BELOW FOR THE OUTPUT
# IT WAS MANUALLY CURRATED FROM MULTIPLE LINK RESULTS
sra_pubmed_links <- entrez_link(dbfrom = 'sra',
                                db = 'pubmed',
                                web_history = sra_search$web_history)


publication_summaries <- entrez_summary(db = 'pubmed',
                                        id = sra_pubmed_links$links$sra_pubmed)

sra_pubmed_links$links$sra_pubmed %>% length()

tibble(pubmed_uid = map_chr(publication_summaries, ~.x$uid),
       lichen = NA_character_,
       field_collected = NA_character_,
       lichen_cleaning = NA_character_,
       authors = map_chr(publication_summaries, ~.x$authors$name %>% paste(collapse = ', ')),
       first_author = str_extract(authors, '^\\w+(?= )'),
       pubdate = map_chr(publication_summaries, ~.x$pubdate),
       year = str_extract(pubdate, '^\\d+(?= )'),
       title = map_chr(publication_summaries, ~.x$title)) %>% 
  dplyr::select(pubmed_uid, lichen, field_collected, first_author, title, lichen_cleaning, authors, pubdate)

# This is a data sheet for each publication, checking to see if:
# - the study actually looked at lichen metagenomes (lichen)
# - the lichen was collected from the field, and not cultured (field_collected)
publications <- read_csv('data/sra_linked_publications.csv')

# now to connect this back with our SRA records

# get SRA uids for the publications that meet our criteria
pubmed_sra_links <- publications %>% 
  filter(lichen == 'yes' & field_collected == 'yes') %>% 
  pull(pubmed_uid) %>% 
  entrez_link(dbfrom = 'pubmed',
              db = 'sra',
              id = .)


# Check SRA records are lichen --------------------------------------------

# Not all SRA records are lichen.
# create a .csv to check if the genus (if listed) could potentially be a lichen
sra_records %>% 
  mutate(publication_check = if_else(is.element(sra_uid,
                                                pubmed_sra_links$links$pubmed_sra),
                                     'confirmed',
                                     'unconfirmed')) %>% 
  filter(!str_detect(sra_id, paste(supp_mater$sra_id, collapse = '|')),
         !str_detect(lichen, paste(supp_mater$lichen %>% 
                                     str_extract('^\\w+(?= )') %>% 
                                     unique(), 
                                   collapse = '|'))) %>% 
  filter(publication_check == 'unconfirmed') %>% 
  pull(lichen) %>% str_extract('^\\w+(?= )') %>% unique() %>% sort() %>% 
  tibble(genus = .)
# add the below code to the above pipes to generate your own list to check
# write.csv('data/is_lichen_genus_check_pezizomycotina.csv', row.names = FALSE, na = '')



# lichen_check is the list of genera I've manually checked for being potentially a lichen
lichen_check <- read.csv('data/is_lichen_genus_check_pezizomycotina_complete.csv')

# check it. There should be a Geoglossum record present, it is NOT a lichen
lichen_check %>%  filter(str_detect(genus, 'Geogloss'))

str(lichen_check)



# Check SRA records associated with lichen metagenomes --------------------


# For our SRA records, we now have the following:
# - If there is a connected publication associated with the record
# - If that record informs us of how the lichen metagenome was generated
# - If the genus associated with the SRA record is for a lichen (and not a moth)

# Now, to generate a list of remaining records to manually check. 
# Need to confirm the lichen name (or mycobiont) and how it was collected
sra_records %>% 
  mutate(publication_check = if_else(is.element(sra_uid, pubmed_sra_links$links$pubmed_sra),
                                     'confirmed',
                                     'unconfirmed')) %>%
  filter(!str_detect(lichen, 
                     paste(lichen_check %>% 
                             filter(is_lichen == 'no') %>%
                             pull(genus),
                           collapse = '|')),
         !str_detect(sra_id, paste(supp_mater$sra_id, collapse = '|')),
         publication_check == 'unconfirmed',
         lichen != 'lichen metagenome') %>% 
  arrange(produced_by)
# Add the below line if to write out a .csv for manual checking
# write.csv('data/unconfirmed_lichen_metagenomes_pezizomycotina.csv', row.names = FALSE, na = '')

# we don't need to read it in yet. We'll only be using it in one step later



# Adding publications on cutoff date --------------------------------------

# Some records were submitted pre-October, 2021, but published post-October, 2021

# goign to write a short function to pull that information in
# manually checking bioprojects for information about reads
pull_records <- function(bioproject_id){
  entrez_search(db = 'bioproject', term = as.character(bioproject_id))$ids %>% 
    entrez_summary(db = 'bioproject', id = .) %>% 
    extract_from_esummary(., elements = 'uid') %>% 
    entrez_link(dbfrom = 'bioproject', db = 'sra', id = .) %>% 
    first() %>% first() %>% 
    entrez_summary(db = 'sra', id = .) %>% 
    tibble(sra_uid = extract_from_esummary(., elements = 'uid'),
           expxml_text = extract_from_esummary(., elements = 'expxml'),
           runs_text = extract_from_esummary(., elements = 'runs'),
           lichen = str_extract(expxml_text, '(?<=ScientificName=\\").+(?=\\"/><Sample)'),
           date = extract_from_esummary(., elements = 'createdate'),
           sra_id =str_extract(runs_text, pattern = '(?<=acc=\\").+(?=\\" total_spots)'),
           bioproject = str_extract(expxml_text, '(?<=Bioproject>).+(?=</Bioproject)'),
           produced_by = str_extract(expxml_text, '(?<=center_name=\\").+(?=\\" contact_name)') %>% 
             if_else(is.na(.), 
                     str_extract(expxml_text, '(?<=contact_name=\\").+(?=\\" lab_name)'),
                     .)) %>%
    dplyr::select(sra_uid, lichen, produced_by, date, sra_id, bioproject) %>% 
    unique()
}




# Publications missing
# Keuler et al. 2020 - Genome-scale data reveal the role of hybridization in lichen-forming fungi

pull_records('PRJNA576709')

# Vecherskii et al. 2022.- Metagenomes of lichens Solorina crocea and Peltigera canina
# NOTE: Housed in two differen bioprojects
pull_records('PRJNA756680')
pull_records('PRJNA756777')

# Vijayakumar et al. 2022 - Metagenomic analysis of lichen-associated bacterial community profiling in Rocella montagnei
# pull_records('PRJNA741221') # AMPLICONS



# Combining information ---------------------------------------------------

metagenome_check <- read_csv('data/unconfirmed_lichen_metagenomes_pezizomycotina_complete.csv')
metagenome_check

# What we already have
supp_mater

# Get a vector of known lichen genera (manually confirmed and those included already)
lichen_genera <- c(supp_mater %>% 
                     mutate(genus = str_extract(lichen, '^\\w+(?= )')) %>% 
                     pull(genus),
                   lichen_check %>% 
                     filter(is_lichen == 'yes') %>% 
                     pull(genus)) %>% unique()

# Get vector of sra id's confirmed to be metagenomes
metagenome_bad_sra <- metagenome_check %>% 
  filter(manual_check != 'confirmed') %>% 
  filter(!str_detect(sra_id, paste(supp_mater$sra_id, collapse = '|'))) %>% 
  pull(sra_id) %>% unique()


lichen_genera

# add in the records found in the previous step

sra_additional <- rbind(
  pull_records('PRJNA576709'),
  pull_records('PRJNA756680'),
  pull_records('PRJNA756777')
) %>% 
  dplyr::select(-date)


# Build full dataset

# get a full list
sra_new_pezizomycotina <- sra_records %>% 
  mutate(genus = str_extract(lichen, '^\\w+(?=\\s|$)')) %>%
  filter(str_detect(genus, paste(lichen_genera, collapse = '|')) &
           !str_detect(sra_id, paste(metagenome_bad_sra, collapse = '|'))) %>% 
  dplyr::select(-genus) %>% 
  rbind(sra_additional) %>% 
  filter(!str_detect(sra_id, paste(supp_mater$sra_id, collapse = '|'))) %>% 
  arrange(produced_by, sra_id) %>% 
  unique() %>% 
  left_join(., metagenome_check %>% 
              dplyr::select(sra_uid, publication_name) %>% 
              mutate(sra_uid = as.character(sra_uid))) 


# somoe of the records need to be filled in manually
sra_new_pezizomycotina %>% 
  filter(is.na(publication_name)) %>% 
  pull(bioproject) %>% 
  unique()

# Adding in the additional publications
# remove some problem genomes (Resl et al. 2022 included cultures)
sra_new_pezizomycotina <- sra_new_pezizomycotina %>% 
  rows_update(read.csv('data/publication_list.csv'),
              by = 'bioproject') %>% 
  filter(!str_detect(lichen, 'Sarea'),
         !str_detect(sra_id, 'SRR10277352|SRR18162812|SRR18162823|SRR18162821|SRR18162819|SRR18162818|SRR18162817|SRR18162815|SRR18162814'))



# List of SRA records to search -------------------------------------------

sra_new_pezizomycotina %>% 
  pull(sra_id) %>% 
  cat(file = 'output/sra_accession_for_sra_download.txt',
      sep = '\n')


write.csv(sra_new_pezizomycotina,
          'output/sra_records_pezizomycotina2.csv', 
          row.names = FALSE,
          na = '')
