# adapted from Gulnara Tagridzhanova's script
# removed setting up the GTDB file and instead read in one already made
# changed directory structure to match my directory structure

# Note, not checking for installation, because this should have been done in
# the 'idtaxa_compiling_db.R' script.


# idTaxa - script that assigned taxonomic positions to bacterial rDNA (assemblies and reads) based on GTDB
library(tidyverse)
library(DECIPHER)
source('scripts/utils.R')
library(R.utils)

# reference DB to compare our reads to
trainingSet <- readRDS('data/db/gtdb_trainingset.rds')

#function to process fasta file
process_fasta<-function(file){
  dna<- readDNAStringSet(file)
  ids <- IdTaxa(dna, trainingSet, strand="both",threshold=50)
  plot(ids, trainingSet)
  assignment <- sapply(ids,
                       function(x)
                         paste(x$taxon,
                               collapse=";"))
  assignment <-data.frame(assignment)
  assignment$file<-file
  return(assignment)
}

#make list of files
# filenames <- list.files(path = "results/metaxa_bacteria/")
# filenames <- filenames[grepl("bacteria.fasta",filenames)]
# filenames_full<-paste0("results/metaxa_bacteria/",filenames)

# temporary for additional records
filenames <- list.files(path = "results/metaxa_bacteria")
filenames <- filenames[grepl("bacteria.fasta",filenames)]
filenames_full<-paste0("results/metaxa_bacteria/",filenames)

#remove empty files
nlines<-sapply(filenames_full,countLines)
filenames_nonempty<-filenames_full[nlines>0]

#reads
filenames_reads<-filenames_nonempty
lreads_1<-process_fasta(filenames_reads[1])
write.table(lreads_1,"results/idtaxa/idtaxa_reads_1.txt",sep="\t",quote = F, row.names = F)

lreads_2<-process_fasta(filenames_reads[2])
write.table(lreads_2,"results/idtaxa/idtaxa_reads_2.txt",sep="\t",quote = F, row.names = F)


lreads_3<-lapply(filenames_reads[3:6],process_fasta)
df_reads_3<-do.call(rbind,lreads_3)
write.table(df_reads_3,"results/idtaxa/idtaxa_reads_3.txt",sep="\t",quote = F, row.names = F)

lreads_4<-lapply(filenames_reads[7:15],process_fasta)
df_reads_4<-do.call(rbind,lreads_4)
write.table(df_reads_4,"results/idtaxa/idtaxa_reads_4.txt",sep="\t",quote = F, row.names = F)

lreads_5<-lapply(filenames_reads[16:100],process_fasta)
df_reads_5<-do.call(rbind,lreads_5)
write.table(df_reads_5,"results/idtaxa/idtaxa_reads_5.txt",sep="\t",quote = F, row.names = F)


lreads_6<-lapply(filenames_reads[101:length(filenames_reads)],process_fasta)
df_reads_6<-do.call(rbind,lreads_6)
write.table(df_reads_6,"results/idtaxa/idtaxa_reads_6.txt",sep="\t",quote = F, row.names = F)


#combine all together
df_reads_all<-rbind(lreads_1,lreads_2,df_reads_3,df_reads_4,df_reads_5,df_reads_6)
df_reads_all$metagenome<-str_replace(df_reads_all$file,
                                     "results/metaxa_bacteria/","") %>% 
  str_replace(".bacteria.fasta","")

write.table(df_reads_all,"results/idtaxa/idtaxa_reads_all.txt",sep="\t",quote = F, row.names = F)

