# idTaxa - script that assigned taxonomic positions to bacterial rDNA (assemblies and reads) based on GTDB
#setwd("~/Documents/gulya/coverage")
library(tidyverse)
library(DECIPHER)
source('code/utils.R')
library(R.utils)

#set up database
seqs <- readDNAStringSet("data/db/gtdb-sbdi-sativa.r06rs202.assignTaxonomy.fna")
taxid <- NULL
seqs <- OrientNucleotides(seqs)

# obtain the taxonomic assignments
names(seqs)<-paste("Root;",names(seqs),sep="")
groups <- names(seqs) # sequence names
# assume the taxonomy begins with 'Root;'
groups <- gsub("(.*)(Root;)", "\\2", groups) # extract the group label
groupCounts <- table(groups)
u_groups <- names(groupCounts) # unique groups
length(u_groups) # number of groups


maxGroupSize <- 10 # max sequences per label (>= 1)
remove <- logical(length(seqs))
for (i in which(groupCounts > maxGroupSize)) {
  index <- which(groups==u_groups[i])
  keep <- sample(length(index),
                 maxGroupSize)
  remove[index[-keep]] <- TRUE
}
sum(remove)

maxIterations <- 3 # must be >= 1
allowGroupRemoval <- FALSE
probSeqsPrev <- integer() # suspected problem sequences from prior iteration
for (i in seq_len(maxIterations)) {
  cat("Training iteration: ", i, "\n", sep="")
  # train the classifier
  trainingSet <- LearnTaxa(seqs[!remove],
                           names(seqs)[!remove],
                           taxid)
  # look for problem sequences
  probSeqs <- trainingSet$problemSequences$Index
  if (length(probSeqs)==0) {
    cat("No problem sequences remaining.\n")
    break
  } else if (length(probSeqs)==length(probSeqsPrev) &&
             all(probSeqsPrev==probSeqs)) {
    cat("Iterations converged.\n")
    break
  }
  if (i==maxIterations)
    break
  probSeqsPrev <- probSeqs
  # remove any problem sequences
  index <- which(!remove)[probSeqs]
  remove[index] <- TRUE # remove all problem sequences
  if (!allowGroupRemoval) {
    # replace any removed groups
    missing <- !(u_groups %in% groups[!remove])
    missing <- u_groups[missing]
    if (length(missing) > 0) {
      index <- index[groups[index] %in% missing]
      remove[index] <- FALSE # don't remove
    }
  }
}
sum(remove) # total number of sequences eliminated
length(probSeqs) # 

plot(trainingSet)

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
filenames <- list.files(path = "analysis/03_metagenome_reanalysis/")
filenames <- filenames[grepl("bacteria.fasta",filenames)]
filenames_full<-paste0("analysis/03_metagenome_reanalysis/",filenames)

#remove empty files
nlines<-sapply(filenames_full,countLines)
filenames_nonempty<-filenames_full[nlines>0]

#process assembly files
filenames_assembly<-filenames_nonempty[grepl("assembly_",filenames_nonempty)]

l<-lapply(filenames_assembly,process_fasta)
df_assembly <- do.call(rbind,l)
df_assembly$metagenome<-str_replace(df_assembly$fil,"analysis/03_metagenome_reanalysis/assembly_","") %>% str_replace(".bacteria.fasta","")
write.table(df_assembly,"analysis/03_metagenome_reanalysis/idtaxa_assemblies.txt",sep="\t",quote = F, row.names = F)



#reads
filenames_reads<-filenames_nonempty[!(filenames_nonempty %in% filenames_assembly)]
lreads_1<-process_fasta(filenames_reads[1])
write.table(lreads_1,"analysis/03_metagenome_reanalysis/idtaxa_reads_1.txt",sep="\t",quote = F, row.names = F)

lreads_2<-process_fasta(filenames_reads[2])
write.table(lreads_2,"analysis/03_metagenome_reanalysis/idtaxa_reads_2.txt",sep="\t",quote = F, row.names = F)


lreads_3<-lapply(filenames_reads[3:6],process_fasta)
df_reads_3<-do.call(rbind,lreads_3)
write.table(df_reads_3,"analysis/03_metagenome_reanalysis/idtaxa_reads_3.txt",sep="\t",quote = F, row.names = F)

lreads_4<-lapply(filenames_reads[7:15],process_fasta)
df_reads_4<-do.call(rbind,lreads_4)
write.table(df_reads_4,"analysis/03_metagenome_reanalysis/idtaxa_reads_4.txt",sep="\t",quote = F, row.names = F)

lreads_5<-lapply(filenames_reads[16:100],process_fasta)
df_reads_5<-do.call(rbind,lreads_5)
write.table(df_reads_5,"analysis/03_metagenome_reanalysis/idtaxa_reads_5.txt",sep="\t",quote = F, row.names = F)


lreads_6<-lapply(filenames_reads[101:length(filenames_reads)],process_fasta)
df_reads_6<-do.call(rbind,lreads_6)
write.table(df_reads_6,"analysis/03_metagenome_reanalysis/idtaxa_reads_6.txt",sep="\t",quote = F, row.names = F)


#combine all together
df_reads_all<-rbind(lreads_1,lreads_2,df_reads_3,df_reads_4,df_reads_5,df_reads_6)
df_reads_all$metagenome<-str_replace(df_reads_all$fil,"analysis/03_metagenome_reanalysis/","") %>% str_replace(".bacteria.fasta","")

write.table(df_reads_all,"analysis/03_metagenome_reanalysis/idtaxa_reads.txt",sep="\t",quote = F, row.names = F)







