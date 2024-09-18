# This script compiles the database used for IDTAXA. 
# This step only needs to be run once. 


# Package names
packages <- c("tidyverse",
              "DECIPHER",
              "R.utils")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


# idTaxa - script that assigned taxonomic positions to bacterial rDNA (assemblies and reads) based on GTDB
#setwd("~/Documents/gulya/coverage")
library(tidyverse)
library(DECIPHER)
source('scripts/utils.R')
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

# Save this for later, so we can load it back in! 
saveRDS(trainingSet, file = 'data/db/gtdb_trainingset.rds')




