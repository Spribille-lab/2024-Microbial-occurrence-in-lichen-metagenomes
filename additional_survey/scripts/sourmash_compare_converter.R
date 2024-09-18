# load in tidyverse library
library(tidyverse)

# read in .csv file
sourmash_compare <- read.csv('results/sourmash_compare_matrix.csv') %>% 
  mutate(names = names(.)) %>%  # create a column for metagenome names
  column_to_rownames(var = 'names') %>%  # move these names to rownames
  as.matrix() %>% # convert the .csv 'matrix' to an actual matrix
  as.table() %>% # turn it in to a table
  as.data.frame() %>% # now a data.frame. Yes... all this is necessary
  rename('metagenome_1' = 'Var1', # rename the first column to 'metagenome_1'
         'metagenome_2' = 'Var2', # rename second column to 'metagenome_2'
         'sourmash_score' = 'Freq') %>% # rename the final colum to 'sourmash_score'
  filter(metagenome_1 != metagenome_2) %>% # remove same-same comparisons
  arrange(desc(sourmash_score)) # arrange by descending sourmash score

# write it out so you have a copy! 
write.csv(sourmash_compare, 
          file = 'results/sourmash_compare_results.csv', 
          row.names = FALSE)

# Display the top 20 comparisons.
print('Here are the top 20 sourmash scores')

head(sourmash_compare, 20)
