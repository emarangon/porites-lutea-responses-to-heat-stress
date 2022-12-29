#####################################################################################################
################### plutea - heat stress experiment - host transcriptomics analyses #################
#####################################################################################################

##EMMA MARANGON


##########################################################
#################### PRE-PROCESSING ######################
##########################################################

###loading libraries
library(tidyverse) #data wrangling (includes purrr, ggplot2, dplyr, tibble -> I need them all)
library(edgeR)
library(limma)
library(KOGMWU) #for kog.mwu test


###import files
files <- dir(path = "./13_mapping_plutea_host_results", pattern = "*plutea_host_relaxed_ReadsPerGene_StrandnessReverse.out.tab", full.names = T) 

host_counts1 <- files %>% purrr::map(read_tsv,  col_names = TRUE ) %>%  #map (function of purrr, which is part of tidyverse) returns a list
  purrr::reduce(cbind) #Reduce a list to a single value by iteratively applying a binary function; reduce() is an operation that combines the elements of a vector into a single value # The cbind function is used to combine vectors, matrices and/or data frames by columns.

host_counts1 %>% head() #now I have one column per sample named gene_id (all identical across samples) and one column per sample with counts (sample_counts). 
#I rename only the first column (any gene_id column would be ok as all identical) and then remove all the 'gene-id' columns
names(host_counts1) <- gsub(x = names(host_counts1), pattern = "_", replacement = ".") #repace _ with .
host_counts1 %>% head() #it worked
colnames(host_counts1)[1] <- "gene_plutea_host_1234567890" #I rename first column
host_counts1 %>% head() #it worked
host_counts1 <- host_counts1 %>%dplyr::select(-ends_with("id")) #I remove all gene_id columns
host_counts1 %>% head() #it worked
names(host_counts1) <-  gsub('.{11}$', '', names(host_counts1)) #I remove last 11 characters from each colname
host_counts1 %>% head() #it worked
