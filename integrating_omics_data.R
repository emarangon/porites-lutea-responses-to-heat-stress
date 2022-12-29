

##############################################################################################
################### plutea - heat stress experiment - integrating omics data #################
##############################################################################################

##author EMMA MARANGON


##########################################################
#################### PRE-PROCESSING ######################
##########################################################


###loading libraries

library(tidyverse) #includes ggplot2, dplyr, tibble, tidyR
library(phyloseq) #microbial analyses
library(edgeR) #transcriptomic analyses
library(limma)  #transcriptomic analyses
library(mixOmics) #intergrating omics data


#############################################################
#################### 16S data (microbes) ####################

############################################################
### loading data ####

##1 ASV TABLE
otu_table = read.csv("./heat-feature-table_R.txt", sep = '\t',  dec = ".", check.names = FALSE, row.names=1)

##2 TAXA TABLE
otu_matrix = read.csv("./heat-taxonomy_R.tsv", sep = '\t', header=T, row.names=1)
TAXA_TABLE_split<-otu_matrix %>% separate(Taxon, c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep=";", remove=TRUE)
TAXA_TABLE_matrix<-as.matrix(TAXA_TABLE_split)

##3 METADATA
metadata = read.csv("./heat-metadata.txt", sep = '\t', row.names=1)
metadata1 <- metadata %>% unite ("time.genotype", timePoint, sampleType, genotype, sep=".", remove=FALSE, na.rm = TRUE)
metadata2 <- metadata1 %>% unite ("TreatTime", treatment, timePoint, sep=".", remove=FALSE, na.rm = TRUE) %>%
  unite ("time.genotype.treat", treatment, timePoint, sampleType, genotype, sep=".", remove=FALSE, na.rm = TRUE) 
  
### import ###
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table (TAXA_TABLE_matrix)
meta = sample_data(metadata2)
phy_tree = read_tree("heat-tree.nwk")  

### merging ###
phyloseq_merged = phyloseq(OTU, TAX, meta, phy_tree)


############################################################
### filtering ####

### filtering eukaryota, mitochondria, chloroplasts ###
phyloseq_merged_filtered <- phyloseq_merged %>% subset_taxa(Domain != "D_0__Eukaryota" & Family  != "D_4__Mitochondria" & Order   != "D_3__Chloroplast")

### filtering out samples not belonging to the plutea experiment
phyloseq_merged_filtered_Porites <- phyloseq_merged_filtered %>% subset_samples(experiment == "Porites")

### removing contaminants (see plutea_16S_microbes.R)
contaminants_only <- c("fd98346394a5c79e554003012cb33826") #based on results decontam
AllTaxa = taxa_names(phyloseq_merged_filtered_Porites)
AllTaxa <- AllTaxa[!(AllTaxa %in% contaminants_only)]
phyloseq_merged_filtered_Porites_NOcontaminants = prune_taxa(AllTaxa, phyloseq_merged_filtered_Porites)

### more filtering (keeping only samples used also for transcriptomics)
phyloseq_merged_filtered_Porites_NOcontaminants_noP2 <- subset_samples(phyloseq_merged_filtered_Porites_NOcontaminants, genotype != "P2")
phyloseq_merged_filtered_Porites_noP2_mixomics <- phyloseq_merged_filtered_Porites_NOcontaminants_noP2 %>% subset_samples(sampleType != "food" & sampleType != "seawater" & sampleType != "blank" & sampleType != "negative" & 
                                                                                            TreatTime != "baseline.baseline" & timePoint != "T1"  & timePoint != "T5" & timePoint != "T3" &
                                                                                            genotype != "P2") 
phyloseq_merged_filtered_Porites_noP2_mixomics2 <- subset_samples(phyloseq_merged_filtered_Porites_noP2_mixomics, 
                                                                             fragID == "P1.12" | fragID == "P1.13" | fragID == "P1.19" | fragID == "P1.21" | 
                                                                               fragID == "P1.24" | fragID == "P1.25" | fragID == "P1.26" | fragID == "P1.28" |
                                                                               fragID == "P1.31" | fragID == "P1.36" | fragID == "P1.6 " | fragID == "P1.8 " | 
                                                                               fragID == "P3.10" | fragID == "P3.12" | fragID == "P3.18" | fragID == "P3.2 " | 
                                                                               fragID == "P3.21" | fragID == "P3.30" | fragID == "P3.31" | fragID == "P3.32" | 
                                                                               fragID == "P3.5 " | fragID == "P3.6 " | fragID == "P3.8 " | fragID == "P3.33" |
                                                                               fragID == "P4.14" | fragID == "P4.15" | fragID == "P4.17" | fragID == "P4.18" | 
                                                                               fragID == "P4.21" | fragID == "P4.24" | fragID == "P4.25" | fragID == "P4.3 " | 
                                                                               fragID == "P4.7 " | fragID == "P4.8 " | fragID == "P4.9 " | fragID == "P4.31")                                                                                            
phyloseq_merged_filtered_Porites_noP2_mixomics3 <- subset_samples( phyloseq_merged_filtered_Porites_noP2_mixomics2, 
                                  fragID != "P4.31" & fragID != "P3.33" & fragID != "P3.18" & fragID != "P1.19" & fragID != "P1.36") #removing samples not included in final anayses host or symb transcriptomics


### removing low abundance ASVs
low_ab <- phyloseq::genefilter_sample(phyloseq_merged_filtered_Porites_noP2_mixomics3, filterfun_sample(function(x) {x / sum(x)} > 1e-5))
phyloseq_merged_filtered_Porites_noP2_mixomics4 <- phyloseq::prune_taxa(low_ab, phyloseq_merged_filtered_Porites_noP2_mixomics3)
phyloseq_merged_filtered_Porites_noP2_mixomics4 = prune_taxa(taxa_sums(phyloseq_merged_filtered_Porites_noP2_mixomics4) > 0, phyloseq_merged_filtered_Porites_noP2_mixomics4) #remove unobserved ASVs (sum 0 across all samples)
#from 22612 to 3916 ASVs

### removing singletons  
phyloseq_merged_filtered_Porites_noP2_mixomics4_NOsingletons <- prune_taxa(taxa_sums(phyloseq_merged_filtered_Porites_noP2_mixomics4) >=2, phyloseq_merged_filtered_Porites_noP2_mixomics4)
phyloseq_merged_filtered_Porites_noP2_mixomics_final = prune_taxa(taxa_sums(phyloseq_merged_filtered_Porites_noP2_mixomics4_NOsingletons) > 0, phyloseq_merged_filtered_Porites_noP2_mixomics4_NOsingletons) #remove unobserved ASVs (sum 0 across all samples)
#no change


############################################################
### importing data into mixomics (following http://mixomics.org/mixmc/mixmc-preprocessing/)

taxo <- tax_table(phyloseq_merged_filtered_Porites_noP2_mixomics_final) # extraction of the taxonomy
meta.data <- phyloseq_merged_filtered_Porites_noP2_mixomics_final@sam_data # extraction of the metadata
data.raw <- t(otu_table(phyloseq_merged_filtered_Porites_noP2_mixomics_final)) # extract OTU table from phyloseq object # samples should be in row and variables in column
data.offset <- data.raw+1 ## STEP 1: OFFSET

pca.result <- pca(data.offset, logratio = 'CLR') # undergo PCA after CLR transformation
plotIndiv(pca.result, group = meta.data$TreatTime, legend=TRUE, legend.title = 'Treatment and Time',
          pch = as.factor(meta.data$genotype), legend.title.pch = 'Genotype',
          ind.names = FALSE,
          title = 'PCA plutea 16S, PCA Comps 1&2') 

data.offset_transf <- logratio.transfo(data.offset, logratio = 'CLR') ##STEP 2: TRANSFORMATION
rownames(data.offset_transf)

