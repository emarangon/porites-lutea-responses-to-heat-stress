

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
library (DESeq2) #for normalization transcriptomics for mixomics
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



########################################################################
#################### host RNA-seq (transcriptomics) ####################

############################################################
### loading data ####

###import files

files <- dir(path = "./13_mapping_plutea_host_results", pattern = "*plutea_host_relaxed_ReadsPerGene_StrandnessReverse.out.tab", full.names = T) 

host_counts1 <- files %>% purrr::map(read_tsv,  col_names = TRUE ) %>%
  purrr::reduce(cbind)

host_counts1 %>% head()
names(host_counts1) <- gsub(x = names(host_counts1), pattern = "_", replacement = ".") 
colnames(host_counts1)[1] <- "gene_plutea_host_1234567890" 
host_counts1 <- host_counts1 %>%dplyr::select(-ends_with("id")) 
names(host_counts1) <-  gsub('.{11}$', '', names(host_counts1)) 
host_counts1 %>% head() 


###filtering samples

remove <- "P4.31" # because of QC
host_counts2 <- host_counts1[, !(names(host_counts1) %in% remove)] #remove sample
remove <- "P3.33" # because of QC
host_counts3 <- host_counts2[, !(names(host_counts2) %in% remove)] #remove sample
remove <- "P1.36" # because identified as outlier (WGCNA clustering)
host_counts4 <- host_counts3[, !(names(host_counts3) %in% remove)] #remove sample
remove <- "P1.19" # because identified as outlier (WGCNA clustering)
host_counts5 <- host_counts4[, !(names(host_counts4) %in% remove)] #remove sample
remove <- "P3.18" # because identified as outlier (WGCNA clustering) in Symbiodiniaceae samples
host_counts <- host_counts5[, !(names(host_counts5) %in% remove)] #remove sample
head (host_counts)
#I'll use 31 samples (36-5) for downstream analyses


###convert into DEGList object

host_samples <- read_tsv("./SampleInfo_plutea.txt") #metadata
host_samples1<- subset (host_samples, sample_id != "P4.31" & sample_id != "P3.33" & sample_id != "P1.19" & sample_id != "P1.36" & sample_id != "P3.18")
gn <- host_counts %>% dplyr::select (-"gene_plutea_host")
table(colnames(gn) == host_samples1$sample_id) #match

host_counts_matrix <- host_counts %>%
  dplyr::select(-gene_plutea_host) %>% 
  as.matrix()
row.names(host_counts_matrix) <- host_counts$gene_plutea_host 
head (host_counts_matrix)

host_DGE2 <- DGEList(host_counts_matrix, samples = host_samples1) #I have a tot of 31,126 genes (some have zero counts tho)


###filtering genes

table(rowSums(host_DGE2$counts==0)==31) #how many genes have zero counts across all 31 samples? 642
host_filt <- filterByExpr(host_DGE2, design = model.matrix(~TreatTime + genotype, 
                                                      data = host_DGE2$samples), min.count = 20)
host_filtered1 <- host_DGE2[host_filt, , keep.lib.sizes = F]
host_filtered1 #24,171 genes after filtering
table(rowSums(host_filtered1$counts==0)==31) #how many genes have zero counts across all 31 samples? zero


###normalization using DESeq2

