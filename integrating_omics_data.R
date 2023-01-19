

######################################################################################################################
################### plutea - heat stress experiment - integrating omics data using DIABLO (mixOmics) #################
######################################################################################################################

##author EMMA MARANGON


##########################################################
#################### PRE-PROCESSING ######################
##########################################################


###loading libraries

library(phyloseq) #microbial analyses
library(plyr) #calculate means etc for metadata
library(tidyverse) #includes ggplot2, dplyr, tibble, tidyR 
library(mixOmics)#intergrating omics data
library(edgeR)#transcriptomic analyses
library(limma)#transcriptomic analyses
library (DESeq2) #for normalization for WGCNA


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

### more filtering (keeping only samples used for transcriptomics)
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
remove <- "P1.36" # because identified as outlier (WGCNA clustering) in host samples
host_counts4 <- host_counts3[, !(names(host_counts3) %in% remove)] #remove sample
remove <- "P1.19" # because identified as outlier (WGCNA clustering) in host samples
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


#### physiology

pr_porites <- read.csv("./PR_Porites.csv", header = TRUE)
pr_porites1 <- dplyr::filter (pr_porites, TimePoint =="T0" | TimePoint =="T2" | TimePoint =="T4") #keep data points used for transcriptomics
pr_porites2 <- pr_porites1 %>% dplyr::select (2, 9, 18, 19, 21) #select columns
names(pr_porites2)[names(pr_porites2) == 'Sample.ID'] <- 'sample_id' #rename column Sample.ID
bleaching_porites <- read.csv("./BleachingScores_porites.csv", header = TRUE)
bleaching_porites1 <- dplyr::filter (bleaching_porites, TimePoint =="T0" | TimePoint =="T2" | TimePoint =="T4")
bleaching_porites2 <- bleaching_porites1 %>% dplyr::select (3, 13, 32)
names(bleaching_porites2)[names(bleaching_porites2) == 'Sample.ID'] <- 'sample_id' #rename column Sample.ID
Light_1_filtered_sd <- read.csv("./Light_porites_filteringSD.csv", header = TRUE, sep=",")
Light_2_filtered_sd <- unite(Light_1_filtered_sd, ID2, TimePoint:Light_Dark, remove=FALSE)
Light_2_filtered_sd$ID2 <- as.factor(Light_2_filtered_sd$ID2)
Light_3_filtered_sd <- unite(Light_2_filtered_sd, ID3, ID2:Sample_ID, remove=FALSE)
Light_3_filtered_sd$ID3 <- as.factor(Light_3_filtered_sd$ID3)
Light_4_filtered_sd <- unite(Light_3_filtered_sd, ID4, ID3:Tank, remove=FALSE)
Light_4_filtered_sd$ID4 <- as.factor(Light_4_filtered_sd$ID4)
Light_5_filtered_sd <- unite(Light_4_filtered_sd, ID5, ID4:Treatment, remove=FALSE)
Light_5_filtered_sd$ID5 <- as.factor(Light_5_filtered_sd$ID5)
Light_6_filtered_sd <- ddply(Light_5_filtered_sd,~ID5,summarise,meanY=mean(Y)) #to calculate Y means
Light_7_filtered_sd <-Light_6_filtered_sd %>%
  separate (ID5, c("TimePoint", "Light_Dark", "Sample_ID", "Tank", "Treatment"), "_") #to separate ID5 into multiple columns
Light_7_filtered_sd$Tank <- as.factor(Light_7_filtered_sd$Tank)
Light_7_filtered_sd$Light_Dark <- factor(Light_7_filtered_sd$Light_Dark)
Light_7_filtered_sd$TimePoint <- factor(Light_7_filtered_sd$TimePoint)
Light_7_filtered_sd$Treatment <- factor(Light_7_filtered_sd$Treatment)
Light_7_filtered_sd$Sample_ID <- factor(Light_7_filtered_sd$Sample_ID)
Light_7_filtered_sd$Parent = Light_7_filtered_sd$Sample_ID #copy column Sample_ID in an identical column called Parent
Light_7_filtered_sd %>% mutate(Parent = substr(Parent, 1, 2)) -> photefficiency_porites #edit Parent column (I keep only the first two values)
names(photefficiency_porites)[names(photefficiency_porites)=="meanY"] <- "PhotochemicalEfficiency" #rename 'mean' into 'PhotochemicalEfficiency'
summary(photefficiency_porites)
photefficiency_porites1 <- dplyr::filter (photefficiency_porites, TimePoint =="T0" | TimePoint =="T2" | TimePoint =="T4") 
photefficiency_porites2 <- photefficiency_porites1 %>% dplyr::select (1, 3, 6)
names(photefficiency_porites2)[names(photefficiency_porites2) == 'Sample_ID'] <- 'sample_id' #rename column Sample.ID
meta_1 <- host_filtered1$samples
meta_2<-  dplyr::left_join(meta_1, pr_porites2, by=c('sample_id'='sample_id', 'time'='TimePoint')) #keep only samples used also for transcriptomics
meta_3 <-  dplyr::left_join (meta_2, bleaching_porites2, by=c('sample_id'='sample_id', 'time'='TimePoint'))
meta_df <-  dplyr::left_join (meta_3, photefficiency_porites2, by=c('sample_id'='sample_id', 'time'='TimePoint'))
head(meta_df)

physio.data <- meta_df %>% column_to_rownames("sample_id") %>% dplyr::select (Respiration, Photosynthesis, RatioPR, BleachigScore, PhotochemicalEfficiency, genotype, TreatTime, DHW, stress)
physio.data <- as.data.frame(physio.data)
dim(physio.data)
head(physio.data)
rownames (physio.data)
rownames(physio.data) <- c(220,240,312, 297, 
                           184, 190, 286, 234,
                           135, 210, 266, 216,
                           263, 200, 142, 152,
                           231, 179, 140, 284,
                           180, 136, 177, 211,
                           315, 173, 199, 178,
                           202, 318, 242) #need to rename

order <- c("135", "136", "140", "142", "152", "173", "177", "178", "179", "180", "184", "190", "199", "200", "202", "210", "211", "216", "220", "231", "234", "240", "242", "263", "266", "284", "286", "297", "312", "315", "318")
physio.data_df <- as.data.frame(physio.data) %>% rownames_to_column("sample")
physio.data_df$sample <- as.factor(physio.data_df$sample)
YwFactors <- physio.data_df[match(order, physio.data_df$sample),] %>% remove_rownames() %>% column_to_rownames("sample") #reordering datasets


###normalization using DESeq2

host_for_deseq <- host_filtered1$counts
host_de_input = as.matrix(host_for_deseq)
meta_df$TreatTime <- factor(meta_df$TreatTime)
levels(meta_df$TreatTime)
all.equal(colnames(host_de_input), meta_df$sample_id) # they match
host_dds <- DESeqDataSetFromMatrix(round(host_de_input),
                              meta_df,
                              design = ~ TreatTime + genotype) 
host_dds_norm <- vst(host_dds, blind=FALSE)


#####I keep only genes with KEGG annotations

host_annotations <- read.csv("./eggNOG_output_plutea_annotation_ALLgo.tsv", sep = '\t', header=TRUE)
host_annotations$query <- gsub (".m1", "", as.character(host_annotations$query)) #to make data match
host_annotations$query <- gsub (".model.", ".TU.", as.character(host_annotations$query)) #to make data match
host_annotations_KEGG = host_annotations[,c(1, 12)] #keep only gene and ko columns
dim (host_annotations_KEGG)
host_genes_to_keep <- host_annotations_KEGG[!grepl("-", host_annotations_KEGG$KEGG_ko),] #remove unknown pathways
host_m <- assay (host_dds_norm) 
dim (host_m)
host_d <- as.data.frame(host_m) %>% rownames_to_column("query")
X_host_df_FINAL_only_annotated <- semi_join(host_d, host_genes_to_keep, by = "query") #filterig out
dim (X_host_df_FINAL_only_annotated)
X_host <- X_host_df_FINAL_only_annotated %>% column_to_rownames("query") %>% as.matrix() %>% t


### filtering samples

rownames (X_host)
rownames(X_host) <- c(220,240,312, 297, 
                      184, 190, 286, 234,
                      135, 210, 266, 216,
                      263, 200, 142, 152,
                      231, 179, 140, 284,
                      180, 136, 177, 211,
                      315, 173, 199, 178,
                      202, 318, 242) 
order <- c("135", "136", "140", "142", "152", "173", "177", "178", "179", "180", "184", "190", "199", "200", "202", "210", "211", "216", "220", "231", "234", "240", "242", "263", "266", "284", "286", "297", "312", "315", "318")
X_host_df <- as.data.frame(X_host) %>% rownames_to_column("sample")
X_host_df$sample <- as.factor(X_host_df$sample)
X_host_df_final <- X_host_df[match(order, X_host_df$sample),] %>% remove_rownames() %>% column_to_rownames("sample") #reordering datasets
dim(X_host_df_final) #12,435 genes




###################################################################################
#################### Symbiodiniaceae RNA-seq (transcriptomics) ####################

############################################################
### loading data ####

###import files

files_symb <- dir(path = "./13_mapping_plutea_UQ_cladocopium_goreaui_results", pattern = "*plutea_UQ_cladocopium_goreaui_relaxed_ReadsPerGene_StrandnessReverse.out.tab", full.names = T) 
cladocopium_counts1 <- files_symb %>% purrr::map(read_tsv,  col_names = TRUE ) %>%
  purrr::reduce(cbind)
names(cladocopium_counts1) <- gsub(x = names(cladocopium_counts1), pattern = "_", replacement = ".")
colnames(cladocopium_counts1)[1] <- "gene_plutea_cladocopium_1234567890" 
cladocopium_counts1 <- cladocopium_counts1 %>%dplyr::select(-ends_with("id")) 
names(cladocopium_counts1) <-  gsub('.{11}$', '', names(cladocopium_counts1)) 


###filtering samples

remove <- "P4.31" # because of QC
cladocopium_counts2 <- cladocopium_counts1[, !(names(cladocopium_counts1) %in% remove)] #remove sample
remove <- "P3.33" # because of QC
cladocopium_counts3 <- cladocopium_counts2[, !(names(cladocopium_counts2) %in% remove)] #remove sample
remove <- "P1.36" # because identified as outlier (WGCNA clustering) in host samples
cladocopium_counts4 <- cladocopium_counts3[, !(names(cladocopium_counts3) %in% remove)] #remove sample
remove <- "P1.19"  # because identified as outlier (WGCNA clustering) in host samples
cladocopium_counts5 <- cladocopium_counts4[, !(names(cladocopium_counts4) %in% remove)] #remove sample
remove <- "P3.18" # because identified as outlier (WGCNA clustering) in Symbiodiniaceae samples
cladocopium_counts <- cladocopium_counts5[, !(names(cladocopium_counts5) %in% remove)] #remove sample
head (cladocopium_counts)
#I'll use 31 samples (36-5) for downstream analyses


###convert into DEGList object

cladocopium_samples <- read_tsv("./SampleInfo_plutea.txt") #metadata
cladocopium_samples1<- subset (cladocopium_samples, sample_id != "P4.31" & sample_id != "P3.33" & sample_id != "P3.18" & sample_id != "P1.19" & sample_id != "P1.36")
gn <- cladocopium_counts %>% dplyr::select (-"gene_plutea_cladocopium")
table(colnames(gn) == cladocopium_samples1$sample_id) #they match

cladocopium_counts_matrix <- cladocopium_counts %>%
  dplyr::select(-gene_plutea_cladocopium) %>% #remove gene_plutea_cladocopium
  as.matrix()
row.names(cladocopium_counts_matrix) <- cladocopium_counts$gene_plutea_cladocopium #add gene_plutea_cladocopium
head (cladocopium_counts_matrix)

cladocopium_DGE2 <- DGEList(cladocopium_counts_matrix, samples = cladocopium_samples1) #I have a tot of 39,006 genes (some have zero counts tho)


###filtering genes

table(rowSums(cladocopium_DGE2$counts==0)==31) #how many genes have zero counts across all 33 samples? 9836
filt <- filterByExpr(cladocopium_DGE2, design = model.matrix(~TreatTime + genotype, 
                                                             data = cladocopium_DGE2$samples), min.count = 20)
cladocopium_filtered1 <- cladocopium_DGE2[filt, , keep.lib.sizes = F]
cladocopium_filtered1 #23,810 genes after filtering
table(rowSums(cladocopium_filtered1$counts==0)==31) #how many genes have zero counts across all 31 samples? zero


###normalization using DESeq2

cladocopium_for_deseq <- cladocopium_filtered1$counts
cladocopium_de_input = as.matrix(cladocopium_for_deseq)
meta_df$TreatTime <- factor(meta_df$TreatTime)
levels(meta_df$TreatTime)
all.equal(colnames(cladocopium_de_input), meta_df$sample_id) #they match
cladocopium_dds <- DESeqDataSetFromMatrix(round(cladocopium_de_input),
                                   meta_df,
                                   design = ~ TreatTime + genotype)
cladocopium_dds_norm <- vst(cladocopium_dds, blind = FALSE)


#####I keep only genes with KEGG annotations

cladocopium_annotations <- read.csv("./eggNOG_output_cladocopium_annotation_ALLgo.tsv", sep = '\t', header=TRUE)
cladocopium_annotations$query <- gsub (".mRNA1", "", as.character(cladocopium_annotations$query)) #to make data match
cladocopium_annotations_KEGG = cladocopium_annotations[,c(1, 12)] #keep only gen and ko columns
dim(cladocopium_annotations_KEGG)
cladocopium_genes_to_keep <- cladocopium_annotations_KEGG[!grepl("-", cladocopium_annotations_KEGG$KEGG_ko),] #remove unknown pathways
dim(cladocopium_genes_to_keep)
cladocopium_m <- assay (cladocopium_dds_norm) 
dim (cladocopium_m)
cladocopium_d <- as.data.frame(cladocopium_m) %>% rownames_to_column("query")
dim(cladocopium_d)
X_cladocopium_df_FINAL_only_annotated <- semi_join(cladocopium_d, cladocopium_genes_to_keep, by = "query") #filterig out
dim (X_cladocopium_df_FINAL_only_annotated)
X_cladocopium <- X_cladocopium_df_FINAL_only_annotated %>% column_to_rownames("query") %>% as.matrix() %>% t
rownames (X_cladocopium)
rownames(X_cladocopium) <- c(220,240,312, 297, 
                      184, 190, 286, 234,
                      135, 210, 266, 216,
                      263, 200, 142, 152,
                      231, 179, 140, 284,
                      180, 136, 177, 211,
                      315, 173, 199, 178,
                      202, 318, 242) 
order <- c("135", "136", "140", "142", "152", "173", "177", "178", "179", "180", "184", "190", "199", "200", "202", "210", "211", "216", "220", "231", "234", "240", "242", "263", "266", "284", "286", "297", "312", "315", "318")
X_cladocopium_df <- as.data.frame(X_cladocopium) %>% rownames_to_column("sample")
X_cladocopium_df$sample <- as.factor(X_cladocopium_df$sample)
X_cladocopium_df_final <- X_cladocopium_df[match(order, X_cladocopium_df$sample),] %>% remove_rownames() %>% column_to_rownames("sample") #reordering datasets
dim(X_cladocopium_df_final) #7,627 genes



##############################################################################################
########################################## DIABLO  ###########################################
##############################################################################################
#following http://mixomics.org/mixdiablo/diablo-tcga-case-study/


#I want to apply multilevel analyses to take into account the effect of Parental colony (i.e. genotype)
YwFactors$genotype = as.factor(YwFactors$genotype)
design <- data.frame(genotype = YwFactors$genotype)
Xw_rna_symb<-withinVariation(X_cladocopium_df_final,design=design)
Xw_rna_host<-withinVariation(X_host_df_final,design=design)
Xw_amplicon<-withinVariation(data.offset_transf,design=design)

data = list (rna_symb = Xw_rna_symb, 
             rna_host = Xw_rna_host, 
             amplicon = Xw_amplicon
             )
lapply(data, dim)

YwFactors$stress = as.factor(YwFactors$stress)
Y = YwFactors$stress 
summary(Y)


#### Pairwise PLS Comparisons - generate three pairwise PLS models to understand the correlation between the datasets

list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(25, 25)

pls1 <- spls(data[["rna_symb"]], data[["rna_host"]], 
             keepX = list.keepX, keepY = list.keepY) 
pls2 <- spls(data[["rna_symb"]], data[["amplicon"]], 
             keepX = list.keepX, keepY = list.keepY)
pls3 <- spls(data[["rna_host"]], data[["amplicon"]], 
             keepX = list.keepX, keepY = list.keepY)

plotVar(pls1, cutoff = 0.5, title = "(a) rna_symb vs rna_host", 
        legend = c("rna_symb", "rna_host"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

plotVar(pls2, cutoff = 0.5, title = "(b) rna_symb vs amplicon", 
        legend = c("rna_symb", "amplicon"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

plotVar(pls3, cutoff = 0.5, title = "(c) rna_host vs amplicon", 
        legend = c("rna_host", "amplicon"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

cor(pls1$variates$X, pls1$variates$Y) 
cor(pls2$variates$X, pls2$variates$Y) 
cor(pls3$variates$X, pls3$variates$Y) 
#data sets are highly correlated


####initial DIABLO model 

design = matrix(0.6, ncol = length(data), nrow = length(data),
                dimnames = list(names(data), names(data)))
diag(design) = 0 # set diagonal to 0s
design

basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 5, design = design)

perf.diablo = perf(basic.diablo.model, validation = 'Mfold',
                   folds = 5, nrepeat = 100) 
plot(perf.diablo)

# set the optimal ncomp value
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "mahalanobis.dist"] 
ncomp #n=2 comp are optimal
perf.diablo$choice.ncomp$WeightedVote

# tuning the number of features
#set grid (the model was first run with bigger grids; based on that results, the following values were chosen)
test.keepX = list (rna_symb = c(5:9, seq(10, 50, 3)), 
                   rna_host = c(5:9, seq(10, 50, 3)),
                   amplicon = c(5:9, seq(10, 50, 3)))

tune.diablo = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, test.keepX = test.keepX, design = design, validation = 'Mfold', folds = 5, nrepeat = 50, dist = "mahalanobis.dist")
save (tune.diablo, file="./tune.diablo.R.data")

head (tune.diablo$error.rate)
head (tune.diablo$error.rate.class)
list.keepX = tune.diablo$choice.keepX 
list.keepX

#final model
final.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                                  keepX = list.keepX, design = design)
save (final.diablo.model, file="./final.diablo.R.data")
final.diablo.model$design

#plots
plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO',
          col.per.group = c("#ECBE92", "#F7794D", "#B3B4B4"),
          pch.size = 2)

circosPlot(final.diablo.model, cutoff = 0.8, line = TRUE,
           color.blocks= c('#0E5E6F', '#BAD1C2','#F2DEBA'),
           color.cor = c("#EEE3CB", "#7895B2"), 
           color.Y = c("#ECBE92", "#F7794D", "#B3B4B4"),
           size.variables = 0.2, size.labels = 1, size.legend = 0.7,
           linkWidth = 1)

pdf(file = "loading_diablo_comp1.pdf", height = 50, width = 50)
loading_diablo_comp1 <- plotLoadings(final.diablo.model, comp = 1, contrib = 'max', method = 'median')
dev.off()

pdf(file = "loading_diablo_comp2.pdf", height = 50, width = 50)
loading_diablo_comp2 <- plotLoadings(final.diablo.model, comp = 2, contrib = 'max', method = 'median')
dev.off()
