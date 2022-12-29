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

host_counts1 <- files %>% purrr::map(read_tsv,  col_names = TRUE ) %>%
  purrr::reduce(cbind)

host_counts1 %>% head() #now I have one column per sample named gene_id (all identical across samples) and one column per sample with counts (sample_counts). 
#I rename only the first column (any gene_id column would be ok as all identical) and then remove all the 'gene-id' columns
names(host_counts1) <- gsub(x = names(host_counts1), pattern = "_", replacement = ".") #repace _ with .
colnames(host_counts1)[1] <- "gene_plutea_host_1234567890" #I rename first column
host_counts1 <- host_counts1 %>%dplyr::select(-ends_with("id")) #I remove all gene_id columns
names(host_counts1) <-  gsub('.{11}$', '', names(host_counts1)) #I remove last 11 characters from each colname
host_counts1 %>% head() #it worked

###filtering samples
remove <- "P4.31" # because of QC
host_counts2 <- host_counts1[, !(names(host_counts1) %in% remove)] #remove sample
remove <- "P3.33" # because of QC
host_counts3 <- host_counts2[, !(names(host_counts2) %in% remove)] #remove sample
remove <- "P1.36" # because identified as outlier (WGCNA clustering)
host_counts4 <- host_counts3[, !(names(host_counts3) %in% remove)] #remove sample
remove <- "P1.19" # because identified as outlier (WGCNA clustering)
host_counts <- host_counts4[, !(names(host_counts4) %in% remove)] #remove sample
head (host_counts)
#I'll use 32 samples (36-4) for downstream analyses


###convert into DEGList object

host_samples <- read_tsv("./SampleInfo_plutea.txt") #metadata
host_samples1<- subset (host_samples, sample_id != "P4.31" & sample_id != "P3.33" & sample_id != "P1.19" & sample_id != "P1.36")  #filtering samples (see above)
#it is very important to make sure the metadata is in the same order as the column names of the counts table !!!
gn <- host_counts %>% dplyr::select (-"gene_plutea_host") #Just removing first gene colum to be able to compare samples in the next step
table(colnames(gn) == host_samples1$sample_id) #check order corresponds -> it does!

host_counts_matrix <- host_counts %>%
  dplyr::select(-gene_plutea_host) %>% 
  as.matrix()
row.names(host_counts_matrix) <- host_counts$gene_plutea_host #add gene_plutea_host
head (host_counts_matrix)

host_DGE2 <- DGEList(host_counts_matrix, samples = host_samples1) #I have a tot of 31,126 genes (some have zero counts tho)


###filtering genes

table(rowSums(host_DGE2$counts==0)==32) #how many genes have zero counts across all 32 samples? 625
filt <- filterByExpr(host_DGE2, design = model.matrix(~TreatTime + genotype, 
                                                     data = host_DGE2$samples), min.count = 20) 
host_filtered1 <- host_DGE2[filt, , keep.lib.sizes = F]
dim(host_filtered1) #24,207 genes, 32 samples
table(rowSums(host_filtered1$counts==0)==32) #how many genes have zero counts across all 32 samples?zero

#To check filter cut off before and after
par(mfrow = c(2, 2)) # Create a 2 x 2 plotting matrix
mean_log_cpm3 <- aveLogCPM(host_DGE2$counts)
filter_threshold <- 0.5
hist(mean_log_cpm3)
abline(v = filter_threshold)
qqnorm(mean_log_cpm3)
abline(h = filter_threshold)
# 
mean_log_cpm4 <- aveLogCPM(host_filtered1$counts)
filter_threshold <- 0.5
hist(mean_log_cpm4)
abline(v = filter_threshold)
qqnorm(mean_log_cpm4)
abline(h = filter_threshold)
#improved a lot!
par(mfrow = c(1, 1)) #back to stnadard visualization


###normalization
host_filtered2 <- calcNormFactors(host_filtered1,method = "TMM") #calcNormFactors doesn’t normalize the data, it just calculates normalization factors for use downstream !!
host_filtered2$samples$norm.factors #all scaling factors are relatively close to 1 -> the effect of TMM-normalisation is mild for my data
table(host_samples1$sample_id == colnames(host_filtered2)) #double check data are still matching
host_filtered2_cpm <- cpm(host_filtered2, log=TRUE)  #converting into log-cpm





#######################################################################
#################### DEGs analyses (limma-vooom) ######################
#######################################################################

treatment <- host_filtered2$samples$treatment
treatment <- factor(treatment)
time <- host_filtered2$samples$time
time <- factor(time)
group <- interaction(treatment, time)
group
group <- factor(group)
levels (group)
genotype <- host_filtered2$samples$genotype
genotype <- factor(genotype)
mm <- model.matrix(~0 + group)
v <- voom(host_filtered2, mm, plot = T)
cor <- duplicateCorrelation(v, mm, block=genotype)
cor$consensus.correlation #0.4885872
fit <- lmFit(object=v, design=mm, 
             block=genotype, correlation=cor$consensus.correlation)
head(coef(fit))
contr <- makeContrasts((
                       groupHeat.T2 - (groupHeat.T0 + groupAmbient.T0 + groupAmbient.T2)/3,
                       groupHeat.T4 - (groupHeat.T0 + groupAmbient.T0 + groupAmbient.T2 + groupAmbient.T4)/4,
                       levels = colnames(coef(fit)))
contr
fit.cont <- contrasts.fit(fit, contr)
e.fit.cont <- eBayes(fit.cont)
de <- decideTests(e.fit.cont)
summary(de) 
#0 down and 0 up in low relative to no stress
#48 down and 103 up in moderate relative to no stress
                 
                       
                       

############################################################################
#################### KOG enrichment analyses (KOGMWU) ######################
############################################################################
                       
host_annotations <- read.csv("./eggNOG_output_plutea_annotation_ALLgo.tsv", sep = '\t', header=TRUE) #annotation reference genome using EggNOG mapper
host_annotations_KOG = host_annotations[,c(1, 7)] #keep only gene and KOG                       
host_annotations_KOG$query <- gsub (".m1", "", as.character(host_annotations_KOG$query))  #editing gene names to match my reads
host_annotations_KOG$query <- gsub (".model.", ".TU.", as.character(host_annotations_KOG$query))  #editing gene names to match my reads
host_annotations_KOG_noNA <- host_annotations_KOG[!grepl("-", host_annotations_KOG$COG_category),] #removing the unassigned
                       
#I want to have one KOG per gene per row
host_annotations_KOG_multiple <- host_annotations_KOG_noNA %>% separate (COG_category, c("kog","B", "C", "D", "E"), sep=cumsum(c(1,1,1,1)))                       
host_annotations_KOG_multiple <- host_annotations_KOG_multiple %>% mutate_all(na_if,"")
first = host_annotations_KOG_multiple[,c(1,2)]
second = host_annotations_KOG_multiple[,c(1,3)]
second <- na.omit(second) #removing NA lines
third = host_annotations_KOG_multiple[,c(1,4)]
third <- na.omit(third) #removing NA lines
fourth = host_annotations_KOG_multiple[,c(1,5)]
fourth <- na.omit(fourth) #removing NA lines
fifth = host_annotations_KOG_multiple[,c(1,6)]
fifth <- na.omit(fifth) #removing NA lines
colnames(second) <- c("query", "kog") #need to be same colname for dplyr::bind_rows step
colnames(third) <- c("query", "kog")
colnames(fourth) <- c("query", "kog")
colnames(fifth) <- c("query", "kog")
combined12 <- dplyr::bind_rows(first, second) #multiple rows per gene if >1 KOG
combined123 <- dplyr::bind_rows(combined12, third)
combined1234 <- dplyr::bind_rows(combined12, fourth)
host_gene2kog <- dplyr::bind_rows(combined1234, fifth) #Two-column dataframe of gene annotations: gene id, KOG class
head (host_gene2kog)                       
                    
                       
################################
##### host comparisons
                       
treatment <- host_filtered2$samples$treatment
treatment <- factor(treatment)
time <- host_filtered2$samples$time
time <- factor(time)
group <- interaction(treatment, time)
group
group <- factor(group)
genotype <- host_filtered2$samples$genotype
genotype <- factor(genotype)
mm <- model.matrix(~0 + group)
v <- voom(host_filtered2, mm, plot = T)
cor <- duplicateCorrelation(v, mm, block=genotype)                       
fit <- lmFit(object=v, design=mm, 
             block=genotype, correlation=cor$consensus.correlation)
head(coef(fit))
                          
                       
#### HOST no stress vs low stress
                       
contr <- makeContrasts(groupHeat.T2 - (groupHeat.T0 + groupAmbient.T0 + groupAmbient.T2)/3,
                       levels = colnames(coef(fit)))
fit.cont <- contrasts.fit(fit, contr)
e.fit.cont <- eBayes(fit.cont)
top.table <- topTable(e.fit.cont, sort.by = "P", n = Inf)
keep <- tibble::rownames_to_column(top.table, "genes") 
host_only_degs_NOvsLOW_block_logF = keep[,c(1,2)] 
head (host_only_degs_NOvsLOW_block_logF)

host_kogmwu_data = host_only_degs_NOvsLOW_block_logF
row.names(host_kogmwu_data)=NULL 
kog_NOvsLOW_host = kog.mwu(host_kogmwu_data, host_gene2kog, Alternative = "t")
kog_NOvsLOW_host                       
                       
  
#### HOST no stress vs moderate stress                       
                       
contr <- makeContrasts(groupHeat.T4 - (groupHeat.T0 + groupAmbient.T0 + groupAmbient.T2 + groupAmbient.T4)/4,
                       levels = colnames(coef(fit)))
fit.cont <- contrasts.fit(fit, contr)
e.fit.cont <- eBayes(fit.cont)
top.table <- topTable(e.fit.cont, sort.by = "P", n = Inf)
keep <- tibble::rownames_to_column(top.table, "genes") 
host_only_degs_NOvsMODERATE_block_logF = keep[,c(1,2)]
head (host_only_degs_NOvsMODERATE_block_logF) 

host_kogmwu_data = host_only_degs_NOvsMODERATE_block_logF
row.names(host_kogmwu_data)=NULL 
kog_NOvsMODERATE_host = kog.mwu(host_kogmwu_data, host_gene2kog, Alternative = "t")
kog_NOvsMODERATE_host                       
             
                       
################################
##### Symbiodiniaceae comparisons
                       
load ("./cladocopium_filtered2.R.data")
load ("./cladocopium_gene2kog.R.data")                      
                       
treatment_s <- cladocopium_filtered2$samples$treatment
treatment_s <- factor(treatment_s)
time_s <- cladocopium_filtered2$samples$time
time_s <- factor(time_s)
group_s <- interaction(treatment_s, time_s)
group_s
group_s <- factor(group_s)
genotype_s <- cladocopium_filtered2$samples$genotype
genotype_s <- factor(genotype_s)
mm_s <- model.matrix(~0 + group_s)
v_s <- voom(cladocopium_filtered2, mm_s, plot = T)
cor_s <- duplicateCorrelation(v_s, mm_s, block=genotype_s)
fit_s <- lmFit(object=v_s, design=mm_s, 
             block=genotype_s, correlation=cor_s$consensus.correlation)
head(coef(fit_s))                       
       
                       
#### SYMB no stress vs low stress
                       
contr_s <- makeContrasts(group_sHeat.T2 - (group_sHeat.T0 + group_sAmbient.T0 + group_sAmbient.T2)/3,
                       levels = colnames(coef(fit_s)))
fit.cont_s <- contrasts.fit(fit_s, contr_s)
e.fit.cont_s <- eBayes(fit.cont_s)
top.table_s <- topTable(e.fit.cont_s, sort.by = "P", n = Inf)
keep_s <- tibble::rownames_to_column(top.table_s, "genes") 
cladocopium_only_degs_NOvsLOW_block_logF = keep_s[,c(1,2)] 
head (cladocopium_only_degs_NOvsLOW_block_logF)

cladocopium_kogmwu_data = cladocopium_only_degs_NOvsLOW_block_logF
row.names(cladocopium_kogmwu_data)=NULL 
kog_NOvsLOW_s = kog.mwu(cladocopium_kogmwu_data, cladocopium_gene2kog, Alternative = "t")
kog_NOvsLOW_s                       
                       
                       
#### SYMB no stress vs moderate stress 
                       
contr_s <- makeContrasts(group_sHeat.T4 - (group_sHeat.T0 + group_sAmbient.T0 + group_sAmbient.T2 + group_sAmbient.T4)/4,
                       levels = colnames(coef(fit_s)))
fit.cont_s <- contrasts.fit(fit_s, contr_s)
e.fit.cont_s <- eBayes(fit.cont_s)
top.table_s <- topTable(e.fit.cont_s, sort.by = "P", n = Inf)
keep_s <- tibble::rownames_to_column(top.table_s, "genes")
cladocopium_only_degs_NOvsMODERATE_block_logF = keep_s[,c(1,2)]
head (cladocopium_only_degs_NOvsMODERATE_block_logF)

cladocopium_kogmwu_data = cladocopium_only_degs_NOvsMODERATE_block_logF
row.names(cladocopium_kogmwu_data)=NULL 
kog_NOvsMODERATE_s = kog.mwu(cladocopium_kogmwu_data, cladocopium_gene2kog, Alternative = "t")
kog_NOvsMODERATE_s                      
                       
                       
################################
##### KOG plots

### heatmap                       
#compiling a table of delta-ranks to compare the results
ktable1=makeDeltaRanksTable(list("host NO-LOW"=kog_NOvsLOW_host, "host NO-MODERATE"=kog_NOvsMODERATE_host, "symb NO-LOW"=kog_NOvsLOW_s, "symb NO-MODERATE"=kog_NOvsMODERATE_s))
ktable <- ktable1 %>% subset(rownames(ktable1) !="S") %>% droplevels() #remove genes associated with Function unknown (S)                       
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation") # heatmap with hierarchical clustering trees  
#I want to replace symbol KOG catgories with full class names                      
kog_map <- read.csv("./kog_classes.txt", sep = '\t', header=FALSE)                      
ktable2 <- ktable%>% rownames_to_column("V1")
for_map <- left_join(ktable2, kog_map, by = "V1")
rownames(ktable) == for_map$V1 #check they match                       
pheatmap(as.matrix(ktable),clustering_distance_cols="correlation", labels_row=for_map$V2) #with KOG class full names         
                       
### correlations plot
corrPlot(x="symb NO-LOW",y="symb NO-MODERATE",ktable)
corrPlot(x="host NO-LOW",y="host NO-MODERATE",ktable)                       
                       
                       

                       
                       
###########################################################################
#################### GO enrichment analyses (GO_MWU) ######################
###########################################################################
#following https://github.com/z0on/GO_MWU/blob/master/GO_MWU.R 
                       
host_annotations <- read.csv("./eggNOG_output_plutea_annotation_ALLgo.tsv", sep = '\t', header=TRUE)
host_annotations_GO = host_annotations[,c(1, 10)] #keep only gene and GO   
host_annotations_GO$GOs <- gsub(',', ';', host_annotations_GO$GOs) #for GO_MWU, string of concatenated GO terms separates by semicolons
host_annotations_GO$GOs <- gsub('-', 'unknown', host_annotations_GO$GOs) #for GO_MWU, the genes without annotation should be called "unknown", if you want to analyze these too
host_annotations_GO$query <- gsub (".m1", "", as.character(host_annotations_GO$query)) #editing gene names to match my reads      
host_annotations_GO$query <- gsub (".model.", ".TU.", as.character(host_annotations_GO$query)) #editing gene names to match my reads                
                       
                       
#################                     
###LogFC between no stress and moderate                      
                       
                       
