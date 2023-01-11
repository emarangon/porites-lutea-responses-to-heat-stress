#####################################################################################################
################### plutea - heat stress experiment - host transcriptomics analyses #################
#####################################################################################################

##author EMMA MARANGON


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
contr <- makeContrasts(groupHeat.T4 - (groupHeat.T0 + groupAmbient.T0 + groupAmbient.T2 + groupAmbient.T4)/4,
                       levels = colnames(coef(fit)))
fit.cont <- contrasts.fit(fit, contr)
e.fit.cont <- eBayes(fit.cont)
top.table <- topTable(e.fit.cont, sort.by = "P", n = Inf)
keep <- tibble::rownames_to_column(top.table, "genes")
host_only_degs_HT4vsAll_block_logF = keep[,c(1,2)]
#save host_only_degs_HT4vsAll_block_logF as 'host_genes_HT4vsAll_block_logFC.csv' 
     
my_GO_host <- left_join (host_only_degs_HT4vsAll_block_logF, host_annotations_GO, by=c("genes"="query"))
my_GO_host$GOs[is.na(my_GO_host$GOs)] <- "unknown"
dim(my_GO_host)
my_GO_host = my_GO_host[,c(1, 3)]
head (my_GO_host)
#save my_GO_host as 'host_annotations_genome_GO.tab'
                       
                       
################# BP                     
input="host_genes_HT4vsAll_block_logFC.csv"
goAnnotations="host_annotations_genome_GO.tab"
goDatabase="gene_ontology_OBO_June22.obo"
goDivision="BP" #biological processes
source("gomwu.functions.R")                        
                       
# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # perl path usr/bin/perl
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.
                       
# ----------- Plotting results

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  #absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.01, # FDR cutoff to print in regular (not italic) font.
                  level3=0.001, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") 
)  

# text representation of results, with actual adjusted p-values
results[[1]]
BP_results <- results[[1]]
write.csv(BP_results, "BP_results_NOvsMODRATE.csv")
                       
                       
################# MF                         
input="host_genes_HT4vsAll_block_logFC.csv"
goAnnotations="host_annotations_genome_GO.tab"
goDatabase="gene_ontology_OBO_June22.obo"
goDivision="MF" #Molecular Function
source("gomwu.functions.R")
                       
# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # perl path usr/bin/perl
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# ----------- Plotting results

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  #absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=1, # un-remark this if you are using log2-fold changes 
                  level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.01, # FDR cutoff to print in regular (not italic) font.
                  level3=0.001, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") 
)

# text representation of results, with actual adjusted p-values
results[[1]]
MF_results <- results[[1]]
write.csv(MF_results, "MF_results_NOvsMODRATE.csv")
                       

                       
#################                     
###LogFC between no stress and low
                       
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
contr <- makeContrasts(groupHeat.T2 - (groupHeat.T0 + groupAmbient.T0 + groupAmbient.T2)/3,
                       levels = colnames(coef(fit)))
fit.cont <- contrasts.fit(fit, contr)
e.fit.cont <- eBayes(fit.cont)
top.table <- topTable(e.fit.cont, sort.by = "P", n = Inf)
keep <- tibble::rownames_to_column(top.table, "genes")
host_only_degs_HT2vsAll_block_logF = keep[,c(1,2)]                    
#save host_only_degs_HT2vsAll_block_logF as 'host_genes_HT2vsAll_block_logFC.csv'                    
                       
my_GO_host <- left_join (host_only_degs_HT2vsAll_block_logF, host_annotations_GO, by=c("genes"="query"))
my_GO_host$GOs[is.na(my_GO_host$GOs)] <- "unknown"
dim(my_GO_host)
my_GO_host = my_GO_host[,c(1, 3)]
head (my_GO_host)
#save my_GO_host as 'host_annotations_genome_GO.tab'
                       
                       
################# BP
input="host_genes_HT2vsAll_block_logFC.csv"
goAnnotations="host_annotations_genome_GO.tab"
goDatabase="gene_ontology_OBO_June22.obo"
goDivision="BP" #biological processes
source("gomwu.functions.R")
                       
# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # perl path usr/bin/perl
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# ----------- Plotting results

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  #absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.01, # FDR cutoff to print in regular (not italic) font.
                  level3=0.001, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)

# text representation of results, with actual adjusted p-values
results[[1]]
BP_results <- results[[1]]
write.csv(BP_results, "BP_results_graph_NOvsLOW.csv")
                       
                       
################# MF
input="host_genes_HT2vsAll_block_logFC.csv"
goAnnotations="host_annotations_genome_GO.tab"
goDatabase="gene_ontology_OBO_June22.obo"
goDivision="MF" #Molecular Function
source("gomwu.functions.R")
                       
# ------------- Calculating stats
# It might take a few minutes for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.

gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="perl", # perl path usr/bin/perl
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=5,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #	Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #	Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.

# ----------- Plotting results

quartz()
results=gomwuPlot(input,goAnnotations,goDivision,
                  #absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.01, # FDR cutoff to print in regular (not italic) font.
                  level3=0.001, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral")
)

# text representation of results, with actual adjusted p-values
results[[1]]
MF_results <- results[[1]]
write.csv(MF_results, "MF_results_graph_NOvsLOW.csv")
                       
                
                       
                       
###########################################################################
#################### mixOmics analyses ######################
###########################################################################
  
###importing physiology measures
                       
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
                       
                       
###normalization for mixOmics                       
                       
host_for_deseq <- host_filtered1$counts
host_de_input = as.matrix(host_for_deseq)
meta_df$stress <- factor (meta_df$stress)
levels(meta_df$stress)     
all.equal(colnames(host_de_input), meta_df$sample_id) #they match
host_dds <- DESeqDataSetFromMatrix(round(host_de_input), #to convert values to integers, returns a matrix
                              meta_df,
                              design = ~1) #I am not specyfing a model here
class(host_dds) #DESeqDataSet object
host_dds_norm <- vst(host_dds) #normalize and transform
X_host <- assay (host_dds_norm) %>% t
dim(X_host)
X <- X_host                       
                       
                       
###PCA
                       
tune_pca<- tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)                       
plot(tune_pca) #3 components seem to explain most of the variance from the data
pca.result <- pca(X, ncomp = 3, center = TRUE, scale = FALSE)
pca.result #results
Colors <- c(
  "#ECBE92", "#F7794D", "#B3B4B4")
plotIndiv(pca.result, 
          group = meta_df$stress, legend=TRUE, legend.title = 'Heat stress', 
          pch = as.factor(meta_df$genotype), legend.title.pch = 'Parental colony',
          title = 'PCA plutea host sample plot', style = 'ggplot2',
          ellipse = FALSE, col.per.group = Colors)                        
       
                       
###multilevel PCA (to account for effect of parental colony)

design <- data.frame(genotype = meta_df$genotype) #repeated measures are accounted for (parental colony in my case)
pca.multilevel.result <- pca(X, ncomp = 3, scale = FALSE, center = TRUE, 
                             multilevel = design)
pca.multilevel.result
plotIndiv(pca.multilevel.result, 
          group = meta_df$stress, legend=TRUE, legend.title = 'Heat stress', 
          pch = as.factor(meta_df$genotype), legend.title.pch = 'Parental colony',
          title = 'MULTILEVEL PCA plutea host sample plot', col.per.group = Colors)                       
                       
                       
###multilevel sPlS-DA (following http://mixomics.org/case-studies/multilevel-vac18-case-study/)  

dim(X)
Y <- as.factor(meta_df$stress)
summary(Y)  
host.splsda.stress <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance                       
perf.splsda.host.stress <- perf(host.splsda.stress, validation = "Mfold", 
                             folds = 5, nrepeat = 100,
                             progressBar = TRUE, auc = TRUE) 
plot(perf.splsda.host.stress, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")                      
perf.splsda.host.stress$choice.ncomp #6 components
perf.splsda.host.stress$error.rate.class$mahalanobis.dist  #BER around 0.45 for 6 comp                      
list.keepX <- c(1:10,  seq(10, 200, 10)) #I chose this grid after multiple trials starting with 1-300                      
tune.splsda.host.stress <- tune.splsda(X, Y, ncomp=6, # calculate for first 6 components (based on optimal see above)
                                    validation = 'Mfold',
                                    folds = 5, nrepeat = 100, # use repeated cross-validation
                                    dist = 'mahalanobis.dist', # use mahalanobis.dist measure
                                    measure = "BER", #use  balanced error rate of dist measure instead of overall because I have unbalanced designed (suggested in mixomics book)
                                    test.keepX = list.keepX, #number of variables to select on each component
                                    progressBar = TRUE)                       
head (tune.splsda.host.stress$error.rate)
head (tune.splsda.host.stress$error.rate.class) 
plot(tune.splsda.host.stress, col = color.jet(6)) #error rate around 0.45 for 1 comp!                    
tune.splsda.host.stress$choice.ncomp$ncomp #1 comp
tune.splsda.host.stress$choice.keepX #4
optimal.keepX <- tune.splsda.host.stress$choice.keepX[1:optimal.ncomp]
final.multilevel.splsda.host.stress <- splsda(X, Y, ncomp = 2,  #I choose 2 comp instead of optimal (1) only for visuaization purposes
                                              keepX = optimal.keepX,
                                           multilevel = design)                      
Colors <- c(
  "#ECBE92", "#F7794D", "#B3B4B4")
plotIndiv(final.multilevel.splsda.host.stress, group = meta_df$stress, 
          comp = c(1,2),
          ind.names = FALSE,
          legend = TRUE, legend.title = 'Heat stress',
          ellipse = TRUE, col.per.group = Colors, 
          title = 'Sample Plot of multilevel sPLS-DA on host data, comp 1 & 2',
          pch = as.factor(meta_df$genotype), legend.title.pch = 'Parental colony')   
pdf(file = "multilevel.splsda.host.stress.pdf", height = 10, width = 10)
cim.final.multilevel.splsda.host.stress <- cim(final.multilevel.splsda.host.stress, margins = c(10, 10), 
                                            row.names = paste(meta_df$TreatTime, meta_df$sample_id, sep = "_"),
                                            row.sideColors = color.mixo(Y),
                                            col.names = TRUE, legend=list(legend = levels(Y), 
                                                                          title = "Heat stress", cex = 0.8), 
                                            comp=1) #optimal component number
dev.off()                       
                       
                       
                       
