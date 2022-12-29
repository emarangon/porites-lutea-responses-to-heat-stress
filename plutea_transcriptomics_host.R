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
                       
