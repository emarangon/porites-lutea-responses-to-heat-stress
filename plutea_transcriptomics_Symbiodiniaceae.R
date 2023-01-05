################################################################################################################
################### plutea - heat stress experiment - Symbiodiniaceae transcriptomics analyses #################
################################################################################################################

##author EMMA MARANGON


##########################################################
#################### PRE-PROCESSING ######################
##########################################################

###loading libraries




###import files

files <- dir(path = "./13_mapping_plutea_UQ_cladocopium_goreaui_results", pattern = "*plutea_UQ_cladocopium_goreaui_relaxed_ReadsPerGene_StrandnessReverse.out.tab", full.names = T) 

cladocopium_counts1 <- files %>% purrr::map(read_tsv,  col_names = TRUE ) %>%
  purrr::reduce(cbind)

cladocopium_counts1 %>% head() #now I have one column per sample named gene_id (all identical across samples) and one column per sample with counts (sample_counts). 
#I rename only the first column (any gene_id column would be ok as all identical) and then remove all the 'gene-id' columns
names(cladocopium_counts1) <- gsub(x = names(cladocopium_counts1), pattern = "_", replacement = ".") #repace _ with .
colnames(cladocopium_counts1)[1] <- "gene_plutea_cladocopium_1234567890" #I rename first column
cladocopium_counts1 <- cladocopium_counts1 %>%dplyr::select(-ends_with("id")) #I remove all gene_id columns
names(cladocopium_counts1) <-  gsub('.{11}$', '', names(cladocopium_counts1)) #I remove last 11 characters from each colname
cladocopium_counts1 %>% head() #it worked


###filtering samples

remove <- "P4.31" # because of QC
cladocopium_counts2 <- cladocopium_counts1[, !(names(cladocopium_counts1) %in% remove)] #remove sample
remove <- "P3.33" # because of QC
cladocopium_counts3 <- cladocopium_counts2[, !(names(cladocopium_counts2) %in% remove)] #remove sample
remove <- "P3.18" # because identified as outlier (WGCNA clustering)
cladocopium_counts <- cladocopium_counts3[, !(names(cladocopium_counts3) %in% remove)] #remove sample
head (cladocopium_counts)
#I'll use 33 samples (36-3) for downstream analyses


###convert into DEGList object

cladocopium_samples <- read_tsv("./SampleInfo_plutea.txt") #metadata
cladocopium_samples1<- subset (cladocopium_samples, sample_id != "P4.31" & sample_id != "P3.33" & sample_id != "P3.18") #filtering out samples (see above)
#it is very important to make sure the metadata is in the same order as the column names of the counts table !!!
gn <- cladocopium_counts %>% dplyr::select (-"gene_plutea_cladocopium") #Just removing first gene column to be able to compare samples in the next step
table(colnames(gn) == cladocopium_samples1$sample_id) #check order corresponds -> it does!

cladocopium_counts_matrix <- cladocopium_counts %>%
  dplyr::select(-gene_plutea_cladocopium) %>%
  as.matrix()
row.names(cladocopium_counts_matrix) <- cladocopium_counts$gene_plutea_cladocopium #add gene_plutea_cladocopium
head (cladocopium_counts_matrix)

cladocopium_DGE2 <- DGEList(cladocopium_counts_matrix, samples = cladocopium_samples1) #I have a tot of 39,006 genes (some have zero counts tho)


###filtering genes
table(rowSums(cladocopium_DGE2$counts==0)==33) #how many genes have zero counts across all 33 samples? 9795
filt <- filterByExpr(cladocopium_DGE2, design = model.matrix(~TreatTime + genotype, 
                                                      data = cladocopium_DGE2$samples), min.count = 20) 
cladocopium_filtered1 <- cladocopium_DGE2[filt, , keep.lib.sizes = F]
dim(cladocopium_filtered1) #23,345 genes, 33 samples
table(rowSums(cladocopium_filtered1$counts==0)==33) #how many genes have zero counts across all 33 samples? zero

#To check filter cut off before and after
par(mfrow = c(2, 2)) # Create a 2 x 2 plotting matrix
mean_log_cpm3 <- aveLogCPM(cladocopium_DGE2$counts)
filter_threshold <- 0.5
hist(mean_log_cpm3)
abline(v = filter_threshold)
qqnorm(mean_log_cpm3)
abline(h = filter_threshold)
# 
mean_log_cpm4 <- aveLogCPM(cladocopium_filtered1$counts)
filter_threshold <- 0.5
hist(mean_log_cpm4)
abline(v = filter_threshold)
qqnorm(mean_log_cpm4)
abline(h = filter_threshold)
#improved a lot!
par(mfrow = c(1, 1)) #back to standard visualization


###normalization
cladocopium_filtered2 <- calcNormFactors(cladocopium_filtered1, method = "TMM") #calcNormFactors doesn’t normalize the data, it just calculates normalization factors for use downstream !!
cladocopium_filtered2$samples$norm.factors #all scaling factors are relatively close to 1 -> the effect of TMM-normalisation is mild for my data
table(cladocopium_samples1$sample_id == colnames(cladocopium_filtered2)) #double check data are still matching
cladocopium_filtered2_cpm <- cpm(cladocopium_filtered2, log=TRUE) #converting into log-cpm



#######################################################################
#################### DEGs analyses (limma-vooom) ######################
#######################################################################

treatment <- cladocopium_filtered2$samples$treatment
treatment <- factor(treatment)
time <- cladocopium_filtered2$samples$time
time <- factor(time)
group <- interaction(treatment, time)
group
group <- factor(group)
levels (group)
genotype <- cladocopium_filtered2$samples$genotype
genotype <- factor(genotype)
mm <- model.matrix(~0 + group)
v <- voom(cladocopium_filtered2, mm, plot = T)
cor <- duplicateCorrelation(v, mm, block=genotype)
cor$consensus.correlation #0.1692757
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
#9 down and 18 up in low relative to no stress
#1353 down and 1841 up in moderate relative to no stress
                       


############################################################################
#################### KOG enrichment analyses (KOGMWU) ######################
############################################################################
                       
cladocopium_annotations <- read.csv("./eggNOG_output_cladocopium_annotation_ALLgo.tsv", sep = '\t', header=TRUE) #annotation reference genome using EggNOG mapper
cladocopium_annotations_KOG = cladocopium_annotations[,c(1, 7)] #keep only gene and KOG                    
cladocopium_annotations_KOG$query <- gsub (".mRNA1", "", as.character(cladocopium_annotations_KOG$query)) #editing gene names to match my reads
cladocopium_annotations_KOG_noNA <- cladocopium_annotations_KOG[!grepl("-", cladocopium_annotations_KOG$COG_category),] #removing the unassigned

#I want to have one KOG per gene per row
cladocopium_annotations_KOG_noNA_multiple <- cladocopium_annotations_KOG_noNA %>% separate (COG_category, c("kog","B", "C", "D", "E"), sep=cumsum(c(1,1,1,1,1)))
cladocopium_annotations_KOG_noNA_multiple <- cladocopium_annotations_KOG_noNA_multiple %>% mutate_all(na_if,"")
first = cladocopium_annotations_KOG_noNA_multiple[,c(1,2)]
second = cladocopium_annotations_KOG_noNA_multiple[,c(1,3)]
second <- na.omit(second) #removing NA lines
third = cladocopium_annotations_KOG_noNA_multiple[,c(1,4)]
third <- na.omit(third) #removing NA lines
fourth = cladocopium_annotations_KOG_noNA_multiple[,c(1,5)]
fourth <- na.omit(fourth) #removing NA lines
fifth = cladocopium_annotations_KOG_noNA_multiple[,c(1,6)]
fifth <- na.omit(fifth) #removing NA lines
colnames(second) <- c("query", "kog") #need to be same colname for dplyr::bind_rows step
colnames(third) <- c("query", "kog")
colnames(fourth) <- c("query", "kog")
colnames(fifth) <- c("query", "kog")
combined12 <- dplyr::bind_rows(first, second) #multiple rows per gene if >1 KOG
combined123 <- dplyr::bind_rows(combined12, third)
combined1234 <- dplyr::bind_rows(combined123, fourth)
cladocopium_gene2kog <- dplyr::bind_rows(combined1234, fifth) #Two-column dataframe of gene annotations: gene id, KOG class
head(cladocopium_gene2kog)                       

save(cladocopium_filtered2, file = "./cladocopium_filtered2.R.data")
save(cladocopium_gene2kog, file = "./cladocopium_gene2kog.R.data")   
                       
#CONTINUATION KOG ANALYSES IS IN plutea_transcriptomics_host.R
                       
 
                       
###########################################################################
#################### GO enrichment analyses (GO_MWU) ######################
###########################################################################
#following https://github.com/z0on/GO_MWU/blob/master/GO_MWU.R                       
                       
cladocopium_annotations <- read.csv("./eggNOG_output_cladocopium_annotation_ALLgo.tsv", sep = '\t', header=TRUE)
cladocopium_annotations_GO = cladocopium_annotations[,c(1, 10)] #keep only gene and GO        
cladocopium_annotations_GO$GOs <- gsub(',', ';', cladocopium_annotations_GO$GOs) #for GO_MWU, string of concatenated GO terms separates by semicolons                       
cladocopium_annotations_GO$GOs <- gsub('-', 'unknown', cladocopium_annotations_GO$GOs) #for GO_MWU, the genes without annotation should be called "unknown", if you want to analyze these too
cladocopium_annotations_GO$query <- gsub (".mRNA1", "", as.character(cladocopium_annotations_GO$query)) #editing gene names to match my reads      

#################                     
###LogFC between no stress and moderate
treatment <- cladocopium_filtered2$samples$treatment
treatment <- factor(treatment)
time <- cladocopium_filtered2$samples$time
time <- factor(time)
group <- interaction(treatment, time)
group
group <- factor(group)
genotype <- cladocopium_filtered2$samples$genotype
genotype <- factor(genotype)
mm <- model.matrix(~0 + group)
v <- voom(cladocopium_filtered2, mm, plot = T)
cor <- duplicateCorrelation(v, mm, block=genotype)
fit <- lmFit(object=v, design=mm, 
             block=genotype, correlation=cor$consensus.correlation)
head(coef(fit))
contr <- makeContrasts(groupHeat.T4 - (groupHeat.T0 + groupAmbient.T0 + groupAmbient.T2 + groupAmbient.T4)/4,
                       levels = colnames(coef(fit)))
contr
fit.cont <- contrasts.fit(fit, contr)
e.fit.cont <- eBayes(fit.cont)
top.table <- topTable(e.fit.cont, sort.by = "P", n = Inf)
keep <- tibble::rownames_to_column(top.table, "genes")
cladocopium_genes_HT4vsAll_block_logF = keep[,c(1,2)] #keep only first and 2nd column
#save as 'cladocopium_genes_HT4vsAll_block_logFC.csv'                        

                       
################# BP              
input="cladocopium_genes_HT4vsAll_block_logFC.csv"
goAnnotations="cladocopium_annotations_genome_GO.tab"
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
write.csv(BP_results, "BP_graph_NOvsMODRATE.csv")

                       
                       
################# MF              
input="cladocopium_genes_HT4vsAll_block_logFC.csv"
goAnnotations="cladocopium_annotations_genome_GO.tab"
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
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)

# text representation of results, with actual adjusted p-values
results[[1]]
MF_results <- results[[1]]
write.csv(MF_results, "MF_graph_NOvsMODRATE.csv")                       
                       

                       
#################                     
###LogFC between no stress and low

treatment <- cladocopium_filtered2$samples$treatment
treatment <- factor(treatment)
time <- cladocopium_filtered2$samples$time
time <- factor(time)
group <- interaction(treatment, time)
group
group <- factor(group)
genotype <- cladocopium_filtered2$samples$genotype
genotype <- factor(genotype)
mm <- model.matrix(~0 + group)
v <- voom(cladocopium_filtered2, mm, plot = T)
cor <- duplicateCorrelation(v, mm, block=genotype)
fit <- lmFit(object=v, design=mm, 
             block=genotype, correlation=cor$consensus.correlation)
head(coef(fit))
contr <- makeContrasts(groupHeat.T2 - (groupHeat.T0 + groupAmbient.T0 + groupAmbient.T2)/3,
                       levels = colnames(coef(fit)))
contr
fit.cont <- contrasts.fit(fit, contr)
e.fit.cont <- eBayes(fit.cont)
top.table <- topTable(e.fit.cont, sort.by = "P", n = Inf)
keep <- tibble::rownames_to_column(top.table, "genes") 
cladocopium_genes_HT2vsAll_block_logF = keep[,c(1,2)] #keep only first and 2nd column
head (cladocopium_genes_HT2vsAll_block_logF)
#save as 'cladocopium_genes_HT2vsAll_block_logFC.csv' 
  
                       
################# BP
input="cladocopium_genes_HT2vsAll_block_logFC.csv"
goAnnotations="cladocopium_annotations_genome_GO.tab"
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
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)

# text representation of results, with actual adjusted p-values
results[[1]]
BP_results <- results[[1]]
write.csv(BP_results, "BP_graph_NOvsLOW.csv")     
                       
                       
################# MF
input="cladocopium_genes_HT2vsAll_block_logFC.csv"
goAnnotations="cladocopium_annotations_genome_GO.tab"
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
                  colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)

# text representation of results, with actual adjusted p-values
results[[1]]
MF_results <- results[[1]]
write.csv(MF_results, "MF_graph_NOvsLOW.csv")                       
                       
                       

                       
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
meta_1 <- cladocopium_filtered1$samples
meta_2<-  dplyr::left_join(meta_1, pr_porites2, by=c('sample_id'='sample_id', 'time'='TimePoint')) #keep only samples used also for transcriptomics
meta_3 <-  dplyr::left_join (meta_2, bleaching_porites2, by=c('sample_id'='sample_id', 'time'='TimePoint'))
meta_df <-  dplyr::left_join (meta_3, photefficiency_porites2, by=c('sample_id'='sample_id', 'time'='TimePoint'))
head(meta_df)                       
                       
                       
###normalization for mixOmics
                       
cladocopium_for_deseq <- cladocopium_filtered1$counts
cladocopium_de_input = as.matrix(cladocopium_for_deseq)
meta_df$stress <- factor (meta_df$stress)
levels(meta_df$stress)
all.equal(colnames(cladocopium_de_input), meta_df$sample_id) #check they match
cladocopium_dds <- DESeqDataSetFromMatrix(round(cladocopium_de_input), #to convert values to integers, returns a matrix
                                          meta_df,
                                          design = ~1) #I am not specyfing a model here
class(cladocopium_dds) #DESeqDataSet object
cladocopium_dds_norm <- vst(cladocopium_dds)

X_cladocopium <- assay (cladocopium_dds_norm) %>% t
dim(X_cladocopium)
X <- X_cladocopium                      

 
###PCA
                       
tune_pca<- tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
plot(tune_pca) #4 components seem to explain most of the variance from the data
pca.result <- pca(X, ncomp = 4, center = TRUE, scale = FALSE)
pca.result #results
Colors <- c(
  "#ECBE92", "#F7794D", "#B3B4B4")
plotIndiv(pca.result, 
          group = meta_df$stress, legend=TRUE, legend.title = 'Heat stress', 
          pch = as.factor(meta_df$genotype), legend.title.pch = 'Parental colony',
          title = 'PCA plutea cldocopium sample plot', style = 'ggplot2',
          ellipse = FALSE, col.per.group = Colors)                        
     
                       
###multilevel PCA (to account for effect of parental colony)
  
design <- data.frame(sample = meta_df$genotype) #repeated measures are accounted for (parental colony in my case)
pca.multilevel.result <- pca(X, ncomp = 4, scale = FALSE, center = TRUE, 
                             multilevel = design)
pca.multilevel.result
plotIndiv(pca.multilevel.result, 
          group = meta_df$stress, legend=TRUE, legend.title = 'Heat stress', 
          pch = as.factor(meta_df$genotype), legend.title.pch = 'Parental colony',
          title = 'MULTILEVEL PCA plutea cladocopium sample plot', col.per.group = Colors) #plot the samples                      
                       
                       
###multilevel sPlS-DA (following http://mixomics.org/case-studies/multilevel-vac18-case-study/)  
                       
dim(X)
Y <- as.factor(meta_df$stress)
summary(Y) 
cladocopium.splsda.stress <- splsda(X, Y, ncomp = 10)  # set ncomp to 10 for performance       
perf.cladocopium.splsda.stress <- perf(cladocopium.splsda.stress, validation = "Mfold", 
  folds = 5, nrepeat = 100, 
  progressBar = FALSE, auc = TRUE)
plot(perf.cladocopium.splsda.stress, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")                       
perf.cladocopium.splsda.stress$choice.ncomp #1 component                       
perf.cladocopium.splsda.stress$error.rate.class$mahalanobis.dist #BER around 0.2 for 1 comp
list.keepX <- c(1:10,  seq(100, 150, 3))  #chose this grad after multiple trials starting with 1-300
tune.splsda.cladocopium.stress <- tune.splsda(X, Y, ncomp = 2, # calculate for 2 components (1 is optimal - see above - but still trying with 2 here - following mixomics book)
   validation = 'Mfold',
  folds = 5, nrepeat = 100, # use repeated cross-validation
  dist = 'mahalanobis.dist', # use mahalanobis.dist measure
 measure = "BER", # use  balanced error rate of dist measure instead of overall because I have unbalanced designed (suggested in mixomics book)
  test.keepX = list.keepX #number of variables to select on each component
)                       
head (tune.splsda.cladocopium.stress$error.rate)
head (tune.splsda.cladocopium.stress$error.rate.class)
plot(tune.splsda.cladocopium.stress, col = color.jet(2)) #error rate around 0.25               
tune.splsda.cladocopium.stress$choice.ncomp$ncomp #1
tune.splsda.cladocopium.stress$choice.keepX #142
final.multilevel.splsda.cladocopium.stress <- splsda(X, Y, ncomp = 2, #I choose 2 comp instead of optimal (1) only for visuaization purposes
                                              keepX = c(142, 142), #kept same number as per first component as suggeated in book when optimal comp =1
                                              multilevel = design)                       
Colors <- c(
  "#ECBE92", "#F7794D", "#B3B4B4")
plotIndiv(final.multilevel.splsda.cladocopium.stress, group = meta_df$stress, 
          ind.names = FALSE, 
          legend = TRUE, legend.title = 'Heat stress',
          pch = as.factor(meta_df$genotype), legend.title.pch = 'Parental colony',
          ellipse = TRUE, col.per.group = Colors,
          title = 'Sample Plot of multilevel sPLS-DA on plutea cladocopium data')                       
 pdf(file = "cimCladocopium_splsda_stress.pdf", height = 50, width = 50)
cim.final.multilevel.splsda.cladocopium.stress <- cim(final.multilevel.splsda.cladocopium.stress,margins = c(10, 6), 
                                                   row.names = paste(meta_df$DHW, meta_df$sample_id, sep = "_"),
                                                   row.sideColors = color.mixo(Y),
                                                   col.names = TRUE, legend=list(legend = levels(Y), 
                                                                                 title = "DHW", cex = 0.8),
                                                   comp=1) #ncomp=1 is the optimal
dev.off()                         
                       
