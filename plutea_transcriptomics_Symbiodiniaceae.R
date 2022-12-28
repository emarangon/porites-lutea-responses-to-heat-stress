################################################################################################################
################### plutea - heat stress experiment - Symbiodiniaceae transcriptomics analyses #################
################################################################################################################

##EMMA MARANGON


##########################################################
#################### PRE-PROCESSING ######################
##########################################################

###loading libraries




###import files

files <- dir(path = "13_mapping_plutea_UQ_cladocopium_goreaui_results", pattern = "*plutea_UQ_cladocopium_goreaui_relaxed_ReadsPerGene_StrandnessReverse.out.tab", full.names = T) 

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
remove <- "P3.18" # because of outlier (WGCNA)
cladocopium_counts <- cladocopium_counts3[, !(names(cladocopium_counts3) %in% remove)] #remove sample
head (cladocopium_counts)
#I'll use 33 samples (36-3) for downstream analyses


###convert into DEGList object

cladocopium_samples <- read_tsv("SampleInfo_plutea.txt") #metadata
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



