####################################################################################
################### plutea - heat stress experiment - 16S analyses #################
####################################################################################

##EMMA MARANGON


##########################################################
#################### PRE-PROCESSING ######################
##########################################################

#loading libraries

library(phyloseq) #microbial analyses
library(tidyverse) #includes ggplot2, dplyr, tibble, tidyR
library(RColorBrewer) #to change colours graphs
library(scales) #to convert into percentage on graph
library (decontam) #identify contaminants
library(vegan) #adonis
library(RVAideMemoire) #post hocs adonis
library (DESeq2) #deseq analyses
library(glmmTMB) #alpha diversity stats
library(emmeans) #alpha diversity stats
library(DHARMa) #checking model assumptions for stats alpha diversity
library(forcats)



############################################################
### load data ####

##1 OTU TABLE (actually ASV)
otu_table = read.csv("heat-feature-table_R.txt", sep = '\t',  dec = ".", check.names = FALSE, row.names=1) #check.names = FALSE to avoid having X be automatically placed in front of each sampleID
head (otu_table)
str (otu_table)

##2 TAXA TABLE
otu_matrix = read.csv("heat-taxonomy_R.tsv", sep = '\t', header=T, row.names=1) #I have separated each taxa level in the csv
head (otu_matrix)
# I split column Taxon in 7 columns and discard the extra 7 columns (Domain,Phylum,Class,Order,Family,Genus,Species) with no information; I keep the Taxon column
TAXA_TABLE_split<-otu_matrix %>% separate(Taxon, c("Domain","Phylum","Class","Order","Family","Genus","Species"), sep=";", remove=TRUE) #I want to remove the old taxonomic column i.e. remove=TRUE) otherwise this column is not going to make work properly tax_glom and potentially other phyloseq functions
TAXA_TABLE_matrix<-as.matrix(TAXA_TABLE_split)
head (TAXA_TABLE_matrix)
str (TAXA_TABLE_matrix)

##3 METADATA
metadata = read.csv("heat-metadata.txt", sep = '\t', row.names=1) #row.names=1 -> R takes the first column of your dataset and uses it as the rownames of your data frame. 
head (metadata)
metadata1 <- metadata %>% unite ("time.genotype", timePoint, sampleType, genotype, sep=".", remove=FALSE, na.rm = TRUE) # I add a column as combination of time and genotype to the metadata
head (metadata1)
metadata2 <- metadata1 %>% unite ("TreatTime", treatment, timePoint, sep=".", remove=FALSE, na.rm = TRUE) %>%
  unite ("time.genotype.treat", treatment, timePoint, sampleType, genotype, sep=".", remove=FALSE, na.rm = TRUE) 
head (metadata2)

### import ###
OTU = otu_table(otu_table, taxa_are_rows = TRUE) #table is oriented with taxa as rows
TAX = tax_table (TAXA_TABLE_matrix)
meta = sample_data(metadata2)
phy_tree = read_tree("heat-tree.nwk") #rooted tree

### merging ###
phyloseq_merged = phyloseq(OTU, TAX, meta, phy_tree)
phyloseq_merged


############################################################
### filtering ####

### filtering eukaryota, mitochondria, chloroplasts ###
phyloseq_merged_filtered <- phyloseq_merged %>% subset_taxa(Domain != "D_0__Eukaryota" & Family  != "D_4__Mitochondria" & Order   != "D_3__Chloroplast")

### filtering out samples not belonging to the plutea experiment
phyloseq_merged_filtered_Porites <- phyloseq_merged_filtered %>% subset_samples(experiment == "Porites")

### data info (including blanks, negatives, contminants) ###
sum(sample_sums(phyloseq_merged_filtered_Porites)) # 3157960 reads
summary(sample_sums(phyloseq_merged_filtered_Porites)) # Between 7 and 84813 reads per sample, mean = 16887
phyloseq_merged_filtered_Porites # 22613 ASVs for 187 samples


############################################################
### identifying contaminants using decontam ###

##### 1. inspect library size
#following tutorial https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
head(sample_data(phyloseq_merged_filtered_Porites))
df <- as.data.frame(sample_data(phyloseq_merged_filtered_Porites))
df$LibrarySize <- sample_sums(phyloseq_merged_filtered_Porites)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=sampleType)) + geom_point()

##### 2. identify contaminants - prevalence
sample_data(phyloseq_merged_filtered_Porites)$is.neg <- sample_data(phyloseq_merged_filtered_Porites)$sampleType == "blank"
contamdf.prev_0.5 <- isContaminant(phyloseq_merged_filtered_Porites, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev_0.5$contaminant) #FALSE 22612 #TRUE 1 -> 1 ASVs has been identified as contaminants
# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(phyloseq_merged_filtered_Porites, function(abund) {1*(abund>0)})
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sampleType == "blank", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sampleType != "blank", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Blanks)") + ylab("Prevalence (All Samples but blanks)")

#### 3. removing contaminants
#based on https://github.com/joey711/phyloseq/issues/652
phyloseq_merged_filtered_Porites_NoBlanks <- subset_samples(phyloseq_merged_filtered_Porites, sampleType!="blank")
phyloseq_merged_filtered_Porites_NoBlanksNonegatives <- subset_samples(phyloseq_merged_filtered_Porites_NoBlanks, sampleType!="negative") #no blanks no negatives
taxa_are_rows(phyloseq_merged_filtered_Porites_NoBlanksNonegatives) #TRUE
contaminants <- subset(contamdf.prev_0.5, contaminant == "TRUE") #I am using the strict threshold
rownames(contaminants)
# fd98346394a5c79e554003012cb33826 is the contaminant
contaminants_only <- c("fd98346394a5c79e554003012cb33826")
AllTaxa = taxa_names(phyloseq_merged_filtered_Porites_NoBlanksNonegatives)
AllTaxa <- AllTaxa[!(AllTaxa %in% contaminants_only)]
phyloseq_merged_filtered_Porites_NoBNC = prune_taxa(AllTaxa, phyloseq_merged_filtered_Porites_NoBlanksNonegatives) #no blanks no negatives no contaminants
phyloseq_merged_filtered_Porites_NoBNC
                                 
### data info (excluding blanks, negatives, contaminants) ###
summary(sample_sums(phyloseq_merged_filtered_Porites_NoBNC)) # Between 1117 and 84479 reads per sample, mean = 17540

                                 
############################################################
### more filtering ### 

### I remove one genotype as it has been identified as Porites lobata ###
phyloseq_merged_filtered_Porites_NoBNC_noP2 <- subset_samples(phyloseq_merged_filtered_Porites_NoBNC, genotype != "P2")
sum(sample_sums(phyloseq_merged_filtered_Porites_NoBNC_noP2)) # 2738889 reads
summary(sample_sums(phyloseq_merged_filtered_Porites_NoBNC_noP2)) # Between 2435 and 84479 reads per sample, mean = 18889
 
### I remove low abundance ASVs based on nseq (<5400 reads) ###                                 
phyloseq_merged_porites_noP2_filtered <- prune_samples(sample_sums(phyloseq_merged_filtered_Porites_NoBNC_noP2) > 5400, phyloseq_merged_filtered_Porites_NoBNC_noP2) 
#4 samples were removed (I have now 141 samples)
summary(sample_sums(phyloseq_merged_porites_noP2_filtered)) # Between 5411 and 84479 reads per sample, mean = 19315                                
                                 
### I remove singletons ### 
phyloseq_merged_porites_noP2_FINAL <- prune_taxa(taxa_sums(phyloseq_merged_porites_noP2_filtered) >=2, phyloseq_merged_porites_noP2_filtered)
phyloseq_merged_porites_noP2_FINAL = prune_taxa(taxa_sums(phyloseq_merged_porites_noP2_FINAL) > 0, phyloseq_merged_porites_noP2_FINAL) #remove unobserved ASVs (sum 0 across all samples)
phyloseq_merged_porites_noP2_FINAL #from 22612 ASVs to 12558                                 
                                 
### data info (FINAL DATASET) ###
phyloseq_merged_porites_noP2_FINAL # 12558 ASVs
sum(sample_sums(phyloseq_merged_porites_noP2_FINAL)) # 2723455 reads 
summary(sample_sums(phyloseq_merged_porites_noP2_FINAL)) # Between 5411 and 84479 reads per sample, mean = 19315       
         
                                 
############################################################
### rarefaction curves ###        

#seawater                                 
phyloseq_merged_porites_noP2_FINAL_sw <- subset_samples(phyloseq_merged_porites_noP2_FINAL, sampleType == "seawater")
rarecurve(t(otu_table(phyloseq_merged_porites_noP2_FINAL_sw)), step=50, cex=0.5, label = FALSE, xlab="Sequencing depth", col = "#00BFC3")

#coral
phyloseq_merged_porites_noP2_FINAL_porites <- subset_samples(phyloseq_merged_porites_noP2_FINAL, sampleType == "porites")
rarecurve(t(otu_table(phyloseq_merged_porites_noP2_FINAL_porites)), step=50, cex=0.5, label = FALSE, xlab="Sequencing depth", col = "#F8766C")

#feed                                 
phyloseq_merged_porites_noP2_FINAL_food <- subset_samples(phyloseq_merged_porites_noP2_FINAL, sampleType == "food")
rarecurve(t(otu_table(phyloseq_merged_porites_noP2_FINAL_food)), step=50, cex=0.5, label = FALSE, xlab="Sequencing depth", col = "#7CAD00")

                                 
############################################################
### normalization ###   
                                 
#proportions
phyloseq_merged_porites_noP2_FINALabundances = transform_sample_counts(phyloseq_merged_porites_noP2_FINAL, function(x){x / sum(x)}) # The first argument to this function is the phyloseq object you want to transform, and the second argument is an R function that defines the transformation (the counts of each sample will be transformed INDIVIDUALLY). 
table(otu_table(phyloseq_merged_porites_noP2_FINALabundances)) #just to check it worked 
phyloseq_merged_porites_noP2_FINALabundances_cutoff <- filter_taxa(phyloseq_merged_porites_noP2_FINALabundances, function(x) mean(x) > 1e-5, TRUE) #keeping only ASvs with mean > 0.00001; from 12558 to 2758 taxa (keeping 22% only) #https://joey711.github.io/phyloseq/preprocess.html
phyloseq_merged_porites_noP2_FINALabundances_cutoff = prune_taxa(taxa_sums(phyloseq_merged_porites_noP2_FINALabundances_cutoff) > 0, phyloseq_merged_porites_noP2_FINALabundances_cutoff) #to make sure no ASVs with sum 0 across all samples 
phyloseq_merged_porites_noP2_FINALabundances_cutoff #2758 ASVs
                                 
#rarefaction                                 
phyloseq_merged_porites_noP2_FINAL_noZero <- prune_taxa(taxa_sums(phyloseq_merged_porites_noP2_FINAL) > 0, phyloseq_merged_porites_noP2_FINAL) #prune ASVs that are not present in at least one sample 
phyloseq_merged_porites_noP2_FINALrarefied = rarefy_even_depth(phyloseq_merged_porites_noP2_FINAL_noZero, rngseed=1, sample.size = min(sample_sums(phyloseq_merged_porites_noP2_FINAL_noZero))) # RAREFYING; `set.seed(1)` was used to initialize repeatable random subsampling
phyloseq_merged_porites_noP2_FINALrarefied = prune_taxa(taxa_sums(phyloseq_merged_porites_noP2_FINALrarefied) > 0, phyloseq_merged_porites_noP2_FINALrarefied)  #remove unobserved ASVs (sum 0 across all samples)
                                                                   
      
                                                                   
                                                                   
                                 
###############################################################
#################### BETA DIVERSITY ANALYSES ##################
###############################################################
                                                                   
#all samples (coral, seawater, feed)                                                                   
All_noP2_sqrt <- transform_sample_counts(phyloseq_merged_porites_noP2_FINALabundances_cutoff, function (x) sqrt(x)) #sqrt transformation

#coral samples (NO seawater NO feed)   
Tissue_noP2 <- subset_samples(phyloseq_merged_porites_noP2_FINALabundances_cutoff, sampleType == "porites")
Tissue_noP2_sqrt <- transform_sample_counts(Tissue_noP2, function (x) sqrt(x)) #sqrt
   
                                            
                                            
#### NMDS of all samples (coral, seawater, feed)
                                                                  
Colors <- c(
  "#7CAC00", "#F8766B", "#00BFC2")
ordinate(All_noP2_sqrt, "NMDS", "bray") %>% 
  plot_ordination(All_noP2_sqrt, ., color = "sampleType", shape = "sampleType", title = "nmds_brays_porites") + 
  theme_bw() + geom_point(size = 3) + 
  scale_color_manual(values = Colors) + 
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

#stress
nmds_sampletype_stress <- ordinate(All_noP2_sqrt, "NMDS", "bray") #for next step
cat("Stress:", nmds_sampletype_stress$stress, fill=TRUE) #to measure Stress
#if stress <0.1 great; 0.1-0.2 good; 0.2-0.3, acceptable (treat with some caution), >0.3 unreliable
       
                                            
                                            
#### NMDS of CORAL samples (incuding baseline) - main figure
                                            
Colors <- c(
  "#B3B4B4", "#F3E0B5", "#EB9623", "#B84040", "#B3B4B4")
ordinate(Tissue_noP2_sqrt, "NMDS", "bray") %>% 
  plot_ordination(Tissue_noP2_sqrt, ., color = "DHW", shape = "treatment", title = "nmds_brays_porites") + 
  theme_bw() + geom_point(size = 4) + 
  scale_color_manual(values = Colors) + 
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14)
  ) 

#stress
nmds_tissue_stress <- ordinate(Tissue_noP2_sqrt, "NMDS", "bray") #for next step
cat("Stress:", nmds_tissue_stress$stress, fill=TRUE) #to measure Stress
                                            

                                            
#### NMDS of CORAL samples (incuding baseline; Parental colonies separated)
                                            
Colors <- c(
  "#B3B4B4", "#F3E0B5", "#EB9623", "#B84040", "#B3B4B4")
ordinate(Tissue_noP2_sqrt, "NMDS", "bray") %>% 
  plot_ordination(Tissue_noP2_sqrt, ., color = "DHW", shape = "treatment", title = "nmds_brays_porites_genotypes") + 
  theme_bw() + geom_point(size = 4) + 
  scale_color_manual(values = Colors) + 
 facet_wrap(~genotype) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14)
    )                                            
                                            
                                            
#### BARPLOT at the family level, coral samples only (baseline, ambient, heat) - main figure
Porites_Fam <- tax_glom(phyloseq_merged_porites_noP2_FINALabundances_cutoff, taxrank = 'Family') #agglomerate data (merging taxa of the same Family)
OnlyTissue <- subset_samples(Porites_Fam, sampleType == "porites") #only coral samples
Porites_fam_alll <- subset_samples(OnlyTissue, treatment == "heat" | treatment == "ambient" | treatment == "baseline")
porites_fam_all<- prune_taxa(taxa_sums(Porites_fam_alll) > 0, Porites_fam_alll) #to delete taxa which abundance sum is 0 
porites_Gfam_sample_all <- merge_samples(porites_fam_all, "TreatTime", fun=mean)  #I create a new group named Sample
porites_Gfam_sample_all_abund = transform_sample_counts(porites_Gfam_sample_all, function(x){x / sum(x)}) 
porites_Gfam_sample_all_abund_melt <- psmelt(porites_Gfam_sample_all_abund)
porites_Gfam_sample_all_abund_melt$Family <- as.character(porites_Gfam_sample_all_abund_melt$Family) # convert Family to a character vector from a factor
porites_Gfam_sample_all_abund_melt$Family[porites_Gfam_sample_all_abund_melt$Abundance < 0.05] <- "Other" #rename classes with < 5% abundance
porites_Gfam_sample_all_abund_melt2 <- porites_Gfam_sample_all_abund_melt %>% 
  separate(Sample, c("treat", "timepoint"), remove=FALSE) #re-creating treatment variable to use in the bar plot
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(17, "Set3"))(nb.cols)
p3 <- ggplot(data=porites_Gfam_sample_all_abund_melt2, aes(x=Sample, y=Abundance, fill=Family))
p3 + geom_bar(aes(), stat="identity", position="stack") +  
  guides(fill=guide_legend(nrow=5)) + labs(title = "porites_Gfam_abund ") +
  scale_fill_manual(values = mycolors) +
 theme_bw() +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )  +
  theme(legend.position="bottom", axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  labs (y = "Mean relative abundance (%)") + #to change label y axis
  scale_y_continuous(labels=percent) #to change y axis tick marks to percentage (need library scales)
                                
                                            
#### BARPLOT at the family level, feed and seawater (ambient, heat) samples only                                           
NOTissue <- subset_samples(Porites_Fam, sampleType == "food" | sampleType == "seawater") #only feed and seawater samples
porites_fam_all<- prune_taxa(taxa_sums(NOTissue) > 0, NOTissue) #to delete taxa which abundance sum is 0 
porites_Gfam_sample_all <- merge_samples(porites_fam_all, "TreatTime", fun=mean)  #I create a new group named Sample 
porites_Gfam_sample_all_abund = transform_sample_counts(porites_Gfam_sample_all, function(x){x / sum(x)})
porites_Gfam_sample_all_abund_melt <- psmelt(porites_Gfam_sample_all_abund)
porites_Gfam_sample_all_abund_melt$Family <- as.character(porites_Gfam_sample_all_abund_melt$Family) # convert Family to a character vector from a factor
porites_Gfam_sample_all_abund_melt$Family[porites_Gfam_sample_all_abund_melt$Abundance < 0.05] <- "Other" #rename classes with < 5% abundance
porites_Gfam_sample_all_abund_melt2 <- porites_Gfam_sample_all_abund_melt %>% 
  separate(Sample, c("treat", "timepoint"), remove=FALSE) #re-creating treatment variable to use in the bar plot
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(17, "Set3"))(nb.cols)
p3 <- ggplot(data=porites_Gfam_sample_all_abund_melt2, aes(x=Sample, y=Abundance, fill=Family))
p3 + geom_bar(aes(), stat="identity", position="stack") +  
  guides(fill=guide_legend(nrow=5)) + labs(title = "porites_Gfam_abund ") +
  scale_fill_manual(values = mycolors) +
  theme_bw() +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )  +
  theme(legend.position="bottom", axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  labs (y = "Mean relative abundance (%)") + #to change label y axis
  scale_y_continuous(labels=percent) #to change y axis tick marks to percentage (need library scales)
         
                                            
#### Boxplot of dominant families
phyloseq::psmelt(OnlyTissue) %>%
  subset (TreatTime != "no.algae" & TreatTime != "no.rotifers") %>%
  subset(Family == "D_4__Alteromonadaceae" | Family == "D_4__Arenicellaceae" | Family == "D_4__Cellvibrionaceae" | Family == "D_4__Endozoicomonadaceae" | Family == "D_4__Flavobacteriaceae" | Family == "D_4__Kiloniellaceae" | Family == "D_4__Nitrincolaceae" | Family == "D_4__Nitrosopumilaceae" | Family == "D_4__Nodosilineaceae" | Family == "D_4__Rhodobacteraceae") %>%
  mutate(name2 = fct_relevel(TreatTime, "baseline.baseline", "ambient.T0", "ambient.T1", "ambient.T2","ambient.T3", "ambient.T4", "ambient.T5", "heat.T0", "heat.T1", "heat.T2", "heat.T3", "heat.T4", "heat.T5")) %>% #reorder x axis
  ggplot(data = ., aes(x = name2, y = Abundance,  fill=treatment)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(height = 0, width = .2, size=0.8) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ Family, scales = "free") +
  theme_bw() +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 25),
    axis.text = element_text(size = 25),
    strip.text.x = element_text(size = 30) #title facets
  )  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(
    values=c("#B3B4B4","#00AFBF", "#B13636"), #colour for treatment
                    guide=guide_legend(reverse=TRUE)) +
  labs (y = "Mean relative abundance (%)") +
  scale_y_continuous(labels=percent) #to change y axis tick marks to percentage (need library scales)
                                            
                                            
#### STATS - adonis - coral vs feed vs seawater microbiota
porites_adonis = as (sample_data(All_noP2_sqrt), "data.frame")
porites_adonis$sampleType = as.factor(porites_adonis$sampleType)
porites_d = phyloseq::distance(All_noP2_sqrt,'bray') 
Adonis_porites <-adonis2(porites_d ~ sampleType, data=porites_adonis,  permutations = 10000, method = "bray")
Adonis_porites 

#dispersion
beta <- betadisper(porites_d, sample_data(All_noP2_sqrt)$sampleType)
boxplot (beta)
disper.test = permutest(beta, permutations =10000)
disper.test # p<0.05 - not ideal
TukeyHSD (beta)
                                            
# post hoc tests
testing_sampleType = pairwise.perm.manova(porites_d, sample_data(All_noP2_sqrt)$sampleType,
                                          nperm=10000, p.method = "BH")
testing_sampleType$p.value
                  
                                            
#### STATS - adonis - coral microbiota: effects of time and treatment                                      
Tissue_noP2_sqrt_noBaseline <- subset_samples (Tissue_noP2_sqrt, treatment != "baseline")
porites_adonis = as (sample_data(Tissue_noP2_sqrt_noBaseline), "data.frame")
porites_adonis$treatment = as.factor(porites_adonis$treatment)
porites_adonis$timePoint = as.factor(porites_adonis$timePoint)
porites_d = phyloseq::distance(Tissue_noP2_sqrt_noBaseline,'bray') 
Adonis_porites <-adonis2(porites_d ~ timePoint*treatment + genotype + tank, data=porites_adonis,  permutations = 10000, method = "bray")
Adonis_porites 
                                            
#dispersion time
beta <- betadisper(porites_d, sample_data(Tissue_noP2_sqrt_noBaseline)$timePoint)
plot(beta)
boxplot (beta)
disper.test = permutest(beta, permutations =10000)
disper.test # OK  
#dispersion treatment                                            
beta <- betadisper(porites_d, sample_data(Tissue_noP2_sqrt_noBaseline)$treatment)
plot(beta)
boxplot (beta)
disper.test = permutest(beta, permutations =10000)
disper.test # OK   
#dispersion parental colony                                               
beta <- betadisper(porites_d, sample_data(Tissue_noP2_sqrt_noBaseline)$genotype)
plot(beta)
boxplot (beta)
disper.test = permutest(beta, permutations =10000)
disper.test # OK    
                                            
# post hoc tests
testing_genotype = pairwise.perm.manova(porites_d, sample_data(Tissue_noP2_sqrt_noBaseline)$genotype,
                                          nperm=10000, p.method = "BH")
testing_genotype$p.value #each parental colony distinct microbiome                                            
sample_data(Tissue_noP2_sqrt_noBaseline)$TimeTemp <-interaction(sample_data(Tissue_noP2_sqrt_noBaseline)$timePoint, sample_data(Tissue_noP2_sqrt_noBaseline)$treatment)
testing_treattime = pairwise.perm.manova(porites_d, sample_data(Tissue_noP2_sqrt_noBaseline)$TimeTemp,
                                          nperm=10000, p.method = "BH")
testing_treattime$p.value                                            
                                            
                                            
####################################################
################## DESEQ ANALYSES ##################
####################################################
phyloseq_merged_porites_noP2_FINAL # I need count data for DESeq (no proportions)
phyloseq_merged_porites_noP2_FINAL_porites <- subset_samples(phyloseq_merged_porites_noP2_FINAL, sampleType == "porites" & treatment != "baseline") #only coral and no baseline samples
phyloseq_merged_porites_noP2_FINAL_porites_noZero <- prune_taxa(taxa_sums(phyloseq_merged_porites_noP2_FINAL_porites) > 0, phyloseq_merged_porites_noP2_FINAL_porites)                                            
sample_data(phyloseq_merged_porites_noP2_FINAL_porites_noZero)$treatment <- as.factor(sample_data(phyloseq_merged_porites_noP2_FINAL_porites_noZero)$treatment) 
sample_data(phyloseq_merged_porites_noP2_FINAL_porites_noZero)$genotype <- as.factor(sample_data(phyloseq_merged_porites_noP2_FINAL_porites_noZero)$genotype)                                            

############### ambient T0 vs heat T0
Filtered <- subset_samples(phyloseq_merged_porites_noP2_FINAL_porites_noZero, TreatTime == "ambient.T0" | TreatTime == "heat.T0")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q1_dds = phyloseq_to_deseq2(Filtered, ~treatment + genotype) #The function phyloseq_to_deseq2 converts your phyloseq-format microbiome data into a DESeqDataSet with dispersions estimated, using the experimental design formula
levels(Q1_dds$treatment) #double check levels as expected
Q1_dds$treatment<-relevel(Q1_dds$treatment, ref="ambient") #I want this as reference
levels(Q1_dds$treatment) #double check it worked (the first listed is the reference)
Q1_dds = DESeq(Q1_dds, test="Wald", fitType="parametric", sfType="poscounts") #run DESEQ function; the default multiple-inference correction is Benjamini-Hochberg, and occurs within the DESeq function. # sfType="poscounts" deals with 0 
head(Q1_dds)
colData(Q1_dds)
resultsNames (Q1_dds) # lists variables of the final model
alpha = 0.01
res_Q1 <- results(Q1_dds, contrast=c("treatment","heat", "ambient"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) #the last one is my reference
head(results(Q1_dds, tidy=TRUE)) #let's look at the results table
summary (res_Q1) #summary of results
res_Q1_order <- res_Q1[order(res_Q1$padj),] #order according to padj 
res_Q1_order = res_Q1_order[order(res_Q1_order$padj, na.last=NA), ]
res_Q1_order_sig = res_Q1_order[(res_Q1_order$padj < alpha), ] #only significant
res_Q1_order_sig #1 ASVs; timePoint heat vs ambient because ambient is the reference
res_Q1_table_sign = cbind(as(res_Q1_order_sig, "data.frame"), as(tax_table(Filtered)[rownames(res_Q1_order_sig), ], "matrix")) #create dataframe with taxonomy

#bubbleplot
filtering1 <- subset_samples(phyloseq_merged_porites_noP2_FINALabundances_cutoff, sampleType == "porites")
table(sample_data(filtering1)$sampleType)
filtering <- subset_samples(filtering1, TreatTime == "ambient.T0" | TreatTime == "heat.T0")
table(sample_data(filtering)$TreatTime)
table(sample_data(filtering)$sampleType)
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("e597cf7906dfdbf035303bf06825ae3c"))
physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Class, Family,  Genus, OTU, sep=",", remove=FALSE, na.rm = TRUE) #I create OTU_fam column
bubbleplot <- ggplot(physeq.subset_m, aes(x = sample, y = OTU_fam, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,20)) +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  facet_grid(~ treatment, scales = "free_x", space = "free_x")
bubbleplot
###no ASVs present in at least 50% of samples within each category


############### ambient T1 vs heat T1
Filtered <- subset_samples(phyloseq_merged_porites_noP2_FINAL_porites_noZero, TreatTime == "ambient.T1" | TreatTime == "heat.T1")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q1_dds = phyloseq_to_deseq2(Filtered, ~treatment + genotype) 
levels(Q1_dds$treatment) 
Q1_dds$treatment<-relevel(Q1_dds$treatment, ref="ambient")
levels(Q1_dds$treatment) 
Q1_dds = DESeq(Q1_dds, test="Wald", fitType="parametric", sfType="poscounts")
head(Q1_dds)
colData(Q1_dds)
resultsNames (Q1_dds)
alpha = 0.01
res_Q1 <- results(Q1_dds, contrast=c("treatment","heat", "ambient"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE)
head(results(Q1_dds, tidy=TRUE))
summary (res_Q1)
res_Q1_order <- res_Q1[order(res_Q1$padj),] 
res_Q1_order = res_Q1_order[order(res_Q1_order$padj, na.last=NA), ]
res_Q1_order_sig = res_Q1_order[(res_Q1_order$padj < alpha), ]
res_Q1_order_sig #5 ASVs
res_Q1_table_sign = cbind(as(res_Q1_order_sig, "data.frame"), as(tax_table(Filtered)[rownames(res_Q1_order_sig), ], "matrix"))

#bubbleplot
filtering1 <- subset_samples(phyloseq_merged_porites_noP2_FINALabundances_cutoff, sampleType == "porites")
table(sample_data(filtering1)$sampleType)
filtering <- subset_samples(filtering1, TreatTime == "ambient.T1" | TreatTime == "heat.T1")
table(sample_data(filtering)$TreatTime)
table(sample_data(filtering)$sampleType)
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("ec706a9c4b9122236e52c792866606cd", "62a42b9060da2e9b281a191e7db519a1", "32bc52e08cc32a1708110af78c22e199",
                                                                              "069f8e6a2a401cce7b6a1292b832c3cc", "e445cce7b3065ed923ad91b18656b594"
))
physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Class, Family,  Genus, OTU, sep=",", remove=FALSE, na.rm = TRUE) 
bubbleplot <- ggplot(physeq.subset_m, aes(x = sample, y = OTU_fam, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,20)) +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  facet_grid(~ treatment, scales = "free_x", space = "free_x")
bubbleplot
###no ASVs present in at least 50% of samples within each category


############### ambient T2 vs heat T2
Filtered <- subset_samples(phyloseq_merged_porites_noP2_FINAL_porites_noZero, TreatTime == "ambient.T2" | TreatTime == "heat.T2")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q1_dds = phyloseq_to_deseq2(Filtered, ~treatment + genotype) 
levels(Q1_dds$treatment) 
Q1_dds$treatment<-relevel(Q1_dds$treatment, ref="ambient") 
levels(Q1_dds$treatment)
Q1_dds = DESeq(Q1_dds, test="Wald", fitType="parametric", sfType="poscounts") 
head(Q1_dds)
colData(Q1_dds)
resultsNames (Q1_dds)
alpha = 0.01
res_Q1 <- results(Q1_dds, contrast=c("treatment","heat", "ambient"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
head(results(Q1_dds, tidy=TRUE)) 
summary (res_Q1) 
res_Q1_order <- res_Q1[order(res_Q1$padj),] 
res_Q1_order = res_Q1_order[order(res_Q1_order$padj, na.last=NA), ]
res_Q1_order_sig = res_Q1_order[(res_Q1_order$padj < alpha), ] 
res_Q1_order_sig #5 ASVs
res_Q1_table_sign = cbind(as(res_Q1_order_sig, "data.frame"), as(tax_table(Filtered)[rownames(res_Q1_order_sig), ], "matrix"))

#bubbleplot
filtering1 <- subset_samples(phyloseq_merged_porites_noP2_FINALabundances_cutoff, sampleType == "porites")
table(sample_data(filtering1)$sampleType)
filtering <- subset_samples(filtering1, TreatTime == "ambient.T2" | TreatTime == "heat.T2")
table(sample_data(filtering)$TreatTime)
table(sample_data(filtering)$sampleType)
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("0b50994d139f92bd27db38f91b5023b7", "19ecabf0b95b8ccb6d75b3ae95178b4e", "d8876415575b0e723bf588954e25ef74",
                                                                              "9fad0e97961eaed64ac123e21461a9aa", "6b9aacca12c4be537807db9f2524fabc"
))
physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Class, Family,  Genus, OTU, sep=",", remove=FALSE, na.rm = TRUE) #I create OTU_fam column
bubbleplot <- ggplot(physeq.subset_m, aes(x = sample, y = OTU_fam, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,20)) +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  facet_grid(~ treatment, scales = "free_x", space = "free_x")
bubbleplot
###no ASVs present in at least 50% of samples within each category
           
                                            
############### ambient T3 vs heat T3
Filtered <- subset_samples(phyloseq_merged_porites_noP2_FINAL_porites_noZero, TreatTime == "ambient.T3" | TreatTime == "heat.T3")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q1_dds = phyloseq_to_deseq2(Filtered, ~treatment + genotype)
levels(Q1_dds$treatment) 
Q1_dds$treatment<-relevel(Q1_dds$treatment, ref="ambient") 
levels(Q1_dds$treatment)
Q1_dds = DESeq(Q1_dds, test="Wald", fitType="parametric", sfType="poscounts")
head(Q1_dds)
colData(Q1_dds)
resultsNames (Q1_dds)
alpha = 0.01
res_Q1 <- results(Q1_dds, contrast=c("treatment","heat", "ambient"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE)
head(results(Q1_dds, tidy=TRUE)) 
summary (res_Q1) 
res_Q1_order <- res_Q1[order(res_Q1$padj),] 
res_Q1_order = res_Q1_order[order(res_Q1_order$padj, na.last=NA), ]
res_Q1_order_sig = res_Q1_order[(res_Q1_order$padj < alpha), ] 
res_Q1_order_sig #5 ASVs
res_Q1_table_sign = cbind(as(res_Q1_order_sig, "data.frame"), as(tax_table(Filtered)[rownames(res_Q1_order_sig), ], "matrix")) 

#bubbleplot
filtering1 <- subset_samples(phyloseq_merged_porites_noP2_FINALabundances_cutoff, sampleType == "porites")
table(sample_data(filtering1)$sampleType)
filtering <- subset_samples(filtering1, TreatTime == "ambient.T3" | TreatTime == "heat.T3")
table(sample_data(filtering)$TreatTime)
table(sample_data(filtering)$sampleType)
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("ead3872a010eab5f7baf2848b98dcee5", "0b50994d139f92bd27db38f91b5023b7", "ab42a580e6e224c7c5205d383a3444df",
                                                                              "7479a25bf010ea5650cc0a934ac5cabf", "a5f2fa7eee368d4fb7dd9ef2134edf63"))
physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Class, Family,  Genus, OTU, sep=",", remove=FALSE, na.rm = TRUE) #I create OTU_fam column
bubbleplot <- ggplot(physeq.subset_m, aes(x = sample, y = OTU_fam, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,20)) +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  facet_grid(~ treatment, scales = "free_x", space = "free_x")
bubbleplot
### 1 ASV present in at least 50% of samples within each category: ead3872a010eab5f7baf2848b98dcee5


############### ambient T4 vs heat T4
Filtered <- subset_samples(phyloseq_merged_porites_noP2_FINAL_porites_noZero, TreatTime == "ambient.T4" | TreatTime == "heat.T4")
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q1_dds = phyloseq_to_deseq2(Filtered, ~treatment + genotype) 
levels(Q1_dds$treatment) 
Q1_dds$treatment<-relevel(Q1_dds$treatment, ref="ambient")
levels(Q1_dds$treatment) 
Q1_dds = DESeq(Q1_dds, test="Wald", fitType="parametric", sfType="poscounts")
head(Q1_dds)
colData(Q1_dds)
resultsNames (Q1_dds)
alpha = 0.01
res_Q1 <- results(Q1_dds, contrast=c("treatment","heat", "ambient"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE)
head(results(Q1_dds, tidy=TRUE)) 
summary (res_Q1)
res_Q1_order <- res_Q1[order(res_Q1$padj),]  
res_Q1_order = res_Q1_order[order(res_Q1_order$padj, na.last=NA), ]
res_Q1_order_sig = res_Q1_order[(res_Q1_order$padj < alpha), ]
res_Q1_order_sig #20 ASVs
res_Q1_table_sign = cbind(as(res_Q1_order_sig, "data.frame"), as(tax_table(Filtered)[rownames(res_Q1_order_sig), ], "matrix")) 

#bubbleplot
filtering1 <- subset_samples(phyloseq_merged_porites_noP2_FINALabundances_cutoff, sampleType == "porites")
table(sample_data(filtering1)$sampleType)
filtering <- subset_samples(filtering1, TreatTime == "ambient.T4" | TreatTime == "heat.T4")
table(sample_data(filtering)$TreatTime)
table(sample_data(filtering)$sampleType)
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("55a3a7dff7f10c5f630c326355572bdf", "e7fc92ade9eef09a7181469aa0f6575f", "effa3245acfbf1689906feb50c144f66",
                                                                              "74eaf8d4a7bbc53947a8b4c3b00cc2ec", "bced101120be21e9dbcd6afdd6cadf29", "c19c9bdafbf13a2e2350768e51743378",
                                                                              "1b8de0faefa8f6891626d91028667b2a", "bda24e4d0e58912087ac259cbc334192", "643f37196749c6ec3ec5915e26bce7b8",
                                                                              "40374b248a7e70f57e32b66900f6d546", "db928247870bf0627a0a0c7cd382795c", "b8807fbf53c538f2671645208c798f11",
                                                                              "849c000404a68b97a3b9e0f8fabf2e42", "04f652b10a1f4b4a0326bdc516a883c4", "dffda5880f1a50ba3287c598412b341e",
                                                                              "0d52758540e7b906a8414b7d90b9f727", "1dcbb8169ffbd18969c5dd4fea6bedae", "7f450be96e9d000e93b555fb75e9de6e",
                                                                              "b65598390a16ce109ce4732249d5b59f", "694df3c7f8b6b66c922ed51a965d75d0"
))

physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Class, Family,  Genus, OTU, sep=",", remove=FALSE, na.rm = TRUE) #I create OTU_fam column
bubbleplot <- ggplot(physeq.subset_m, aes(x = sample, y = OTU_fam, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,20)) +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  facet_grid(~ treatment, scales = "free_x", space = "free_x")
bubbleplot
## 3 ASVs differntially abundant in at least 50 % samples: b65598390a16ce109ce4732249d5b59f, 694df3c7f8b6b66c922ed51a965d75d0, 7f450be96e9d000e93b555fb75e9de6e


############### ambient T5 vs heat T5
Filtered <- subset_samples(phyloseq_merged_porites_noP2_FINAL_porites_noZero, TreatTime == "ambient.T5" | TreatTime == "heat.T5")
table(sample_data(Filtered)$sampleType)
Filtered <- prune_taxa(taxa_sums(Filtered) > 0, Filtered)
Q1_dds = phyloseq_to_deseq2(Filtered, ~treatment + genotype)
levels(Q1_dds$treatment) 
Q1_dds$treatment<-relevel(Q1_dds$treatment, ref="ambient")
levels(Q1_dds$treatment) 
Q1_dds = DESeq(Q1_dds, test="Wald", fitType="parametric", sfType="poscounts")
head(Q1_dds)
colData(Q1_dds)
resultsNames (Q1_dds) 
alpha = 0.01
res_Q1 <- results(Q1_dds, contrast=c("treatment","heat", "ambient"), independentFiltering = TRUE, alpha=alpha, pAdjustMethod = "BH", parallel=TRUE) 
head(results(Q1_dds, tidy=TRUE)) 
summary (res_Q1) 
res_Q1_order <- res_Q1[order(res_Q1$padj),] 
res_Q1_order = res_Q1_order[order(res_Q1_order$padj, na.last=NA), ]
res_Q1_order_sig = res_Q1_order[(res_Q1_order$padj < alpha), ]
res_Q1_order_sig #26 ASVs
res_Q1_table_sign = cbind(as(res_Q1_order_sig, "data.frame"), as(tax_table(Filtered)[rownames(res_Q1_order_sig), ], "matrix"))

#bubbleplot
filtering1 <- subset_samples(phyloseq_merged_porites_noP2_FINALabundances_cutoff, sampleType == "porites")
table(sample_data(filtering1)$sampleType)
filtering <- subset_samples(filtering1, TreatTime == "ambient.T5" | TreatTime == "heat.T5")
table(sample_data(filtering)$TreatTime)
table(sample_data(filtering)$sampleType)
physeq.subset <- subset_taxa(filtering, rownames(tax_table(filtering)) %in% c("753bf49fb914820d19755bec2fd3bd73", "b65598390a16ce109ce4732249d5b59f", "b8807fbf53c538f2671645208c798f11",
                                                                              "849c000404a68b97a3b9e0f8fabf2e42", "760e446f1f1a17e507022634b1bd3f4b", "a08e03bc999c94b24f2512cb2c85677c",
                                                                              "fa3586c81b9c2c808496a040fd751ffa", "27c54877d2a5b58d5b569609cecfcd51", "5aaf9af5f2a07af52e4da802fe1c0b89",
                                                                              "f7f6e9c18ba50fcc4d2e245728062e7c", "cd98aa9fdffe73513428a05e68fdf1da", "4865c6bd53b4c415ae9c9e676cde06be",
                                                                              "694df3c7f8b6b66c922ed51a965d75d0", "6255e5c6afad6ba4213e20224eb4aa15", "598ffdfead8438e27976c14c36098607",
                                                                              "74820d2c3ba5a931dc083dd159abd8bb", "ea9ab787d580ae677de7a846c16a0466", "aab1ddf395b392bc5414a73e7f2e85a1", 
                                                                              "5eab1ddfbbe3615ea476c7e94dabf486", "48a6a0e363dfb3823afdb092a0bb834d", "4f3c49786ad830a7b0b107c93704f4ce",
                                                                              "cdd5d101fdefb45735fd15533ef52d98", "ae3dedc407ae084178fce8f6e060ae13", "36db701374d60d600bcfee2828c30b1c",
                                                                              "eb1516031724c4bc8a71287dfda3aa75", "ab5ecef1c061a5a36d43e89a7fe19fdc"
))
physeq.subset <- prune_taxa(taxa_sums(physeq.subset) > 0, physeq.subset)
physeq.subset_m <- psmelt (physeq.subset)
physeq.subset_m <- physeq.subset_m %>% unite("OTU_fam", Class, Family,  Genus, OTU, sep=",", remove=FALSE, na.rm = TRUE)
bubbleplot <- ggplot(physeq.subset_m, aes(x = sample, y = OTU_fam, color = Class)) + 
  geom_point(aes(size = ifelse(Abundance==0, NA, Abundance))) + 
  scale_size(range = c(3,20)) +
  theme_bw() +
  theme(axis.title.y=element_blank(),  axis.ticks.y=element_blank(),
        legend.text=element_text(size=20), legend.title=element_text(size=20), 
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=22), axis.text.y=element_text(size=22, face="bold")) +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  facet_grid(~ treatment, scales = "free_x", space = "free_x")
bubbleplot
# 15 ASVs differentially abudant in at least 50% samples: b8807fbf53c538f2671645208c798f11, ab5ecef1c061a5a36d43e89a7fe19fdc, 694df3c7f8b6b66c922ed51a965d75d0, aab1ddf395b392bc5414a73e7f2e85a1, ea9ab787d580ae677de7a846c16a0466
#ae3dedc407ae084178fce8f6e060ae13, eb1516031724c4bc8a71287dfda3aa75, 6255e5c6afad6ba4213e20224eb4aa15, 598ffdfead8438e27976c14c36098607, 48a6a0e363dfb3823afdb092a0bb834d, cdd5d101fdefb45735fd15533ef52d98,
#5eab1ddfbbe3615ea476c7e94dabf486, 4f3c49786ad830a7b0b107c93704f4ce, 74820d2c3ba5a931dc083dd159abd8bb, 36db701374d60d600bcfee2828c30b1c

                                            
############### summarizing differentially abundant ASVs in plots

#main
phyloseq::psmelt(phyloseq_merged_porites_noP2_FINALabundances_cutoff) %>%
  subset (TreatTime != "no.algae" & TreatTime != "no.rotifers" & TreatTime != "baseline.baseline") %>% #I have removed baseline because not in deseq analyses
  subset (sampleType == "porites") %>%
  subset (OTU == "694df3c7f8b6b66c922ed51a965d75d0" | OTU == "b65598390a16ce109ce4732249d5b59f"| OTU == "48a6a0e363dfb3823afdb092a0bb834d" |
            
            OTU == "ea9ab787d580ae677de7a846c16a0466" | OTU == "6255e5c6afad6ba4213e20224eb4aa15"| OTU == "7f450be96e9d000e93b555fb75e9de6e"
          ) %>%
  ggplot(data = ., aes(x = TreatTime, y = Abundance,  fill=treatment)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(height = 0, width = .2, size=0.8) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +
  theme_bw() +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme (strip.text = element_text(size = 6, margin = margin())) + 
  scale_fill_manual(values=c("#B4B4B5", "#B13636"), 
                    guide=guide_legend(reverse=TRUE)) +
  labs (y = "Mean relative abundance (%)") +
  scale_y_continuous(labels=percent) 
                                            
#supplementary figure
phyloseq::psmelt(phyloseq_merged_porites_noP2_FINALabundances_cutoff) %>%
  subset (TreatTime != "no.algae" & TreatTime != "no.rotifers" & TreatTime != "baseline.baseline") %>% #I have removed baseline because not in deseq analyses
  subset (sampleType == "porites") %>%
  subset (OTU == "ead3872a010eab5f7baf2848b98dcee5" | OTU == "b8807fbf53c538f2671645208c798f11"| OTU == "ab5ecef1c061a5a36d43e89a7fe19fdc"| 
            OTU ==  "aab1ddf395b392bc5414a73e7f2e85a1"| OTU == "ae3dedc407ae084178fce8f6e060ae13"| OTU == "eb1516031724c4bc8a71287dfda3aa75"| 
            OTU == "598ffdfead8438e27976c14c36098607"| OTU == "cdd5d101fdefb45735fd15533ef52d98"| OTU == "5eab1ddfbbe3615ea476c7e94dabf486"| 
            OTU =="4f3c49786ad830a7b0b107c93704f4ce"| OTU == "74820d2c3ba5a931dc083dd159abd8bb"| OTU == "36db701374d60d600bcfee2828c30b1c") %>%
  ggplot(data = ., aes(x = TreatTime, y = Abundance,  fill=treatment)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(height = 0, width = .2, size=0.8) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free") +
  theme_bw() +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme (strip.text = element_text(size = 6, margin = margin())) + 
  scale_fill_manual(values=c("#B4B4B5", "#B13636"), 
                    guide=guide_legend(reverse=TRUE)) +
  labs (y = "Mean relative abundance (%)") +
  scale_y_continuous(labels=percent) 
                                           
                                            
                                            
################################################################
#################### ALPHA DIVERSITY ANALYSES ##################
################################################################   
                                            
#I use rarefied data for Shannon diversity 
rarefied_Porites_noP2_OnlyTissue <- subset_samples(phyloseq_merged_porites_noP2_FINALrarefied, sampleType == "porites")

#alpha diversity plot                                            
plot_richness(rarefied_Porites_noP2_OnlyTissue, x="timePoint", color = "treatment", measures=c("Shannon")) + 
  geom_boxplot(alpha=0.3, colour = "black", lwd=0.2) + labs(title = "Shannon") +
  theme( axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + aes (fill=treatment) +
  theme_bw() +
  theme( 
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  scale_color_manual(values = Colors) +
  scale_fill_manual (values = Colors)                                            
                                            
#alpha diversity stats                                            
rich = estimate_richness(rarefied_Porites_noP2_OnlyTissue, split=TRUE)
rich<-as.data.frame(rich)
RICHNESSVER<-tibble::rownames_to_column(as.data.frame(rich), var="sampleID")
RICHNESSVER$sampleID = gsub("X", "", RICHNESSVER$sampleID)
metadata3 <- tibble::rownames_to_column(as.data.frame(metadata2), var="sampleID")
richness_table_ver <- inner_join(RICHNESSVER, metadata3, by="sampleID")
porites_richness <- as.data.frame(richness_table_ver)                                            
richness_table_porites_NObaseline <- subset (porites_richness, DHW != "no") #remove baseline samples
richness_table_porites_NObaseline$timePoint <- as.factor(richness_table_porites_NObaseline$timePoint)
richness_table_porites_NObaseline$treatment <- as.factor(richness_table_porites_NObaseline$treatment)
richness_table_porites_NObaseline$genotype <- as.factor(richness_table_porites_NObaseline$genotype)
model_Shannon_porites <- glmmTMB(Shannon ~ timePoint*treatment+genotype + (1|tank), data = richness_table_porites_NObaseline)
S.resid <- model_Shannon_porites %>% simulateResiduals(plot=TRUE, integerResponse = TRUE) #good!
model_Shannon_porites %>% summary ()
model_Shannon_porites %>% emmeans (~treatment|timePoint) %>% pairs %>%rbind(adjust='sidak')                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
