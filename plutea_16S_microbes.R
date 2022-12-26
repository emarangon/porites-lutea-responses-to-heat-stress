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
         
                                            
######### boxplot of dominant families
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
                                            
                                            
#### STATS - adonis
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
                                            
