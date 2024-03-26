
######################################################################################################################
################### plutea - heat stress experiment - physiology #################
######################################################################################################################

##author EMMA MARANGON


###loading libraries

library(dplyr) #data wrangling
library(ggplot2) #graphs
library(Rmisc) #calculation mean SE etc
library (glmmTMB)  #models
library (DHARMa) #model assumption
library(MuMIn)  #r.squaredGLMM
library(emmeans) #post hocs
library (tidyr) 



##################################################################################
#################### Photosynthesis, Respiration and PR ratio ####################
##################################################################################

###preprocessig
All <- read.csv("./PR_Porites.csv", header = TRUE)
All <- All %>% 
  mutate(Parent = factor(Parent),
         Tank = factor(Tank),
         Treatment = factor(Treatment),
         TimePoint = factor(TimePoint),
         Sample.ID=factor(Sample.ID)) 
All$Time <- as.numeric(All$Time) #time is continuous
All$temp_celsius <- as.numeric(All$temp_celsius)
class(All$Time)
class(All$temp_celsius)
class(All$Tank)
class (All$TimePoint)
All <- subset(All, Respiration != "NA" ) %>% droplevels()
All <- subset(All, Photosynthesis != "NA" )  %>% droplevels()
All <- subset(All, RatioPR != "NA" )  %>% droplevels()
summary(All)


###data exploration
ggplot(All, aes(y=Photosynthesis, x=TimePoint, color=Treatment))+
  geom_boxplot()+ 
  theme_classic() +
  geom_point(position=position_dodge(width=0.75),aes(group=Treatment))

ggplot(All, aes(y=Respiration, x=TimePoint, color=Treatment))+
  geom_boxplot()+ 
  theme_classic() +
  geom_point(position=position_dodge(width=0.75),aes(group=Treatment))

ggplot(All, aes(y=RatioPR, x=TimePoint, color=Treatment))+
  geom_boxplot()+ #points
  theme_classic() +
  geom_point(position=position_dodge(width=0.75),aes(group=Treatment))


###final plots

#Photosynthesis
tgc <- summarySE(All, measurevar="Photosynthesis", groupvars=c("Time", "Treatment"))
tgc

ggplot(tgc, aes(Time, Photosynthesis, colour = Treatment)) +
  geom_line() +
  scale_x_continuous(breaks=c(0, 17, 25, 32, 39, 46, 53)) + 
    theme_classic() + geom_line(size=0.7)+ #line width
scale_color_manual(values=c("#B4B4B4", "#AC664B")) + 
geom_point (aes(x = Time, y = Photosynthesis, fill=Treatment),
                        shape=21) +
  scale_fill_manual(values=c("#B4B4B4", "#AC664B")) +
  geom_errorbar(aes(ymin=Photosynthesis-se, ymax=Photosynthesis+se), 
                width=.4, lwd=.7,  
                position=position_dodge(.05)) +
  ggtitle("Porites Net Photosynthesis") 

#Respiration
tgc <- summarySE(All, measurevar="Respiration", groupvars=c("Time", "Treatment"))
tgc

ggplot(tgc, aes(Time, Respiration, colour = Treatment)) +
  geom_line() +
  scale_x_continuous(breaks=c(0, 17, 25, 32, 39, 46, 53)) +
    theme_classic() + geom_line(size=0.7)+ 
scale_color_manual(values=c("#B4B4B4", "#AC664B")) +
geom_point (aes(x = Time, y = Respiration, fill=Treatment),
                        shape=21) +
  scale_fill_manual(values=c("#B4B4B4", "#AC664B")) +
  geom_errorbar(aes(ymin=Respiration-se, ymax=Respiration+se), 
                width=.4, lwd=.7,  
                position=position_dodge(.05)) +
  ggtitle("Porites Respiration") 

#PRratio
tgc <- summarySE(All, measurevar="RatioPR", groupvars=c("Time", "Treatment"))
tgc

ggplot(tgc, aes(Time, RatioPR, colour = Treatment)) +
  geom_line() +
  scale_x_continuous(breaks=c(0, 17, 25, 32, 39, 46, 53)) +
    theme_classic() + geom_line(size=0.7)+ #line width
scale_color_manual(values=c("#B4B4B4", "#AC664B")) +
geom_point (aes(x = Time, y = RatioPR, fill=Treatment),
                        shape=21) +
  scale_fill_manual(values=c("#B4B4B4", "#AC664B")) +
  geom_errorbar(aes(ymin=RatioPR-se, ymax=RatioPR+se), 
                width=.4, lwd=.7,  
                position=position_dodge(.05)) +
  ggtitle("Porites P/R ratio")



###statistics (final models)

#Photosynthesis
All_posP <- All %>% dplyr::mutate(pos_photosynthesis = (Photosynthesis+0.009)) #to have all positive values
All_posP_sqrt <- All_posP %>% dplyr::mutate(sqrt_pos_photosynthesis = sqrt(pos_photosynthesis)) #sqrt transformation
Pmod11 = glmmTMB(sqrt_pos_photosynthesis ~ Parent+TimePoint*Treatment + 
                  (1|Tank/Sample.ID), data = All_posP_sqrt, family='gaussian', REML = TRUE)
P.resid <- Pmod11 %>% simulateResiduals(plot =TRUE) #assumptions met!
Pmod11 %>% summary () #interaciton sign
Pmod11 %>% r.squaredGLMM() #our model explains 78% of variability.
Pmod11 %>% emmeans (~Treatment|TimePoint) %>% pairs %>%rbind(adjust='sidak') 

#Respiration
All_abs_respiration <- All %>% dplyr::mutate(abs_respiration = abs(Respiration))
Rmod7 = glmmTMB(abs_respiration ~ Parent+TimePoint*Treatment + 
                  (1|Tank/Sample.ID), data = All_abs_respiration, family='Gamma', REML = TRUE)
R.resid <- Rmod7 %>% simulateResiduals(plot =TRUE)  #assumptions met!
Rmod7 %>% summary () #interaction time point treatment sign
Rmod7 %>% emmeans (~Treatment|TimePoint) %>% pairs %>%rbind(adjust='sidak')

#PR ratio
PRmod5 = glmmTMB(RatioPR ~ Parent+Treatment*TimePoint + 
                  (1|Tank/Sample.ID), data = All, family='Gamma', REML = TRUE)
PR.resid <- PRmod5 %>% simulateResiduals(plot =TRUE) #assumptions met!
PRmod5 %>% summary () #Treatment*TimePoint sign
PRmod5 %>% emmeans (~Treatment|TimePoint) %>% pairs %>%rbind(adjust='sidak')




#######################################################
#################### Health scores ####################
#######################################################

###preprocessing
All <- read.csv("./BleachingScores_porites.csv", header = TRUE)
All <- All %>% 
  mutate(Parent = factor(Parent),
         Tank = factor(Tank),
         Treatment = factor(Treatment),
         TimePoint = factor(TimePoint),
         dhw = factor(dhw),
         Sample.ID = factor (Sample.ID))
All$Time <- as.numeric(All$Time) #time is continuous 
All$temp_celsius <- as.numeric(All$temp_celsius)
summary(All)
class(All$Time)
class(All$temp_celsius)
class(All$Tank)


###data exploration
ggplot(All, aes(y=HealthScore, x=TimePoint, color=Treatment))+
  geom_boxplot()+ 
  theme_classic() +
  geom_point(position=position_dodge(width=0.75),aes(group=Treatment))


###final plot
tgc <- summarySE(All, measurevar="HealthScore", groupvars=c("TimePoint", "Treatment"))
tgc

tgc2 <- tgc %>% #add a Time column
  mutate(
    Time = case_when(
      TimePoint == "T0" ~ "0",
      TimePoint == "T1" ~ "19",
      TimePoint == "T2" ~ "27",
      TimePoint == "T3" ~ "34",
      TimePoint == "T4" ~ "41",
      TimePoint == "T5" ~ "56"
    )
  )
tgc2$Time <- as.numeric(tgc2$Time)

ggplot(tgc2, aes(Time, HealthScore, colour = Treatment)) +
  geom_line() +
  scale_x_continuous(breaks=c(0, 19, 27, 34, 41, 56)) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
    theme_classic() + geom_line(size=0.7)+ 
scale_color_manual(values=c("#B4B4B4", "#AC664B")) + 
geom_point (aes(x = Time, y = HealthScore, fill=Treatment),
                        shape=21) + #to include points on lines
  scale_fill_manual(values=c("#B4B4B4", "#AC664B")) + 
  geom_errorbar(aes(ymin=HealthScore-se, ymax=HealthScore+se), 
                width=.4, lwd=.7,  
                position=position_dodge(.05)) +
  ggtitle("Porites Health score")


###statistics (final model)
Bmod5 = glmmTMB(HealthScore ~ Parent+TimePoint*Treatment + 
                  (1|Tank/Sample.ID), data = All, family='gaussian', REML = TRUE)
B.resid <- Bmod5 %>% simulateResiduals(plot =TRUE) #assumptions met!
Bmod5 %>% summary ()
Bmod5 %>% r.squaredGLMM() #our model explains 55% of variability
Bmod5 %>% emmeans (~Treatment|TimePoint) %>% pairs %>%rbind(adjust='sidak')



############################################################################
#################### Photochemical effective efficiency ####################
############################################################################

###preprocessing
Light_1_filtered_sd <- read.csv("./PhotochemicalEfficiency_porites.csv", header = TRUE, sep=",")
Light_2_filtered_sd <- unite(Light_1_filtered_sd, ID2, TimePoint:Light_Dark, remove=FALSE)
Light_2_filtered_sd$ID2 <- as.factor(Light_2_filtered_sd$ID2) 

Light_3_filtered_sd <- unite(Light_2_filtered_sd, ID3, ID2:Sample_ID, remove=FALSE)
Light_3_filtered_sd$ID3 <- as.factor(Light_3_filtered_sd$ID3)

Light_4_filtered_sd <- unite(Light_3_filtered_sd, ID4, ID3:Tank, remove=FALSE)
Light_4_filtered_sd$ID4 <- as.factor(Light_4_filtered_sd$ID4)

Light_5_filtered_sd <- unite(Light_4_filtered_sd, ID5, ID4:Treatment, remove=FALSE)
Light_5_filtered_sd$ID5 <- as.factor(Light_5_filtered_sd$ID5)

Light_6_filtered_sd <- ddply(Light_5_filtered_sd,~ID5,summarise,meanY=mean(Y))

Light_7_filtered_sd <-Light_6_filtered_sd %>%
  separate (ID5, c("TimePoint", "Light_Dark", "Sample_ID", "Tank", "Treatment"), "_")
Light_7_filtered_sd$Tank <- as.factor(Light_7_filtered_sd$Tank)
Light_7_filtered_sd$Light_Dark <- factor(Light_7_filtered_sd$Light_Dark)
Light_7_filtered_sd$TimePoint <- factor(Light_7_filtered_sd$TimePoint)
Light_7_filtered_sd$Treatment <- factor(Light_7_filtered_sd$Treatment)
Light_7_filtered_sd$Sample_ID <- factor(Light_7_filtered_sd$Sample_ID)
summary (Light_7_filtered_sd)

Light_7_filtered_sd$Parent = Light_7_filtered_sd$Sample_ID
Light_7_filtered_sd %>% mutate(Parent = substr(Parent, 1, 2)) -> Light_8_filtered_sd
names(Light_8_filtered_sd)[names(Light_8_filtered_sd)=="meanY"] <- "PhotochemicalEfficiency"

Light_filtered_sd <- Light_8_filtered_sd %>% 
  mutate(Parent = factor(Parent),
         Tank = factor(Tank),
         Treatment = factor(Treatment),
         TimePoint = factor (TimePoint),
         Sample_ID = factor(Sample_ID))
summary(Light_filtered_sd)



###final plot
Light_filtered_sd_graph <- Light_filtered_sd %>% 
  mutate(
    Time = case_when(
      TimePoint == "T0" ~ "0",
      TimePoint == "T1" ~ "19",
      TimePoint == "T2" ~ "27",
      TimePoint == "T3" ~ "34",
      TimePoint == "T4" ~ "41",
      TimePoint == "T5" ~ "56"
    )
  )

tgc <- summarySE(Light_filtered_sd_graph, measurevar="PhotochemicalEfficiency", groupvars=c("Time", "Treatment"))
tgc

tgc$Time = as.numeric(tgc$Time) #need to use this variabe as numeric if I want to add geom_line (lines) to graph

ggplot(tgc, aes(Time, PhotochemicalEfficiency, colour = Treatment)) +
  geom_line() +
  scale_x_continuous(breaks=c(0, 19, 27, 34, 41, 56)) + 
    theme_classic() + geom_line(size=0.7)+ 
  scale_color_manual(values=c("#B4B4B4", "#AC664B")) + 
  geom_point (aes(x = Time, y = PhotochemicalEfficiency, fill=Treatment),
                        shape=21) + 
  scale_fill_manual(values=c("#B4B4B4", "#AC664B")) + 
  geom_errorbar(aes(ymin=PhotochemicalEfficiency-se, ymax=PhotochemicalEfficiency+se), 
                width=.4, lwd=.7,  
                position=position_dodge(.05)) +
  ggtitle("Porites Photochemical Effective Efficiency") 


###statistics (final model)
Light_filtered_sd_sqrt <- Light_filtered_sd %>% dplyr::mutate(sqrt_efficiency = sqrt(PhotochemicalEfficiency))
Emod6 = glmmTMB(sqrt_efficiency ~ Parent+Treatment*TimePoint + 
                  (1|Tank/Sample_ID), data = Light_filtered_sd_sqrt, beta_family (link='logit'), REML = TRUE)

B.resid <- Emod6 %>% simulateResiduals(plot=TRUE, integerResponse = TRUE) #acceptable
Emod6 %>% summary ()
Emod6 %>% emmeans (~Treatment|TimePoint) %>% pairs %>%rbind(adjust='sidak')



