
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
All_noP2 <- subset (All, Parent !="2") %>% droplevels # I remove P2 because it was identified as Porites lobata
summary(All_noP2)

###data exploration
ggplot(All_noP2, aes(y=Photosynthesis, x=TimePoint, color=Treatment))+
  geom_boxplot()+ 
  theme_classic() +
  geom_point(position=position_dodge(width=0.75),aes(group=Treatment))

ggplot(All_noP2, aes(y=Respiration, x=TimePoint, color=Treatment))+
  geom_boxplot()+ 
  theme_classic() +
  geom_point(position=position_dodge(width=0.75),aes(group=Treatment))

ggplot(All_noP2, aes(y=RatioPR, x=TimePoint, color=Treatment))+
  geom_boxplot()+ #points
  theme_classic() +
  geom_point(position=position_dodge(width=0.75),aes(group=Treatment))


###final plots

#Photosynthesis
tgc <- summarySE(All_noP2, measurevar="Photosynthesis", groupvars=c("Time", "Treatment"))
tgc

ggplot(tgc, aes(Time, Photosynthesis, colour = Treatment)) +
  geom_line() +
  scale_x_continuous(breaks=c(0, 17, 25, 32, 39, 46, 53)) + #to set tick marks on x axis
    theme_classic() + geom_line(size=0.7)+ #line width
scale_color_manual(values=c("#B4B4B4", "#AC664B")) + #colour lines
geom_point (aes(x = Time, y = Photosynthesis, fill=Treatment),
                        shape=21) + #to include points on lines
  scale_fill_manual(values=c("#B4B4B4", "#AC664B")) + #colour points
  geom_errorbar(aes(ymin=Photosynthesis-se, ymax=Photosynthesis+se), 
                width=.4, lwd=.7,  
                position=position_dodge(.05)) +
  ggtitle("Porites Photosynthesis") 

#Respiration
tgc <- summarySE(All_noP2, measurevar="Respiration", groupvars=c("Time", "Treatment"))
tgc

ggplot(tgc, aes(Time, Respiration, colour = Treatment)) +
  geom_line() +
  scale_x_continuous(breaks=c(0, 17, 25, 32, 39, 46, 53)) +
    theme_classic() + geom_line(size=0.7)+ #line width
scale_color_manual(values=c("#B4B4B4", "#AC664B")) +
geom_point (aes(x = Time, y = Respiration, fill=Treatment),
                        shape=21) +
  scale_fill_manual(values=c("#B4B4B4", "#AC664B")) +
  geom_errorbar(aes(ymin=Respiration-se, ymax=Respiration+se), 
                width=.4, lwd=.7,  
                position=position_dodge(.05)) +
  ggtitle("Porites Respiration") 

#PRratio
tgc <- summarySE(All_noP2, measurevar="RatioPR", groupvars=c("Time", "Treatment"))
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
  ggtitle("Porites RatioPR")



###statistics (final models)

#Photosynthesis
All_noP2_posP <- All_noP2 %>% dplyr::mutate(pos_photosynthesis = (Photosynthesis+0.009)) #to have all positive values
All_noP2_posP_sqrt <- All_noP2_posP %>% dplyr::mutate(sqrt_pos_photosynthesis = sqrt(pos_photosynthesis)) #sqrt transformation
Pmod11 = glmmTMB(sqrt_pos_photosynthesis ~ Parent+TimePoint*Treatment + 
                  (1|Tank/Sample.ID), data = All_noP2_posP_sqrt, family='gaussian', REML = TRUE)
P.resid <- Pmod11 %>% simulateResiduals(plot =TRUE) #assumptions met!
Pmod11 %>% summary () #interaciton sign
Pmod11 %>% r.squaredGLMM() #our model explains 78% of variability.
Pmod11 %>% emmeans (~Treatment|TimePoint) %>% pairs %>%rbind(adjust='sidak') 

#Respiration
All_noP2_abs_respiration <- All_noP2 %>% dplyr::mutate(sqrt_abs_respiration = sqrt(abs_respiration))
Rmod7 = glmmTMB(abs_respiration ~ Parent+TimePoint*Treatment + 
                  (1|Tank/Sample.ID), data = All_noP2_abs_respiration, family='Gamma', REML = TRUE)
R.resid <- Rmod7 %>% simulateResiduals(plot =TRUE)  #assumptions met!
Rmod7 %>% summary () #interaction time point treatment sign
Rmod7 %>% emmeans (~Treatment|TimePoint) %>% pairs %>%rbind(adjust='sidak')

#PR ratio
PRmod5 = glmmTMB(RatioPR ~ Parent+Treatment*TimePoint + 
                  (1|Tank/Sample.ID), data = All_noP2, family='Gamma', REML = TRUE)
PR.resid <- PRmod5 %>% simulateResiduals(plot =TRUE) #assumptions met!
PRmod5 %>% summary () #Treatment*TimePoint sign
PRmod5 %>% emmeans (~Treatment|TimePoint) %>% pairs %>%rbind(adjust='sidak')
