setwd("C:/Users/davin/Documents/PhD/Results/Westerns/")
library(readxl)
library(tidyverse)
library(reshape2)
library(openxlsx)
library(directlabels)
library(gtools)
library(gridExtra)
library(ggpubr)
library(dplyr)
library(tidyr)
library(plotrix)
library(extrafont)
library(rstatix)
library(grafify)
library(viridis)

####plotting###
fold <- read.xlsx("foldonly.xlsx")
fold

fold_mu<- fold %>% rowwise() %>% mutate(avg=mean(c(r1,r2,r3)), se=std.error(c(r1,r2,r3))) #GOOD
fold_mu

#order samples by order of sds page
ord<-c("OST1-MUT1", "MFA-MFA", "HSP150-HSP150", "HSP150-MFA", "HSP150-MUT1", "PIR1-PIR1", "PIR1-MFA","PIR1-MUT1", "BP")
fold_mu$sample <-factor(fold_mu$sample, levels=unique(ord))

long <- melt(fold_mu)
long

#tukey and anova
tuk.test <- function(t){
  tuk <- TukeyHSD(t) #perform TukeyHSD test on anova results
  df <- tuk$sample %>% as.data.frame() #keep sample column, change results to data frame 
  fp <- filter(df, df$`p adj`<0.05) #keep only significant values <0.05
  ##  out <- fp[grep("BP", rownames(fp)), ] #keep only rows with BP in
  return(fp)
}

long #CONTAINS AVG AND SE 
aovout <- aov(value~sample, data=long) #WRONG
summary(aovout) #f=3.637 p=0.0034 ** so continue to tukey

tuk1<-tuk.test(aovout) #TUKEY
tuk1 #                      diff       lwr        upr       p adj
##BP-HSP150-HSP150 -1.0482035 -1.860088 -0.2363188 0.004006943
##BP-PIR1-PIR1     -0.9962936 -1.808178 -0.1844090 0.007204617

#student's t-test version
fonly<- melt(fold)
fonly
test1 <- fonly %>% t_test(value ~ sample) %>% add_significance()
test1 %>% filter(p.adj<0.05)

?t_test
#pairwise t-test? not sure which is better...
#maybe don't do fold change? otherwise can't test against mfa-mfa

#aov and tukey with fonly
aovout3 <- aov(value~sample, data=fonly)
summary(aovout3) #F=70.75 p=5.08e-12 *** so continue to tukey
tuk3<-tuk.test(aovout3) #TUKEY
tuk3 #bare results, tuk1 is definitely wrong as it contains mus already which is messing up the analysis




#####one colour plot#### - not in thesis
ggplot(fold_mu, aes(x = sample, y = avg)) +
  geom_col(position = position_dodge(width = 0.5), color = "black", size = 0.8, fill="#909fc1", alpha=1) +
  geom_point(aes(y = r1), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r2), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r3), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 0.2) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, by = 0.25), expand = expansion(mult = c(0, 0.05))) +
  ##scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0, 1.75)) +
  xlab("Secretion Signal") +
  ylab("Fold change") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        text = element_text(family = "Calibri Light"),
        ##legend.title = element_blank(),
        axis.line = element_line(),
        axis.ticks = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")



#####colours based on pre-signal, and with legend for pro signal#### - not in thesis
ggplot(fold_mu, aes(x = sample, y = avg, fill=sample)) +
  geom_col(position = position_dodge(width = 0.5), color = "black", size = 0.8) +
  #changing x-axis labels from ord to split pre and pro signal schemes
  scale_x_discrete(labels = c("OST1-MUT1" = "MUT1", "MFA-MFA"="MFA", 
                              "HSP150-HSP150" = "HSP150", "HSP150-MFA" = "MFA", 
                              "HSP150-MUT1" = "MUT1", 
                              "PIR1-PIR1"="PIR1", "PIR1-MFA"="MFA", "PIR1-MUT1"="MUT1", "BP"="BP")) +
  scale_fill_manual(values=c("#CCB7AE","#D6CFCB","#706677","#706677","#706677","#a7b0cb","#a7b0cb","#a7b0cb","white")) +
  geom_point(aes(y = r1), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r2), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r3), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 0.2) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, by = 0.25), expand = expansion(mult = c(0, 0.05))) +
  ##scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0, 1.75)) +
  xlab("Secretion Signal") +
  ylab("Fold change") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        text = element_text(family = "Calibri Light"),
        legend.position = "right",
        axis.line = element_line(),
        axis.ticks = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")


###raw values plot###
raw <- read.xlsx("rawonly.xlsx")
raw

raw_mu<- raw %>% rowwise() %>% mutate(avg=mean(c(r1,r2,r3)), se=std.error(c(r1,r2,r3))) #GOOD
raw_mu

#order samples by order of sds page
##ord<-c("OST1-MUT1", "MFA-MFA", "HSP150-HSP150", "HSP150-MFA", "HSP150-MUT1", "PIR1-PIR1", "PIR1-MFA","PIR1-MUT1", "BP")
raw_mu$sample <-factor(raw_mu$sample, levels=unique(ord))

rawl <- melt(raw_mu)
rawl #has avg and se


#plot raw values - not in thesis
rplot <- ggplot(raw_mu, aes(x = sample, y = avg, fill=sample)) +
  geom_col(position = position_dodge(width = 0.5), color = "black", size = 0.8) +
  #changing x-axis labels from ord to split pre and pro signal schemes
  scale_x_discrete(labels = c("OST1-MUT1" = "MUT1", "MFA-MFA"="MFA", 
                              "HSP150-HSP150" = "HSP150", "HSP150-MFA" = "MFA", 
                              "HSP150-MUT1" = "MUT1", 
                              "PIR1-PIR1"="PIR1", "PIR1-MFA"="MFA", "PIR1-MUT1"="MUT1", "BP"="BP")) +
  scale_fill_manual(values=c("#CCB7AE","#D6CFCB","#706677","#706677","#706677","#a7b0cb","#a7b0cb","#a7b0cb","white")) +
  geom_point(aes(y = r1), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r2), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r3), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 0.2) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 5000), expand = expansion(mult = c(0, 0.05))) +
  ##scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0, 1.75)) +
  xlab("Secretion Signal") +
  ylab("Pixel Coverage (Area)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        text = element_text(family = "Calibri Light"),
        legend.position = "right",
        axis.line = element_line(),
        axis.ticks = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1)) 
rplot

#ANOVAs and Tukey tests
ronly
aovout2 <- aov(value~sample, data=ronly)
summary(aovout2) #F=33.21  p=3.08e-09 *** so continue to TUKEY tests

tuk2<-tuk.test(aovout2) #TUKEY
tuk2
write.csv(tuk2,"pixel_aov.csv")



#plot sig for mfa vs everything

#START HERE 10/10/23
#plot not in thesis
rplot + geom_signif(comparisons = list(c("MFA-MFA", "HSP150-HSP150")),  #mfa-mfa vs new signals with increased expr
            annotations="*", y_position = 25000, tip_length = 0.0, vjust = 0.5) +
         geom_signif(comparisons = list(c("OST1-MUT1", "HSP150-HSP150")),  #published ost1 vs new signals 
              annotations="***", y_position = 27500, tip_length = 0.0, vjust = 0.5) +
         geom_signif(comparisons = list(c("OST1-MUT1", "PIR1-PIR1")),  
              annotations="**", y_position = 28500, tip_length = 0.0, vjust = 0.5) +
          geom_signif(comparisons = list(c("MFA-MFA", "OST1-MUT1")),  
              annotations="ns", y_position = 26500, tip_length = 0.0, vjust = 0.1, textsize = 3) 
      
  

tuk2
TukeyHSD(aovout2)

#####Delete - significance testing for increase only####
# Perform an independent two-sample t-test
ronly
raw
MFA_MFA <- raw[raw$sample == "MFA-MFA", c("r1", "r2", "r3")]
MFA_MFA
HSP150_HSP150<- raw[raw$sample == "HSP150-HSP150", c("r1", "r2", "r3")]
HSP150_HSP150
##tr <- t(raw)
t_MH <- t.test(MFA_MFA, HSP150_HSP150)

# Print the result
print(t_MH) #t = -4.1728, df = 3.1541, p-value = 0.02271


#for PIR1
raw
PIR1_PIR1<- raw[raw$sample == "PIR1-PIR1", c("r1", "r2", "r3")]
PIR1_PIR1
##tr <- t(raw)
t_MP <- t.test(MFA_MFA, PIR1_PIR1)
print(t_MP) #t = -2.8924, df = 2.8675, p-value = 0.0664 - same ns result as ANOVA/TUKEY


####Delete - BP NORM####
#using values that are normalised to the negative control (BP)
bpn <- read.xlsx("bpnorm.xlsx")
bpn

bpn_mu<- bpn %>% rowwise() %>% mutate(avg=mean(c(r1,r2,r3)), se=std.error(c(r1,r2,r3))) #GOOD
bpn_mu

#order samples by order of sds page
##ord<-c("OST1-MUT1", "MFA-MFA", "HSP150-HSP150", "HSP150-MFA", "HSP150-MUT1", "PIR1-PIR1", "PIR1-MFA","PIR1-MUT1", "BP")

bpn_mu$sample <-factor(bpn_mu$sample, levels=unique(ord))
bpn_mu

bplot <- ggplot(bpn_mu, aes(x = sample, y = avg, fill=sample)) +
  geom_col(position = position_dodge(width = 0.5), color = "black", size = 0.8) +
  #changing x-axis labels from ord to split pre and pro signal schemes
  scale_x_discrete(labels = c("OST1-MUT1" = "MUT1", "MFA-MFA"="MFA", 
                              "HSP150-HSP150" = "HSP150", "HSP150-MFA" = "MFA", 
                              "HSP150-MUT1" = "MUT1", 
                              "PIR1-PIR1"="PIR1", "PIR1-MFA"="MFA", "PIR1-MUT1"="MUT1")) +
  scale_fill_manual(values=c("#CCB7AE","#D6CFCB","#706677","#706677","#706677","#a7b0cb","#a7b0cb","#a7b0cb")) +
  geom_point(aes(y = r1), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r2), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r3), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 0.2) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 5000), expand = expansion(mult = c(0, 0.05))) +
  ##scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0, 1.75)) +
  xlab("Secretion Signal") +
  ylab("Pixel Coverage (Area)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        text = element_text(family = "Calibri Light"),
        legend.position = "right",
        axis.line = element_line(),
        axis.ticks = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1)) 

bplot


#sig tests
bpnl <- melt(bpn) #long format for anova
bpnl
aovout4 <- aov(value~sample, data=bpnl)
summary(aovout4) #F=12.27  p=2.26e-05 ***, going ahead with tukey

tuk4<-tuk.test(aovout4) #TUKEY
tuk4
write.csv(tuk4,"pixel_aov4.csv")


bplot + geom_signif(comparisons = list(c("MFA-MFA", "HSP150-HSP150")),  
                    annotations="*", y_position = 25000, tip_length = 0, vjust = 0.5)






####splitting plots####
#native signals
nat <- read.xlsx("raw_native.xlsx")
nat

nat_mu<- nat %>% rowwise() %>% mutate(avg=mean(c(r1,r2,r3)), se=std.error(c(r1,r2,r3))) #GOOD
nat_mu

#order samples by order of sds page
nat_mu$sample <-factor(nat_mu$sample, levels=unique(ord))

nplot <- ggplot(nat_mu, aes(x = sample, y = avg, fill=sample)) +
  geom_col(position = position_dodge(width = 0.5), color = "black", size = 0.8) +
  #changing x-axis labels from ord to split pre and pro signal schemes
  scale_x_discrete(labels = c("OST1-MUT1" = "OST1-MUT1", 
                              "MFA-MFA"="MFA", 
                              "HSP150-HSP150" = "HSP150", 
                              "PIR1-PIR1"="PIR1", 
                              "BP"="BP")) +
  scale_fill_manual(values=c("#D6CFCB","#CCB7AE","#9A91A1","#BFC5D9","white")) +
  geom_point(aes(y = r1), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r2), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r3), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 0.2) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 5000), expand = expansion(mult = c(0, 0.05))) +
  ##scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0, 1.75)) +
  xlab("Secretion Signal") +
  ylab("Pixel Coverage (Area)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        text = element_text(family = "Calibri Light"),
        legend.position = "none",
        axis.line = element_line(),
        axis.ticks = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1)) 

A<-nplot + geom_signif(comparisons = list(c("MFA-MFA", "HSP150-HSP150")),  #mfa-mfa vs new signals with increased expr
                    annotations="*", y_position = 25100, tip_length = 0.0, vjust = 0.5) +
  geom_signif(comparisons = list(c("OST1-MUT1", "HSP150-HSP150")),  #published ost1 vs new signals 
              annotations="***", y_position = 27100, tip_length = 0.0, vjust = 0.5) +
  geom_signif(comparisons = list(c("OST1-MUT1", "PIR1-PIR1")),  
              annotations="**", y_position = 28100, tip_length = 0.0, vjust = 0.5) 




#MODIFIFIED SIGNALS
mod <- read.xlsx("raw_mod.xlsx")
mod

mod_mu<- mod %>% rowwise() %>% mutate(avg=mean(c(r1,r2,r3)), se=std.error(c(r1,r2,r3))) #GOOD
mod_mu

#order samples by order of sds page
mod_mu$sample <-factor(mod_mu$sample, levels=unique(ord))

mplot <- ggplot(mod_mu, aes(x = sample, y = avg, fill=sample)) +
  geom_col(position = position_dodge(width = 0.5), color = "black", size = 0.8) +
  #changing x-axis labels from ord to split pre and pro signal schemes
  scale_x_discrete(labels = c("HSP150-HSP150" = "HSP150", "HSP150-MFA" = "MFA","HSP150-MUT1" = "MUT1", 
                              "PIR1-PIR1"="PIR1", "PIR1-MFA"="MFA", "PIR1-MUT1"="MUT1", 
                              "BP"="BP")) +
  scale_fill_manual(values=c("#9A91A1","#9A91A1","#9A91A1","#BFC5D9","#BFC5D9","#BFC5D9","white")) +
  geom_point(aes(y = r1), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r2), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r3), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 0.2) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 5000), expand = expansion(mult = c(0, 0.05))) +
  ##scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0, 1.75)) +
  xlab("Secretion Signal") +
  ylab("Pixel Coverage (Area)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        text = element_text(family = "Calibri Light"),
        legend.position = "none",
        axis.line = element_line(),
        axis.ticks = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1)) 

mplot


B<-mplot + geom_signif(comparisons = list(c("HSP150-HSP150", "HSP150-MFA")),  #hsp150 vs modified versions
                    annotations="***", y_position = 25100, tip_length = 0.0, vjust = 0.5) +
  geom_signif(comparisons = list(c("HSP150-HSP150", "HSP150-MUT1")),  #published ost1 vs new signals 
              annotations="**", y_position = 26100, tip_length = 0.0, vjust = 0.5) +
  geom_signif(comparisons = list(c("PIR1-PIR1", "PIR1-MFA")),  
              annotations="***", y_position = 25100, tip_length = 0.0, vjust = 0.5) +
  geom_signif(comparisons = list(c("PIR1-PIR1", "PIR1-MUT1")),  
              annotations="***", y_position = 26100, tip_length = 0.0, vjust = 0.5)

ggarrange(A,B, ncol = 2, common.legend = FALSE, align = "h")

A


####colour palette change####
nplot <- ggplot(nat_mu, aes(x = sample, y = avg, fill=sample)) +
  geom_col(position = position_dodge(width = 0.5, preserve = 'single'), color = "black", size = 0.8) +
  #changing x-axis labels from ord to split pre and pro signal schemes
  scale_x_discrete(labels = c("OST1-MUT1" = "OST1-MUT1", 
                              "MFA-MFA"="MFA", 
                              "HSP150-HSP150" = "HSP150", 
                              "PIR1-PIR1"="PIR1", 
                              "BP"="BP")) +
  scale_fill_manual(values=c("#E4F5E0","#A9DBB9","#7CA5B8","#466286","white")) +
  geom_point(aes(y = r1), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r2), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r3), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 0.2) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 5000), expand = expansion(mult = c(0, 0.05))) +
  ##scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0, 1.75)) +
  xlab("Secretion Signal") +
  ylab("Pixel Coverage (Area)") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        text = element_text(family = "Calibri Light"),
        legend.position = "none",
        axis.line = element_line(),
        axis.ticks = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1)) 

C <-nplot + geom_signif(comparisons = list(c("MFA-MFA", "HSP150-HSP150")),  #mfa-mfa vs new signals with increased expr
                       annotations="*", y_position = 25600, tip_length = 0.0, vjust = 0.5) +
  geom_signif(comparisons = list(c("OST1-MUT1", "HSP150-HSP150")),  #published ost1 vs new signals 
              annotations="***", y_position = 27100, tip_length = 0.0, vjust = 0.5) +
  geom_signif(comparisons = list(c("OST1-MUT1", "PIR1-PIR1")),  
              annotations="**", y_position = 28600, tip_length = 0.0, vjust = 0.5) +
  coord_cartesian(ylim = c(0, 30000))  # Match y-axis limits to Plot D



 mplot <- ggplot(mod_mu, aes(x = sample, y = avg, fill=sample)) +
  geom_col(position = position_dodge(width = 0.5, preserve = 'single'), color = "black", size = 0.8) +
  #changing x-axis labels from ord to split pre and pro signal schemes
  scale_x_discrete(labels = c("HSP150-HSP150" = "HSP150", "HSP150-MFA" = "MFA","HSP150-MUT1" = "MUT1", 
                              "PIR1-PIR1"="PIR1", "PIR1-MFA"="MFA", "PIR1-MUT1"="MUT1", 
                              "BP"="BP")) +
  scale_fill_manual(values=c("#7CA5B8","#7CA5B8","#7CA5B8","#466286","#466286","#466286","white")) +
  geom_point(aes(y = r1), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r2), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_point(aes(y = r3), shape = 1, position = position_dodge(width = 0.5), size = 1.5) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se), width = 0.2) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 5000), expand = expansion(mult = c(0, 0.05))) +
  ##scale_fill_viridis_c(option = "rocket", direction = -1, limits = c(0, 1.75)) +
  xlab("Pro Signal") +
  ylab("Pixel Coverage (Area)") +
  labs(fill="Pre Signal") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        text = element_text(family = "Calibri Light"),
        legend.position = "none",
        axis.line = element_line(),
        axis.ticks = element_line(),
        axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1)) 

D <- mplot + geom_signif(comparisons = list(c("HSP150-HSP150", "HSP150-MFA")),  #hsp150 vs modified versions
                       annotations="***", y_position = 25600, tip_length = 0.0, vjust = 0.5) +
  geom_signif(comparisons = list(c("HSP150-HSP150", "HSP150-MUT1")),  #published ost1 vs new signals 
              annotations="**", y_position = 27100, tip_length = 0.0, vjust = 0.5) +
  geom_signif(comparisons = list(c("PIR1-PIR1", "PIR1-MFA")),  
              annotations="***", y_position = 25600, tip_length = 0.0, vjust = 0.5) +
  geom_signif(comparisons = list(c("PIR1-PIR1", "PIR1-MUT1")),  
              annotations="***", y_position = 27100, tip_length = 0.0, vjust = 0.5) +
  coord_cartesian(ylim = c(0, 30000))  # Match y-axis limits to Plot C


D


ggarrange(C,D, ncol = 2, align = "h") #in thesis
library(patchwork)
combined_plots <- C & D / plot_layout(guides = "collect")
combined_plots





