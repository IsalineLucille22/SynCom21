#Load libraries
library(reshape2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(vegan)
library(ggpubr)
library(readxl)
library(grid)
library(Hmisc)
library(FME)
library(deSolve)
library(MASS)
library(permute)


#Set working directory
setwd("/Users/isalinelucille-guex/Documents/CoCulture_Soil/Data/")

#Load data 
# absabund <- read.table("Updated_phase1_cfu&absabund_syn21&20.csv", header= TRUE, sep =",") #abundance data for all samples
absabund = read_excel("S20_S20Sim_rel_abund.xlsx", sheet = 5, range = "A1:AA169") #78, 48, 96, 97, 192

#Get mean of each replicate
absabund.m <- melt(absabund, id.vars = c(1:6))

absabund.m$Group <- paste0(absabund.m$condition, "_", absabund.m$time.c, "_", absabund.m$variable)

mean_abs <- absabund.m %>% 
  group_by(Group) %>%
  summarise(Mean = mean(value, na.rm = T))

#Get per species effect

abs_s20 <- absabund.m[absabund.m$condition =="S21",]
abs_s21 <- absabund.m[absabund.m$condition =="SAPost",]

abs_s20$Type <- paste0(abs_s20$Group,"_",abs_s20$replicate) 
abs_s20$Type <- gsub("S21","T",abs_s20$Type)
abs_s20$value_20 <- abs_s20$value
abs_s20$value <- NULL

abs_s21$Type <- paste0(abs_s21$Group,"_",abs_s21$replicate) 
abs_s21$Type <- gsub("SAPost","T",abs_s21$Type)
abs_s21$value_21 <- abs_s21$value
abs_s21$value <- NULL

merg_abs <- merge(abs_s20,abs_s21, by=c("Type"))
merg_abs$difference <- merg_abs$value_20 - merg_abs$value_21
merg_abs$variable.y <- as.factor(merg_abs$variable.y)
merg_abs$variable.y <- ordered(merg_abs$variable.y, levels=c("Bradyrhizobium", "Burkholderia", "Caulobacter", "Cellulomonas", "Chitinophaga",
                                                             "Cohnella", "Curtobacterium", "Devosia", "Flavobacterium", "Luteibacter", "Lysobacter",
                                                             "Mesorhizobium", "Microbacterium", "Mucilaginibacter", "Phenylobacterium", "Pseudomonas1_R2",
                                                             "Pseudomonas2_R2", "Rahnella", "Rhodococcus", "Tardiphaga", "Variovorax"))


myCol = c("#b35806", "#e08214", "#d53e4f", "#b2182b", "#d6604d", 
          "#c51b7d","#de77ae", "#f1b6da", "#fdae61", "#fee090", 
          "#a6d96a", "#5aae61", "#01665e","#35978f","#1b7837", "#c2a5cf", "#9970ab", 
          "#762a83",  "#80cdc1", "#c7eae5",
          "#2166ac", "#4393c3")

# Convert time.c.y to a factor with sorted levels
merg_abs$time.c.y <- factor(merg_abs$time.c.y, levels = c("0","1","3","7","10","21"))


#6x20
ggplot(merg_abs, aes(x=variable.y, y=difference, fill=variable.y))+
  geom_bar(stat = "summary", fun = "mean", width = 0.7)+
  geom_point(size=2, alpha=0.9)+
  coord_flip()+
  facet_wrap(~time.c.y, nrow=1)+
  scale_fill_manual(values=c(myCol))+
  ylim(-3.8e-04,3.8e-04)+
  theme_bw()

