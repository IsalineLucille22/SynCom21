#This script analyses SynCom Phase1 dynamics (Syn20 & Syn21) in order to investigate the effect of Lysobacter, as well as generaly assembly dynamics

#Load libraries
library(reshape2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(grid)
# install.packages("Hmisc")
library(Hmisc)
library(FME)
library(deSolve)
library(MASS)
library(permute)
# install.packages("vegan", type = "binary")
library(vegan)

#Set working directory
# setwd("/Users/isalinelucille-guex/Documents/CoCulture_Soil/Data/")
setwd("/Users/iguex/Documents/CoCulture_Soil/Data")

#Load data 
# relabund = read_excel("S20_S20Sim_rel_abund.xlsx", sheet = 6, range = "A1:AA85") #78,97,72,48. Put 6 to compute absolute abundances
relabund = read_excel("S20_S20Sim_rel_abund.xlsx", sheet = 1, range = "A1:AA78") #78, 85,97,72,48. Put 6 to compute absolute abundances
#Get mean of each replicate
relabund.m <- melt(relabund, id.vars = c(1:6))
relabund.m$Group <- paste0(relabund.m$condition, "_", relabund.m$time.c, "_", relabund.m$variable)
mean_rel <- relabund.m %>% 
  group_by(Group) %>%
  summarise(Mean = mean(value, na.rm = T))

#Format tables to add back metadat
mean_rel[c("condition","time.c","variable")] <- str_split_fixed(mean_rel$Group, "_", 3)

mean_rel$time.c <- as.factor(mean_rel$time.c)
mean_rel$time.c <- ordered(mean_rel$time.c, levels =c("0","1","3","7","10","21"))

###################################################################

#Plot stacked barplot showing relative abundances of each strain (Lysobacter in light green)
myCol = c("#b35806", "#e08214", "#d53e4f", "#b2182b", "#d6604d", 
          "#c51b7d","#de77ae", "#f1b6da", "#fdae61", "#fee090", 
          "#a6d96a", "#5aae61", "#01665e","#35978f","#1b7837", "#c2a5cf", "#9970ab", 
          "#762a83",  "#80cdc1", "#c7eae5",
          "#2166ac", "#4393c3")


#10x9
ggplot(mean_rel, aes(x = time.c, y = Mean, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal()+
  scale_fill_manual(values=c(myCol))+
  facet_wrap(~condition, scales = "free_x")


###################################################################

#Alpha diversity calculation: shannon takes into account evennes and richness for 
#diversity calculation; richness only counts amount of species detected

str(relabund)
#only numeric data can be used as input (no metadata)
richness <- as.data.frame(specnumber(relabund[,c(7:27)]))
temp_df <- as.data.frame(relabund[,c(7:27)])
shannondiv <- as.data.frame(diversity(temp_df, index="shannon"))
#merge data back
diversity <- cbind(shannondiv,richness)
colnames(diversity) <- c("shannon", "richness")
#add back metadat
diversity_meta <- cbind(relabund[,c(1:6)],diversity)

time_levels <- c("0","1","3","7","10","21")

# Convert time.c to a factor with specified levels
diversity_meta$time.c <- factor(diversity_meta$time.c, levels = time_levels)

# Plot data: shannon
ggplot(diversity_meta, aes(x = time.c, y = shannon, fill = condition)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_point(position = position_dodge(width = 0.8), size = 3, shape = 21) +
  scale_fill_manual(values = c("steelblue2", "maroon2", "darkseagreen3", "darkorange")) +
  theme_bw()

#Plot data:richness
#add line
diversity_meta$Group <- paste0(diversity_meta$condition, "_", diversity_meta$time.c)
rich_avg <- diversity_meta %>% 
  group_by(Group) %>%
  summarise(Mean = mean(richness, na.rm = T))
rich_avg[c("condition","time.c")] <- str_split_fixed(rich_avg$Group, "_", 2)
rich_avg$time <- as.numeric(rich_avg$time.c)
diversity_meta$time <- as.numeric(diversity_meta$time.c)

ggplot(diversity_meta, aes(x = time, y = richness, fill = condition)) +
  geom_point(position = position_jitter(width = 0.15), size = 3, shape = 21) +
  scale_fill_manual(values=c("steelblue2", "maroon2", "darkseagreen3", "darkorange"))+
  geom_line(data=rich_avg, aes(x = time, y = Mean, color = condition))+
  scale_color_manual(values=c("steelblue2", "maroon2", "darkseagreen3", "darkorange"))+
  ylim(15,22)+
  theme_bw()


###################################################################

#Beta diversity 
#NMDS can gewnerate different result each time, set seed
set.seed(234)

#For calculations metadata and numerical data need to be separated
relabund$Group <- paste0(relabund$condition, "_", relabund$time.c)
names <- relabund[,c(1:6,28)]
numbers <- relabund[,c(7:27)]
# num_list <- lapply(numbers, function(x) as.numeric(x))

#Calculate bray curtis distances
bray_curtis_dist <- vegdist(numbers, method = "bray")
#Run NMDS
nmds_result <- metaMDS(bray_curtis_dist)

#Check NMDs parameters
goodness(nmds_result) # Produces a results of test statistics for goodness of fit for each point
stressplot(nmds_result) # Produces a Shepards diagram
#As each time the result will vary slightly, I prefer to save at this point
saveRDS(nmds_result, "./NMDS_result.RDS")
saveRDS(bray_curtis_dist, "./bray_curtis_dist.RDS")

#Extract scores and points to calculate community centroids (these will be plotted as larger dots for
# easier visualisation)
#Get points
z.points <- data.frame(nmds_result$points)

z.points <- data.frame(z.points,
                       relabund$sample.id, relabund$condition, relabund$time.c, relabund$Group)
#Get scores
scrs <- scores(nmds_result, display = "sites", "species")
#Get centroid
cent <- aggregate(scrs ~ Group, data = names, FUN = "mean")
z.points$Group <- z.points$relabund.Group
merged <- merge(z.points, cent, by= c("Group"))

nmds_col <- c("lightskyblue1","#87CEFA","royalblue3", "#00008B","steelblue2","royalblue1",
              "pink", "#FF69B4","#8B008B","darkorchid4","hotpink2", "#C71585",
              "darkseagreen1", "darkseagreen2", "darkolivegreen", "darkgreen", "darkseagreen3", "darkseagreen4",
              "darkorange", "darkorange2", "chocolate4", "brown", "darkorange3", "darkorange4" )

ggplot(data = z.points, aes(x = MDS1, y = MDS2, color=relabund.Group)) +
  theme_bw() +
  geom_point(shape=16,size=4, alpha=0.7)+
  geom_point(data=cent, aes(x=NMDS1, y=NMDS2, color=Group), size= 5)+
  geom_segment(data = merged, 
               aes(xend = NMDS1, yend = NMDS2), 
               size = 0.5, 
               linetype = "dashed", 
               alpha = 0.5)+
  scale_color_manual(values=c(nmds_col))

