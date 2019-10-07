################# Packages needed ################# 

library('vegan')
library(dplyr)
library(car)
library(psych)
library(tidyverse)
library(reshape)

# for plots : alpha diversity, PCoA,
library(ggplot2)
library(ggpubr)

# set to directory with all files
setwd("~/Dropbox/projects/Finch1/manuscripts/Finch1_sampleComparison/dataAnalysis_clean/to publish")
################################~#################################### 

################# relative abundance setup #################
# table from taxa-barplots.qzv, exported .csv at class level from https://view.qiime2.org/
# read in table without metadata at the end, taxa are columns
abun <- (read.csv("level-3.csv",
                  header = TRUE,
                  row.names = 1))[1:18]

  # check which have unassigned and rename manually otherwise sub will make them blank
colnames(abun)
colnames(abun)[9] <- "Epsilonproteobacteria"
colnames(abun)[17] <- "Unassigned Bacteria"
colnames(abun)[18] <- "Unassigned"
colnames(abun) <- sub(".*_", "", colnames(abun))

# Change to proportions
for (i in 1:(length(rownames(abun)))) {
  total <- sum(abun[i,])
  for (j in 1:(length(colnames(abun)))) {
  abun[i,j]  <- (abun[i,j])/total
  }
}
# Check that rows add up to 1
sum(abun[5,])


# find out which are <1% for all samples
abun_maxmin <- NULL
for (i in 1:(length(colnames(abun)))) {
  abun_maxmin <- bind_rows(abun_maxmin, 
                           data.frame("class" = colnames(abun)[i],
                                      "min" = min(abun[,i]),
                                      "max" = max(abun[,i])))
}
cbind(abun_maxmin$class, colnames(abun)) # make sure taxa in same order still
onePerc <- which(abun_maxmin$max>=0.01) #list of column >1%
gr1 <- colnames(abun)[onePerc] # list of taxa >1%, for key later

# pulling only columns >1%
abun <- abun[,onePerc]

# make rownames 1st row (can't have duplicate row name, important for later)
abun <- cbind(rownames(abun), data.frame(abun, row.names=NULL))

# melt for easy plotting
abun <- melt(abun)

# replace sample IDs with only sample type
abun[,1] <- sub("T[0-9]*.pv", "proventriculus", abun[,1])
abun[,1] <- sub("T[0-9]*.sm", "small int.", abun[,1])
abun[,1] <- sub("T[0-9]*.lg", "large int.", abun[,1])
abun[,1] <- sub("T[0-9]*.cs", "cloacal swab", abun[,1])
abun[,1] <- sub("T[0-9]*.f", "feces", abun[,1])

colnames(abun) <- c("X", "variable", "value")


################# relative abundance plots - col. & g.s. #################

### RELATIVE ABUNDANCE PLOT - COLOR

ggplot(abun,aes(x = X, y = value,fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") +
  theme_classic() +
  scale_fill_manual(name = "Class", 
                    values = c("#953272",
                               "#1E2085",
                               "#005B94",
                               "#9CCADE",
                               "#00AD00",
                               "#B2FF2E",
                               "#DB6A00",
                               "#D12600")) +
  scale_x_discrete(limits=c("proventriculus",
                            "small int.",
                            "large int.",
                            "cloacal swab",
                            "feces")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15,
                                face="bold")) +
  theme(legend.text = element_text(size=14))+
  theme(legend.background = element_rect(size=0.5, 
                                         linetype="solid", 
                                         colour ="black")) +  
  labs(x="Sample Type",
       y="") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title=element_text(size=14), 
        legend.text=element_text(size=10)) 



### RELATIVE ABUNDANCE PLOT - GRAYSCALE

ggplot(abun,aes(x = X, y = value,fill = variable)) + 
  geom_bar(position = "fill",stat = "identity") +
  theme_classic() +
  scale_fill_grey(name = "Class", start = 0.8, end = 0) +
  scale_x_discrete(limits=c("proventriculus",
                            "small int.",
                            "large int.",
                            "cloacal swab",
                            "feces")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=15,
                                face="bold")) +
  theme(legend.text = element_text(size=14))+
  theme(legend.background = element_rect(size=0.5, 
                                         linetype="solid", 
                                         colour ="black")) +  
  labs(x="Sample Type",
       y="") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.title=element_text(size=14), 
        legend.text=element_text(size=10)) 

################################~#################################### 

################# alpha diversity file setup #################

adiv <- read.csv("alpha-diversity-shannon.csv")

levels(adiv$sample)

  # reorder so they make more sense on figures
adiv$sample <- ordered(adiv$sample,
                       levels = c("pv", "sm",  "lg", "cs", "f"))

################# alpha diversity ANOVA & post-hoc #################

describeBy(adiv, group = "sample", mat = FALSE)
adiv.aov <- aov(shannon ~ sample, data = adiv)

# Checking ANOVA test assumptions
  # 1. Normality
res <- adiv.aov$residuals
shapiro.test(res)
  # 2. Homogeneity of variances
bartlett.test(shannon ~ sample, data = adiv)
plot(adiv.aov, 2)       

summary(adiv.aov)
TukeyHSD(adiv.aov)


################# alpha diversity figure - color & grayscale #################

# significantly different relationships from tukey post-hoc
my_comparisons = list( c("pv", "sm"),
                       c("pv", "lg"),
                       c("pv", "cs"),
                       c("pv", "f"))

# COLOR VERSION
ggplot(adiv, aes(x = sample,
                 y = shannon,
                 fill = sample)) +
  geom_boxplot() + 
  theme_classic() +
  scale_x_discrete(limits=c("pv",
                            "sm",
                            "lg",
                            "cs",
                            "f"),
                   labels=c("pv" = "proventriculus",
                            "sm" = "small int",
                            "lg" = "large int",
                            "cs" = "cloaca",
                            "f" = "feces")) +  
  scale_fill_manual(values = c("#D12600", "#DB6A00", "#B2FF2E", "#00AD00", "#005B94")) +
  theme(axis.text.x = element_text(angle = 30,
                                   vjust = .5)) +
  theme(legend.position = "none") +  
  labs(x="Sample Type",
       y="Shannon Index") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=20,
                                face="bold")) +
  stat_compare_means(label = "p.signif", 
                     method = "t.test",
                     comparisons = my_comparisons,
                     hide.ns = TRUE)    

# GRAYSCALE VERSION
ggplot(adiv, aes(x = sample,
                 y = shannon,
                 fill = sample)) +
  geom_boxplot() + 
  theme_classic() +
  scale_fill_grey(start = 0, end = 0.9) +
  scale_x_discrete(limits=c("pv",
                            "sm",
                            "lg",
                            "cs",
                            "f"),
                   labels=c("pv" = "proventriculus",
                            "sm" = "small int",
                            "lg" = "large int",
                            "cs" = "cloaca",
                            "f" = "feces")) +
  theme(axis.text.x = element_text(angle = 30,
                                   vjust = .5)) +
  theme(legend.position = "none") +  
  labs(x="Sample Type",
       y="Shannon Index") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=20,
                                face="bold")) +
  stat_compare_means(label = "p.signif", 
                     method = "t.test",
                     comparisons = my_comparisons,
                     hide.ns = TRUE)    


#############################~#################################### 

################# adonis file setup #################

# unweighted unifrac distance matrix
uwU <- read.csv("unweighted-distance-matrix.csv",
               sep=",",
               header=TRUE)
# weighted unifrac distance matrix
wU <- read.csv("weighted-distance-matrix.csv",
              sep=",",
              header=TRUE)
# map/metadata file
  # row names should be same as distance matricies, and in same order
map <- read.table("map2.txt", 
                 sep="\t", 
                 header=TRUE, 
                 row.names = 1)

#check all are in same order
cbind(colnames(uwU[-c(1)]), 
      colnames(wU[-c(1)]), 
      row.names(map))


uwU <- uwU[-c(1)] # make square
uwU <- data.matrix(uwU) #distance matrix format
wU <- wU[-c(1)] # make square
wU <- data.matrix(wU) #distance matrix format


################# setup for sub-comparisons #################

  # pulling specific sets - weighted uniFrac

fec_lg_wU_rows <- (map$sample=="f"|
                     map$sample=="lg")
fec_lg_wU_map <- map[fec_lg_wU_rows,]
fec_lg_wU <- data.matrix(wU[fec_lg_wU_rows,fec_lg_wU_rows])

fec_sm_wU_rows <- (map$sample=="f"|
                     map$sample=="sm")
fec_sm_wU_map <- map[fec_sm_wU_rows,]
fec_sm_wU <- data.matrix(wU[fec_sm_wU_rows,fec_sm_wU_rows])

fec_pv_wU_rows <- (map$sample=="f"|
                     map$sample=="pv")
fec_pv_wU_map <- map[fec_pv_wU_rows,]
fec_pv_wU <- data.matrix(wU[fec_pv_wU_rows,fec_pv_wU_rows])

fec_all_wU_rows <- (map$sample!="cs")
fec_all_wU_map <- map[fec_all_wU_rows,]
fec_all_wU <- data.matrix(wU[fec_all_wU_rows,fec_all_wU_rows])

clo_lg_wU_rows <- (map$sample=="cs"|
                     map$sample=="lg")
clo_lg_wU_map <- map[clo_lg_wU_rows,]
clo_lg_wU <- data.matrix(wU[clo_lg_wU_rows,clo_lg_wU_rows])

clo_sm_wU_rows <- (map$sample=="cs"|
                     map$sample=="sm")
clo_sm_wU_map <- map[clo_sm_wU_rows,]
clo_sm_wU <- data.matrix(wU[clo_sm_wU_rows,clo_sm_wU_rows])

clo_pv_wU_rows <- (map$sample=="cs"|
                     map$sample=="pv")
clo_pv_wU_map <- map[clo_pv_wU_rows,]
clo_pv_wU <- data.matrix(wU[clo_pv_wU_rows,clo_pv_wU_rows])

clo_all_wU_rows <- (map$sample!="f")
clo_all_wU_map <- map[clo_all_wU_rows,]
clo_all_wU <- data.matrix(wU[clo_all_wU_rows,clo_all_wU_rows])

  # pulling specific sets - unweighted

fec_lg_uwU_rows <- (map$sample=="f"|
                      map$sample=="lg")
fec_lg_uwU_map <- map[fec_lg_uwU_rows,]
fec_lg_uwU <- data.matrix(uwU[fec_lg_uwU_rows,fec_lg_uwU_rows])

fec_sm_uwU_rows <- (map$sample=="f"|
                      map$sample=="sm")
fec_sm_uwU_map <- map[fec_sm_uwU_rows,]
fec_sm_uwU <- data.matrix(uwU[fec_sm_uwU_rows,fec_sm_uwU_rows])

fec_pv_uwU_rows <- (map$sample=="f"|
                      map$sample=="pv")
fec_pv_uwU_map <- map[fec_pv_uwU_rows,]
fec_pv_uwU <- data.matrix(uwU[fec_pv_uwU_rows,fec_pv_uwU_rows])

fec_all_uwU_rows <- (map$sample!="cs")
fec_all_uwU_map <- map[fec_all_uwU_rows,]
fec_all_uwU <- data.matrix(uwU[fec_all_uwU_rows,fec_all_uwU_rows])

clo_lg_uwU_rows <- (map$sample=="cs"|
                      map$sample=="lg")
clo_lg_uwU_map <- map[clo_lg_uwU_rows,]
clo_lg_uwU <- data.matrix(uwU[clo_lg_uwU_rows,clo_lg_uwU_rows])

clo_sm_uwU_rows <- (map$sample=="cs"|
                      map$sample=="sm")
clo_sm_uwU_map <- map[clo_sm_uwU_rows,]
clo_sm_uwU <- data.matrix(uwU[clo_sm_uwU_rows,clo_sm_uwU_rows])

clo_pv_uwU_rows <- (map$sample=="cs"|
                      map$sample=="pv")
clo_pv_uwU_map <- map[clo_pv_uwU_rows,]
clo_pv_uwU <- data.matrix(uwU[clo_pv_uwU_rows,clo_pv_uwU_rows])

clo_all_uwU_rows <- (map$sample!="f")
clo_all_uwU_map <- map[clo_all_uwU_rows,]
clo_all_uwU <- data.matrix(uwU[clo_all_uwU_rows,clo_all_uwU_rows])

################# adonis comparisons #################
# UNWEIGHTED

set.seed(1990)
adonis(clo_lg_uwU ~ location, clo_lg_uwU_map, permutations = 999)

set.seed(1990)
adonis(clo_sm_uwU ~ location, clo_sm_uwU_map, permutations = 999)

set.seed(1990)
adonis(clo_pv_uwU ~ location, clo_pv_uwU_map, permutations = 999)

set.seed(1990)
adonis(clo_all_uwU ~ location, clo_all_uwU_map, permutations = 999)

set.seed(1990)
adonis(fec_lg_uwU ~ location, fec_lg_uwU_map, permutations = 999)

set.seed(1990)
adonis(fec_sm_uwU ~ location, fec_sm_uwU_map, permutations = 999)

set.seed(1990)
adonis(fec_pv_uwU ~ location, fec_pv_uwU_map, permutations = 999)

set.seed(1990)
adonis(fec_all_uwU ~ location, fec_all_uwU_map, permutations = 999)

# WEIGHTED


set.seed(1990)
adonis(clo_lg_wU ~ location, clo_lg_wU_map, permutations = 999)

set.seed(1990)
adonis(clo_sm_wU ~ location, clo_sm_wU_map, permutations = 999)

set.seed(1990)
adonis(clo_pv_wU ~ location, clo_pv_wU_map, permutations = 999)

set.seed(1990)
adonis(clo_all_wU ~ location, clo_all_wU_map, permutations = 999)

set.seed(1990)
adonis(fec_lg_wU ~ location, fec_lg_wU_map, permutations = 999)

set.seed(1990)
adonis(fec_sm_wU ~ location, fec_sm_wU_map, permutations = 999)

set.seed(1990)
adonis(fec_pv_wU ~ location, fec_pv_wU_map, permutations = 999)

set.seed(1990)
adonis(fec_all_wU ~ location, fec_all_wU_map, permutations = 999)


################################~#################################### 

################# PCoA file setup #################

unweighted_pc=read.table("unweighted_unifrac_pc.txt",header=TRUE)
samp_unweight <- (unweighted_pc$samp)
#ordering for key
levels(unweighted_pc$samp)
unweighted_pc$samp_unweight <- factor(unweighted_pc$samp, 
                                      levels = c("proventriculus", "small", "large", "cloaca", "feces"))

weighted_pc=read.table("weighted_unifrac_pc.txt",header=TRUE)
samp_weight <- (weighted_pc$samp)
#ordering for key
levels(weighted_pc$samp)
weighted_pc$samp_weight <- factor(weighted_pc$samp, 
                                  levels = c("proventriculus", "small", "large", "cloaca", "feces"))

################# Color PCoA plots #################
 
# color unweighted uniFrac plot
unweight <- ggplot(unweighted_pc, 
       aes(PCoA1, PCoA2, 
           shape = factor(samp_unweight),
           colour = factor(samp_unweight))) +
  geom_point(size = 3) + 
  scale_shape_manual(values=c(15, 18, 16, 2, 1)) +
  scale_color_manual(values = c("#D12600", "#DB6A00", "#B2FF2E", 
                                "#00AD00", "#005B94"),
                     guide = guide_legend(override.aes = list(linetype = c(rep("blank"))))) +
  labs(x="PCoA1 (34.17%)",
       y="PCoA2 (26.86%)",
       color=NULL,
       shape=NULL) + 
  theme_classic(base_size = 18) +
  ggtitle("Unweighted Unifrac") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(type='t',
               level = .90)

# color weighted uniFrac plot
weight <- ggplot(weighted_pc, 
       aes(PCoA1, PCoA2, 
           shape = factor(samp_weight),
           colour = factor(samp_weight))) +
  geom_point(size = 3) + 
  scale_shape_manual(values=c(15, 18, 16, 2, 1)) +
  scale_color_manual(values = c("#D12600", "#DB6A00", "#B2FF2E", 
                                "#00AD00", "#005B94"),
                     guide = guide_legend(override.aes = list(linetype = c(rep("blank"))))) +
  labs(x="PCoA1 (58.56%)",
       y="PCoA2 (19.53%)",
       color=NULL,
       shape=NULL) + 
  theme_classic(base_size = 18) +
  ggtitle("Weighted Unifrac") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(type='t',
               level = .90)

# combine plots & view
ggarrange(unweight, 
          weight, 
          labels = c("a)", "b)"),
          common.legend = TRUE,
          legend="bottom")



################# Grayscale PCoA plots #################

# grayscale unweighted uniFrac plot 
unweight <- ggplot(unweighted_pc, 
       aes(PCoA1, PCoA2,
           shape = factor(samp_unweight),
           colour = factor(samp_unweight))) +
  scale_colour_grey(start = 0.85, end = 0,
                    guide = guide_legend(override.aes = list(linetype = c(rep("blank"))))) +
  geom_point(size = 4) +
  scale_shape_manual(values=c(15, 18, 16, 2, 1)) +
  labs(x="PCoA1 (34.17%)",
       y="PCoA2 (26.86%)",
       color=NULL,
       shape=NULL) +
  theme_classic(base_size = 18) +
  ggtitle("Unweighted Unifrac") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(type='t',
               level = .90)

# grayscale weighted uniFrac plot 
weight <- ggplot(weighted_pc, 
       aes(PCoA1, PCoA2,
           shape = factor(samp_weight),
           colour = factor(samp_weight))) +
  scale_colour_grey(start = 0.85, end = 0,
                    guide = guide_legend(override.aes = list(linetype = c(rep("blank"))))) +
  geom_point(size = 4) +
  scale_shape_manual(values=c(15, 18, 16, 2, 1)) +
  labs(x="PCoA1 (58.56%)",
       y="PCoA2 (19.53%)",
       color=NULL,
       shape=NULL) + 
  theme_classic(base_size = 18) +
  ggtitle("Weighted Unifrac") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(type='t',
               level = .90)

# combine plots & view
ggarrange(unweight, 
          weight, 
          labels = c("a)", "b)"),
          common.legend = TRUE,
          legend="bottom")


################################~#################################### 

################# Setup - paired t-tests ###############

uwU = read.csv("unweighted-distance-matrix.csv",
               row.names = 1)

wU = read.csv("weighted-distance-matrix.csv",
               row.names = 1)

# check files in same order
cbind(rownames(map),rownames(uwU),rownames(wU))

# specifying which birds for functions (all with feces)
fec_birds <- c("T4","T6","T8","T32","T35","T37","T38","T39","T44")

# specifying which birds for functions (all with cloacalSwab)
clo_birds <- c("T4","T8","T32","T37","T39")

#dataframe where all the results will go
ttestResults <- (matrix(c("group1", "group2", "distance", "t", "df", "p"), nrow = 1))

################# Within vs between - paired t-tests - unweighted uniFrac ###############

# function to pull correct pairwise distances - WITHIN BIRD VS BETWEEN BIRDS
# (non-lethal sample compared to own bird's lg intestine, 
# vs non-lethal sample compared to average of all other birds' large intestine)
prws_dsts <- function(bird = NA, dat = NA, map = NA, refloc = NA, altloc = NA){
  # pull the parts of the matrix we want
  # focal bird reference location
  focalsamp <- which(map$bird==bird & map$location==refloc)
  # focal bird alternate location
  altsamebird <- which(map$location==altloc & map$bird==bird)
  # alternate location, all other birds
  othersamp <- which(map$location==altloc & map$bird!=bird)
  
  # pulling single point (distance of self comparison for refloc vs altloc)
  self <- dat[focalsamp,altsamebird]
  # pulling all rows/columns from dm for among and averaging, excluding focal bird's altloc
  # (distance of among comparison for focal refloc vs all other bird's altloc)
  # averaging first row without focal lgint
  other <- mean(as.matrix(dat[c(focalsamp,othersamp),c(focalsamp,othersamp)])[1,-1])
  
  return(c(self,other))
}

## WITHIN BIRD VS BETWEEN BIRDS - FECES & LG INTESTINE ##

# make array of all within (column 1) and all between (column 2) using function
large_f_uw <- array(dim=c(length(fec_birds),2))
for(i in (1:length(fec_birds))){
  large_f_uw[i,] <- prws_dsts(fec_birds[i],dat=uwU,map=map,refloc="large_int",altloc="fecal")
}
# paired t-test of column 1 and column 2 of large_f_uw
tRes <- t.test(large_f_uw[,1],large_f_uw[,2],paired=TRUE)

ttestResults <- rbind(ttestResults,
                      c("fec-lg within", 
                        "fec-lg between", 
                        "unweighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))


## WITHIN BIRD VS BETWEEN BIRDS - CLOACAL SWAB & LG INTESTINE ##

# make array of all within (column 1) and all between (column 2) using function
large_c_uw <- array(dim=c(length(clo_birds),2))
for(i in (1:length(clo_birds))){
  large_c_uw[i,] <- prws_dsts(clo_birds[i],dat=uwU,map=map,refloc="large_int",altloc="cloacal_swab")
}
# paired t-test of column 1 and column 2 of large_c_uw
tRes <- t.test(large_c_uw[,1],large_c_uw[,2],paired=TRUE)
as.character(tRes[[9]])
ttestResults <- rbind(ttestResults,
                      c("clo-lg within", 
                        "clo-lg between", 
                        "unweighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))


################# Within vs between - paired t-tests - weighted uniFrac ###############

## FECES & LG INTESTINE ##

# make array of all within (column 1) and all between (column 2) using function
large_f_w <- array(dim=c(length(fec_birds),2))
for(i in (1:length(fec_birds))){
  large_f_w[i,] <- prws_dsts(fec_birds[i],dat=wU,map=map,refloc="large_int",altloc="fecal")
}
# paired t-test of column 1 and column 2 of large_f_w
tRes <- t.test(large_f_w[,1],large_f_w[,2],paired=TRUE)

ttestResults <- rbind(ttestResults,
                      c("fecal-large int within", 
                        "fecal-large int between", 
                        "weighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))



## CLOACAL SWAB & LG INTESTINE ##

# make array of all within (column 1) and all between (column 2) using function
large_c_w <- array(dim=c(length(clo_birds),2))
for(i in (1:length(clo_birds))){
  large_c_w[i,] <- prws_dsts(clo_birds[i],dat=wU,map=map,refloc="large_int",altloc="cloacal_swab")
}
# paired t-test of column 1 and column 2 of large_c_w
tRes <- t.test(large_c_w[,1],large_c_w[,2],paired=TRUE)

ttestResults <- rbind(ttestResults,
                      c("cloacal swab-large int within", 
                        "cloacal swab-large int between", 
                        "weighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))



################# Lethal sample comparisons - paired t-tests - unweighted uniFrac #######

# Same function as before,
# but  comparing distances from non-lethal to lethal 1 & 2
alltoall <- function(dat=NA, map=NA, refloc="non-lethal", altloc1="lethal-1", altloc2="lethal-2"){
  library(reshape)
  focalsamp <- which(map$location==refloc)
  othersamp1 <- which(map$location==altloc1)
  othersamp2 <- which(map$location==altloc2)
  
  leth1 <- (melt(as.matrix(dat[c(focalsamp),c(othersamp1)])))[,3]
  leth2 <- (melt(as.matrix(dat[c(focalsamp),c(othersamp2)])))[,3]
  
  n <- max(length(leth1), length(leth2))
  length(leth1) <- n                      
  length(leth2) <- n
  
  return(cbind(leth1,leth2))
}


## DIST TO SM. INT. vs DIST TO LG. INT. - FECES ##

all_a <- alltoall(dat=uwU, map=map, refloc="fecal", altloc1="small_int", altloc2="large_int")
#run paired t-test of column 1 and column 2 of all_a
tRes <- t.test(all_a[,1],all_a[,2])

ttestResults <- rbind(ttestResults,
                      c("fecal vs small int", 
                        "fecal vs large int", 
                        "unweighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))


## DIST TO SM. INT. vs DIST TO LG. INT. - CLOACAL SWAB ##

all_a <- alltoall(dat=uwU, map=map, refloc="cloacal_swab", altloc1="small_int", altloc2="large_int")
#run paired t-test of column 1 and column 2 of all_a
tRes <- t.test(all_a[,1],all_a[,2])

ttestResults <- rbind(ttestResults,
                      c("cloacal swab vs small int", 
                        "cloacal swab vs large int", 
                        "unweighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))

## DIST TO PROV. vs DIST TO LG. INT. - FECES ##

all_a <- alltoall(dat=uwU, map=map, refloc="fecal", altloc1="proventriculus", altloc2="large_int")
#run paired t-test of column 1 and column 2 of all_a
tRes <- t.test(all_a[,1],all_a[,2])

ttestResults <- rbind(ttestResults,
                      c("fecal vs proventriculus", 
                        "fecal vs large int", 
                        "unweighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))


## DIST TO PROV. vs DIST TO LG. INT. - CLOACAL SWAB ##

all_a <- alltoall(dat=uwU, map=map, refloc="cloacal_swab", altloc1="proventriculus", altloc2="large_int")
#run paired t-test of column 1 and column 2 of all_a
tRes <- t.test(all_a[,1],all_a[,2])

ttestResults <- rbind(ttestResults,
                      c("cloacal swab vs proventriculus", 
                        "cloacal swab vs large int", 
                        "unweighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))


## DIST TO PROV. vs DIST TO SM. INT. - FECES ##

all_a <- alltoall(dat=uwU, map=map, refloc="fecal", altloc1="proventriculus", altloc2="small_int")
#run paired t-test of column 1 and column 2 of all_a
tRes <- t.test(all_a[,1],all_a[,2])

ttestResults <- rbind(ttestResults,
                      c("fecal vs proventriculus", 
                        "fecal vs small int", 
                        "unweighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))


## DIST TO PROV. vs DIST TO SM. INT. - CLOACAL SWAB ##

all_a <- alltoall(dat=uwU, map=map, refloc="cloacal_swab", altloc1="proventriculus", altloc2="small_int")
#run paired t-test of column 1 and column 2 of all_a
tRes <- t.test(all_a[,1],all_a[,2])

ttestResults <- rbind(ttestResults,
                      c("cloacal swab vs proventriculus", 
                        "cloacal swab vs small int", 
                        "unweighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))


################# Lethal sample comparisons - paired t-tests - weighted uniFrac #######

## DIST TO SM. INT. vs DIST TO LG. INT. - FECES ##

all_a <- alltoall(dat=wU, map=map, refloc="fecal", altloc1="small_int", altloc2="large_int")
#run paired t-test of column 1 and column 2 of all_a
tRes <- t.test(all_a[,1],all_a[,2])

ttestResults <- rbind(ttestResults,
                      c("feces vs small int", 
                        "feces vs large int", 
                        "weighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))


## DIST TO SM. INT. vs DIST TO LG. INT. - CLOACAL SWAB ##

all_a <- alltoall(dat=wU, map=map, refloc="cloacal_swab", altloc1="small_int", altloc2="large_int")
#run paired t-test of column 1 and column 2 of all_a
tRes <- t.test(all_a[,1],all_a[,2])

ttestResults <- rbind(ttestResults,
                      c("cloacal swab vs small int", 
                        "cloacal swab vs large int", 
                        "weighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))

## DIST TO PROV. vs DIST TO LG. INT. - FECES ##

all_a <- alltoall(dat=wU, map=map, refloc="fecal", altloc1="proventriculus", altloc2="large_int")
#run paired t-test of column 1 and column 2 of all_a
tRes <- t.test(all_a[,1],all_a[,2])

ttestResults <- rbind(ttestResults,
                      c("feces vs proventriculus", 
                        "feces vs large int", 
                        "weighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))


## DIST TO PROV. vs DIST TO LG. INT. - CLOACAL SWAB ##

all_a <- alltoall(dat=wU, map=map, refloc="cloacal_swab", altloc1="proventriculus", altloc2="large_int")
#run paired t-test of column 1 and column 2 of all_a
tRes <- t.test(all_a[,1],all_a[,2])

ttestResults <- rbind(ttestResults,
                      c("cloacal swab vs proventriculus", 
                        "cloacal swab vs large int", 
                        "weighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))


## DIST TO PROV. vs DIST TO SM. INT. - FECES ##

all_a <- alltoall(dat=wU, map=map, refloc="fecal", altloc1="proventriculus", altloc2="small_int")
#run paired t-test of column 1 and column 2 of all_a
tRes <- t.test(all_a[,1],all_a[,2])

ttestResults <- rbind(ttestResults,
                      c("feces vs proventriculus", 
                        "feces vs small int", 
                        "weighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))


## DIST TO PROV. vs DIST TO SM. INT. - CLOACAL SWAB ##

all_a <- alltoall(dat=wU, map=map, refloc="cloacal_swab", altloc1="proventriculus", altloc2="small_int")
#run paired t-test of column 1 and column 2 of all_a
tRes <- t.test(all_a[,1],all_a[,2])

ttestResults <- rbind(ttestResults,
                      c("cloacal swab vs proventriculus", 
                        "cloacal swab vs small int", 
                        "weighted unifrac", 
                        as.numeric(tRes[[1]]), 
                        as.numeric(tRes[[2]]), 
                        as.numeric(tRes[[3]])))

      ##### SAVE ALL t-test RESULTS IN .csv
  write.csv(ttestResults, file = "ttestResults.csv")

###############################~##################################### 

################# all sample comparison boxplot - col. & g.s.################# 

uwU <- read.csv("unweighted-distance-matrix.csv",
                sep=",",
                header=TRUE)
# transform distance matrix to distances in one column and pairs in other columns
uwU <- melt(uwU)

# replace all sample IDs with only sample type (unweighted)
uwU[,1] <- sub("T[0-9]*.pv", "proventriculus", uwU[,1])
uwU[,1] <- sub("T[0-9]*.sm", "small", uwU[,1])
uwU[,1] <- sub("T[0-9]*.lg", "large", uwU[,1])
uwU[,1] <- sub("T[0-9]*.cs", "cloacal swab", uwU[,1])
uwU[,1] <- sub("T[0-9]*.f", "feces", uwU[,1])
uwU[,2] <- sub("T[0-9]*.pv", "proventriculus", uwU[,2])
uwU[,2] <- sub("T[0-9]*.sm", "small", uwU[,2])
uwU[,2] <- sub("T[0-9]*.lg", "large", uwU[,2])
uwU[,2] <- sub("T[0-9]*.cs", "cloacal swab", uwU[,2])
uwU[,2] <- sub("T[0-9]*.f", "feces", uwU[,2])
uwU <- cbind(uwU,"Unweighted UniFrac")

wU <- read.csv("weighted-distance-matrix.csv",
                sep=",",
                header=TRUE)
# transform distance matrix to distances in one column and pairs in other columns
wU <- melt(wU)

# replace all sample IDs with only sample type (weighted)
wU[,1] <- sub("T[0-9]*.pv", "proventriculus", wU[,1])
wU[,1] <- sub("T[0-9]*.sm", "small", wU[,1])
wU[,1] <- sub("T[0-9]*.lg", "large", wU[,1])
wU[,1] <- sub("T[0-9]*.cs", "cloacal swab", wU[,1])
wU[,1] <- sub("T[0-9]*.f", "feces", wU[,1])
wU[,2] <- sub("T[0-9]*.pv", "proventriculus", wU[,2])
wU[,2] <- sub("T[0-9]*.sm", "small", wU[,2])
wU[,2] <- sub("T[0-9]*.lg", "large", wU[,2])
wU[,2] <- sub("T[0-9]*.cs", "cloacal swab", wU[,2])
wU[,2] <- sub("T[0-9]*.f", "feces", wU[,2])
wU <- cbind(wU,"Weighted UniFrac")

# give both dataframes same/desired column names so we can combine
colnames(wU) <- c("focal", "compar", "values", "betaDiv")
colnames(uwU) <- c("focal", "compar", "values", "betaDiv")
# combine both dataframes
U <- rbind(uwU, wU)
# only rows with non-lethal samples as focal sample
U <- U[c(which(U$focal=="feces"|U$focal=="cloacal swab")),]
# only rows with lethal samples as comparison samples
U <- U[c(which(U$compar=="proventriculus"|U$compar=="small"|U$compar=="large")),]
# make factor and order levels so the colors don't get messed up when x-axis is plotted properly
U$compar <- factor(U$compar, levels = c("proventriculus", "small", "large"))


# plot all comparisons in 2x2 grid - color
ggplot(U, 
       aes(x = compar, y = values)) + 
  xlab("") +
  ylab("Distance") +
  # ylim(0,1.25) +
  geom_boxplot() +  
  scale_fill_manual(values = c("#D12600", "#DB6A00", "#B2FF2E")) +
  geom_boxplot(outlier.size = .1,
               outlier.shape = 20,
               aes(fill = compar)) + 
  theme_classic() +
  facet_grid(betaDiv ~ focal, scales="free") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("small", "large"),
                                        c("proventriculus", "small"),
                                        c("proventriculus", "large")),
                     label = "p.signif",
                     step_increase = .2) +
  theme(legend.position="none", 
        strip.text.x = element_text(size = 14),
        axis.text.x=element_text(size=18, 
                                 angle = 45,
                                 hjust=1), 
        strip.text.y = element_text(size = 20),
        axis.text.y=element_text(size=14),
        axis.title.y = element_text(size = 24))


# plot all comparisons in 2x2 grid - grayscale
ggplot(U, 
       aes(x = compar, y = values)) + 
  xlab("") +
  ylab("Distance") +
  # ylim(0,1.25) +
  geom_boxplot() +  
  scale_fill_grey(start = 0.3, end = .9) +
  geom_boxplot(outlier.size = .1,
               outlier.shape = 20,
               aes(fill = compar)) + 
  theme_classic() +
  facet_grid(betaDiv ~ focal, scales="free") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("small", "large"),
                                        c("proventriculus", "small"),
                                        c("proventriculus", "large")),
                     label = "p.signif",
                     step_increase = .2) +
  theme(legend.position="none", 
        strip.text.x = element_text(size = 14),
        axis.text.x=element_text(size=18, 
                                 angle = 45,
                                 hjust=1), 
        strip.text.y = element_text(size = 20),
        axis.text.y=element_text(size=14),
        axis.title.y = element_text(size = 24))





################# individual specificity boxplot - col. & g.s. #################
library(reshape)
# setup: 
  # using objects generated above for t-tests for individual specificity
large_c_uw <- cbind("cloacal swab", "Unweighted UniFrac", large_c_uw)
large_c_w <- cbind("cloacal swab", "Weighted UniFrac", large_c_w)
large_f_uw <- cbind("feces", "Unweighted UniFrac", large_f_uw)
large_f_w <- cbind("feces", "Weighted UniFrac", large_f_w)

all <- rbind(large_c_uw, large_c_w, large_f_uw, large_f_w)
colnames(all) <- c("sample", "betaDiv", "self", "other")

#give each distance value it's own row, with which pairwise comparison in own columns
all <- reshape::melt(as.data.frame(all), 
                     id = (c("betaDiv", "sample")))
# colnames(all) <- c("betaDiv", "sample", "variable", "values")
class(all$value)
# for some reason the values are a factor, which we cannot convert directly to a number so we must first make it characters then numbers. 
all$value <- as.numeric(as.character(all$value))



# boxplot - color
ggplot(all, 
       aes(x= all$variable, y = all$value)) + 
  xlab("") +
  ylab("Distance to large intestines") +
  geom_boxplot() +  
  geom_boxplot(outlier.size = .1,
               outlier.shape = 20,
               aes(fill = variable)) + 
  scale_fill_manual(values = c("#00AD00", "#005B94")) +
  theme_classic() +
  facet_grid(betaDiv ~ sample, scales = "free") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("self", "other")),
                     label = "p.signif",
                     paired = TRUE) +
  theme(legend.position="none", 
        strip.text.x = element_text(size = 14),
        axis.text.x=element_text(size=18, 
                                 angle = 45,
                                 hjust=1), 
        strip.text.y = element_text(size = 20),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=20))



# boxplot - grayscale
ggplot(all, 
       aes(x= all$variable, y = all$value)) + 
  xlab("") +
  ylab("Distance to large intestines") +
  geom_boxplot() +  
  geom_boxplot(outlier.size = .1,
               outlier.shape = 20,
               aes(fill = variable)) + 
  scale_fill_grey(start = 0.3, end = .9) +
  theme_classic() +
  facet_grid(betaDiv ~ sample, scales = "free") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("self", "other")),
                     label = "p.signif",
                     paired = TRUE) +
  theme(legend.position="none", 
        strip.text.x = element_text(size = 14),
        axis.text.x=element_text(size=18, 
                                 angle = 45,
                                 hjust=1), 
        strip.text.y = element_text(size = 20),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=20))


