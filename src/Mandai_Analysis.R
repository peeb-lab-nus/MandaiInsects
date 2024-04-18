

rm(list = ls())

library(sf)
library(sp)
library(dplyr)
library(ggplot2)
library(readxl)
library(reshape2)
library(vegan)
library(FD)
library(picante)
library(viridis)
library(plyr)
library(iNEXT)
library(tidyverse)
library(terra)

if(Sys.info()[["user"]] == "junyinglim"){
  main.dir <- "~/Dropbox/Projects/MandaiInsects"
  data.dir <- file.path(main.dir, "data")  
  fig.dir <- file.path(main.dir, "figures")  
}

## GRAPHICAL PARAMETERS ====================
berkcol3 = c("#003262", "#3B7EA1", "#9BBEA9", "#00B0DA", "#00A598",
             "#006D2C", "#CFDD45", "#859438", "#FDB515", "#FD8D3C",
             "#ED4E33", "#C4820E", "#D9661F", "#6C3302")

## IMPORT DATA ====================
MandaiTree <- vect(file.path(data.dir, "TreesMPD_202001.shp"))
MandaiCoords <- read.csv(file.path(data.dir, "MandaiCoordinates.csv"), header = TRUE)
MandaiInsects <- read_excel(file.path(data.dir, "Mandai_FullSeq_May2019only_Clustered_JP06Mar2024.xlsx"),
                            sheet = 2)

# Only analysing data from May
table(MandaiInsects$`Coll. Date`)
table(MandaiInsects[,c("trapid", "Coll. Date")])
table(MandaiInsects$`Trap Type`) 

## PROCESSING DATA ====================
MandaiCoords$Trap <- paste0("MIS-", MandaiCoords$Trap)
# Only analyse
MandaiInsects_MT <- MandaiInsects %>% subset(`Trap Type` == "MT")

table(MandaiInsects_MT[,c("trapid", "Coll. Date")]) # only 6 sites

# Calculate summary statistics for each 
MandaiInsects_MT_abund <- ddply(.data = MandaiInsects_MT,
                                .variables = .(trapid),
                                .fun = summarise,
                                insect_div = length(unique(`3pClusterID`)),
                                insect_abund = length(trapid))

# Create site x OTU matrix
insect_mat <- acast(formula = trapid ~ `3pClusterID`, data = MandaiInsects_MT,
                    fun.aggregate = length, fill = 0)

## SAMPLE RAREFACTION ====================
raremax <- min(rowSums(insect_mat))

MandaiInsects_MT_div <- data.frame(trapid = rownames(insect_mat),
                                   insect_div_rarefy = rarefy(insect_mat, sample = raremax))

insect_rarefy_curve <- rarecurve(insect_mat, sample = raremax, tidy = TRUE)

insect_rarefy_curve_plot <- ggplot(data = insect_rarefy_curve) + 
  geom_path(aes(y = Species, x = Sample, color = Site)) + 
  scale_colour_manual(values = berkcol3[c(1,3, 5, 7, 9, 13)]) +
  labs(x = "Number of sequenced individuals", y = "Rarefied OTU richness") +
  coord_fixed()
ggsave(insect_rarefy_curve_plot, 
       filename = file.path(fig.dir, "insect_rarefy_curve.png"),
       width = 5, height = 4)

insect_rarefy_bar_plot <- ggplot(data = MandaiInsects_MT_div) + 
  geom_bar(aes(y = insect_div_rarefy, x = trapid, fill = trapid), stat = "identity") +
  scale_fill_manual(values = berkcol3[c(1,3, 5, 7, 9, 13)]) + 
  labs(x = "Trap", y = "Rarefied OTU richness") 

ggsave(insect_rarefy_bar_plot, 
       filename = file.path(fig.dir, "insect_rarefy_bar.png"),
       width = 5, height = 5)

## CALCULATE TREE DIVERSITY ====================
MandaiBuffer <- MandaiCoords %>% 
  subset(Trap %in% unique(MandaiInsects_MT$trapid)) %>%
  terra::vect(geom = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84") %>%
  terra::buffer(width = 20)

MandaiTreeClean <- MandaiTree %>%
  as.data.frame() %>%
  subset(DBH_M >= 0.3) %>%
  subset(!(Lat == 0 | Long == 0)) %>%
  terra::vect(geom = c("Long", "Lat"))

MandaiTreeBufferIntersect <- terra::relate(MandaiTreeClean, MandaiBuffer, relation = "within")

# Create list of trees in buffer zones around each trap
TrapTreeBuffer <- list()
for(i in 1:ncol(MandaiTreeBufferIntersect)){
  # 
  TrapTreeBuffer[[i]] <- as.data.frame(MandaiTreeClean)[MandaiTreeBufferIntersect[,i],]
  TrapTreeBuffer[[i]]$trapid <- MandaiBuffer$Trap[i]
}
TrapTreeBuffer <- TrapTreeBuffer %>% 
  do.call("rbind", .)

# Calculate summary statistics of trees in buffer zones
MandaiTreeSummary <- ddply(.data = TrapTreeBuffer,
                           .variables = .(trapid),
                           .fun = summarise,
                           sd_dbh = sd(DBH_M, na.rm = TRUE),
                           sd_height = sd(Height_M, na.rm = TRUE),
                           max_height = quantile(Height_M, 0.95),
                           n_trees = length(Species),
                           tree_div = length(unique(Species)))

MandaiTreeMat <- acast(formula = trapid ~ Species, data = TrapTreeBuffer, fill = 0)
treemin <- min(rowSums(MandaiTreeMat))

MandaiTree_Rarefy <- data.frame(trapid = rownames(MandaiTreeMat),
                                tree_div_rarefy = rarefy(MandaiTreeMat, sample = treemin))

MandaiDataFinal <- MandaiInsects_MT_abund %>%
  left_join(MandaiInsects_MT_div) %>%
  left_join(MandaiTreeSummary) %>%
  left_join(MandaiTree_Rarefy)

ggplot(data = MandaiDataFinal) + geom_point(aes(y = insect_div_rarefy, x = sd_dbh))
ggplot(data = MandaiDataFinal) + geom_point(aes(y = insect_div_rarefy, x = max_height))
ggplot(data = MandaiDataFinal) + geom_point(aes(y = insect_div_rarefy, x = sd_height))

summary(lm(insect_div_rarefy~ tree_div_rarefy, data = MandaiDataFinal))
#summary(lm(insect_div_rarefy~ tree_div, data = MandaiDataFinal))
summary(lm(insect_div_rarefy~ sd_height, data = MandaiDataFinal))
summary(lm(insect_div_rarefy~ max_height, data = MandaiDataFinal))
summary(lm(scale(insect_div_rarefy)~ scale(n_trees), data = MandaiDataFinal))



target.col <- c("tree_div_rarefy", "sd_height", "max_height", "n_trees")

betas <- vector()
confint_lwr <- vector()
confint_upp <- vector()
for(i in target.col){
  mod <- lm(formula = paste0("scale(insect_div_rarefy) ~ scale(", i, ")"), data = MandaiDataFinal)
  betas[i] <- coefficients(mod)[2]
  confint_lwr[i] <- confint(mod)[2,1]
  confint_upp[i] <- confint(mod)[2,2]
}

res <- data.frame("Variable" = target.col,
                  "Coefficient" = betas,
                  "Confint_upper" = confint_upp,
                  "Confint_lower" = confint_lwr)

# Relationship between rarefied insect diversity and plant and vegetation variables
insect_plant_relate_plot <- ggplot(data = res) + geom_point(aes(y = Variable, x = Coefficient)) +
  geom_errorbarh(aes(xmin = Confint_lower, xmax = Confint_upper, y = Variable), linewidth = 0.5, height = 0) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", colour = "red", size = 1) +
  scale_x_continuous(breaks = round(seq(-1.8, 1.8, 0.6), 1), limits = c(-1.8, 1.8))

ggsave(insect_plant_relate_plot,
       filename = file.path(fig.dir, "insect_plant_relate.png"),
       width = 5, height = 5)

## CALCULATE DISSIMILARITIES ====================
# Define function to calculate dissimilarities
calc_dissim <- function(x){
  vegdist(x, method = "bray") %>% 
    as.matrix() %>%
    melt()
}


# Calculate raw insect dissimilarities
insect_raw_dist <- calc_dissim(insect_mat)

insect_raw_diss_plot <- ggplot(data = insect_raw_dist) + 
  geom_tile(aes(y = Var1, x = Var2, fill = value)) + 
  scale_fill_viridis() +
  coord_fixed() + 
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(legend.position = "bottom")

ggsave(insect_raw_diss_plot, 
       filename = file.path(fig.dir, "insect_raw_diss.png"),
       width = 5, height = 5)

# Perform dissimilarities of 100 rarefied insect matrices
insect_mat_random <- list()
for(i in 1:100){
  insect_mat_random[[i]] <- rrarefy(insect_mat, sample = raremax)  
}

insect_random_dist <- lapply(insect_mat_random, FUN = calc_dissim) %>%
  do.call("rbind", .) %>%
  ddply(.variables = .(Var1, Var2), .fun = summarise, mean_insect_dissim = mean(value))

insect_rarefy_diss_plot <- ggplot(data = insect_random_dist) + 
  geom_tile(aes(y = Var1, x = Var2, fill = mean_insect_dissim)) + 
  scale_fill_viridis() +
  coord_fixed() +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  theme(legend.position = "bottom")

ggsave(insect_rarefy_diss_plot, 
       filename = file.path(fig.dir, "insect_rarefy_diss.png"),
       width = 5, height = 5)

# Calculate dissimilarities of 100 veg matrices
tree_mat_random <- list()
for(i in 1:100){
  tree_mat_random[[i]] <- rrarefy(MandaiTreeMat, sample = treemin)
}

tree_random_dist <- lapply(tree_mat_random, FUN = calc_dissim) %>%
  do.call("rbind", .) %>%
  ddply(.variables = .(Var1, Var2), .fun = summarise, mean_plant_dissim = mean(value))


rownames(MandaiDataFinal) <- MandaiDataFinal$trapid
veg_dist_list <- list()
for(i in target.col){
  veg_dist_list[[i]] <- dist(MandaiDataFinal[,i, drop = FALSE]) %>% 
    as.matrix() %>%
    melt(value.name = i)
}
veg_dist <- Reduce("left_join",veg_dist_list)

rownames(MandaiCoords) <- MandaiCoords$Trap

geog_dist <- MandaiCoords %>% 
  subset(Trap %in% unique(MandaiInsects_MT$trapid)) %>%
  dplyr::select(c("Latitude", "Longitude")) %>%
  dist() %>%
  as.matrix() %>%
  melt(value.name = "geog_dist")

MandaiDataCompFinal <- veg_dist %>%
  left_join(insect_random_dist, by = c("Var1", "Var2")) %>%
  left_join(tree_random_dist, by = c("Var1", "Var2")) %>%
  left_join(geog_dist, by = c("Var1", "Var2"))


insect_plant_dissim <- ggplot(data = MandaiDataCompFinal) + 
  geom_point(aes(y = mean_insect_dissim, x = mean_plant_dissim)) +
  labs(x = "Tree community dissimilarity", y = "Insect community dissimilarity")
ggsave(insect_plant_dissim,
       filename = file.path(fig.dir, "insect_plant_dissim.png"),
       width = 5, height = 5)

mantel(xdis = acast(Var1 ~ Var2, data = insect_random_dist),
       ydis = acast(Var1 ~ Var2, data = tree_random_dist))
mantel(xdis = acast(Var1 ~ Var2, data = insect_random_dist),
       ydis = acast(Var1 ~ Var2, data = geog_dist))
mantel.partial(xdis = acast(Var1 ~ Var2, data = insect_random_dist),
               ydis = acast(Var1 ~ Var2, data = tree_random_dist),
               zdis = acast(Var1 ~ Var2, data = geog_dist))


insect_geog_dissim <- ggplot(data = MandaiDataCompFinal) + 
  geom_point(aes(y = mean_insect_dissim, x = geog_dist)) +
  labs(x = "Geographic distance", y = "Insect community dissimilarity")
ggsave(insect_geog_dissim,
       filename = file.path(fig.dir, "insect_geog_dissim.png"),
       width = 5, height = 5)

mod <- lm(mean_insect_dissim ~ mean_plant_dissim + geog_dist, data = MandaiDataCompFinal)
crPlot(mod, variable = "mean_plant_dissim")
library(car)

library(cowplot)
plot_grid(insect_plant_dissim, insect_geog_dissim)


# ##only keeping traps with more than 100 specimens
# MandaiInsectsMay <- subset(MandaiInsectsMay, TrapID %in% c("L02","L03", "L04","L05","L06" ,"L07","L08","L10"))
# table(MandaiInsectsMay$TrapID)
# 
# 
# 
# 
# ##Malaise Trap only Data
# MandaiInsectsMayMT <- subset(MandaiInsectsMay, `Trap Type` %in% "MT")
# table(MandaiInsectsMayMT$TrapID)
# 
# MayTraps <- unique(MandaiInsectsMay$TrapID)
# MayTrapsMT <- unique(MandaiInsectsMayMT$TrapID)
# 
# #subset out coordinates for traps in each DF
# MandaiCoordsMay <- subset(MandaiCoords, Trap %in% MayTraps)
# MandaiCoordsMayMT <- subset(MandaiCoords, Trap %in% MayTrapsMT)
# 
# ##Creating a DF with only the coordinates of the trees and renaming columns 
# MandaiTree2 <- MandaiTree[,c(3,4)] 
# MandaiTree2 <- as.data.frame(MandaiTree2)
# MandaiTree2 <- MandaiTree2[,c(1,2)] 
# MandaiTree2$Latitude <- MandaiTree2$Lat
# MandaiTree2$Longitude <- MandaiTree2$Long
# MandaiTree2 <- MandaiTree2[,c(3,4)] 
# 
# ###define parameters and create empty vectors
# Trap = vector()
# #20m plot
# width = 0.00018
# Veg_df <- list()
# Latitude <- vector()
# Longitude <- vector()
# 
# ### getting tree data around each insect trap 
# for (i in 1:length(MandaiCoords$Trap)){
#   
#   MandaiTree3 <- vector()
#   MandaiTree3 <- as.data.frame(MandaiTree3)
#   Trap[i] = MandaiCoords$Trap[i]
#   Latitude[i] <- MandaiCoords$Latitude[i]
#   Longitude[i] <- MandaiCoords$Longitude[i]
#   new_row <- cbind(Latitude[i],Longitude[i])
#   colnames(new_row) <- c("Latitude","Longitude")
#   MandaiTree3 <- rbind(new_row,MandaiTree2)
#   Dist <- dist(MandaiTree3, method = "euclidean")
#   # Distances from the first point to all other points
#   distances_from_first_point <- as.matrix(Dist)[1, ]  
#   # Remove the distance of the first point to itself (which is 0)
#   distances_from_first_point <- distances_from_first_point[-1]
#   Treeindex <- which(distances_from_first_point <= width)
#   ##Selecting trees from main DF which are within 20m of the malaise trap
#   selected_trees <- MandaiTree[Treeindex, ]
#   Veg_df[[as.character(Trap[i])]] <- selected_trees
#   
# }
# 
# MandaiTrees <- do.call(rbind, Veg_df)
# MandaiTrees$Trap <- rownames(MandaiTrees)
# MandaiTrees$Trap <- gsub("\\.\\d+", "", MandaiTrees$Trap)
# rownames(MandaiTrees) <- NULL
# MandaiTrees <- st_drop_geometry(MandaiTrees)
# MandaiTrees <- as.data.frame(MandaiTrees) 
# 
# MandaiTrees<- subset(MandaiTrees, Trap!= "T02")
# MandaiTrees <- subset(MandaiTrees, Trap != "T03")
# 
# MandaiTrees$TrapID <- MandaiTrees$Trap
# MandaiTrees<- MandaiTrees[,-c(10)]
# 
# table(MandaiTrees$TrapID)
# ##Checking for any duplicates 
# length(unique(MandaiTrees$TreeID))
# duplicated(MandaiTrees,by = TreeID)
# 
# 
# ###Cleaning up Tree Data
# Smalltrees <- subset(MandaiTrees, DBH_M < 0.3)
# BigTrees <- subset(MandaiTrees, DBH_M >=0.3)
# 
# 
# BigTreesMay <- subset(BigTrees, TrapID %in% MayTraps)
# BigTreesMayMT <- subset(BigTrees, TrapID %in% MayTrapsMT)
# 
# ###Structural stats for trees 
# BigTreesDF <- BigTrees %>%
#   group_by(TrapID) %>%
#   group_map(~ data.frame(
#     TrapID = .y,
#     MinHeight95 = quantile(.x$Height_M, 0.05),
#     MaxHeight95 = quantile(.x$Height_M, 0.95),
#     SD_Height = sd(.x$Height_M)
#     )) %>%
#   bind_rows()
# row.names(BigTreesDF) <- NULL
# 
# 
# BigTreesMayDF <- subset(BigTreesDF, TrapID %in% MayTraps)
# BigTreesMayMTDF <- subset(BigTreesDF, TrapID %in% MayTrapsMT)
# 
# ######
# MandaiInsects$TrapID <- as.factor(MandaiInsects$TrapID)
# MandaiInsectsMay$TrapID <- as.factor(MandaiInsectsMay$TrapID)
# MandaiInsectsMayMT$TrapID <- as.factor(MandaiInsectsMayMT$TrapID)
# 
# #########
# #Matrices 
# 
# BigTreesMatrix <- acast(data = BigTrees,
#                            formula = TrapID~Species,
#                            fill = 0,
#                            fun.aggregate = length,
#                            value.var = "Species")
# 
# BigTreesMatrixMay <- acast(data = BigTreesMay,
#                            formula = TrapID~Species,
#                            fill = 0,
#                            fun.aggregate = length,
#                            value.var = "Species")
# 
# BigTreesMatrixMayMT <- acast(data = BigTreesMayMT,
#                            formula = TrapID~Species,
#                            fill = 0,
#                            fun.aggregate = length,
#                            value.var = "Species")
# 
# MandaiInsectsMatrix <- acast(data = MandaiInsects,
#                              formula = TrapID~Cluster,
#                              fill = 0,
#                              fun.aggregate = length,
#                              value.var = "Cluster")
# 
# MandaiInsectsMatrixMay <- acast(data = MandaiInsectsMay,
#                              formula = TrapID~Cluster,
#                              fill = 0,
#                              fun.aggregate = length,
#                              value.var = "Cluster")
# 
# MandaiInsectsMatrixMayMT <- acast(data = MandaiInsectsMayMT,
#                              formula = TrapID~Cluster,
#                              fill = 0,
#                              fun.aggregate = length,
#                              value.var = "Cluster")
# 
# ########
# #Rarefaction
# ###I've reused the SR and RarefiedSR names for each time I run the code since it keeps the colnames in each DF consistent 
# 
# SR <- specnumber(BigTreesMatrix)
# Raremax <- min(rowSums(BigTreesMatrix))
# RarefiedSR <- rarefy(BigTreesMatrix, Raremax)
# plot(SR,RarefiedSR,xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main = "Mandai Trees")
# print(SR)
# 
# 
# ###Adding the stats to the tree data
# SR <- as.data.frame(SR)
# RarefiedSR <- as.data.frame(RarefiedSR)
# BigTreesDF <- cbind(BigTreesDF,SR,RarefiedSR)
# 
# ##For tree data in May
# 
# SR <- specnumber(BigTreesMatrixMay)
# Raremax <- min(rowSums(BigTreesMatrixMay))
# RarefiedSR <- rarefy(BigTreesMatrixMay, Raremax)
# plot(SR,RarefiedSR,xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main = "Mandai Trees May")
# 
# SR <- as.data.frame(SR)
# RarefiedSR <- as.data.frame(RarefiedSR)
# BigTreesMayDF <- cbind(BigTreesMayDF,SR,RarefiedSR)
# 
# ###For tree data in May (MT data only)
# SR <- specnumber(BigTreesMatrixMayMT)
# Raremax <- min(rowSums(BigTreesMatrixMayMT))
# RarefiedSR <- rarefy(BigTreesMatrixMayMT, Raremax)
# plot(SR,RarefiedSR,xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main = "Mandai Trees May (MT)")
# 
# SR <- as.data.frame(SR)
# RarefiedSR <- as.data.frame(RarefiedSR)
# BigTreesMayMTDF <- cbind(BigTreesMayMTDF,SR,RarefiedSR)
# 
# ###For insect data 
# InsectSR <- specnumber(MandaiInsectsMatrix)
# Raremax <- min(rowSums(MandaiInsectsMatrix))
# RarefiedInsectSR <- rarefy(MandaiInsectsMatrix, Raremax)
# plot(InsectSR,RarefiedInsectSR,xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main = "Mandai Insects")
# 
# InsectSR <- as.data.frame(InsectSR)
# RarefiedInsectSR <- as.data.frame(RarefiedInsectSR)
# BigTreesDF <- cbind(BigTreesDF,InsectSR,RarefiedInsectSR)
# 
# ###For insect data (May)
# 
# InsectSR <- specnumber(MandaiInsectsMatrixMay)
# Raremax <- min(rowSums(MandaiInsectsMatrixMay))
# RarefiedInsectSR <- rarefy(MandaiInsectsMatrixMay, Raremax)
# plot(InsectSR,RarefiedInsectSR,xlab = "Observed No. of Species", ylab = "Rarefied No. of Species", main = "Mandai Insects May")
# 
# InsectSR <- as.data.frame(InsectSR)
# RarefiedInsectSR <- as.data.frame(RarefiedInsectSR)
# BigTreesMayDF <- cbind(BigTreesMayDF,InsectSR,RarefiedInsectSR)
# 
# ###For insect data (May and MT only)
# 
# InsectSR <- specnumber(MandaiInsectsMatrixMayMT)
# Raremax <- min(rowSums(MandaiInsectsMatrixMayMT))
# RarefiedInsectSR <- rarefy(MandaiInsectsMatrixMayMT, Raremax)
# plot(InsectSR,RarefiedInsectSR,xlab = "Observed No. of Species", ylab = "Rarefied No. of Species",main = "Mandai Insects May (MT)")
# 
# ###Adding the InsectSRats to the tree data
# InsectSR <- as.data.frame(InsectSR)
# RarefiedInsectSR <- as.data.frame(RarefiedInsectSR)
# BigTreesMayMTDF <- cbind(BigTreesMayMTDF,InsectSR,RarefiedInsectSR)
# 
# 






