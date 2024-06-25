library(tidyverse)
library(Rfast)
library(flexclust)
setwd("set-pathway")


###############################################################
# Set parameters, modeling approach, and import datasets
ROI = "specify-site-name"
KSSL_num = 500


# Extract SOC calibration set
###############################################################
# Import calibration (KSSL) and validation (Year 1 samples) datasets
Cali<- read.table('./KSSL.f.2022.csv', comment.char ="", quote = "\"", header = T, sep = ",")
Vali<- read.table('./CSV-with-SOC-and-site-covariates-filename.csv',comment.char ="", quote = "\"", header = T, sep = ",")
Cali$Dataset <- "Cali"
Vali$Dataset <- "Vali"
# Refine KSSL dataset
Cali <- Cali[which((Cali$Med_depth <120)),]
# Use the same unit for calibration and validation datasets, if not yet done
Vali$ppt <- Vali$ppt/20
Vali$tmin <- Vali$tmin/20
Vali$tmax <- Vali$tmax/20
Vali$VPD <- Vali$VPD/20
Vali$GPP <- Vali$GPP/20
Vali$AS <- sin(pi*(Vali$AS)/360)
Vali$SOC <- Vali$SOC/10
# First reorganize the datasets in order to combine them for PCA analysis
Vali_reorg <- Vali[,c("Dataset","X.1","Sample_ID","Num","Med_depth","Measured_SOC","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                      "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE")]
Cali_reorg <- cbind(Cali[,c("Dataset")],Cali[,2:dim(Cali)[2]-1])
Cali_reorg <- na.omit(Cali_reorg) # Remove data points outside of western U.S. in order to use RAP covariates, if not yet done in previous steps
colnames(Cali_reorg)[1] <- "Dataset"
Cali_reorg <- Cali_reorg[,c("Dataset","X","sampleids","lay_id","Med_depth","oc_usda",
                            "BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD","LULC",
                            "Tree","GPP","EVI","AS","EL","SL","TWI",
                            "AFGC","BG","LTR","PFGC","SHR","TREE")]
colnames(Vali_reorg) <- colnames(Cali_reorg)
# Select KSSL calibration dataset based on similarity of covariate space
combo <- rbind(Cali_reorg, Vali_reorg) # Combine calibration and prediction/validation sets
combo.p <- combo[,c(7:14,16:28)] # Exclude LULC (categorical) covariate in the model
combo.p.pca = prcomp(combo.p, scale. = TRUE) # Use scaled values to run PCA
var_explained <- combo.p.pca$sdev^2/sum(combo.p.pca$sdev^2)
scores.p = combo.p.pca$x
s.p <- scores.p [,1:(max(which(var_explained > 0.05)))] # Retain PCs explained >5% 
combo.new <- cbind(combo,s.p)
combo.new.cali.p <- combo.new[which(combo.new$Dataset == "Cali"),c("Dataset","sampleids",colnames(s.p))]
combo.new.vali.p <- combo.new[which(combo.new$Dataset == "Vali"),c("Dataset","sampleids",colnames(s.p))]
combo.new.cali.p.ID <- tibble::rowid_to_column(combo.new.cali.p, "new_ID")
combo.new.vali.p.ID <- tibble::rowid_to_column(combo.new.vali.p, "new_ID")
cali_ID <- combo.new.cali.p.ID[,c("new_ID","sampleids")]
vali_ID <- combo.new.vali.p.ID[,c("new_ID","sampleids")]
# Use distance calculation to identify similar KSSL samples
r.p <- dist2(combo.new.vali.p.ID[,4:(dim(s.p)[2]+3)],combo.new.cali.p.ID[,4:(dim(s.p)[2]+3)],method = "euclidean")
threshold.i <- Rfast::rownth(as.matrix(r.p), rep(KSSL_num, nrow(r.p)))
r.i.cutoff <- cbind(r.p,threshold.i)
r.i.cutoff.zero <- sweep(r.i.cutoff, 1, threshold.i)
cali.i <- which(r.i.cutoff.zero[,1:(dim(r.i.cutoff.zero)[2]-1)] <= 0,arr.ind=T)
cali.reorg.i <- data.frame(cali.i)
colnames(cali.reorg.i)[1:2] <- c("ID","ID.cali")
calibration.i <- merge(cali.reorg.i,cali_ID, by.x = "ID.cali", by.y = "new_ID")
calibration.i <- merge(calibration.i,vali_ID,  by.x = "ID", by.y = "new_ID")
calibration.i <- data.frame(calibration.i[,3])
colnames(calibration.i)[1] <- c("sampleids")
calibration.site <- data.frame(unique(calibration.i$sampleids))
colnames(calibration.site)[1] <- c("sampleids")
dim(calibration.site) # Find out the number of unique KSSL samples retained
write.csv(calibration.site,paste0("./SOC/",ROI,"/",ROI,"_calibration_site_0_",KSSL_num,".csv"))


# Extract BD calibration set
###############################################################
# Import calibration (KSSL) and validation (Year 1 samples) datasets
Cali<- read.table('./KSSL.f.w.BD.csv', comment.char ="", quote = "\"", header = T, sep = ",")
Vali<- read.table('./CSV-with-BD-and-site-covariates-filename.csv',comment.char ="", quote = "\"", header = T, sep = ",")
Cali$Dataset <- "Cali"
Vali$Dataset <- "Vali"
# Refine KSSL dataset
Cali <- Cali[which((Cali$Med_depth <120)),]
# Use the same unit for calibration and validation datasets, if not yet done
Vali$ppt <- Vali$ppt/20
Vali$tmin <- Vali$tmin/20
Vali$tmax <- Vali$tmax/20
Vali$VPD <- Vali$VPD/20
Vali$GPP <- Vali$GPP/20
Vali$AS <- sin(pi*(Vali$AS)/360)
Vali$SOC <- Vali$SOC/10
# First reorganize the datasets in order to combine them for PCA analysis
Vali_reorg <- Vali[,c("Dataset","X.1","Sample_ID","Num","Med_depth","Measured_SOC","Measured_BD","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                      "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE")]
Cali_reorg <- cbind(Cali[,c("Dataset")],Cali[,2:dim(Cali)[2]-1])
Cali_reorg <- na.omit(Cali_reorg) # Remove data points outside of western U.S. in order to use RAP covariates, if not yet done in previous steps
colnames(Cali_reorg)[1] <- "Dataset"
Cali_reorg <- Cali_reorg[,c("Dataset","X","sampleids","lay_id","Med_depth","oc_usda","bd_kssl",
                            "BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD","LULC",
                            "Tree","GPP","EVI","AS","EL","SL","TWI",
                            "AFGC","BG","LTR","PFGC","SHR","TREE")]
colnames(Vali_reorg) <- colnames(Cali_reorg)
# Select KSSL calibration dataset based on similarity of covariate space
combo <- rbind(Cali_reorg, Vali_reorg) # Combine calibration and prediction/validation sets
combo.p <- combo[,c(8:15,17:29)] # Exclude LULC (categorical) covariate in the model
combo.p.pca = prcomp(combo.p, scale. = TRUE) # Use scaled values to run PCA
var_explained <- combo.p.pca$sdev^2/sum(combo.p.pca$sdev^2)
scores.p = combo.p.pca$x
s.p <- scores.p [,1:(max(which(var_explained > 0.05)))] # Retain PCs explained >5% 
combo.new <- cbind(combo,s.p)
combo.new.cali.p <- combo.new[which(combo.new$Dataset == "Cali"),c("Dataset","sampleids",colnames(s.p))]
combo.new.vali.p <- combo.new[which(combo.new$Dataset == "Vali"),c("Dataset","sampleids",colnames(s.p))]
combo.new.cali.p.ID <- tibble::rowid_to_column(combo.new.cali.p, "new_ID")
combo.new.vali.p.ID <- tibble::rowid_to_column(combo.new.vali.p, "new_ID")
cali_ID <- combo.new.cali.p.ID[,c("new_ID","sampleids")]
vali_ID <- combo.new.vali.p.ID[,c("new_ID","sampleids")]
# Use distance calculation to identify similar KSSL samples
r.p <- dist2(combo.new.vali.p.ID[,4:(dim(s.p)[2]+3)],combo.new.cali.p.ID[,4:(dim(s.p)[2]+3)],method = "euclidean")
threshold.i <- Rfast::rownth(as.matrix(r.p), rep(KSSL_num, nrow(r.p)))
r.i.cutoff <- cbind(r.p,threshold.i)
r.i.cutoff.zero <- sweep(r.i.cutoff, 1, threshold.i)
cali.i <- which(r.i.cutoff.zero[,1:(dim(r.i.cutoff.zero)[2]-1)] <= 0,arr.ind=T)
cali.reorg.i <- data.frame(cali.i)
colnames(cali.reorg.i)[1:2] <- c("ID","ID.cali")
calibration.i <- merge(cali.reorg.i,cali_ID, by.x = "ID.cali", by.y = "new_ID")
calibration.i <- merge(calibration.i,vali_ID,  by.x = "ID", by.y = "new_ID")
calibration.i <- data.frame(calibration.i[,3])
colnames(calibration.i)[1] <- c("sampleids")
calibration.site <- data.frame(unique(calibration.i$sampleids))
colnames(calibration.site)[1] <- c("sampleids")
dim(calibration.site) # Find out the number of unique KSSL samples retained
write.csv(calibration.site,paste0("./BD/",ROI,"/",ROI,"_calibration_site_0_",KSSL_num,".csv"))


# Extract SOCS calibration set
###############################################################
# Import calibration (KSSL) and validation (Year 1 samples) datasets
Cali<- read.table('./KSSL.f.SOCS.csv', comment.char ="", quote = "\"", header = T, sep = ",")
Vali<- read.table('./CSV-with-SOCS-and-site-covariates-filename.csv',comment.char ="", quote = "\"", header = T, sep = ",")
Cali$Dataset <- "Cali"
Vali$Dataset <- "Vali"
# Refine KSSL dataset
Cali <- Cali[which((Cali$Med_depth <120)),]
# Use the same unit for calibration and validation datasets, if not yet done
Vali$ppt <- Vali$ppt/20
Vali$tmin <- Vali$tmin/20
Vali$tmax <- Vali$tmax/20
Vali$VPD <- Vali$VPD/20
Vali$GPP <- Vali$GPP/20
Vali$AS <- sin(pi*(Vali$AS)/360)
Vali$SOC <- Vali$SOC/10
# First reorganize the datasets in order to combine them for PCA analysis
Vali_reorg <- Vali[,c("Dataset","X","Sample_ID","Num","Depth","SOCS_sum","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                      "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE")]
Cali_reorg <- cbind(Cali[,c("Dataset")],Cali[,2:dim(Cali)[2]-1])
Cali_reorg <- na.omit(Cali_reorg) # Remove data points outside of western U.S. in order to use RAP covariates, if not yet done in previous steps
colnames(Cali_reorg)[1] <- "Dataset"
Cali_reorg <- Cali_reorg[,c("Dataset","X.1","sampleids","lay_id","Med_depth","socs_kssl",
                            "BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD","LULC",
                            "Tree","GPP","EVI","AS","EL","SL","TWI",
                            "AFGC","BG","LTR","PFGC","SHR","TREE")]
colnames(Vali_reorg) <- colnames(Cali_reorg)
# Select KSSL calibration dataset based on similarity of covariate space
combo <- rbind(Cali_reorg, Vali_reorg) # Combine calibration and prediction/validation sets
combo.p <- combo[,c(7:14,16:28)]# Exclude LULC (categorical) covariate in the model
combo.p.pca = prcomp(combo.p, scale. = TRUE) # Use scaled values to run PCA
var_explained <- combo.p.pca$sdev^2/sum(combo.p.pca$sdev^2)
scores.p = combo.p.pca$x
s.p <- scores.p [,1:(max(which(var_explained > 0.05)))] # Retain PCs explained >5% 
combo.new <- cbind(combo,s.p)
combo.new.cali.p <- combo.new[which(combo.new$Dataset == "Cali"),c("Dataset","sampleids",colnames(s.p))]
combo.new.vali.p <- combo.new[which(combo.new$Dataset == "Vali"),c("Dataset","sampleids",colnames(s.p))]
combo.new.cali.p.ID <- tibble::rowid_to_column(combo.new.cali.p, "new_ID")
combo.new.vali.p.ID <- tibble::rowid_to_column(combo.new.vali.p, "new_ID")
cali_ID <- combo.new.cali.p.ID[,c("new_ID","sampleids")]
vali_ID <- combo.new.vali.p.ID[,c("new_ID","sampleids")]
# Use distance calculation to identify similar KSSL samples
r.p <- dist2(combo.new.vali.p.ID[,4:(dim(s.p)[2]+3)],combo.new.cali.p.ID[,4:(dim(s.p)[2]+3)],method = "euclidean")
threshold.i <- Rfast::rownth(as.matrix(r.p), rep(KSSL_num, nrow(r.p)))
r.i.cutoff <- cbind(r.p,threshold.i)
r.i.cutoff.zero <- sweep(r.i.cutoff, 1, threshold.i)
cali.i <- which(r.i.cutoff.zero[,1:(dim(r.i.cutoff.zero)[2]-1)] <= 0,arr.ind=T)
cali.reorg.i <- data.frame(cali.i)
colnames(cali.reorg.i)[1:2] <- c("ID","ID.cali")
calibration.i <- merge(cali.reorg.i,cali_ID, by.x = "ID.cali", by.y = "new_ID")
calibration.i <- merge(calibration.i,vali_ID,  by.x = "ID", by.y = "new_ID")
calibration.i <- data.frame(calibration.i[,3])
colnames(calibration.i)[1] <- c("sampleids")
calibration.site <- data.frame(unique(calibration.i$sampleids))
colnames(calibration.site)[1] <- c("sampleids")
dim(calibration.site) # Find out the number of unique KSSL samples retained
write.csv(calibration.site,paste0("./SOCS/",ROI,"/",ROI,"_calibration_site_0_",KSSL_num,".csv"))