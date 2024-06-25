library(caret)
library(rgdal)
library(flexclust)
library(tidyverse)
library(Rfast)
setwd("set-pathway")


###############################################################
# Define parameters, site, and import datasets
# Select site and desired parameters
ROI = "specify-site-name"
KSSL_num = 500 
# Import KSSL dataset for calibration
Cali<- read.table('./KSSL.f.SOCS.csv', comment.char ="", quote = "\"", header = T, sep = ",")
Cali$Dataset <- "Cali"
Cali$x <- 0
Cali$y <- 0
Cali <- na.omit(Cali) # Only use western US samples that have RAP results
Cali_reorg <- Cali[,c("Dataset","sampleids","x","y",
                      "BD","SOC","Clay","Sand",
                      "ppt","tmin","tmax","VPD",
                      "Tree","GPP","EVI","AS","EL","SL","TWI",
                      "AFGC","BG","LTR","PFGC","SHR","TREE")]
# Import raster image of target site with extracted covariates
predict.raster <- paste0('./',ROI,'_cov_all_RAP.tif')
x <- new("GDALReadOnlyDataset", predict.raster)
width <- dim(x)[2]
height <- dim(x)[1]
imagedata <- data.frame(getRasterTable(x))
names(imagedata) <- c('x','y','EL','SL','AS','TWI','ppt','tmin','tmax','VPD',
                      'BD1','BD2','BD3','BD4','BD5','Clay1','Clay2','Clay3','Clay4','Clay5',
                      'Sand1','Sand2','Sand3','Sand4','Sand5','SOC2','SOC5','SOC2','SOC5','SOC5',
                      'Tree','GPP','EVI','AFGC','PFGC','TREE','SHR','LTR','BG',
                      'Med_depth')
# Select the depth corresponding to med_depth  = 15
colnames(imagedata)[colnames(imagedata) %in% c("BD3","Clay3","Sand3","SOC3")] <- c("BD","Clay","Sand","SOC")
# Scale covariates to match with KSSL data
imagedata$SOC <- imagedata$SOC/10
imagedata$ppt <- imagedata$ppt/20
imagedata$tmin <- imagedata$tmin/20
imagedata$tmax <- imagedata$tmax/20
imagedata$VPD <- imagedata$VPD/20
imagedata$GPP <- imagedata$GPP/20
imagedata$AS <- sin(pi*(imagedata$AS)/360)
imagedata <- na.omit(imagedata)
# Reorganize the grid dataset
imagedata$Dataset <- "Pred"
pred.ID <- tibble::rowid_to_column(imagedata, "ID")
imagedata <- cbind(pred.ID[,1], imagedata)
colnames(imagedata)[1] <- c("pred.ID")
Pred_reorg <- imagedata[,c("Dataset","pred.ID","x","y",
                           "BD","SOC","Clay","Sand",
                           "ppt","tmin","tmax","VPD",
                           "Tree","GPP","EVI","AS","EL","SL","TWI",
                           "AFGC","BG","LTR","PFGC","SHR","TREE")]
colnames(Pred_reorg) <- colnames(Cali_reorg)
combo <- rbind(Cali_reorg, Pred_reorg) 


###############################################################
# Carry out PCA and distance calculation to select similar samples in KSSL for subsequent calibration
combo.p <- combo[,c(5:25)]  # Or specify different column number or names to match with covariate datasets
combo.p.pca = prcomp(combo.p, scale. = TRUE)
var_explained <- combo.p.pca$sdev^2/sum(combo.p.pca$sdev^2)
scores.p = combo.p.pca$x
s.p <- scores.p [,1:(max(which(var_explained > 0.05)))] # Retain important PCs 
combo.new <- cbind(combo,s.p)
combo.new.cali.p <- combo.new[which(combo.new$Dataset == "Cali"),c("Dataset","sampleids",colnames(s.p))]
combo.new.pred.p <- combo.new[which(combo.new$Dataset == "Pred"),c("Dataset","sampleids",colnames(s.p))]
combo.new.cali.p.ID <- tibble::rowid_to_column(combo.new.cali.p, "new_ID")
combo.new.pred.p.ID <- tibble::rowid_to_column(combo.new.pred.p, "new_ID")
cali_ID <- combo.new.cali.p.ID[,c("new_ID","sampleids")]
pred_ID <- combo.new.pred.p.ID[,c("new_ID","sampleids")]
r.p <- dist2(combo.new.pred.p.ID[,4:(dim(s.p)[2]+3)],combo.new.cali.p.ID[,4:(dim(s.p)[2]+3)],method = "euclidean")
rm("combo.p.pca","s.p","scores.p","combo.p","combo.new")
threshold.i <- Rfast::rownth(as.matrix(r.p), rep(KSSL_num, nrow(r.p)))
r.i.cutoff <- cbind(r.p,threshold.i)
rm("r.p","combo.new.cali.p")
r.i.cutoff.zero <- sweep(r.i.cutoff, 1, threshold.i)
rm("r.i.cutoff","threshold.i")
cali.i <- which(r.i.cutoff.zero[,1:(dim(r.i.cutoff.zero)[2]-1)] <= 0,arr.ind=T)
# Generate pixel-based KSSL calibration data ID and export
cali.reorg.i <- data.frame(cali.i)
colnames(cali.reorg.i)[1:2] <- c("ID","ID.cali")
calibration.i <- merge(cali.reorg.i,cali_ID, by.x = "ID.cali", by.y = "new_ID")
calibration.i <- merge(calibration.i,pred_ID,  by.x = "ID", by.y = "new_ID")
calibration.i <- cbind(calibration.i[,4],calibration.i[,3])
colnames(calibration.i)[1:2] <- c("predid","sampleids")
write.csv(calibration.i,paste0("./",ROI,"_calibration_",KSSL_num,"_pred_set_SOCS.csv"))
# Generate site-based KSSL calibration data ID and export
cali.reorg.i <- data.frame(cali.i)
colnames(cali.reorg.i)[1:2] <- c("ID","ID.cali")
calibration.i <- merge(cali.reorg.i,cali_ID, by.x = "ID.cali", by.y = "new_ID")
calibration.i <- merge(calibration.i,pred_ID,  by.x = "ID", by.y = "new_ID")
calibration.i <- data.frame(calibration.i[,3])
colnames(calibration.i)[1] <- c("sampleids")
calibration.site <- data.frame(unique(calibration.i$sampleids))
colnames(calibration.site)[1] <- c("sampleids")
write.csv(calibration.i,paste0("./",ROI,"_calibration_site_",KSSL_num,"_pred_set_SOCS.csv"))
