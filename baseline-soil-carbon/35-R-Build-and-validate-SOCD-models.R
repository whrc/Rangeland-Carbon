library(tidyverse)
library(randomForest)
library(clhs)
setwd("set-pathway")


###############################################################
# Define parameters, site, and import datasets
# In the paper, spiking samples were from Year 2, calibration samples were from KSSL, and validation samples were from Year 1. 
ROI = "specify-site-name"
spike_num = 36 # Choose number of spiking cores
KSSL_num = 500 # Choose number of calibration samples
model_scheme = 1 # 0: Site-based model; 1: Pixel-based model
if (model_scheme == 0) {scheme = "site"} else {scheme = "pixel"}

# Import files for model calibration and validation
KSSL <- read.table('./KSSL.f.SOCSd.csv', comment.char ="", quote = "\"", header = T, sep = ",")
Vali <- read.table(paste0('./',ROI,'CSV-with-measured-SOCS-filename.csv')
                   , comment.char ="", quote = "\"", header = T, sep = ",")
Vali$density <- Vali$Measured_SOC * Vali$BD_m / 100
if (model_scheme == 0){
  CaliAll <- read.table(paste0('./SOCD/', ROI,'/', ROI,'_calibration_site_0_', KSSL_num,'.csv')
                        , comment.char ="", quote = "\"", header = T, sep = ",")
} else {
  CaliAll <- read.table(paste0('./SOCD/', ROI,'/', ROI,'_calibration_0_', KSSL_num,'.csv')
                        , comment.char ="", quote = "\"", header = T, sep = ",")
}
Spike_all <- read.table(paste0('./',ROI,'CSV-with-local-spiking-samples-joined-with-covariates.csv') 
                        , comment.char ="", quote = "\"", header = T, sep = ",")
Spike_all$density <- Spike_all$SOC_m * Spike_all$BD_m / 100
total_spike = length(unique(Spike_all$Point)) 
total_spike # Total number of available cores for spiking


###############################################################
# Get subset of spike samples
ncomp <- length(unique(Spike_all$Cluster))
cluster <- Spike_all %>% group_by(Cluster) %>% summarise (perc = n()/nrow(Spike_all))
Spike <- data.frame()
for(i in 1:ncomp){
  Spike_i <- Spike_all[which((Spike_all$Cluster==i)&(Spike_all$Med_depth < 10)),]
  size_i = round(cluster$perc[i] * spike_num)
  if (size_i == 0){size_i = 1}
  if (size_i < dim(Spike_i)[1]) {
    sample_indices <- clhs(Spike_i[,XX:XX], size_i, # Specify column names or numbers corresponding to covariates
                           progress = TRUE, iter = 500)
    thecores <- unique(Spike_i[sample_indices,]$Point)
    thesamples <- Spike_all[which(Spike_all$Point %in% thecores),]} else {
      thesamples <- Spike_all[which(Spike_all$Cluster==i),]}
  Spike <- rbind(Spike,thesamples)
}
length(unique(Spike$Point)) # Number of local cores for spiking
dim(Spike)[1] # Number of local samples for spiking
length(unique(Vali$Num)) # Number of local cores for validation
dim(Vali)[1] # Number of local samples for validation


#######################################################################
# Merge datasets to get KSSL calibration data for each point in the pixel-based model
Cali <- merge(CaliAll,KSSL,by = "sampleids")
if (model_scheme == 1) { # Pixel-based model
  Cali.i <- within(Cali, Cali_ID <- match(predid, unique(predid)))
  ID <- Cali.i[,c("predid","Cali_ID")] 
  ID <- unique(ID)
  Vali.i <- merge(Vali,ID,by.x = "Sample_ID",by.y="predid")
  Cali.i <- Cali.i[,c("sampleids","lay_id","Med_depth","density","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                      "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE","Cali_ID")]
  Spike <- Spike[,c("Point","X","Depth","density","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                    "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE")]
  Spike$Cali_ID = 0 # Placeholder for spiking samples
  colnames(Spike) <- colnames(Cali.i)}


#######################################################################
# Run model with spiking: combine KSSL and local samples
if (model_scheme == 1) { # Pixel-based model
  imp.cov <- data.frame()
  preds.i <- data.frame()
  Vali.individual.i <- data.frame()
  rmse.f = 0
  r2.f = 0
  mbe.f = 0
  for (i in 1:length(unique(Cali.i$Cali_ID))){
    Cali.KSSL <- Cali.i[which(Cali.i$Cali_ID == i),]
    Cali.individual <- rbind(Cali.KSSL,Spike)
    Vali.individual <- Vali.i[which(Vali.i$Cali_ID == i),]
    rf.i <- randomForest(density ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                           Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                         data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    r2 <- summary(lm(rf.i$predicted ~ Cali.individual$density))$adj.r.squared
    r2.f <- r2 + r2.f
    resids <- rf.i$predicted - Cali.individual$density
    rmse <- sqrt(mean(resids^2))
    rmse.f <- rmse + rmse.f
    mbe <- mean(resids)
    mbe.f <- mbe + mbe.f
    imp <- importance(rf.i, type = 1)[,1]
    imp.cov <- rbind(imp.cov, imp)
    preds <- predict(rf.i, newdata = Vali.individual)
    preds.i <- rbind(preds.i,preds) 
    Vali.individual.i <- rbind(Vali.individual.i,Vali.individual[,c("Num","Med_depth","density")])
  }
  r2.c <- r2.f/length(unique(Vali.i$Cali_ID))
  rmse.c <- rmse.f/length(unique(Vali.i$Cali_ID))
  mbe.c <- mbe.f/length(unique(Vali.i$Cali_ID))
  colnames(imp.cov) <- names(imp)
  imp.cov <- sapply(imp.cov, mean)
  write.csv(imp.cov,paste0("./SOCD/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_",spike_num,"_imp.csv"))
  dotchart(sort(imp.cov), xlab="Variable importance (%)")
  colnames(preds.i) <- c("preds.i")
  colnames(Vali.individual.i) <- c("Core","Med_depth","vali.i")
  vali.result <- cbind(preds.i,Vali.individual.i)
  r2.v <- summary(lm(vali.result$preds.i~vali.result$vali.i))
  rmse.v <- sqrt(mean((vali.result$preds.i-vali.result$vali.i)^2))
  mbe.v <- mean(vali.result$preds.i-vali.result$vali.i)
  plot(vali.result$vali.i, vali.result$preds.i, 
       main = "Measured vs. predicted SOC density of the validation set \n
       (KSSL spiked with local samples for individual-based model)", cex.main = 0.95,
       xlab = "Measured SOCD (g/cm3)",
       ylab = "Predicted SOCD (g/cm3)",
       xlim = c(0,0.1), # Specify range for SOCD visualization
       ylim = c(0,0.1))
  abline(coef = c(0,1))
  vali.result <- vali.result[,c("Core","Med_depth","vali.i","preds.i")]
  write.csv(vali.result, paste0("./SOCD/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_",spike_num,"_pred.csv"))
}
if (model_scheme == 0) { # Site-based model
  Spike <- Spike[,c("Sample_ID","Med_depth","density","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                    "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE")]
  Cali <- Cali[,c("sampleids","Med_depth","density","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                  "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE")]
  colnames(Spike) <- colnames(Cali)
  Site_Cali <- rbind(Cali,Spike)
  rf <- randomForest(density ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                       Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                     data = Site_Cali, importance = TRUE, ntree=500, mtry=3)
  r2.c <- summary(lm(rf$predicted ~ Site_Cali$density))$adj.r.squared
  r2.c
  resids.c <- rf$predicted - Site_Cali$density
  rmse.c <- sqrt(mean(resids.c^2))
  rmse.c
  mbe.c <- mean(rf$predicted - Site_Cali$density)
  mbe.c
  imp <- importance(rf, type = 1)[,1]
  write.csv(imp,paste0("./SOCD/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_",spike_num,"_imp.csv"))
  dotchart(sort(imp), xlab="Variable importance (%)")
  preds <- predict(rf, newdata = Vali)
  r2.v <- summary(lm(preds ~ Vali$density))$adj.r.squared
  resids.v <- preds - Vali$density
  rmse.v <- sqrt(mean(resids.v^2))
  mbe.v <- mean(preds - Vali$density)
  plot(Vali$density, preds, 
       main = "Measured vs. predicted SOC density of the validation set \n
     (KSSL spiked with local samples for site-based model)", cex.main = 0.95,
       xlab = "Measured SOCD (g/cm3)",
       ylab = "Predicted SOCD (g/cm3)",
       xlim = c(0,0.1), # Specify range for SOCD visualization
       ylim = c(0,0.1))
  abline(coef = c(0,1))
  vali.result <- cbind(Vali[,c("Num","Med_depth","SOCD")],preds)
  colnames(vali.result) <- c("Core","Med_depth","vali.i","preds.i")
  write.csv(vali.result, paste0("./SOCD/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_",spike_num,"_pred.csv"))
}
r2.c
rmse.c
mbe.c 
r2.v
rmse.v
mbe.v


#################################################################################
# Run model without spiking: KSSL only 
if (model_scheme == 1) { # Pixel-based model
  imp.cov <- data.frame()
  preds.i <- data.frame()
  Vali.individual.i <- data.frame()
  rmse.f = 0
  r2.f = 0
  mbe.f = 0
  for (i in 1:length(unique(Cali.i$Cali_ID))){
    Cali.individual <- Cali.i[which(Cali.i$Cali_ID == i),]
    Vali.individual <- Vali.i[which(Vali.i$Cali_ID == i),]
    rf.i <- randomForest(density ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                           Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                         data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    r2 <- summary(lm(rf.i$predicted ~ Cali.individual$density))$adj.r.squared
    r2.f <- r2 + r2.f
    resids <- rf.i$predicted - Cali.individual$density
    rmse <- sqrt(mean(resids^2))
    rmse.f <- rmse + rmse.f
    mbe <- mean(resids)
    mbe.f <- mbe + mbe.f
    imp <- importance(rf.i, type = 1)[,1]
    imp.cov <- rbind(imp.cov, imp)
    preds <- predict(rf.i, newdata = Vali.individual)
    preds.i <- rbind(preds.i,preds) 
    Vali.individual.i <- rbind(Vali.individual.i,Vali.individual[,c("Num","Med_depth","density")])
  }
  r2.c <- r2.f/length(unique(Vali.i$Cali_ID))
  rmse.c <- rmse.f/length(unique(Vali.i$Cali_ID))
  mbe.c <- mbe.f/length(unique(Vali.i$Cali_ID))
  colnames(imp.cov) <- names(imp)
  imp.cov <- sapply(imp.cov, mean)
  write.csv(imp.cov,paste0("./SOCD/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_0_imp.csv"))
  dotchart(sort(imp.cov), xlab="Variable importance (%)")
  colnames(preds.i) <- c("preds.i")
  colnames(Vali.individual.i) <- c("Core","Med_depth","vali.i")
  vali.result <- cbind(preds.i,Vali.individual.i)
  r2.v <- summary(lm(vali.result$preds.i~vali.result$vali.i))
  rmse.v <- sqrt(mean((vali.result$preds.i-vali.result$vali.i)^2))
  mbe.v <- mean(vali.result$preds.i-vali.result$vali.i)
  plot(vali.result$vali.i, vali.result$preds.i, 
       main = "Measured vs. predicted SOC density of the validation set \n
       (KSSL only individual-based model)", cex.main = 0.95,
       xlab = "Measured SOCD (g/cm3)",
       ylab = "Predicted SOCD (g/cm3)",
       xlim = c(0,0.1), # Specify range for SOCD visualization
       ylim = c(0,0.1))
  abline(coef = c(0,1))
  vali.result <- vali.result[,c("Core","Med_depth","vali.i","preds.i")]
  write.csv(vali.result, paste0("./SOCD/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_0_pred.csv"))
}
if (model_scheme == 0) { # Site-based model
  rf <- randomForest(density ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                       Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                     data = Cali, importance = TRUE, ntree=500, mtry=3)
  r2.c <- summary(lm(rf$predicted ~ Cali$density))$adj.r.squared
  r2.c
  resids.c <- rf$predicted - Cali$density
  rmse.c <- sqrt(mean(resids.c^2))
  rmse.c
  mbe.c <- mean(rf$predicted - Cali$density)
  mbe.c
  imp <- importance(rf, type = 1)[,1]
  write.csv(imp,paste0("./SOCD/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_0_imp.csv"))
  dotchart(sort(imp), xlab="Variable importance (%)")
  preds <- predict(rf, newdata = Vali)
  r2.v <- summary(lm(preds ~ Vali$density))$adj.r.squared
  resids.v <- preds - Vali$density
  rmse.v <- sqrt(mean(resids.v^2))
  mbe.v <- mean(preds - Vali$density)
  plot(Vali$density, preds, 
       main = "Measured vs. predicted SOC density of the validation set \n
     (KSSL only site-based model)", cex.main = 0.95,
       xlab = "Measured SOCD (g/cm3)",
       ylab = "Predicted SOCD (g/cm3)",
       xlim = c(0,0.1), # Specify range for SOCD visualization
       ylim = c(0,0.1))
  abline(coef = c(0,1))
  vali.result <- cbind(Vali[,c("Num","Med_depth","density")],preds)
  colnames(vali.result) <- c("Core","Med_depth","vali.i","preds.i")
  write.csv(vali.result, paste0("./SOCD/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_0_pred.csv"))
}
r2.c
rmse.c
mbe.c 
r2.v
rmse.v
mbe.v 


#################################################################################
# Run model with local samples only
rf <- randomForest(density ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                     Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                   data = Spike, importance = TRUE, ntree=500, mtry=3)
r2.c <- summary(lm(rf$predicted ~ Spike$density))$adj.r.squared
r2.c
resids.c <- rf$predicted - Spike$density
rmse.c <- sqrt(mean(resids.c^2))
rmse.c
mbe.c <- mean(rf$predicted - Spike$density)
mbe.c
imp <- importance(rf, type = 1)[,1]
write.csv(imp, paste0("2022_calibration/Results/SOCD/",ROI,"/",ROI,"_",scheme,"_0_",spike_num,"_imp.csv"))
dotchart(sort(imp), xlab="Variable importance (%)")
preds <- predict(rf, newdata = Vali)
r2.v <- summary(lm(preds ~ Vali$density))$adj.r.squared
r2.v
resids.v <- preds - Vali$density
rmse.v <- sqrt(mean(resids.v^2))
rmse.v
mbe.v <- mean(preds - Vali$density)
mbe.v
plot(Vali$density, preds, 
     main = "Measured vs. predicted SOC density of the validation set \n
     (local only model)", cex.main = 0.95,
     xlab = "Measured SOCD (g/cm3)",
     ylab = "Predicted SOCD (g/cm3)",
     xlim = c(0,0.1), # Specify range for SOCD visualization
     ylim = c(0,0.1))
abline(coef = c(0,1))
vali.result <- cbind(Vali[,c("Num","Med_depth","density")],preds)
colnames(vali.result) <- c("Core","Med_depth","vali.i","preds.i")
write.csv(vali.result, paste0("./SOCD/",ROI,"/",ROI,"_0_",spike_num,"_pred.csv"))
