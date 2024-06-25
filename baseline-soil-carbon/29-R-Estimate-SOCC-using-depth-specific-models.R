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
# Import files for model calibration and validation
KSSL <- read.table('./KSSL.f.2022.csv', comment.char ="", quote = "\"", header = T, sep = ",")
Vali <- read.table(paste0('./',ROI,'_vali_0.csv') # Choose all samples from Year 1 for validation
                   , comment.char ="", quote = "\"", header = T, sep = ",")
Vali1 <- Vali[which(Vali$Med_depth <15),]
Vali2 <- Vali[which((Vali$Med_depth >=15)&(Vali$Med_depth <30)),]
Vali3 <- Vali[which(Vali$Med_depth >30),]
CaliAll1 <- read.table(paste0('./SOC/D1/', ROI,'/', ROI,'_calibration_0_', KSSL_num,'.csv')
                       , comment.char ="", quote = "\"", header = T, sep = ",")
CaliAll2 <- read.table(paste0('./SOC/D2/', ROI,'/', ROI,'_calibration_0_', KSSL_num,'.csv')
                       , comment.char ="", quote = "\"", header = T, sep = ",")
CaliAll3 <- read.table(paste0('./SOC/D3/', ROI,'/', ROI,'_calibration_0_', KSSL_num,'.csv')
                       , comment.char ="", quote = "\"", header = T, sep = ",")
Spike_all <- read.table(paste0('./',ROI,'CSV-with-local-spiking-samples-joined-with-covariates.csv'), 
                        comment.char ="", quote = "\"", header = T, sep = ",")
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
# Merge datasets to get KSSL calibration data for each point in the pixel-based model by depth
Cali1 <- merge(CaliAll1,KSSL,by = "sampleids")
Cali2 <- merge(CaliAll2,KSSL,by = "sampleids")
Cali3 <- merge(CaliAll3,KSSL,by = "sampleids")
Cali.i1 <- within(Cali1, Cali_ID1 <- match(predid, unique(predid)))
ID1 <- Cali.i1[,c("predid","Cali_ID1")] 
ID1 <- unique(ID1) 
Cali.i2 <- within(Cali2, Cali_ID2 <- match(predid, unique(predid)))
ID2 <- Cali.i2[,c("predid","Cali_ID2")] 
ID2 <- unique(ID2) 
Cali.i3 <- within(Cali3, Cali_ID3 <- match(predid, unique(predid)))
ID3 <- Cali.i3[,c("predid","Cali_ID3")] 
ID3 <- unique(ID3) 
Vali.i1 <- merge(Vali1,ID1,by.x = "Sample_ID",by.y="predid")
Vali.i2 <- merge(Vali2,ID2,by.x = "Sample_ID",by.y="predid")
Vali.i3 <- merge(Vali3,ID3,by.x = "Sample_ID",by.y="predid")
Cali.i1 <- Cali.i1[,c("sampleids","lay_id","Med_depth","oc_usda","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                      "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE","Cali_ID1")]
Cali.i2 <- Cali.i1[,c("sampleids","lay_id","Med_depth","oc_usda","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                      "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE","Cali_ID2")]
Cali.i3 <- Cali.i1[,c("sampleids","lay_id","Med_depth","oc_usda","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                      "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE","Cali_ID3")]
Spike <- Spike[,c("Sample_ID","X","Med_depth","SOC_m","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                    "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE")]
Spike1 <- Spike[which(Spike$Med_depth <15),]
Spike2 <- Spike[which((Spike$Med_depth >=15)&(Spike$Med_depth <30)),]
Spike3 <- Spike[which(Spike$Med_depth >30),]
Spike1$Cali_ID1 = 0 
Spike2$Cali_ID2 = 0 
Spike3$Cali_ID3 = 0 
colnames(Spike1) <- colnames(Cali.i1)
colnames(Spike2) <- colnames(Cali.i2)
colnames(Spike3) <- colnames(Cali.i3)


# Run model with spiking by depth: combine KSSL and local samples
#######################################################################
# Top depth
imp.cov <- data.frame()
preds.i <- data.frame()
Vali.individual.i <- data.frame()
rmse.f = 0
r2.f = 0
mbe.f = 0
for (i in 1:length(unique(Cali.i1$Cali_ID1))){
    Cali.KSSL <- Cali.i1[which(Cali.i1$Cali_ID1 == i),]
    Cali.individual <- rbind(Cali.KSSL,Spike1)
    Vali.individual <- Vali.i1[which(Vali.i1$Cali_ID1 == i),]
    rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                           Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                         data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    r2 <- summary(lm(rf.i$predicted ~ Cali.individual$oc_usda))$adj.r.squared
    r2.f <- r2 + r2.f
    resids <- rf.i$predicted - Cali.individual$oc_usda
    rmse <- sqrt(mean(resids^2))
    rmse.f <- rmse + rmse.f
    mbe <- mean(resids)
    mbe.f <- mbe + mbe.f
    imp <- importance(rf.i, type = 1)[,1]
    imp.cov <- rbind(imp.cov, imp)
    preds <- predict(rf.i, newdata = Vali.individual)
    preds.i <- rbind(preds.i,preds) 
    Vali.individual.i <- rbind(Vali.individual.i,Vali.individual[,c("Num","Med_depth","Measured_SOC")])
}
r2.c <- r2.f/length(unique(Vali.i1$Cali_ID1))
rmse.c <- rmse.f/length(unique(Vali.i1$Cali_ID1))
mbe.c <- mbe.f/length(unique(Vali.i1$Cali_ID1))
colnames(imp.cov) <- names(imp)
imp.cov <- sapply(imp.cov, mean)
write.csv(imp.cov,paste0("./SOC/Imp3/D1/",ROI,"/",ROI,"_",KSSL_num,"_",spike_num,"_imp.csv"))
dotchart(sort(imp.cov), xlab="Variable importance (%)")
colnames(preds.i) <- c("preds.i")
colnames(Vali.individual.i) <- c("Core","Med_depth","vali.i")
vali.result <- cbind(preds.i,Vali.individual.i)
r2.v <- summary(lm(vali.result$preds.i~vali.result$vali.i))
rmse.v <- sqrt(mean((vali.result$preds.i-vali.result$vali.i)^2))
mbe.v <- mean(vali.result$preds.i-vali.result$vali.i)
plot(vali.result$vali.i, vali.result$preds.i, 
       main = "Measured vs. predicted SOC of the validation set \n
       (KSSL spiked with local samples for individual-based model)", cex.main = 0.95,
       xlab = "Measured SOC (%)",
       ylab = "Predicted SOC (%)",
       xlim = c(0,20), # Specify range for SOC visualization
       ylim = c(0,20))
abline(coef = c(0,1))
vali.result <- vali.result[,c("Core","Med_depth","vali.i","preds.i")]
r2.c
rmse.c
mbe.c 
r2.v
rmse.v
mbe.v
write.csv(vali.result, paste0("./SOC/Imp3/D1/",ROI,"/",ROI,"_",KSSL_num,"_",spike_num,"_pred.csv"))
# Medium depth
imp.cov <- data.frame()
preds.i <- data.frame()
Vali.individual.i <- data.frame()
rmse.f = 0
r2.f = 0
mbe.f = 0
for (i in 1:length(unique(Cali.i2$Cali_ID2))){
  Cali.KSSL <- Cali.i2[which(Cali.i2$Cali_ID2 == i),]
  Cali.individual <- rbind(Cali.KSSL,Spike2)
  Vali.individual <- Vali.i2[which(Vali.i2$Cali_ID2 == i),]
  rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                         Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                       data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
  r2 <- summary(lm(rf.i$predicted ~ Cali.individual$oc_usda))$adj.r.squared
  r2.f <- r2 + r2.f
  resids <- rf.i$predicted - Cali.individual$oc_usda
  rmse <- sqrt(mean(resids^2))
  rmse.f <- rmse + rmse.f
  mbe <- mean(resids)
  mbe.f <- mbe + mbe.f
  imp <- importance(rf.i, type = 1)[,1]
  imp.cov <- rbind(imp.cov, imp)
  preds <- predict(rf.i, newdata = Vali.individual)
  preds.i <- rbind(preds.i,preds) 
  Vali.individual.i <- rbind(Vali.individual.i,Vali.individual[,c("Num","Med_depth","Measured_SOC")])
}
r2.c <- r2.f/length(unique(Vali.i2$Cali_ID2))
rmse.c <- rmse.f/length(unique(Vali.i2$Cali_ID2))
mbe.c <- mbe.f/length(unique(Vali.i2$Cali_ID2))
colnames(imp.cov) <- names(imp)
imp.cov <- sapply(imp.cov, mean)
write.csv(imp.cov,paste0("./SOC/Imp3/D2/",ROI,"/",ROI,"_",KSSL_num,"_",spike_num,"_imp.csv"))
dotchart(sort(imp.cov), xlab="Variable importance (%)")
colnames(preds.i) <- c("preds.i")
colnames(Vali.individual.i) <- c("Core","Med_depth","vali.i")
vali.result <- cbind(preds.i,Vali.individual.i)
r2.v <- summary(lm(vali.result$preds.i~vali.result$vali.i))
rmse.v <- sqrt(mean((vali.result$preds.i-vali.result$vali.i)^2))
mbe.v <- mean(vali.result$preds.i-vali.result$vali.i)
plot(vali.result$vali.i, vali.result$preds.i, 
     main = "Measured vs. predicted SOC of the validation set \n
       (KSSL spiked with local samples for individual-based model)", cex.main = 0.95,
     xlab = "Measured SOC (%)",
     ylab = "Predicted SOC (%)",
     xlim = c(0,20), # Specify range for SOC visualization
     ylim = c(0,20))
abline(coef = c(0,1))
vali.result <- vali.result[,c("Core","Med_depth","vali.i","preds.i")]
r2.c
rmse.c
mbe.c 
r2.v
rmse.v
mbe.v
write.csv(vali.result, paste0("./SOC/Imp3/D2/",ROI,"/",ROI,"_",KSSL_num,"_",spike_num,"_pred.csv"))
# Bottom depth
imp.cov <- data.frame()
preds.i <- data.frame()
Vali.individual.i <- data.frame()
rmse.f = 0
r2.f = 0
mbe.f = 0
for (i in 1:length(unique(Cali.i3$Cali_ID3))){
  Cali.KSSL <- Cali.i3[which(Cali.i3$Cali_ID3 == i),]
  Cali.individual <- rbind(Cali.KSSL,Spike3)
  Vali.individual <- Vali.i3[which(Vali.i3$Cali_ID3 == i),]
  rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                         Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                       data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
  r2 <- summary(lm(rf.i$predicted ~ Cali.individual$oc_usda))$adj.r.squared
  r2.f <- r2 + r2.f
  resids <- rf.i$predicted - Cali.individual$oc_usda
  rmse <- sqrt(mean(resids^2))
  rmse.f <- rmse + rmse.f
  mbe <- mean(resids)
  mbe.f <- mbe + mbe.f
  imp <- importance(rf.i, type = 1)[,1]
  imp.cov <- rbind(imp.cov, imp)
  preds <- predict(rf.i, newdata = Vali.individual)
  preds.i <- rbind(preds.i,preds) 
  Vali.individual.i <- rbind(Vali.individual.i,Vali.individual[,c("Num","Med_depth","Measured_SOC")])
}
r2.c <- r2.f/length(unique(Vali.i3$Cali_ID3))
rmse.c <- rmse.f/length(unique(Vali.i3$Cali_ID3))
mbe.c <- mbe.f/length(unique(Vali.i3$Cali_ID3))
colnames(imp.cov) <- names(imp)
imp.cov <- sapply(imp.cov, mean)
write.csv(imp.cov,paste0("./SOC/Imp3/D3/",ROI,"/",ROI,"_",KSSL_num,"_",spike_num,"_imp.csv"))
dotchart(sort(imp.cov), xlab="Variable importance (%)")
colnames(preds.i) <- c("preds.i")
colnames(Vali.individual.i) <- c("Core","Med_depth","vali.i")
vali.result <- cbind(preds.i,Vali.individual.i)
r2.v <- summary(lm(vali.result$preds.i~vali.result$vali.i))
rmse.v <- sqrt(mean((vali.result$preds.i-vali.result$vali.i)^2))
mbe.v <- mean(vali.result$preds.i-vali.result$vali.i)
plot(vali.result$vali.i, vali.result$preds.i, 
     main = "Measured vs. predicted SOC of the validation set \n
       (KSSL spiked with local samples for individual-based model)", cex.main = 0.95,
     xlab = "Measured SOC (%)",
     ylab = "Predicted SOC (%)",
     xlim = c(0,20), # Specify range for SOC visualization
     ylim = c(0,20))
abline(coef = c(0,1))
vali.result <- vali.result[,c("Core","Med_depth","vali.i","preds.i")]
r2.c
rmse.c
mbe.c 
r2.v
rmse.v
mbe.v
write.csv(vali.result, paste0("./SOC/Imp3/D3/",ROI,"/",ROI,"_",KSSL_num,"_",spike_num,"_pred.csv"))


#######################################################################
# SOC modeling results by combining different depth layers
KSSL_local1 <- read.table(paste0("./SOC/Imp3/D1/",ROI,"/",ROI,"_",KSSL_num,"_",spike_num,"_pred.csv"),
                          comment.char ="", quote = "\"", header = T, sep = ",")
KSSL_local2 <- read.table(paste0("./SOC/Imp3/D2/",ROI,"/",ROI,"_",KSSL_num,"_",spike_num,"_pred.csv"),
                          comment.char ="", quote = "\"", header = T, sep = ",")
KSSL_local3 <- read.table(paste0("./SOC/Imp3/D3/",ROI,"/",ROI,"_",KSSL_num,"_",spike_num,"_pred.csv"),
                          comment.char ="", quote = "\"", header = T, sep = ",")
vali.result <- rbind(KSSL_local1,KSSL_local2,KSSL_local3)
r2.v <- summary(lm(vali.result$preds.i~vali.result$vali.i))
rmse.v <- sqrt(mean((vali.result$preds.i-vali.result$vali.i)^2))
mbe.v <- mean(vali.result$preds.i-vali.result$vali.i)
r2.v
rmse.v
mbe.v
plot(vali.result$vali.i, vali.result$preds.i,
     main = "Measured vs. predicted SOC of the validation set \n
     (KSSL spiked with local samples for individual-based model)", cex.main = 0.95,
     xlab = "Measured SOC (%)",
     ylab = "Predicted SOC (%)",
     xlim = c(0,20), # Specify range for SOC visualization
     ylim = c(0,20))
abline(coef = c(0,1))