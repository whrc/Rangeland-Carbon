library(tidyverse)
library(randomForest)
library(clhs)
setwd("set-pathway")


###############################################################
# Define parameters, site, and import datasets
ROI = "specify-site-name"
spike_num = 36 # Choose number of spiking cores
KSSL_num = 500 # Choose number of calibration samples
# Import files for model calibration and validation
KSSL <- read.table('./KSSL.f.2022.csv', comment.char ="", quote = "\"", header = T, sep = ",")
Vali <- read.table(paste0('./',ROI,'_vali_0.csv') # Choose all samples from Year 1 for validation
                   , comment.char ="", quote = "\"", header = T, sep = ",")
CaliAll <- read.table(paste0('./SOC/', ROI,'/', ROI,'_calibration_0_', KSSL_num,'.csv')
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
# Merge datasets to get KSSL calibration data for each point in the pixel-based model
Cali <- merge(CaliAll,KSSL,by = "sampleids")
Cali.i <- within(Cali, Cali_ID <- match(predid, unique(predid)))
ID <- Cali.i[,c("predid","Cali_ID")] 
ID <- unique(ID) 
Vali.i <- merge(Vali,ID,by.x = "Sample_ID",by.y="predid")
Cali.i <- Cali.i[,c("sampleids","lay_id","Med_depth","oc_usda","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                      "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE","Cali_ID")]
Spike <- Spike[,c("Sample_ID","X","Med_depth","SOC_m","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                    "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE")]
Spike$Cali_ID = 0 # Placeholder for spiking samples
colnames(Spike) <- colnames(Cali.i)


#######################################################################
# Run model with spiking: combine KSSL and local samples
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
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
    #                        Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M1 - without depth
    rf.i <- randomForest(oc_usda ~ BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                           Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                         data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M2 - without RAP
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
    #                        Tree + GPP + EVI + AS + EL + SL + TWI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M3 - without climate
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC +
    #                        Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M4 - without topography
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
    #                        Tree + GPP + EVI + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M5 - without RS
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
    #                        AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M6 - without soil
    # rf.i <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD +
    #                        Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M7 - depth and soil
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M8 - depth and climate
    # rf.i <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M9 - depth and RS
    # rf.i <- randomForest(oc_usda ~ Med_depth + Tree + GPP + EVI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M10 - depth and topography
    # rf.i <- randomForest(oc_usda ~ Med_depth + AS + EL + SL + TWI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M11 - depth and RAP
    # rf.i <- randomForest(oc_usda ~ Med_depth + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M12 - depth, soil, and climate
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M13 - depth, soil, and RS
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + Tree + GPP + EVI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M14 - depth, soil, and RAP
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M15 - depth, soil, and topography
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + AS + EL + SL + TWI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M16 - depth, climate, and RS
    # rf.i <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD + Tree + GPP + EVI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M17 - depth, climate, and RAP
    # rf.i <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M18 - depth, climate, and topography
    # rf.i <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD + AS + EL + SL + TWI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M19 - depth, RS, and RAP
    # rf.i <- randomForest(oc_usda ~ Med_depth + Tree + GPP + EVI + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M20 - depth, RS, and topography
    # rf.i <- randomForest(oc_usda ~ Med_depth + Tree + GPP + EVI + AS + EL + SL + TWI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M21 - depth, RAP, and topography
    # rf.i <- randomForest(oc_usda ~ Med_depth + AFGC + BG + LTR + PFGC + SHR + TREE + AS + EL + SL + TWI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
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
r2.c <- r2.f/length(unique(Vali.i$Cali_ID))
rmse.c <- rmse.f/length(unique(Vali.i$Cali_ID))
mbe.c <- mbe.f/length(unique(Vali.i$Cali_ID))
colnames(imp.cov) <- names(imp)
imp.cov <- sapply(imp.cov, mean)
write.csv(imp.cov,paste0("./SOC/Imp1/",ROI,"/",ROI,"_",KSSL_num,"_",spike_num,"_imp.csv"))
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
write.csv(vali.result, paste0("./SOC/Imp1/",ROI,"/",ROI,"_",KSSL_num,"_",spike_num,"_pred.csv"))


#################################################################################
# Run model without spiking: KSSL only 
imp.cov <- data.frame()
preds.i <- data.frame()
Vali.individual.i <- data.frame()
rmse.f = 0
r2.f = 0
mbe.f = 0
for (i in 1:length(unique(Cali.i$Cali_ID))){
    Cali.individual <- Cali.i[which(Cali.i$Cali_ID == i),]
    Vali.individual <- Vali.i[which(Vali.i$Cali_ID == i),]
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
    #                        Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M1 - without depth
    rf.i <- randomForest(oc_usda ~ BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                           Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                         data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M2 - without RAP
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
    #                        Tree + GPP + EVI + AS + EL + SL + TWI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M3 - without climate
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC +
    #                        Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M4 - without topography
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
    #                        Tree + GPP + EVI + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M5 - without RS
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
    #                        AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M6 - without soil
    # rf.i <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD +
    #                        Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M7 - depth and soil
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M8 - depth and climate
    # rf.i <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M9 - depth and RS
    # rf.i <- randomForest(oc_usda ~ Med_depth + Tree + GPP + EVI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M10 - depth and topography
    # rf.i <- randomForest(oc_usda ~ Med_depth + AS + EL + SL + TWI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M11 - depth and RAP
    # rf.i <- randomForest(oc_usda ~ Med_depth + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M12 - depth, soil, and climate
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M13 - depth, soil, and RS
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + Tree + GPP + EVI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M14 - depth, soil, and RAP
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M15 - depth, soil, and topography
    # rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + AS + EL + SL + TWI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M16 - depth, climate, and RS
    # rf.i <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD + Tree + GPP + EVI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M17 - depth, climate, and RAP
    # rf.i <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M18 - depth, climate, and topography
    # rf.i <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD + AS + EL + SL + TWI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M19 - depth, RS, and RAP
    # rf.i <- randomForest(oc_usda ~ Med_depth + Tree + GPP + EVI + AFGC + BG + LTR + PFGC + SHR + TREE,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M20 - depth, RS, and topography
    # rf.i <- randomForest(oc_usda ~ Med_depth + Tree + GPP + EVI + AS + EL + SL + TWI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    # M21 - depth, RAP, and topography
    # rf.i <- randomForest(oc_usda ~ Med_depth + AFGC + BG + LTR + PFGC + SHR + TREE + AS + EL + SL + TWI,
    #                      data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
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
r2.c <- r2.f/length(unique(Vali.i$Cali_ID))
rmse.c <- rmse.f/length(unique(Vali.i$Cali_ID))
mbe.c <- mbe.f/length(unique(Vali.i$Cali_ID))
colnames(imp.cov) <- names(imp)
imp.cov <- sapply(imp.cov, mean)
write.csv(imp.cov,paste0("./SOC/Imp1/",ROI,"/",ROI,"_",KSSL_num,"_0_imp.csv"))
dotchart(sort(imp.cov), xlab="Variable importance (%)")
colnames(preds.i) <- c("preds.i")
colnames(Vali.individual.i) <- c("Core","Med_depth","vali.i")
vali.result <- cbind(preds.i,Vali.individual.i)
r2.v <- summary(lm(vali.result$preds.i~vali.result$vali.i))
rmse.v <- sqrt(mean((vali.result$preds.i-vali.result$vali.i)^2))
mbe.v <- mean(vali.result$preds.i-vali.result$vali.i)
plot(vali.result$vali.i, vali.result$preds.i, 
       main = "Measured vs. predicted SOC of the validation set \n
       (KSSL only individual-based model)", cex.main = 0.95,
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
write.csv(vali.result, paste0("./SOC/Imp1/",ROI,"/",ROI,"_",KSSL_num,"_0_pred.csv"))


#################################################################################
# Run model with local samples only
# rf <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
#                      Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
#                    data = Spike, importance = TRUE, ntree=500, mtry=3)
# M1 - without depth
rf <- randomForest(oc_usda ~ BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                       Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                     data = Spike, importance = TRUE, ntree=500, mtry=3)
# M2 - without RAP
# rf <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
#                        Tree + GPP + EVI + AS + EL + SL + TWI,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M3 - without climate
# rf <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC +
#                        Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M4 - without topography
# rf <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
#                        Tree + GPP + EVI + AFGC + BG + LTR + PFGC + SHR + TREE,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M5 - without RS
# rf <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
#                        AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M6 - without soil
# rf <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD +
#                        Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M7 - depth and soil
# rf <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M8 - depth and climate
# rf <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M9 - depth and RS
# rf <- randomForest(oc_usda ~ Med_depth + Tree + GPP + EVI,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M10 - depth and topography
# rf <- randomForest(oc_usda ~ Med_depth + AS + EL + SL + TWI,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M11 - depth and RAP
# rf <- randomForest(oc_usda ~ Med_depth + AFGC + BG + LTR + PFGC + SHR + TREE,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M12 - depth, soil, and climate
# rf <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M13 - depth, soil, and RS
# rf <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + Tree + GPP + EVI,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M14 - depth, soil, and RAP
# rf <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + AFGC + BG + LTR + PFGC + SHR + TREE,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M15 - depth, soil, and topography
# rf <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + AS + EL + SL + TWI,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M16 - depth, climate, and RS
# rf <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD + Tree + GPP + EVI,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M17 - depth, climate, and RAP
# rf <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD + AFGC + BG + LTR + PFGC + SHR + TREE,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M18 - depth, climate, and topography
# rf <- randomForest(oc_usda ~ Med_depth + ppt + tmin + tmax + VPD + AS + EL + SL + TWI,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M19 - depth, RS, and RAP
# rf <- randomForest(oc_usda ~ Med_depth + Tree + GPP + EVI + AFGC + BG + LTR + PFGC + SHR + TREE,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M20 - depth, RS, and topography
# rf <- randomForest(oc_usda ~ Med_depth + Tree + GPP + EVI + AS + EL + SL + TWI,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
# M21 - depth, RAP, and topography
# rf <- randomForest(oc_usda ~ Med_depth + AFGC + BG + LTR + PFGC + SHR + TREE + AS + EL + SL + TWI,
#                      data = Spike, importance = TRUE, ntree=500, mtry=3)
r2.c <- summary(lm(rf$predicted ~ Spike$oc_usda))$adj.r.squared
r2.c
resids.c <- rf$predicted - Spike$oc_usda
rmse.c <- sqrt(mean(resids.c^2))
rmse.c
mbe.c <- mean(rf$predicted - Spike$oc_usda)
mbe.c
imp <- importance(rf, type = 1)[,1]
write.csv(imp, paste0("./SOC/Imp1/",ROI,"/",ROI,"_0_",spike_num,"_imp.csv"))
dotchart(sort(imp), xlab="Variable importance (%)")
preds <- predict(rf, newdata = Vali)
r2.v <- summary(lm(preds ~ Vali$Measured_SOC))$adj.r.squared
r2.v
resids.v <- preds - Vali$Measured_SOC
rmse.v <- sqrt(mean(resids.v^2))
rmse.v
mbe.v <- mean(preds - Vali$Measured_SOC)
mbe.v
plot(Vali$Measured_SOC, preds, 
     main = "Measured vs. predicted SOC of the validation set \n
     (local only model)", cex.main = 0.95,
     xlab = "Measured SOC (%)",
     ylab = "Predicted SOC (%)",
     xlim = c(0,20), # Specify range for SOC visualization
     ylim = c(0,20))
abline(coef = c(0,1))
vali.result <- cbind(Vali[,c("Num","Med_depth","Measured_SOC")],preds)
colnames(vali.result) <- c("Core","Med_depth","vali.i","preds.i")
write.csv(vali.result, paste0("./SOC/Imp1/",ROI,"/",ROI,"_0_",spike_num,"_pred.csv"))
