library(randomForest)
library(outliers)
library(caret)
library(tidyverse)
setwd("set-pathway")


###############################################################
# Set parameters, modeling approach, and import datasets
KSSL_num_site = 500 
KSSL_num_pixel = 500 
#0: Calibration datasets selected with site covariates only;
#1: Calibration datasets selected with MIR and site covariates;
model_scheme = 1 
# Import calibration and validation datasets
CaliAll<- read.table('./KSSL.f.2022.csv', comment.char ="", quote = "\"", header = T, sep = ",")
ValiAll<- read.table('./CSV-with-site-covariates-filename.csv',comment.char ="", quote = "\"", header = T, sep = ",")
ValiAll$Measured_SOC <- ValiAll$Measured_SOC/10 # Convert to %, if not yet converted in the spreadsheet
CaliAll <- na.omit(CaliAll) # Remove data points outside of western U.S. in order to use RAP covariates, if not yet done in previous steps
if (model_scheme == 0){
  # Site-based calibration datasets selected with site covariates only  
  site.calibration <- read.table(paste0("./site_calibration_0_",KSSL_num_site,"_site_cov.csv"),
                                comment.char ="", quote = "\"", header = T, sep = ",")
} else {
  # Site-based calibration datasets selected with MIR and covariates including RAP
  site.calibration <- read.table(paste0("./site_calibration_0_",KSSL_num_site,"_site_MIR_cov.csv"),
                                comment.char ="", quote = "\"", header = T, sep = ",")}

if (model_scheme == 0){
  # Pixel-based calibration datasets selected with site covariates only    
  site.calibration.i <- read.table(paste0("./site_calibration_0_",KSSL_num_pixel,"_pixel_cov.csv"),
                                  comment.char ="", quote = "\"", header = T, sep = ",")
} else {
  # Pixel-based calibration datasets selected with MIR and covariates including RAP
  site.calibration.i <- read.table(paste0("./site_calibration_0_",KSSL_num_pixel,"_pixel_MIR_cov.csv"),
                                  comment.char ="", quote = "\"", header = T, sep = ",")
}
CaliAll.site <- merge(site.calibration,CaliAll,by = "sampleids")
CaliAll.site.i <- merge(site.calibration.i,CaliAll,by = "sampleids")


###############################################################
# Run RF for site-based models
# Find optimal parameters for RF
set.seed(123)
control <- trainControl(method="repeatedcv", number=10, repeats=3, search="random")
metric <- "Rsquared"
mtry <- sqrt(ncol(CaliAll.site))
rf_site <- train(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                data = CaliAll.site, method="rf", metric=metric, tuneLength=15, trControl=control)
print(rf_site)
plot(rf_site)
# Use optimized parameters to run the model at the site level
rf_site <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                         Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE, 
                       data = CaliAll.site, importance = TRUE, ntree=500, mtry=3)
# Calculate error metrics for the calibration set
summary(lm(rf_site$predicted ~ CaliAll.site$oc_usda))
resids <- rf_site$predicted - CaliAll.site$oc_usda
rmse <- sqrt(mean(resids^2))
rmse
mbe <- mean(resids)
mbe
plot(CaliAll.site$oc_usda, rf_site$predicted, 
     main = "Measured vs. predicted SOC of the calibration set", 
     cex.main = 0.95,
     xlab = "Measured SOC (%)",
     ylab = "Predicted SOC (%)",
     xlim = c(0,10), # Set to fit data range for SOC%
     ylim = c(0,10))
abline(coef = c(0,1))
# Generate covariate importance based on the calibration model
impToPlot <- importance(rf_site, type = 1)
dotchart(sort(impToPlot[,1]), xlim=c(0,60), xlab="Variable importance (%)")
impToPlot[,1]
if (model_scheme == 0){
  write.csv(impToPlot,paste0("./site_calibration_0_",KSSL_num_site,"_site_cov_imp.csv"))
} else {
  write.csv(impToPlot,paste0("./site_calibration_0_",KSSL_num_site,"_site_MIR_cov_imp.csv"))
}  
# Prediction on the validation set
length(unique(CaliAll.site$sampleids))
preds.site <- predict(rf_site, newdata = ValiAll)
summary(lm(preds.site~ValiAll$Measured_SOC))
sqrt(mean((preds.site-ValiAll$Measured_SOC)^2))
mean(preds.site-ValiAll$Measured_SOC)
plot(ValiAll$Measured_SOC, preds.site, 
     main = "Measured vs. predicted SOC of the validation set", cex.main = 0.95,
     xlab = "Measured SOC (%)",
     ylab = "Predicted SOC (%)",
     xlim = c(0.5,10),
     ylim = c(0.5,10))
abline(coef = c(0,1))
site_result <- data.frame(ValiAll$Num,ValiAll$whrc_id,ValiAll$Med_depth,preds.site,ValiAll$Measured_SOC)
colnames(site_result) <- c("Num","whrc_id","Med_depth","Predicted_SOC","Measured_SOC")
if (model_scheme == 0){
  write.csv(site_result,paste0("./site_",KSSL_num_site,"_site_cov_pred.csv"))
} else {
  write.csv(site_result,paste0("./site_",KSSL_num_site,"_site_MIR_cov_pred.csv"))
}  


###############################################################
# Run RF for individual-based models
colnames(CaliAll.site.i)[XX] <- c("predid") # Rename column XX if needed to keep the column names consistent
CaliAll.site.i <- within(CaliAll.site.i, Cali_ID <- match(predid, unique(predid)))
CaliAll.site.i.reorg <-CaliAll.site.i[,c("predid","Cali_ID")] 
colnames(CaliAll.site.i.reorg)[1] <- "whrc_id"
ValiAll.site.i <- merge(ValiAll,CaliAll.site.i.reorg,by = "whrc_id")
ValiAll.site.i <- ValiAll.site.i %>% distinct(whrc_id, .keep_all = TRUE)
length(unique(CaliAll.site.i$sampleids))
# Initialize variables for calculating average
imp.cov <- data.frame()
preds.i <- data.frame()
Vali.individual.i <- data.frame()
rmse.f = 0
r2.f = 0
mbe.f = 0
# Loop the RF model
for (i in 1:length(unique(CaliAll.site.i$predid))){
  Cali.individual <- CaliAll.site.i[which(CaliAll.site.i$Cali_ID == i),]
  Vali.individual <- ValiAll.site.i[which(ValiAll.site.i$Cali_ID == i),]
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
  Vali.individual.i <- rbind(Vali.individual.i,Vali.individual[,c("Num","whrc_id","Med_depth","Measured_SOC")])
}
# Report averaged results from model calibration
r2.c <- r2.f/length(unique(CaliAll.site.i$predid))
rmse.c <- rmse.f/length(unique(CaliAll.site.i$predid))
mbe.c <- mbe.f/length(CaliAll.site.i$predid)
r2.c
rmse.c
mbe.c
colnames(imp.cov) <- names(imp)
imp.cov <- sapply(imp.cov, mean)
if (model_scheme == 0){
  write.csv(imp.cov,paste0("./site_calibration_0_",KSSL_num_pixel,"_pixel_cov_imp.csv"))
} else {
  write.csv(imp.cov,paste0("./site_calibration_0_",KSSL_num_pixel,"_pixel_MIR_cov_imp.csv"))
}  
dotchart(sort(imp.cov), xlim=c(0,60), xlab="Variable importance (%)")
colnames(preds.i) <- c("preds.i")
colnames(Vali.individual.i) <- c("Num","whrc_id","Med_depth","vali.i")
vali.result <- cbind(Vali.individual.i$Num,Vali.individual.i$whrc_id,Vali.individual.i$Med_depth,preds.i,Vali.individual.i$vali.i)
colnames(vali.result) <- c("Num","whrc_id","Med_depth","Predicted_SOC","Measured_SOC")
summary(lm(vali.result$Predicted_SOC~vali.result$Measured_SOC))
sqrt(mean((vali.result$Predicted_SOC-vali.result$Measured_SOC)^2))
mean(vali.result$Predicted_SOC-vali.result$Measured_SOC)
plot(vali.result$Measured_SOC, vali.result$Predicted_SOC, 
     main = "Measured vs. predicted SOC of the validation set (individual-based model)", cex.main = 0.95,
     xlab = "Measured SOC (%)",
     ylab = "Predicted SOC (%)",
     xlim = c(0.5,10),
     ylim = c(0.5,10))
abline(coef = c(0,1))
if (model_scheme == 0){
  write.csv(site_result,paste0("./site_",KSSL_num_site,"_pixel_cov_pred.csv"))
} else {
  write.csv(vali.result,paste0("./site_",KSSL_num_pixel,"_pixel_MIR_cov_pred.csv"))
}  

