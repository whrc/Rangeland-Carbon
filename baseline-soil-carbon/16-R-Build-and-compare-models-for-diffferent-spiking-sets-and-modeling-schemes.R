library(tidyverse)
library(randomForest)
setwd("C:/Users/yxia/Desktop/Optimized_field_sampling")
setwd("set-pathway")


###############################################################
# Define parameters, site, and import datasets
ROI = "specify-site-name"
spike_num = 25
KSSL_num = 500
# Import files for model calibration and validation
KSSL <- read.table('KSSL/KSSL.f.csv', comment.char ="", quote = "\"", header = T, sep = ",")
CaliAll<- read.table(paste0("./",ROI,"_calibration_",spike_num,"_",KSSL_num,".csv"), 
                       comment.char ="", quote = "\"", header = T, sep = ",") # Get ID for calibration dataset
BD <- BD[,c("Sample_ID","BD")] # Import BD dataset
colnames(BD)[2] <- "Measured_BD"
if (spike_num == 0){
  Spike <- list()
  spikelist <- list()
} else {
  Spike <- read.table(paste0("./",ROI,"_spike_",spike_num,".csv"), 
                      comment.char ="", quote = "\"", header = T, sep = ",")
  spikelist <- unique(Spike$Num)
} # Get local spiking dataset
Vali<- read.table("CSV-corresponding-to-site-specific-validation-dataset-not-used-for-spiking.csv", 
                    comment.char ="", quote = "\"", header = T, sep = ",") 
Cali <- merge(CaliAll,KSSL,by = "sampleids")
Cali.i <- within(Cali, Cali_ID <- match(predid, unique(predid)))
ID <- Cali.i[,c("predid","Cali_ID")] 
ID <- unique(ID)
Vali.i <- merge(Vali,ID,by.x = "Sample_ID",by.y="predid")
Vali.i <- merge(Vali.i, BD, by="Sample_ID")
Cali.i <- Cali.i[,c("sampleids","lay_id","Med_depth","oc_usda","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                    "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE","Cali_ID")]
dim(Spike)[1] # Num of spiking samples
dim(Vali)[1] # Num of validation samples
# Reorganize spiking samples for subsequent combination
if (spike_num != 0)
{Spike <- Spike[,c("Sample_ID","Num","Med_depth","Measured_SOC","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                   "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE")]
Spike$Cali_ID <- 1 # Placeholder since all spike samples will be used
colnames(Spike) <- colnames(Cali.i)}


###############################################################
# Run model with or without spiking samples
imp.cov <- data.frame()
preds.i <- data.frame()
Vali.individual.i <- data.frame()
rmse.f = 0
r2.f = 0
for (i in 1:length(unique(Cali.i$Cali_ID))){
  Cali.KSSL <- Cali.i[which(Cali.i$Cali_ID == i),]
  if (spike_num != 0) {Cali.individual <- rbind(Cali.KSSL,Spike)} # Model with spiking samples
  else {Cali.individual <- Cali.KSSL} # Model without use of spiking samples
  Vali.individual <- Vali.i[which(Vali.i$Cali_ID == i),]
  rf.i <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                         Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                       data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
  r2 <- summary(lm(rf.i$predicted ~ Cali.individual$oc_usda))$adj.r.squared
  r2.f <- r2 + r2.f
  resids <- rf.i$predicted - Cali.individual$oc_usda
  rmse <- sqrt(mean(resids^2))
  rmse.f <- rmse + rmse.f
  imp <- importance(rf.i, type = 1)[,1]
  imp.cov <- rbind(imp.cov, imp)
  preds <- predict(rf.i, newdata = Vali.individual)
  preds.i <- rbind(preds.i,preds) 
  Vali.individual.i <- rbind(Vali.individual.i,Vali.individual[,c("Num","Med_depth","Measured_SOC","Measured_BD")])
}
r2.f <- r2.f/length(unique(Vali.i$Cali_ID))
rmse.f <- rmse.f/length(unique(Vali.i$Cali_ID))
r2.f
rmse.f
colnames(imp.cov) <- names(imp)
imp.cov <- sapply(imp.cov, mean)
dotchart(sort(imp.cov), xlab="Variable importance (%)")
colnames(preds.i) <- c("preds.i")
colnames(Vali.individual.i) <- c("Core","Med_depth","vali.i","BD")
vali.result <- cbind(preds.i,Vali.individual.i)
summary(lm(vali.result$preds.i~vali.result$vali.i))
sqrt(mean((vali.result$preds.i-vali.result$vali.i)^2))
plot(vali.result$vali.i, vali.result$preds.i, 
     main = "Measured vs. predicted SOC of the validation set (individual-based model)", cex.main = 0.95,
     xlab = "Measured SOC (%)",
     ylab = "Predicted SOC (%)",
     xlim = c(0,20), #Define SOC% range for visualization
     ylim = c(0,20))
abline(coef = c(0,1))
# Summarize results based on soil cores
vali.result$BD <- as.numeric(vali.result$BD)
vali.result$SOCS_pred <- vali.result$preds.i*15*vali.result$BD*1.102311
vali.result$SOCS_measured <- vali.result$vali.i*15*vali.result$BD*1.102311
vali_stock_pred <- aggregate(vali.result$SOCS_pred, by=list(vali.result$Core), FUN=sum)
vali_stock_measured <- aggregate(vali.result$SOCS_measured, by=list(vali.result$Core), FUN=sum)
vali_stock <- data.frame(vali_stock_pred,vali_stock_measured)[,c(1,2,4)]
colnames(vali_stock) <- c("Core","pred","measured")
summary(lm(vali_stock$pred~vali_stock$measured))
sqrt(mean((vali_stock$pred-vali_stock$measured)^2))
plot(vali_stock$measured, vali_stock$pred, 
     main = "Measured vs. predicted SOC stocks of the validation set (individual-based model)", cex.main = 0.95,
     xlab = "Measured SOC stocks (t/ha)",
     ylab = "Predicted SOC stocks (t/ha)",
     xlim = c(0,200),
     ylim = c(0,200))
abline(coef = c(0,1))


###############################################################
# Run local model without using KSSL samples
rf <- randomForest(oc_usda ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                     Tree + GPP + EVI + AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                   data = Spike, importance = TRUE, ntree=100, mtry=3)
r2.c <- summary(lm(rf$predicted ~ Spike$oc_usda))$adj.r.squared
r2.c
resids.c <- rf$predicted - Spike$oc_usda
rmse.c <- sqrt(mean(resids.c^2))
rmse.c
imp <- importance(rf.i, type = 1)[,1]
dotchart(sort(imp), xlab="Variable importance (%)")
preds <- predict(rf, newdata = Vali.i)
r2.v <- summary(lm(preds ~ Vali.i$Measured_SOC))$adj.r.squared
r2.v
resids.v <- preds - Vali.i$Measured_SOC
rmse.v <- sqrt(mean(resids.v^2))
rmse.v