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
BD <- read.table("CSV-with-measured-BD-filename.csv", comment.char ="", quote = "\"", header = T, sep = ",")
BD <- BD[,c("Sample_ID","BD")]
colnames(BD)[2] <- "BD_m"
KSSL <- read.table('./KSSL.f.w.BD.csv', comment.char ="", quote = "\"", header = T, sep = ",")
Vali <- read.table(paste0('./',ROI,'_vali_0.csv') # Choose all samples from Year 1 for validation
                   , comment.char ="", quote = "\"", header = T, sep = ",")
Vali <- merge(Vali, BD, by="Sample_ID")
CaliAll <- read.table(paste0('./BD/BS/', ROI,'/', ROI,'_calibration_0_', KSSL_num,'.csv')
                        , comment.char ="", quote = "\"", header = T, sep = ",")
Spike_all <- read.table(paste0('./',ROI,'CSV-with-local-spiking-samples-joined-with-covariates.csv.csv'), comment.char ="", quote = "\"", header = T, sep = ",")
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
Cali.i <- Cali.i[,c("sampleids","lay_id","Med_depth","bd_kssl","BD",
                      "SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                      "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","Cali_ID")]
Spike <- Spike[,c("Sample_ID","X","Med_depth","BD_m","BD",
                    "SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                    "LULC","Tree","GPP","EVI","AS","EL","SL","TWI")]
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
    rf.i <- randomForest(bd_kssl ~ Med_depth + BD + Clay + Sand + SOC + ppt + tmin + tmax + VPD +
                           Tree + GPP + EVI + AS + EL + SL + TWI,
                         data = Cali.individual, importance = TRUE, ntree=500, mtry=3)
    r2 <- summary(lm(rf.i$predicted ~ Cali.individual$bd_kssl))$adj.r.squared
    r2.f <- r2 + r2.f
    resids <- rf.i$predicted - Cali.individual$bd_kssl
    rmse <- sqrt(mean(resids^2))
    rmse.f <- rmse + rmse.f
    mbe <- mean(resids)
    mbe.f <- mbe + mbe.f
    imp <- importance(rf.i, type = 1)[,1]
    imp.cov <- rbind(imp.cov, imp)
    preds <- predict(rf.i, newdata = Vali.individual)
    preds.i <- rbind(preds.i,preds) 
    Vali.individual.i <- rbind(Vali.individual.i,Vali.individual[,c("Num","Med_depth","BD_m")])
}
r2.c <- r2.f/length(unique(Vali.i$Cali_ID))
rmse.c <- rmse.f/length(unique(Vali.i$Cali_ID))
mbe.c <- mbe.f/length(unique(Vali.i$Cali_ID))
colnames(imp.cov) <- names(imp)
imp.cov <- sapply(imp.cov, mean)
write.csv(imp.cov,paste0("2022_calibration/Results/BD/Imp2/",ROI,"/",ROI,"_",KSSL_num,"_",spike_num,"_imp.csv"))
dotchart(sort(imp.cov), xlab="Variable importance (%)")
colnames(preds.i) <- c("preds.i")
colnames(Vali.individual.i) <- c("Core","Med_depth","vali.i")
vali.result <- cbind(preds.i,Vali.individual.i)
r2.v <- summary(lm(vali.result$preds.i~vali.result$vali.i))
rmse.v <- sqrt(mean((vali.result$preds.i-vali.result$vali.i)^2))
mbe.v <- mean(vali.result$preds.i-vali.result$vali.i)
plot(vali.result$vali.i, vali.result$preds.i, 
       main = "Measured vs. predicted BD of the validation set \n
       (KSSL spiked with local samples for individual-based model)", cex.main = 0.95,
       xlab = "Measured BD (g/cm3)",
       ylab = "Predicted BD (g/cm3)",
       xlim = c(0,2), # Specify range for BD visualization
       ylim = c(0,2))
abline(coef = c(0,1))
vali.result <- vali.result[,c("Core","Med_depth","vali.i","preds.i")]
r2.c
rmse.c
mbe.c 
r2.v
rmse.v
mbe.v
write.csv(vali.result, paste0("./BD/Imp2/",ROI,"/",ROI,"_",KSSL_num,"_",spike_num,".csv"))
