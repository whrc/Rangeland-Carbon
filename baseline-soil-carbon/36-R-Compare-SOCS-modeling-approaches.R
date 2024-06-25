library(dplyr)
setwd("set-pathway")


###############################################################
# Define parameters, site, and import datasets
ROI = "specify-site-name"
spike_num = 36 # Choose number of spiking cores
KSSL_num = 500 # Choose number of calibration samples
model_scheme = 1 # 0: Site-based model; 1: Pixel-based model
if (model_scheme == 0) {scheme = "site"} else {scheme = "pixel"}
# Import files for summarizing modeling results
Vali_SOCS <- read.table(paste0('./',ROI,'CSV-with-measured-SOCS-filename-in-Year1.csv')
                   , comment.char ="", quote = "\"", header = T, sep = ",")
SpikeSOCS <- read.table(paste0('./',ROI,'CSV-with-measured-SOCS-filename-in-Year2.csv'), 
                        comment.char ="", quote = "\"", header = T, sep = ",")
Cluster <- unique(Vali_SOCS[,c("Num","Cluster")])
KSSL_local_soc <- read.table(paste0("./SOC/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_",spike_num,".csv"),
                             comment.char ="", quote = "\"", header = T, sep = ",")
KSSL_local_bd <- read.table(paste0("./BD/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_",spike_num,".csv"),
                            comment.char ="", quote = "\"", header = T, sep = ",")
KSSL_local_socs <- read.table(paste0("./SOCS/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_",spike_num,".csv"),
                              comment.char ="", quote = "\"", header = T, sep = ",")
KSSL_local_socd <- read.table(paste0("./SOCD/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_",spike_num,".csv"),
                              comment.char ="", quote = "\"", header = T, sep = ",")
KSSL_local_soc <- merge(KSSL_local_soc,Cluster,by.x = "Core", by.y = "Num")
KSSL_local_bd <- merge(KSSL_local_bd,Cluster,by.x = "Core", by.y = "Num")
KSSL_local_socs <- merge(KSSL_local_socs,Cluster,by.x = "Core", by.y = "Num")
KSSL_local_socd <- merge(KSSL_local_socd,Cluster,by.x = "Core", by.y = "Num")


###############################################################
# Method 1 for modeling SOC stocks directly
r2.v <- summary(lm(KSSL_local_socs$preds.i~KSSL_local_socs$vali.i))
rmse.v <- sqrt(mean((KSSL_local_socs$preds.i-KSSL_local_socs$vali.i)^2))
mbe.v <- mean(KSSL_local_socs$preds.i-KSSL_local_socs$vali.i)
r2.v
rmse.v
mbe.v
plot(KSSL_local_socs$vali.i, KSSL_local_socs$preds.i, 
     main = "Measured vs. predicted SOC stocks (0-30 cm) of the validation set \n
     (Modeled with KSSL-local individual-based model for SOC stocks)", cex.main = 0.95,
     xlab = "Measured SOCS (g/m2)",
     ylab = "Predicted SOCS (g/m2)",
     xlim = c(0,20000), # Specify range for SOCS visualization
     ylim = c(0,20000))
abline(coef = c(0,1))


###############################################################
# Method 2 for calculating SOCS: model BD and SOC separately
KSSL_local <- merge(KSSL_local_soc,KSSL_local_bd,by=c("Cluster","Core","Med_depth","X"))
KSSL_local$vali.socs <- KSSL_local$vali.i * KSSL_local$vali.i.y * 100 * 15
KSSL_local$pred.socs <- KSSL_local$preds.i * KSSL_local$preds.i.y * 100 * 15
KSSL_local_30 <- KSSL_local[which(KSSL_local$Med_depth <= 15),]
M2_SOCS <- KSSL_local_30 %>% group_by(Core) %>% 
  summarize(mean_vali = sum(vali.socs),mean_pred = sum(pred.socs),
            sd_vali = sd(vali.socs), sd_pred = sd(pred.socs))
M2_SOCS <- M2_SOCS[which(!is.na(M2_SOCS$sd_vali)),]
r2.M2 <- summary(lm(M2_SOCS$mean_vali ~ M2_SOCS$mean_pred))$r.squared
resids.M2 <- M2_SOCS$mean_pred - M2_SOCS$mean_vali
rmse.M2 <- sqrt(mean(resids.M2^2))
mbe.M2 <- mean(M2_SOCS$mean_pred - M2_SOCS$mean_vali)
plot(M2_SOCS$mean_vali, M2_SOCS$mean_pred, 
     main = "Measured vs. predicted SOC stocks (0-30 cm) of the validation set \n
     (KSSL-local individual-based models built separately for BD and SOC%)", cex.main = 0.95,
     xlab = "Measured SOCS (g/m2)",
     ylab = "Predicted SOCS (g/m2)",
     xlim = c(0,20000), # Specify range for SOCS visualization
     ylim = c(0,20000))
abline(coef = c(0,1))
r2.M2 
rmse.M2
mbe.M2 
write.csv(M2_SOCS, paste0("./SOCS/M2/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_",spike_num,".csv"))


###############################################################
# Method 3 for calculating SOCS: model SOCD before integrating results to SOCD
KSSL_local_30_SOCD <- KSSL_local_socd[which(KSSL_local_socd$Med_depth <= 15),]
M3_SOCS <- KSSL_local_30_SOCD %>% group_by(Core) %>% 
  summarize(mean_vali = mean(vali.i),mean_pred = mean(preds.i),
            sd_vali = sd(vali.i), sd_pred = sd(preds.i))
M3_SOCS$mean_vali <- M3_SOCS$mean_vali*30*10000 
M3_SOCS$mean_pred <- M3_SOCS$mean_pred*30*10000 
M3_SOCS <- M3_SOCS[which(!is.na(M3_SOCS$sd_vali)),]
r2.M3 <- summary(lm(M3_SOCS$mean_vali ~ M3_SOCS$mean_pred))$r.squared
resids.M3 <- M3_SOCS$mean_pred - M3_SOCS$mean_vali
rmse.M3 <- sqrt(mean(resids.M3^2))
mbe.M3 <- mean(M3_SOCS$mean_pred - M3_SOCS$mean_vali)
plot(M3_SOCS$mean_vali, M3_SOCS$mean_pred, 
     main = "Measured vs. predicted SOC stocks (0-30 cm) of the validation set \n
     (Modeled with KSSL-local individual-based model for SOC density)", cex.main = 0.95,
     xlab = "Measured SOCS (g/m2)",
     ylab = "Predicted SOCS (g/m2)",
     xlim = c(0,20000), # Specify range for SOCS visualization
     ylim = c(0,20000))
abline(coef = c(0,1))
r2.M3 
rmse.M3
mbe.M3 
write.csv(M3_SOCS, paste0("./SOCS/M3/",ROI,"/",ROI,"_",scheme,"_",KSSL_num,"_",spike_num,".csv"))


###############################################################
# Calculate cluster-based results
M1_cluster <- merge(M1_SOCS,Cluster,by.x = "Core", by.y = "Num")
M1_cluster %>% group_by(Cluster) %>% 
  summarize(mean_vali = mean(mean_vali),mean_pred = mean(mean_pred))
mean(M1_cluster$mean_vali)
sd(M1_cluster$mean_vali)
mean(M1_cluster$mean_pred)
sd(M1_cluster$mean_pred)
M2_cluster <- merge(M2_SOCS,Cluster,by.x = "Core", by.y = "Num")
M2_cluster %>% group_by(Cluster) %>% 
  summarize(mean_vali = mean(mean_vali),mean_pred = mean(mean_pred))
mean(M2_cluster$mean_pred)
sd(M2_cluster$mean_pred)
M3_cluster <- merge(KSSL_local_socs,Cluster,by.x = "Core", by.y = "Num") 
M3_cluster %>% group_by(Cluster) %>% 
  summarize(mean_vali = mean(vali.i),mean_pred = mean(preds.i))
mean(M3_cluster$preds.i)
sd(M3_cluster$preds.i)