library(dplyr)
setwd("set-pathway")


###############################################################
# In the paper, validation samples were from Year 1. 
ROI = "specify-site-name"
spike_num = 36 # Choose number of spiking cores
KSSL_num = 500 # Choose number of calibration samples
# Import files for model calibration and validation
KSSL <- read.table('./KSSL.f.SOCS.csv', comment.char ="", quote = "\"", header = T, sep = ",")
Vali_SOCS <- read.table(paste0('./',ROI,'CSV-with-measured-SOCS-filename.csv')
                   , comment.char ="", quote = "\"", header = T, sep = ",")
#Get SoilGrid+ performance
summary(lm(Vali_SOCS$SOCS_sum ~ Vali_SOCS$SoilGrid))
resids <- Vali_SOCS$SOCS_sum - Vali_SOCS$SoilGrid
rmse <- sqrt(mean(resids^2))
mbe <- mean(Vali_SOCS$SoilGrid -Vali_SOCS$SOCS_sum)
rmse
mbe