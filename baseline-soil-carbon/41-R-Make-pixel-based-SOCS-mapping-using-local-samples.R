library(randomForest)
library(caret)
library(rgdal)
library(flexclust)
library(tidyverse)
library(raster)
setwd("set-pathway")


###############################################################
# Define site and import datasets
# Select site 
ROI = "specify-site-name"
# Import raster image of target site with extracted covariates
predict.raster <- paste0('./',ROI,'_cov_all_RAP.tif')
x <- new("GDALReadOnlyDataset", predict.raster)
width <- dim(x)[2]
height <- dim(x)[1]
imagedata <- data.frame(getRasterTable(x))
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
# Match with pixel ID
pred.ID <- tibble::rowid_to_column(imagedata, "ID")
imagedata <- cbind(pred.ID[,1], imagedata)
colnames(imagedata)[1] <- c("pred.ID")
# Import local samples
Spike <- read.table(paste0('./',ROI,'CSV-with-local-spiking-samples-joined-with-covariates.csv') 
                    , comment.char ="", quote = "\"", header = T, sep = ",")
length(unique(Spike$Point))
Spike <- Spike[,c("Point","X","Depth","SOCS_sum","BD3","SOC3","Clay3","Sand3","ppt","tmin","tmax","VPD",
                  "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE")]
colnames(Spike) <- c("sampleids","lay_id","Med_depth","socs_kssl","BD","SOC","Clay","Sand","ppt","tmin","tmax","VPD",
                     "LULC","Tree","GPP","EVI","AS","EL","SL","TWI","AFGC","BG","LTR","PFGC","SHR","TREE")


###############################################################
# Run the RF models with local samples to generate estimates
rf.p <- randomForest(socs_kssl ~ BD + Clay + Sand + SOC + 
                       ppt + tmin + tmax + VPD + Tree + GPP + EVI +
                       AS + EL + SL + TWI + AFGC + BG + LTR + PFGC + SHR + TREE,
                     data = Spike, importance = TRUE, ntree=500, mtry=3)
r2 <- summary(lm(rf.p$predicted ~ Spike$socs_kssl))$adj.r.squared
resids <- rf.p$predicted - Spike$socs_kssl
rmse <- sqrt(mean(resids^2))
mbe <- mean(resids)
imp <- importance(rf.p, type = 1)[,1]
preds <- predict(rf.p, newdata = imagedata)
dotchart(sort(imp), xlab="Variable importance (%)")
write.csv(imp.cov,paste0("./",ROI,"_",KSSL_num,"_pred_SOCS_local_imp.csv"))
preds <- data.frame(preds)
preds.i.reorg <- cbind(pred.ID[,1:3], preds[,1])
colnames(preds.i.reorg)[4] <- c("Pred_SOCS")
write.csv(preds.i.reorg,paste0("./",ROI,"_",KSSL_num,"_pred_SOCS_local_pixel.csv"))
# Create maps
res <- c(0.000269494585235856472, 0.000269494585235856472) # Define resolution 
x <- raster(xmn= XXX, xmx= XXX, 
            ymn= XXX, ymx= XXX, res=res,
            crs="+proj=longlat +datum=WGS84") # Define XXX with boundary parameters
SOCS_map <- rasterize(preds.i.reorg[, c('x', 'y')], x,
                      preds.i.reorg[, 'Pred_SOCS'], fun=mean, na.rm = TRUE)
plot(SOCS_map)
filepath <- paste0("./",ROI,"_",KSSL_num,"_predicted_SOCS_local")
writeRaster(SOCS_map, filepath,format="GTiff",overwrite=TRUE)