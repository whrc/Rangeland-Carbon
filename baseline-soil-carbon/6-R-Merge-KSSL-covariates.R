setwd("set-pathway")


###############################################################
# Import and merge KSSL datasets
KSSL.cov <- read.table('KSSL/KSSL_cov.txt', header = TRUE)
drop <- c("BD","SOC","clay","sand")
KSSL.cov = KSSL.cov[,!(names(KSSL.cov) %in% drop)]
KSSL.cov2022 <- read.table('KSSL/KSSL_cov_2022.txt', header = TRUE, sep=",")
KSSL.cov2022 <- KSSL.cov2022[,c("sampleids","BD","SOC","Clay","Sand")]
KSSL_cov <- merge(KSSL.cov,KSSL.cov2022,by="sampleids")
KSSL.cov2 <- read.table('KSSL/KSSL_cov_V2.csv', comment.char="", quote = "\"", header=T, sep=",")
KSSL.cov2 <- unique(KSSL.cov2)
KSSL.RAP <- read.table('KSSL/KSSL_RAP.txt', header = TRUE)
KSSL.all <- merge(KSSL_cov,KSSL.cov2,by=c("lay_id","sampleids","Med_depth"))
KSSL.all <- na.omit(KSSL.all)


###############################################################
# Process and refine covariate datasets
KSSL.all$ppt <- KSSL.all$ppt20/20
KSSL.all$tmin <- KSSL.all$tmin20/20
KSSL.all$tmax <- KSSL.all$tmax20/20
KSSL.all$VPD <- KSSL.all$VPD20/20
KSSL.all$EVI <- KSSL.all$EVI20
KSSL.all$GPP <- KSSL.all$GPP20
KSSL.all$Tree <- KSSL.all$Tree20
KSSL.all$AS <- sin(pi*(KSSL.all$AS)/360) 
drop <- c("oc_usda.y","ppt20","tmin20","tmax20","T","VPD20","Tree20","GPP20","EVI20",
          "coords","FIPS","Low_depth","Upper_depth")
KSSL.refined = KSSL.all[,!(names(KSSL.all) %in% drop)]
colnames(KSSL.refined)[XX] <- "oc_usda" # Specify column number XX that corresponds to SOC measures
KSSL.refined <- KSSL.refined[,c(XX)] # Specify column numbers corresponding to the needed covariates
KSSL.refined <- merge(KSSL.refined,KSSL.RAP,by="sampleids", all= TRUE) # Include all records with or without RAP for now
KSSL.refined <- KSSL.refined[which(!is.na(KSSL.refined$ppt)),]
KSSL <- KSSL.refined[which(KSSL.refined$oc_usda >0),]
write.csv(KSSL,"KSSL/KSSL.f.2022.csv") 