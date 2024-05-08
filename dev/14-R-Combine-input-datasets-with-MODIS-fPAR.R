library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)
library(stringr)
library(zoo) 
library(R.utils)

# Merge GPP datasets
setwd("C:/Users/username/path-to-the-folder/GPP")
data_cov <- read.table('final_extracted_GPP_covariates.csv', header=T, sep = ",")
data <- read.table('final_measured_GPP.csv', header=T, sep = ",")
data_cov$Date <- as.Date(data_cov$Date,  format = "%Y-%m-%d")
data$Date <- as.Date(data$Date,  format = "%Y-%m-%d")
data_cov$Site <- sapply(strsplit(data_cov$Site, ".csv"),`[`, 1)
data$Site <- sapply(strsplit(data$Site, ".csv"),`[`, 1)
data <- merge(data,data_cov,by = c("Site","Date"))
drop <- c("X.x","X.y")
data = data[,!(names(data) %in% drop)]
data <- data[which(!is.na(data$fPAR)),]
data <- data[which(!is.na(data$GPP_measured)),]

# Divide MODIS GPP dataset into different groups
listG <- c("LS1","Wjs","xAE","xCP","xDC","xKA","xKZ","xNG","xWD",
           "AR1","AR2","Aud","BMM","BRG","Bkg","Ctn","FPe",
           "IB2","KM2","KM3","KM4","RFW","Ro4","Var") # Add a list of sites for the specific vegetation type
listH <- c("Snf","Tx1","Tx2")
listS <- c("LS2","xJR","xMB","xNQ","xYE","Jo1","Rls",
           "Rms","Rwe","Rwf","Rws","SRM","Seg","Ses","Wkg")
listT <- c("Mpj","Ton","xCL","xSJ","Fwf")
data.G <- data[data$Site %in% listG,]
summary(as.factor(data.G$Site))
data.H <- data[data$Site %in% listH,]
summary(as.factor(data.H$Site))
data.S <- data[data$Site %in% listS,]
summary(as.factor(data.S$Site))
data.T <- data[data$Site %in% listT,]
summary(as.factor(data.T$Site))
write.csv(data.G, "C:/Users/username/path-to-the-folder/GPP/data.G.csv")
write.csv(data.H, "C:/Users/username/path-to-the-folder/GPP/data.H.csv")
write.csv(data.S, "C:/Users/username/path-to-the-folder/GPP/data.S.csv")
write.csv(data.T, "C:/Users/username/path-to-the-folder/GPP/data.T.csv")

# Merge GPP datasets
setwd("C:/Users/username/path-to-the-folder/NEE")
data_cov <- read.table('final_extracted_NEE_covariates.csv', header=T, sep = ",")
data_cov$Site <- sapply(strsplit(data_cov$Site, ".csv"),`[`, 1)
data <- read.table('final_measured_NEE.csv', header=T, sep = ",")
data$Site <- sapply(strsplit(data$Site, ".csv"),`[`, 1)
data$Date <- as.Date(data$Date,  format = "%Y-%m-%d")
data_cov$Date <- as.Date(data_cov$date,  format = "%Y-%m-%d")
data_cov <- data_cov[which(data_cov$Date >= '2003-01-01'),]
data_cov$Site <- str_replace(data_cov$Site, "_2022", "") # If data is extracted separately for a year due to covariate dataset being unavailable
summary(as.factor(data_cov$Site))
data <- merge(data,data_cov,by = c("Site","Date"), all = TRUE) 
drop <- c("X.x","X.y","date")
data = data[,!(names(data) %in% drop)]
data$site_num = as.numeric(as.factor(with(data,paste(Site))))
site <- unique(data[c("site_num","Site")])
data$DOY <- yday(data$Date)
data$day <- difftime(data$Date,"2003-01-01",units = "days")
data$day <- round(as.numeric(sub(" .*", "", data$day)))

# Optional gap-filling of dataset if not done so previously; Change time-step of the input dataset if needed
data_f = data.frame()
for (i in 1:length(unique(data$site_num))){
  data_i <- data[which(data$site_num == i),]
  x<- data_i$fPAR
  for (j in 1:7060){# Total number of days for the time period 
    lagx = dplyr::lag(x)
    leadx = dplyr::lead(x)
    x = ifelse(is.na(x),(leadx+lagx)/2, x)
    x = ifelse(is.na(x), leadx, x)
    x = ifelse(is.na(x), lagx, x)
    if (sum(is.na(x))==0) {break}
  }
  data_i$fPAR = x
  data_f = rbind(data_f,data_i)
}
data_f = data_f[which((data_f$DOY %% 4 == 1)&(data_f$DOY < 365)),]
data_f$Year <- as.numeric(substr(data_f$Date,1,4))
data_f$POY = round((data_f$DOY+2)/4)
date <- expand_grid(Year= 2003:2022, POY = 1:91, value = "NA")
date <- date[which((date$Year <= 2021)|(date$POY <= 31)),]
date$Date <- strptime(paste(date$Year, date$POY*4-3), format="%Y %j")
colnames(date)[1:3] <- c("Year","POY","Value") 
fill <- function(x) {
  for (i in 1:1760){# total number of 4-days for the time period #1729 if 2022 is not considered
    lagx = dplyr::lag(x)
    leadx = dplyr::lead(x)
    x = ifelse(is.na(x),(leadx+lagx)/2, x)
    x = ifelse(is.na(x), leadx, x)
    x = ifelse(is.na(x), lagx, x)
    if (sum(is.na(x))==0) {break}
  }
  return (x)
}
data_ff = data.frame()
for (i in 1:length(unique(data_f$site_num))){
  data_i <- data_f[which(data_f$site_num == i),]
  data_i <- merge(data_i,date,by=c("Year","POY"),all=TRUE)
  data_i$SWC1 = fill(data_i$SWC1)
  data_i$SWC2 = fill(data_i$SWC2)
  data_i$Clay = fill(data_i$Clay)
  data_i$SW_IN = fill(data_i$SW_IN)
  data_i$SW_IN_NLDAS = fill(data_i$SW_IN_NLDAS)
  data_i$TA = fill(data_i$TA)
  data_i$TA_min = fill(data_i$TA_min)
  data_i$TS = fill(data_i$TS)
  data_i$VPD = fill(data_i$VPD)
  data_i$fPAR = fill(data_i$fPAR)
  data_i$ppt = fill(data_i$ppt)
  data_i$site_num = i
  data_i$Site = site[which(site$site_num == i),]$Site
  data_ff = rbind(data_ff,data_i)
}
data_ff[which(is.na(data_ff$DOY)),]$Date.y = strptime(paste(data_ff[which(is.na(data_ff$DOY)),]$Year, 
                                                            data_ff[which(is.na(data_ff$DOY)),]$POY*4-3), format="%Y %j")
data_ff[which(is.na(data_ff$day)),]$day = difftime(data_ff[which(is.na(data_ff$day)),]$Date.y,"2003-01-01",units = "days")
data_ff$day <- round(as.numeric(sub(" .*", "", data_ff$day)))
drop <- c("DOY","Value","site_num", "Date.x")
data_ff = data_ff[,!(names(data_ff) %in% drop)]
colnames(data_ff)[dim(data_ff)[2]] <- "Date"
data_ff <- data_ff[which(data_ff$day >=0),]

# Combine all sites and then divide them by PFTs
summary(as.factor(data_ff$Site))
data0 <- read.table('final_measured_NEE.csv', header=T, sep = ",")
summary(as.factor(data0$Site))
listG <- c("LS1","Wjs","SCg","SdH","xAE","xCP","xDC","xKA","xKZ","xNG","xWD",
           "AR1","AR2","ARb","ARc","Aud","BMM","BRG","Bkg","Ctn","Dia","FPe",
           "Hn2","IB2","KFS","KLS","KM2","KM3","KM4","Kon","RFW","Ro4","Shd",
           "Var")
listH <- c("A32","Snd","Tx1","Tx2")
listS <- c("LS2","xJR","xMB","xNQ","xYE","Jo1","Rls","Hn1",
           "Rms","Rwe","Rwf","Rws","SRM","Seg","Ses","Wkg","Wdn","Cop")
listT <- c("Mpj","Ton","xCL","xSJ","Fwf","CZ1")
data.G <- data_ff[data_ff$Site %in% listG,]
summary(as.factor(data.G$Site))
data.H <- data_ff[data_ff$Site %in% listH,]
summary(as.factor(data.H$Site))
data.S <- data_ff[data_ff$Site %in% listS,]
summary(as.factor(data.S$Site))
data.T <- data_ff[data_ff$Site %in% listT,]
summary(as.factor(data.T$Site))
write.csv(data.G, "C:/Users/username/path-to-the-folder/NEE/data.G.csv")
write.csv(data.H, "C:/Users/username/path-to-the-folder/NEE/data.H.csv")
write.csv(data.S, "C:/Users/username/path-to-the-folder/NEE/data.S.csv")
write.csv(data.T, "C:/Users/username/path-to-the-folder/NEE/data.T.csv")