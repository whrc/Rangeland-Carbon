library(ncdf4) 
library(lubridate)

# Read nc data file
setwd("C:/Users/username/path-to-the-file")
data <- nc_open('filename.nc')

# Extract flux variables
names(data$var) # Get variables of interest, below are just some examples
lon <- ncvar_get(data, "longitude")
lat <- ncvar_get(data, "latitude")
dom <- ncvar_get(data, "dom")
month <- ncvar_get(data, "month")
doy <- ncvar_get(data, "start_doy")
year <- ncvar_get(data, "year")
NEE <- ncvar_get(data, "NEE")
GPP <- ncvar_get(data, "GPP")
RECO <- ncvar_get(data, "Reco")
NEE_uncertainty <- ncvar_get(data, "NEE_Un_tot")
GPP_uncertainty <- ncvar_get(data, "GPP_Un_tot")
RECO_uncertainty <- ncvar_get(data, "RECO_Un_tot")
flux <- data.frame(lon,lat,year,month,doy,dom,NEE,GPP,RECO,NEE_uncertainty,GPP_uncertainty,RECO_uncertainty)
flux$Date <- with(flux,ymd(paste(year,month,dom,sep=" ")))
data <- flux[,c("Date","NEE","GPP","RECO")]
data <- data %>%
  mutate_at(vars(contains("GPP")), ~replace(., .< 0, NA))
data <- data %>%
  mutate_at(vars(contains("RECO")), ~replace(., .< 0, NA)) 
data$NEE <- -data$NEE
data$RECO <- -data$RECO
colnames(data)[colnames(data) == "NEE"] <- "NEE_PI_F_avg"
colnames(data)[colnames(data) == "RECO"] <- "RECO_PI_F_avg"
colnames(data)[colnames(data) == "GPP"] <- "GPP_PI_F_avg"
data0 <- data

# Extract ancillary variables
data <- nc_open('filename-for-ancillary.nc')
names(data$var)
# Please Check if labels are correct for the variables
lon <- ncvar_get(data, "lon")
lat <- ncvar_get(data, "lat")
time <- ncvar_get(data, "timestp")
TA <- ncvar_get(data, "Tair") # Air temperature, K, convert to C to be consistent with other Ameriflux sites
SH <- ncvar_get(data, "Qair") # Specific humidity, kg/kg 
WS <- ncvar_get(data, "Wind") # Wind speed, m/s
ppt <- ncvar_get(data, "Rainf") # precipitation (kg/m2/s), convert to mm to be consistent with Ameriflux data
PA <- ncvar_get(data, "Psurf") # Air pressure, Pa
SW <- ncvar_get(data, "SWdown") # Shortwave radiation (W/m2)
LW <- ncvar_get(data, "LWdown") # longwave radiation (W/m2)
ancillary <- data.frame(time,SW,ppt,PA,TA,SH,WS)
ancillary$RH <- 0.263*ancillary$SH*ancillary$PA/exp((17.67*(ancillary$TA-273.15))/(ancillary$TA-29.65))
ancillary$TA <- ancillary$TA - 273.15
ancillary$ppt <- ancillary$ppt*24*60*60
ancillary$PA <- ancillary$PA/1000
# DOY and Year not provided in this example so they need to be assigned manually
ancillary$day <- ceiling((ancillary$time+1)/48)
a1<- ancillary[which(ancillary$day <= 366),]
a1$Year <- "2004"
a2<- ancillary[which((ancillary$day > 366)&(ancillary$day <= 731)),]
a2$Year <- "2005"
a2$day <- a2$day - 366
a3<- ancillary[which((ancillary$day > 731)&(ancillary$day <= 1096)),]
a3$Year <- "2006"
a3$day <- a3$day - 731
a4<- ancillary[which(ancillary$day > 1096),]
a4$Year <- "2007"
a4$day <- a4$day - 1096
ancillary <- rbind(a1,a2,a3,a4)
ancillary$Date <- strptime(paste(ancillary$Year, ancillary$day), format="%Y %j")
rm("a1","a2","a3","a4","data")
data <- ancillary[,c("Date","SW","ppt","PA","TA","RH","WS")]
colnames(data)[2:3] <- c("SW_IN","P")
rm("ancillary")

# Calculate averages and standard deviations
mt_list_a <- c("TA","PA","RH","WS")
mt_list_s <- c("P","SW_IN")
all_list = c(mt_list_a,mt_list_s)
for (i in 2:(dim(data)[2])){
  if (colnames(data)[i] %in% mt_list_a){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_avg")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := mean(!!sym(prop), na.rm=TRUE)))
  }
}
for (i in 2:(dim(data)[2])){
  if (colnames(data)[i] %in% mt_list_s){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_sum")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := sum(!!sym(prop), na.rm=TRUE)))
  }
}
for (i in 2:(dim(data)[2])){
  if (colnames(data)[i] %in% all_list){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_std")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := sd(!!sym(prop), na.rm=TRUE)))
  }
}

# The rest of the codes can be adapted from step 2 after combining these files with original Ameriflux files

# An alternative example of site file in the matlab format
library(R.matlab)
library(dplyr)
data <- readMat('sitename/filename.mat')
data <- data$DATA
time <- data[,1] 
WD <- data[,3]
WS <- data[,4]
TA <- data[,5]
RH <- data[,9]
SR <- data[,13]
NDVI <- data[,21]
SWC1 <- data[,23] 
SWC2 <- data[,24] 
SWC3 <- data[,25] 
SWC4 <- data[,26] 
TS1 <- data[,27] 
TS2 <- data[,28] 
TS3 <- data[,29] 
TS4 <- data[,30] 
TS5 <- data[,31] 
NEE <- data[,19] # umol CO2/m2/s
flux <- data.frame(time,WD,WS,TA,RH,SR,NDVI,SWC1,SWC2,SWC3,SWC4,TS1,TS2,TS3,TS4,TS5,NEE)
flux$time <- as.integer(flux$time)
flux <- aggregate(flux[,2:dim(flux)[2]],list(flux$time),mean)
flux <- transform(flux,SWC = rowMeans(flux[,8:11], na.rm = TRUE))
flux$SWC <- flux$SWC*100 # Unit conversion for soil moisture data if needed
flux$Date <- as.Date("2006-01-01") + flux$Group.1 # Customize start date for the dataset
flux <- flux[,c("Date","WD","WS","TA","RH","SWC","TS1","TS2","TS3","TS4","TS5","NEE")]
flux <- flux %>% mutate_all(~ifelse(is.nan(.), NA, .))
colnames(flux)[2:12] <- c("WD_avg","WS_avg","TA_avg","RH_avg","SWC_avg","TS1_avg","TS2_avg","TS3_avg","TS4_avg","TS5_avg","NEE_PI_F_avg")
flux$Date <- as.Date(flux$Date, format = "%Y-%m-%d")
flux$NEE_PI_F_avg <- -flux$NEE_PI_F_avg