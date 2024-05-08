library(xts)
library(dygraphs)
library(zoo) 
library(R.utils)
setwd("C:/Users/username/path-to-the-flux-tower-folder")
data <- read.table('sitename/sitename_summary.csv', header=T, sep = ",") # Read file with daily summaries
drop <- c("X")
data$Date <- as.Date(data$Date, format = "%Y-%m-%d")

# Check data trends/spikes with meteorological data if needed
mt_list_a = c("TA","RH")
mt_avg_list <- paste(mt_list_a, "avg", sep = "_")
comb_all_mt_a <- data[,c('Date',mt_avg_list)]
comb_all_mt_a_reorg <- xts(x = comb_all_mt_a[,-1], order.by = as.POSIXct(comb_all_mt_a$Date))
plot <- dygraph(comb_all_mt_a_reorg, main = "Visualization for meterology") %>%
  dySeries("RH_avg", label = "Relative humidity", axis = 'y') %>%
  dySeries("TA_avg", label = "air temp", axis = 'y2') %>%
  dyHighlight(highlightSeriesOpts = list(strokeWidth = 1)) %>%
  dyAxis("y", label = "RH (%)", valueRange = c(20, 110)) %>%
  dyAxis("y2", label = "Temperature (C)", valueRange = c(-20, 40)) %>%
  dyRangeSelector() %>%
  dyLegend(width = 650) %>%
  dyOptions(colors = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(10))
print(plot)

#set parameters and range for data filtering
window_size <- 7 # Window size can be changed, here we used 15 days
start_m <- 4 # Start month can be customized 
end_m <- 10 # End month can be customized 
#Date_min <- "2009-01-01"
#data <- data[which(data$Date>=Date_min),] 
data_reorg <- xts(x = data[,c("GPP_PI_F_avg")], order.by = data$Date)
data_GPP <- data[c("Date","GPP_PI_F_avg")]
data_GPP$Date <- as.Date(data_GPP$Date, format = "%Y-%m-%d")
data_GPP$Month <- as.numeric(substr(data_GPP$Date, 6,7))

# Fill in winter month with 0 and use linear interpolation for the growing season
data_GPP[which((is.na(data_GPP$GPP_PI_F_avg))&((data_GPP$Month <=start_m)|(data_GPP$Month >=end_m))),]$GPP_PI_F_avg <- 0
data_GPP$GPP_PI_F_avg <- na.approx(data_GPP$GPP_PI_F_avg)
moving_median <- rollmedian(data_GPP$GPP_PI_F_avg, k = 2*window_size+1) 
moving_median <- insert(moving_median, ats=1:window_size, values=moving_median[1])
moving_median <- insert(moving_median, ats=length(moving_median)+1, values=rep(moving_median[length(moving_median)],times = window_size))
data_GPP$ref <- moving_median
data_GPP$diff <- abs(data_GPP$GPP_PI_F_avg-data_GPP$ref)
threashold <- 2*sd(data[which(!is.na(data$GPP_PI_F_avg)),]$GPP_PI_F_avg) # Customize threshold
date_list <- data_GPP[which((data_GPP$diff > threashold)&((data_GPP$Month <=start_m)|(data_GPP$Month >=end_m))),]$Date
data[which(data$Date %in% date_list),]$GPP_PI_F_avg <- NA
data_reorg2 <- xts(x = data[,c("GPP_PI_F_avg")], order.by = data$Date)
plot2 <- dygraph(data_reorg2) %>%
  dySeries("V1", label = "GPP") %>%
  dyAxis("y", label = "µmol CO<sub>2</sub> m<sup>-2</sup> s<sup>-1</sup>", valueRange = c(0, 15)) %>%
  print(plot2)
data <- within (data, NEE_PI_F_avg[is.na(GPP_PI_F_avg)] <- NA)
data <- within (data, RECO_PI_F_avg[is.na(GPP_PI_F_avg)] <- NA)

# Same process for RECO outlier removal
data_RECO <- data[c("Date","RECO_PI_F_avg")]
data_RECO$Date <- as.Date(data_RECO$Date, format = "%Y-%m-%d")
data_RECO$Month <- as.numeric(substr(data_RECO$Date, 6,7))
data_RECO[which((is.na(data_RECO$RECO_PI_F_avg))&((data_RECO$Month <=start_m)|(data_RECO$Month >=end_m))),]$RECO_PI_F_avg <- 0
data_RECO$RECO_PI_F_avg <- na.approx(data_RECO$RECO_PI_F_avg)
moving_median2 <- rollmedian(data_RECO$RECO_PI_F_avg, k = 2*window_size+1) 
moving_median2 <- insert(moving_median2, ats=1:window_size, values=moving_median2[1])
moving_median2 <- insert(moving_median2, ats=length(moving_median2)+1, values=rep(moving_median2[length(moving_median2)],times = window_size))
data_RECO$ref <- moving_median2
data_RECO$diff <- abs(data_RECO$RECO_PI_F_avg-data_RECO$ref)
threashold2 <- 2*sd(data[which(!is.na(data$RECO_PI_F_avg)),]$RECO_PI_F_avg)
date_list2 <- data_RECO[which((data_RECO$diff > threashold2)&((data_RECO$Month <=start_m)|(data_RECO$Month >=end_m))),]$Date
data[which(data$Date %in% date_list2),]$RECO_PI_F_avg <- NA
data <- within (data, NEE_PI_F_avg[is.na(RECO_PI_F_avg)] <- NA)
data <- within (data, GPP_PI_F_avg[is.na(RECO_PI_F_avg)] <- NA)

# Take out filtered data and visualize
flux_list <- c("GPP","NEE","RECO")
flux_avg_list <- paste(flux_list, "PI_F_avg", sep = "_")
data_flux <- data[,c('Date',flux_avg_list)]
data_flux_reorg <- xts(x = data_flux[,-1], order.by = data_flux$Date)
flux_plot <- dygraph(data_flux_reorg, main = "Visualization for flux measures") %>%
  dyHighlight(highlightSeriesOpts = list(strokeWidth = 1)) %>%
  dySeries("GPP_PI_F_avg", label = "GPP") %>%
  dySeries("NEE_PI_F_avg", label = "NEE") %>%
  dySeries("RECO_PI_F_avg", label = "RECO") %>%
  dyAxis("y", label = "µmol CO<sub>2</sub> m<sup>-2</sup> s<sup>-1</sup>", valueRange = c(-20, 30)) %>%
  dyRangeSelector() %>%
  dyLegend(width = 650) %>%
  dyOptions(colors = RColorBrewer::brewer.pal(3, "Set1"))
print(flux_plot)

# Export quality controlled datasets
data <- data[which(!is.na(data$GPP_PI_F_avg)),]
data = data[,!(names(data) %in% drop)]
write.csv(data,"username/path-to-the-flux-tower-folder/Final_GPP/sitename.csv")

