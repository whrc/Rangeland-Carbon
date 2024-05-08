library(xts)
library(dygraphs)
library(zoo) 
library(R.utils)
setwd("C:/Users/username/path-to-the-flux-tower-folder")
data <- read.table('sitename/sitename_summary.csv', header=T, sep = ",") # Read file with daily summaries
drop <- c("X")
data$Date <- as.Date(data$Date, format = "%Y-%m-%d")

#set parameters and range for data filtering
window_size <- 7 # Window size can be changed, here we used 15 days
start_m <- 4 # Customize range for data filling 
end_m <- 10 
s_m = 4 # Customize range for cutting spikes
e_m = 10

# QC NEE data
data_NEE <- data[c("Date","NEE_PI_F_avg")]
data_NEE$Date <- as.Date(data_NEE$Date, format = "%Y-%m-%d")
data_NEE$Month <- as.numeric(substr(data_NEE$Date, 6,7))
# Fill in winter month with 0 and use linear interpolation for the growing season
data_NEE[which((is.na(data_NEE$NEE_PI_F_avg))&((data_NEE$Month <=start_m)|(data_NEE$Month >=end_m))),]$NEE_PI_F_avg <- 0
data_NEE$NEE_PI_F_avg <- na.approx(data_NEE$NEE_PI_F_avg)
moving_median <- rollmedian(data_NEE$NEE_PI_F_avg, k = 2*window_size+1) 
moving_median <- insert(moving_median, ats=1:window_size, values=moving_median[1])
moving_median <- insert(moving_median, ats=length(moving_median)+1, values=rep(moving_median[length(moving_median)],times = window_size))
data_NEE$ref <- moving_median
data_NEE$diff <- abs(data_NEE$NEE_PI_F_avg-data_NEE$ref)
threashold <- 2*sd(data[which(!is.na(data$NEE_PI_F_avg)),]$NEE_PI_F_avg) # Customize threshold
date_list <- data_NEE[which((data_NEE$diff > threashold)&((data_NEE$Month <=s_m)|(data_NEE$Month >=e_m))),]$Date
data[which(data$Date %in% date_list),]$NEE_PI_F_avg <- NA
data_reorg2 <- xts(x = data[,c("NEE_PI_F_avg")], order.by = data$Date)
plot2 <- dygraph(data_reorg2) %>%
  dySeries("V1", label = "NEE") %>%
  dyAxis("y", label = "µmol CO<sub>2</sub> m<sup>-2</sup> s<sup>-1</sup>", valueRange = c(-10, 20)) %>%
  print(plot2)

# Export QC results
data <- data[which(!is.na(data$NEE_PI_F_avg)),]
drop <- c("X")
data = data[,!(names(data) %in% drop)]
write.csv(data,"C:/Users/username/path-to-the-flux-tower-folder/Final_NEE/sitename.csv")