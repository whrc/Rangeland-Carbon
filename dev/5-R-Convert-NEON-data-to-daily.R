library(tidyverse)
site = "sitename"
abbv = "sitename-starts-with-x"
num_a = 3 # Change parameter based on file name
num_b = 5 # Change parameter based on file name
setwd(paste0("C:/Users/username/path-to-NEON-data-summary/",site))
data <- read.table(paste0('AMF_US-',abbv,'_BASE_HH_',num_a,'-',num_b,'.csv'), skip = 2, header=T, sep = ",")
data1 <- read.table('flux_summary.csv', header=T, sep = ",")[,-1]
data2 <- read.table('met_summary.csv', header=T, sep = ",")[,-1]
data1$Date <- as.Date(data1$Date, format = "%Y-%m-%d")
data2$Date <- as.Date(data2$Date, format = "%Y-%m-%d")
data2$PA <- data2$PA/1000 # Change unit if needed

# Process Ameriflux data for soil temp, mois, and met data
data[,(dim(data)[2]+1)] <- round(data[,1]/1e8)
data[,(dim(data)[2])+1] <- round((data[,1] - data[,(dim(data)[2])]*1e8)/1e6)
data[,(dim(data)[2])+1] <- round((data[,1] - data[,(dim(data)[2]-1)]*1e8 - data[,(dim(data)[2])]*1e6)/1e4)
data[,(dim(data)[2])+1] <- round((data[,1] - data[,(dim(data)[2]-2)]*1e8 - data[,(dim(data)[2]-1)]*1e6 - data[,(dim(data)[2])]*1e4)/1e2)
data[,(dim(data)[2])+1] <- (data[,1] - data[,(dim(data)[2]-3)]*1e8 - data[,(dim(data)[2]-2)]*1e6 - data[,(dim(data)[2]-1)]*1e4 - data[,(dim(data)[2])]*1e2)
colnames(data)[(dim(data)[2]-4):dim(data)[2]] <- c("Year","Month","Day","Hour","Minute")
data[,(dim(data)[2]+1)] <- as.Date(with(data, paste(Year, Month, Day, sep="-")), "%Y-%m-%d")
colnames(data)[(dim(data)[2])] <- c("Date")
names(data)
colnames(data)[colnames(data) == "VPD_PI"] <- "VPD"

# Treat flux outliers
data1$RECO <- data1$GPP + data1$NEE
colnames(data1)[colnames(data1) == "GPP"] <- "GPP_PI_F"
colnames(data1)[colnames(data1) == "NEE"] <- "NEE_PI_F"
colnames(data1)[colnames(data1) == "RECO"] <- "RECO_PI_F"
data1.GPP.check <- data1[which(data1$GPP_PI_F != -9999),]
data1.GPP.neg <- data1.GPP.check[which(data1.GPP.check$GPP_PI_F < 0),]
summary(data1.GPP.neg$GPP_PI_F)
data1.GPP.pos <- data1[which(data1$GPP_PI_F > 0),]
threshold.GPP <- -sd(data1.GPP.pos$GPP_PI_F)
data1.RECO.check <- data1[which(data1$RECO_PI_F != -9999),]
data1.RECO.neg <- data1.RECO.check[which(data1.RECO.check$RECO_PI_F < 0),]
summary(data1.RECO.neg$RECO_PI_F)
data1.RECO.pos <- data1[which(data1$RECO_PI_F > 0),]
threshold.RECO <- -sd(data1.RECO.pos$RECO_PI_F)

# Replace missing data and outliers with NA
data1 <- data1 %>%
  mutate(GPP_PI_F = replace(GPP_PI_F, ((GPP_PI_F == -9999)|(GPP_PI_F <= threshold.GPP)), NA))
data1 <- data1 %>%
  mutate(RECO_PI_F = replace(RECO_PI_F, ((RECO_PI_F == -9999)|(RECO_PI_F <= threshold.RECO)), NA))
data <- data %>%
  mutate_at(vars(contains("SWC")), ~na_if(.,-9999))
data <- data %>%
  mutate_at(vars(contains("SWC")), ~replace(., .< 0, NA))
mois_list <- ls(name = data, pattern = "SWC")
mois_list
data <- data %>%
  mutate_at(vars(contains("TS")), ~na_if(.,-9999))
temp_list <- ls(name = data, pattern = "TS")
temp_list
mt_list <- ls(name = data, pattern = "(WD)|(VPD)|(ALB)")
mt_list
data <- data %>%
  mutate_at(vars(any_of(mt_list)), ~na_if(.,-9999))

# Summarize availability, average, and std for ancillary data
list = c(mois_list,temp_list, mt_list)
for (i in 3:(dim(data)[2]-6)){
  if (colnames(data)[i] %in% list){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_availability")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := 1-sum(is.na(!!sym(prop))/length(!!sym(prop)))))
  }
}
for (i in 3:(dim(data)[2]-6)){
  if (colnames(data)[i] %in% list){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_avg")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := mean(!!sym(prop), na.rm=TRUE)))
  }
}
for (i in 3:(dim(data)[2]-6)){
  if (colnames(data)[i] %in% list){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_std")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := sd(!!sym(prop), na.rm=TRUE)))
  }
}

# Summarize availability, average, and std for flux data
data <- data1
flux_list <- c("GPP_PI_F","RECO_PI_F","NEE_PI_F")
for (i in 2:(dim(data)[2])){
  if (colnames(data)[i] %in% flux_list){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_availability")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := 1-sum(is.na(!!sym(prop))/length(!!sym(prop)))))
  }
}
for (i in 2:(dim(data)[2])){
  if (colnames(data)[i] %in% flux_list){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_avg")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := mean(!!sym(prop), na.rm=TRUE)))
  }
}
for (i in 2:(dim(data)[2])){
  if (colnames(data)[i] %in% flux_list){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_std")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := sd(!!sym(prop), na.rm=TRUE)))
  }
}

# Summarize availability, average, and std for met data
data <- data2
mt_list_a <- c("RH","WS","PA","TA")
mt_list_s <- c("ppt","SW_IN")
met_list <- list(mt_list_a,mt_list_s)
for (i in 2:(dim(data)[2])){
  if (colnames(data)[i] %in% met_list){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_availability")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := 1-sum(is.na(!!sym(prop))/length(!!sym(prop)))))
  }
}
for (i in 2:(dim(data)[2])){
  if (colnames(data)[i] %in% mt_list_a){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_avg")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := mean(!!sym(prop), na.rm=TRUE)))
  }
}
for (i in 3:(dim(data)[2]-6)){
  if (colnames(data)[i] %in% mt_list_s){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_sum")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := sum(!!sym(prop), na.rm=TRUE)))
  }
}
for (i in 2:(dim(data)[2])){
  if (colnames(data)[i] %in% met_list){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_std")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := sd(!!sym(prop), na.rm=TRUE)))
  }
}

# Merge files 
comb <- reduce(mget(ls(pattern="availability")), full_join)
comb <- comb %>%
  mutate_at(vars(contains("availability")),
            function(x) (x*100))
comb_avg <- reduce(mget(ls(pattern="_avg")), full_join)
comb_sum <- reduce(mget(ls(pattern="_sum")), full_join)
comb_std <- reduce(mget(ls(pattern="_std")), full_join)
comb$Date <- as.Date(comb$Date, format = "%Y-%m-%d")
comb_avg$Date <- as.Date(comb_avg$Date, format = "%Y-%m-%d")
comb_sum$Date <- as.Date(comb_sum$Date, format = "%Y-%m-%d")
comb_std$Date <- as.Date(comb_std$Date, format = "%Y-%m-%d")
comb_all <- merge(merge(merge(comb, comb_avg, by = "Date", all = TRUE),
                        comb_sum, by = "Date", all = TRUE),
                  comb_std, by = "Date", all = TRUE)
comb_all <- cbind(comb_all[,1],
                  comb_all[,-1] %>%
                    mutate_all(~ifelse(is.nan(.), NA, .)))
colnames(comb_all)[1] <-"Date"
rm(list = ls(pattern="_availability"))
rm(list = ls(pattern="_avg"))
rm(list = ls(pattern="_std"))
rm(list = ls(pattern="_sum"))
names(comb) <- gsub(x = names(comb), pattern = " _availability", replacement = "_availability")

# Exclude outliers or records with low availability
names(comb_all) <- gsub(x = names(comb_all), pattern = " _availability", replacement = "_availability")
comb_all <- comb_all %>% 
  mutate (Outlier_GPP_RECO = if_else (((GPP_PI_F_avg < 0)&(!is.na(GPP_PI_F_avg)))|((RECO_PI_F_avg < 0)&(!is.na(RECO_PI_F_avg))), "Y","N")) 
comb_all <- comb_all %>%
  mutate_at(vars(contains("GPP")), ~replace(., .< 0, NA))
comb_all <- comb_all %>%
  mutate_at(vars(contains("RECO")), ~replace(., .< 0, NA)) 
comb_all <- comb_all %>% 
  mutate (Low_availability = if_else ((GPP_PI_F_availability < 40)|(RECO_PI_F_availability < 70), "Y","N")) 
comb_all <- within (comb_all, GPP_PI_F_avg[GPP_PI_F_availability < 40] <- NA)
comb_all <- within (comb_all, RECO_PI_F_avg[RECO_PI_F_availability < 70] <- NA)
comb_all <- comb_all %>% select(-contains("availability"))
comb_all$NEE_PI_F_avg <- -comb_all$NEE_PI_F_avg
comb_all$RECO_PI_F_avg <- -comb_all$RECO_PI_F_avg

# Visualization
library("dygraphs")
library("xts")
library("htmlwidgets")
flux_list <- c("GPP","NEE","RECO")
flux_avg_list <- paste(flux_list, "PI_F_avg", sep = "_")
data_flux <- comb_all[,c('Date',flux_avg_list)]
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

#Export files
write.csv(comb,paste0(site,"_availability.csv"))
write.csv(comb_all,paste0(site,"_summary.csv"))

