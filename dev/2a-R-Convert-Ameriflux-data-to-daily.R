library(tidyverse)
library(lubridate)
setwd("C:/Users/username/path-to-the-flux-tower-dataset")

# Data from this site is supplied by PI so is already quality-controlled and gap-filled
data <- read.table('sitename/filename.Daily.csv', header=T, sep = ",") # PI supplied file already at daily time step, if applicable
data$date <- as.Date(data$date, format = "%Y-%m-%d")
data$e_f <- data$e_f*10 # Convert pressure from kPa to hPa, if needed
data$GPP_g <- data$GPP_g/12*10^6/3600/24 # Convert from gC/day to umol CO2/s, if needed
data$Reco_g <- data$Reco_g/12*10^6/3600/24
data$NEE_PI_F <- data$GPP_g - data$Reco_g
data_reorg <- data[,c("date","GPP_g","NEE_PI_F","Reco_g","Ta_f","e_f","P_f","SW_In_f","SWC_f")]

# We need to combine PI-supplied data with Ameriflux deposit for several missing covariates
data <- read.table('sitename/AMF_US-sitename_BASE_HH_parameters.csv', skip = 2, header=T, sep = ",") # Import standard Ameriflux files
data[,(dim(data)[2]+1)] <- round(data[,1]/1e8)
data[,(dim(data)[2])+1] <- round((data[,1] - data[,(dim(data)[2])]*1e8)/1e6)
data[,(dim(data)[2])+1] <- round((data[,1] - data[,(dim(data)[2]-1)]*1e8 - data[,(dim(data)[2])]*1e6)/1e4)
data[,(dim(data)[2])+1] <- round((data[,1] - data[,(dim(data)[2]-2)]*1e8 - data[,(dim(data)[2]-1)]*1e6 - data[,(dim(data)[2])]*1e4)/1e2)
data[,(dim(data)[2])+1] <- (data[,1] - data[,(dim(data)[2]-3)]*1e8 - data[,(dim(data)[2]-2)]*1e6 - data[,(dim(data)[2]-1)]*1e4 - data[,(dim(data)[2])]*1e2)
colnames(data)[(dim(data)[2]-4):dim(data)[2]] <- c("Year","Month","Day","Hour","Minute")
data[,(dim(data)[2]+1)] <- as.Date(with(data, paste(Year, Month, Day, sep="-")), "%Y-%m-%d") # Convert time into a standard format
colnames(data)[(dim(data)[2])] <- c("Date")
names(data)
colnames(data)[colnames(data) == "RH_1_1_1"] <- "RH" # Change column names, if needed for follow-up processing
data <- data %>%
  mutate_at(vars(contains("TS")), ~na_if(.,-9999)) # Change -9999 values into NA. An example is given for soil temperature, but depending on the variables available this can be applied to more columns
temp_list <- ls(name = data, pattern = "TS")
temp_list 
mt_list = c("PA","ALB","RH","WS","WD") # Depending on the available variable, this list can be longer
data <- data %>%
  mutate_at(vars(any_of(mt_list)), ~na_if(.,-9999))
list = c(temp_list, mt_list) # Some sites have soil moisture data under the 'SWC' group, which can be added here if applicable 
for (i in 3:(dim(data)[2]-6)){
  if (colnames(data)[i] %in% list){
    prop <- colnames(data)[i]
    var_name <- paste0(prop,"_availability")
    assign(var_name, data %>% group_by(Date) %>%
             summarise(!!var_name := 1-sum(is.na(!!sym(prop))/length(!!sym(prop)))))
  }
}
comb <- reduce(mget(ls(pattern="availability")), full_join)
comb <- comb %>%
  mutate_at(vars(contains("availability")),
            function(x) (x*100))
# Export data availability summary file
write.csv(comb,"sitename/sitename_availability.csv")

# Combine refined data with the availability file;
var_list <- c('Date',temp_list,mt_list)
data_refined <- data[,c(colnames(data) %in% var_list)]
mt_list_a = c("PA","RH","WS","WD","ALB")
#mt_list_a <- ls(name = data, pattern = "(TA)|(PA)|(RH)|(WS)|(WD)|(VPD)|(ALB)") # For sites with more complicated column names
#mt_list_s <- c("P","SW_IN") # Some sites have variables that needs to be summed up, such as precipitation
var_list_a <- c(temp_list,mt_list_a)
for (i in 2:(dim(data_refined)[2])){
  if (colnames(data_refined)[i] %in% var_list_a){
    prop <- colnames(data_refined)[i]
    var_name <- paste0(prop,"_avg") # Calculate average values
    assign(var_name, data_refined %>% group_by(Date) %>%
             summarise(!!var_name := mean(!!sym(prop), na.rm=TRUE)))
  }
}
for (i in 2:(dim(data_refined)[2])){
  if (colnames(data_refined)[i] %in% var_list_a){
    prop <- colnames(data_refined)[i]
    var_name <- paste0(prop,"_std") # Calculate standard deviations
    assign(var_name, data_refined %>% group_by(Date) %>%
             summarise(!!var_name := sd(!!sym(prop), na.rm=TRUE)))
  }
}
# Combine files and use the same format for NAs
comb_avg <- reduce(mget(ls(pattern="_avg")), full_join)
comb_std <- reduce(mget(ls(pattern="_std")), full_join)
comb_all <- merge(comb_avg, comb_std, by = "Date", all = TRUE)
comb_all <- cbind(comb_all[,1],
                  comb_all[,-1] %>%
                    mutate_all(~ifelse(is.nan(.), NA, .)))
colnames(comb_all)[1] <-"Date"

# Filter by date before matching the two datasets
maxd <- max(data_reorg$date)
mind <- min(data_reorg$date)
comb_all <- comb_all[comb_all$Date >= mind & comb_all$Date <= maxd,]
comb_all <- cbind(data_reorg,comb_all)
comb_all <- comb_all[,-which(names(comb_all) %in% c("Date"))]
colnames(comb_all)[c(1:9)] <- c("Date","GPP_PI_F_avg","NEE_PI_F_avg","RECO_PI_F_avg","TA_avg","VPD_avg","P_sum","SW_IN_sum","SWC_avg")
comb_all$RECO_PI_F_avg <- -comb_all$RECO_PI_F_avg # Change the sign if needed
temp_list <- c("TS_4_1_1_avg","TS_PI_1_avg","TS_PI_2_avg","TS_PI_3_avg") #Customize with relevant column names
comb_all$TS_avg <- rowMeans(comb_all[,temp_list], na.rm=TRUE)
write.csv(comb_all,"sitename/sitename_summary.csv")

# Alternative sites where more variables need to be summarized
#Calculate availability and average files
# a_list = c(flux_list,temp_list,mois_list,mt_list_a)
# all_list = c(flux_list,temp_list,mois_list,mt_list_a,mt_list_s)
# for (i in 3:(dim(data)[2]-6)){
#   if (colnames(data)[i] %in% all_list){
#     prop <- colnames(data)[i]
#     var_name <- paste(prop,"_availability")
#     assign(var_name, data %>% group_by(Date) %>%
#              summarise(!!var_name := 1-sum(is.na(!!sym(prop))/length(!!sym(prop)))))
#   }
# }
# for (i in 3:(dim(data)[2]-6)){
#   if (colnames(data)[i] %in% a_list){
#     prop <- colnames(data)[i]
#     var_name <- paste0(prop,"_avg")
#     assign(var_name, data %>% group_by(Date) %>%
#              summarise(!!var_name := mean(!!sym(prop), na.rm=TRUE)))
#   }
# }
# for (i in 3:(dim(data)[2]-6)){
#   if (colnames(data)[i] %in% mt_list_s){
#     prop <- colnames(data)[i]
#     var_name <- paste0(prop,"_sum")
#     assign(var_name, data %>% group_by(Date) %>%
#              summarise(!!var_name := sum(!!sym(prop), na.rm=TRUE)))
#   }
# }
# for (i in 3:(dim(data)[2]-6)){
#   if (colnames(data)[i] %in% all_list){
#     prop <- colnames(data)[i]
#     var_name <- paste0(prop,"_std")
#     assign(var_name, data %>% group_by(Date) %>%
#              summarise(!!var_name := sd(!!sym(prop), na.rm=TRUE)))
#   }
# }
# comb <- reduce(mget(ls(pattern="availability")), full_join)
# comb <- comb %>%
#   mutate_at(vars(contains("availability")),
#             function(x) (x*100))
# comb_avg <- reduce(mget(ls(pattern="_avg")), full_join)
# comb_sum <- reduce(mget(ls(pattern="_sum")), full_join)
# comb_std <- reduce(mget(ls(pattern="_std")), full_join)
# comb$Date <- as.Date(comb$Date, format = "%Y-%m-%d")
# comb_avg$Date <- as.Date(comb_avg$Date, format = "%Y-%m-%d")
# comb_sum$Date <- as.Date(comb_sum$Date, format = "%Y-%m-%d")
# comb_std$Date <- as.Date(comb_std$Date, format = "%Y-%m-%d")
# comb_all <- merge(merge(merge(comb, comb_avg, by = "Date", all = TRUE),
#                         comb_sum, by = "Date", all = TRUE),
#                   comb_std, by = "Date", all = TRUE)
# comb_all <- cbind(comb_all[,1],
#                   comb_all[,-1] %>%
#                     mutate_all(~ifelse(is.nan(.), NA, .)))
# colnames(comb_all)[1] <-"Date"
# rm(list = ls(pattern="_availability"))
# rm(list = ls(pattern="_avg"))
# rm(list = ls(pattern="_std"))
# rm(list = ls(pattern="_sum"))

# Alternative sites might need to have outlier data removed
# names(comb_all) <- gsub(x = names(comb_all), pattern = " _availability", replacement = "_availability")
# comb_all <- comb_all %>% 
#   mutate (Outlier_GPP_RECO = if_else (((GPP_PI_F_avg < 0)&(!is.na(GPP_PI_F_avg)))|((RECO_PI_F_avg < 0)&(!is.na(RECO_PI_F_avg))), "Y","N")) 
# comb_all <- comb_all %>%
#   mutate_at(vars(contains("GPP")), ~replace(., .< 0, NA))
# comb_all <- comb_all %>%
#   mutate_at(vars(contains("RECO")), ~replace(., .< 0, NA)) 
# comb_all <- comb_all %>% 
#   mutate (Low_availability = if_else ((GPP_PI_F_availability < 40)|(RECO_PI_F_availability < 70), "Y","N")) 
# comb_all <- within (comb_all, GPP_PI_F_avg[GPP_PI_F_availability < 40] <- NA)
# comb_all <- within (comb_all, RECO_PI_F_avg[RECO_PI_F_availability < 70] <- NA)
# comb_all <- comb_all %>% select(-contains("availability"))
# comb_all$RECO_PI_F_avg <- -comb_all$RECO_PI_F_avg
# comb_all$NEE_PI_F_avg <- comb_all$GPP_PI_F_avg + comb_all$RECO_PI_F_avg

# Visualization with interactive charts in R
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