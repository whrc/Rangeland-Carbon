library(data.table)
setwd("C:/Users/username/path-to-the-modeled-results-data-folder")

# Import datasets and process to the format needed for comparison
comp_all <- read.table('file-with-modeled-versus-measured-results-for-the-vegetation-type.csv', header=T, sep = ",")
colnames(comp_all)[3] <- "NEE"
comp_all$Year <- as.numeric(tstrsplit(comp_all$Date, "-")[[1]])
comp_all$Month <- as.numeric(tstrsplit(comp_all$Date, "-")[[2]])
comp_all$site_num = as.numeric(as.factor(with(comp_all,paste(Site))))
replace.values <- function(search, replace, x){
  return(replace[ match(x, search) ])
} # Only needed if more than one file (different time periods) correspond to a site
comp_all <- do.call("rbind", by(comp_all, comp_all[,"Site"], function(x){
  ux <- unique(x[,"Year"])
  return(cbind(x, subgroup=replace.values(ux, 1:length(ux), x[,"Year"])))
}))
all_year_NEE <- aggregate(NEE~Site+Year,data = comp_all, mean)
all_year_NEE_measured <- aggregate(NEE_measured~Site+Year,data = comp_all, mean)
all_year <- merge(all_year_NEE,all_year_NEE_measured, by = c("Site","Year"))
all_year$NEE <- all_year$NEE*365
all_year$NEE_measured <- all_year$NEE_measured*365

# Calculate results for site-year combinations
lr_RCTMs_year <- lm(NEE_measured ~ NEE, data = all_year)
summary(lr_RCTMs_year)
RMSE <- sqrt(mean((all_year$NEE_measured - all_year$NEE)^2))
RMSE
all_year$site_num = as.numeric(as.factor(with(all_year,paste(Site))))
RMSE_I <- data.frame()
r2_I <- data.frame()
MBE_I <- data.frame()
site_I <- data.frame()
site_num_I <- data.frame()
for (i in 1:length(unique(all_year$Site))){
  all_year_i <- all_year[which(all_year$site_num == i),]
  if((dim(all_year_i)[1] > 2)&(sd(all_year_i$NEE_measured)> 0)){
    r2 <- cor(all_year_i$NEE_measured, all_year_i$NEE)^2
    mbe <- mean(all_year_i$NEE_measured - all_year_i$NEE)
    RMSE <- sqrt(mean((all_year_i$NEE_measured - all_year_i$NEE)^2))
    r2_I <- rbind(r2_I, r2)
    MBE_I <- rbind(MBE_I, mbe)
    RMSE_I <- rbind(RMSE_I, RMSE)
    site_I <- rbind(site_I, unique(all_year_i$Site))
    site_num_I <- rbind(site_num_I, i)}
}
individual_result <- data.frame(site_I, site_num_I,r2_I,RMSE_I,MBE_I)
colnames(individual_result)[3:5] <- c("R2","RMSE","MBE")
colMeans(individual_result[3:5])
sd(individual_result$R2)
sd(individual_result$RMSE)
sd(individual_result$MBE)

# Calculate monthly results
all_month_NEE <- aggregate(NEE~Site+Year+Month,data = comp_all, mean)
all_month_NEE_measured <- aggregate(NEE_measured~Site+Year+Month,data = comp_all, mean)
all_month <- merge(all_month_NEE,all_month_NEE_measured, by = c("Site","Year","Month"))
all_month$NEE <- all_month$NEE*30
all_month$NEE_measured <- all_month$NEE_measured*30
lr_RCTMs_month <- lm(NEE_measured ~ NEE, data = all_month)
summary(lr_RCTMs_month)
RMSE <- sqrt(mean((all_month$NEE_measured - all_month$NEE)^2))
RMSE
all_month$site_num = as.numeric(as.factor(with(all_month,paste(Site))))
RMSE_I <- data.frame()
r2_I <- data.frame()
MBE_I <- data.frame()
site_I <- data.frame()
site_num_I <- data.frame()
for (i in 1:length(unique(all_month$Site))){
  all_month_i <- all_month[which(all_month$site_num == i),]
  if((dim(all_month_i)[1] > 2)&(sd(all_month_i$NEE_measured)> 0)){
    r2 <- cor(all_month_i$NEE_measured, all_month_i$NEE)^2
    mbe <- mean(all_month_i$NEE_measured - all_month_i$NEE)
    RMSE <- sqrt(mean((all_month_i$NEE_measured - all_month_i$NEE)^2))
    r2_I <- rbind(r2_I, r2)
    MBE_I <- rbind(MBE_I, mbe)
    RMSE_I <- rbind(RMSE_I, RMSE)
    site_I <- rbind(site_I, unique(all_month_i$Site))
    site_num_I <- rbind(site_num_I, i)}
}
individual_result <- data.frame(site_I,site_num_I,r2_I,RMSE_I,MBE_I)
colnames(individual_result) <- c("Site","Site_num","R2","RMSE","MBE")
colMeans(individual_result[3:5])
sd(individual_result$R2)
sd(individual_result$RMSE)
sd(individual_result$MBE)

# Calculate seasonal results
comparison_s1 <- comp_all[which((comp_all$Month == 12)|(comp_all$Month == 1)|(comp_all$Month == 2)),]
comparison_s2 <- comp_all[which((comp_all$Month == 3)|(comp_all$Month == 4)|(comp_all$Month == 5)),]
comparison_s3 <- comp_all[which((comp_all$Month == 6)|(comp_all$Month == 7)|(comp_all$Month == 8)),]
comparison_s4 <- comp_all[which((comp_all$Month == 9)|(comp_all$Month == 10)|(comp_all$Month == 11)),]
RCTM_s1_NEE <- aggregate(NEE~Site+Year,data = comparison_s1, mean)
RCTM_s1_NEE_measured <- aggregate(NEE_measured~Site+Year,data = comparison_s1, mean)
RCTM_s2_NEE <- aggregate(NEE~Site+Year,data = comparison_s2, mean)
RCTM_s2_NEE_measured <- aggregate(NEE_measured~Site+Year,data = comparison_s2, mean)
RCTM_s3_NEE <- aggregate(NEE~Site+Year,data = comparison_s3, mean)
RCTM_s3_NEE_measured <- aggregate(NEE_measured~Site+Year,data = comparison_s3, mean)
RCTM_s4_NEE <- aggregate(NEE~Site+Year,data = comparison_s4, mean)
RCTM_s4_NEE_measured <- aggregate(NEE_measured~Site+Year,data = comparison_s4, mean)
comparison_gs <- comp_all[which((comp_all$Month == 4)|(comp_all$Month == 5)|(comp_all$Month == 6)
                                |(comp_all$Month == 7)|(comp_all$Month == 8)|(comp_all$Month == 9)
                                |(comp_all$Month == 10)),]
RCTM_gs_NEE <- aggregate(NEE~Site+Year,data = comparison_gs, mean)
RCTM_gs_NEE_measured <- aggregate(NEE_measured~Site+Year,data = comparison_gs, mean)

# S1
RCTM_year1 <- merge(RCTM_s1_NEE,RCTM_s1_NEE_measured, by = c("Site","Year"))
RCTM_year1$NEE <- RCTM_year1$NEE*91
RCTM_year1$NEE_measured <- RCTM_year1$NEE_measured*91
lr_RCTM_year1 <- lm(NEE_measured ~ NEE, data = RCTM_year1)
summary(lr_RCTM_year1)
RMSE <- sqrt(mean((RCTM_year1$NEE_measured - RCTM_year1$NEE)^2))
RMSE
RCTM_year1$site_num = as.numeric(as.factor(with(RCTM_year1,paste(Site))))
RMSE_I <- data.frame()
r2_I <- data.frame()
MBE_I <- data.frame()
site_I <- data.frame()
site_num_I <- data.frame()
for (i in 1:length(unique(RCTM_year1$Site))){
  RCTM_year1_i <- RCTM_year1[which(RCTM_year1$site_num == i),]
  if((dim(RCTM_year1_i)[1] > 2)&(sd(RCTM_year1_i$NEE_measured)> 0)){
    r2 <- cor(RCTM_year1_i$NEE_measured, RCTM_year1_i$NEE)^2
    mbe <- mean(RCTM_year1_i$NEE_measured - RCTM_year1_i$NEE)
    RMSE <- sqrt(mean((RCTM_year1_i$NEE_measured - RCTM_year1_i$NEE)^2))
    r2_I <- rbind(r2_I, r2)
    MBE_I <- rbind(MBE_I, mbe)
    RMSE_I <- rbind(RMSE_I, RMSE)
    site_I <- rbind(site_I, unique(RCTM_year1_i$Site))
    site_num_I <- rbind(site_num_I, i)}
}
individual_result <- data.frame(site_I,site_num_I,r2_I,RMSE_I,MBE_I)
colnames(individual_result) <- c("Site","Site_num","R2","RMSE","MBE")
colMeans(individual_result[3:5])
sd(individual_result$R2)
sd(individual_result$RMSE)
sd(individual_result$MBE)

# S2
RCTM_year2 <- merge(RCTM_s2_NEE,RCTM_s2_NEE_measured, by = c("Site","Year"))
RCTM_year2$NEE <- RCTM_year2$NEE*91
RCTM_year2$NEE_measured <- RCTM_year2$NEE_measured*91
lr_RCTM_year2 <- lm(NEE_measured ~ NEE, data = RCTM_year2)
summary(lr_RCTM_year2)
RMSE <- sqrt(mean((RCTM_year2$NEE_measured - RCTM_year2$NEE)^2))
RMSE
RCTM_year2$site_num = as.numeric(as.factor(with(RCTM_year2,paste(Site))))
RMSE_I <- data.frame()
r2_I <- data.frame()
MBE_I <- data.frame()
site_I <- data.frame()
site_num_I <- data.frame()
for (i in 1:length(unique(RCTM_year2$Site))){
  RCTM_year2_i <- RCTM_year2[which(RCTM_year2$site_num == i),]
  if((dim(RCTM_year2_i)[1] > 2)&(sd(RCTM_year2_i$NEE_measured)> 0)){
    r2 <- cor(RCTM_year2_i$NEE_measured, RCTM_year2_i$NEE)^2
    mbe <- mean(RCTM_year2_i$NEE_measured - RCTM_year2_i$NEE)
    RMSE <- sqrt(mean((RCTM_year2_i$NEE_measured - RCTM_year2_i$NEE)^2))
    r2_I <- rbind(r2_I, r2)
    MBE_I <- rbind(MBE_I, mbe)
    RMSE_I <- rbind(RMSE_I, RMSE)
    site_I <- rbind(site_I, unique(RCTM_year2_i$Site))
    site_num_I <- rbind(site_num_I, i)}
}
individual_result <- data.frame(site_I,site_num_I,r2_I,RMSE_I,MBE_I)
colnames(individual_result) <- c("Site","Site_num","R2","RMSE","MBE")
colMeans(individual_result[3:5])
sd(individual_result$R2)
sd(individual_result$RMSE)
sd(individual_result$MBE)

# S3
RCTM_year3 <- merge(RCTM_s3_NEE,RCTM_s3_NEE_measured, by = c("Site","Year"))
RCTM_year3$NEE <- RCTM_year3$NEE*91
RCTM_year3$NEE_measured <- RCTM_year3$NEE_measured*91
lr_RCTM_year3 <- lm(NEE_measured ~ NEE, data = RCTM_year3)
summary(lr_RCTM_year3)
RMSE <- sqrt(mean((RCTM_year3$NEE_measured - RCTM_year3$NEE)^2))
RMSE
RCTM_year3$site_num = as.numeric(as.factor(with(RCTM_year3,paste(Site))))
RMSE_I <- data.frame()
r2_I <- data.frame()
MBE_I <- data.frame()
site_I <- data.frame()
site_num_I <- data.frame()
for (i in 1:length(unique(RCTM_year3$Site))){
  RCTM_year3_i <- RCTM_year3[which(RCTM_year3$site_num == i),]
  if((dim(RCTM_year3_i)[1] > 2)&(sd(RCTM_year3_i$NEE_measured)> 0)){
    r2 <- cor(RCTM_year3_i$NEE_measured, RCTM_year3_i$NEE)^2
    mbe <- mean(RCTM_year3_i$NEE_measured - RCTM_year3_i$NEE)
    RMSE <- sqrt(mean((RCTM_year3_i$NEE_measured - RCTM_year3_i$NEE)^2))
    r2_I <- rbind(r2_I, r2)
    MBE_I <- rbind(MBE_I, mbe)
    RMSE_I <- rbind(RMSE_I, RMSE)
    site_I <- rbind(site_I, unique(RCTM_year3_i$Site))
    site_num_I <- rbind(site_num_I, i)}
}
individual_result <- data.frame(site_I,site_num_I,r2_I,RMSE_I,MBE_I)
colnames(individual_result) <- c("Site","Site_num","R2","RMSE","MBE")
colMeans(individual_result[3:5])
sd(individual_result$R2)
sd(individual_result$RMSE)
sd(individual_result$MBE)

# S4
RCTM_year4 <- merge(RCTM_s4_NEE,RCTM_s4_NEE_measured, by = c("Site","Year"))
RCTM_year4$NEE <- RCTM_year4$NEE*91
RCTM_year4$NEE_measured <- RCTM_year4$NEE_measured*91
lr_RCTM_year4 <- lm(NEE_measured ~ NEE, data = RCTM_year4)
summary(lr_RCTM_year4)
RMSE <- sqrt(mean((RCTM_year4$NEE_measured - RCTM_year4$NEE)^2))
RMSE
RCTM_year4$site_num = as.numeric(as.factor(with(RCTM_year4,paste(Site))))
RMSE_I <- data.frame()
r2_I <- data.frame()
MBE_I <- data.frame()
site_I <- data.frame()
site_num_I <- data.frame()
for (i in 1:length(unique(RCTM_year4$Site))){
  RCTM_year4_i <- RCTM_year4[which(RCTM_year4$site_num == i),]
  if((dim(RCTM_year4_i)[1] > 2)&(sd(RCTM_year4_i$NEE_measured)> 0)){
    r2 <- cor(RCTM_year4_i$NEE_measured, RCTM_year4_i$NEE)^2
    mbe <- mean(RCTM_year4_i$NEE_measured - RCTM_year4_i$NEE)
    RMSE <- sqrt(mean((RCTM_year4_i$NEE_measured - RCTM_year4_i$NEE)^2))
    r2_I <- rbind(r2_I, r2)
    MBE_I <- rbind(MBE_I, mbe)
    RMSE_I <- rbind(RMSE_I, RMSE)
    site_I <- rbind(site_I, unique(RCTM_year4_i$Site))
    site_num_I <- rbind(site_num_I, i)}
}
individual_result <- data.frame(site_I,site_num_I,r2_I,RMSE_I,MBE_I)
colnames(individual_result) <- c("Site","Site_num","R2","RMSE","MBE")
colMeans(individual_result[3:5])
sd(individual_result$R2)
sd(individual_result$RMSE)
sd(individual_result$MBE)
