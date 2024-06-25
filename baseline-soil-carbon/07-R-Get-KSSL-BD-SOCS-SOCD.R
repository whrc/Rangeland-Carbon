library(tidyverse)
library(data.table)
setwd("set-pathway")


###############################################################
# Import KSSL and KSSL covariate datasets
KSSL<- read.table('./KSSL.f.2022.csv', comment.char ="", quote = "\"", header = T, sep = ",")
KSSL_BD<- read.table('./KSSL_BD.csv', comment.char ="", quote = "\"", header = T, sep = ",")
KSSL_layer<- read.table('./KSSL.CONUS.V2.csv', comment.char ="", quote = "\"", header = T, sep = ",")


###############################################################
# Merge and process datasets to get KSSL BD dataset
KSSL_all <- merge(KSSL,KSSL_layer,by ="lay_id")
# Check and remove duplicated records if any
check <- data.frame(KSSL_all[duplicated(KSSL_all$X),])
KSSL_all <- KSSL_all[!duplicated(KSSL_all$X), ]
drop <- c("oc_usda.y","med_dp")
KSSL_all = KSSL_all[,!(names(KSSL_all) %in% drop)]
colnames(KSSL_all)[colnames(KSSL_all) == "oc_usda.x"] <- "oc_usda"
# Retain only needed columns
drop <- c("BD_method","BD_core","BD_clod")
KSSL_BD = KSSL_BD[,!(names(KSSL_BD) %in% drop)]
colnames(KSSL_BD)[colnames(KSSL_BD) == "BD"] <- "bd_kssl"
KSSL_all <- merge(KSSL_all, KSSL_BD, by.x="sampleids", by.y="sample_id")
drop <- c("Num.y")
KSSL_all = KSSL_all[,!(names(KSSL_all) %in% drop)]
colnames(KSSL_all)[colnames(KSSL_all) == "Num.x"] <- "Num"
write.csv(KSSL_all,"./KSSL.f.w.BD.csv")


###############################################################
# Merge datasets to get KSSL SOC density
KSSL_all$density <- KSSL_all$oc_usda * KSSL_all$bd_kssl *100 # Calculate SOC density in g/cm3
write.csv(KSSL_all,"./KSSL.f.SOCSd.csv")


###############################################################
# Merge and process datasets to get KSSL SOC stocks for 0-30 cm
KSSL_ID<- read.table('.CSV-with-ID-filename.csv', comment.char ="", quote = "\"", header = T, sep = ",")
colnames(KSSL_ID)[XX] <- "sampleids" # Specify XX with column number corresponds to sample id
KSSL_30 <- KSSL_all[which(KSSL_all$layer_upper < 30),]
KSSL_30 <- KSSL_30[which((KSSL_30$layer_upper < 25)|(KSSL_30$layer_lower < 45)),]
KSSL_30 <- merge(KSSL_30,KSSL_ID,by = c("lay_id","sampleids"))
KSSL_30 <- KSSL_30 %>% arrange(pedon_id, Med_depth) # Rearrangement is necessary to assign the correct layer num
KSSL_30$lay_count <- ave(KSSL_30$layer_upper, KSSL_30$pedon_id, FUN = seq_along)
KSSL_30_count <- KSSL_30 %>% count(pedon_id)
KSSL_30 <- merge(KSSL_30,KSSL_30_count,by="pedon_id")
KSSL_30_1 <- KSSL_30[which(KSSL_30$n == 1),] # Pedons with only one layer
KSSL_30_1$SOCS <- KSSL_30_1$bd_kssl * KSSL_30_1$oc_usda * 30 * 100 
KSSL_30_2 <- KSSL_30[which(KSSL_30$n == 2),] # Pedons with two layers
pedon_2 <- KSSL_30_2[,c("pedon_id","lay_count","layer_upper","layer_lower")]
pedon_2 <- reshape(pedon_2,direction = "wide",idvar = "pedon_id",timevar = "lay_count")
pedon_2$layer_med <- (pedon_2$layer_lower.1 + pedon_2$layer_upper.2)/2
pedon_2 <- pedon_2[,c("pedon_id","layer_med")]
KSSL_30_2 <- merge(KSSL_30_2,pedon_2,by="pedon_id")
KSSL_30_2_1 <- KSSL_30_2[which(KSSL_30_2$lay_count == 1),]
KSSL_30_2_2 <- KSSL_30_2[which(KSSL_30_2$lay_count == 2),]
KSSL_30_2_1$SOCS <- KSSL_30_2_1$bd_kssl * KSSL_30_2_1$oc_usda * KSSL_30_2_1$layer_med * 100
KSSL_30_2_2$SOCS <- KSSL_30_2_2$bd_kssl * KSSL_30_2_2$oc_usda * (30 - KSSL_30_2_2$layer_med) * 100
KSSL_30_2 <- rbind(KSSL_30_2_1,KSSL_30_2_2)
drop <- c("layer_med")
KSSL_30_2 = KSSL_30_2[,!(names(KSSL_30_2) %in% drop)]
KSSL_30_3 <- KSSL_30[which(KSSL_30$n == 3),] # Pedons with three layers
pedon_3 <- KSSL_30_3[,c("pedon_id","lay_count","layer_upper","layer_lower")]
pedon_3 <- reshape(pedon_3,direction = "wide",idvar = "pedon_id",timevar = "lay_count")
pedon_3$layer_med1 <- (pedon_3$layer_lower.1 + pedon_3$layer_upper.2)/2
pedon_3$layer_med2 <- (pedon_3$layer_lower.2 + pedon_3$layer_upper.3)/2
pedon_3 <- pedon_3[,c("pedon_id","layer_med1","layer_med2")]
KSSL_30_3 <- merge(KSSL_30_3,pedon_3,by="pedon_id")
KSSL_30_3_1 <- KSSL_30_3[which(KSSL_30_3$lay_count == 1),]
KSSL_30_3_2 <- KSSL_30_3[which(KSSL_30_3$lay_count == 2),]
KSSL_30_3_3 <- KSSL_30_3[which(KSSL_30_3$lay_count == 3),]
KSSL_30_3_1$SOCS <- KSSL_30_3_1$bd_kssl * KSSL_30_3_1$oc_usda * KSSL_30_3_1$layer_med1 * 100
KSSL_30_3_2$SOCS <- KSSL_30_3_2$bd_kssl * KSSL_30_3_2$oc_usda * (KSSL_30_3_2$layer_med2 - KSSL_30_3_2$layer_med1) * 100
KSSL_30_3_3$SOCS <- KSSL_30_3_3$bd_kssl * KSSL_30_3_3$oc_usda * (30 - KSSL_30_3_3$layer_med2) * 100
KSSL_30_3 <- rbind(KSSL_30_3_1,KSSL_30_3_2,KSSL_30_3_3)
drop <- c("layer_med1","layer_med2")
KSSL_30_3 = KSSL_30_3[,!(names(KSSL_30_3) %in% drop)]
KSSL_30_4 <- KSSL_30[which(KSSL_30$n == 4),] # Pedons with four layers
pedon_4 <- KSSL_30_4[,c("pedon_id","lay_count","layer_upper","layer_lower")]
pedon_4 <- reshape(pedon_4,direction = "wide",idvar = "pedon_id",timevar = "lay_count")
pedon_4$layer_med1 <- (pedon_4$layer_lower.1 + pedon_4$layer_upper.2)/2
pedon_4$layer_med2 <- (pedon_4$layer_lower.2 + pedon_4$layer_upper.3)/2
pedon_4$layer_med3 <- (pedon_4$layer_lower.3 + pedon_4$layer_upper.4)/2
pedon_4 <- pedon_4[,c("pedon_id","layer_med1","layer_med2","layer_med3")]
KSSL_30_4 <- merge(KSSL_30_4,pedon_4,by="pedon_id")
KSSL_30_4_1 <- KSSL_30_4[which(KSSL_30_4$lay_count == 1),]
KSSL_30_4_2 <- KSSL_30_4[which(KSSL_30_4$lay_count == 2),]
KSSL_30_4_3 <- KSSL_30_4[which(KSSL_30_4$lay_count == 3),]
KSSL_30_4_4 <- KSSL_30_4[which(KSSL_30_4$lay_count == 4),]
KSSL_30_4_1$SOCS <- KSSL_30_4_1$bd_kssl * KSSL_30_4_1$oc_usda * KSSL_30_4_1$layer_med1 * 100
KSSL_30_4_2$SOCS <- KSSL_30_4_2$bd_kssl * KSSL_30_4_2$oc_usda * (KSSL_30_4_2$layer_med2 - KSSL_30_4_2$layer_med1) * 100
KSSL_30_4_3$SOCS <- KSSL_30_4_3$bd_kssl * KSSL_30_4_3$oc_usda * (KSSL_30_4_3$layer_med3 - KSSL_30_4_3$layer_med2) * 100
KSSL_30_4_4$SOCS <- KSSL_30_4_4$bd_kssl * KSSL_30_4_4$oc_usda * (30 - KSSL_30_4_4$layer_med3) * 100
KSSL_30_4 <- rbind(KSSL_30_4_1,KSSL_30_4_2,KSSL_30_4_3,KSSL_30_4_4)
drop <- c("layer_med1","layer_med2","layer_med3")
KSSL_30_4 = KSSL_30_4[,!(names(KSSL_30_4) %in% drop)]
KSSL_30_5 <- KSSL_30[which(KSSL_30$n == 5),] # Pedons with five layers
pedon_5 <- KSSL_30_5[,c("pedon_id","lay_count","layer_upper","layer_lower")]
pedon_5 <- reshape(pedon_5,direction = "wide",idvar = "pedon_id",timevar = "lay_count")
pedon_5$layer_med1 <- (pedon_5$layer_lower.1 + pedon_5$layer_upper.2)/2
pedon_5$layer_med2 <- (pedon_5$layer_lower.2 + pedon_5$layer_upper.3)/2
pedon_5$layer_med3 <- (pedon_5$layer_lower.3 + pedon_5$layer_upper.4)/2
pedon_5$layer_med4 <- (pedon_5$layer_lower.4 + pedon_5$layer_upper.5)/2
pedon_5 <- pedon_5[,c("pedon_id","layer_med1","layer_med2","layer_med3","layer_med4")]
KSSL_30_5 <- merge(KSSL_30_5,pedon_5,by="pedon_id")
KSSL_30_5_1 <- KSSL_30_5[which(KSSL_30_5$lay_count == 1),]
KSSL_30_5_2 <- KSSL_30_5[which(KSSL_30_5$lay_count == 2),]
KSSL_30_5_3 <- KSSL_30_5[which(KSSL_30_5$lay_count == 3),]
KSSL_30_5_4 <- KSSL_30_5[which(KSSL_30_5$lay_count == 4),]
KSSL_30_5_5 <- KSSL_30_5[which(KSSL_30_5$lay_count == 5),]
KSSL_30_5_1$SOCS <- KSSL_30_5_1$bd_kssl * KSSL_30_5_1$oc_usda * KSSL_30_5_1$layer_med1 * 100
KSSL_30_5_2$SOCS <- KSSL_30_5_2$bd_kssl * KSSL_30_5_2$oc_usda * (KSSL_30_5_2$layer_med2 - KSSL_30_5_2$layer_med1) * 100
KSSL_30_5_3$SOCS <- KSSL_30_5_3$bd_kssl * KSSL_30_5_3$oc_usda * (KSSL_30_5_3$layer_med3 - KSSL_30_5_3$layer_med2) * 100
KSSL_30_5_4$SOCS <- KSSL_30_5_4$bd_kssl * KSSL_30_5_4$oc_usda * (KSSL_30_5_4$layer_med4 - KSSL_30_5_4$layer_med3) * 100
KSSL_30_5_5$SOCS <- KSSL_30_5_5$bd_kssl * KSSL_30_5_5$oc_usda * (30 - KSSL_30_5_5$layer_med4) * 100
KSSL_30_5 <- rbind(KSSL_30_5_1,KSSL_30_5_2,KSSL_30_5_3,KSSL_30_5_4,KSSL_30_5_5)
drop <- c("layer_med1","layer_med2","layer_med3","layer_med4")
KSSL_30_5 = KSSL_30_5[,!(names(KSSL_30_5) %in% drop)]
KSSL_30_6 <- KSSL_30[which(KSSL_30$n == 6),] # Pedons with six layers
KSSL_30_6_1 <- KSSL_30_6[which(KSSL_30_6$lay_count < 6),]
KSSL_30_6_2 <- KSSL_30_6[which(KSSL_30_6$lay_count == 6),]
KSSL_30_6_1$SOCS <- KSSL_30_6_1$bd_kssl * KSSL_30_6_1$oc_usda * (KSSL_30_6_1$layer_lower - KSSL_30_6_1$layer_upper) * 100
KSSL_30_6_2$SOCS <- KSSL_30_6_2$bd_kssl * KSSL_30_6_2$oc_usda * (30 - KSSL_30_6_2$layer_upper) * 100
KSSL_30_6 <- rbind(KSSL_30_6_1,KSSL_30_6_2)
KSSL_30_8 <- KSSL_30[which(KSSL_30$n == 8),] # Pedons with eight layers
KSSL_30_8 <- KSSL_30_8[-c(2,4,6,8),] # Remove repeated pedons manually
KSSL_30_8 <- KSSL_30_8[-c(5,7,8,10),]
KSSL_30_8 <- KSSL_30_8[-c(10,12,14,16),]
KSSL_30_8_1 <- KSSL_30_8[which(KSSL_30_8$layer_lower < 30),]
KSSL_30_8_2 <- KSSL_30_8[which(KSSL_30_8$layer_lower > 30),]
KSSL_30_8_1$SOCS <- KSSL_30_8_1$bd_kssl * KSSL_30_8_1$oc_usda * (KSSL_30_8_1$layer_lower - KSSL_30_8_1$layer_upper) * 100
KSSL_30_8_2$SOCS <- KSSL_30_8_2$bd_kssl * KSSL_30_8_2$oc_usda * (30 - KSSL_30_8_2$layer_upper) * 100
KSSL_30_8 <- rbind(KSSL_30_8_1,KSSL_30_8_2)
KSSL_30_10 <- KSSL_30[which(KSSL_30$n == 10),] # Pedons with ten layers
KSSL_30_10 <- KSSL_30_10[-c(2,4,6,8,10,12,14,16,18,20),] # Remove repeated pedons manually
KSSL_30_10_1 <- KSSL_30_10[which(KSSL_30_10$layer_lower < 30),]
KSSL_30_10_2 <- KSSL_30_10[which(KSSL_30_10$layer_lower > 30),]
KSSL_30_10_1$SOCS <- KSSL_30_10_1$bd_kssl * KSSL_30_10_1$oc_usda * (KSSL_30_10_1$layer_lower - KSSL_30_10_1$layer_upper) * 100
KSSL_30_10_2$SOCS <- KSSL_30_10_2$bd_kssl * KSSL_30_10_2$oc_usda * (30 - KSSL_30_10_2$layer_upper) * 100
KSSL_30_10 <- rbind(KSSL_30_10_1,KSSL_30_10_2)
# Merge and get final SOCS cores joined with covariates
KSSL_f <- rbind(KSSL_30_1,KSSL_30_2,KSSL_30_3,KSSL_30_4,KSSL_30_5,KSSL_30_6,KSSL_30_8,KSSL_30_10)
SOCS <- KSSL_f %>% group_by(pedon_id) %>% summarise(socs_kssl = sum(SOCS))
KSSL.f <- KSSL_f[which(KSSL_f$lay_count == 1),]
KSSL.f <- merge(KSSL.f,SOCS,by="pedon_id")
data0<- read.table('./KSSL.cov.csv',comment.char ="", quote = "\"", header = T, sep = ",") #Join data with soil properties at 15 cm depth for the cores
data0 <- data0[,c("smp_id","BD3","SOC3","Clay3","Sand3")]
KSSL_f <- merge(KSSL.f,data0,by.x ="sampleids",by.y ="smp_id")
drop <- c("lay_count","n","SOCS","BD","SOC","Clay","Sand")
KSSL_f = KSSL_f[,!(names(KSSL_f) %in% drop)]
colnames(KSSL_f)[colnames(KSSL_f) == "BD3"] <- "BD"
colnames(KSSL_f)[colnames(KSSL_f) == "SOC3"] <- "SOC"
colnames(KSSL_f)[colnames(KSSL_f) == "Clay3"] <- "Clay"
colnames(KSSL_f)[colnames(KSSL_f) == "Sand3"] <- "Sand"
KSSL_f$SOC <- KSSL_f$SOC/10 # Change unit to %
write.csv(KSSL_f,"./KSSL.f.SOCS.csv")

