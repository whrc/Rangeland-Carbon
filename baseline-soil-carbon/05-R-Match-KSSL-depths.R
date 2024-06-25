setwd("set-pathway")


###############################################################
# Import extracted V1 KSSL covariate datasets
data0 <- read.table('./KSSL_cov.csv', comment.char ="", quote = "\"", header = T, sep = ",")
data0$SOC <- data0$SOC1/10 # Change unit to %
data0$BD <- data0$BD1
data0$Clay <- data0$Clay1
data0$Sand <- data0$Sand1
data0[((data0$Med_depth >= 2.5)&(data0$Med_depth < 7.5)),]$SOC <- data0[((data0$Med_depth >= 2.5)&(data0$Med_depth < 7.5)),]$SOC2/10
data0[((data0$Med_depth >= 2.5)&(data0$Med_depth < 7.5)),]$BD <- data0[((data0$Med_depth >= 2.5)&(data0$Med_depth < 7.5)),]$BD2
data0[((data0$Med_depth >= 2.5)&(data0$Med_depth < 7.5)),]$Clay <- data0[((data0$Med_depth >= 2.5)&(data0$Med_depth < 7.5)),]$Clay2
data0[((data0$Med_depth >= 2.5)&(data0$Med_depth < 7.5)),]$Sand <- data0[((data0$Med_depth >= 2.5)&(data0$Med_depth < 7.5)),]$Sand2
data0[((data0$Med_depth >= 7.5)&(data0$Med_depth < 22.5)),]$SOC <- data0[((data0$Med_depth >= 7.5)&(data0$Med_depth < 22.5)),]$SOC3/10
data0[((data0$Med_depth >= 7.5)&(data0$Med_depth < 22.5)),]$BD <- data0[((data0$Med_depth >= 7.5)&(data0$Med_depth < 22.5)),]$BD3
data0[((data0$Med_depth >= 7.5)&(data0$Med_depth < 22.5)),]$Clay <- data0[((data0$Med_depth >= 7.5)&(data0$Med_depth < 22.5)),]$Clay3
data0[((data0$Med_depth >= 7.5)&(data0$Med_depth < 22.5)),]$Sand <- data0[((data0$Med_depth >= 7.5)&(data0$Med_depth < 22.5)),]$Sand3 
data0[((data0$Med_depth >= 22.5)&(data0$Med_depth < 37.5)),]$SOC <- data0[((data0$Med_depth >= 22.5)&(data0$Med_depth < 37.5)),]$SOC4/10
data0[((data0$Med_depth >= 22.5)&(data0$Med_depth < 37.5)),]$BD <- data0[((data0$Med_depth >= 22.5)&(data0$Med_depth < 37.5)),]$BD4
data0[((data0$Med_depth >= 22.5)&(data0$Med_depth < 37.5)),]$Clay <- data0[((data0$Med_depth >= 22.5)&(data0$Med_depth < 37.5)),]$Clay4
data0[((data0$Med_depth >= 22.5)&(data0$Med_depth < 37.5)),]$Sand <- data0[((data0$Med_depth >= 22.5)&(data0$Med_depth < 37.5)),]$Sand4
data0[((data0$Med_depth >= 37.5)&(data0$Med_depth < 82.5)),]$SOC <- data0[((data0$Med_depth >= 37.5)&(data0$Med_depth < 82.5)),]$SOC5/10
data0[((data0$Med_depth >= 37.5)&(data0$Med_depth < 82.5)),]$BD <- data0[((data0$Med_depth >= 37.5)&(data0$Med_depth < 82.5)),]$BD5
data0[((data0$Med_depth >= 37.5)&(data0$Med_depth < 82.5)),]$Clay <- data0[((data0$Med_depth >= 37.5)&(data0$Med_depth < 82.5)),]$Clay5
data0[((data0$Med_depth >= 37.5)&(data0$Med_depth < 82.5)),]$Sand <- data0[((data0$Med_depth >= 37.5)&(data0$Med_depth < 82.5)),]$Sand5
data0[((data0$Med_depth >= 82.5)&(data0$Med_depth < 117.5)),]$SOC <- data0[((data0$Med_depth >= 82.5)&(data0$Med_depth < 117.5)),]$SOC6/10
data0[((data0$Med_depth >= 82.5)&(data0$Med_depth < 117.5)),]$BD <- data0[((data0$Med_depth >= 82.5)&(data0$Med_depth < 117.5)),]$BD6
data0[((data0$Med_depth >= 82.5)&(data0$Med_depth < 117.5)),]$Clay <- data0[((data0$Med_depth >= 82.5)&(data0$Med_depth < 117.5)),]$Clay6
data0[((data0$Med_depth >= 82.5)&(data0$Med_depth < 117.5)),]$Sand <- data0[((data0$Med_depth >= 82.5)&(data0$Med_depth < 117.5)),]$Sand6
data0[(data0$Med_depth >= 117.5),]$SOC <- data0[(data0$Med_depth >= 117.5),]$SOC7/10
data0[(data0$Med_depth >= 117.5),]$BD <- data0[(data0$Med_depth >= 117.5),]$BD7
data0[(data0$Med_depth >= 117.5),]$Clay <- data0[(data0$Med_depth >= 117.5),]$Clay7
data0[(data0$Med_depth >= 117.5),]$Sand <- data0[(data0$Med_depth >= 117.5),]$Sand7
drop <- c("FIPS","SOC1","SOC2","SOC3","SOC4","SOC5","SOC6","SOC7","BD1","BD2","BD3","BD4","BD5","BD6","BD7",
          "Clay1","Clay2","Clay3","Clay4","Clay5","Clay6","Clay7","Sand1","Sand2","Sand3","Sand4","Sand5","Sand6","Sand7")
data0 = data0[,!(names(data0) %in% drop)]
colnames(data0)[2] <- "sampleids"
write.csv(data0, "./KSSL_cov_2022.txt")