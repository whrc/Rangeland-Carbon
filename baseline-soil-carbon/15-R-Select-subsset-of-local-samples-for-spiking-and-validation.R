library(tidyverse)
library(ggmap)
library(clhs)
setwd("set-pathway")


###############################################################
# Define parameters, site, and import datasets
# Select site and desired parameters
ROI = "specify-site-name"
sample_num = 25 # Define number of samples desired for testing out the spiking approach
ncomp = 3 # Change cluster number here
# Import site measured dataset and join with covariate datasets
all <- read.table(paste0("./",ROI,"CSV-with-covariates-filename.csv"), 
                  comment.char ="", quote = "\"", header = T, sep = ",")
measured <- read.table("./CSV-with-measurement-filename.csv", 
                  comment.char ="", quote = "\"", header = T, sep = ",")
# Convert scale to match with KSSL dataset, if not yet converted
all$ppt <- all$ppt/20
all$tmin <- all$tmin/20
all$tmax <- all$tmax/20
all$VPD <- all$VPD/20
all$GPP <- all$GPP/20
all$AS <- sin(pi*(all$AS)/360) 
# Keep only soil core samples for clustering analysis
all_cores <- all[!duplicated(all[,c("Field","Num","Cluster","Point")]),]
reduced <- c("Field","Num","Cluster","Point","EL","SL","AS","TWI",
             "ppt","tmin","tmax","VPD","GPP","EVI","Tree",
             "SOC2","SOC5","Clay2","Clay5","Sand2","Sand5","BD2","BD5", # Use two different depths
             "AFGC","PFGC","TREE","SHR","LTR","BG","LULC")
all.reduced <- all_cores[,(names(all_cores) %in% reduced)]
reduce <- c("LULC","Field","Num","Cluster","Point")
subset <- all.reduced[,!(names(all.reduced) %in% reduce)] # Customize and remove unnecessary covarites or columns, if any
df <- data.frame(scale(subset))


###############################################################
# Select spiking samples with clustering analysis + cLHS
# Carry out clustering analysis
set.seed(1)
kmeans_basic <- kmeans(df,centers=ncomp,nstart = 10)
kmeans_basic_table <- data.frame(kmeans_basic$size, kmeans_basic$centers) 
kmeans_basic_df <- data.frame(Cluster = kmeans_basic$cluster, df)
propdf <- data.frame(t(table(kmeans_basic_df$Cluster)))[,2:3]
colnames(propdf) <- c("cluster","counts")
propdf <- propdf[order(propdf$counts),] 
propdf$prop <- propdf$counts/sum(propdf$counts)
propdf$perc <- round(propdf$prop*100) 
geodf <- data.frame(all.reduced,kmeans_basic_df$Cluster)
colnames(geodf)[dim(geodf)[2]] = c("new_cluster")
# Get samples within the cluster using CLHS
samples <- data.frame()
for(i in 1:ncomp){
  oneclust <- geodf[geodf$new_cluster==i,]
  clusterperc <- round(propdf$perc[propdf$cluster==i]*sample_num/100) # Round the data to ensure that the suggested number of samples is an integer
  sample_indices <- clhs(oneclust[,5:28], size = clusterperc, progress = TRUE, iter = 500) # Customize column selection
  thesesamples <- oneclust[sample_indices,]
  samples <- rbind(samples,thesesamples)
}
# Combine and reorganize files
comb <- merge(samples,all[c("Num","Sample_ID","Med_depth","SOC1","SOC3","SOC4","BD1","BD3","BD4",
                            "Clay1","Clay3","Clay4","Sand1","Sand3","Sand4")],by="Num")
spike <- merge(comb,measured[c("Field","Num","Y","X","Med_depth","Measured_SOC")],by=c("Field","Med_depth","Num"))
# Assign soil properties according to the sampling depth
spike$SOC <- spike$SOC1/10 # Convert to % if not yet done
spike$BD <- spike$BD1
spike$Clay <- spike$Clay1
spike$Sand <- spike$Sand1
spike[((spike$Med_depth >= 2.5)&(spike$Med_depth < 7.5)),]$SOC <- spike[((spike$Med_depth >= 2.5)&(spike$Med_depth < 7.5)),]$SOC2/10
spike[((spike$Med_depth >= 2.5)&(spike$Med_depth < 7.5)),]$BD <- spike[((spike$Med_depth >= 2.5)&(spike$Med_depth < 7.5)),]$BD2
spike[((spike$Med_depth >= 2.5)&(spike$Med_depth < 7.5)),]$Clay <- spike[((spike$Med_depth >= 2.5)&(spike$Med_depth < 7.5)),]$Clay2
spike[((spike$Med_depth >= 2.5)&(spike$Med_depth < 7.5)),]$Sand <- spike[((spike$Med_depth >= 2.5)&(spike$Med_depth < 7.5)),]$Sand2
spike[((spike$Med_depth >= 7.5)&(spike$Med_depth < 22.5)),]$SOC <- spike[((spike$Med_depth >= 7.5)&(spike$Med_depth < 22.5)),]$SOC3/10
spike[((spike$Med_depth >= 7.5)&(spike$Med_depth < 22.5)),]$BD <- spike[((spike$Med_depth >= 7.5)&(spike$Med_depth < 22.5)),]$BD3
spike[((spike$Med_depth >= 7.5)&(spike$Med_depth < 22.5)),]$Clay <- spike[((spike$Med_depth >= 7.5)&(spike$Med_depth < 22.5)),]$Clay3
spike[((spike$Med_depth >= 7.5)&(spike$Med_depth < 22.5)),]$Sand <- spike[((spike$Med_depth >= 7.5)&(spike$Med_depth < 22.5)),]$Sand3 
spike[((spike$Med_depth >= 22.5)&(spike$Med_depth < 37.5)),]$SOC <- spike[((spike$Med_depth >= 22.5)&(spike$Med_depth < 37.5)),]$SOC4/10
spike[((spike$Med_depth >= 22.5)&(spike$Med_depth < 37.5)),]$BD <- spike[((spike$Med_depth >= 22.5)&(spike$Med_depth < 37.5)),]$BD4
spike[((spike$Med_depth >= 22.5)&(spike$Med_depth < 37.5)),]$Clay <- spike[((spike$Med_depth >= 22.5)&(spike$Med_depth < 37.5)),]$Clay4
spike[((spike$Med_depth >= 22.5)&(spike$Med_depth < 37.5)),]$Sand <- spike[((spike$Med_depth >= 22.5)&(spike$Med_depth < 37.5)),]$Sand4 
spike[(spike$Med_depth >= 37.5),]$SOC <- spike[(spike$Med_depth >= 37.5),]$SOC5/10
spike[(spike$Med_depth >= 37.5),]$BD <- spike[(spike$Med_depth >= 37.5),]$BD5
spike[(spike$Med_depth >= 37.5),]$Clay <- spike[(spike$Med_depth >= 37.5),]$Clay5
spike[(spike$Med_depth >= 37.5),]$Sand <- spike[(spike$Med_depth >= 37.5),]$Sand5


###############################################################
# Reorganize and export final files of local samples for spiking
drop <- c("SOC1","SOC2","SOC3","SOC4","SOC5","BD1","BD2","BD3","BD4","BD5",
          "Clay1","Clay2","Clay3","Clay4","Clay5","Sand1","Sand2","Sand3","Sand4","Sand5")
spike = spike[,!(names(spike) %in% drop)]
spike <- spike[,c("Field","Y","X","Num","Cluster","Point","new_cluster","Med_depth","Sample_ID","EL","SL","AS","TWI","Tree",
                  "ppt","tmin","tmax","VPD","GPP","EVI","SOC","Clay","Sand","BD","AFGC","PFGC","TREE","SHR","LTR","BG","LULC","Measured_SOC")]
write.csv(spike,paste0("./",ROI,"_spike_",sample_num,".csv"))


###############################################################
# Reorganize and export final files of local samples held-out for validation
# Use the following codes to retain all local samples not used for spiking for validation
# Alternatively, retain the minimum validation dataset that corresponds to the highest number of local samples for spiking
spikelist <- unique(spike$Num)
vali <- all[!all$Num%in%spikelist,]
# Assign soil properties according to the sampling depth
vali$SOC <- vali$SOC1/10 # Convert to % if not yet done
vali$BD <- vali$BD1
vali$Clay <- vali$Clay1
vali$Sand <- vali$Sand1
vali[((vali$Med_depth >= 2.5)&(vali$Med_depth < 7.5)),]$SOC <- vali[((vali$Med_depth >= 2.5)&(vali$Med_depth < 7.5)),]$SOC2/10
vali[((vali$Med_depth >= 2.5)&(vali$Med_depth < 7.5)),]$BD <- vali[((vali$Med_depth >= 2.5)&(vali$Med_depth < 7.5)),]$BD2
vali[((vali$Med_depth >= 2.5)&(vali$Med_depth < 7.5)),]$Clay <- vali[((vali$Med_depth >= 2.5)&(vali$Med_depth < 7.5)),]$Clay2
vali[((vali$Med_depth >= 2.5)&(vali$Med_depth < 7.5)),]$Sand <- vali[((vali$Med_depth >= 2.5)&(vali$Med_depth < 7.5)),]$Sand2
vali[((vali$Med_depth >= 7.5)&(vali$Med_depth < 22.5)),]$SOC <- vali[((vali$Med_depth >= 7.5)&(vali$Med_depth < 22.5)),]$SOC3/10
vali[((vali$Med_depth >= 7.5)&(vali$Med_depth < 22.5)),]$BD <- vali[((vali$Med_depth >= 7.5)&(vali$Med_depth < 22.5)),]$BD3
vali[((vali$Med_depth >= 7.5)&(vali$Med_depth < 22.5)),]$Clay <- vali[((vali$Med_depth >= 7.5)&(vali$Med_depth < 22.5)),]$Clay3
vali[((vali$Med_depth >= 7.5)&(vali$Med_depth < 22.5)),]$Sand <- vali[((vali$Med_depth >= 7.5)&(vali$Med_depth < 22.5)),]$Sand3 
vali[((vali$Med_depth >= 22.5)&(vali$Med_depth < 37.5)),]$SOC <- vali[((vali$Med_depth >= 22.5)&(vali$Med_depth < 37.5)),]$SOC4/10
vali[((vali$Med_depth >= 22.5)&(vali$Med_depth < 37.5)),]$BD <- vali[((vali$Med_depth >= 22.5)&(vali$Med_depth < 37.5)),]$BD4
vali[((vali$Med_depth >= 22.5)&(vali$Med_depth < 37.5)),]$Clay <- vali[((vali$Med_depth >= 22.5)&(vali$Med_depth < 37.5)),]$Clay4
vali[((vali$Med_depth >= 22.5)&(vali$Med_depth < 37.5)),]$Sand <- vali[((vali$Med_depth >= 22.5)&(vali$Med_depth < 37.5)),]$Sand4 
vali[(vali$Med_depth >= 37.5),]$SOC <- vali[(vali$Med_depth >= 37.5),]$SOC5/10
vali[(vali$Med_depth >= 37.5),]$BD <- vali[(vali$Med_depth >= 37.5),]$BD5
vali[(vali$Med_depth >= 37.5),]$Clay <- vali[(vali$Med_depth >= 37.5),]$Clay5
vali[(vali$Med_depth >= 37.5),]$Sand <- vali[(vali$Med_depth >= 37.5),]$Sand5
drop <- c("SOC1","SOC2","SOC3","SOC4","SOC5","BD1","BD2","BD3","BD4","BD5",
          "Clay1","Clay2","Clay3","Clay4","Clay5","Sand1","Sand2","Sand3","Sand4","Sand5")
vali = vali[,!(names(vali) %in% drop)]
vali <- merge(vali,measured[c("Field","Num","Y","X","Med_depth","Measured_SOC")],by=c("Field","Med_depth","Num"))
vali <- vali[,c("Field","Y","X","Num","Cluster","Point","Med_depth","Sample_ID","EL","SL","AS","TWI","Tree",
                "ppt","tmin","tmax","VPD","GPP","EVI","SOC","Clay","Sand","BD","AFGC","PFGC","TREE","SHR","LTR","BG","LULC","Measured_SOC")]
write.csv(vali,paste0("./",ROI,"_vali_",spike_num,".csv"))
