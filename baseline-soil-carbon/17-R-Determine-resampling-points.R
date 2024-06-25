library(raster)
library(ggplot2) 
library(factoextra)
library(ggmap)
library(clhs)
library(rgdal)
library(flexclust)
library(broom)
setwd("set-pathway")


# For fields not restricted by roads
###############################################################
# Define site and parameters
ROI = "specify-site-name"
sample_num = 15
ncomp = 7
exclude_area = "No" # Need to exclude inaccessible area in addition to road restrictions
res <- c(0.000269494585235856472, 0.000269494585235856472) # Define resolution
boundary <- raster(xmn= XXX, xmx= XXX, 
                   ymn= XXX, ymx= XXX, res=res,
                   crs="+proj=longlat +datum=WGS84") # Define XXX with boundary parameters
# Import datasets
predict.raster <- paste0("./",ROI,"_cov_all_RAP.tif")
x <- new("GDALReadOnlyDataset", predict.raster)
width <- dim(x)[2]
height <- dim(x)[1]
imagedata <- data.frame(getRasterTable(x))
# Define and process covariates
names(imagedata) <- c('x','y','EL','SL','AS','TWI',  
                      'ppt','tmin','tmax','VPD',
                      'BD2','BD5','Clay2','Clay5','Sand2','Sand5','SOC2','SOC5',
                      'Tree','GPP','EVI',
                      'AFGC','PFGC','TREE','SHR','LTR','BG','Depth')
# Scale covariates to match with KSSL data
imagedata$SOC2 <- imagedata$SOC2/10
imagedata$SOC5 <- imagedata$SOC5/10
imagedata$ppt <- imagedata$ppt/20
imagedata$tmin <- imagedata$tmin/20
imagedata$tmax <- imagedata$tmax/20
imagedata$VPD <- imagedata$VPD/20
imagedata$GPP <- imagedata$GPP/20
imagedata$AS <- sin(pi*(imagedata$AS)/360)
imagedata <- na.omit(imagedata)
# Match with pixel ID
pred.ID <- tibble::rowid_to_column(imagedata, "ID")
imagedata <- cbind(pred.ID[,1], imagedata)
colnames(imagedata)[1] <- c("pred.ID")
subset <- imagedata[,4:(dim(imagedata)[2]-1)] # Customize column names or numbers to exclude depth, coordinates, and pred ID
df <- data.frame(scale(subset))
# Carry out clustering analysis
factoextra::fviz_nbclust(df, kmeans, nstart=1, method = "wss") + geom_vline(xintercept = ncomp, linetype = 1) 
set.seed(1)
kmeans_basic <- kmeans(df,centers=ncomp,nstart = 10)
kmeans_basic_table <- data.frame(kmeans_basic$size, kmeans_basic$centers) 
kmeans_basic_df <- data.frame(Cluster = kmeans_basic$cluster, df)
propdf <- data.frame(t(table(kmeans_basic_df$Cluster)))[,2:3]
colnames(propdf) <- c("cluster","counts")
propdf <- propdf[order(propdf$counts),] 
propdf$prop <- propdf$counts/sum(propdf$counts)
propdf$perc <- round(propdf$prop*100)
ggplot(data=propdf, aes(x=counts, y=cluster, fill=cluster)) +
  geom_bar(stat="identity") + geom_text(aes(label=counts), hjust=+0.5, size=3.5)
# Visualize cluster distribution in space
geodf <- data.frame(imagedata,kmeans_basic_df$Cluster)
colnames(geodf)[dim(geodf)[2]] = c("cluster")
bbox <- make_bbox(lat = y, lon = x, data = geodf)
bbox
big <- get_map(location = bbox, source = "google", maptype = "satellite")
ggmap(big) + 
  geom_point(data = geodf, mapping = aes(x, y, color = as.factor(cluster)))
# cLHS sampling based on cluster proportions
if (exclude_area == "Yes") {
  whole_list <- c(1:ncomp)
  exclude_cluster_list <- c(1,5,6,8) # Read the map and add cluster number that should be excluded
  retain_cluster_list <- whole_list[!whole_list%in%exclude_cluster_list]
  total_count <- sum(propdf[which(propdf$cluster %in% retain_cluster_list),"counts"])
  propdf$percentage <- propdf$counts/total_count
  reorg_propdf <- propdf[which(propdf$cluster %in% retain_cluster_list),]
} 
samples <- data.frame()
if (exclude_area == "No") {
  for(i in 1:ncomp){
    oneclust <- geodf[geodf$cluster==i,]
    clusterperc <- round(propdf$perc[propdf$cluster==i]*sample_num/100) 
    if (clusterperc >= 0.04) {
      sample_indices <- clhs(oneclust[,4:dim(geodf)[2]-2], size = clusterperc, progress = TRUE, iter = 500)
      thesesamples <- oneclust[sample_indices,]
      samples <- rbind(samples,thesesamples)}}
} else {
  for(i in 1:length(retain_cluster_list)[1]){
    oneclust <- geodf[geodf$cluster==retain_cluster_list[i],]
    clusterperc <- round(reorg_propdf$percentage[reorg_propdf$cluster==retain_cluster_list[i]]*sample_num) 
    if (clusterperc >= 0.04) {
      sample_indices <- clhs(oneclust[,4:dim(geodf)[2]-2], size = clusterperc, progress = TRUE, iter = 500)
      thesesamples <- oneclust[sample_indices,]
      samples <- rbind(samples,thesesamples)}}
}
# Visualize sampling points in space
AOI <- readOGR(paste0("./",ROI,".shp"))
AOI <-spTransform(AOI, CRS("+init=epsg:4326"))
AOI <- tidy(AOI)
samples_geodf <- geodf[row.names(samples),]
ggmap(big)+
  geom_path(data = fortify(AOI),aes(long, lat, group = group),size = 2)+
  geom_point(data = samples_geodf, mapping = aes(x, y, fill = as.factor(cluster)),color="Black",pch=21,size=4)
# Export sampling points
write.csv(samples, paste0("./",ROI,"_",sample_num,"_samples.csv"))
# Export and visualize cluster map
write.csv(geodf[,c("x","y","cluster")], file=paste0("./",ROI,"_clusters.csv"), row.names=FALSE)
cluster.map <- read.table(paste0("./",ROI,"_clusters.csv"),header=T, sep=",")
cluster_map <- rasterize(cluster.map[, c('x', 'y')], boundary, cluster.map[, 'cluster'], fun=mean, na.rm = TRUE)
plot(cluster_map)
writeRaster(cluster_map, paste0("./",ROI,"_clusters.tif"),format="GTiff",overwrite=TRUE)


# For fields restricted by roads
###############################################################
# Define site and parameters
ROI = "specify-site-name"
AOI <- readOGR(paste0("./",ROI,".shp"))  
AOI <-spTransform(AOI, CRS("+init=epsg:4326"))
AOI <- tidy(AOI)
roads <- readOGR(paste0("./",ROI,"_roads.shp"))  
roads <-spTransform(roads, CRS("+init=epsg:4326"))
roads <- tidy(roads)
sample_num = 15
ncomp = 7
exclude_area = "No" # Need to exclude inaccessible area in addition to road restrictions
res <- c(0.000269494585235856472, 0.000269494585235856472) # Define resolution
boundary <- raster(xmn= XXX, xmx= XXX, 
                   ymn= XXX, ymx= XXX, res=res,
                   crs="+proj=longlat +datum=WGS84") # Define XXX with boundary parameters
predict.raster <- paste0("./",ROI,"_cov_all_RAP.tif")
x <- new("GDALReadOnlyDataset", predict.raster)
width <- dim(x)[2]
height <- dim(x)[1]
imagedata <- data.frame(getRasterTable(x))
# Define and process covariates
names(imagedata) <- c('x','y','EL','SL','AS','TWI',  
                      'ppt','tmin','tmax','VPD',
                      'BD2','BD5','Clay2','Clay5','Sand2','Sand5','SOC2','SOC5',
                      'Tree','GPP','EVI',
                      'AFGC','PFGC','TREE','SHR','LTR','BG','Depth')
# Scale covariates to match with KSSL data
imagedata$SOC2 <- imagedata$SOC2/10
imagedata$SOC5 <- imagedata$SOC5/10
imagedata$ppt <- imagedata$ppt/20
imagedata$tmin <- imagedata$tmin/20
imagedata$tmax <- imagedata$tmax/20
imagedata$VPD <- imagedata$VPD/20
imagedata$GPP <- imagedata$GPP/20
imagedata$AS <- sin(pi*(imagedata$AS)/360)
imagedata <- na.omit(imagedata)
# Match with pixel ID
pred.ID <- tibble::rowid_to_column(imagedata, "ID")
imagedata <- cbind(pred.ID[,1], imagedata)
colnames(imagedata)[1] <- c("pred.ID")
imagedata.point <- imagedata[,c(1:3)] 
write.csv(imagedata.point,paste0("./",ROI,"_points.csv"),row.names = FALSE) # Export coords file and then calculate road distance in GIS 
imagedata.point.dist <- read.table(paste0("./",ROI,"_points_w_dist.csv"), comment.char ="", quote = "\"", header = T, sep = ",")
imagedata.point.dist <- imagedata.point.dist[,c("farmid","x","y","NEAR_DIST")]
imagedata.comb <- merge(imagedata,imagedata.point.dist,by=c("farmid"))
imagedata.comb <- subset(imagedata.comb, select = -c(x.y,y.y))
colnames(imagedata.comb)[2:3] <- c("x","y")
subset <- imagedata.comb[,4:(dim(imagedata.comb)[2]-2)] # Customize column names or numbers to exclude depth, coordinates, distance to road, and pred ID
df <- data.frame(scale(subset))
# Carry out clustering analysis
factoextra::fviz_nbclust(df, kmeans, nstart=1, method = "wss") + geom_vline(xintercept = ncomp, linetype = 1) 
set.seed(1)
kmeans_basic <- kmeans(df,centers=ncomp,nstart = 10)
kmeans_basic_table <- data.frame(kmeans_basic$size, kmeans_basic$centers) 
kmeans_basic_df <- data.frame(Cluster = kmeans_basic$cluster, df)
propdf <- data.frame(t(table(kmeans_basic_df$Cluster)))[,2:3]
colnames(propdf) <- c("cluster","counts")
propdf <- propdf[order(propdf$counts),] 
propdf$prop <- propdf$counts/sum(propdf$counts)
propdf$perc <- round(propdf$prop*100)
ggplot(data=propdf, aes(x=counts, y=cluster, fill=cluster)) +
  geom_bar(stat="identity") + geom_text(aes(label=counts), hjust=+0.5, size=3.5)
# Visualize cluster distribution in space
geodf <- data.frame(imagedata.comb,kmeans_basic_df$Cluster)
colnames(geodf)[dim(geodf)[2]] = c("cluster")
bbox <- make_bbox(lat = y, lon = x, data = geodf)
bbox
big <- get_map(location = bbox, source = "google", maptype = "satellite")
ggmap(big) + 
  geom_point(data = geodf, mapping = aes(x, y, color = as.factor(cluster)))
# Sample selection
if (exclude_area == "Yes") {
  whole_list <- c(1:ncomp)
  exclude_cluster_list <- c(1,5,6,8) # Read the map and add cluster number that should be excluded
  retain_cluster_list <- whole_list[!whole_list%in%exclude_cluster_list]
  total_count <- sum(propdf[which(propdf$cluster %in% retain_cluster_list),"counts"])
  propdf$percentage <- propdf$counts/total_count
  reorg_propdf <- propdf[which(propdf$cluster %in% retain_cluster_list),]
} 
samples <- data.frame()
if (exclude_area == "No") {
  for(i in 1:ncomp){
    oneclust <- geodf[geodf$cluster==i,]
    oneclustclose <- oneclust[oneclust$NEAR_DIST<=600,] 
    oneclustclose <- oneclustclose[oneclustclose$NEAR_DIST>=50,]
    clusterperc <- round(propdf$perc[propdf$cluster==i]*sample_num/100) #total % of all cluster is 100 so if choosing 100 samples no scaling is needed
    if (clusterperc >= 0.04) {
      sample_indices <- clhs(oneclustclose[,4:dim(geodf)[2]-2], size = clusterperc, progress = TRUE, iter = 500)
      thesesamples <- oneclustclose[sample_indices,]
      samples <- rbind(samples,thesesamples)}}
} else {
  for(i in 1:length(retain_cluster_list)[1]){
    oneclust <- geodf[geodf$cluster==retain_cluster_list[i],]
    oneclustclose <- oneclust[oneclust$NEAR_DIST<=600,] 
    oneclustclose <- oneclustclose[oneclustclose$NEAR_DIST>=50,]
    clusterperc <- round(reorg_propdf$percentage[reorg_propdf$cluster==retain_cluster_list[i]]*sample_num) 
    if (clusterperc >= 0.04) {
      sample_indices <- clhs(oneclustclose[,4:dim(geodf)[2]-2], size = clusterperc, progress = TRUE, iter = 500)
      thesesamples <- oneclustclose[sample_indices,]
      samples <- rbind(samples,thesesamples)}}
}
# Plot sampling points on the farm boundary and roads
samples_geodf <- geodf[row.names(samples),]
ggmap(big)+
  geom_path(data = fortify(AOI),aes(long, lat, group = group),size = 2)+
  geom_path(data = fortify(roads),aes(long, lat, group = group),size = 1,color="red")+
  geom_point(data = samples_geodf, mapping = aes(x, y, fill = as.factor(cluster)),color="Black",pch=21,size=4)
# Export raster for clusters
write.csv(geodf[,c("x","y","cluster")], file=paste0("./",ROI,"_clusters.csv"), row.names=FALSE)
cluster.map <- read.table(paste0("./",ROI,"_clusters.csv"),header=T, sep=",")
cluster_map <- rasterize(cluster.map[, c('x', 'y')], boundary, cluster.map[, 'cluster'], fun=mean, na.rm = TRUE)
plot(cluster_map)
writeRaster(cluster_map, paste0("./",ROI,"_clusters.tif"),format="GTiff",overwrite=TRUE)
# Export sampling points
write.csv(samples, paste0("./",ROI,"_",sample_num,"_samples.csv"))
