library(readr)
library(ggplot2)
library(factoextra)
library(ggmap)
library(raster)
library(clhs)
library(scales)
library(broom)
library(rgdal)
setwd("set-pathway")


###############################################################
# Import dataset
dataset <- read_csv("./CSV-with-covariates-filename.csv")
dataset$aspect <- sin(pi*(dataset$aspect)/360) # Adjust aspect since 0 and 360 are same
subset <- dataset[,2:9] # Grad columns with covariates
df <- data.frame(scale(subset)) # Scale the data


###############################################################
# Carry out clustering analysis
ncomp = 5 # Define number of clusters
set.seed(1)
kmeans_basic <- kmeans(df,centers=ncomp,nstart = 1)
kmeans_basic_table <- data.frame(kmeans_basic$size, kmeans_basic$centers)
kmeans_basic_df <- data.frame(Cluster = kmeans_basic$cluster, df)
propdf <- data.frame(t(table(kmeans_basic_df$Cluster)))[,2:3] # Get proportions of each cluster
colnames(propdf) <- c("cluster","counts")
propdf <- propdf[order(propdf$counts),] # Show proportion of records within each cluster
propdf$prop <- propdf$counts/sum(propdf$counts)
propdf$perc <- round(propdf$prop*100)
# Visualize proportion of records within each cluster
ggplot(data=propdf, aes(x=counts, y=cluster, fill=cluster)) +
  geom_bar(stat="identity") + geom_text(aes(label=counts), hjust=+0.5, size=3.5)
# Determine number of clusters needed
factoextra::fviz_nbclust(df, kmeans, nstart=1, method = "wss") +
 geom_vline(xintercept = ncomp, linetype = 1)
# Carry out PCA and visualize results
p <- fviz_cluster(kmeans_basic, data = df, geom = c("point"),ellipse.type = "euclid")
pg <- ggplot_build(p)
p
geodf <- data.frame(dataset[,".geo"])
geodf$latlon <- gsub("[\\[\\]]", "", regmatches(geodf$.geo, gregexpr("\\[.*?\\]", geodf$.geo)))
geodf$latitude <- sapply(1:nrow(geodf),
                         function(i){
                           row <- geodf$latlon[i]
                           split <- strsplit(row,",")[[1]][2]
                           sub <- substr(split,1,nchar(split)-1)
                           return(as.numeric(sub))
                         })
geodf$longitude <- sapply(1:nrow(geodf),
                          function(i){
                            row <- geodf$latlon[i]
                            split <- strsplit(row,",")[[1]][1]
                            sub <- substr(split,2,nchar(split))
                            return(as.numeric(sub))
                          })
geodf$cluster <- kmeans_basic_df$Cluster
geo_kmeans_df <- data.frame(geodf,kmeans_basic_df)
bbox <- make_bbox(lat = latitude, lon = longitude, data = geodf)
bbox
big <- get_map(location = bbox, source = "google", maptype = "satellite")
# Visualize clusters
ggmap(big) + 
  geom_point(data = geodf, mapping = aes(x = longitude, y = latitude, color = as.factor(cluster)))
# Get raster file for determined clusters
write.csv(geodf[,c("latitude","longitude","cluster")], file="./define-filename1.csv", row.names=FALSE)
cluster.map <- read.table("./define-filename1.csv",header=T, sep=",")
res <- c(0.000898, 0.000898) # Define resolution
boundary <- raster(xmn= XXX, xmx= XXX, 
                   ymn= XXX, ymx= XXX, res=res,
                   crs="+proj=longlat +datum=WGS84") # Define XXX with boundary parameters
cluster_map <- rasterize(cluster.map[, c('longitude', 'latitude')], boundary, cluster.map[, 'cluster'], fun=mean, na.rm = TRUE)
plot(cluster_map)
writeRaster(cluster_map,"./define-filename2.tif",format="GTiff",overwrite=TRUE)


###############################################################
# Select sampling points
road_distance <- read_csv("./CSV-with-distance-to-roads-filename.csv") # Import distance to road file
road_distance$cluster <- geodf$cluster
finaldf <- data.frame(road_distance,df)
samples <- data.frame()
for(i in 1:ncomp){
  oneclust <- finaldf[finaldf$cluster==i,]
  ggmap(big) + 
    geom_point(data = oneclust, mapping = aes(x = longitude, y = latitude, color = HubDist, fill=HubDist))
  # Retain just samples close to road
  oneclustclose <- oneclust[oneclust$HubDist<=0.373,] # Up to 600 m away
  oneclustclose <- oneclustclose[oneclustclose$HubDist>=0.031,] # At least 50 m away
  ggmap(big) + 
    geom_point(data = oneclustclose, mapping = aes(x = longitude, y = latitude, color = HubDist, fill=HubDist))
  # CLHS Sampling of this group
  clusterperc <- propdf$perc[propdf$cluster==i]
  sample_indices <- clhs(oneclustclose[,6:13], size = clusterperc, progress = TRUE, iter = 500)
  thesesamples <- oneclustclose[sample_indices,]
  ggmap(big) + 
    geom_point(data = samples, mapping = aes(x = longitude, y = latitude, color = HubDist, fill=HubDist))
  samples <- rbind(samples,thesesamples)
}
write.csv(samples, file="./define-filename3.csv", row.names=FALSE)


###############################################################
# Map Samples
roads <- readOGR("./shapefile-with-distance-to-roads.shp")
roads <-spTransform(roads, CRS("+init=epsg:4326"))
roads <- tidy(roads)
ggmap(big)+
  geom_path(data = fortify(roads),aes(long, lat, group = group))+
  geom_point(data = samples, mapping = aes(x = longitude, y = latitude, 
            fill = as.factor(cluster)),color="Black",pch=21)
# Plot in PCA space
colors <- hue_pal()(ncomp)
pcadata <- data.frame(pg[["data"]][[1]])
pcasamples <- pcadata[row.names(samples),]
fviz_cluster(kmeans_basic, data = df, geom = c("point"),ellipse.type = "euclid")+
  geom_point(data=pcasamples, aes(x,y, fill=as.factor(pcasamples$group)), size=2, color="black",pch=21)+
  scale_fill_manual(values = colors)
