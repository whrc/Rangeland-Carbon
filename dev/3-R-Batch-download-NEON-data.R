# Define site and temporal coverage
site <- "sitename"
DateBgn <- as.Date("2018-08-01") # Define start date
DateEnd <- as.Date("2022-04-01") # Define end date

#Define files and pathway for downloading
TypeFile <- c("atm/cdeps", "surf_files", "eval")[1] # Meteorological data
TypeFile2 <- c("atm/cdeps", "surf_files", "eval")[3] # NEE/GPP/RECO
versData <- c("v1","v2")[2] # Define version
versDataSurf <- "v1"
DirDnldBase <- "C:/Users/username/path-to-the-NEON-folder"
dnldDate <- paste0("v",format(Sys.time(), "%Y%m%d"))
DirDnld <- paste0(DirDnldBase,"/",dnldDate,"/",site)

#Data downloading algorithms
if(!dir.exists(DirDnld)) dir.create(DirDnld, recursive = TRUE)
DateSeq <- seq.Date(from = DateBgn,to = DateEnd, by = "month")
PrdWndwDnld <- strftime(DateSeq, format = "%Y-%m")
if(TypeFile == "surf_files"){
  urlDnld <- paste0("https://storage.googleapis.com/neon-ncar/NEON/",TypeFile,"/",versDataSurf,"/",site,"_surfaceData.csv")
  fileDnld <- paste0(DirDnld,"/",site,"_",versData,"_surfaceData.csv")
} else
{
  TypeFileDnld <- unlist(strsplit(TypeFile, "/|_"))[1]
  urlDnld <- paste0("https://storage.googleapis.com/neon-ncar/NEON/",TypeFile,"/",versData,"/",site,"/",site,"_",ifelse(TypeFile == "atm/cdeps",TypeFileDnld,TypeFile),"_",PrdWndwDnld,".nc")
  fileDnld <-  paste0(DirDnld,"/",site,"_",TypeFileDnld,"_",versData,"_",PrdWndwDnld,".nc") 
}
if(TypeFile2 == "surf_files"){
  urlDnld2 <- paste0("https://storage.googleapis.com/neon-ncar/NEON/",TypeFile2,"/",versDataSurf,"/",site,"_surfaceData.csv")
  fileDnld2 <- paste0(DirDnld,"/",site,"_",versData,"_surfaceData.csv")
} else
{
  TypeFileDnld2 <- unlist(strsplit(TypeFile2, "/|_"))[1]
  urlDnld2 <- paste0("https://storage.googleapis.com/neon-ncar/NEON/",TypeFile2,"_files/",versData,"/",site,"/",site,"_",ifelse(TypeFile2 == "atm/cdeps",TypeFileDnld2,TypeFile2),"_",PrdWndwDnld,".nc")
  fileDnld2 <-  paste0(DirDnld,"/",site,"_",TypeFileDnld2,"_",versData,"_",PrdWndwDnld,".nc") 
} 
# Download files
sapply(seq_along(urlDnld), function(x){
  download.file(url = urlDnld[x], destfile = fileDnld[x], mode = 'wb')
})
sapply(seq_along(urlDnld2), function(x){
  download.file(url = urlDnld2[x], destfile = fileDnld2[x], mode = 'wb')
})