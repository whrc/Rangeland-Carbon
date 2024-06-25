###############################################################
# Method 1 - Use soilspec4gg library 
library(mongolite)
library(jsonify)
setwd("set-pathway")

# Import KSSL dataset
soilspec4gg.db = list(
  host = 'api.soilspectroscopy.org',
  name = 'soilspec4gg',
  user = 'soilspec4gg',
  pw = 'soilspec4gg'
)
soilspec4gg.db$url <- paste0(
  'mongodb://', soilspec4gg.db$user, ':', 
  soilspec4gg.db$pw, '@', 
  soilspec4gg.db$host, '/', 
  soilspec4gg.db$name, '?ssl=true'
)
soilspec4gg.init <- function() {
  print('Creating the access for mongodb collections.')
  soilspec4gg.db$collections <<- list(
    soillab = mongo(collection = 'soillab', url = soilspec4gg.db$url, verbose = TRUE),
    soilsite = mongo(collection = 'soilsite', url = soilspec4gg.db$url, verbose = TRUE),
    mir = mongo(collection = 'mir', url = soilspec4gg.db$url, verbose = TRUE),
    visnir = mongo(collection = 'visnir', url = soilspec4gg.db$url, verbose = TRUE)
  )  
}
soilspec4gg.samplesById <- function (ids) {
  print('Accessing mongodb collections.')
  query <- to_json( list( 
    id_layer_local_c = list("$in"=ids) 
  ))
  do.merge <- function(x, y) {
    if (nrow(y) > 0) {
      suppressWarnings(
        merge(x, y, by="id_layer_local_c", all=TRUE),
      )
    } else{
      x
    }
  }
  return(
    Reduce(f = do.merge,
           list(
             soilspec4gg.db$collections$soilsite$find(query = query),
             soilspec4gg.db$collections$soillab$find(query = query),
             soilspec4gg.db$collections$mir$find(query = query),
             soilspec4gg.db$collections$visnir$find(query = query)
           )
    )
  )
}
# Initialize the database to allow extraction
soilspec4gg.init()
# Import data file showing sample IDs needed for extraction, the IDs can be obtained from the KSSL database in the ACCESS format
KSSL.IDS <- read.table('./ID-filename.txt', header = TRUE)
# Choose the corresponding column for lay_ids if multiple columns are included
IDS_KSSL <- unlist(lapply(KSSL.IDS,as.numeric))
soilspec.samples.KSSL = soilspec4gg.samplesById(IDS_KSSL)
dim(soilspec.samples.KSSL)
KSSL.extracted <- soilspec.samples.KSSL[,c(XXX)] #Choose proper column names/numbers for extraction
write.csv(KSSL.extracted,"./define-filename1.csv")


###############################################################
# Method 2 - Use downloaded KSSL database and upload it to Google Cloud
# This approach requires the use of Ubuntu system in order to use the mdbtool
library("tidyverse")
library("Hmisc")
setwd("set-pathway")

# Import Access database
mdb.get("KSSL_DB_Access2000.mdb", tables = TRUE)
kssl.analyte <- mdb.get("KSSL_DB_Access2000.mdb", tables = c("analyte"))
head(kssl.analyte)
# Select soil properties
analyte.selection <- c("db_13b","db_fmstc") #Example of selecting bulk density measured with different methods
analyte.descriptions <- kssl.analyte %>%
  filter(analyte.abbrev %in% analyte.selection) %>%
  select(analyte.id, analyte.name, analyte.abbrev)
analyte.descriptions
analyte.id.selection <- kssl.analyte %>%
  filter(analyte.abbrev %in% analyte.selection) %>%
  pull(analyte.id)
kssl.layer_analyte <- mdb.get("KSSL_DB_Access2000.mdb",tables = c("layer_analyte"))
head(kssl.layer_analyte)
# Get soil data by layer ID
kssl.soildata <- kssl.layer_analyte %>%
  as_tibble() %>%
  filter(analyte.id %in% analyte.id.selection) %>%
  select(lay.id, analyte.id, calc.value) %>%
  arrange(lay.id, analyte.id) %>%
  left_join(analyte.descriptions, by = "analyte.id") %>%
  select(-analyte.id, -analyte.name) %>%
  mutate(lay.id = as.numeric(lay.id),
         analyte.abbrev = as.character(analyte.abbrev),
         calc.value = as.numeric(calc.value)) %>%
  group_by(lay.id, analyte.abbrev) %>%
  summarise(calc.value = mean(calc.value, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = "analyte.abbrev", values_from = "calc.value")
kssl.soildata
# Get soil data by sample id
kssl.sample <- mdb.get("KSSL_DB_Access2000.mdb", tables = c("sample"))
head(kssl.sample)
kssl.sample.ids <- kssl.sample %>%
  select(smp.id, lay.id) %>%
  rename(sample_id = smp.id) %>%
  mutate(sample_id = as.numeric(sample_id),
         lay.id = as.numeric(lay.id)) %>%
  arrange(lay.id)
head(kssl.sample.ids)
kssl.sample.ids <- kssl.sample.ids %>%
  group_by(lay.id) %>%
  summarise(count = n()) %>%
  filter(count == 1) %>%
  select(lay.id)
head(kssl.sample.ids)
kssl.sample.ids <- kssl.sample.ids %>%
  left_join({kssl.sample %>%
      select(smp.id, lay.id) %>%
      rename(sample_id = smp.id) %>%
      mutate(sample_id = as.numeric(sample_id),
             lay.id = as.numeric(lay.id))}, by = "lay.id")
kssl.sample.ids
kssl.sample.ids %>%
  distinct(sample_id) %>%
  nrow()
kssl.soildata
# Join files and export
kssl.soildata <- kssl.soildata %>%
  left_join(kssl.sample.ids, by = "lay.id") %>%
  select(-lay.id) %>%
  select(sample_id, everything()) %>%
  filter(!is.na(sample_id))
kssl.soildata
kssl.soildata <- kssl.soildata %>%
  rename(BD_clod = db_13b, BD_core = db_fmstc)
kssl.soildata
kssl.soildata$BD <- kssl.soildata$BD_core
kssl.soildata[which(is.na(kssl.soildata$BD)),]$BD <- kssl.soildata[which(is.na(kssl.soildata$BD)),]$BD_clod
kssl.soildata$BD_method <- "Core"
kssl.soildata$BD_method[kssl.soildata$BD == kssl.soildata$BD_clod] <- "Clod"
write.csv(kssl.soildata,"./define-filename2.csv")
