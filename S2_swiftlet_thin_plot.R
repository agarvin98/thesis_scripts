## in this script
# 1. create presence and absence subset csv files
# 2. spatially thin each and re-concatenate (if needed)
# 3. primary plot 

# 1. create presence and absence subset csv files
# load libraries
library(remotes)
library(sf)
library(rnaturalearth)
library(auk)
library(lubridate)
library(gridExtra)
library(tidyverse)
library(dplyr)

##establish directories (double check changes here)
localDir = "~/Desktop/swiftlet_data/ebird/zf_data/zf_files/pluglo"
setwd(localDir)
dir.create("thinned_pres", showWarnings = FALSE)
outputDir <- file.path(localDir, "thinned_pres")
outputname <- c("pluglo")
## add in here replace name for more reproducibility?

## start here if need to subset pres/abs --------------------------------------
#load zf ebd
full_zf <- read.csv("gerswi_zf.csv")
##checkout summary and presense/absence data
head(full_zf,20)
count(full_zf, species_observed==TRUE) # take note of dimensions
# create a variable for presence and absences 
# presence
presence_only <- full_zf[full_zf$species_observed == "TRUE",]
head(presence_only,20) #checking
#create presence only file
write.csv(presence_only, "barswa_presence.csv") #change 
presence_only <- read.csv("~/Desktop/swiftlet_data/ebird/zf_data/barswa_presence.csv") #change
dim(presence_only)
# absence
absence_only <- full_zf[full_zf$species_observed == "FALSE",]
head(absence_only,20) #checking
#create absence only file
write.csv(absence_only, "barswa_absence.csv") #change
absence_only <- read.csv("~/Desktop/swiftlet_data/ebird/zf_data/barswa_absence.csv") #change
head(absence_only)
dim(absence_only)


## Start here if already have pres/abs datasets --------------------------------
# 2. spatially thin to address sampling bias and re-concatenate (if needed)
##applying spthin to presence data set and absence data set
library(spThin)
#spatially thin presence files
data <-read.csv("pluglo_pres.csv") ## read in your absence data
head(data)
dim(data)
#thin
thin_data <-
  thin(loc.data = data,
       lat.col = "latitude", long.col = "longitude",
       spec.col = "scientific_name",
       thin.par = 5, reps = 30,
       locs.thinned.list.return = TRUE, max.files = 1,
       out.dir = outputDir, #change
       out.base = outputname, write.files = TRUE, #change
       write.log.file = FALSE)

## check single file
thinned <- read.csv(file.path(outputDir, "pluglo_thin1.csv"))
head(thinned) #checking
dim(thinned) #checking
### if vector was not exhausted can move to plotting below

#loop to partition off thinning if vector gets exhausted -----------------------
#so we dont have to manually filter by 2 degrees latitude
# :) yay 
data<-arrange(data, latitude) ## sort your data by latitude
min(data$latitude) # checking
max(data$latitude) # checking
# create a vector of your latitude range
sections<-data$latitude[seq(1, length(data$latitude), 500)] #you can mess around with this number to change how fine its subsetting your latitude
## add a new element at the end with a number larger than your highest latitude- this will make sure you don't miss data in your loop
sections[length(sections)+1]<-max(data$lat)
sections # check how many elements are in your vector
#output<-vector() #if you are testing out your code, you can initialize this vector
for (i in 1:length(sections)){
  test_data<-filter(data, latitude >= sections[i] & latitude <sections[i+1]) # filters your data based on your latitude ranges
  #output[i]<-max(test_data$latitude) #this is just for testing that your code works
  thin(loc.data = test_data,
       lat.col = "latitude", long.col = "longitude",
       spec.col = "scientific_name",
       thin.par = 5, reps = 30,
       locs.thinned.list.return = TRUE, max.files = 1,
       out.dir = "~/Desktop/swiftlet_data/ebird/zf_data/con",
       out.base = paste0("thinned_barswa",i), write.files = TRUE,
       write.log.file = FALSE)
}
#recombine data back together
#combine thinned datasets for new master thinned dataset
setwd("~/Desktop/swiftlet_data/ebird/zf_data/con")
list_thinned <- list.files("~/Desktop/swiftlet_data/ebird/zf_data/con")
thinned <- lapply(list_thinned, read.csv, header = TRUE)
thinned_swi<-do.call(rbind, thinned)
nrow(thinned_swi)
head(thinned_swi)
write.csv(thinned_swi, "pres_barswa_thinned.csv", row.names = FALSE)
#read concatendated csv back in
swi_rethin <- read.csv("abs_monswi_thin.csv")
#rethin if you used loop to be thorough
rethin <-
  thin(loc.data = swi_rethin,
       lat.col = "latitude", long.col = "longitude",
       spec.col = "scientific_name",
       thin.par = 5, reps = 30,
       locs.thinned.list.return = TRUE, max.files = 1,
       out.dir = "~/Desktop/swiftlet_data/ebird/zf_data/absence",
       out.base = "monswi_thin", write.files = TRUE,
       write.log.file = FALSE)


# 3. primary plots on street map data first ------------------------------------
#install.packages("ecospat")
#install.packages("ape")
#install.package = 'rgdal'
library("rgdal")
library(ape)
library(ecospat)
library(raster)
library(RColorBrewer)
library(maptools)
library(tidyverse)
library(rasterVis)
data("wrld_simpl")
# set options
options(stringsAsFactors = F)         # no automatic data transformation
options("scipen" = 100, "digits" = 4) # suppress math annotation
op <- options(gvis.plot.tag='chart')  # set gViz options
# load package
library(OpenStreetMap)
library(DT)
library(mapproj)
library(RgoogleMaps)
library(scales)
library(rworldmap)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)
library(ggspatial)
library(maptools)
library(leaflet)
library(tmap)
library(here)
library(scales)
library(flextable)
# activate klippy for copy-to-clipboard button
klippy::klippy()

##different map type options
opt <- c("osm", "osm-bw","maptoolkit-topo", "waze", "bing", "stamen-toner", "stamen-terrain", "stamen-watercolor", "osm-german", "osm-wanderreitkarte", "mapbox")
opt2 <- c("esri", "esri-topo", "nps", "apple-iphoto", "skobbler", "hillshade", "opencyclemap", "osm-transport", "osm-public-transport", "osm-bbike", "osm-bbike-german")
opt <- data.frame(opt, opt2) %>%
  dplyr::rename(options = colnames(.)[1],
                more = colnames(.)[2])
opt %>%
  as.data.frame() %>%
  head(15) %>%
  flextable() %>%
  flextable::set_table_properties(width = .5, layout = "autofit") %>%
  flextable::theme_zebra() %>%
  flextable::fontsize(size = 12) %>%
  flextable::fontsize(size = 12, part = "header") %>%
  flextable::align_text_col(align = "center") %>%
  flextable::set_caption(caption = "Map display options.")  %>%
  flextable::border_outer()

# load data
setwd("~/Desktop/swiftlet_data/ebird/zf_data/zf_files/pluglo/")
localDir = "~/Desktop/swiftlet_data/ebird/zf_data/zf_files/pluglo"
RAWpresence <- read.csv(file.path(localDir, "pluglo_pres.csv"))
head(RAWpresence)
THINpresence <- read.csv(file.path(localDir, "thinned_pres/pluglo_thin1.csv"))
head(THINpresence)
### COMPARE against pre thinned data
### just run thinned points after for pres map
##plot on map to see where it is clustering
#south east asia
SEA <- getMap(resolution = "low")
# upper left (18,89) lower right (-21,170)
# plot data on world map
plot(SEA, xlim = c(89, 170), ylim = c(-21, 18), 
     asp = 1, bg = "azure1", border = "darkgrey", 
     col = "wheat2", fill = T)
# add points (PRE-THINNED)
points(RAWpresence$longitude, RAWpresence$latitude,
       # define colors as transparent
       col = rgb(red = 0, green = 0, blue = 1, alpha = 0.3),
       # define text and symbol size
       cex = 1.5, pch = 20)
# add points (THINNED)
points(THINpresence$longitude, THINpresence$latitude,
       # define colors as transparent
       col = rgb(red = 1, green = 0, blue = 0, alpha = 0.3),
       # define text and symbol size
       cex = 1.5, pch = 20)

##save files for future reference

###make sure species_observed column is back on
head(final_combdata)
final_combdata$species_observed <- TRUE ##adding back on species observed column
head(final_combdata)
