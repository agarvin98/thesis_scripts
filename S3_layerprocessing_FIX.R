#swiftlet single species distribution model layer processing
## mask raster 
## crop masked raster
#load packages
detectCores()
library(tidyverse)
library(dismo)
library(dplyr)
library(reshape2)
require(rgdal)
library(rgeos)
require(rgeos)
library(maptools)
library(sp)
library(rasterVis)
library(raster)
library(grid)
library(tmaptools)
library(sf)
library(doParallel)
library(foreach)
data("wrld_simpl")

#check files that were unzipped
hyde_asc <-list.files("~/Desktop/swiftlet_data/enviro_data/CHELSA/clim_bio/", pattern=".tif$", full.names=TRUE)
hyde_asc #should be 561 
#cropping
setwd("~/Desktop/swiftlet_data/enviro_data/CHELSA/clim_bio/")
hyde_layers <-list.files("~/Desktop/swiftlet_data/enviro_data/CHELSA/clim_bio/", full.names=TRUE)
output_dir<-"~/Desktop/swiftlet_data/enviro_data/cropped_chelsa/"

cores<- 6 #can set to more cores, doParallel means run once per core
cl<- makeCluster(cores)
registerDoParallel(cl)

i=1
foreach(i=1:NROW(hyde_layers)) %dopar% {
  library(raster)
  #read in raster
  tempraster <- raster(hyde_layers[i]) #change this to 1 and run each line until you can check the plot before running entire loop
  # ## center over Pacific
  # x1 <- crop(tempraster, extent(-180, 0, -60, 90))
  # x2 <- crop(tempraster, extent(0, 180, -60, 90))   
  # extent(x1) <- c(180, 360, -60, 90)
  # mergeraster <- merge(x1, x2)
  # names(mergeraster)<-names(tempraster)
  # # crop mergeraster numbers subject to change
  sample.range <-extent(82,162,-19, 30) #CHANGE THIS TO EXTENT OF FOCUS
  croppedraster <-crop(tempraster, sample.range, snap='near')
  #plot(croppedraster) #use when first checking extent is properly cropped
  #write cropped raster
  writeRaster(croppedraster, filename=file.path(output_dir, names(croppedraster)), format="ascii", bylayer=TRUE, overwrite=FALSE, suffix=names(hyde_layers))
}
stopCluster(cl)

#cropping
setwd("~/Desktop/swiftlet_data/enviro_data/land_use/")
hyde_layers <-list.files("~/Desktop/swiftlet_data/enviro_data/land_use/", full.names=TRUE)
output_dir<-"~/Desktop/swiftlet_data/enviro_data/cropped_landuse/"

cores<- 6 #can set to more cores, doParallel means run once per core
cl<- makeCluster(cores)
registerDoParallel(cl)

i=1
foreach(i=1:NROW(hyde_layers)) %dopar% {
  library(raster)
  #read in raster
  tempraster <- raster(hyde_layers[i]) #change this to 1 and run each line until you can check the plot before running entire loop
  # ## center over Pacific
  # x1 <- crop(tempraster, extent(-180, 0, -60, 90))
  # x2 <- crop(tempraster, extent(0, 180, -60, 90))   
  # extent(x1) <- c(180, 360, -60, 90)
  # mergeraster <- merge(x1, x2)
  # names(mergeraster)<-names(tempraster)
  # # crop mergeraster numbers subject to change
  sample.range <-extent(82,162,-19, 30) #CHANGE THIS TO EXTENT OF FOCUS
  croppedraster <-crop(tempraster, sample.range, snap='near')
  #plot(croppedraster) #use when first checking extent is properly cropped
  #downscale resolution- comment this out for hyde layers as it will slow the processing down
  #maskedcropped<-resample(croppedraster, resample_layer, method="bilinear")
  #write cropped raster
  writeRaster(croppedraster, filename=file.path(output_dir, names(croppedraster)), format="ascii", bylayer=TRUE, overwrite=FALSE, suffix=names(hyde_layers))
}
stopCluster(cl)

########### creating mask for the raster layers
#crops the environmental predictors tighter to your sighting records
###load unmasked layers
#cropping
setwd("~/Desktop/swiftlet_data/enviro_data/cropped_chelsa/")
hyde_layers <-list.files("~/Desktop/swiftlet_data/enviro_data/cropped_chelsa/", full.names=TRUE)
output_dir<-"~/Desktop/swiftlet_data/enviro_data/masked_chelsa/"


### speed up by cropping, masking and writing re-centered rasters in loop and in parallel
#can resample layers in loop too but commented out to save time
cores<- 6 #can set to more cores (don't want to use all of computer's cores or it will crash)
cl<- makeCluster(cores)
registerDoParallel(cl)

i=1
foreach(i=1:NROW(hyde_layers)) %dopar% {
  library(raster)
  #read in raster
  tempraster <- raster(hyde_layers[i])
  ## center over Pacific
  #x1 <- crop(tempraster, extent(-180, 0, -60, 90))
  #x2 <- crop(tempraster, extent(0, 180, -60, 90))   
  #extent(x1) <- c(180, 360, -60, 90)
  #mergeraster <- merge(x1, x2)
  #names(mergeraster)<-names(tempraster)
  # mask using polygon created above outside of loop
  maskedraster <- mask(x = tempraster, mask = wrld_simpl)
  #plot(tempraster)
  #plot(maskedraster)
  #crop masked raster to study extent
  #sample.range<-extent(82,162,-19, 30)
  #maskedcropped<-crop(maskedraster, sample.range, snap='near')
  # downscale resolution- comment this out for hyde layers as it will slow the processing down
  #maskedcropped<-resample(maskedcropped, resample_layer, method="bilinear")
  #write raster
  writeRaster(maskedraster, filename=file.path(output_dir, names(maskedraster)), format="ascii", bylayer=TRUE, overwrite=TRUE)
}
stopCluster(cl)


########### downscaling raster layers
###load low res model layers
setwd("~/Desktop/swiftlet_data/enviro_data/")
resample_layer <- raster("~/Desktop/swiftlet_data/enviro_data/cropland2017AD.asc")
plot(resample_layer) #check
#crop to extent of layers youre resampling
sample.range<-extent(82,162,-19, 30)
resample_layer<-crop(resample_layer, sample.range, snap='near')
#read in layers to be rescaled
layers <-list.files("~/Desktop/swiftlet_data/enviro_data/present_day/", full.names=TRUE)
output_dir<-"~/Desktop/swiftlet_data/enviro_data/DS_present_day/"


### speed up by cropping, masking and writing re-centered rasters in loop and in parallel
#can resample layers in loop too but commented out to save time
cores<- 6 #can set to more cores (don't want to use all of computer's cores or it will crash)
cl<- makeCluster(cores)
registerDoParallel(cl)

i=1
foreach(i=1:NROW(layers)) %dopar% {
  library(raster)
  #read in raster
  tempraster <- raster(layers[4])
  # downscale resolution- comment this out for hyde layers as it will slow the processing down
  scaled<-resample(tempraster, resample_layer, method="bilinear")
  plot(tempraster) #check first time through and check res
  plot(scaled) #check first time through and check res
  # #write raster
  writeRaster(scaled, filename=file.path(output_dir, names(scaled)), format="ascii", bylayer=TRUE, overwrite=TRUE)
}
stopCluster(cl)


###BUFFERING ENVIRO DATA
##read in environmental variables
#cropping
setwd("~/Desktop/swiftlet_data/enviro_data/present_day/")
envs_files <-list.files("~/Desktop/swiftlet_data/enviro_data/present_day/", full.names=TRUE)
output_dir<-"~/Desktop/swiftlet_data/enviro_data/buffer_presentday/"
##stack environmental variables
envs <- raster::stack(envs_files)
##read in point data of all swiftlets
occs <- read.csv("~/Desktop/swiftlet_data/ebird/zf_data/zf_files/unique/new_swiswa.csv")
table(occs$scientific_name)
table(occs$species_observed)
occs<- occs[occs$scientific_name !="Hirundo rustica",]
occs<- occs[occs$scientific_name !="Hirundo tahitica",]
occs<- occs[occs$species_observed!="FALSE",]
head(occs)
# We'll now experiment with a different spatial R package called sf (simple features).
# Let's make our occs into a sf object -- as the coordinate reference system (crs) for these 
# points is WGS84, a geographic crs (lat/lon) and the same as our envs rasters, we specify it 
# as the RasterStack's crs.
occs.sf <- sf::st_as_sf(occs, coords = c("longitude","latitude"), crs = raster::crs(envs))

# Now, we project our point data to an equal-area projection, which converts our 
# degrees to meters, which is ideal for buffering (the next step). 
# We use the typical Eckert IV projection.
eckertIV <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
occs.sf <- sf::st_transform(occs.sf, crs = eckertIV)

# Buffer all occurrences by 500 km, union the polygons together 
# (for visualization), and convert back to a form that the raster package 
# can use. Finally, we reproject the buffers back to WGS84 (lat/lon).
# We choose 500 km here to avoid sampling the Caribbean islands.
occs.buf <- sf::st_buffer(occs.sf, dist = 500000) %>% 
  sf::st_union() %>% 
  sf::st_sf() %>%
  sf::st_transform(crs = raster::crs(envs))
plot(envs[[1]], main = names(envs)[1])
points(occs)
# To add sf objects to a plot, use add = TRUE
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)
# Crop environmental rasters to match the study extent
envs.bg <- raster::crop(envs[[1]], occs.buf)
# Next, mask the rasters to the shape of the buffers
envs.bg <- raster::mask(envs.bg, occs.buf)
plot(envs.bg[[1]], main = names(envs)[1])
points(occs)
plot(occs.buf, border = "blue", lwd = 3, add = TRUE)

# Assuming 'envs' is a list of rasters and 'occs.buf' is a spatial polygon
for (i in seq_along(envs)) {
  # Crop environmental rasters to match the study extent
  envs.bg <- raster::crop(envs[[i]], occs.buf)
  # Next, mask the rasters to the shape of the buffers
  envs.bg <- raster::mask(envs.bg, occs.buf)
  # Plot the cropped and masked raster
  #plot(envs.bg[[1]], main = names(envs)[i])
  # Save the cropped and masked raster to the output directory
  output_name <- paste0(output_dir, "/", names(envs)[i], "_cropped_masked.asc")
  raster::writeRaster(envs.bg, filename = output_name, format = "ascii")
}

###### overlay raster grid for layers to be put into Arcgis 




