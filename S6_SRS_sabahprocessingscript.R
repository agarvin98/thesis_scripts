### consolidating enmeval models of swiftlets and swallows
##liz edits
library(ENMeval)
library(raster)
library(RColorBrewer)
library(tidyverse)
library(cowplot)
library(maptools)
data("wrld_simpl")

input_dir<-"~/Desktop/swiftlet_data/enviro_data/niche_overlap/all"

## get list of input files
model_files_list<-list.files(path=input_dir , recursive = TRUE  , pattern=".RData" , full.names = TRUE)
model_files_list

## read all models into a list
enmeval_models_list<-list()

for(i in 1:length(model_files_list)){
  # new temporary environment for just the enmeval results
  temp<-new.env()
  # load results into new environment
  load(model_files_list[i], temp)
  # get list of objects in the temporary env (should be only the one enmeval result)
  models<-ls(temp)
  # rename the enmeval model to generic "mod", can use the [[1]] because there is only 1 object in environment
  mod<-get(models[[1]], temp)
  # set up naming for the model
  new_name<-basename(model_files_list[i])
  new_name<-gsub(".RData","",new_name)
  # add model to list with correct name
  enmeval_models_list[[new_name]]<-mod
}

rm(temp)
names(enmeval_models_list)

#calculate max boyce index for each model in the list, output into new list
best_mods_list<-list()

for(i in seq_along(enmeval_models_list)){
  ## pull out the model results
  res<-eval.results(enmeval_models_list[[i]])
  ## make rm numeric so it can be filtered to less than 6 (and not underfit model)
  res$rm<-as.numeric(as.character(res$rm))
  res<-filter(res, rm < 6)
  # get max cbi from results
  opt.boyce<- res%>% filter(cbi.val.avg==max(cbi.val.avg, na.rm=TRUE))
  ## put into list
  best_mods_list[[i]]<-opt.boyce
}

names(best_mods_list)<-names(enmeval_models_list)
names(best_mods_list)

## get the raster layer corresponding to the best model
best_mods_rasters<-list()

for(i in seq_along(best_mods_list)){
  best_mods_rasters[[i]]<-eval.predictions(enmeval_models_list[[i]])[[best_mods_list[[i]]$tune.args]]
}
names(best_mods_rasters)<-names(enmeval_models_list)

## plot all best models, make sure they look reasonable
par(mfrow=c(2,2))
cols<-rev(brewer.pal(10, "Spectral"))
for(i in seq_along(best_mods_rasters)){
  plot(best_mods_rasters[[i]], col = cols, main=names(best_mods_rasters[i]))
}

################## get variable importance for best models (because I was curious)
var_imp_list<-list()

for(i in seq_along(enmeval_models_list)){
  var_imp<-enmeval_models_list[[i]]@variable.importance[[best_mods_list[[i]]$tune.args]]
  var_imp_list[[i]]<-var_imp
}
names(var_imp_list)<-names(enmeval_models_list)
## turn into dataframe for plotting
var_imp_df<-do.call(rbind, var_imp_list)
## add columns for overall model
var_imp_df$model<-gsub("\\.\\d+", "", rownames(var_imp_df))
## split model names so each component is its own column
labels_df<-as.data.frame(str_split_fixed(var_imp_df$model, "_", 3))
names(labels_df)<-c("env_vars","hyp","species")
## add new columns to df
var_imp_df<-as.data.frame(cbind(var_imp_df, labels_df))
head(var_imp_df)

# filter to permutation importance > 3% to make plot easier to interpret
var_imp_df_filter<-filter(var_imp_df, permutation.importance > 3)

## make different color palettes for climate & land use variables
cols_chelsa<-c("#ffeda0", "#ffffcc", brewer.pal(9,"OrRd"))
cols_glc<-rev(c(brewer.pal(6, "BuGn")))

# plot
ggplot(var_imp_df_filter, aes(species, permutation.importance, fill=variable))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = c(cols_chelsa, cols_glc))+
  theme_cowplot()


################## make best models into a raster stack, plot
sabah_stack<-stack(best_mods_rasters)
plot(sabah_stack)

## remove the uniform swiftlet- not present in our study area
best_mods_rasters<-best_mods_rasters[1:6]
sabah_stack<-stack(best_mods_rasters)

## crop to sabah
sabah_stack<-crop(sabah_stack, extent(116.5, 119, 5, 6.5))
plot(sabah_stack)

### create a stack that is the average suitability of each cell across all species
average_suit<-calc(sabah_stack, mean)

breaks <- seq(0, 1, by = 0.1)
plot(average_suit, col = cols, breaks=breaks)

## histogram of suitabilities
breaks <- seq(0, 1, by = 0.2)
hist(values(average_suit), breaks = breaks)

# plot with bins of 0.2 rather than 0.1
breaks <- seq(0, 1, by = 0.2)
plot(average_suit, col = cols, breaks=breaks)

## get values of pixels per bin- can change breaks to change bin size
breaks <- seq(0, 1, by = 0.2)
hist_counts <- hist(values(average_suit), breaks = breaks, plot = FALSE)$counts
## you'll want to randomly sample from these bins
hist_counts

####### total suitability rather than average
## max = 6
total_suit<-calc(sabah_stack, sum)
breaks <- seq(0, 6, by = 0.1)
plot(total_suit, col = cols)

# Specify the file path and name for the output raster file
output_file <- "~/Desktop/swiftlet_data/enviro_data/sabah_stack.asc"
# Use the writeRaster function to save the raster layer
writeRaster(total_suit, filename = output_file, format = "ascii")

##input into arcgis... summarize 1km polygons within the 10km suitability 
##additionally create subsets based off of how close to roads and infrastructure
## input back into here with attribute tables for stratified random sampling

### creating gridded raster 
#cropping
setwd("~/Desktop/swiftlet_data/enviro_data/present_day/")
layers <-list.files("~/Desktop/swiftlet_data/enviro_data/CHELSA/clim_bio", full.names=TRUE)
output_dir<-"~/Desktop/swiftlet_data/enviro_data/gridded_layers/"
library(raster)
library(terra)
temp_raster <- raster(layers[1])
plot(temp_raster, col=c(topo.colors(200)), axes=FALSE, box=FALSE)
cropped_raster<-crop(temp_raster, extent(116, 119, 5, 7))
plot(cropped_raster)
gridded_raster <- rasterToPolygons(cropped_raster)
plot(gridded_raster) #check
shapefile(gridded_raster, filename="~/Desktop/swiftlet_data/enviro_data/gridded_layers/proper_1km_bio1.shp")
raster::writeRaster(gridded_raster, filename=file.path(output_dir,"1km_gridded_bio1"), format="ascii", bylayer=TRUE, overwrite=TRUE)

###loading in layer after summarize within from arcgis
###loading in polygons with summarize_within features on suitability scores from 1 (no) to 4 (high)
#read in table and layer of polygons
library(sf)
library(rnaturalearth)
# Read the shapefile
shapefile <- st_read("~/Desktop/swiftlet_data/enviro_data/1km_sum/1km_sum.shp")
# Print summary of the shapefile
print(shapefile)
# Plot the shapefile
plot(shapefile)
# Get a polygon layer representing land areas
world_land <- ne_countries(returnclass = "sf")
# Perform a spatial overlay to keep only the polygons within land areas
shapefile <- st_make_valid(shapefile)
shapefile <- st_intersection(shapefile, world_land)
# Print summary of the shapefile
print(shapefile)
# Plot the shapefile
plot(shapefile)
# Set the seed for reproducibility
set.seed(123)
# Create an empty list to store the selected polygons
selected_polygons <- list()
# Define the unique suitabilit scores
unique_scores <- unique(shapefile$suitabilit)
# Set the number of polygons to select per score
polygons_per_score <- 20
# Loop through each unique suitabilit score
for (score in unique_scores) {
  # Subset the shapefile to only include polygons with the current score
  subset_shapefile <- shapefile[shapefile$suitabilit == score, ]
  # If there are fewer than 2 polygons with this score, select all of them
  if (nrow(subset_shapefile) <= polygons_per_score) {
    selected_polygons[[as.character(score)]] <- subset_shapefile
  } else {
    # Randomly select 2 polygons without replacement
    selected_indices <- sample(1:nrow(subset_shapefile), polygons_per_score, replace = FALSE)
    selected_polygons[[as.character(score)]] <- subset_shapefile[selected_indices, ]
  }
}
# Combine the selected polygons into a single dataframe
selected_polygons_df <- do.call(rbind, selected_polygons)
# Print the selected polygons
print(selected_polygons_df)
plot(selected_polygons)
# Specify the file path where you want to save the shapefile
output_shapefile <- "selected_polygons.shp"
# Write the selected polygons to a shapefile
st_write(selected_polygons_df, output_shapefile)

##### END OF SCRIPT
##extra liz edits
### overlay map with roads and cities (still working on this)
library("rnaturalearth")
library("rnaturalearthdata")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
