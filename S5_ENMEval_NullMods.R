## new ENMEval models 
install.packages("ecospat")
install.packages("ape")
install.package = 'rgdal'
library("rgdal")
library(ape)
library(ecospat)
library(ENMeval)
library(raster)
library(RColorBrewer)
library(maptools)
library(tidyverse)
library(rasterVis)
data("wrld_simpl")

## load occurrence
input_dir = "~/Desktop/swiftlet_data/ebird/zf_data/zf_files/pluglo"
enviro_dir = "~/Desktop/swiftlet_data/enviro_data/DS_present_day"

## load occurrence data
setwd(input_dir)
occs_data <- read.csv(file.path(input_dir,"thinned_pres/pluglo_thin1.csv"))
head(occs_data,20)
occs_data$species_observed <- TRUE
head(occs_data,20)

#subset for presence
presence <- subset(occs_data, species_observed==TRUE)
dim(presence)
coords<-presence[,2:3]
head(coords)

# load environmental layers -----------------------------------------------
setwd(enviro_dir)
files<-list.files(file.path(enviro_dir), pattern=".asc")
files #35 with lithology
layers<-stack(files)
crs(layers)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
layers #check resolution and crs

## Background points zerofilled -------------------------------------------------
# bg.buff<- subset(welswallow, species_observed==FALSE)
# dim(bg.buff)
# head(bg.buff)
# bg.buff<-bg.buff[,2:3]
# head(bg.buff)

## Random Sampling background points--------------------------------------------
## make sure seed is set for repeatability of background point generation
set.seed(13)
bg.buff <- dismo::randomPoints(layers[[11]], p = coords, n = 20000) %>% as.data.frame()
colnames(bg.buff) <- colnames(coords)
dim(bg.buff)
head(bg.buff)
#20,000 random sample points

# check that every cell is covered by background points
plot(layers[[11]], main = names(layers[[11]]))
points(bg.buff, pch = 20, cex = 0.2)
points(coords, pch=16, cex=0.2, col="red")


# create subsets of layers for different hypotheses ---------------------------
## these are layers identified by select07 method as uncorrelated
## list of layers is saved in that script so input can be read in here

## SELECT 07
pred_sel<-read.table(paste0("~/Desktop/swiftlet_data/enviro_data/pluglo_uncorrelated_var.txt"), header=FALSE)
pred_sel
layers

pred_sel<-layers[[c("CHELSA_bio4_1981","CHELSA_bio2_1981","CHELSA_bio5_1981", "glc_shv10_03",
                    "glc_shv10_05","glc_shv10_04","CHELSA_bio8_1981","CHELSA_bio9_1981","glc_shv10_01","glc_shv10_02",
                    "CHELSA_bio18_1981","glc_shv10_07","glc_shv10_09","glc_shv10_06","CHELSA_bio16_1981","glc_shv10_08")]]

## subset the rasterstack of layers to combinations for each hypothesis
all_predictors<-pred_sel #all layers
human_predictors<- layers[[c("glc_shv10_03","glc_shv10_05","glc_shv10_04","glc_shv10_01","glc_shv10_02",
                            "glc_shv10_07","glc_shv10_09","glc_shv10_06","glc_shv10_08")]] #hyde only
climate_predictors<-layers[[c("CHELSA_bio4_1981","CHELSA_bio2_1981","CHELSA_bio5_1981",
                              "CHELSA_bio8_1981","CHELSA_bio9_1981","CHELSA_bio18_1981","CHELSA_bio16_1981")]] #climate only

#check
names(all_predictors)
names(climate_predictors)
names(human_predictors)

# ##lithology needs to be outlined as a categorical variable 
# layers$lithology_ #check
# layers$lithology_ <- raster::as.factor(layers$lithology_) #categorize
# #ur a queen continue 

#plot
plot(all_predictors[[9]], main=names(all_predictors[[9]]))
points(bg.buff, pch = 20, cex = 0.2)
points(coords, pch=16, cex=0.2, col="red")

## NA values
plot(all_predictors[[11]], main=names(all_predictors[[11]]), colNA="red")
plot(wrld_simpl, add=TRUE)

### running models -------------------------------------------------------
tune.args=list(fc = c("L", "LQ","H", "LQP", "LQHP", "LQHPT"), rm = c(seq(0.5,4,0.5),5:8))

#SEL7
Sel7_all <- ENMevaluate(coords, all_predictors, bg.buff, partitions="block", tune.args=tune.args, algorithm='maxent.jar', parallel=TRUE, numCores = 7)
Sel7_climate <- ENMevaluate(coords, climate_predictors, bg.buff, partitions="block", tune.args=tune.args, algorithm='maxent.jar', parallel=TRUE, numCores = 7)
Sel7_hyde <- ENMevaluate(coords, human_predictors, bg.buff, partitions="block", tune.args=tune.args, algorithm='maxent.jar', parallel=TRUE, numCores = 7)


## save models with closer cropped extent
output_dir<-"~/Desktop/swiftlet_data/enviro_data/SDMs/"
save(Sel7_all, file=paste0(output_dir, "Sel7_all_pluglo.RData"))
save(Sel7_climate, file=paste0(output_dir, "Sel7_climate_pluglo.RData"))
save(Sel7_hyde, file=paste0(output_dir, "Sel7_hyde_pluglo.RData"))

### if reading in from here make sure load your .RData --------------------
load (output_dir, "Sel7_all_pluglo.RData")
#Sel7
model<-Sel7_all #all layers Sel7
#model<-Sel7_climate #climate only Sel7
#model<-Sel7_hyde #hyde only Sel7

# plotting
cols <- colorRampPalette(rev(brewer.pal(10, "Spectral")))
countries <- maps::map("world", plot=FALSE) 
countries <- map2SpatialLines(countries, proj4string = CRS("+proj=longlat"))

##### model selection 
res <- eval.results(model)
model
res
#############max boyce
opt.boyce<- res%>% filter(cbi.val.avg==max(cbi.val.avg))
ENMeval_all_opt.boyce<-opt.boyce
### variable importance plots
mod.seq <- eval.models(model)[[opt.boyce$tune.args]]
mod.seq@lambdas
plot(mod.seq, type = "cloglog") #save as png 

pred.boyce <- eval.predictions(model)[[opt.boyce$tune.args]]
plot(pred.boyce)

plot(wrld_simpl, add=T)
##dimensions to save png: 880, 535
levelplot(pred.boyce, margin=FALSE, main="Plume-toed Swiftlet Select07 All Predictors", col.regions=cols, at=seq(0, 1, len=20), 
          xlab = "Longitude", ylab = "Latitude") 

model
####NEXT STEPS
# null models -------------------------------------------------------------
# We first run the null simulations with 100 iterations to get a reasonable null distribution 
# for comparisons with the empirical values
##load models
load("~/Desktop/swiftlet_data/enviro_data/SDMs/Sel7_all_pluglo.RData")
#run null for each model hypothesis 
### change fc and LQ and rm and change according to model that loads
### save best perform as mod.setting so more generalized/repeatable in future
mod.null <- ENMnulls(Sel7_climate, mod.settings = list(fc = "LQ", rm = 8), no.iter = 100)

# We can inspect the results of each null simulation.
null1<- null.results(mod.null) %>% head()
# And even inspect the results of each partition of each null simulation.
null2<- null.results.partitions(mod.null) %>% head()
# For a summary, we can look at a comparison between the empirical and simulated results.
null3<- null.emp.results(mod.null)
# Finally, we can make plots of the null model results as a histogram.
evalplot.nulls(mod.null, stats = c("auc.val"), plot.type = "histogram")
evalplot.nulls(mod.null, stats = c("cbi.val"), plot.type = "histogram")

###all null results
evalplot.nulls(mod.null, stats = c("or.10p", "auc.val", "cbi.val"), plot.type = "histogram")
# Or we can visualize the results with a violin plot.
evalplot.nulls(mod.null, stats = c("or.10p", "auc.val", "cbi.val"), plot.type = "violin")

#keep files
write.csv(null1, "~/Desktop/null/ENMhyde_null1.csv", 
          row.names = FALSE)

write.csv(null2, "~/Desktop/null/ENMhyde_null2.csv", 
          
          row.names = FALSE)
write.csv(null3, "~/Desktop/null/ENMhyde_null3.csv", 
          row.names = FALSE)

##all done!! congrats b

### niche overlap between models
## take best performance all predictors model and run niche overlap analysis
# load data do each separate before loading the next
library(ENMeval)
library(dplyr)     
library(purrr)     
library(stringr)
load("~/Desktop/swiftlet_data/enviro_data/SDMs/Sel7_climate_pacswa.RData")
#Sel7
model<-Sel7_climate
##### model selection 
res <- eval.results(model)
model
res
##### max boyce
opt.boyce<- res%>% filter(cbi.val.avg==max(cbi.val.avg,na.rm=TRUE))
ENMeval_all_opt.boyce<-opt.boyce
pacswa.model <- eval.predictions(model)[[opt.boyce$tune.args]]
#check
pluglo.model
ebnger.model
pacswa.model
barswa.model
uniswi.model
blnswi.model
monswi.model

### calc.niche.overlap
predictors<- stack(pacswa.model,barswa.model,ebnger.model,blnswi.model,monswi.model,uniswi.model,pluglo.model)
overlap <- calc.niche.overlap(predictors, overlapStat = "D")
overlap

##niche breadth for each
library(ENMeval)
library(ENMTools)
library(terra)
# Convert RasterLayer to SpatRaster
model_spat <- terra::rast(pacswa.model)
# Now you can use raster.breadth() with the converted SpatRaster object
breadth_result <- raster.breadth(model_spat, verbose = FALSE)
breadth_result

###END OF SCRIPT

  