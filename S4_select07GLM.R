##in this script
#use thinned data
#make final thinned dataset with pres and abs
#plot zero filled absences along with presences
#variable correlations

#load packages
library(raster)
library(RColorBrewer)
library(maptools)
library(tidyverse)
library(rasterVis)
library(ecospat)
data("wrld_simpl")
library(doParallel)
library(foreach)
detectCores()
library(tidyverse)

###set directories
enviroDir = "~/Desktop/swiftlet_data/enviro_data/DS_present_day"
localDir = "~/Desktop/swiftlet_data/ebird/zf_data/zf_files/pluglo"
setwd(enviroDir)
outputname <- c("pluglo")

#load in present day layers
files<-list.files(file.path(enviroDir), pattern=".asc")
files #35 with lithology
layers<-stack(files)

#load occurrence data
setwd(localDir)
occs_data <- read.csv(file.path(localDir,"thinned_pres/pluglo_thin1.csv"))
head(occs_data,20)
occs_data$species_observed <- 1 ##adding back on species observed column
#visualise
plot(layers[[1]], main=names(layers[[1]]))
points(occs_data[,2:3], pch=16, cex=.5)

#coords object with pres data
coords<-occs_data[,2:3]
points(coords, cex=0.2)
head(coords,20)
#test coords plot
plot(layers[[1]], main=names(layers[[1]]))
points(coords)
plot(wrld_simpl, add=TRUE) #check alignment

## check resolution and coordinate ref of layers
crs(layers)<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
layers

#### MAKE SURE YOU KNOW WHICH BACKGROUND POINT TESTING/RUNNING
## absense if using zerofilling background --------------------------------
# #subset absence and plot
# absence<-subset(swallow_zf, species_observed==FALSE)
# points(absence[,2:3], pch=16, cex=.2, col="red")

# add random generated absence points -------------------------------------
set.seed(13)
bg.buff <- dismo::randomPoints(layers[[11]], p=coords, 20000) %>% as.data.frame()
colnames(bg.buff) <- colnames(coords)
head(bg.buff)
##check how we good coverage (every cell).
plot(layers[[11]], main = names(layers[[11]]))
points(bg.buff, pch = 20, cex = 0.2)
points(coords, pch=16, cex=0.2, col="red")

## add to data frame with coords
head(coords)
head(bg.buff)
coords$occurrence<-rep(1, length(coords$longitude))
bg.buff$occurrence<-rep(0, length(bg.buff$longitude))

# dataset with both presence and absences ((should make column from above with 0/1))
coords<-rbind(coords, bg.buff)
table(coords$occurrence)
head(coords)
##combined coords and bg.buff = occs_data with presences and absences
head(occs_data)
table(occs_data$species_observed)
table(coords$occurrence)

###### 
# extract layers to set up glms -------------------------------------------
#extract information from masked layers
m<-raster::extract(layers, coords[,1:2])
m
#extract data from both presences and background points
data<-as.data.frame(cbind(coords[,1:2], m))

#summary should have values for each layer
names(data)
tail(data)
summary(data)
#normal to have NAs this is how many that were with welcome swallows
#more NAS because of cropping so close it doesnt include absence points outside extent
# Plot presences in black and absences in red 
plot(layers[[4]], main=names(layers[[4]]))
points(data[coords$occurrence==1,1:2], pch=16, cex=0.2, col="black")
points(data[coords$occurrence==0,1:2], pch=16, cex=0.2, col="red")

# assess correlations among predictor variables ---------------------------
#Correlation among predictors
library(corrplot)
# We first estimate a correlation matrix from the predictors.
#Generally, correlations below |r|<0.7 are considered unproblematic (or below |r|<0.5 as more conservative threshold).
# We use Spearman rank correlation coefficient, as we do not know # whether all variables are normally distributed.
cor_mat_points <- cor(data[,-c(1:2)], method='spearman', use="pairwise.complete.obs")
cor_mat_points
# We can visualise this correlation matrix. For better visibility,
# we plot the correlation coefficients as percentages.
#positive correllations in blue, negative correlations in red
par(mfrow = c(1, 1))
corrplot.mixed(cor_mat_points, tl.pos='lt', tl.cex=0.6, number.cex=0.5, addCoefasPercent=T)

var.imp <- function (predictor, response) {
  AIC(glm(response ~ predictor + I(predictor^2), binomial)) }

##glm usage: 
# glm(formula, family = gaussian, data, weights, subset,
# na.action, start = NULL, etastart, mustart, offset,
# control = list(...), model = TRUE, method = "glm.fit",
# x = FALSE, y = TRUE, singular.ok = TRUE, contrasts = NULL, ...)


#apply this function (above) to all predictor variables in the data
#this includes correlated variables 
aic_imp <- apply(data[,-c(1:2)], 2, var.imp, response=coords$occurrence)
# Sort the predictors; the lower the AIC the better:
sort(aic_imp)

#variable selection: removing highly correlated variables
#need to identify all pairs of variables that have correlation |r|>0.7 and remove the less important variable
#test for pairwise correlation for the variables (if they are above 0.7 might have a problem — 
#find which is more important for explaining distribution of my species)
select07 <- function(predictor_dat, response_dat, cor_mat=NULL, threshold=0.7){
  # Function for calculating AIC - we use univariate GLMs with linear and quadratic term
  var.imp <- function (predictor, response)
  {
    AIC(glm(response ~ predictor + I(predictor^2), binomial))
  }
  # Calculate AIC for all predictor variables
  aic_imp <- apply(predictor_dat, 2, var.imp, response= response_dat)
  # Names of sorted variables
  sort_imp <- names(sort(aic_imp))  
  # Calculate correlation matrix if not provided in function call
  if (is.null(cor_mat)) {
    cor_mat <- cor(predictor_dat, method='spearman')
  }
  # Identifies correlated variable pairs:
  diag(cor_mat)=NA
  pairs <- which(abs(cor_mat)>= threshold, arr.ind=T)  
  # Identify which variables should be excluded
  exclude <- NULL
  for (i in 1:length(sort_imp))
  {
    if ((sort_imp[i] %in% row.names(pairs))& 
        ((sort_imp[i] %in% exclude)==F)) {
      cv <- cor_mat[setdiff(row.names(cor_mat),exclude),sort_imp[i]] 
      cv <- cv[setdiff(names(cv),sort_imp[1:i])]
      exclude <- c(exclude,names(which((abs(cv)>=threshold))))
    } 
  }
  # Select set of weakly correlated predictors:
  pred_sel <- sort_imp[!(sort_imp %in% exclude)]
  # Return list with AIC, correlation matrix, and final predictors:
  return(list(AIC=sort(aic_imp), cor_mat=cor_mat, pred_sel=pred_sel)) 
}

# apply function
var_sel <- select07(predictor_dat=data[,-c(1:2)], response_dat=coords$occurrence,
                    cor_mat=cor_mat_points, threshold=0.7) 
pred_sel <- var_sel$pred_sel
pred_sel

## write to file to import into next set of scripts
write.table(pred_sel, "~/Desktop/swiftlet_data/enviro_data/pluglo_uncorrelated_var.txt", 
            row.names = FALSE, col.names = FALSE)

#yay first set of uncorrelated environmental predictors! 
