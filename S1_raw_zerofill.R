##raw data to zero filled
##load packages
library(remotes)
library(sf)
library(rnaturalearth)
library(auk)
library(lubridate)
library(gridExtra)
library(tidyverse)
library(dplyr)

#make sure you put the .txt files in the same folder
### MAKE EDITS HERE AND REST WILL BE AUTOMATED START
#set ebd path
auk::auk_set_ebd_path("~/Desktop/swiftlet_data/ebird/zf_data/zf_files/gloswi/", overwrite = TRUE)
# resolve namespace conflicts
select <- dplyr::select
#setup data directory
localDir = "~/Desktop/swiftlet_data/ebird/zf_data/zf_files/gloswi/"
setwd(localDir)
dir.create("data", showWarnings = FALSE)
library(auk)
#### if not working make sure you open .tar sampling and .zip occurrence 
#### record to access .txt files
ebd <- auk_ebd("ebd_gloswi1_201501_202012_relJan-2024.txt", 
               file_sampling = "ebd_sampling_relJan-2024.txt")
spname <- c("Collocalia esculenta")
# output files
f_ebd <- file.path("~/Desktop/swiftlet_data/ebird/zf_data/zf_files/gloswi/ebd_gloswi1_201501_202012_relJan-2024.txt")
f_sampling <- file.path("~/Desktop/swiftlet_data/ebird/zf_data/zf_files/gloswi/ebd_sampling_relJan-2024.txt")
# Define the bounding box coordinates
bbox_se_asia_oceania <- c(90,-15,160,25) #c(lng_min, lat_min, lng_max, lat_max) 90,-15,160,25
####AUTOMATED END, RUN REST

#choose filters
ebd_filters <- ebd %>% 
  auk_species(spname) %>% 
  # get_ebird_taxonomy(2020) if names need resolution
  auk_distance(distance = c(0, 5))%>%
  #auk_country(country = c("French Polynesia")) %>%
  auk_date(date = c("2015-01-01", "2020-12-31")) %>% 
  auk_bbox(bbox = c(bbox_se_asia_oceania)) %>%
  auk_protocol(protocol = c("Stationary", "Traveling")) %>% 
  auk_duration(duration = c(0, 300)) %>%
  auk_complete()

#check filters
ebd_filters

# only run if the files don't already exist
if (!file.exists(f_ebd)) {
  auk_filter(ebd_filters, file = f_ebd, file_sampling = f_sampling, overwrite=TRUE)
}

#establish output files, very important lol
tax_sed_filtered<-auk_filter(ebd_filters, file="ebd_filtered.txt", file_sampling="sampling_filtered.txt")
tax_sed_filtered

##importing and zero filling--------------------------------------------------------------------------
ebd_zf <- auk_zerofill(tax_sed_filtered, collapse = TRUE)
ebd_zf
table(ebd_zf$species_observed)

# function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# clean up variables
ebd_zf <- ebd_zf %>%
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X",
                                NA_character_, observation_count,),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling",
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

#accounting for variation in detectability----------------------------------------------------------------
# additional filtering
ebd_zf_filtered <- ebd_zf %>%
  filter(
    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 5,
    #years of data from 2018 onwards
    year>= 1990,
    # 10 or fewer observers
    number_observers <= 10)

ebird <- ebd_zf_filtered %>% 
  select(scientific_name, 
         observation_count, species_observed, 
         country, country_code, latitude, longitude,
         observation_date, year, day_of_year)
head(ebird$country)
setwd(file.path(localDir, "data"))
# Concatenate species name with "_zf.csv"
filename <- paste(spname, "_zf.csv", sep = "")
filename <- gsub(" ", "_", filename) #remove spaces
# Write the data frame to the file
write_csv(ebird, filename, na = "")

###### extra code
#test for unique lat/long points
test <- unique(ebird[,c(6:7)])
dim(test)
head(test)
plot (test)
