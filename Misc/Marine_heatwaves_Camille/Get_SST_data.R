library(tidyverse)
library(raster)

site_data <- read_csv(file = "Data/rls_site_ls.csv") 

colSums(is.na(site_data)) # No missing data

# Get the extent of the region for cropping 
site_data %>% with(c(min(SiteLat),
                   max(SiteLat),
                   min(SiteLong),
                   max(SiteLong)))

ROI <- extent(112 ,157,-25,-8)

# Get the lat and long
xy <- dplyr::select(site_data, SiteLong, SiteLat)

# Load in the files and crop as we go 
dhw_data <- crop(stack(x = "Data/ct5km_dhw-max_v3.1_2016.nc", varname = "degree_heating_week"), ROI)
ssta_data <- crop(stack(x = "Data/ct5km_ssta-max_v3.1_2016.nc", varname = "sea_surface_temperature_anomaly"), ROI)

# Take a look at the data
plot(dhw_data)
plot(ssta_data)

# Save and load in as a tiff
writeRaster(dhw_data, filename = "Data/cropped_dhw_max_2016.tif", format="GTiff", overwrite = TRUE)
writeRaster(ssta_data, filename = "Data/cropped_ssta_max_2016.tif", format="GTiff", overwrite = TRUE)
writeRaster(dhw_data, filename = "Data/cropped_dhw_max_2016.raster", format="raster", overwrite = TRUE)
writeRaster(ssta_data, filename = "Data/cropped_ssta_max_2016.raster", format="raster", overwrite = TRUE)
# cropped_dhw<-raster("Data/cropped_dhw_max_2016.tif") 
# all(as.vector(extract(dhw_data, xy)) - extract(cropped_dhw, xy) < 0.000001) # Proof that they are the same thing
# plot(cropped_dhw)

# Extract the data we need based on lat and long
site_data$maxDHW <- as.vector(extract(dhw_data, xy))
site_data$maxSSTA <- as.vector(extract(ssta_data, xy))

# The First time round - whole file is 600mb so good to make this a little smaller
# climatology_data <- crop(stack(x = "Data/noaa_crw_thermal_history_climatology_v2.1.nc", varname = "clim_monthly"), ROI)
#writeRaster(climatology_data, filename = "Data/cropped_climatology_data.nc", format="CDF", overwrite = TRUE)

# Load the cropped data back in
climatology_data <- stack("Data/cropped_climatology_data.nc")
# all(extract(look, xy) - extract(climatology_data, xy) < 0.00000001)
plot(climatology_data)
# plot(look)

# Get a matrix with one col for each month
climatology_extract <- extract(climatology_data, xy)

# Append some more appropriate column names
colnames(climatology_extract) <- paste(month.abb, "SST_cl", sep = "_")

site_data <- cbind(site_data, climatology_extract)

# Save
write_csv(site_data, path = "Data/sites_with_sst.csv")
