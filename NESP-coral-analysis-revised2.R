
##############################################################
# NESP CORAL ANALYSIS - REVISED - CM 29/11/22                #
##############################################################


# Load libraries ------------
rm(list = ls())

#library(RLSMetrics)
library(RColorBrewer)
library(PBSmapping)
library(goeveg)
library(ggplot2)
library(ggplotify)
library(gridExtra)
library(tidyverse)
library(vegan)
library(betapart)
library(labdsv)
#library(mvpart) #devtools::install_github("cran/mvpart")
library(adespatial)
library(psych)
library(gbm)
library(pdist)
library(ggpubr)
library(sp)
library(sf)
library(tmap)
library(raster)
library(dplyr)
library(rgdal)
library(ggmap)
library(maptools)
library(cowplot)
library(pairwiseAdonis)

# Functions ----------------

addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

summary_indval <- function (object, p = 0.05) 
{   tmp <- data.frame(object$maxcls[object$pval <= p], round(object$indcls[object$pval <= p], 3), object$pval[object$pval <= p])
names(tmp) <- c("cluster", "indicator_value", "probability")
data.frame(tmp[order(tmp$cluster, -tmp$indicator_value),])
}

# Load and filter data --------------

#setwd("/Users/Camille/Dropbox/My documents/Projects/Future Fellowship/RESEARCH/NESP-coral-analysis-revised")
load("NESP-coral-analysis-revised2.RData")
source("brt.functions.R")

# Load covariates
# dhw.raster <- raster::raster("Marine_heatwaves_Camille/Data/cropped_dhw_max_2016.grd")
# ssta.raster <- raster::raster("Marine_heatwaves_Camille/Data/cropped_ssta_max_2016.grd")

SiteCov_msec <- read.csv("MSEC/msec_out_npp_sst_wave.csv")
SiteCov_msec_rd <- read.csv("MSEC/msec_out_ReefLandArea_HumanPop.csv")
SiteCov <- SiteCov_msec %>% inner_join(SiteCov_msec_rd) %>% mutate(SiteCode = NULL, Region = NULL)

aus <- importShapefile("shapefiles/250K_coastline", readDBF=FALSE)
aus2 <- aus %>% dplyr::select(group=PID, POS=POS,long=X,lat=Y)
aus.sf <- st_read("shapefiles/250K_coastline.shp")

reef <- importShapefile("shapefiles/SDE_OWNER_crcgis_land250", readDBF=FALSE)

# Fix shapefile (takes forever!)
shape.clip <- closePolys(clipPolys(reef, xlim=c(142,155), ylim=c(-25,-10)))
reefs <- as.PolySet(subset(shape.clip, PID != 2401), projection = "LL")
reef2 <- reefs %>% dplyr::select(group=PID, POS=POS,long=X,lat=Y)

# Load temperature anomalies
SiteSSTA <- read.csv("Marine_heatwaves_Camille/sites_with_sst.csv")
SiteSSTA$MMM <- apply(subset(SiteSSTA, select = Jan_SST_cl:Dec_SST_cl), 1, max)

SiteSSTA <- SiteSSTA %>% dplyr::select(SiteLat, SiteLong, MMM,
                          maxDHW, maxSSTA, Mean_SST_cl_noaa = Mean_SST_cl, SD_SST_cl_noaa = SD_SST_cl)

# Import new PQ extract and filter with survey list
bent.dat <- read.csv("PQ data/NESP_PQ_endpoint.csv")

bent.dat$survey_date <- as.Date(bent.dat$survey_date, tryFormats = "%d/%m/%Y")
bent.dat$survey_year <- format(bent.dat$survey_date, '%Y')
bent.dat$survey_year <- as.numeric(bent.dat$survey_year)
  
bent.dat$survey_id <- factor(bent.dat$survey_id)
bent.dat$dataset_id <- factor(bent.dat$dataset_id)


# Remove duplicate surveys
dataset_id_duplicates <- c("912346798_RLS Australian Coral Species List_982", 
                           "912346792_RLS Australian Coral Species List_917",
                           "912351889_RLS Catalogue_946", 
                           "912352031_RLS Catalogue_946",
                           "912352656_RLS Catalogue_939",
                           "912352687_RLS Catalogue_928")

bent.dat <- bent.dat %>% filter(!dataset_id %in% dataset_id_duplicates)

# Filter surveys scored by Emre
bent.dat <- bent.dat %>% filter(label_scheme == "RLS Australian Coral Species List")

# Check unique survey_id/dataset_id match (i.e. each transect is scored only once)
bent.meta <- bent.dat %>% group_by(survey_id, dataset_id) %>% summarize() %>%
  group_by(survey_id) %>% summarize(count = n_distinct(dataset_id))

# Load in Emre's groupings
ET_grouping <- read.csv("PQ data/coral_species_grouping ET_Genus GF.csv")
ET_grouping$Genus_GF <- as.character(ET_grouping$Genus_GF)


# Build coral-only dataset (wide format) ------------


coral.dat <- bent.dat %>% filter(RLE_category == "Coral") %>%
              dplyr::select(survey_id, area, location, site_code, site_name, latitude, longitude, depth,
                            survey_date, survey_year, pre_post = is_NESP_pre.post,
                            RLS_category, label, flag_bleached, percent_cover) 

coral.dat$RLS_category_2 <- sub("_Bleached", coral.dat$RLS_category, replacement = "")

# Calculate total live coral cover (total_live) and total bleached coral cover (total_bleached)
coral.total <- coral.dat %>% group_by(survey_id) %>%
                              summarize(total_live = sum(percent_cover),
                                        total_bleached = sum(percent_cover[flag_bleached == TRUE]))


# Match to coral.total and Genus GF and switch to wide format

coral.dat <- coral.dat %>% left_join(coral.total) %>% left_join(ET_grouping)

coral.w.sp <- coral.dat %>% group_by(survey_id, area, location, site_code, site_name, latitude, longitude, depth,
                                   survey_date, survey_year, pre_post, total_live, total_bleached, label) %>%
                          summarize(percent_cover = sum(percent_cover)) %>%
                          pivot_wider(names_from = label, names_sort = T, values_from = percent_cover, values_fill = 0) %>%
                          data.frame()

coral.w.sp[,12:417] <- round(coral.w.sp[,12:417],2)            
table(round(rowSums(coral.w.sp[,14:417])) == round(coral.w.sp$total_live))


coral.w.gff <- coral.dat %>% group_by(survey_id, area, location, site_code, site_name, latitude, longitude, depth,
           survey_date, survey_year, pre_post, total_live, total_bleached, Genus_GF) %>%
          summarize(percent_cover = sum(percent_cover)) %>%
          pivot_wider(names_from = Genus_GF, names_sort = T, values_from = percent_cover, values_fill = 0) %>%
          data.frame()

coral.w.gff[,12:132] <- round(coral.w.gff[,12:132],2)            
table(round(rowSums(coral.w.gff[,14:132])) == round(coral.w.gff$total_live))

# Map coral labels to RLS categories
coral.map <- coral.dat %>% group_by(survey_id, RLS_category_2, label) %>%
                            summarize(percent_cover = sum(percent_cover)) %>%
                            group_by(RLS_category_2, label) %>% 
                            summarize(mean_cover_present = mean(percent_cover[percent_cover > 0]))

# Export for Emre
# write.csv(coral.map, "EmreExports_2/coral_species_mapping_updated.csv", row.names = F)
# write.csv(coral.dat, "EmreExports_2/coral_data_updated.csv", row.names = F)


# Match coral dataset to covariates --------------

# MSEC data
coral.cov <- coral.w.gff %>% dplyr::select(survey_id, survey_year, pre_post, site_code, site_name, latitude, longitude) %>%
                left_join(SiteCov, by = c("latitude" = "SiteLat", "longitude" = "SiteLong")) %>%
                group_by(survey_id, survey_year, pre_post, site_code, site_name, latitude, longitude) %>%
                summarize_all(mean)

# MSEC data: Fill in NAs
coral.cov$SiteLat_rd[is.na(coral.cov$SiteLat_rd)] <- round(coral.cov$latitude[is.na(coral.cov$SiteLat_rd)],1)
coral.cov$SiteLong_rd[is.na(coral.cov$SiteLong_rd)] <- round(coral.cov$longitude[is.na(coral.cov$SiteLong_rd)],1)

coral.cov.na <- coral.cov[is.na(coral.cov$Mean_SST_cl),]
coral.cov.na <- coral.cov.na %>% dplyr::select(survey_id, survey_year, pre_post, site_code, site_name, latitude, longitude, SiteLong_rd, SiteLat_rd) %>%
                      left_join(SiteCov, by = c("SiteLat_rd", "SiteLong_rd")) %>%
                      mutate(SiteLong = NULL, SiteLat = NULL) %>%
                      group_by(survey_id, survey_year, pre_post, site_code, site_name, latitude, longitude, SiteLong_rd, SiteLat_rd) %>%
                      summarize_all(mean)

coral.cov <- rbind(na.omit(coral.cov), coral.cov.na)

# Survey depth
SurveyDepth <- coral.dat %>% group_by(survey_id) %>% summarise(depth = mean(depth))
coral.cov <- coral.cov %>% inner_join(SurveyDepth)


# Thermal stress data: Match SSTA and DHW to each site in "Post" years
coral.cov <- coral.cov %>% left_join(SiteSSTA, by = c("latitude" = "SiteLat", "longitude" = "SiteLong")) %>%
              group_by(survey_id, survey_year, pre_post, site_code, site_name, latitude, longitude) %>%
              summarize_all(mean)

SiteSSTA$SiteLat_rd <- round(SiteSSTA$SiteLat,1)
SiteSSTA$SiteLong_rd <- round(SiteSSTA$SiteLong,1)
Site_SSTA_NA <- SiteSSTA %>% dplyr::select(SiteLat_rd, SiteLong_rd, MMM, maxDHW, maxSSTA, Mean_SST_cl_noaa, SD_SST_cl_noaa)

coral.cov.na <- coral.cov[is.na(coral.cov$MMM),]
coral.cov.na <- coral.cov.na %>% dplyr::select(survey_id:depth) %>%
  left_join(Site_SSTA_NA, by = c("SiteLat_rd", "SiteLong_rd")) %>%
  group_by(survey_id, survey_year, pre_post, site_code, site_name, latitude, longitude, SiteLong_rd, SiteLat_rd) %>%
  summarize_all(mean)

coral.cov <- rbind(na.omit(coral.cov), coral.cov.na)
coral.cov$maxDHW[coral.cov$pre_post == "Pre"] <- 0
coral.cov$maxSSTA[coral.cov$pre_post == "Pre"] <- 0

coral.cov <- coral.cov %>% mutate(Mean_SST_cl_noaa = NULL, SD_SST_cl_noaa = NULL)

coral.cov <- coral.cov[order(coral.cov$survey_id),]

# Calculate sampling interval (in years) and add as an offset
sampl_interval <- coral.dat %>% group_by(site_code) %>% summarize(Year_Pre = mean(survey_year[pre_post == "Pre"]), Year_Post = mean(survey_year[pre_post == "Post"]))
sampl_interval$interval <- with(sampl_interval, Year_Post - Year_Pre)

coral.cov <- coral.cov %>% left_join(sampl_interval)
coral.cov[,c("Year_Pre", "Year_Post", "interval")] <- round(coral.cov[,c("Year_Pre", "Year_Post", "interval")],1)





# Coral MRT for pre-bleaching surveys at the site level ----------------


coral.w.mrt <- coral.w.gff %>% filter(pre_post == "Pre") %>% group_by(site_code) %>%
                            summarize_at(vars("Acanthastrea.Submassive":"Turbinaria.Foliose.plates"), mean) %>% 
                            mutate(site_code = NULL) %>% data.frame()
coral.w.mrt <- coral.w.mrt[,colSums(coral.w.mrt)>0]
coral.w.mrt.hel <- decostand(coral.w.mrt, "hellinger")

coral.cov.mrt <- coral.cov %>% filter(pre_post == "Pre") %>% group_by(site_code, site_name) %>% 
                              summarize_at(vars("latitude":"interval"), mean) %>% data.frame()

coral.cov.mrt.sub <- coral.cov.mrt %>% dplyr::select(Mean_SST_cl, SD_SST_cl, npp_mean, npp_sd, wave_mean, land_area_20km,
                                                    reef_area_20km, pop2015_20km, depth, MMM) %>% data.frame()
                        

coral.mrt <- mvpart(as.matrix(coral.w.mrt.hel)~., coral.cov.mrt.sub, xv=c("pick"), legend=F)
coral.mrt <- mvpart(as.matrix(coral.w.mrt.hel)~., coral.cov.mrt.sub, size=7, legend=F)

par(mai = c(1,1,.1,1))
plot(coral.mrt)
text(coral.mrt, cex=0.8)
print(coral.mrt$cptable)


# Rename clusters
coral.clusters <- coral.mrt$where
for (i in 1:length(table(coral.clusters))) {
  coral.clusters[coral.mrt$where==as.numeric(names(table(coral.mrt$where))[i])] <- i
}
table(coral.clusters)


# Indicators GF
ind.Coral <- indval(coral.w.mrt, coral.clusters, numitr=100)
ind.Coral.tb <- summary_indval(ind.Coral)

all.indval <- data.frame(ind.Coral[["indval"]], ind.Coral[["pval"]])
# write.csv(ind.Coral.tb, "MRT_Indicator_Corals.csv", row.names = T)
# write.csv(all.indval, "MRT_Indicator_Corals_all.csv", row.names = T)

# Indicator species
coral.w.sp.mrt <- coral.w.sp %>% filter(pre_post == "Pre") %>% group_by(site_code) %>%
  summarize_at(vars("Acanthastrea":"Turbinaria.stellulata"), mean) %>% 
  mutate(site_code = NULL) %>% data.frame()
coral.w.sp.mrt <- coral.w.sp.mrt[,colSums(coral.w.sp.mrt)>0]

ind.Coral.sp <- indval(coral.w.sp.mrt, coral.clusters, numitr=100)
ind.Coral.sp.tb <- summary_indval(ind.Coral.sp)

all.indval.sp <- data.frame(ind.Coral.sp[["indval"]], ind.Coral.sp[["pval"]])
# write.csv(ind.Coral.sp.tb, "MRT_Indicator_Corals-revised-Species.csv", row.names = T)
# write.csv(all.indval.sp, "MRT_Indicator_Corals_all-revised-Species.csv", row.names = T)


# Map communities
mrt_coords <- coral.cov.mrt %>% dplyr::select(site_code, longitude, latitude) %>% mutate(cluster = coral.clusters)

mrt_nws <- ggplot() +
  geom_polygon(data=reef2, aes(long, lat, group=group),  fill="darkgray", color="darkgray") +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  #coord_map(xlim=c(118,128), ylim=c(-20,-11)) +
  xlab(expression(paste(Longitude^o, ~'E'))) +
  ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_point(data = mrt_coords, aes(x=longitude, y=latitude, fill = factor(cluster)), size=2, shape=21) +
  scale_fill_brewer(palette = "Paired", name = "Cluster")+
  theme_light() +
  theme(legend.position = "none") +
  #theme(legend.position = c(.85,.7), legend.direction = "vertical")+
  #ggtitle("CategoryName MRT") +
  #coord_cartesian(xlim = c(142, 157), ylim = c(-25,-9)) # GBR
  coord_cartesian(xlim = c(112, 128), ylim = c(-24,-9)) # NWS

mrt_gbr <- ggplot() +
  geom_polygon(data=reef2, aes(long, lat, group=group),  fill="darkgray", color="darkgray") +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  #coord_map(xlim=c(118,128), ylim=c(-20,-11)) +
  xlab(expression(paste(Longitude^o, ~'E'))) +
  ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_point(data = mrt_coords, aes(x=longitude, y=latitude, fill = factor(cluster)), size=2, shape=21) +
  scale_fill_brewer(palette = "Paired", name = "Cluster")+
  theme_light() +
  #theme(legend.position = "none") +
  theme(legend.position = c(.85,.7), legend.direction = "vertical")+
  #ggtitle("CategoryName MRT") +
  coord_cartesian(xlim = c(142, 157), ylim = c(-24,-9)) # GBR
  #coord_cartesian(xlim = c(112, 128), ylim = c(-24,-12)) # NWS

plot_grid(mrt_nws, mrt_gbr, nrow = 1)

ggplot() +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="black", color="black") +
  geom_rect(data = data.frame(xmin = 142, xmax = 157, ymin = -25, ymax = -9), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), col = "red", fill = "transparent", lwd = 3)+
  geom_rect(data = data.frame(xmin = 112, xmax = 128, ymin = -24, ymax = -12), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), col = "red", fill = "transparent", lwd = 3)+
  theme_void()

# Export as Google Earth file for Emre
coordinates(mrt_coords)<-c("longitude","latitude")         
proj4string(mrt_coords)<-CRS("+proj=longlat +datum=WGS84")
#writeOGR(mrt_coords, dsn="mrt_coords.kml", layer= "cluster", driver="KML")

col_kml <- col2kml(RColorBrewer::brewer.pal("Set3", n = 12)[7:1])
col_kml = c("50F0FA14", "50B44614", "5078FF78", "50008214", "508278F0", "501400E6", "5014B4FA")
data(worldgrids_pal)
# See palettes at http://plotkml.r-forge.r-project.org/settings.php

plotKML::kml(mrt_coords,
             file.name    = "EmreExports_2/mrt_coords_GenusGF.kml",
             points_names = mrt_coords$site_code,
             colour    = mrt_coords$cluster,
             colour_scale = rev(worldgrids_pal[[4]]),
             alpha     = .8,
             size      = .5,
             shape     = "http://maps.google.com/mapfiles/kml/pal2/icon18.png")

data_mrt <- data.frame(mrt_coords, coral.w.mrt)
write.csv(data_mrt, "EmreExports_2/coral_data_mrt_GenusGF.csv", row.names = F)
write.csv(coral.cov.mrt, "EmreExports_2/covariate_data_mrt_GenusGF.csv", row.names = F)

# Link to original coral.dat dataset and export for Emre
coral.dat <- coral.dat %>% left_join(mrt_coords)
write.csv(coral.dat, "EmreExports_2/coral_data_updated_with cluster.csv", row.names = F)

# Check distribution of env covariates in each cluster
par(mfrow = c(4,1), mai = c(.1,.5,.2,.1))
boxplot(coral.cov.mrt.sub$npp_mean ~ coral.clusters, xlab = "Cluster", ylab = "npp_mean")
boxplot(coral.cov.mrt.sub$npp_sd ~ coral.clusters, xlab = "Cluster", ylab = "npp_sd")
boxplot(coral.cov.mrt.sub$MMM ~ coral.clusters, xlab = "Cluster", ylab = "MMM")
boxplot(coral.cov.mrt.sub$depth ~ coral.clusters, xlab = "Cluster", ylab = "Depth")

cov.by.cluster <- data.frame(mrt_coords, coral.cov.mrt.sub) %>% 
  dplyr::select(cluster, npp_mean, npp_sd, depth, MMM) %>%
  group_by(cluster) %>%
  summarize_all(list(min = min, max = max, mean = mean, median = median))
  
cov.by.cluster <- data.frame(cluster = cov.by.cluster[,1], cov.by.cluster[,-1][,order(names(cov.by.cluster)[-1], decreasing = T)])

write.csv(cov.by.cluster, "EmreExports_2/Env_covariates_by_cluster.csv", row.names = F)

# Distance-based RDA: Genus GF level ----------------

# Match cluster with covariates and compute site-level averages
coral.cov <- mrt_coords %>% left_join(coral.cov)
coral.cov.site <- coral.cov %>% group_by(site_code, pre_post, longitude, latitude, cluster, site_name) %>%
                                summarize_at(vars("Mean_SST_cl":"interval"), mean)       

coral.cov.site$clust.prepost <- with(coral.cov.site, paste(pre_post, cluster, sep = "_"))

# Calculate T1 as time between pre-surveys and heatwave, T2 as time between heatwave and post-surveys
coral.cov.site$T1 <- 2016 - coral.cov.site$Year_Pre
coral.cov.site$T2 <- coral.cov.site$Year_Post - 2016

par(mfcol = c(1,2))
boxplot(maxSSTA ~ cluster, data = coral.cov.site[coral.cov.site$pre_post == "Post",])
boxplot(maxDHW ~ cluster, data = coral.cov.site[coral.cov.site$pre_post == "Post",])

# Compute site-level benthic covers
coral.w.gff.site <- coral.w.gff %>% filter(site_code %in% mrt_coords$site_code) %>%
                            group_by(site_code, pre_post) %>%
                            summarize_at(vars("Acanthastrea.Submassive":"Turbinaria.Foliose.plates"), mean) %>%
                            ungroup() %>%
                            mutate(site_code = NULL, pre_post = NULL)
  
coral.w.gff.site <- subset(coral.w.gff.site, select=colSums(coral.w.gff.site)>0)

# Remove sites surveyed only once
coral.w.gff.site <- coral.w.gff.site[!is.na(coral.cov.site$interval),]
coral.cov.site <- coral.cov.site[!is.na(coral.cov.site$interval),]

# Export for Emre
write.csv(coral.w.gff.site, "EmreExports_2/coral.data.site.prepost_genus gf.csv", row.names = F)
write.csv(coral.cov.site, "EmreExports_2/covariate.data.site.prepost.csv", row.names = F)

# ADONIS
#adonis <- adonis(sqrt(coral.w.gff.site) ~ cluster * pre_post, data = coral.cov.site, distance = "bray", add = T)
adonis <- adonis2(sqrt(coral.w.gff.site) ~ clust.prepost + offset(T1) + offset (T2), data = coral.cov.site, distance = "bray", add = T)
adonis$aov.tab

adonis_1 <- adonis(sqrt(coral.w.gff.site[coral.cov.site$cluster ==1,]) ~ pre_post, data = coral.cov.site[coral.cov.site$cluster ==1,], distance = "bray", add = T)
adonis_1$aov.tab

adonis_2 <- adonis(sqrt(coral.w.gff.site[coral.cov.site$cluster ==2,]) ~ pre_post, data = coral.cov.site[coral.cov.site$cluster ==2,], distance = "bray", add = T)
adonis_2$aov.tab

adonis_3 <- adonis(sqrt(coral.w.gff.site[coral.cov.site$cluster ==3,]) ~ pre_post, data = coral.cov.site[coral.cov.site$cluster ==3,], distance = "bray", add = T)
adonis_3$aov.tab

adonis_4 <- adonis(sqrt(coral.w.gff.site[coral.cov.site$cluster ==4,]) ~ pre_post, data = coral.cov.site[coral.cov.site$cluster ==4,], distance = "bray", add = T)
adonis_4$aov.tab

adonis_5 <- adonis(sqrt(coral.w.gff.site[coral.cov.site$cluster ==5,]) ~ pre_post, data = coral.cov.site[coral.cov.site$cluster ==5,], distance = "bray", add = T)
adonis_5$aov.tab

adonis_6 <- adonis(sqrt(coral.w.gff.site[coral.cov.site$cluster ==6,]) ~ pre_post, data = coral.cov.site[coral.cov.site$cluster ==6,], distance = "bray", add = T)
adonis_6$aov.tab

adonis_7 <- adonis(sqrt(coral.w.gff.site[coral.cov.site$cluster ==7,]) ~ pre_post, data = coral.cov.site[coral.cov.site$cluster ==7,], distance = "bray", add = T)
adonis_7$aov.tab

pairwise.adonis2(sqrt(coral.w.gff.site) ~ clust.prepost, data = coral.cov.site)


# CAPSCALE

# capscale.preds <- data.frame(clust.bentCat = bentCat_w_site$clust.bentCat, PrePost = bentCat_w_site$PrePost, clust.prepost,
#                              maxSSTA = ssta_site$maxSSTA, maxDHW = ssta_site$maxDHW,
#                              resurveyed = sampl_years$resurveyed, from2016 = sampl_years$from2016)


capscale <- capscale(sqrt(coral.w.gff.site) ~ clust.prepost + Condition(T2) + Condition(interval), data = coral.cov.site, distance = "bray", add = T)

#capscale <- capscale(sqrt(coral.w.gff.site) ~ cluster * pre_post + Condition(interval), data = coral.cov.site, distance = "bray", add = T)
#capscale <- capscale(sqrt(coral.w.gff.site) ~ 1, data = coral.cov.site, distance = "bray", add = T)

capscale
plot(capscale)
anova(capscale, permutations = how(nperm=99), by = "term")


fit.preds <- data.frame(subset(coral.cov.site, select = c(maxDHW, maxSSTA)), LiveCoralCover = rowSums(coral.w.gff.site))
fit <- envfit(capscale, fit.preds, perm = 0, display = "lc", scaling = "sites", na.rm = T)

centroid.scores <- data.frame(scores(capscale, choices = c(1,2))$centroids)
n.clust <- as.numeric(dim(centroid.scores)[1]/2)
col.pal <- rep(brewer.pal(n.clust, "Paired"),2)

dev.off()
pl <- ordiplot(capscale, type = "none", xlim = c(-3,2), ylim = c(-3,2))
points(pl, "sites", pch = 19, col = "lightgrey", cex = .5)
#points(pl, "sites", pch = 21, col = col.pal, cex = .5, alpha = .5)
points(centroid.scores, pch = 19, col = col.pal)
ellipses <- ordiellipse(capscale, coral.cov.site$clust.prepost, kind = "se", draw = "polygon", alpha = .5, col = col.pal, lty = 0)
arrows(x0 = centroid.scores[(n.clust+1):(n.clust*2),1], y0 = centroid.scores[(n.clust+1):(n.clust*2),2], 
       x1 = centroid.scores[1:n.clust,1], y1 = centroid.scores[1:n.clust,2], col = col.pal, length = .1, angle = 20, lwd = 4)
arrows(x0 = centroid.scores[(n.clust+1):(n.clust*2),1], y0 = centroid.scores[(n.clust+1):(n.clust*2),2], 
       x1 = centroid.scores[1:n.clust,1], y1 = centroid.scores[1:n.clust,2], col = "dimgrey", length = .1, angle = 20, lwd = 1)
plot(fit, col = "black")
legend(3.5, 3.5, legend = 1:n.clust, fill = col.pal, title = "Cluster")


#Export site and species scores
sites.scores <- data.frame(site_code = coral.cov.site$site_code, 
                           clust.prepost = coral.cov.site$clust.prepost,
                           cluster = coral.cov.site$cluster,
                           pre_post = coral.cov.site$pre_post,
                           scores(capscale, choices = c(1:2), display = "sites"),
                           scores(pco, choices = c(1:2), display = "sites"))
names(sites.scores)[7:8] <- c("PCO1","PCO2")

# Site scores between CAP1 =0 and 2
sites.scores.sub <- sites.scores %>% filter(CAP1 > 0 & CAP1 < 2 & pre_post == "Pre")
table(sites.scores.sub$cluster)/sum(table(sites.scores.sub$cluster))

# Site scores at CAP2 >1
sites.scores.sub2 <- sites.scores %>% filter(CAP2 > 1 & pre_post == "Pre")
table(sites.scores.sub2$cluster)
table(sites.scores.sub2$cluster)/sum(table(sites.scores.sub2$cluster))


# Species contributions
species.scores <- data.frame(label = names(coral.w.gff.site), scores(capscale, choices = c(1,2), display = "species"))

sp.list <- data.frame(label = factor(ordiselect(coral.w.gff.site, capscale, fitlim = .1))) %>% inner_join(species.scores)

plot(capscale, type="none", choice = c(1,2), display="species", xlim = c(-1.5,3.2))
points(species.scores[,2:3], pch = 19, col = "grey", cex = .5)
points(sp.list[,2:3], pch = "+", col = "red", cex = .75)
#set.seed(314)
ordipointlabel(capscale, display="species", select=sp.list$label, pch = 19, col = "red", cex = .8, add=T)



# Plot PCO
centroid.scores_PCO <- sites.scores %>% group_by(clust.prepost) %>% summarise(PCO1 = mean(PCO1), PCO2 = mean(PCO2)) %>%
  dplyr::select(PCO1,PCO2) %>% data.frame()

pco.fit <- envfit(pco, fit.preds, perm = 99, scaling = "sites", na.rm = T, choices = c(1,2))

pl_pco <- ordiplot(pco, type = "none", xlim = c(-2,2.5), ylim = c(-2.5,3.5))
points(pl_pco, "sites", pch = 19, col = "lightgrey", cex = .5)
ellipses <- ordiellipse(pco, sites.scores$clust.prepost, draw = "polygon", alpha = .5, col = col.pal, lty = 0)
points(centroid.scores_PCO, pch = 19, col = col.pal)
arrows(x0 = centroid.scores_PCO[(n.clust+1):(n.clust*2),1], y0 = centroid.scores_PCO[(n.clust+1):(n.clust*2),2], 
       x1 = centroid.scores_PCO[1:n.clust,1], y1 = centroid.scores_PCO[1:n.clust,2], col = col.pal, length = .1, angle = 20, lwd = 4)
arrows(x0 = centroid.scores_PCO[(n.clust+1):(n.clust*2),1], y0 = centroid.scores_PCO[(n.clust+1):(n.clust*2),2], 
       x1 = centroid.scores_PCO[1:n.clust,1], y1 = centroid.scores_PCO[1:n.clust,2], col = "dimgrey", length = .1, angle = 20, lwd = 1)
plot(pco.fit, col = "black")
legend(3, 2, legend = 1:n.clust, fill = col.pal, title = "Cluster")

# PCO species plot
pco.species.scores <- data.frame(label = names(coral.w.gff.site), scores(pco, choices = c(1,2), display = "species"))
pco.sp.list <- data.frame(label = factor(ordiselect(coral.w.gff.site, pco, fitlim = .1))) %>% inner_join(pco.species.scores)

plot(pco, type="none", choice = c(1,2), display="species") # , xlim = c(-1.5,3.2)
points(pco.species.scores[,2:3], pch = 19, col = "grey", cex = .5)
points(pco.sp.list[,2:3], pch = "+", col = "red", cex = .75)
#set.seed(314)
ordipointlabel(pco, display="species", select=pco.sp.list$label, pch = 19, col = "red", cex = .8, add=T)


# Distance-based RDA: Species level ----------------

# Match cluster with covariates and compute site-level averages
# coral.cov <- mrt_coords %>% left_join(coral.cov)
# coral.cov.site <- coral.cov %>% group_by(site_code, pre_post, longitude, latitude, cluster, site_name) %>%
#   summarize_at(vars("Mean_SST_cl":"interval"), mean)       
# 
# coral.cov.site$clust.prepost <- with(coral.cov.site, paste(pre_post, cluster, sep = "_"))

#par(mfcol = c(1,2))
#boxplot(maxSSTA ~ cluster, data = coral.cov.site[coral.cov.site$pre_post == "Post",])
#boxplot(maxDHW ~ cluster, data = coral.cov.site[coral.cov.site$pre_post == "Post",])

# Compute site-level benthic covers
coral.w.sp.site <- coral.w.sp %>% filter(site_code %in% coral.cov.site$site_code) %>%
  group_by(site_code, pre_post) %>%
  summarize_at(vars("Acanthastrea":"Turbinaria.stellulata"), mean) %>%
  ungroup() %>%
  mutate(site_code = NULL, pre_post = NULL)

coral.w.sp.site <- subset(coral.w.sp.site, select=colSums(coral.w.sp.site)>0)

# Remove sites surveyed only once
#coral.w.sp.site <- coral.w.sp.site[!is.na(coral.cov.site$interval),]
#coral.cov.site <- coral.cov.site[!is.na(coral.cov.site$interval),]

# Export for Emre
write.csv(coral.w.sp.site, "EmreExports_2/coral.data.site.prepost_species.csv", row.names = F)
#write.csv(coral.cov.site, "EmreExports_2/covariate.data.site.prepost.csv", row.names = F)

# CAPSCALE

# capscale.preds <- data.frame(clust.bentCat = bentCat_w_site$clust.bentCat, PrePost = bentCat_w_site$PrePost, clust.prepost,
#                              maxSSTA = ssta_site$maxSSTA, maxDHW = ssta_site$maxDHW,
#                              resurveyed = sampl_years$resurveyed, from2016 = sampl_years$from2016)


capscale_sp <- capscale(sqrt(coral.w.sp.site) ~ clust.prepost + Condition(interval), data = coral.cov.site, distance = "bray", add = T)
capscale2_sp <- capscale(sqrt(coral.w.sp.site) ~ cluster * pre_post + Condition(interval), data = coral.cov.site, distance = "bray", add = T)
pco_sp <- capscale(sqrt(coral.w.sp.site) ~ 1, data = coral.cov.site, distance = "bray", add = T)

capscale_sp
plot(capscale_sp)
anova(capscale_sp, permutations = how(nperm=99), by = "term")

fit.preds <- data.frame(subset(coral.cov.site, select = c(maxDHW, maxSSTA)), LiveCoralCover = rowSums(coral.w.sp.site))
fit <- envfit(capscale_sp, fit.preds, perm = 0, display = "lc", scaling = "sites", na.rm = T)

centroid.scores_sp <- data.frame(scores(capscale_sp, choices = c(1,2))$centroids)
n.clust <- as.numeric(dim(centroid.scores)[1]/2)
col.pal <- rep(brewer.pal(n.clust, "Paired"),2)

pl <- ordiplot(capscale_sp, type = "none", xlim = c(-3,4.5))
points(pl, "sites", pch = 19, col = "lightgrey", cex = .5)
points(centroid.scores_sp, pch = 19, col = col.pal)
ellipses <- ordiellipse(capscale_sp, coral.cov.site$clust.prepost, draw = "polygon", alpha = .5, col = col.pal, lty = 0)
arrows(x0 = centroid.scores_sp[(n.clust+1):(n.clust*2),1], y0 = centroid.scores_sp[(n.clust+1):(n.clust*2),2], 
       x1 = centroid.scores_sp[1:n.clust,1], y1 = centroid.scores_sp[1:n.clust,2], col = col.pal, length = .1, angle = 20, lwd = 4)
arrows(x0 = centroid.scores_sp[(n.clust+1):(n.clust*2),1], y0 = centroid.scores_sp[(n.clust+1):(n.clust*2),2], 
       x1 = centroid.scores_sp[1:n.clust,1], y1 = centroid.scores_sp[1:n.clust,2], col = "dimgrey", length = .1, angle = 20, lwd = 1)
plot(fit, col = "black")
legend(3.5, -1, legend = 1:n.clust, fill = col.pal, title = "Cluster")


#Export site and species scores
sites.scores_sp <- data.frame(site_code = coral.cov.site$site_code, 
                           clust.prepost = coral.cov.site$clust.prepost,
                           cluster = coral.cov.site$cluster,
                           pre_post = coral.cov.site$pre_post,
                           scores(capscale_sp, choices = c(1:2), display = "sites"),
                           scores(pco_sp, choices = c(1:2), display = "sites"))
names(sites.scores_sp)[7:8] <- c("PCO1","PCO2")


# Species contributions
species.scores_sp <- data.frame(label = names(coral.w.sp.site), scores(capscale_sp, choices = c(1,2), display = "species"))

sp.list <- data.frame(label = factor(ordiselect(coral.w.sp.site, capscale_sp, fitlim = .05))) %>% inner_join(species.scores_sp)

plot(capscale_sp, type="none", choice = c(1,2), display="species", xlim = c(-1.5,2.5))
points(species.scores_sp[,2:3], pch = 19, col = "grey", cex = .5)
points(sp.list[,2:3], pch = "+", col = "red", cex = .75)
#set.seed(314)
ordipointlabel(capscale_sp, display="species", select=sp.list$label, pch = 19, col = "red", cex = .8, add=T)



# Plot PCO
centroid.scores_sp_PCO <- sites.scores_sp %>% group_by(clust.prepost) %>% summarise(PCO1 = mean(PCO1), PCO2 = mean(PCO2)) %>%
  dplyr::select(PCO1,PCO2) %>% data.frame()

pco.fit <- envfit(pco, fit.preds, perm = 99, scaling = "sites", na.rm = T, choices = c(1,2))

pl_pco <- ordiplot(pco_sp, type = "none", xlim = c(-2,2.5), ylim = c(-2.5,3.5))
points(pl_pco, "sites", pch = 19, col = "lightgrey", cex = .5)
ellipses <- ordiellipse(pco_sp, sites.scores$clust.prepost, draw = "polygon", alpha = .5, col = col.pal, lty = 0)
points(centroid.scores_sp_PCO, pch = 19, col = col.pal)
arrows(x0 = centroid.scores_sp_PCO[(n.clust+1):(n.clust*2),1], y0 = centroid.scores_sp_PCO[(n.clust+1):(n.clust*2),2], 
       x1 = centroid.scores_sp_PCO[1:n.clust,1], y1 = centroid.scores_sp_PCO[1:n.clust,2], col = col.pal, length = .1, angle = 20, lwd = 4)
arrows(x0 = centroid.scores_sp_PCO[(n.clust+1):(n.clust*2),1], y0 = centroid.scores_sp_PCO[(n.clust+1):(n.clust*2),2], 
       x1 = centroid.scores_sp_PCO[1:n.clust,1], y1 = centroid.scores_sp_PCO[1:n.clust,2], col = "dimgrey", length = .1, angle = 20, lwd = 1)
plot(fit, col = "black")
legend(3, 2, legend = 1:n.clust, fill = col.pal, title = "Cluster")

# PCO species plot
pco.species.scores <- data.frame(label = names(coral.w.sp.site), scores(pco_sp, choices = c(1,2), display = "species"))
pco.sp.list <- data.frame(label = factor(ordiselect(coral.w.sp.site, pco_sp, fitlim = .05))) %>% inner_join(pco.species.scores)

plot(pco_sp, type="none", choice = c(1,2), display="species") # , xlim = c(-1.5,3.2)
points(pco.species.scores[,2:3], pch = 19, col = "grey", cex = .5)
points(pco.sp.list[,2:3], pch = "+", col = "red", cex = .75)
#set.seed(314)
ordipointlabel(pco_sp, display="species", select=pco.sp.list$label, pch = 19, col = "red", cex = .8, add=T)



# Temporal beta-diversity index  --------

coral.w.gff.site.pre <- coral.w.gff.site[coral.cov.site$pre_post == "Pre",] %>% data.frame()
coral.w.gff.site.post <- coral.w.gff.site[coral.cov.site$pre_post == "Post",] %>% data.frame()

coral.cov.site.pre <- coral.cov.site[coral.cov.site$pre_post == "Pre",]
coral.cov.site.post <- coral.cov.site[coral.cov.site$pre_post == "Post",]

coral.cover.pre <- rowSums(coral.w.sp.site.pre)
coral.cover.post <- rowSums(coral.w.sp.site.post)
coral.cover.delta <- coral.cover.post - coral.cover.pre
  
tbi <- TBI(coral.w.gff.site.pre, coral.w.gff.site.post)
tbi.pa <- TBI(coral.w.gff.site.pre, coral.w.gff.site.post, pa.tr = T)

tbi.dat <- data.frame(subset(coral.cov.site.post, select = c(site_code, longitude, latitude, cluster, maxDHW, maxSSTA, interval)),
                      TBI = tbi$TBI, 
                      TBI.p = tbi$p.TBI, 
                      TBI_PA = tbi.pa$TBI, 
                      TBI_PA.p = tbi.pa$p.TBI,
                      coral.cover.pre,
                      coral.cover.post,
                      coral.cover.delta)


tbi.dat <- tbi.dat %>% left_join(subset(sites.scores, select = c(site_code, CAP1, CAP2, PCO1, PCO2), pre_post == "Pre"))

richness.pre <- rowSums(ifelse(coral.w.gff.site.pre > 0, 1, 0))

#tbi.richness <- data.frame(TBI = tbi.dat$TBI, richness.pre)

summary(lm(TBI ~ richness.pre, data = tbi.richness))
plot(TBI ~ richness.pre, data = tbi.richness, xlab = "Coral richness (Pre)", ylab = "TBI", pch = 19)


# Test TBI significance based on permutations

tbi.test <- TBI(coral.w.gff.site.pre, coral.w.gff.site.post, nperm = 999, test.BC = TRUE, test.t.perm = TRUE)

TBI(coral.w.gff.site.pre[coral.cov.site.pre$cluster == 1,], coral.w.gff.site.post[coral.cov.site.post$cluster == 1,], nperm = 999, test.BC = TRUE, test.t.perm = TRUE)$t.test_B.C
TBI(coral.w.gff.site.pre[coral.cov.site.pre$cluster == 2,], coral.w.gff.site.post[coral.cov.site.post$cluster == 2,], nperm = 999, test.BC = TRUE, test.t.perm = TRUE)$t.test_B.C
TBI(coral.w.gff.site.pre[coral.cov.site.pre$cluster == 3,], coral.w.gff.site.post[coral.cov.site.post$cluster == 3,], nperm = 999, test.BC = TRUE, test.t.perm = TRUE)$t.test_B.C
TBI(coral.w.gff.site.pre[coral.cov.site.pre$cluster == 4,], coral.w.gff.site.post[coral.cov.site.post$cluster == 4,], nperm = 999, test.BC = TRUE, test.t.perm = TRUE)$t.test_B.C
TBI(coral.w.gff.site.pre[coral.cov.site.pre$cluster == 5,], coral.w.gff.site.post[coral.cov.site.post$cluster == 5,], nperm = 999, test.BC = TRUE, test.t.perm = TRUE)$t.test_B.C
TBI(coral.w.gff.site.pre[coral.cov.site.pre$cluster == 6,], coral.w.gff.site.post[coral.cov.site.post$cluster == 6,], nperm = 999, test.BC = TRUE, test.t.perm = TRUE)$t.test_B.C
TBI(coral.w.gff.site.pre[coral.cov.site.pre$cluster == 7,], coral.w.gff.site.post[coral.cov.site.post$cluster == 7,], nperm = 999, test.BC = TRUE, test.t.perm = TRUE)$t.test_B.C

# Exploratory analysis: relationships between TBI, delta CC, DHW, PCO1 and 2

tbi.dat2 <- tbi.dat

pairs.panels(tbi.dat2[,-(1:4)])

par(mfcol = c(6,1), mai = c(.1,.8,.2,.1))
boxplot(maxSSTA ~ cluster, data = tbi.dat2, ylab = "maxSSTA")
abline(h = mean(tbi.dat2$maxSSTA), lty = 2, col = "blue")
boxplot(maxDHW ~ cluster, data = tbi.dat2, ylab = "maxDHW")
abline(h = median(tbi.dat2$maxDHW), lty = 2, col = "blue")
boxplot(TBI ~ cluster, data = tbi.dat2, ylab = "TBI")
abline(h = mean(tbi.dat2$TBI), lty = 2, col = "blue")
boxplot(CAP.dist.pre.post ~ cluster, data = tbi.dat.dist, ylab = "CAP_dist")
abline(h = mean(tbi.dat.dist$CAP.dist.pre.post), lty = 2, col = "blue")
boxplot(CAP1 ~ cluster, data = tbi.dat2, ylab = "CAP1")
abline(h = mean(tbi.dat2$CAP1), lty = 2, col = "blue")
boxplot(CAP2 ~ cluster, data = tbi.dat2, ylab = "CAP2")
abline(h = mean(tbi.dat2$CAP2), lty = 2, col = "blue")


kruskal.test(maxSSTA ~ cluster, data = tbi.dat)
kruskal.test(maxDHW ~ cluster, data = tbi.dat)
kruskal.test(TBI ~ cluster, data = tbi.dat)
kruskal.test(CAP1 ~ cluster, data = tbi.dat)
kruskal.test(CAP2 ~ cluster, data = tbi.dat)

ggplot() +
  geom_point(data = tbi.dat, aes(x=maxDHW, y=coral.cover.delta, col = TBI))+
  geom_hline(yintercept = 0, col = "red")+
  facet_grid(cluster~.)

g.mm <- ggplot(data = tbi.dat, aes(x = abs(coral.cover.delta), y = TBI, group = cluster)) +
  facet_wrap(~cluster, ncol = 4) +
  geom_point() +
  stat_smooth(method = "lm", se = FALSE) +
  stat_cor(method="pearson", mapping = aes(group = cluster), label.y = 1)
g.mm

summary(lm(coral.cover.delta ~ maxDHW*PCO1 + maxDHW*PCO2, data = tbi.dat))
summary(lm(TBI ~ maxDHW*PCO1 + maxDHW*PCO2, data = tbi.dat))

summary(lm(coral.cover.delta ~ maxDHW*cluster, data = tbi.dat))
summary(lm(TBI ~ maxDHW*cluster, data = tbi.dat))


# Check that TBI correlates with pre/post distance on dbRDA

CAP.sites.scores.pre <- sites.scores %>% filter(pre_post == "Pre") %>% dplyr::select(CAP1, CAP2) %>% as.matrix()
CAP.sites.scores.post <- sites.scores %>% filter(pre_post == "Post") %>% dplyr::select(CAP1, CAP2) %>% as.matrix()

CAP.d <- pdist(CAP.sites.scores.pre, CAP.sites.scores.post)
CAP.dist.pre.post <- diag(as.matrix(CAP.d))
par(mfcol = c(1,1), mai = c(1,1,1,1))
plot(CAP.dist.pre.post, tbi.dat$TBI)

PCO.sites.scores.pre <- sites.scores %>% filter(pre_post == "Pre") %>% dplyr::select(PCO1, PCO2) %>% as.matrix()
PCO.sites.scores.post <- sites.scores %>% filter(pre_post == "Post") %>% dplyr::select(PCO1, PCO2) %>% as.matrix()

PCO.d <- pdist(PCO.sites.scores.pre, PCO.sites.scores.post)
PCO.dist.pre.post <- diag(as.matrix(PCO.d))
par(mfcol = c(1,1), mai = c(1,1,1,1))
plot(PCO.dist.pre.post, tbi.dat$TBI)

tbi.dat.dist <- data.frame(tbi.dat, CAP.dist.pre.post, PCO.dist.pre.post)

ggplot(data = tbi.dat.dist, aes(x=CAP.dist.pre.post, y=TBI))+
  geom_point()+
  stat_smooth(method = "lm", se = FALSE) +
  stat_cor(method="pearson")

# ggplot(data = tbi.dat.dist, aes(x=PCO.dist.pre.post, y=TBI))+
#   geom_point()+
#   stat_smooth(method = "lm", se = FALSE) +
#   stat_cor(method="pearson")
# 
# ggplot(data = tbi.dat.dist, aes(x=coral.cover.delta, y=TBI))+
#   geom_point()+
#   stat_smooth(method = "lm", se = FALSE) +
#   stat_cor(method="pearson")

ggplot(data = tbi.dat.dist, aes(x=abs(coral.cover.delta), y=TBI))+
  geom_point()+
  stat_smooth(method = "lm", se = FALSE) +
  stat_cor(method="pearson")

summary(lm(TBI ~ CAP.dist.pre.post, data = tbi.dat.dist))
summary(lm(TBI ~ PCO.dist.pre.post, data = tbi.dat.dist))
summary(lm(TBI ~ coral.cover.delta, data = tbi.dat.dist))

# Stats of TBI and CAP.dist.pre.post by cluster

stats <- tbi.dat.dist %>% group_by(cluster) %>% summarize(mean_dist = mean(CAP.dist.pre.post),
                                                          sd_dist = sd(CAP.dist.pre.post),
                                                          mean_TBI = mean(TBI),
                                                          sd_TBI = sd(TBI))

write.csv(data.frame(stats), file = "cluster_stats_TBI.csv", row.names = F)

# Test which species changed over time in each cluster - Problem with SIMPER, fix with t-tests? -------


for (i in c(1:7)) {
  
  t <- tpaired.krandtest(sqrt(coral.w.sp.site.pre[tbi.dat$cluster == i,]), sqrt(coral.w.sp.site.post[tbi.dat$cluster == i,]))
  #t <- tpaired.krandtest(sqrt(mat1[bent_w.site.1$clust.bent == i,]), sqrt(mat2[bent_w.site.2$clust.bent == i,]))
  
  t$t.tests[,c(1,2,6)] <- t$t.tests[,c(1,2,6)]*(-1)
  names(t$t.tests)[1] <- "mean(T2-T1)"
  names(t$t.tests)[6] <- "sign(T2-T1)"
  
  #View(t$t.tests)
  print(paste("Clust No", i, sep="_"))
  print(mean(tbi.dat$coral.cover.delta[tbi.dat$cluster == i]))
  print(sd(tbi.dat$coral.cover.delta[tbi.dat$cluster == i]))
  print(paste("N = ", dim(tbi.dat[tbi.dat$cluster == i,])[1], sep = ""))
  print(tbi.dat$site_code[tbi.dat$cluster == i])
  
  # SIMPER Pre vs. post 2016 (and reformat output)
  simper.dat <- sqrt(rbind(coral.w.sp.site.pre[tbi.dat$cluster == i,], coral.w.sp.site.post[tbi.dat$cluster == i,]))
  simper.dat <- subset(simper.dat, select=colSums(simper.dat)!=0)
  simper.group <- c(rep("Pre", dim(coral.w.sp.site.pre[tbi.dat$cluster == i,])[1]), rep("Post", dim(coral.w.sp.site.post[tbi.dat$cluster == i,])[1]))
  simper <- vegan::simper(as.matrix(simper.dat), group = simper.group)
  simper_species <- data.frame(Cat = row.names(summary(simper, ordered = T)$Pre_Post), 
                               summary(simper, ordered = T)$Pre_Post) #%>% 
  #filter(cumsum < .9)
  
  simper_species[,-1] <- round(simper_species[,-1], 3)
  simper_species$meanPre <- round((simper_species$ava),2)
  simper_species$meanPost <- round((simper_species$avb),2)
  simper_species$meanDiff <- simper_species$meanPost - simper_species$meanPre
  
  # Filter SIMPER species with significant change over time
  simper_species_sub <- simper_species %>% filter(Cat %in% row.names(t$t.tests)[t$t.tests$p.perm < .05]) %>%
    dplyr::select(CategoryName = Cat, average_contrib = average, sd_contrib = sd, meanPre, meanPost, meanDiff)
  
  # Add pre-disturbance indicators
  simper_species_sub <- data.frame(simper_species_sub[order(simper_species_sub$meanDiff, decreasing = T),], pre.Indic = 0)
  simper_species_sub$pre.Indic <- ifelse(simper_species_sub$CategoryName %in% row.names(ind.Coral.tb)[ind.Coral.tb$cluster == i], 1, 0)
  
  write.csv(simper_species_sub, paste("SIMPER-2/simper_cluster_", i, ".csv", sep = ""), row.names = F)
}

# Export pre-disturbance indicators
write.csv(ind.Coral.tb, "SIMPER/indicator_coral_categories.csv")

# Test which GFF changed over time in each cluster --------

coral.w.gff.site.pre <- coral.w.gff.site[coral.cov.site$pre_post == "Pre",] %>% data.frame()
coral.w.gff.site.post <- coral.w.gff.site[coral.cov.site$pre_post == "Post",] %>% data.frame()

for (i in c(1:7)) {
  
  t <- tpaired.krandtest(sqrt(coral.w.gff.site.pre[tbi.dat$cluster == i,]), sqrt(coral.w.gff.site.post[tbi.dat$cluster == i,]))
  #t <- tpaired.krandtest(sqrt(mat1[bent_w.site.1$clust.bent == i,]), sqrt(mat2[bent_w.site.2$clust.bent == i,]))
  
  t$t.tests[,c(1,2,6)] <- t$t.tests[,c(1,2,6)]*(-1)
  names(t$t.tests)[1] <- "mean(T2-T1)"
  names(t$t.tests)[6] <- "sign(T2-T1)"
  
  ttest.res <- data.frame(GenusGF = row.names(t$t.tests), t$t.tests) %>% filter(p.perm < 0.05) %>% data.frame()
  ttest.res <- ttest.res[order(ttest.res$p.perm, decreasing = F),]

  write.csv(ttest.res, paste("Ttests-2/ttest_cluster_", i, ".csv", sep = ""), row.names = F)
}

# Global t-test

t <- tpaired.krandtest(sqrt(coral.w.gff.site.pre), sqrt(coral.w.gff.site.post))

t$t.tests[,c(1,2,6)] <- t$t.tests[,c(1,2,6)]*(-1)
names(t$t.tests)[1] <- "mean(T2-T1)"
names(t$t.tests)[6] <- "sign(T2-T1)"

ttest.res <- data.frame(GenusGF = row.names(t$t.tests), t$t.tests) %>% filter(p.perm < 0.05) %>% data.frame()
ttest.res <- ttest.res[order(ttest.res$p.perm, decreasing = F),]

write.csv(ttest.res, "Ttests-2/ttest_global.csv", row.names = F)

# Coral cover change without zeros
coral.w.gff.site.pre.NA <- coral.w.gff.site.pre
coral.w.gff.site.pre.NA[coral.w.gff.site.pre == 0] <- NA

coral.w.gff.site.post.NA <- coral.w.gff.site.post
coral.w.gff.site.post.NA[coral.w.gff.site.post == 0] <- NA

mean.pre <- colMeans(coral.w.gff.site.pre.NA, na.rm = T)
mean.post <- colMeans(coral.w.gff.site.post.NA, na.rm = T)

range(mean.post - mean.pre, na.rm = T)
diff <- data.frame(cover.pre = mean.pre,
                   cover.post = mean.post,
                   diff = mean.post - mean.pre)

diff.255 <- data.frame(cover.pre = t(coral.w.gff.site.pre[coral.cov.site.pre$site_code == "GBR255",]),
                       cover.post = t(coral.w.gff.site.post[coral.cov.site.post$site_code == "GBR255",]))

# Maximum change in cover
ttest.res <- read.csv("Ttests-2/ttest_global.csv")

diff.all <- coral.w.gff.site.post - coral.w.gff.site.pre

diff.min <- diff.all %>% summarize_all(min)
diff.max <- diff.all %>% summarize_all(max)
diff.mean <- diff.all %>% summarize_all(mean)

diff <- data.frame(min = t(diff.min), max = t(diff.max), mean = t(diff.mean))

diff.signif <- diff[row.names(diff) %in% ttest.res$GenusGF,]
diff.signif <- data.frame(GenusGF = row.names(diff.signif), diff.signif)

ttest.res2 <- ttest.res %>% left_join(diff.signif)

write.csv(ttest.res2, "Ttests-2/ttest_global_withMinMax.csv", row.names = F)
write.csv(diff, "Ttests-2/diff_all.csv", row.names = T)

# Test if indicator Genus GF lost their indicator value after heatwave --------

coral.w.gff.site.post <- coral.w.gff.site[coral.cov.site$pre_post == "Post",] %>% data.frame()

ind.Coral.post <- indval(subset(coral.w.gff.site.post, select = colSums(coral.w.gff.site.post)>0), coral.cov.site.post$cluster, numitr=100)
ind.Coral.post.tb <- summary_indval(ind.Coral.post)
#write.csv(ind.Coral.post.tb, "MRT_Indicator_Corals_post.csv", row.names = T)

table(row.names(ind.Coral.tb) %in% row.names(ind.Coral.post.tb))
row.names(ind.Coral.tb)[!row.names(ind.Coral.tb) %in% row.names(ind.Coral.post.tb)]

row.names(ind.Coral.post.tb)[!row.names(ind.Coral.post.tb) %in% row.names(ind.Coral.post.tb)]





# Boosted regression tree: relative importance of initial community composition (CAP1, CAP2) vs thermal stress (DHW, SST) in explaining TBI ------

HumanPop <- coral.cov.site %>% group_by(site_code) %>% summarize(HumanPop = mean(pop2015_20km))
#MPA <- bent_CM %>% group_by(SiteCode, MPA) %>% summarise()

tbi.dat <- tbi.dat %>% left_join(HumanPop) #%>% inner_join(MPA)
tbi.dat$HumanPop <- factor(ifelse(tbi.dat$HumanPop > 0, 1, 0))

tbi.dat$cluster <- factor(tbi.dat$cluster)

# GBM on TBI
gbm.TBI <- gbm.step(data = tbi.dat,
                    gbm.x = c(4,5,6,17,18,19),
                    gbm.y = 8,
                    family = "gaussian",
                    tree.complexity = 2,
                    learning.rate = 0.001,
                    bag.fraction = 0.7)

gbm.TBI.CAP <- gbm.step(data = tbi.dat,
                    gbm.x = c(4,5,6,15,16),
                    gbm.y = 8,
                    family = "gaussian",
                    tree.complexity = 2,
                    learning.rate = 0.001,
                    bag.fraction = 0.7)

# Explained deviance = 39.4%%
par(mai = c(.5,.6,.1,.2))
gbm.plot(gbm.TBI.CAP, write.title=F, n.plots = 5, plot.layout = c(5,1))
summary(gbm.TBI.CAP)


dev.off()
qqnorm(gbm.TBI.CAP$residuals, pch=19)
qqline(gbm.TBI.CAP$residuals)
# find.int <- gbm.interactions(gbm.TBI)
# find.int$rank.list
# gbm.perspec(gbm.TBI,4,3)

# GBM on TBI_PA
gbm.TBI.PA <- gbm.step(data = tbi.dat,
                       gbm.x = c(4,5,6,15,16,19),
                       gbm.y = 10,
                       family = "gaussian",
                       tree.complexity = 2,
                       learning.rate = 0.001,
                       bag.fraction = 0.7)

# Explained deviance = 30.4%%
par(mai = c(.6,.6,.2,.6))
gbm.plot(gbm.TBI.PA, write.title=F, plot.layout = c(6,1))
summary(gbm.TBI.PA)


qqnorm(gbm.TBI.PA$residuals, pch=19)
qqline(gbm.TBI.PA$residuals)
# find.int <- gbm.interactions(gbm.TBI)
# find.int$rank.list
# gbm.perspec(gbm.TBI,4,3)

# GBM on HC
gbm.HC <- gbm.step(data = tbi.dat,
                   gbm.x = c(4,5,6,15,16),
                   gbm.y = 14,
                   family = "gaussian",
                   tree.complexity = 2,
                   learning.rate = 0.001,
                   bag.fraction = 0.7)

# Explained deviance = 44.5%
par(mai = c(.5,.6,.1,.2))
gbm.plot(gbm.HC, write.title=F, n.plots = 5, plot.layout = c(5,1))
summary(gbm.HC)

dev.off()
qqnorm(gbm.HC$residuals, pch=19)
qqline(gbm.HC$residuals)

# find.int <- gbm.interactions(gbm.HC)
# find.int$rank.list
# gbm.perspec(gbm.HC, 3, 2, x.range = c(0,25), y.range = c(0,7), z.range = c(-25,5))
# gbm.HC.simpl <- gbm.simplify(gbm.HC, n.drops = 3)


# Paired boxplots of HCC pre. vs post for each cluster --------


HCC.pre.post <- data.frame(subset(coral.cov.site, select = c(site_code, cluster, pre_post)), LiveCoralCover = rowSums(coral.w.gff.site))
HCC.pre.post$PrePost <- factor(HCC.pre.post$pre_post, levels = c("Pre", "Post"))
HCC.pre.post$cluster <- factor(HCC.pre.post$cluster)

g <- ggplot(HCC.pre.post, aes(PrePost, LiveCoralCover, factor(cluster))) +
  facet_grid(.~ factor(cluster)) +
  geom_boxplot(aes(colour = factor(cluster)), width=0.3, size=1.5, fatten=1.2) +
  geom_point(aes(colour=factor(cluster)), size=2, alpha=0.5) +
  geom_line(aes(x = PrePost, y = LiveCoralCover, group = site_code, colour=factor(cluster)), linetype="11", alpha = .5) +
  geom_boxplot(aes(colour = factor(cluster)), width=0.3, size=1.5, fatten=1.5) +
  scale_color_brewer(palette = "Paired")+
  theme_classic() + theme(legend.position = "none")

g

HCC.pre.post.summary <- HCC.pre.post %>% mutate(pre_post = NULL) %>% pivot_wider(names_from = PrePost, values_from = LiveCoralCover) %>%
  mutate(Diff = Post - Pre) %>%
  group_by(cluster) %>% summarize(mean.Diff = mean(Diff, na.rm = T), sd.Diff = sd(Diff, na.rm = T))

HCC.pre.post.test <- HCC.pre.post %>% mutate(pre_post = NULL) %>% pivot_wider(names_from = PrePost, values_from = LiveCoralCover) %>%
  mutate(Diff = Post - Pre) %>%
  group_by(cluster) %>% summarize(P = t.test(Diff)$p.value)

tbi.dat %>% group_by(cluster) %>% summarize(mean_tbi = mean(TBI))


# List of surveys (ProcB revisions) ----

coral.survey.ls <- bent.dat %>% dplyr::select(survey_id, site_code, site_name, latitude, longitude,
                                              survey_date, program, database, is_NESP_pre.post, survey_year) %>%
                                filter(site_code %in% coral.cov.site.post$site_code) %>%
                                group_by(survey_id, site_code, site_name, latitude, longitude,
                                         survey_date, program, database, is_NESP_pre.post, survey_year) %>%
                                summarize()

write.csv(coral.survey.ls, "coral_survey_ls.csv", row.names = F)

original.survey.ls <- read.csv("Ningaloo_NW_GBR_SurveyList.csv")

discarded.survey.ls <- original.survey.ls[!(original.survey.ls$SurveyID %in% coral.survey.ls$survey_id),]

table(discarded.survey.ls$SurveyID %in% bent.dat$survey_id) #403 out of 485 missing surveys from original list are not in PQ data

# New Analysis: biplot of species-level changes in cover and no. occupied sites ----------

ini.no.sites <- coral.w.gff.site.pre
ini.no.sites[coral.w.gff.site.pre > 0] <- 1
ini.no.sites <- colSums(ini.no.sites)

fin.no.sites <- coral.w.gff.site.post
fin.no.sites[coral.w.gff.site.post > 0] <- 1
fin.no.sites <- colSums(fin.no.sites)

delta.no.sites <- fin.no.sites - ini.no.sites

delta.cover <- coral.w.gff.site.post - coral.w.gff.site.pre
delta.cover.no.zero <- delta.cover
delta.cover.no.zero[delta.cover == 0] <- NA
delta.cover <- colMeans(delta.cover.no.zero, na.rm = T)

data.win.los <- data.frame(ini.no.sites, delta.no.sites, delta.cover)

plot.win.los <- ggplot() +
  geom_point(data = data.win.los, aes(x=delta.no.sites, y=delta.cover, size = ini.no.sites), col = "blue") +
  geom_text(data = data.win.los, aes(x=delta.no.sites, y=delta.cover, label=row.names(data.win.los))) +
  theme_light() +
  ylim(-2,2) + xlim(-30,30) +
  geom_hline(yintercept = 0, col = "black")+
  geom_vline(xintercept = 0)



