
##############################################################
# NESP CORAL ANALYSIS - REVISED - CM 29/11/21                #
##############################################################

# Load libraries ------------
rm(list = ls())

library(RLSMetrics)
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
library(mvpart)
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
library(tidyverse)

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

setwd("/Users/Camille/Dropbox/My documents/Projects/Future Fellowship/RESEARCH/NESP-coral-analysis-revised")
load("NESP-coral-analysis-revised.RData")
source("brt.functions.R")

# Load covariates
#dhw.raster <- raster::raster("Marine_heatwaves_Camille/Data/cropped_dhw_max_2016.grd")
#ssta.raster <- raster::raster("Marine_heatwaves_Camille/Data/cropped_ssta_max_2016.grd")

SiteCov_msec <- read.csv("MSEC/msec_out_npp_sst_wave.csv")
SiteCov_msec_rd <- read.csv("MSEC/msec_out_ReefLandArea_HumanPop.csv")
SiteCov <- SiteCov_msec %>% inner_join(SiteCov_msec_rd) %>% mutate(SiteCode = NULL, Region = NULL)

aus <- importShapefile("shapefiles/250K_coastline", readDBF=FALSE)
aus2 <- aus %>% dplyr::select(group=PID, POS=POS,long=X,lat=Y)
aus.sf <- st_read("shapefiles/250K_coastline.shp")

reef <- importShapefile("shapefiles/SDE_OWNER_crcgis_land250", readDBF=FALSE)

# Fix shapefile (takes forever!)
#shape.clip <- closePolys(clipPolys(reef, xlim=c(142,155), ylim=c(-25,-10)))
#reefs <- as.PolySet(subset(shape.clip, PID != 2401), projection = "LL")
#reef2 <- reefs %>% dplyr::select(group=PID, POS=POS,long=X,lat=Y)

# Load temperature anomalies
SiteSSTA <- read.csv("Marine_heatwaves_Camille/sites_with_sst.csv")
SiteSSTA$MMM <- apply(subset(SiteSSTA, select = Jan_SST_cl:Dec_SST_cl), 1, max)

SiteSSTA <- SiteSSTA %>% dplyr::select(SiteLat, SiteLong, MMM,
                          maxDHW, maxSSTA, Mean_SST_cl_noaa = Mean_SST_cl, SD_SST_cl_noaa = SD_SST_cl)

# Load benthic data

# Compile survey list from old extract and export
bent.dat.old <- read.csv("NESP-Extract_2021-11-09.csv")

bent.dat.old$survey_date <- as.Date(bent.dat.old$survey_date, tryFormats = "%d/%m/%y")
bent.dat.old$survey_year <- format(bent.dat.old$survey_date, '%Y')

bent.dat.old$survey_id <- factor(bent.dat.old$survey_id)
bent.dat.old$dataset_id <- factor(bent.dat.old$dataset_id)

bent.dat.survey.ls <- bent.dat.old %>% group_by(survey_id) %>% summarize()
write.csv(bent.dat.survey.ls, "bent.dat.survey.ls.csv", row.names = F)

# Import new extract and filter with survey list
bent.dat <- read.csv("PQ_FullRes_infilled.csv")

bent.dat$survey_date <- as.Date(bent.dat$survey_date, tryFormats = "%d/%m/%y")
bent.dat$survey_year <- format(bent.dat$survey_date, '%Y')

bent.dat$survey_id <- factor(bent.dat$survey_id)
bent.dat$dataset_id <- factor(bent.dat$dataset_id)

bent.dat <- bent.dat %>% filter(survey_id %in% bent.dat.survey.ls$survey_id)


# Remove duplicate surveys
dataset_id_duplicates <- c("912346798_RLS Australian Coral Species List_982", 
                           "912351889_RLS Catalogue_946", 
                           "912352031_RLS Catalogue_946",
                           "912352656_RLS Catalogue_939",
                           "912352687_RLS Catalogue_928")

bent.dat <- bent.dat %>% filter(!dataset_id %in% dataset_id_duplicates)

# Add 'unfinished' missing info
# bent.dat$is_finished[!(bent.dat$is_finished %in% c('yes', 'no'))] <- 'no'


# Metadata --------------

# Check unique survey_id/dataset_id match (i.e. each transect is scored only once)
bent.meta <- bent.dat %>% group_by(survey_id, dataset_id) %>% summarize() %>%
                          group_by(survey_id) %>% summarize(count = n_distinct(dataset_id))

# Check number of finishes vs unfinished surveys
bent.meta <- bent.dat %>% group_by(survey_id, is_finished, is.infilled.in.previous.extract) %>% summarize()


# Check labeller bias in non-coral benthic scores --------------

bent.dat.check <- bent.dat %>% dplyr::select(survey_id, location, site_code, site_name, latitude, longitude, survey_date, survey_year,
                                             annotation_labeller_username, label, flag_dead, flag_bleached, RLS_category, RLE_category,
                                             percent_cover, is_finished, is_NESP_pre.post, is.infilled.in.previous.extract) %>%
                                filter(is_finished == "yes" & RLE_category %in% c("0", "Kelp"))

# Check number of labellers per survey_id
bent.dat.check2 <- bent.dat.check %>% group_by(survey_id, annotation_labeller_username) %>% summarize() %>%
                                      group_by(survey_id) %>% summarize(n_labellers = n_distinct(annotation_labeller_username))

# Wide format for single-labeller surveys only
bent.dat.check.wide <- bent.dat.check %>% filter(survey_id %in% bent.dat.check2$survey_id[bent.dat.check2$n_labellers == 1]) %>%
                        group_by(survey_id, location, survey_year, is_NESP_pre.post, labeller = annotation_labeller_username, label) %>%
                        summarize(percent_cover = sum(percent_cover)) %>%
                        pivot_wider(names_from = label, values_from = percent_cover, values_fill = 0) %>% data.frame

# Divide by total cover of non-coral categories
bent.dat.check.wide$total <- rowSums(base::subset(bent.dat.check.wide, select = Bare.Rock:Colonial.Anemones..Zoanthids.and.Corallimorphs))

bent.dat.check.wide[, 6:44] <- bent.dat.check.wide[, 6:44]/bent.dat.check.wide$total

bent.dat.check.wide$total <- rowSums(base::subset(bent.dat.check.wide, select = Bare.Rock:Colonial.Anemones..Zoanthids.and.Corallimorphs))

bent.dat.check.wide <- bent.dat.check.wide %>% filter(!is.na(total))


# CAPSCALE testing for labeller effect

capscale.dat <- sqrt(subset(bent.dat.check.wide, select = Bare.Rock:Colonial.Anemones..Zoanthids.and.Corallimorphs))

capscale <- capscale(capscale.dat ~ labeller + location + survey_year + is_NESP_pre.post, data = bent.dat.check.wide, distance = "bray", add = T)

capscale
plot(capscale)
anova(capscale, permutations = how(nperm=  999), by = "term")

boxplot(Turfing.algae...2.cm.high.algal.sediment.mat.on.rock. ~ labeller, data = bent.dat.check.wide)

# Check which benthic categories are sensitive to labeller effect: Kruskal-Wallis test + boxplot
kw.pvalues <- data.frame(label = names(bent.dat.check.wide)[6:44], kw.pvalues = NA)

for (i in 6:44) {
  kw.test <- kruskal.test(bent.dat.check.wide[,i] ~ bent.dat.check.wide$labeller)
  kw.pvalues$kw.pvalues[i-5] <- ifelse(kw.test$p.value < 0.001, 0, 0.05)
}

par(mfcol = c(5,2), mai = c(.3,.6,.2,.2))
boxplot(Bare.Rock ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Bare.Rock), col = "red")

boxplot(Coral.rubble ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Coral.rubble), col = "red")

boxplot(Coral.rubble.with.turf.encrusting.algae ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Coral.rubble.with.turf.encrusting.algae), col = "red")

boxplot(Turfing.algae...2.cm.high.algal.sediment.mat.on.rock. ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Turfing.algae...2.cm.high.algal.sediment.mat.on.rock.), col = "red")

boxplot(Crustose.coralline.algae ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Crustose.coralline.algae), col = "red")

boxplot(Sponges..encrusting. ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Sponges..encrusting.), col = "red")

boxplot(Encrusting.leathery.algae ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Encrusting.leathery.algae), col = "red")

boxplot(Geniculate.coralline.algae ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Geniculate.coralline.algae), col = "red")

boxplot(Slime..not.trapping.sediment. ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Slime..not.trapping.sediment.), col = "red")

boxplot(Filamentous.brown.algae_epiphyte ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Filamentous.brown.algae_epiphyte), col = "red")




# Build coral-only dataset (wide format) ------------


coral.dat$survey_year <- as.numeric(coral.dat$survey_year)

coral.dat <- bent.dat %>% filter(RLE_category == "Live Coral") %>%
              dplyr::select(survey_id, area, location, site_code, site_name, latitude, longitude, depth,
                            survey_date, survey_year, pre_post = is_NESP_pre.post,
                            RLS_category, label, flag_bleached, percent_cover) 

coral.dat$RLS_category_2 <- sub("_Bleached", coral.dat$RLS_category, replacement = "")

# Calculate total live coral cover (total_live) and total bleached coral cover (total_bleached)
coral.total <- coral.dat %>% group_by(survey_id) %>%
                              summarize(total_live = sum(percent_cover),
                                        total_bleached = sum(percent_cover[flag_bleached == TRUE]))

# Match to coral.total and switch to wide format
coral.w <- coral.dat %>% left_join(coral.total) %>% 
                          group_by(survey_id, area, location, site_code, site_name, latitude, longitude, depth,
                                   survey_date, survey_year, pre_post, total_live, total_bleached, label) %>%
                          summarize(percent_cover = sum(percent_cover)) %>%
                          pivot_wider(names_from = label, values_from = percent_cover, values_fill = 0) %>%
                          data.frame()

coral.w[,12:413] <- round(coral.w[,12:413],2)            
table(round(rowSums(coral.w[,14:413])) == round(coral.w$total_live))

# Map coral labels to RLS categories
coral.map <- coral.dat %>% group_by(survey_id, RLS_category_2, label) %>%
                            summarize(percent_cover = sum(percent_cover)) %>%
                            group_by(RLS_category_2, label) %>% 
                            summarize(mean_cover_present = mean(percent_cover[percent_cover > 0]))


# Match coral dataset to covariates --------------

# MSEC data
coral.cov <- coral.w %>% dplyr::select(survey_id, survey_year, pre_post, site_code, site_name, latitude, longitude) %>%
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

# Lump coral categories based on frequency and mean cover where present (i.e. reduce number of labels from 400 to 106) ---------

# Calculate species frequencies and match to coral mapping
spp.freq <- coral.dat %>% group_by(survey_id, label) %>% summarize(cover = sum(percent_cover))
spp.freq$count <- ifelse(spp.freq$cover > 0, 1, 0)
spp.freq <- spp.freq %>% group_by(label) %>% summarize(freq = sum(count))
spp.freq <- spp.freq %>% left_join(coral.map) %>% dplyr::select(RLS_category_2, label, freq, mean_cover_present)

# Frequency rank
spp.freq <- spp.freq[order(spp.freq$freq, decreasing = T),]
spp.freq$rank.freq <- 1:dim(spp.freq)[1]

# Abundance (cover) rank - conditional on presence
spp.freq <- spp.freq[order(spp.freq$mean_cover_present, decreasing = T),]
spp.freq$rank.cover <- 1:dim(spp.freq)[1]

# Plots
with(spp.freq, plot(rank.freq, rank.cover)) # congruence between ranking systems
with(spp.freq, plot(rank.freq, freq))   # rank frequency curve
with(spp.freq, plot(rank.cover, log10(mean_cover_present+1)))   # rank abundance curve

# 50 species are present on at least 5% of all transects
# 57 species have at least 3% cover where present

# Lumping less common or frequent species within RLS_category_2
spp.freq$new.label <- with(spp.freq, ifelse(rank.cover < 57 | rank.freq < 50, label, RLS_category_2))
n_distinct(spp.freq$new.label) # 106 unique new coral categories (reduced from 400)

with(spp.freq, table(label == RLS_category_2)) # 15 corals were already assigned to their broader category 
with(spp.freq, table(new.label == label)) # 103 corals kept their original label
with(spp.freq, table(new.label == RLS_category_2)) # 312 additional corals now assigned their broader category

with(spp.freq, table(new.label == RLS_category_2 & new.label == label))

# Export for Emre
write.csv(spp.freq, "EmreExports/coral_species_mapping.csv", row.names = F)


# Rebuild coral.w with lumped categories
new.label <- spp.freq %>% dplyr::select(label, new.label)

coral.w <- coral.dat %>% left_join(coral.total) %>% left_join(new.label) %>%
  group_by(survey_id, area, location, site_code, site_name, latitude, longitude, depth,
           survey_date, survey_year, pre_post, total_live, total_bleached, new.label) %>%
  summarize(percent_cover = sum(percent_cover)) %>%
  pivot_wider(names_from = new.label, values_from = percent_cover, values_fill = 0) %>%
  data.frame()


# Coral MRT for pre-bleaching surveys at the site level ----------------


coral.w.mrt <- coral.w %>% filter(pre_post == "Pre") %>% group_by(site_code) %>%
                            summarize_at(vars("Branching.Acropora":"Acropora.sukarnoi"), mean) %>% 
                            mutate(site_code = NULL) %>% data.frame()
coral.w.mrt <- coral.w.mrt[,colSums(coral.w.mrt)>0]
coral.w.mrt.hel <- decostand(coral.w.mrt, "hellinger")

coral.cov.mrt <- coral.cov %>% filter(pre_post == "Pre") %>% group_by(site_code, site_name) %>% 
                              summarize_at(vars("latitude":"interval"), mean) %>% data.frame()

coral.cov.mrt.sub <- coral.cov.mrt %>% dplyr::select(Mean_SST_cl, SD_SST_cl, npp_mean, npp_sd, wave_mean, land_area_20km,
                                                    reef_area_20km, pop2015_20km, depth, MMM) %>% data.frame()
                        

coral.mrt <- mvpart(as.matrix(coral.w.mrt.hel)~., coral.cov.mrt.sub, xv=c("pick"), legend=F)
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


# Indicators species
ind.Coral <- indval(coral.w.mrt, coral.clusters, numitr=100)
ind.Coral.tb <- summary_indval(ind.Coral)

all.indval <- data.frame(ind.Coral[["indval"]], ind.Coral[["pval"]])
#write.csv(ind.Coral.tb, "MRT_Indicator_Corals.csv", row.names = T)
#write.csv(all.indval, "MRT_Indicator_Corals_all.csv", row.names = T)


# Map communities
coords <- coral.cov.mrt %>% dplyr::select(site_code, longitude, latitude) %>% mutate(cluster = coral.clusters)

ggplot() +
  geom_polygon(data=reef2, aes(long, lat, group=group),  fill="darkgray", color="darkgray") +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  #coord_map(xlim=c(118,128), ylim=c(-20,-11)) +
  xlab(expression(paste(Longitude^o, ~'E'))) +
  ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_point(data = coords, aes(x=longitude, y=latitude, fill = factor(cluster)), size=4, shape=21) +
  scale_fill_brewer(palette = "Paired", name = "Cluster")+
  theme_light() +
  theme(legend.position = "none") +
  #theme(legend.position = c(.85,.7), legend.direction = "vertical")+
  #ggtitle("CategoryName MRT") +
  #coord_cartesian(xlim = c(142, 157), ylim = c(-25,-9)) # GBR
  coord_cartesian(xlim = c(112, 128), ylim = c(-24,-12)) # NWS

ggplot() +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="black", color="black") +
  geom_rect(data = data.frame(xmin = 142, xmax = 157, ymin = -25, ymax = -9), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), col = "red", fill = "transparent", lwd = 3)+
  geom_rect(data = data.frame(xmin = 112, xmax = 128, ymin = -24, ymax = -12), aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), col = "red", fill = "transparent", lwd = 3)+
  theme_void()

# Export for Emre

data_mrt <- data.frame(coords, coral.w.mrt)
write.csv(data_mrt, "EmreExports/coral_data_mrt.csv", row.names = F)
write.csv(coral.cov.mrt, "EmreExports/covariate_data_mrt.csv", row.names = F)

# Distance-based RDA ----------------

# Match cluster with covariates and compute site-level averages
coral.cov <- coords %>% right_join(coral.cov)
coral.cov.site <- coral.cov %>% group_by(site_code, pre_post, longitude, latitude, cluster, site_name) %>%
                                summarize_at(vars("Mean_SST_cl":"interval"), mean)       

coral.cov.site$clust.prepost <- with(coral.cov.site, paste(pre_post, cluster, sep = "_"))

boxplot(maxSSTA ~ cluster, data = coral.cov.site[coral.cov.site$pre_post == "Post",])
boxplot(maxDHW ~ cluster, data = coral.cov.site[coral.cov.site$pre_post == "Post",])

# Compute site-level benthic covers
coral.w.site <- coral.w %>% group_by(site_code, pre_post) %>%
                            summarize_at(vars("Branching.Acropora":"Acropora.sukarnoi"), mean) %>%
                            ungroup() %>%
                            mutate(site_code = NULL, pre_post = NULL)
  
coral.w.site <- subset(coral.w.site, select=colSums(coral.w.site)>0)

# Remove sites surveyed only once
coral.w.site <- coral.w.site[!is.na(coral.cov.site$interval),]
coral.cov.site <- coral.cov.site[!is.na(coral.cov.site$interval),]

# Export for Emre
write.csv(coral.w.site, "EmreExports/coral.data.site.prepost.csv", row.names = F)
write.csv(coral.cov.site, "EmreExports/covariate.data.site.prepost.csv", row.names = F)

# CAPSCALE

# capscale.preds <- data.frame(clust.bentCat = bentCat_w_site$clust.bentCat, PrePost = bentCat_w_site$PrePost, clust.prepost,
#                              maxSSTA = ssta_site$maxSSTA, maxDHW = ssta_site$maxDHW,
#                              resurveyed = sampl_years$resurveyed, from2016 = sampl_years$from2016)


capscale <- capscale(sqrt(coral.w.site) ~ clust.prepost + Condition(interval), data = coral.cov.site, distance = "bray", add = T)
pco <- capscale(sqrt(coral.w.site) ~ 1, data = coral.cov.site, distance = "bray", add = T)

capscale
plot(capscale)
anova(capscale, permutations = how(nperm=99), by = "term")

fit.preds <- data.frame(subset(coral.cov.site, select = c(maxDHW, maxSSTA)), LiveCoralCover = rowSums(coral.w.site))
fit <- envfit(capscale, fit.preds, perm = 0, display = "lc", scaling = "sites", na.rm = T)

centroid.scores <- data.frame(scores(capscale, choices = c(1,2))$centroids)
n.clust <- as.numeric(dim(centroid.scores)[1]/2)
col.pal <- rep(brewer.pal(n.clust, "Paired"),2)

pl <- ordiplot(capscale, type = "none", xlim = c(-4,4))
points(pl, "sites", pch = 19, col = "lightgrey", cex = .5)
points(centroid.scores, pch = 19, col = col.pal)
ellipses <- ordiellipse(capscale, coral.cov.site$clust.prepost, draw = "polygon", alpha = .5, col = col.pal, lty = 0)
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


# Species contributions
species.scores <- data.frame(label = names(coral.w.site), scores(capscale, choices = c(1,2), display = "species"))

sp.list <- data.frame(label = factor(ordiselect(coral.w.site, capscale, fitlim = .1))) %>% inner_join(species.scores)

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
plot(fit, col = "black")
legend(3, 2, legend = 1:n.clust, fill = col.pal, title = "Cluster")

# PCO species plot
pco.species.scores <- data.frame(label = names(coral.w.site), scores(pco, choices = c(1,2), display = "species"))
pco.sp.list <- data.frame(label = factor(ordiselect(coral.w.site, pco, fitlim = .1))) %>% inner_join(pco.species.scores)

plot(pco, type="none", choice = c(1,2), display="species") # , xlim = c(-1.5,3.2)
points(pco.species.scores[,2:3], pch = 19, col = "grey", cex = .5)
points(pco.sp.list[,2:3], pch = "+", col = "red", cex = .75)
#set.seed(314)
ordipointlabel(pco, display="species", select=pco.sp.list$label, pch = 19, col = "red", cex = .8, add=T)


# Temporal beta-diversity index  --------

coral.w.site.pre <- coral.w.site[coral.cov.site$pre_post == "Pre",] %>% data.frame()
coral.w.site.post <- coral.w.site[coral.cov.site$pre_post == "Post",] %>% data.frame()

coral.cov.site.pre <- coral.cov.site[coral.cov.site$pre_post == "Pre",]
coral.cov.site.post <- coral.cov.site[coral.cov.site$pre_post == "Post",]

coral.cover.pre <- rowSums(coral.w.site.pre)
coral.cover.post <- rowSums(coral.w.site.post)
coral.cover.delta <- coral.cover.post - coral.cover.pre
  
tbi <- TBI(coral.w.site.pre, coral.w.site.post)
tbi.pa <- TBI(coral.w.site.pre, coral.w.site.post, pa.tr = T)

tbi.dat <- data.frame(subset(coral.cov.site.post, select = c(site_code, longitude, latitude, cluster, maxDHW, maxSSTA, interval)),
                      TBI = tbi$TBI, 
                      TBI.p = tbi$p.TBI, 
                      TBI_PA = tbi.pa$TBI, 
                      TBI_PA.p = tbi.pa$p.TBI,
                      coral.cover.pre,
                      coral.cover.post,
                      coral.cover.delta)


tbi.dat <- tbi.dat %>% left_join(subset(sites.scores, select = c(site_code, CAP1, CAP2, PCO1, PCO2), pre_post == "Pre"))

# Exploratory analysis: relationships between TBI, delta CC, DHW, PCO1 and 2
pairs.panels(tbi.dat[,-(1:4)])

par(mfcol = c(4,1), mai = c(.5,.7,0,.1))
boxplot(maxDHW ~ cluster, data = tbi.dat, ylab = "DHW")
boxplot(TBI ~ cluster, data = tbi.dat, ylab = "TBI")
boxplot(TBI_PA ~ cluster, data = tbi.dat, ylab = "TBI_PA")
boxplot(coral.cover.delta ~ cluster, data = tbi.dat, ylab = "Change in live coral cover")
abline(h = 0, lty = 2, col = "red")

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

# Map TBI
ggplot() +
  geom_polygon(data=aus2, aes(long, lat, group=group), fill="lightgray", color="darkgray") +
  #coord_map(xlim=c(118,128), ylim=c(-20,-11)) +
  xlab(expression(paste(Longitude^o, ~'E'))) +
  ylab(expression(paste(Latitude^o, ~'S'))) +
  geom_point(data = tbi.dat, aes(x=longitude, y=latitude, col = TBI), size=2, shape=19) +
  scale_colour_distiller(type = "div", name = "TBI")+
  theme_light() +
  theme(legend.position = "right", legend.direction = "vertical")

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

ggplot(data = tbi.dat.dist, aes(x=PCO.dist.pre.post, y=TBI))+
  geom_point()+
  stat_smooth(method = "lm", se = FALSE) +
  stat_cor(method="pearson")

ggplot(data = tbi.dat.dist, aes(x=coral.cover.delta, y=TBI))+
  geom_point()+
  stat_smooth(method = "lm", se = FALSE) +
  stat_cor(method="pearson")


summary(lm(TBI ~ CAP.dist.pre.post, data = tbi.dat.dist))
summary(lm(TBI ~ PCO.dist.pre.post, data = tbi.dat.dist))
summary(lm(TBI ~ coral.cover.delta, data = tbi.dat.dist))

# Test which species changed over time in each cluster -------

for (i in c(1:7)) {
  
  t <- tpaired.krandtest(sqrt(coral.w.site.pre[tbi.dat$cluster == i,]), sqrt(coral.w.site.post[tbi.dat$cluster == i,]))
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
  simper.dat <- sqrt(rbind(coral.w.site.pre[tbi.dat$cluster == i,], coral.w.site.post[tbi.dat$cluster == i,]))
  simper <- simper(simper.dat, group = c(rep("Pre", dim(coral.w.site.pre[tbi.dat$cluster == i,])[1]), rep("Post", dim(coral.w.site.post[tbi.dat$cluster == i,])[1])))
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
  
  write.csv(simper_species_sub, paste("SIMPER/simper_cluster_", i, ".csv", sep = ""), row.names = F)
}

# Export pre-disturbance indicators
write.csv(ind.Coral.tb, "SIMPER/indicator_coral_categories.csv")

# Test if indicator taxa lost their indicator value after heatwave --------


ind.Coral.post <- indval(subset(coral.w.site.post, select = colSums(coral.w.site.post)>0), coral.cov.site.post$cluster, numitr=100)
ind.Coral.post.tb <- summary_indval(ind.Coral.post)
#write.csv(ind.Coral.post.tb, "MRT_Indicator_Corals_post.csv", row.names = T)

table(row.names(ind.Coral.tb) %in% row.names(ind.Coral.post.tb))
row.names(ind.Coral.tb)[!row.names(ind.Coral.tb) %in% row.names(ind.Coral.post.tb)]





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
                    gbm.x = c(4,5,6,15,16,19),
                    gbm.y = 8,
                    family = "gaussian",
                    tree.complexity = 2,
                    learning.rate = 0.001,
                    bag.fraction = 0.7)

# Explained deviance = 39.4%%
par(mai = c(.6,.6,.2,.2))
gbm.plot(gbm.TBI, write.title=F, n.plots = 5, plot.layout = c(6,1))
summary(gbm.TBI)


qqnorm(gbm.TBI$residuals, pch=19)
qqline(gbm.TBI$residuals)
# find.int <- gbm.interactions(gbm.TBI)
# find.int$rank.list
# gbm.perspec(gbm.TBI,4,3)

# GBM on TBI_PA
gbm.TBI.PA <- gbm.step(data = tbi.dat,
                       gbm.x = c(4,5,6,17,18,19),
                       gbm.y = 10,
                       family = "gaussian",
                       tree.complexity = 2,
                       learning.rate = 0.001,
                       bag.fraction = 0.7)

# Explained deviance = 30.0%%
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
                   gbm.x = c(4,5,6,17,18,19),
                   gbm.y = 14,
                   family = "gaussian",
                   tree.complexity = 2,
                   learning.rate = 0.001,
                   bag.fraction = 0.7)

# Explained deviance = 46.5%
par(mai = c(.6,.6,.2,.2))
gbm.plot(gbm.HC, write.title=F, n.plots = 5, plot.layout = c(6,1))
summary(gbm.HC)

qqnorm(gbm.HC$residuals, pch=19)
qqline(gbm.HC$residuals)

# find.int <- gbm.interactions(gbm.HC)
# find.int$rank.list
gbm.perspec(gbm.HC, 3, 2, x.range = c(0,25), y.range = c(0,7), z.range = c(-25,5))

gbm.HC.simpl <- gbm.simplify(gbm.HC, n.drops = 3)



# Paired boxplots of HCC pre. vs post for each cluster --------


HCC.pre.post <- data.frame(subset(coral.cov.site, select = c(site_code, cluster, pre_post)), LiveCoralCover = rowSums(coral.w.site))
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


HCC.pre.post.summary <- HCC.pre.post %>% mutate(pre_post = NULL) %>% pivot_wider(names_from = PrePost, values_from = LiveCoralCover) %>%
  mutate(Diff = Post - Pre) %>%
  group_by(cluster) %>% summarize(mean.Diff = mean(Diff, na.rm = T), sd.Diff = sd(Diff, na.rm = T))

HCC.pre.post.test <- HCC.pre.post %>% mutate(pre_post = NULL) %>% pivot_wider(names_from = PrePost, values_from = LiveCoralCover) %>%
  mutate(Diff = Post - Pre) %>%
  group_by(cluster) %>% summarize(P = t.test(Diff)$p.value)



