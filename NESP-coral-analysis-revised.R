
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
bent.dat <- read.csv("NESP-Extract_2021-11-09.csv")

bent.dat$survey_date <- as.Date(bent.dat$survey_date, tryFormats = "%d/%m/%y")
bent.dat$survey_year <- format(bent.dat$survey_date, '%Y')

bent.dat$survey_id <- factor(bent.dat$survey_id)
bent.dat$dataset_id <- factor(bent.dat$dataset_id)


# Remove duplicate surveys
dataset_id_duplicates <- c("912346798_RLS Australian Coral Species List_982", 
                           "912351889_RLS Catalogue_946", 
                           "912352031_RLS Catalogue_946",
                           "912352656_RLS Catalogue_939",
                           "912352687_RLS Catalogue_928")

bent.dat <- bent.dat %>% filter(!dataset_id %in% dataset_id_duplicates)

# Add 'unfinished' missing info
bent.dat$is_finished[!(bent.dat$is_finished %in% c('yes', 'no'))] <- 'no'


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



# Check labeller bias using updated extract (for Lizzie) ---------



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
#write.csv(ind.Coral.tb, "MRT_Indicator_Corals.csv", row.names = T)


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
                            mutate(site_code = NULL, pre_post = NULL) %>%
                            sqrt()
  
coral.w.site <- subset(coral.w.site, select=colSums(coral.w.site)>0)

# Remove sites surveyed only once
coral.w.site <- coral.w.site[!is.na(coral.cov.site$interval),]
coral.cov.site <- coral.cov.site[!is.na(coral.cov.site$interval),]

# CAPSCALE

# capscale.preds <- data.frame(clust.bentCat = bentCat_w_site$clust.bentCat, PrePost = bentCat_w_site$PrePost, clust.prepost,
#                              maxSSTA = ssta_site$maxSSTA, maxDHW = ssta_site$maxDHW,
#                              resurveyed = sampl_years$resurveyed, from2016 = sampl_years$from2016)


capscale <- capscale(coral.w.site ~ clust.prepost + Condition(interval), data = coral.cov.site, distance = "bray", add = T)
pco <- capscale(coral.w.site ~ 1, data = coral.cov.site, distance = "bray", add = T)

capscale
plot(capscale)
anova(capscale, permutations = how(nperm=99), by = "term")

fit.preds <- subset(coral.cov.site, select = c(maxDHW, maxSSTA))
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
