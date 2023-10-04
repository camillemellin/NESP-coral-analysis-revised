
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
#load("NESP-coral-analysis-revised.RData")
#source("brt.functions.R")

# Load covariates
#dhw.raster <- raster::raster("Marine_heatwaves_Camille/Data/cropped_dhw_max_2016.grd")
#ssta.raster <- raster::raster("Marine_heatwaves_Camille/Data/cropped_ssta_max_2016.grd")

# Load benthic data
bent.dat <- read.csv("NESP-Extract-updated_2021-12-03.csv")

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
                                             percent_cover, is_finished, is_NESP_pre.post, is.infilled.in.previous.extract, 
                                             main.labeller = finished.data...main.benthic.labeller) %>%
                                filter(is_finished == "yes" & RLE_category %in% c("0", "Kelp"))

# Check number of labellers per survey_id
bent.dat.check2 <- bent.dat.check %>% group_by(survey_id, main.labeller) %>% summarize() %>%
                                      group_by(survey_id) %>% summarize(n_labellers = n_distinct(main.labeller))

# Wide format
bent.dat.check.wide <- bent.dat.check %>% 
                        group_by(survey_id, location, survey_year, is_NESP_pre.post, labeller = main.labeller, label) %>%
                        summarize(percent_cover = sum(percent_cover)) %>%
                        pivot_wider(names_from = label, values_from = percent_cover, values_fill = 0) %>% data.frame

# Calculate relative cover (divide by total cover of non-coral categories)
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
kw.pvalues <- data.frame(label = names(bent.dat.check.wide)[6:44], kw.pvalues = NA, chi2 = NA)

for (i in 6:44) {
  kw.test <- kruskal.test(bent.dat.check.wide[,i] ~ bent.dat.check.wide$labeller)
  kw.pvalues$kw.pvalues[i-5] <- ifelse(kw.test$p.value < 0.001, 0, 0.05)
  kw.pvalues$chi2[i-5] <- kw.test$statistic
}

par(mfcol = c(5,2), mai = c(.3,.6,.2,.2))
boxplot(Bare.Rock ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Bare.Rock), col = "red")

boxplot(Sand ~ labeller, data = bent.dat.check.wide) 
abline(h = mean(bent.dat.check.wide$Sand), col = "red")

boxplot(Coral.rubble.with.turf.encrusting.algae ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Coral.rubble.with.turf.encrusting.algae), col = "red")

boxplot(Turfing.algae...2.cm.high.algal.sediment.mat.on.rock. ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Turfing.algae...2.cm.high.algal.sediment.mat.on.rock.), col = "red")

boxplot(Filamentous.green.algae_epiphyte ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Filamentous.green.algae_epiphyte), col = "red")

boxplot(Sponges..massive. ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Sponges..massive.), col = "red")

boxplot(Sponges..encrusting. ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Sponges..encrusting.), col = "red")

boxplot(Encrusting.leathery.algae ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Encrusting.leathery.algae), col = "red")

boxplot(Geniculate.coralline.algae ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Geniculate.coralline.algae), col = "red")

boxplot(Slime..not.trapping.sediment. ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Slime..not.trapping.sediment.), col = "red")

boxplot(Filamentous.red.algae_epiphyte ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Filamentous.red.algae_epiphyte), col = "red")

boxplot(Small..2cm.foliose.algal.cover..not.trapping.sediment. ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Small..2cm.foliose.algal.cover..not.trapping.sediment.), col = "red")

boxplot(Medium.foliose.green.algae ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Medium.foliose.green.algae), col = "red")

boxplot(Green.calcified.algae..Halimeda. ~ labeller, data = bent.dat.check.wide)
abline(h = mean(bent.dat.check.wide$Green.calcified.algae..Halimeda.), col = "red")



