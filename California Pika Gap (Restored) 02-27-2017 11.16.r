### Ochotona princeps - Spatially-varying importance of variables
### Adam B. Smith | 2017-02

# source('C:/ecology/Drive/Research/Iconic Species/Analysis - California Gap/California Pika Gap.r')

rm(list=ls())
memory.limit(memory.limit() * 2^30)

library(checkpoint)
checkpoint(
	snapshotDate = '2016-12-01',
	project = 'C:/ecology/Drive/Research/Iconic Species',
	checkpointLocation = 'C:/ecology/Drive/Research/Iconic Species/California Pika Gap'
)

options(stringsAsFactors=FALSE)
gc()


### CONTENTS ###
### libraries, variables, and functions ###

######################
### generalization ###
######################

#########################
### collate presences ###
#########################


###########################################
### libraries, variables, and functions ###
###########################################

workDir <- 'C:/ecology/Drive/Research/Iconic Species/'

library(scales)
library(sp)
library(rgdal)
library(mboost)
library(geosphere)
library(cluster)
library(rJava)
library(fossil)
library(rgeos)

source('C:/ecology/Drive/R/!Omnibus.r')
source('C:/ecology/Drive/R/Graphics/Graphics - Spoke Plot.r')
source('C:/ecology/Drive/R/SDM/SDM - Train SDMs.r')
source('C:/ecology/Drive/R/SDM/SDM - Evaluate Models.r')
source('C:/ecology/Drive/R/SDM/SDM - Collate Data for a Species and Create Partitions.r')
source('C:/ecology/Drive/R/Center and Standardize Predictors.r')

say('#########################')
say('### collate presences ###')
say('#########################')

pres <- readRDS(paste0(workDir, 'Species Records - Pika/!Collated Data 2016-06-30 1256/02 Ochotona princeps - Usable - Presences & Non-detections 1895-2015 - PRISM & DayMet Climate Data Extracted.rds'))

## subset to presences AND years 1981-2015 AND providers who agreed to allow their data to be used and publically-available databases
pres <- pres[
	pres$origRecentPikaOrSignsObserved & pres$obsYear >= 1981 & pres$obsYear <= 2015 & (
		pres$xxxx == 'xxxx' |
		pres$xxxx == 'xxxx' |
		pres$xxxx == 'xxxx' |
		pres$xxxx == 'xxxx' |
		pres$xxxx == 'xxxx' |
		pres$xxxx == 'xxxx' |
		pres$xxxx == 'xxxx' |
		pres$xxxx == 'xxxx' |
	), ]

## subset to presences in Sierra Nevada EPALevel III ecoregion
epa3 <- readOGR('', '', verbose=FALSE)

pres$epa3 <- extract(epa3, cbind(pres$longWgs84, pres$lastWgs84))$xxxxxxxxx
pres <- pres[pres$epa3 == 'xxxx', ]

## save
saveRDS(pres, paste0(workDir, '/Analysis - California Gap/Presences for California Gap.rds'))







