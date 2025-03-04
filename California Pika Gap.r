### OCHOTONA PRINCEPS - SPATIALLY-VARYING IMPORTANCE OF VARIABLES
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2017-02
###
### source('C:/Ecology/Research/Pikas - California Gap (Erik Beever et al)/pika_californiaGap/California Pika Gap.r')
### source('E:/Adam/Research/Pikas - California Gap (Erik Beever et al)/pika_californiaGap/California Pika Gap.r')
###
### CONTENTS ###
### libraries, variables, and functions ###
### define study region ###
### process predictor rasters ###
### collate TRAINING detections and non-detections ###
### collate TEST detections and non-detections ###
### construct KDE on TRAINING occurrences to identify gap ###

### map of gap sampling in focal region ###
### map of gap sampling across Sierras ###
### map of gap sampling across Sierras sans survey plots ###

### assign weights to TRAINING sites based on spatial autocorrelation ###
### calculate spatial autocorrelation within TEST and between TEST and TRAINING survey sites ###
### plot TRAINING sites scaled by SAC-based weight ###
### plot TEST sites scaled by SAC-based weight ###
### generate background sites and calculate PCA ###
### PCA biplot ###
### match predictors with TRAINING and TEST sites ###
### train ENMs ###

### assess predictor importance ###
### write ENM rasters ###
### extract predictions to TEST survey sites ###
### calculate performance statistics against TEST survey sites ###
### boxplot of distribution of predictions at TEST sites ###
### PCA of background and TEST site classes ###
### visual comparison of background, TRAINING, and TEST sites by predictor ###
### map of predictions ###
### map of incorrect predictions ###
### confusion matrix of predictions ###

### make KML of test sites with predictions ###
### calculate climate means for gap area ###

### graph of elevation of calibration and survey sites by type ###
### ancillary analyses of sites ###
### measure gap area ###
### compare inter-point distances between sets of sites ###
### pairwise distances between surveyed test sites with detections ###

### response curve plots ###

###########################################
### libraries, variables, and functions ###
###########################################

	rm(list=ls())
	gc()

	set.seed(1234567890)
	
	# working drive
	drive <- 'C:/Ecology/'
	# drive <- 'E:/Adam/'

	# PRISM
	# prismDrive <- 'F:'
	prismDrive <- 'I:'
	
	# TerraClimate
	tcDrive <- 'F:'

	setwd(paste0(drive, '/Research/Pikas - California Gap (Erik Beever et al)'))

	prismCrs <- '+proj=longlat +datum=NAD83 +no_defs'

	###############################
	### libraries and functions ###
	###############################

		library(cowplot)
		library(data.table)
		library(dismo)
		library(enmSdm)
		library(gbm)
		library(rgdal)
		library(mgcv)
		# library(fasterRaster)
		library(geosphere)
		library(cluster)
		library(ggplot2)
		# library(fossil)
		library(party)
		library(patchwork)
		library(rgeos)
		library(raster)
		# library(rJava)
		library(scales)
		library(sf)
		library(sp)
		library(spatialEco)
		library(terra)

		library(openxlsx)
		
		library(omnibus)
		library(enmSdm)
		library(legendary)
		library(statisfactory)
		
		# rasterOptions(format='GTiff', overwrite=TRUE, tmpdir=paste0('C:/!Scratch/_raster'))
		
		### load predictor rasters
		loadPredRasters <- function() {
			
			env <- stack(c(
				listFiles('./Predictors/Derived', pattern='.tif'),
				'./Predictors/elevPrism_studyRegionExtent.tif'
			))
			
			names(env)[nlayers(env)] <- 'elevation'
			env
			
		}

		# create stack of environmental predictors
		stackEnv <- function() {
		
			out <- raster::stack(c(listFiles('./Data/Climate - Derived', pattern='.tif'), './Data/NDVI 1990-2010/ndvi.tif'))
			out <- subset(out, predictors)
			out
		
		}

		# get vector of PCs to use for predictors
		whichPcs <- function(pca, thold=0.90) {
			
			# pca		object of class "princomp"
			# thold		cumulative variance explaned must be <= this for PC to be used
			
			varExplained <- pca$sdev^2 / sum(pca$sdev^2)
			cumVarExplained <- cumsum(varExplained)
			
			usePcs <- which(cumVarExplained <= thold)
			usePcs
		
		}

		# ensemble prediction across ENMs
		predictEnsemble <- function(data, brt=NULL, gam=NULL, glm=NULL) {
			
			# data		data frame
			# brt, gam, glm	model objects, declare to speed up
			
			preds <- matrix(NA, ncol=length(algos), nrow=nrow(data))
			
			if (is.null(brt)) { load('./ENMs/ENM BRT.rda'); brt <- model }
			if (is.null(gam)) { load('./ENMs/ENM GAM.rda'); gam <- model }
			if (is.null(glm)) { load('./ENMs/ENM GLM.rda'); glm <- model }
		
			for (i in seq_along(algos)) {
				algo <- algos[i]
				model <- get(algo)
				preds[ , i] <- predict(model, data, type='response', n.trees=model$gbm.call$n.trees)
			}
			
			preds <- rowMeans(preds)
			preds
		
		}

	#################
	### variables ###
	#################
	
		# GRASS directory
		# grassDir <- c('C:/OSGeo4W64/', 'grass-7.4.1', 'osgeo4W')
		grassDir <- 'C:/Program Files/GRASS GIS 8.3' # Windows
	
		# days of year
		data(doyLeap, doyNonLeap, package='omnibus')
		
		# long/lat
		ll <- c('longWgs84', 'latWgs84')
		
		# distance by which to buffer Sierra Nevada EPA Level III ecoregion to create study region (in km)
		studyRegionBuffer <- 70
		
		# predictors
		# predictors <- c('eastness', 'northness', 'acuteCold_C', 'chronicHeat_C', 'growSeasonWaterDef_mm', 'summerSrad', 'summerNightTemp_C', 'swe_mm', 'ndvi', 'vpdmax')
		predictors <- c('eastness', 'northness', 'acuteCold_C', 'chronicHeat_C', 'growSeasonWaterDef_mm', 'summerSrad', 'summerNightTemp_C', 'swe_mm', 'ndvi')
		
		# length of distance bins for calculating characteristic scale of clustering of survey sites (in meters)
		sacInterval <- 20000
		
		# gray color palette
		grays <- paste0('gray', 0:100)
	
		# ENM algorithms
		algos <- c('brt', 'glm', 'gam')
	
		years <- 2010:2019
		say('Using ', max(years) - min(years) + 1, '-yr window for climate layers:  ', min(years), ' through ', max(years), '.')
	
	###############
	### options ###
	###############
	
		options(stringsAsFactors=FALSE)
		wopt <- list(filetype='GTiff')
		# rasterOptions(format='GTiff', overwrite=TRUE)

	#########################
	### project structure ###
	#########################
	
		dirCreate('./Data')
		dirCreate('./Figures & Tables')
		dirCreate('./Study Region')
		dirCreate('./KDE')
		dirCreate('./ENMs')

# say('###########################')
# say('### define study region ###')
# say('###########################')

	# ## study region extent derived from modified EPA Level III "Sierra Nevada" ecoregion plus 70-km buffer
	# say('study region')
		
		# epa3 <- shapefile(paste0(drive, '/Research Done/Iconic Species/Extents_Masks_Maps/EcoRegions/!OmernikEcoregions/us_eco_l3SarrREV'))

		# snEa <- epa3[epa3$L3_KEY == '5  Sierra Nevada', ]
		# snBuff <- gBuffer(snEa, width=1000 * studyRegionBuffer)

		# snBuff <- sp::spTransform(snBuff, CRS(prismCrs))
		# sn <- sp::spTransform(snEa, CRS(prismCrs))

		# shapefile(snBuff, './Study Region/epa3_sierraNevadaSarrModified_200kmBuffer', overwrite=TRUE)
		# shapefile(sn, './Study Region/epa3_sierraNevadaSarrModified', overwrite=TRUE)
		
		# rm(epa3); gc()
		
	# ### PRISM elevation
	# say('PRISM elevation')
	
		# elevPrism <- raster(paste0(drive, '/Ecology/Climate/PRISM/PRISM_us_dem_800m.tif'))
		# elevPrism <- crop(elevPrism, snBuff)
		# names(elevPrism) <- 'elevation_prism_m'
		
		# dirCreate('./Data/Topography - PRISM')
		# writeRaster(elevPrism, './Data/Topography - PRISM/elevation_prism_m_studyRegion')
		
	# ### mask
	# say('mask')
	
		# maskPrism <- 1 + elevPrism * 0
		# maskPrism <- rasterize(snBuff, maskPrism)
		# names(maskPrism) <- 'mask'
		# writeRaster(maskPrism, paste0('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus', studyRegionBuffer, 'kmBuffer'), datatype='INT1U')

		# rm(maskPrism); gc()
		
	# ### SRTM elevation
	# say('SRTM elevation')
	
		# srtm1 <- rast('F:/ecology/Topography/SRTM - CGIAR/Tiles/cut_n30w120.tif')
		# srtm2 <- rast('F:/ecology/Topography/SRTM - CGIAR/Tiles/cut_n30w150.tif')
	
		# snBuffWgs84 <- sp::spTransform(snBuff, getCRS('wgs84', TRUE))
		# srtm1 <- crop(srtm1, snBuff)
		# srtm2 <- crop(srtm2, snBuff)
		
		# srtmSn <- mosaic(srtm1, srtm2)
		
		# names(srtmSn) <- 'elevation_srtm_m'
	
		# dirCreate('./Data/Topography - SRTM')
		# writeRaster(srtmSn, './Data/Topography - SRTM/elevation_srtm_m.tif', datatype='INT2S', overwrite=TRUE)
		
		# rm(srtmSn); gc()
		
	# ### SRTM hill shade
	# say('SRTM hill shade')
	
		# snBuffLarger <- gBuffer(snEa, width=6 * 1000 * studyRegionBuffer)
		# snBuffLargerWgs84 <- sp::spTransform(snBuffLarger, getCRS('wgs84', TRUE))

		# srtm1 <- rast('F:/ecology/Topography/SRTM - CGIAR/Tiles/cut_n30w120.tif')
		# srtm2 <- rast('F:/ecology/Topography/SRTM - CGIAR/Tiles/cut_n30w150.tif')
	
		# srtm1 <- crop(srtm1, snBuffLargerWgs84)
		# srtm2 <- crop(srtm2, snBuffLargerWgs84)

		# expandedSrtm <- mosaic(srtm1, srtm2)

		# slope <- terrain(expandedSrtm, 'slope', unit='rad')
		# aspect <- terrain(expandedSrtm, 'aspect', unit='rad')
		
		# names(slope) <- 'slope_srtm_rad'
		# names(aspect) <- 'aspect_srtm_rad'
		
		# hs <- shade(slope, aspect, direction=45)
		# hs <- stretch(hs, 0, 1)
		# hs <- round(hs, 5)
		# names(hs) <- 'hillshade_srtm'
		# writeRaster(hs, './Data/Topography - SRTM/hillshade_srtm.tif', datatype='FLT4S', overwrite=TRUE)
		
		# hsEa <- project(hs, getCRS('climateNA'))
		# names(hsEa) <- 'hillshade_srtm_ea'
		# writeRaster(hsEa, './Data/Topography - SRTM/hillshade_srtm_ea.tif', datatype='FLT4S', overwrite=TRUE)

	# ## political entities
	# say('political geography')
	
		# # California state
		# gadm <- raster::getData('GADM', level=1, country='USA', path='C:/ecology/!Scratch')
		# west1 <- gadm[gadm$NAME_1 %in% c('California', 'Nevada', 'Oregon'), ]
		# save(west1, file='./Study Region/GADM California, Nevada, Oregon States.rda')
		
		# # counties
		# gadm <- raster::getData('GADM', level=2, country='USA', path='C:/ecology/!Scratch')
		# west2 <- gadm[gadm$NAME_1 %in% c('California', 'Nevada', 'Oregon'), ]
		# plumas <- gadm[gadm$NAME_1 == 'California' & gadm$NAME_2 == 'Plumas', ]
		# sierra <- gadm[gadm$NAME_1 == 'California' & gadm$NAME_2 == 'Sierra', ]
		# washoe <- gadm[gadm$NAME_1 == 'Nevada' & gadm$NAME_2 == 'Washoe', ]
		
		# nearGap <- gTouches(plumas, west2, byid=TRUE)
		# gapCounties <- west2[c(nearGap), ]
		# gapCounties <- rbind(plumas, sierra, washoe, gapCounties)
		
		# save(west2, file='./Study Region/GADM California, Nevada, Oregon Counties.rda')
		# save(plumas, file='./Study Region/GADM Plumas County.rda')
		# save(gapCounties, file='./Study Region/GADM GAP Counties.rda')

		# # plumas <- sp::spTransform(plumas, getCRS('climateNA', TRUE))
		# # plumasBuffer <- gBuffer(plumas, width=35000)
		# # plumasBuffer <- sp::spTransform(plumasBuffer, CRS(projection(gadm)))
		
		# # save(plumasBuffer, file='./Study Region/GADM Plumas County + 35-km Buffer.rda')

say('#################################')
say('### process predictor rasters ###')
say('#################################')

	### calculate raster representing mean monthly values for each variable across a set time period
	### from these calculate specific predictor variable rasters

	# calculate time period for climate layers
	load('./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.rda')
	
	png(paste0('./Figures & Tables/Year of Observations.png'), height=900, width=900, res=200)
		hist(surveys$obsYear, breaks=(min(surveys$obsYear) - 1):(max(surveys$obsYear) + 1), main='Observations', xlab='Year')
	dev.off()

	mask <- terra::rast(paste0('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus', studyRegionBuffer, 'kmBuffer.tif'))

	data('doyNonLeap', package='omnibus')

	# calculate monthly means of max/min temperature and ppt from PRISM and TerraClimate
	for (variable in c('tmin', 'tmax', 'srad', 'def', 'swe', 'vpdmax', )) { # ALL
	# for (variable in c('tmin', 'tmax', 'vpdmax')) { # PRISM variables
	# for (variable in c('srad', 'def', 'swe')) { # TerraClimate variables
		
		dirCreate('./Data/Climate - Monthly Means/', variable)

		months <- if (variable %in% c('tmax', 'srad', 'def', 'vpdmax')) {
			6:9
		} else if (variable %in% c( 'tmin', 'swe')) {
			1:12
		}

		for (month in months) {
		
			say(variable, ' | month ', month, ' | year', post=0)
		
			for (year in years) {
		
				say(year, post=ifelse(year == max(years), 1, 0))

				# mid-day of this month
				midMonthDoy <- doyNonLeap[15, paste0('month', month)]
			
				if (variable %in% c('ppt', 'tmin', 'tmax', 'vpdmax')) {
					thisMonthYear <- terra::rast(paste0(prismDrive, '/ecology/Climate/PRISM/working/an81/' , variable, '/monthly/', year, '/prism_', variable, '_us_30s_' , year, prefix(month, 2), '.tif'))
				} else if (variable %in% c('pet', 'srad', 'swe')) {
					thisMonthYear <- terra::rast(paste0(tcDrive, '/ecology/Climate/TerraClimate/ORIGINALS/', variable, '/TerraClimate_', variable, '_', year, '.nc'))
					thisMonthYear <- thisMonthYear[[month]]
				}
				
				thisMonthYear <- terra::resample(thisMonthYear, mask)
				thisMonthYear <- terra::crop(thisMonthYear, mask)
				
				thisMonth <- if (!exists('thisMonth', inherits=FALSE)) {
					thisMonthYear
				} else {
					c(thisMonth, thisMonthYear)
				}
				
			} # next year
			
			# calculate mean for this month across years
			meanForMonth <- mean(thisMonth)
			names(meanForMonth) <- variable
			
			terra::writeRaster(meanForMonth, paste0('./Data/Climate - Monthly Means/', variable, '/', variable, '_month', prefix(month, 2), '_mean', min(years), 'to', max(years), '.tif'), overwrite=TRUE, wopt=c(wopt, names=variable))
			
			rm(thisMonth, meanForMonth); gc()
			
		} # next month
		
	} # next variable

	dirCreate('./Data/Climate - Derived')
	
	say('ASPECT')
	
		elev <- terra::rast(paste0(drive, '/Ecology/Climate/PRISM/PRISM_us_dem_800m.tif'))
		aspect <- terra::terrain(elev, v='aspect')
		aspect <- terra::crop(aspect, mask)
		aspect <- pi * aspect / 180
		northness <- sin(aspect)
		eastness <- cos(aspect)
		names(northness) <- 'northness'
		names(eastness) <- 'eastness'
		
		terra::writeRaster(northness, './Data/Climate - Derived/northness.tif', overwrite=TRUE, wopt=c(wopt, names='northness'))
		terra::writeRaster(eastness, './Data/Climate - Derived/eastness.tif', overwrite=TRUE, wopt=c(wopt, names='eastness'))
		
	say('CHRONIC HEAT')
	
		rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/tmax', pattern='.tif'))
		chronicHeat <- mean(rasts)
		terra::writeRaster(chronicHeat, './Data/Climate - Derived/chronicHeat_C.tif', overwrite=TRUE, wopt=c(wopt, names='chronicHeat_C'))
		
	say('SUMMER NIGHTTIME HEAT')
	
		rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/tmin', pattern='.tif'))
		rasts <- rasts[[6:9]]
		summerNightHeat <- mean(rasts)
		terra::writeRaster(summerNightHeat, './Data/Climate - Derived/summerNightTemp_C.tif', overwrite=TRUE, wopt=c(wopt, names='summerNightTemp_C'))
		
	say('ACUTE COLD')
	
		rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/tmin', pattern='.tif'))
		acuteCold <- min(rasts)
		terra::writeRaster(acuteCold, './Data/Climate - Derived/acuteCold_C.tif', overwrite=TRUE, wopt=c(wopt, names='acuteCold_C'))
		
	say('VPDMAX')
	
		rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/vpdmax', pattern='.tif'))
		vpdmax_haPa <- mean(rasts)
		terra::writeRaster(vpdmax_haPa, './Data/Climate - Derived/vpdmax_haPa.tif', overwrite=TRUE, wopt=c(wopt, names='vpdmax_haPa'))
		
	say('SWE')

		rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/swe', pattern='.tif'))
		swe <- sum(rasts)
		terra::writeRaster(swe, './Data/Climate - Derived/swe_mm.tif', overwrite=TRUE, wopt=c(wopt, names='swe_mm'))
		
	say('SUMMER SRAD')

		rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/srad', pattern='.tif'))
		srad <- sum(rasts)
		terra::writeRaster(srad, './Data/Climate - Derived/summerSrad.tif', overwrite=TRUE, wopt=c(wopt, names='summerSrad'))

	say('GSWD')

		rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/def', pattern='.tif'))
		gswb <- sum(rasts)
		terra::writeRaster(gswb, './Data/Climate - Derived/growSeasonWaterDef_mm.tif', overwrite=TRUE, wopt=c(wopt, names='growSeasonWaterDef_mm'))
	
	say('NDVI')
		
		ndvi <- terra::rast('./Data/NDVI 1990-2010/Raw/NDVI.tif')
		ndvi <- mean(ndvi)
		ndvi <- terra::project(ndvi, mask)
		ndvi <- crop(ndvi, mask)
		terra::writeRaster(ndvi, paste0('./Data/NDVI 1990-2010/ndvi.tif'), overwrite=TRUE, wopt=c(wopt, names='ndvi'))

	# ##########################
	# say('PLOT ALL PREDICTORS')
	
	# pres <- read.csv('./Data/Occurrences/Test Surveys 02 Cleaned.csv')

	# load('./Study Region/GADM California, Nevada, Oregon States.rda')
	# load('./Study Region/GADM California, Nevada, Oregon Counties.rda')
	# load('./Study Region/GADM Plumas County.rda')
	
	# west1 <- vect(west1)
	# west2 <- vect(west2)
	# plumas <- vect(plumas)
	
	# studyRegion <- vect('./Study Region/epa3_sierraNevadaSarrModified_200kmBuffer.shp')

	# period <- paste0('(', min(years), '-', max(years), ')')
	
	# # get rasters
	# rasts <- rast(listFiles('./Data/Climate - Derived', pattern='.tif'))
	# ndvi <- rast('./Data/NDVI 1990-2010/ndvi.tif')
	# elevPrism <- rast('./Data/Topography - PRISM/elevation_prism_m_studyRegion.tif')
	# rasts <- c(rasts, ndvi, elevPrism)
	
	# grays <- paste0('gray', 0:100)

	# say('study region')
	# png(paste0('./Figures & Tables/Predictor Maps - Study Region.png'), width=1920, height=1080, res=200)
	
		# par(mfrow=c(2, 5), oma=c(0, 0, 0, 2), mar=c(0, 0, 0, 0), mgp=c(1, 0.2, 0), cex.main=0.6, cex.lab=0.5, cex.axis=0.5, tck=-0.02)

		# overlays <- function() {
		
			# plot(west2, border='gray50', lwd=0.3, add=TRUE)
			# plot(west1, border='gray20', lwd=0.7, add=TRUE)
			# plot(plumas, border='chartreuse', lwd=0.7, add=TRUE)
			# plot(studyRegion, border='chartreuse', lwd=1, add=TRUE)
			# points(pres$longWgs84, pres$latWgs84, pch=3, cex=0.4, col='yellow')
			
		# }
		
		# plot(studyRegion, main=paste0('Elevation\n', period)); plot(rasts[['elevation_prism_m_studyRegion']], col=grays, add=TRUE); overlays()
		# plot(studyRegion, main=paste0('Northness\n', period)); plot(rasts[['northness']], col=grays, add=TRUE); overlays()
		# plot(studyRegion, main=paste0('GS Chronic Heat\n', period)); plot(rasts[['chronicHeat_C']], col=grays, add=TRUE); overlays()
		# plot(studyRegion, main=paste0('Maximum Vapor Pressure Deficit\n', period)); plot(rasts[['vpdmax_haPa']], col=grays, add=TRUE); overlays()
		# plot(studyRegion, main=paste0('Acute Cold\n', period)); plot(rasts[['acuteCold_C']], col=grays, add=TRUE); overlays()
		# plot(studyRegion, main=paste0('GS Night Heat\n', period)); plot(rasts[['summerNightTemp_C']], col=grays, add=TRUE); overlays()
		# plot(studyRegion, main=paste0('SWE\n', period)); plot(rasts[['swe_mm']], col=grays, add=TRUE); overlays()
		# plot(studyRegion, main=paste0('GS Atmospheric Water Deficit\n', period)); plot(rasts[['growSeasonWaterDef_mm']], col=grays, add=TRUE); overlays()
		# plot(studyRegion, main=paste0('GS Solar Radiation\n', period)); plot(rasts[['summerSrad']], col=grays, add=TRUE); overlays()
		# plot(studyRegion, main='NDVI\n(1990-2010)'); plot(rasts[['ndvi']], col=grays, add=TRUE); overlays()
		
	# dev.off()
	
	# say('gap')
	# plumasBuff <- buffer(plumas, 30000)
	# rastsCrop <- crop(rasts, plumasBuff)
	
	# png(paste0('./Figures & Tables/Predictor Maps - Gap.png'), width=1920, height=700, res=200)
	
		# par(mfrow=c(2, 5), oma=c(0, 0, 0, 2), mar=c(0, 0, 0, 0), mgp=c(1, 0.2, 0), cex.main=0.6, cex.lab=0.5, cex.axis=0.5, tck=-0.02)

		# overlays <- function() {
		
			# plot(west2, border='gray50', lwd=0.3, add=TRUE)
			# plot(west1, border='gray20', lwd=0.7, add=TRUE)
			# plot(plumas, border='chartreuse', lwd=1, add=TRUE)
			# plot(studyRegion, border='chartreuse', lwd=1, add=TRUE)
			# points(pres$longWgs84, pres$latWgs84, pch=3, cex=0.8, col='red')
			
		# }
		
		# plot(plumasBuff, border=NA, main=paste0('Elevation\n', period)); plot(rastsCrop[['elevation_prism_m_studyRegion']], col=grays, add=TRUE); overlays()
		# plot(plumasBuff, border=NA, main=paste0('Northness\n', period)); plot(rastsCrop[['northness']], col=grays, add=TRUE); overlays()
		# plot(plumasBuff, border=NA, main=paste0('GS Chronic Heat\n', period)); plot(rastsCrop[['chronicHeat_C']], col=grays, add=TRUE); overlays()
		# plot(plumasBuff, border=NA, main=paste0('Maximum Vapor Pressure Deficit\n', period)); plot(rastsCrop[['vpdmax_haPa']], col=grays, add=TRUE); overlays()
		# plot(plumasBuff, border=NA, main=paste0('Acute Cold\n', period)); plot(rastsCrop[['acuteCold_C']], col=grays, add=TRUE); overlays()
		# plot(plumasBuff, border=NA, main=paste0('GS Night Heat\n', period)); plot(rastsCrop[['summerNightTemp_C']], col=grays, add=TRUE); overlays()
		# plot(plumasBuff, border=NA, main=paste0('SWE\n', period)); plot(rastsCrop[['swe_mm']], col=grays, add=TRUE); overlays()
		# plot(plumasBuff, border=NA, main=paste0('GS Atmospheric Water Deficit\n', period)); plot(rastsCrop[['growSeasonWaterDef_mm']], col=grays, add=TRUE); overlays()
		# plot(plumasBuff, border=NA, main=paste0('GS Solar Radiation\n', period)); plot(rastsCrop[['summerSrad']], col=grays, add=TRUE); overlays()
		# plot(plumasBuff, border=NA, main='NDVI\n(1990-2010)'); plot(rastsCrop[['ndvi']], col=grays, add=TRUE); overlays()
		
	# dev.off()

# say('######################################################')
# say('### collate TRAINING detections and non-detections ###')
# say('######################################################')	
	
	# dirCreate('./Data/Occurrences')
	# surveys <- readRDS(paste0(drive, '/Research Done/Iconic Species/Species Records - Pika/!Collated Data 2016-06-30 1256/00 Pika - Cleaned Using R - Usable Records for All Ochotona.rds'))

	# mask <- raster(paste0('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus', studyRegionBuffer, 'kmBuffer.tif'))
	
	# ## subset to presences AND years 1981-2015 AND providers who agreed to allow their data to be used and publicly-available databases
	# surveys <- surveys[
		# !is.na(surveys$origRecentPikaOrSignsObserved) &
		# surveys$origRecentPikaOrSignsObserved & surveys$obsYear >= 1981 & surveys$obsYear <= 2015 & (
			# surveys$contact == 'NA' |
			# surveys$contact == 'Adam Smith;' |
			# surveys$contact == 'Chris Ray (Colorado University, Boulder)' |
			# surveys$contact == 'Erik Beever, Chris Ray' |
			# surveys$contact == 'Erik Beever' |
			# surveys$contact == 'Erik Beever;' |
			# surveys$contact == 'Ken Goehring;' |
			# surveys$contact == 'Joseph Stewart' |
			# surveys$contact == 'Joseph Stewart, David Wright, Chris Curlis' |
			# surveys$contact == 'Jason Brewer & Mary Flores' |
			# surveys$contact == 'Clint Epps, Jessica Castillo' |
			# surveys$contact == 'Jess Castillo, Clint Epps' |
			# surveys$contact == 'Connie Millar'
		# ), ]

	# ## subset to presences in study region
	# studyRegion <- raster(paste0('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus', studyRegionBuffer, 'kmBuffer.tif'))
	# inStudyRegion <- raster::extract(mask, surveys[ , ll])

	# surveys$origUsable <- as.logical(surveys$origUsable)
	# surveys <- surveys[surveys$origUsable, ]
	
	# surveys <- surveys[!is.na(inStudyRegion), ] # just in Sierra Nevada + buffer study region
	# surveys <- surveys[!(surveys$origRecencyOfSighting %in% 'years to decades (old pellets only)'), ]

	# save(surveys, file='./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.rda')

# say('##################################################')
# say('### collate TEST detections and non-detections ###')
# say('##################################################')
	
	# ### clean
	# ok <- file.copy('./Data/Occurrences/Gap Surveys - Raw/TR_Edits_NC hole pika   locations__updated16Sept2016_KBO__Status+EvidenceTypeAdded__FINAL-SEND....xlsx', './Data/Occurrences/Test Surveys 00.xlsx', overwrite=TRUE)
	
	# say('Raw test occurrence data file copied: ', ok)
	
	# testSurveys <- readxl::read_excel('./Data/Occurrences/Test Surveys 00.xlsx', range='All patches, sorted Geographica!A3:L264')
	
		# names(testSurveys)[names(testSurveys)=='Map Label in CalTopo.com'] <- 'mapLabel'
		# names(testSurveys)[names(testSurveys)=='...2'] <- 'origOccStatus'
		# names(testSurveys)[names(testSurveys)=='Latitude (dec deg)'] <- 'latWgs84'
		# names(testSurveys)[names(testSurveys)=='Longitude (dec deg)'] <- 'longWgs84'
		# names(testSurveys)[names(testSurveys)=='Zone (in relation to predicted pika occupancy map): red is highest likelihood, then orange, yellow, grey, and white'] <- 'predZone'
		# names(testSurveys)[names(testSurveys)=='Size       (in acres)'] <- 'patchSize_ac'
		# names(testSurveys)[names(testSurveys)=='Elevation (in feet)'] <- 'elev_ft'
		# names(testSurveys)[names(testSurveys)=='Talus Quality/Likelihood of Pika'] <- 'talusQual'
		# names(testSurveys)[names(testSurveys)=='...9'] <- 'empty1'
		# names(testSurveys)[names(testSurveys)=='...10'] <- 'latitude_dms'
		# names(testSurveys)[names(testSurveys)=='...11'] <- 'longitude_dms'
		# names(testSurveys)[names(testSurveys)=='Type of evidence (see Table 3 for classes)'] <- 'typeOfEvidence'
	
	# ### site-level corrections
	# this <- which(testSurveys$mapLabel == 'Grey 13+' & is.na(testSurveys$longWgs84) & is.na(testSurveys$latWgs84))
	# testSurveys$longWgs84[this] <- -121.3745
	# testSurveys$latWgs84[this] <- 40.2122
	# testSurveys$elev_ft[this] <- 2241 * 3.28084 # meters to feet
	
	# # get sites that 1) were surveyed; 2) were potentially suitable (talus-wise) for pikas
	# testSurveys <- testSurveys[!is.na(testSurveys$origOccStatus), ]
	# testSurveys <- testSurveys[testSurveys$origOccStatus %in% c('Evidence of only past pika occupancy', 'Sampled, but no pika evidence det\'d', 'Current pika occupancy (unequivocal)'), ]

	# testSurveys$status <- ifelse(
		# testSurveys$origOccStatus == 'Current pika occupancy (unequivocal)', '2 detected',
		# ifelse(testSurveys$origOccStatus == 'Evidence of only past pika occupancy', '1 recent absence',
		# ifelse(testSurveys$origOccStatus == 'Sampled, but no pika evidence det\'d', '0 long absence', NA)
	# ))
	
	# if (any(is.na(testSurveys$status))) stop('Site failed to be assigned an occupancy status.')

	# testSurveys$longWgs84 <- as.numeric(testSurveys$longWgs84)
	# testSurveys$latWgs84 <- as.numeric(testSurveys$latWgs84)
	
	# testSurveys <- testSurveys[complete.cases(testSurveys[ , ll]), ]
	
	# ### plot
	# load('./Study Region/GADM California, Nevada, Oregon Counties.rda')
	# load('./Study Region/GADM California, Nevada, Oregon States.rda')
	
	# testSurveysSp <- SpatialPointsDataFrame(testSurveys[ , ll], data=testSurveys, proj4=getCRS('wgs84', TRUE))
	
	# plot(testSurveysSp, col='white')
	# plot(west2, border='gray', add=TRUE)
	# plot(west1, add=TRUE)
	# points(testSurveysSp, pch=1)
	
	# write.csv(testSurveys, './Data/Occurrences/Test Surveys 02 Cleaned.csv', row.names=FALSE)
	
# say('#############################################################')
# say('### construct KDE on TRAINING occurrences to identify gap ###')
# say('#############################################################')

	# # mask
	# mask <- rast(paste0('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus', studyRegionBuffer, 'kmBuffer.tif'))
	# mask <- project(mask, getCRS('albersNA'))

	# # occurrences
	# load('./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.rda')
	# pres <- surveys[surveys$origRecentPikaOrSignsObserved, ]
	# presSp <- vect(pres, geom = ll, getCRS('wgs84'))
	# presSpEa <- project(presSp, mask)

	# # cropping occurrences to vicinity around gap because the long north-south extent of the Sierras makes the bandwidth biased latitudinally
	# load('./Study Region/GADM GAP Counties.rda')
	# gapCounties <- vect(gapCounties)
	# gapCounties <- buffer(gapCounties, 50000)
	# gapCounties <- project(gapCounties, getCRS('albersNA'))
	# extent <- as.vector(ext(gapCounties))
	
	# coords <- crds(presSpEa)
	# presSpEa <- presSpEa[coords[ , 'y'] <= extent[['ymax']] & coords[ , 'y'] >= extent[['ymin']], ]
	# coords <- crds(presSpEa)

	# ### KDE
	# # # Gaussian bandwidth
	# # bw <- 0.5 * sum(apply(crds(presSpEa), 2, sd)) * nrow(presSpEa)^(-1/6)
	
	# # Epanechnikov bandwidth
	# # bw <- 1.77 * 0.5 * sum(apply(crds(presSpEa), 2, sd, na.rm=T)) * nrow(presSp)^(-1/6)
	# # bw <- 1.77 * 0.5 * sum(apply(coords, 2, sd, na.rm=T)) * nrow(coords)^(-1/6)
	
	# kde <- ks::kde(coords)
	
	# longLat <- longLatRasts(mask)
	# names(longLat) <- c('x', 'y')
	# longLat <- as.data.frame(longLat, cells = TRUE)
	# kde <- predict(kde, x = longLat[ , c('x', 'y')])

	# kde <- enmSdmX::setValueByCell(mask, val = kde, cell = longLat$cell)
	
	# writeRaster(kde, './KDE/kde.tif', overwrite = TRUE)
	
# say('###########################################')
# say('### map of gap sampling in focal region ###')
# say('###########################################')

	# # generalize
	# focusBuff <- 25 # buffer size around test sites for generating focus of plot (in km)

	# # spatial data
	# load('./Study Region/GADM California, Nevada, Oregon States.rda')
	# load('./Study Region/GADM California, Nevada, Oregon Counties.rda')
	# load('./Study Region/GADM Plumas County.rda')
	# load('./Study Region/GADM GAP Counties.rda')
	# hs <- raster('./Data/Topography - SRTM/hillshade_srtm_ea.tif')
	
	# gapCounties <- sp::spTransform(gapCounties, getCRS('climateNA', TRUE))
	# plumas <- sp::spTransform(plumas, getCRS('climateNA', TRUE))
	# west1 <- sp::spTransform(west1, getCRS('climateNA', TRUE))
	# west2 <- sp::spTransform(west2, getCRS('climateNA', TRUE))

	# nam1 <- vect(paste0(drive, 'Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector Albers.gpkg'))
	# nam2 <- vect(paste0(drive, 'Research Data/GADM/Version 4.1/High Res North America Level 2 sans Great Lakes SpatVector Albers.gpkg'))
	# nam1 <- nam1[nam1$NAME_1 %in% c('California', 'Arizona', 'New Mexico', 'Utah', 'Oregon', 'Washington', 'Baja California', 'Baja California Sur', 'Idaho', 'Nevada', 'British Columbia', 'Sonora', 'Chihuahua', 'Montana', 'Sinola'), ]
	# nam2 <- nam2[nam2$NAME_1 %in% c('California'), ]

	# ### lakes
	# lakes <- shapefile(paste0(drive, '/Research Data/North America Rivers and Lakes USGS/hydrography_p_lakes_v2.shp'))
	# for (i in 1:20) {
		# maxs <- which.max(lakes$Shape_Leng)
		# lakes <- lakes[-maxs, ]
	# }
	# lakes <- sp::spTransform(lakes, getCRS('climateNA', TRUE))

	# # training occurrences
	# load('./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.rda')
	# trainPres <- surveys[surveys$origRecentPikaOrSignsObserved, ]
	# trainPres$origUsable <- as.logical(trainPres$origUsable)
	# trainPres <- trainPres[trainPres$origUsable, ]
	
	# trainPres <- trainPres[-which(trainPres$dataName == 'Connie Millar - Subset A - Pika Present' & trainPres$origLocality == 'Sierra Nevada; Humphreys Basin; Goethe rock glacier; at rg outlet stream; Fresno; CA' & trainPres$origElev == 3523 & trainPres$origNotes != 'RGC; granitic;' & trainPres$origDateOfObs == '8/16/2011'), ]

	# trainPresSp <- SpatialPointsDataFrame(trainPres[ , ll], data=trainPres, proj4=getCRS('wgs84', TRUE))
	# trainPresSpEa <- sp::spTransform(trainPresSp, getCRS('climateNA', TRUE))

	# # test occurrences
	# testSurveys <- read.csv('./Data/Occurrences/Test Surveys 02 Cleaned.csv', stringsAsFactors=FALSE, fileEncoding='latin1')
	# testSurveysSp <- SpatialPointsDataFrame(testSurveys[ , ll], data=testSurveys, proj4=getCRS('nad83', TRUE))
	# testSurveysSpEa <- sp::spTransform(testSurveysSp, getCRS('climateNA', TRUE))
	
	# # plot focus
	# focus <- gBuffer(testSurveysSpEa, width=1000 * focusBuff)
	
	# ### crop hillshade
	# plot(focus)
	# usr <- par('usr')
	# dev.off()

	# ext <- extent(usr)
	# ext <- as(ext, 'SpatialPolygons')
	# projection(ext) <- getCRS('climateNA')

	# hs <- crop(hs, ext)
	# hs <- hs - minValue(hs)
	# hs <- hs / maxValue(hs)
	
	# # kde
	# kde <- raster('./KDE/kde.tif')
	# beginCluster(4)
		# kde <- projectRaster(kde, crs=getCRS('climateNA'))
	# endCluster()

	# kdeVals <- extract(kde, trainPresSpEa)
	# quants <- quantile(kdeVals[kdeVals > 0], c(0.01, 0.025), na.rm = TRUE)

	# # kdeVals <- extract(kde, trainPresSpEa)
	# # quants <- quantile(kdeVals, c(0.05, 0.1))
	# breaks <- c(0, min(kdeVals[kdeVals > 0]), quants, 1)
	# kdeClass <- cut(kde, breaks=breaks)

	# # colors
	# kdeCols <- c(NA, 'deepskyblue', 'dodgerblue1', 'dodgerblue4')
	# for (i in seq_along(kdeCols)) kdeCols[i] <- alpha(kdeCols[i], 0.5)
	
	# ### shapes and colors for points
	# trainPresCol <- 'black'
	# testPresCol <- 'darkgreen'
	# # testShortTermAbsCol <- 'darkorange4'
	# testShortTermAbsCol <- 'yellow'
	# testLongTermAbsCol <- 'darkred'
	
	# trainPresPch <- 3
	# testPresPch <- 1
	# testShortTermPch <- 5
	# testLongTermPch <- 6

	# # plot
	# png('./Figures & Tables/Gap Sampling.png', width=1200, height=1200, res=300)
		
		# par(mar=c(2, 1, 1, 1), cex.axis=0.4, mgp=c(3, 0, 0.2))
		
		# plot(focus, border='white')
		
		# usr <- par('usr')
		# xs <- pretty(c(usr[1], usr[2]))
		# ys <- pretty(c(usr[3], usr[4]))
		# axis(1, at=xs, tck=0.01, labels=xs, col=NA, col.ticks='black')
		# axis(2, at=ys, tck=0.01, labels=ys, col=NA, col.ticks='black')
		
		# plot(hs, col=grays, legend=FALSE, add=TRUE)
		# plot(lakes, col = 'gray60', border = NA, add = TRUE)
		# plot(kdeClass, col=kdeCols, legend=FALSE, add=TRUE)
		# plot(west2, border='gray40', add=TRUE)
		# plot(west1, border='gray40', lwd=2, add=TRUE)

		# testSurveysSpEa <- testSurveysSpEa[order(testSurveysSpEa$status), ]
		
		# cols <- rep(NA, nrow(testSurveys))
		# cols[testSurveysSpEa$status == '0 long absence'] <- testLongTermAbsCol
		# cols[testSurveysSpEa$status == '1 recent absence'] <- testShortTermAbsCol
		# cols[testSurveysSpEa$status == '2 detected'] <- testPresCol
		
		# pchs <- rep(NA, nrow(testSurveys))
		# pchs[testSurveysSpEa$status == '0 long absence'] <- testLongTermPch
		# pchs[testSurveysSpEa$status == '1 recent absence'] <- testShortTermPch
		# pchs[testSurveysSpEa$status == '2 detected'] <- testPresPch
		
		# points(testSurveysSpEa, pch=pchs, col=cols, cex=0.6)
		# points(trainPresSpEa, pch=trainPresPch, cex=0.5, col=trainPresCol)
		# box()
		
		# # legend
		# legendBreaks(
			# 'bottomleft',
			# inset=0.01,
			# height=0.44,
			# width=0.27,
			# title='Previously-known\noccurrences',
			# titleAdj=c(0.5, 0.92),
			# col=kdeCols,
			# adjX=c(0.05, 0.225),
			# adjY=c(0.38, 0.85),
			# labels=c('\U2265min presence', paste0('\U2265', '1st percentile'), paste0('  \U2265', '2.5th percentile')),
			# labAdjX=0.58,
			# cex=0.54,
			# boxBg=alpha('white', 0.8)
		# )
		
		# usr <- par('usr')
		# width <- usr[2] - usr[1]
		# height <- usr[4] - usr[3]
		# x <- usr[1] + 0.01 * width
		# y <- usr[3] + 0.17 * height
		
		# legend(
			# x,
			# y,
			# legend=c('Previously-known occ.', 'No evidence', 'Formerly occupied', 'Currently occupied'),
			# pch=c(trainPresPch, testLongTermPch, testShortTermPch, testPresPch),
			# col=c(trainPresCol, testLongTermAbsCol, testShortTermAbsCol, testPresCol),
			# bty='n',
			# title='Surveys',
			# cex=0.45,
			# pt.cex=0.7
		# )
		
		# # scale bar
		# size <- 50000 # length of scale bar in meters
		# x <- usr[2] - size - 0.02 * width
		# x <- c(x, x + size)
		# y <- usr[3] + rep(0.02 * height, 2)
		
		# lines(x, y, lwd=4, xpd=NA, col='black', lend=1)
		
		# x <- mean(x)
		# y <- usr[3] + 0.05 * height
		# text(x, y[1], labels=paste(size / 1000, 'km'), cex=0.5)

		# title(sub=date(), cex.sub=0.3, outer=TRUE, line=-1)
		
	# dev.off()

						
	# # # # plot
	# # # png('./Figures & Tables/Gap Sampling Plumas Emphasized with No Survey Sites.png', width=1200, height=1200, res=300)
		
		# # # par(mar=c(2, 1, 1, 1), cex.axis=0.4, mgp=c(3, 0, 0.2))
		
		# # # plot(focus, border='white')
		
		# # # usr <- par('usr')
		# # # xs <- pretty(c(usr[1], usr[2]))
		# # # ys <- pretty(c(usr[3], usr[4]))
		# # # axis(1, at=xs, tck=0.01, labels=xs, col=NA, col.ticks='black')
		# # # axis(2, at=ys, tck=0.01, labels=ys, col=NA, col.ticks='black')
		
		# # # plot(hs, col=grays, legend=FALSE, add=TRUE)
		# # # plot(lakes, col = 'gray60', border = NA, add = TRUE)
		# # # plot(kdeClass, col=kdeCols, legend=FALSE, add=TRUE)
		# # # plot(west2, border='gray40', add=TRUE)
		# # # plot(west1, border='gray40', lwd=2, add=TRUE)
		# # # plot(west2[west2$NAME_2 == 'Plumas', ], border='gray40', lwd=3, add=TRUE)

		# # # testSurveysSpEa <- testSurveysSpEa[order(testSurveysSpEa$status), ]
		
		# # # # cols <- rep(NA, nrow(testSurveys))
		# # # # cols[testSurveysSpEa$status == '0 long absence'] <- testLongTermAbsCol
		# # # # cols[testSurveysSpEa$status == '1 recent absence'] <- testShortTermAbsCol
		# # # # cols[testSurveysSpEa$status == '2 detected'] <- testPresCol
		
		# # # # pchs <- rep(NA, nrow(testSurveys))
		# # # # pchs[testSurveysSpEa$status == '0 long absence'] <- testLongTermPch
		# # # # pchs[testSurveysSpEa$status == '1 recent absence'] <- testShortTermPch
		# # # # pchs[testSurveysSpEa$status == '2 detected'] <- testPresPch
		
		# # # points(trainPresSpEa, pch=trainPresPch, cex=0.5, col=trainPresCol)
		# # # # points(testSurveysSpEa, pch=pchs, col=cols, cex=0.6)
		# # # box()
		
		# # # # legend
		# # # legendBreaks(
			# # # 'bottomleft',
			# # # inset=0.01,
			# # # height=0.44,
			# # # width=0.27,
			# # # title='Training occurrence\ndensity',
			# # # titleAdj=c(0.5, 0.92),
			# # # col=kdeCols,
			# # # adjX=c(0.05, 0.225),
			# # # adjY=c(0.38, 0.85),
			# # # labels=c('\U2265min presence', paste0('\U2265', '5th percentile'), paste0('\U2265', '10th percentile')),
			# # # labAdjX=0.52,
			# # # cex=0.42,
			# # # boxBg=alpha('white', 0.8)
		# # # )
		
		# # # usr <- par('usr')
		# # # width <- usr[2] - usr[1]
		# # # height <- usr[4] - usr[3]
		# # # x <- usr[1] + 0.01 * width
		# # # y <- usr[3] + 0.17 * height
		
		# # # legend(
			# # # x,
			# # # y,
			# # # legend=c('Training presence', 'Long-term test absence', 'Recent test absence', 'Test presence'),
			# # # pch=c(trainPresPch, testLongTermPch, testShortTermPch, testPresPch),
			# # # col=c(trainPresCol, testLongTermAbsCol, testShortTermAbsCol, testPresCol),
			# # # bty='n',
			# # # title='Surveys',
			# # # cex=0.45,
			# # # pt.cex=0.7
		# # # )
		
		# # # # scale bar
		# # # size <- 50000 # length of scale bar in meters
		# # # x <- usr[2] - size - 0.02 * width
		# # # x <- c(x, x + size)
		# # # y <- usr[3] + rep(0.02 * height, 2)
		
		# # # lines(x, y, lwd=4, xpd=NA, col='black', lend=1)
		
		# # # x <- mean(x)
		# # # y <- usr[3] + 0.05 * height
		# # # text(x, y[1], labels=paste(size / 1000, 'km'), cex=0.5)

		# # # title(sub=date(), cex.sub=0.3, outer=TRUE, line=-1)
		
	# # # dev.off()
	
	# iucn <- vect(paste0(drive, '/Research Data/IUCN Range Maps/Mammals 2022-05-23/MAMMALS_TERRESTRIAL_ONLY.shp'))
	# iucn <- iucn[iucn$binomial == 'Ochotona princeps', ]
	# iucn <- project(iucn, nam1)

	# # plot
	# png('./Figures & Tables/Gap Sampling Plumas Emphasized with Range Maps & No Survey Sites.png', width=1200, height=1200, res=300)
		
		# par(mar=c(2, 1, 1, 1), cex.axis=0.4, mgp=c(3, 0, 0.2), fig=c(0, 1, 0, 1))
		
		# plot(focus, border='white')
		
		# usr <- par('usr')
		# xs <- pretty(c(usr[1], usr[2]))
		# ys <- pretty(c(usr[3], usr[4]))
		# axis(1, at=xs, tck=0.01, labels=xs, col=NA, col.ticks='black')
		# axis(2, at=ys, tck=0.01, labels=ys, col=NA, col.ticks='black')
		
		# plot(hs, col=grays, legend=FALSE, add=TRUE)
		# plot(lakes, col = 'gray60', border = NA, add = TRUE)

		# iucnSp <- as(iucn, 'Spatial')
		# iucnSp <- sp::spTransform(iucnSp, CRS(proj4string(west1)))
		# iucnSp <- gBuffer(iucnSp, width=20*1000)
		# iucnSp <- gBuffer(iucnSp, width=-20*1000)
		# plot(iucnSp, col=alpha('forestgreen', 0.5), add=TRUE)

		# testSurveysSpEa <- testSurveysSpEa[order(testSurveysSpEa$status), ]
		
		# plot(kdeClass, col=kdeCols, legend=FALSE, add=TRUE)
		# plot(west2, border='gray40', add=TRUE)
		# plot(west1, border='gray40', lwd=2, add=TRUE)
		# plot(west2[west2$NAME_2 == 'Plumas', ], border='gray40', lwd=3, add=TRUE)

		# # cols <- rep(NA, nrow(testSurveys))
		# # cols[testSurveysSpEa$status == '0 long absence'] <- testLongTermAbsCol
		# # cols[testSurveysSpEa$status == '1 recent absence'] <- testShortTermAbsCol
		# # cols[testSurveysSpEa$status == '2 detected'] <- testPresCol
		
		# # pchs <- rep(NA, nrow(testSurveys))
		# # pchs[testSurveysSpEa$status == '0 long absence'] <- testLongTermPch
		# # pchs[testSurveysSpEa$status == '1 recent absence'] <- testShortTermPch
		# # pchs[testSurveysSpEa$status == '2 detected'] <- testPresPch
		
		# points(trainPresSpEa, pch=trainPresPch, cex=0.5, col=trainPresCol)
		# # points(testSurveysSpEa, pch=pchs, col=cols, cex=0.6)
		# box()
		
		# # legend
		# legendBreaks(
			# 'bottomleft',
			# inset=0.01,
			# height=0.39,
			# width=0.30,
			# title='Density of previously-\nknown occurrences',
			# titleAdj=c(0.5, 0.9),
			# col=kdeCols,
			# adjX=c(0.08, 0.225),
			# adjY=c(0.35, 0.77),
			# labels=c('\U2265min presence', paste0('\U2265', '1st percentile'), paste0('  \U2265', '2.5th percentile')),
			# labAdjX=0.58,
			# cex=0.57,
			# boxBg=alpha('white', 0.8)
		# )
		
		# usr <- par('usr')
		# width <- usr[2] - usr[1]
		# height <- usr[4] - usr[3]
		# x <- usr[1] + 0.01 * width
		# y <- usr[3] + 0.15 * height
		
		# legend(
			# x,
			# y,
			# legend=c('Previously-known\noccurrence', 'IUCN range'),
			# pch=c(trainPresPch, NA),
			# col=c(trainPresCol, 'black'),
			# fill=c(NA, 'darkolivegreen3'),
			# border=c(NA, 'black'),
			# bty='n',
			# cex=0.57,
			# pt.cex=0.7
		# )
		
		# # scale bar
		# size <- 50000 # length of scale bar in meters
		# x <- usr[2] - size - 0.02 * width
		# x <- c(x, x + size)
		# y <- usr[3] + rep(0.02 * height, 2)
		
		# lines(x, y, lwd=4, xpd=NA, col='black', lend=1)
		
		# x <- mean(x)
		# y <- usr[3] + 0.05 * height
		# text(x, y[1], labels=paste(size / 1000, 'km'), cex=0.7)

		# title(sub=date(), cex.sub=0.3, outer=TRUE, line=-1)
		
		# # inset
		# insetNam <- nam1[nam1$NAME_1 %in% c('California', 'Oregon', 'Baja California', 'Baja California Sur', 'Washington'), ]
		# insetNam <- ext(insetNam)
		# insetNam <- vect(insetNam, crs=crs(nam1))
		# insetNam <- buffer(insetNam, 100 * 1000)
		# insetNam <- ext(insetNam)
		# insetNam <- vect(insetNam, crs=crs(nam1))

		# par(fig = c(0.6, 1, 0.0, 0.65), bg='white', new=TRUE)
		# plot(insetNam, col='white', axes=FALSE, bty='o')
		# plot(crop(nam1, insetNam), col='gray80', lwd=0.1, add=TRUE)
		# plot(crop(iucn, insetNam), col='forestgreen', lwd=0.2, add=TRUE)
		
		# counties <- nam2[nam2$NAME_2 %in% gapCounties$NAME_2, ]
		# counties <- ext(counties)
		# foc <- vect(counties, crs=crs(nam2))
		# plot(foc, lwd=2, add=TRUE)

		
	# dev.off()

# say('########################################## NOT REDONE in V2!!!')
# say('### map of gap sampling across Sierras ### NOT REDONE in V2!!!')
# say('########################################## NOT REDONE in V2!!!')

	# # generalize
	# focusBuff <- 1 # buffer size around test sites for generating focus of plot (in km)

	# # spatial data
	# load('./Study Region/GADM California, Nevada, Oregon States.rda')
	# load('./Study Region/GADM California, Nevada, Oregon Counties.rda')
	# load('./Study Region/GADM Plumas County.rda')
	# load('./Study Region/GADM GAP Counties.rda')
	# hs <- raster('./Data/Topography - SRTM/hillshade_srtm_ea.tif')
	
	# gapCounties <- sp::spTransform(gapCounties, getCRS('climateNA', TRUE))
	# plumas <- sp::spTransform(plumas, getCRS('climateNA', TRUE))
	# west1 <- sp::spTransform(west1, getCRS('climateNA', TRUE))
	# west2 <- sp::spTransform(west2, getCRS('climateNA', TRUE))

	# # training occurrences
	# load('./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.rda')
	# trainPres <- surveys[surveys$origRecentPikaOrSignsObserved, ]
	# trainPresSp <- SpatialPointsDataFrame(trainPres[ , ll], data=trainPres, proj4=getCRS('wgs84', TRUE))
	# trainPresSpEa <- sp::spTransform(trainPresSp, getCRS('climateNA', TRUE))

	# # test occurrences
	# testSurveys <- read.csv('./Data/Occurrences/Test Surveys 02 Cleaned.csv')
	# testSurveysSp <- SpatialPointsDataFrame(testSurveys[ , ll], data=testSurveys, proj4=getCRS('nad83', TRUE))
	# testSurveysSpEa <- sp::spTransform(testSurveysSp, getCRS('climateNA', TRUE))
	
	# # plot focus
	# focus <- gBuffer(trainPresSpEa, width=1000 * focusBuff)
	
	# ### crop hillshade
	# plot(focus)
	# usr <- par('usr')
	# dev.off()

	# ext <- extent(usr)
	# ext <- as(ext, 'SpatialPolygons')
	# projection(ext) <- getCRS('climateNA')

	# hs <- crop(hs, ext)
	
	# # kde
	# kde <- raster('./KDE/kde.tif')
	# kde <- projectRaster(kde, crs=getCRS('climateNA'))
	# kdeVals <- extract(kde, trainPresSpEa)
	# quants <- quantile(kdeVals, c(0.05, 0.1))
	# breaks <- c(0, min(kdeVals), quants, 1)
	# kdeClass <- cut(kde, breaks=breaks)

	# # colors
	# kdeCols <- c(NA, 'deepskyblue', 'dodgerblue1', 'dodgerblue4')
	# for (i in seq_along(kdeCols)) kdeCols[i] <- alpha(kdeCols[i], 0.5)
	
	# ### shapes and colors for points
	# trainPresCol <- 'black'
	# testPresCol <- 'darkgreen'
	# testShortTermAbsCol <- 'darkorange4'
	# testLongTermAbsCol <- 'darkred'
	
	# trainPresPch <- 3
	# testPresPch <- 1
	# testShortTermPch <- 5
	# testLongTermPch <- 6

	# # plot
	# png('./Figures & Tables/Gap Sampling - Sierra Nevada.png', width=1200, height=1200, res=300)
		
		# par(mar=c(2, 1, 1, 1), cex.axis=0.4, mgp=c(3, 0, 0.2))
		
		# plot(focus, border='white')
		
		# usr <- par('usr')
		# xs <- pretty(c(usr[1], usr[2]))
		# ys <- pretty(c(usr[3], usr[4]))
		# axis(1, at=xs, tck=0.01, labels=xs, col=NA, col.ticks='black')
		# axis(2, at=ys, tck=0.01, labels=ys, col=NA, col.ticks='black')
		
		# plot(hs, col=grays, legend=FALSE, add=TRUE)
		# plot(kdeClass, col=kdeCols, legend=FALSE, add=TRUE)
		# plot(west2, border='gray40', lwd=0.8, add=TRUE)
		# plot(west1, border='gray40', lwd=1.6, add=TRUE)

		# testSurveysSpEa <- testSurveysSpEa[order(testSurveysSpEa$status), ]
		
		# cols <- rep(NA, nrow(testSurveys))
		# cols[testSurveysSpEa$status == '0 long absence'] <- testLongTermAbsCol
		# cols[testSurveysSpEa$status == '1 recent absence'] <- testShortTermAbsCol
		# cols[testSurveysSpEa$status == '2 detected'] <- testPresCol
		
		# pchs <- rep(NA, nrow(testSurveys))
		# pchs[testSurveysSpEa$status == '0 long absence'] <- testLongTermPch
		# pchs[testSurveysSpEa$status == '1 recent absence'] <- testShortTermPch
		# pchs[testSurveysSpEa$status == '2 detected'] <- testPresPch
		
		# points(trainPresSpEa, pch=trainPresPch, cex=0.5, col=trainPresCol)
		# points(testSurveysSpEa, pch=pchs, col=cols, cex=0.6)
		# box()
		
		# # legend
		# legendBreaks(
			# 'bottomleft',
			# inset=0.01,
			# height=0.44,
			# width=0.27,
			# title='Training occurrence\ndensity',
			# titleAdj=c(0.5, 0.92),
			# col=kdeCols,
			# adjX=c(0.05, 0.225),
			# adjY=c(0.38, 0.85),
			# labels=c('\U2265min presence', paste0('\U2265', '5th percentile'), paste0('\U2265', '10th percentile')),
			# labAdjX=0.52,
			# cex=0.42,
			# boxBg=alpha('white', 0.8)
		# )
		
		# usr <- par('usr')
		# width <- usr[2] - usr[1]
		# height <- usr[4] - usr[3]
		# x <- usr[1] + 0.01 * width
		# y <- usr[3] + 0.17 * height
		
		# legend(
			# x,
			# y,
			# legend=c('Training presence', 'Long-term test absence', 'Recent test absence', 'Test presence'),
			# pch=c(trainPresPch, testLongTermPch, testShortTermPch, testPresPch),
			# col=c(trainPresCol, testLongTermAbsCol, testShortTermAbsCol, testPresCol),
			# bty='n',
			# title='Surveys',
			# cex=0.45,
			# pt.cex=0.7
		# )
		
		# # scale bar
		# size <- 50000 # length of scale bar in meters
		# x <- usr[2] - size - 0.02 * width
		# x <- c(x, x + size)
		# y <- usr[3] + rep(0.02 * height, 2)
		
		# lines(x, y, lwd=4, xpd=NA, col='black', lend=1)
		
		# x <- mean(x)
		# y <- usr[3] + 0.05 * height
		# text(x, y[1], labels=paste(size / 1000, 'km'), cex=0.5)

		# title(sub=date(), cex.sub=0.3, outer=TRUE, line=-1)
		
	# dev.off()

# say('############################################################')
# say('### map of gap sampling across Sierras sans survey plots ###')
# say('############################################################')

	# # generalize
	# focusBuff <- 1 # buffer size around test sites for generating focus of plot (in km)

	# # spatial data
	# load('./Study Region/GADM California, Nevada, Oregon States.rda')
	# load('./Study Region/GADM California, Nevada, Oregon Counties.rda')
	# load('./Study Region/GADM Plumas County.rda')
	# load('./Study Region/GADM GAP Counties.rda')
	# hs <- raster('./Data/Topography - SRTM/hillshade_srtm_ea.tif')
	
	# gapCounties <- sp::spTransform(gapCounties, getCRS('climateNA', TRUE))
	# plumas <- sp::spTransform(plumas, getCRS('climateNA', TRUE))
	# west1 <- sp::spTransform(west1, getCRS('climateNA', TRUE))
	# west2 <- sp::spTransform(west2, getCRS('climateNA', TRUE))

	# # training occurrences
	# load('./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.rda')
	# trainPres <- surveys[surveys$origRecentPikaOrSignsObserved, ]
	# trainPresSp <- SpatialPointsDataFrame(trainPres[ , ll], data=trainPres, proj4=getCRS('wgs84', TRUE))
	# trainPresSpEa <- sp::spTransform(trainPresSp, getCRS('climateNA', TRUE))

	# # plot focus
	# focus <- gBuffer(trainPresSpEa, width=1000 * focusBuff)
	
	# ### crop hillshade
	# plot(focus)
	# usr <- par('usr')
	# dev.off()

	# ext <- extent(usr)
	# ext <- as(ext, 'SpatialPolygons')
	# projection(ext) <- getCRS('climateNA')

	# hs <- crop(hs, ext)
	
	# # kde
	# kde <- raster('./KDE/kde.tif')
	# kde <- projectRaster(kde, crs=getCRS('climateNA'))
	# kdeVals <- extract(kde, trainPresSpEa)
	# quants <- quantile(kdeVals, c(0.05, 0.1))
	# breaks <- c(0, min(kdeVals), quants, 1)
	# kdeClass <- cut(kde, breaks=breaks)

	# # colors
	# kdeCols <- c(NA, 'deepskyblue', 'dodgerblue1', 'dodgerblue4')
	# for (i in seq_along(kdeCols)) kdeCols[i] <- alpha(kdeCols[i], 0.5)
	
	# ### shapes and colors for points
	# trainPresCol <- 'black'
	# trainPresPch <- 3

	# # plot
	# png('./Figures & Tables/Gap Sampling - Sierra Nevada Sans Survey Plots.png', width=1200, height=1200, res=300)
		
		# par(mar=c(2, 1, 1, 1), cex.axis=0.4, mgp=c(3, 0, 0.2))
		
		# plot(focus, border='white')
		
		# usr <- par('usr')
		# xs <- pretty(c(usr[1], usr[2]))
		# ys <- pretty(c(usr[3], usr[4]))
		# axis(1, at=xs, tck=0.01, labels=xs, col=NA, col.ticks='black')
		# axis(2, at=ys, tck=0.01, labels=ys, col=NA, col.ticks='black')
		
		# plot(hs, col=grays, legend=FALSE, add=TRUE)
		# plot(kdeClass, col=kdeCols, legend=FALSE, add=TRUE)
		# plot(west2, border='gray40', lwd=0.8, add=TRUE)
		# plot(west1, border='gray40', lwd=1.6, add=TRUE)

		# points(trainPresSpEa, pch=trainPresPch, cex=0.5, col=trainPresCol)
		# box()
		
		# # legend
		# legendBreaks(
			# 'bottomleft',
			# inset=0.01,
			# height=0.24,
			# width=0.24,
			# title='Training occurrence\ndensity',
			# titleAdj=c(0.5, 0.87),
			# col=kdeCols,
			# adjX=c(0.05, 0.225),
			# adjY=c(0.05, 0.80),
			# labels=c('\U2265min presence', paste0('\U2265', '5th percentile'), paste0('\U2265', '10th percentile')),
			# labAdjX=0.59,
			# cex=0.42,
			# boxBg=alpha('white', 0.8)
		# )
		
		# usr <- par('usr')
		# width <- usr[2] - usr[1]
		# height <- usr[4] - usr[3]

		# # scale bar
		# size <- 50000 # length of scale bar in meters
		# x <- usr[2] - size - 0.02 * width
		# x <- c(x, x + size)
		# y <- usr[3] + rep(0.02 * height, 2)
		
		# lines(x, y, lwd=4, xpd=NA, col='black', lend=1)
		
		# x <- mean(x)
		# y <- usr[3] + 0.05 * height
		# text(x, y[1], labels=paste(size / 1000, 'km'), cex=0.5)

		# title(sub=date(), cex.sub=0.3, outer=TRUE, line=-1)
		
	# dev.off()

# say('#########################################################################')
# say('### assign weights to TRAINING sites based on spatial autocorrelation ###')
# say('#########################################################################')

	# say('I want to use inverse-p weighting applied to survey sites to correct for spatial sampling bias. I will assume that two survey sites at the exact same location each have a weight of 1/2, three each have a weight of 1/3, and so on. At the other end of the spectrum survey sites that are far enough away should each have a weight of 1. I will define "far enough away" as the distance at which the proportion of pairwise observed distances falls below the upper 95th quantile of the distribution of pairwise distances from randomly located sites (one-tail). I will draw a number of randomly located sites in each iteration so it is equal to the total number of survey sites.', breaks=80)

	# ### generalization
	# minSurveys <- 0 # minimum number of survey sites by one provider to use in SAC assessment
	
	# ### mask
	# mask <- raster(paste0('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus', studyRegionBuffer, 'kmBuffer.tif'))

	# ### occurrences
	# load('./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.rda')
	# trainPres <- surveys[surveys$origRecentPikaOrSignsObserved, ]
	
	# trainPres$provider <- NA
	# trainPres$provider[trainPres$contact %in% c('Chris Ray (Colorado University, Boulder)')] <- 'Chris Ray'
	# trainPres$provider[trainPres$contact %in% c('Clint Epps, Jessica Castillo', 'Jess Castillo, Clint Epps')] <- 'Clint Epps & Jessica Castillo'
	# trainPres$provider[trainPres$contact %in% c('Connie Millar')] <- 'Connie Millar'
	# trainPres$provider[trainPres$contact %in% c('Erik Beever', 'Erik Beever;')] <- 'Erik Beever'
	# trainPres$provider[trainPres$contact %in% c('Erik Beever, Chris Ray')] <- 'Erik Beever & Chris Ray'
	# trainPres$provider[trainPres$contact %in% c('Adam Smith;')] <- 'online database'
	# trainPres$provider[trainPres$contact %in% c('Jason Brewer & Mary Flores')] <- 'Jason Brewer & Mary Flores'
	# trainPres$provider[trainPres$contact %in% c('Ken Goehring;')] <- 'Ken Goehring'
	# trainPres$provider[trainPres$contact %in% c('Joseph Stewart')] <- 'Joseph Stewart'
	# trainPres$provider[trainPres$contact %in% c('Joseph Stewart, David Wright, Chris Curlis')] <- 'Joseph Stewart, David Wright, Chris Curlis'

	# providers <- sort(unique(trainPres$provider))
	
	# ### generate breaks for distance bins
	# # maximum possible distance (across diagonal of study region mask)
	# ext <- extent(mask)
	# ext <- as(ext, 'SpatialPolygons')
	# projection(ext) <- projection(mask)
	# ext <- sp::spTransform(ext, getCRS('climateNA', TRUE))
	# ext <- extent(ext)
	# corners <- cbind(x=c(ext@xmin, ext@xmax), y=c(ext@ymin, ext@ymax))
	# corners <- SpatialPoints(corners, getCRS('climateNA', TRUE))
	# corners <- sp::spTransform(corners, CRS(projection(mask)))
	# maxDist <- max(distm(corners))
	
	# maxDist <- sacInterval * ceiling(maxDist / sacInterval)
	
	# # breaks
	# breaks <- data.frame(
		# lower=seq(0, maxDist - sacInterval, by=sacInterval / 2),
		# upper=seq(sacInterval, maxDist, by=sacInterval / 2)
	# )
	
	# breaks <- as.matrix(breaks)
	
	# dirCreate('./Figures & Tables/Spatial Autocorrelation between Survey Sites')
	
	# # stores minimum distance of SAC for each provider
	# sacDist <- data.frame()

	# ### evaluate SAC by provider
	# ### assign weights
	
	# trainPres$weight <- NA
	
	# for (provider in providers) {

		# say(provider, pre=2)
	
		# # provider dists
		# presProviderIndex <- which(trainPres$provider == provider)
		# presProvider <- trainPres[presProviderIndex, ]
		
		# if (nrow(presProvider) >= minSurveys) {
			
			# pts <- presProvider[ , ll]
			
			# nullCorrTrain <- spatialCorrForPoints(
				# pts=pts,
				# rast=mask,
				# breaks=breaks,
				# iters=100,
				# verbose=FALSE
			# )

			# # remember minimum SAC distance
			# thisSacDist <- spatialCorrForPointsSummary(nullCorrTrain)
			# sacDist <- rbind(sacDist, data.frame(provider=provider, sacDist_m=thisSacDist))
			
			# # plot pairwise distance distribution
			# png(paste0('./Figures & Tables/Spatial Autocorrelation between Survey Sites/Training Sites - ', provider, '.png'), width=1200, height=800, res=200)
				# main <- paste0(provider, '\nMinimum distance of SAC: ', thisSacDist / 1000, ' km')
				# spatialCorrForPointsPlot(x=nullCorrTrain, rescale=0.001, xlab='Distance (km)', main=main, cex=0.8)
			# dev.off()
		
			# # assign weights
			# weight <- spatialCorrForPointsWeight(x=nullCorrTrain, pts=pts, verbose=FALSE)
			# trainPres$weight[presProviderIndex] <- weight
		
		# } # if adequate number of occurrences for this provider
		
	# } # next provider

	# write.csv(sacDist, './Figures & Tables/Spatial Autocorrelation between Survey Sites/!Characteristic Clustering Distance for Training Sites.csv', row.names=FALSE)
	
	# save(trainPres, file='./Data/Occurrences/Training Surveys 02 Pika - Assigned Weights Based on Pairwise Distance Distribution.rda')

# say('################################################################################################')
# say('### calculate spatial autocorrelation within test and between TEST and TRAINING survey sites ###')
# say('################################################################################################')

	# say('I want to use inverse-p weighting applied to survey testSites to correct for spatial sampling bias. I will assume that two survey testSites at the exact same location each have a weight of 1/2, three each have a weight of 1/3, and so on. At the other end of the spectrum survey testSites that are far enough away should each have a weight of 1. I will define "far enough away" as the distance at which the proportion of pairwise observed distances falls below the upper 95th quantile of the distribution of pairwise distances from randomly located testSites (one-tail). I will draw a number of randomly located testSites in each iteration so it is equal to the total number of survey testSites.', breaks=80)

	# ### test occurrences
	# testSites <- fread('./Data/Occurrences/Test Surveys 02 Cleaned.csv')
	# testSites <- as.data.frame(testSites)
	# testSitesSp <- SpatialPoints(testSites[ , ll], getCRS('wgs84', TRUE))

	# ### mask
	# mask <- raster(paste0('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus', studyRegionBuffer, 'kmBuffer.tif'))
	
	# # restrict mask to buffer around test test sites within max nearest-neighbor distance between test testSites
	# sitesDist <- distm(testSitesSp)
	# diag(sitesDist) <- NA
	# nearestNeighDist <- apply(sitesDist, 1, min, na.rm=TRUE)
	# maxNearestNeighDist <- max(nearestNeighDist)
	
	# testSitesSpEa <- sp::spTransform(testSitesSp, getCRS('climateNA', TRUE))
	# buff <- gBuffer(testSitesSpEa, width=maxNearestNeighDist)
	# buff <- sp::spTransform(buff, CRS(projection(mask)))
	# mask <- crop(mask, buff)
	# buffRast <- rasterize(buff, mask)
	# mask <- buffRast * mask

	# ### get only training sites within the buffer
	# load('./Data/Occurrences/Training Surveys 02 Pika - Assigned Weights Based on Pairwise Distance Distribution.rda')
	# trainPresSp <- SpatialPoints(trainPres[ , ll], CRS(projection(mask)))
	# trainPresSpEa <- sp::spTransform(trainPresSp, getCRS('climateNA', TRUE))
	
	# inBuff <- over(trainPresSp, buff)
	# trainPres <- trainPres[!is.na(inBuff), ]
	
	# ### generate breaks for distance bins
	# # maximum possible distance (across diagonal of study region mask)
	# ext <- extent(mask)
	# ext <- as(ext, 'SpatialPolygons')
	# projection(ext) <- projection(mask)
	# ext <- sp::spTransform(ext, getCRS('climateNA', TRUE))
	# ext <- extent(ext)
	# corners <- cbind(x=c(ext@xmin, ext@xmax), y=c(ext@ymin, ext@ymax))
	# corners <- SpatialPoints(corners, getCRS('climateNA', TRUE))
	# corners <- sp::spTransform(corners, CRS(projection(mask)))
	# maxDist <- max(distm(corners))

	# maxDist <- sacInterval * ceiling(maxDist / sacInterval)
	
	# # breaks
	# breaks <- data.frame(
		# lower=seq(0, maxDist - sacInterval, by=sacInterval / 2),
		# upper=seq(sacInterval, maxDist, by=sacInterval / 2)
	# )
	
	# breaks <- as.matrix(breaks)
	
	# midDists <- rowMeans(breaks[ , c('lower', 'upper')])

	# dirCreate('./Figures & Tables/Spatial Autocorrelation between Survey Sites')

	# ### SAC within TEST surveys
	# nullCorrTest <- spatialCorrForPoints(
		# pts=testSites[ , ll],
		# rast=mask,
		# breaks=breaks,
		# iters=100,
		# verbose=FALSE
	# )

	# # assign weights
	# weightTest <- spatialCorrForPointsWeight(x=nullCorrTest, pts=testSites[ , ll], verbose=FALSE)
	# testSites$weightTestVsTest <- weightTest

	# sacDistTest <- spatialCorrForPointsSummary(nullCorrTest, verbose=FALSE)

	# # plot
	# png(paste0('./Figures & Tables/Spatial Autocorrelation between Survey Sites/Test Sites.png'), width=1200, height=800, res=200)
		# main <- paste0('Test vs Test Sites\nMinimum distance of SAC: ', 0.001 * sacDistTest, ' km')
		# spatialCorrForPointsPlot(x=nullCorrTest, rescale=0.001, xlab='Distance (km)', main=main, cex=0.8)
	# dev.off()

	# ### SAC within TRAINING surveys
	# nullCorrTestTrain <- spatialCorrForPoints(
		# pts=testSites[ , ll],
		# fixed=trainPres[ , ll],
		# rast=mask,
		# breaks=breaks,
		# iters=100,
		# verbose=FALSE
	# )

	# # assign weights
	# weightTestTrain <- spatialCorrForPointsWeight(x=nullCorrTest, pts=testSites[ , ll], fixed=trainPres[ , ll], verbose=FALSE)
	# testSites$weightTestVsTrain <- weightTestTrain

	# sacDistTestTrain <- spatialCorrForPointsSummary(nullCorrTestTrain, verbose=FALSE)

	# # plot
	# png(paste0('./Figures & Tables/Spatial Autocorrelation between Survey Sites/Test vs Train Sites.png'), width=1200, height=800, res=200)
		# main <- paste0('Test vs Training Sites\nMinimum distance of SAC: ', 0.001 * sacDistTestTrain, ' km')
		# spatialCorrForPointsPlot(x=nullCorrTestTrain, rescale=0.001, xlab='Distance (km)', main=main, cex=0.8)
	# dev.off()

	# out <- data.frame(
		# surveys=c('test', 'test vs train fixed'),
		# sacDist_m=c(sacDistTest, sacDistTestTrain)
	# )
	
	# write.csv(out, file='./Figures & Tables/Spatial Autocorrelation between Survey Sites/!Characteristic Clustering Distance for Test vs Test and Test vs Training Sites.csv')

	# write.csv(testSites, './Data/Occurrences/Test Surveys 03 Assigned Weights Based on Pairwise Spatial Autocorrelation.csv', row.names=FALSE)
	
# say('######################################################')
# say('### plot TRAINING sites scaled by SAC-based weight ###')
# say('######################################################')

	# # spatial data
	# study <- shapefile('./Study Region/epa3_sierraNevadaSarrModified_200kmBuffer.shp')
	# sn <- shapefile('./Study Region/epa3_sierraNevadaSarrModified.shp')
	# load('./Study Region/GADM California, Nevada, Oregon States.rda')
	# load('./Study Region/GADM California, Nevada, Oregon Counties.rda')
	
	# # training occurrences
	# load('./Data/Occurrences/Training Surveys 02 Pika - Assigned Weights Based on Pairwise Distance Distribution.rda')
	# sacDist <- read.csv('./Figures & Tables/Spatial Autocorrelation between Survey Sites/!Characteristic Clustering Distance for Training Sites.csv')
	
	# trainPresSp <- SpatialPointsDataFrame(trainPres[ , ll], data=trainPres, proj4=CRS(prismCrs))

	# providers <- sacDist$provider
	
	# png('./Figures & Tables/Spatial Autocorrelation between Survey Sites/!Maps of Weighting by Provider.png', width=1600, height=800, res=200)
		
		# par(mfrow=c(2, 5), mar=c(1, 1, 4, 1))
		
		# for (countProv in seq_along(providers)) {
			
			# provider <- providers[countProv]
			
			# # add buffer to points to increase focal area
			# pts <- trainPresSp[trainPres$provider == provider, ]
			# ptsEa <- sp::spTransform(pts, getCRS('climateNA'))
			# ptsEaBuff <- gBuffer(ptsEa, width=35000)
			# ptsBuff <- sp::spTransform(ptsEaBuff, CRS(projection(pts)))

			# # plot background
			# plot(ptsBuff, border='white', main=provider)
			# plot(west2, border='gray70', add=TRUE)
			# plot(west1, border='black', add=TRUE)

			# # points
			# alpha <- pts$weight
			# alpha <- alpha / max(trainPres$weight)
			# col <- alpha('darkblue', alpha)
			# points(pts[ , ll], bg=col, col=alpha('darkblue', 0.5), pch=21, cex=1.2)

			# # scale circle
			# usr <- par('usr')
			
			# provSacDist <- sacDist$sacDist[sacDist[ , 1] == provider]

			# center <- cbind(mean(usr[1:2]), mean(usr[3:4]))
			# center <- SpatialPoints(center, getCRS('wgs84', TRUE))
			# centerEa <- sp::spTransform(center, getCRS('climateNA', TRUE))
			# centerEaBuff <- gBuffer(centerEa, width=0.5 * provSacDist)
			# centerBuff <- sp::spTransform(centerEaBuff, getCRS('wgs84', TRUE))
			
			# plot(centerBuff, lty='dashed', add=TRUE)
			
		# }
		
		# title(sub=date(), outer=TRUE, line=-1, cex.sub=0.4)
		
	# dev.off()

# say('##################################################')
# say('### plot TEST sites scaled by SAC-based weight ###')
# say('##################################################')

	# # spatial data
	# study <- shapefile('./Study Region/epa3_sierraNevadaSarrModified_200kmBuffer.shp')
	# sn <- shapefile('./Study Region/epa3_sierraNevadaSarrModified.shp')
	# load('./Study Region/GADM California, Nevada, Oregon States.rda')
	# load('./Study Region/GADM California, Nevada, Oregon Counties.rda')
	
	# # test occurrences
	# sites <- fread('./Data/Occurrences/Test Surveys 03 Assigned Weights Based on Pairwise Spatial Autocorrelation.csv')
	# sites <- as.data.frame(sites)
	# sacDist <- read.csv('./Figures & Tables/Spatial Autocorrelation between Survey Sites/!Characteristic Clustering Distance for Test vs Test and Test vs Training Sites.csv')
	
	# sitesSp <- SpatialPointsDataFrame(sites[ , ll], data=sites, proj4=CRS(prismCrs))

	# hs <- raster('./Data/Topography - SRTM/hillshade_srtm.tif')
	
	# # training occurrences
	# load('./Data/Occurrences/Training Surveys 02 Pika - Assigned Weights Based on Pairwise Distance Distribution.rda')
	
	# png('./Figures & Tables/Spatial Autocorrelation between Survey Sites/!Maps of Weighting for Test Sites.png', width=1600, height=800, res=200)
		
		# par(mfrow=c(1, 2), mar=c(1, 1, 4, 1))
		
		# for (set in c('test', 'test vs train fixed')) {
			
			# weightCol <- if (set == 'test') {
				# 'weightTestVsTest'
			# } else if (set == 'test vs train fixed') {
				# 'weightTestVsTrain'
			# }
			
			# # add buffer to points to increase focal area
			# sitesSpEa <- sp::spTransform(sitesSp, getCRS('climateNA'))
			# sitesEaBuff <- gBuffer(sitesSpEa, width=25000)
			# sitesBuff <- sp::spTransform(sitesEaBuff, CRS(projection(sitesSp)))

			# # plot background
			# thisSacDist <- sacDist$sacDist_m[sacDist$surveys == set]

			# plot(sitesBuff, border='white', main=capIt(set))
			# plot(hs, col=grays, legend=FALSE, add=TRUE)
			# plot(west2, border='gray90', add=TRUE)
			# plot(west1, border='black', add=TRUE)

			# # training sites
			# points(trainPres[ , ll], pch=0, cex=1)
			
			# # test points
			# alpha <- sites[ , weightCol]
			# alpha <- alpha / max(alpha)
			# col <- alpha('darkblue', alpha)
			# points(sites[ , ll], bg=col, col=alpha('darkblue', 0.5), pch=21, cex=1.4)

			# # scale circle
			# usr <- par('usr')
			
			# center <- cbind(mean(usr[1:2]), mean(usr[3:4]))
			# center <- SpatialPoints(center, getCRS('wgs84', TRUE))
			# centerEa <- sp::spTransform(center, getCRS('climateNA', TRUE))
			# centerEaBuff <- gBuffer(centerEa, width=0.5 * thisSacDist)
			# centerBuff <- sp::spTransform(centerEaBuff, getCRS('wgs84', TRUE))
			
			# plot(centerBuff, add=TRUE)
			
			# legend('bottomleft', legend=c('Training site', 'Test site'), pch=c(0, 21), pt.bg=c(NA, alpha('darkblue', 0.3)), col=c('black', alpha('darkblue', 0.3)), bty='n')
			
		# }
		
		# title(sub=date(), outer=TRUE, line=-1, cex.sub=0.4)
		
	# dev.off()
		
# say('###################################################')
# say('### generate background sites and calculate PCA ###')
# say('###################################################')

	# dirCreate('./ENMs/Background Sites')

	# mask <- raster('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus70kmBuffer.tif')
	# set.seed(123)
	# bgSites <- randomPoints(mask, 11000)
	
	# env <- stackEnv()

	# bgEnv <- extract(env, bgSites)
	
	# colnames(bgSites) <- c('longWgs84', 'latWgs84')
	
	# bg <- cbind(bgSites, bgEnv)
	# bg <- as.data.frame(bg)
	
	# nas <- naRows(bg)
	# if (length(nas) > 0) bg <- bg[-nas, ]
	# bg <- bg[1:10000, ]
	
	# rownames(bg) <- 1:10000
	
	# # calculate PCA
	# pca <- princomp(bg[ , names(env)], cor=TRUE)
	# pcs <- pca$scores
	# colnames(pcs) <- paste0('pc', 1:ncol(bgEnv))
	# bg <- cbind(bg, pcs)
	
	# save(bg, file='./ENMs/Background Sites/Random Background Sites from across Study Region Mask.rda')
	# save(pca, file='./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.rda')

# say('##################')
# say('### PCA biplot ###')
# say('##################')

	# # # pcaBiplot <- function(x, pcs=c(1, 2)) {
		# # # data <- as.data.frame(x$scores)
		# # # plot <- ggplot(data, aes_string(x=data[ , pcs[1]], y=data[ , pcs[1]])) +
		# # # # + geom_text(alpha=.4, size=3, aes(label=obsnames)) +
		# # # geom_hline(yintercept=0, size=0.2) + geom_vline(xintercept=0, size=0.2)
		
		# # # datapc <- as.data.frame(varnames=rownames(x$loadings), x$rotation)
		# # # mult <- min(
			# # # (max(data[, pcs[2]]) - min(data[, pcs[1]])/(max(datapc[,y])-min(datapc[,y]))),
			# # # (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
			# # # )
		# # # datapc <- transform(datapc,
				# # # v1 = .7 * mult * (get(x)),
				# # # v2 = .7 * mult * (get(y))
				# # # )
		# # # plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
		# # # plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
		# # # plot
	# # # }

	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.rda')
	
	# totalVar <- sum(pca$sdev^2)
	# varPc1 <- pca$sdev[1]^2 / totalVar
	# varPc2 <- pca$sdev[2]^2 / totalVar
	
	# xlab <- paste0('PC1 (', sprintf('%.1f', 100 * varPc1), '%)')
	# ylab <- paste0('PC2 (', sprintf('%.1f', 100 * varPc2), '%)')
	
	# png('./Figures & Tables/Background Environment as PCA.png', width=1000, height=1000, res=200)
	
		# # par(mar=c(4, 4, 4, 4), cex.lab=0.6, cex.axis=0.6)
		
		# # plot(pca$scores[ , 1:2], type='p', pch=16, col=alpha('darkblue', 0.05), xlab=xlab, ylab=ylab)
		
		# # stretch <- 3
		# # for (i in 1:nrow(pca$loadings)) {
		
			# # arrows(x0=0, y0=0, x1=stretch * pca$loadings[i, 1], y1=stretch * pca$loadings[i, 2], angle=15, length=0.07)
			# # text(stretch * pca$loadings[i, 1], stretch * pca$loadings[i, 2], labels=rownames(pca$loadings)[i])
			
		# # }
		
		# biplot(pca, xlab=xlab, ylab=ylab, pch=1, xpd=NA)
		
	# dev.off()
	
	# sink('./Figures & Tables/PCA.txt', split=TRUE)
	
		# print(pca$loadings)
		# say('')
		# print(summary(pca))
		
	# sink()
	
# say('#####################################################')
# say('### match predictors with TRAINING and TEST sites ###')
# say('#####################################################')

	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.rda')

	# load('./Data/Occurrences/Training Surveys 02 Pika - Assigned Weights Based on Pairwise Distance Distribution.rda')
	# testSurveys <- fread('./Data/Occurrences/Test Surveys 03 Assigned Weights Based on Pairwise Spatial Autocorrelation.csv')
	# testSurveys <- as.data.frame(testSurveys)
	
	# env <- stackEnv()
	
	# trainEnv <- extract(env, trainPres[ , ll])
	# testEnv <- extract(env, testSurveys[ , ll])
	
	# trainPres <- cbind(trainPres, trainEnv)
	# testSurveys <- cbind(testSurveys, testEnv)
	
	# trainPc <- predict(pca, trainEnv)
	# testPc <- predict(pca, testEnv)
	
	# colnames(trainPc) <- paste0('pc', 1:ncol(trainEnv))
	# colnames(testPc) <- paste0('pc', 1:ncol(testEnv))
	
	# trainPres <- cbind(trainPres, trainPc)
	# testSurveys <- cbind(testSurveys, testPc)
	
	# save(trainPres, file='./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.rda')
	# write.csv(testSurveys, './Data/Occurrences/Test Surveys 04 Extracted Environmental Values with PCA.csv', row.names=FALSE)
	
# say('##################')
# say('### train ENMs ###')
# say('##################')	

	# # training presences and background
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.rda')
	# load('./ENMs/Background Sites/Random Background Sites from across Study Region Mask.rda')
	
	# # PCA
	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.rda')
	
	# bg$weight <- 1
	
	# pcs <- paste0('pc', whichPcs(pca))
	
	# ### collate training sites
	# train <- rbind(
		# trainPres[ , c(pcs, 'weight')],
		# bg[ , c(pcs, 'weight')]
	# )
	
	# presBg <- data.frame(presBg = c(rep(1, nrow(trainPres)), rep(0, nrow(bg))))
	# train <- insertCol(presBg, into=train, at=1)

	# # equalize total weight of presences and background sites
	# train$weight[train$presBg == 1] <- train$weight[train$presBg == 1] / max(train$weight[train$presBg == 1])
	# presWeight <- sum(train$weight[train$presBg == 1])
	# bgWeight <- sum(train$weight[train$presBg == 0])
	
	# ratioPresBg <- presWeight / bgWeight
	
	# if (ratioPresBg > 1) {
		# train$weight[train$presBg == 0] <- ratioPresBg * train$weight[train$presBg == 0]
	# } else {
		# train$weight[train$presBg == 1] <- train$weight[train$presBg == 1] / ratioPresBg
	# }

	# train$weight <- train$weight / max(train$weight)
	
	# ### BRT
	# set.seed(123)
	# say('brt')
	# model <- enmSdmX::trainBRT(train, resp='presBg', preds=pcs, w='weight', verbose=TRUE)
	# save(model, file='./ENMs/ENM BRT.rda')
	# rm(model); gc()
	
	# ### GLM
	# say('glm')
	# model <- enmSdmX::trainGLM(train, resp='presBg', preds=pcs, w='weight', interaction = TRUE, quadratic = TRUE, verbose=TRUE)
	# save(model, file='./ENMs/ENM GLM.rda')
	# rm(model); gc()
	
	# ### GAM
	# say('gam')
	# model <- enmSdm::trainGam(train, resp='presBg', preds=pcs, w='weight', verbose=TRUE)
	# save(model, file='./ENMs/ENM GAM.rda')
	# rm(model); gc()

# say('###################################')
# say('### assess predictor importance ###')
# say('###################################')

	# # PCA
	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.rda')
	# pcs <- paste0('pc', whichPcs(pca))
	
	# # training presences and background
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.rda')
	# load('./ENMs/Background Sites/Random Background Sites from across Study Region Mask.rda')
	
	# # PCA
	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.rda')
	
	# bg$weight <- 1
	
	# niter <- 100
	
	# pcs <- paste0('pc', whichPcs(pca))
	
	# ### collate training sites
	# train <- rbind(
		# trainPres[ , pcs],
		# bg[ , pcs]
	# )

	# load(paste0('./ENMs/ENM BRT.rda'))
	# brt <- model
				
	# load(paste0('./ENMs/ENM GAM.rda'))
	# gam <- model
				
	# load(paste0('./ENMs/ENM GLM.rda'))
	# glm <- model
	
	# rm(model); gc()
				
	# obsPred <- predictEnsemble(train, brt=brt, gam=gam, glm=glm)

	# varImp <- matrix(NA, ncol=length(pcs), nrow=niter)
	# colnames(varImp) <- pcs
	
	# set.seed(123)
	# for (pc in seq_along(pcs)) {
	
		# say(paste0('pc', pc))
	
		# corTest <- rep(NA, niter)
		# for (iter in 1:niter) {

			# trainRand <- train
			# trainRand[ , paste0('pc', pc)] <- sample(trainRand[ , paste0('pc', pc)], nrow(trainRand))
			# randPred <- predictEnsemble(trainRand, brt=brt, gam=gam, glm=glm)
			
			# corTest[iter] <- cor(obsPred, randPred)
			
		# }

		# corTest <- 1 - corTest
		# varImp[ , paste0('pc', pc)] <- corTest
		
	# } # next PC
	
	# write.csv(varImp, './Figures & Tables/Predictor Importance in Ensemble Model.csv')

# say('#########################')
# say('### write ENM rasters ###')
# say('#########################')	

	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.rda')

	# # stack of PC rasters
	# env <- stackEnv()
	
	# envPc <- predict(env, pca, index=1:nlayers(env))
	# names(envPc) <- paste0('pc', 1:nlayers(envPc))
	# envPc <- subset(envPc, whichPcs(pca))

	# # predict
	# if (exists('predStack')) rm(predStack)
	
	# for (algo in toupper(algos)) {
	
		# say(algo)
	
		# load(paste0('./ENMs/ENM ', algo, '.rda'))
		# if (algo == 'CRF') {
			# envPcDf <- as.data.frame(envPc)

			# crfPred <- predict(model, newdata=envPcDf, type='prob')
			# out <- envPc[[1]]
			# out[] <- crfPred
			# rm(envPcDf, crfPred); gc()
		# } else {
			# out <- predict(envPc, model, type='response', n.trees=model$gbm.call$n.trees)
		# }
		
		# names(out) <- tolower(algo)
		
		# predStack <- if (!exists('predStack')) {
			# out
		# } else {
			# stack(predStack, out)
		# }
	
		# rm(model); gc()
	
	# }
	
	# pred <- mean(predStack)
	
	# pred <- 1000 * pred
	# pred <- round(pred)
	# names(pred) <- 'ensemblePrediction'
	
	# predStack <- 1000 * predStack
	# predStack <- round(predStack)
	# names(predStack) <- algos
	
	# writeRaster(pred, './ENMs/predictionEnsemble', datatype='INT2U', overwrite=TRUE)
	# writeRaster(predStack, './ENMs/predictionByAlgorithm', datatype='INT2U', overwrite=TRUE)

# say('################################################')
# say('### extract predictions to TEST survey sites ###')
# say('################################################')	

	# ensPred <- raster('./ENMs/predictionEnsemble.tif')
	# algosPred <- stack('./ENMs/predictionByAlgorithm.tif')
	# names(algosPred) <- algos
	
	# testSurveys <- fread('./Data/Occurrences/Test Surveys 04 Extracted Environmental Values with PCA.csv')
	# testSurveys <- as.data.frame(testSurveys)
	
	# predictEnsemble <- extract(ensPred, testSurveys[ , ll]) / 1000
	# predictAlgos <- extract(algosPred, testSurveys[ , ll]) / 1000
	
	# predicts <- cbind(predictEnsemble, predictAlgos)
	
	# testSurveys <- cbind(testSurveys, predicts)
	
	# write.csv(testSurveys, './Data/Occurrences/Test Surveys 05 Extracted Predictions.csv', row.names=FALSE)
	
# say('##################################################################')
# say('### calculate performance statistics against TEST survey sites ###')
# say('##################################################################')	

	# testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# testSurveys <- as.data.frame(testSurveys)
	# testSurveys$weight <- pmin(testSurveys$weightTestVsTest, testSurveys$weightTestVsTrain)
	# longTermAbs <- testSurveys[testSurveys$status=='0 long absence', c('predictEnsemble', 'weight')]
	# recentAbs <- testSurveys[testSurveys$status=='1 recent absence' , c('predictEnsemble', 'weight')]
	# pres <- testSurveys[testSurveys$status=='2 detected' , c('predictEnsemble', 'weight')]

	# aucs <- aucMultiWeighted(longTermAbs, recentAbs, pres)

	# numLongTermAbs <- nrow(longTermAbs)
	# numRecentAbs <- nrow(recentAbs)
	# numPres <- nrow(pres)
		
	# sink('./Figures & Tables/Model Performance - AUC.txt', split=TRUE)
		# say('AUC:')
		
		# say('Multivariate: ........................ ', sprintf('%.3f', aucs[['multivariate']]))
		# say('Presences over recent absences: ...... ', sprintf('%.3f', aucs[['case3_over_case2']]))
		# say('Presences over long-term absences: ... ', sprintf('%.3f', aucs[['case3_over_case1']]))
		# say('Recent over long-term absences: ...... ', sprintf('%.3f', aucs[['case2_over_case1']]), post=2)
		
		# say('Number of presences: ..................', numPres)
		# say('Number of recent absences: ............', numRecentAbs)
		# say('Number of long-term absences: .........', numLongTermAbs)
		
		# say(date())
	# sink()
	
# say('############################################################')
# say('### boxplot of distribution of predictions at TEST sites ###')
# say('############################################################')

	# testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# testSurveys <- as.data.frame(testSurveys)
	
	# testSurveys$weight <- pmin(testSurveys$weightTestVsTest, testSurveys$weightTestVsTrain)
	
	# png('./Figures & Tables/Predicted Suitability by Test Site Class - Unweighted.png', width=1200, height=1200, res=300)
		# par(mar=c(3, 4, 2, 1) + 0.1, cex.lab=0.8, cex.axis=0.8, cex.sub=0.4)
		# boxplot(predictEnsemble ~ status, data=testSurveys, names=c('No\nevidence', 'Formerly\noccupied', 'Currently\noccupied'), ylim=c(0, 1), ylab='Predicted suitability')
		# title(sub=date(), line=2)
	# dev.off()
	
	# png('./Figures & Tables/Predicted Suitability by Test Site Class - Weighted.png', width=1200, height=1200, res=300)
		# par(mar=c(3, 4, 2, 1) + 0.1, cex.lab=0.8, cex.axis=0.8, cex.sub=0.4)
		# boxplot(predictEnsemble * weight ~ status, data=testSurveys, names=c('No\nevidence', 'Formerly\noccupied', 'Currently\noccupied'), ylab='Predicted suitability (weighted)')
		# title(sub=date(), line=2)
	# dev.off()
	
# say('###############################################')
# say('### PCA of background and TEST site classes ###')
# say('###############################################')
	
	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.rda')
	
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.rda')
	# trainPres$origUsable <- as.logical(trainPres$origUsable)
	# trainPres <- trainPres[trainPres$origUsable, ]
	
	# # remove erroneous point in desert
	# trainPres <- trainPres[-which(trainPres$dataName == 'Connie Millar - Subset A - Pika Present' & trainPres$origLocality == 'Sierra Nevada; Humphreys Basin; Goethe rock glacier; at rg outlet stream; Fresno; CA' & trainPres$origElev == 3523 & trainPres$origNotes != 'RGC; granitic;' & trainPres$origDateOfObs == '8/16/2011'), ]
	
	# testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# testSurveys <- as.data.frame(testSurveys)
	
	# pcs <- c('pc1', 'pc2')
	# cex <- 0.4

	# ### loadings
	# loads <- pca$loadings
	# loads <- as.data.frame(loads)
	# loads <- loads[complete.cases(loads), ]
	# loads <- as.data.frame(loads)
	# names(loads) <- paste0('PC', 1:9)
	# rownames(loads) <- c('eastness', 'northness', 'acute cold', 'chronic heat', 'water deficit', 'solar radiation', 'summer respite', 'SWE', 'NDVI')
	
	# xlab <- paste0('PC1 (', round(100 * pca$sdev[1]^2 / sum(pca$sdev^2), 1), '%)')
	# ylab <- paste0('PC2 (', round(100 * pca$sdev[2]^2 / sum(pca$sdev^2), 1), '%)')
	# png('./Figures & Tables/PCA with Test Classes.png', width=1800, height=1800, res=600)

		# par(cex.axis=0.5, cex.lab=0.6, mar=rep(1.8, 4), oma=rep(0, 4), lwd=0.6, bty='n', mgp=c(0.7, -0, 0), tck=-0.01)
		
		# smoothScatter(pca$scores[ , 1:2], pch=16, nrpoints=0, xlab=xlab, ylab=ylab)
		# points(trainPres[ , pcs], pch=3, col=alpha('black', 0.4), cex=0.5 * cex)
		
		# pres <- testSurveys[testSurveys$status == '2 detected', ]
		# recentAbs <- testSurveys[testSurveys$status == '1 recent absence', ]
		# longTermAbs <- testSurveys[testSurveys$status == '0 long absence', ]
		
		# points(longTermAbs[ , pcs], pch=25, bg='red', cex=cex)
		# points(recentAbs[ , pcs], pch=23, bg='yellow', cex=cex)
		# points(pres[ , pcs], pch=21, bg='green', cex=cex)

		# mult <- 4
		# x0 <- 10.3
		# y0 <- -3
		# for (i in 1:9) {
		
			# arrows(x0, y0, x0 + loads[i, 1] * mult, y0 + loads[i, 2] * mult, angle=20, length=0.07, xpd=NA)
			
			# label <- rownames(loads)[i]
			# if (label == 'summer respite') {
				# ydelta <- + 0.2
			# } else {
				# ydelta <- 0
			# }
			
			# text(x0 + loads[i, 1] * mult, y0 + loads[i, 2] * mult + ydelta, labels = label, cex=0.5, xpd=NA)
		# }

		# legend('topleft', inset=c(0.01, 0.1), legend=c('Background', 'Previously-known occurrence', 'Currently occupied', 'Formerly occupied', 'No evidence'), pch=c(NA, 3, 21, 23, 25), fill=c(blues9[6], NA, NA, NA, NA), col=c(NA, 'black', 'black', 'black', 'black'), pt.bg = c(NA, NA, 'green', 'yellow', 'red'), border=c('black', NA, NA, NA, NA), bty='n', bg=NA, cex=cex, xpd=NA)


	# dev.off()	

# say('##############################################################################')
# say('### visual comparison of background, TRAINING, and TEST sites by predictor ###')
# say('##############################################################################')	
	
	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.rda')
	# pcs <- paste0('pc', whichPcs(pca))

	# load('./ENMs/Background Sites/Random Background Sites from across Study Region Mask.rda')
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.rda')
	# trainPres$origUsable <- as.logical(trainPres$origUsable)
	# trainPres <- trainPres[trainPres$origUsable, ]

	# testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# testSurveys <- as.data.frame(testSurveys)
	
	# for (pc in pcs) {
		
		# png(paste0('./Figures & Tables/', toupper(pc), '.png'), width=1000, height=800, res=300)
			
			# par(mgp=c(1.6, 0.4, 0), mar=c(3, 4, 2, 1) + 0.1, cex.main=0.8, cex.lab=0.6, cex.axis=0.6, tck=-0.02)
			
			# thisBg <- bg[ , pc]
			# thisTrainPres <- trainPres[ , pc]
			# longTermAbs <- testSurveys[testSurveys$status=='0 long absence', pc]
			# recentAbs <- testSurveys[testSurveys$status=='1 recent absence', pc]
			# detected <- testSurveys[testSurveys$status=='2 detected', pc]
			
			# thisBg <- sort(thisBg)
			# thisTrainPres <- sort(thisTrainPres)
			# longTermAbs <- sort(longTermAbs)
			# recentAbs <- sort(recentAbs)
			# detected <- sort(detected)
			
			# bgCs <- cumsum(seq_along(thisBg)) / sum(seq_along(thisBg))
			# trainPresCs <- cumsum(seq_along(thisTrainPres)) / sum(seq_along(thisTrainPres))
			# longTermAbsCs <- cumsum(seq_along(longTermAbs)) / sum(seq_along(longTermAbs))
			# recentAbsCs <- cumsum(seq_along(recentAbs)) / sum(seq_along(recentAbs))
			# detectedCs <- cumsum(seq_along(detected)) / sum(seq_along(detected))
		
			# allVals <- c(thisBg, thisTrainPres, longTermAbs, recentAbs, detected)
			# xlim <- range(allVals)
		
			# plot(thisBg, bgCs, xlab=toupper(pc), ylab='Cumulative proportion', type='l', col='blue', lwd=2, main=toupper(pc), xlim=xlim)
			# lines(thisTrainPres, trainPresCs, col='black', lwd=2)
			# lines(longTermAbs, longTermAbsCs, col='red', lwd=2)
			# lines(recentAbs, recentAbsCs, col='yellow', lwd=2)
			# lines(detected, detectedCs, col='chartreuse', lwd=2)
			
			# legend('topleft', inset=0.01, legend=c('Background', 'Training presences', 'Test presences', 'Recent absences', 'Long-term absences'), col=c('blue', 'black', 'chartreuse', 'yellow', 'red'), lwd=2, bty='n', cex=0.4)
		
		# dev.off()
	
	# }
	
# say('##########################')
# say('### map of predictions ###')
# say('##########################')

	# # generalize
	# focusBuff <- 25 # buffer size around test sites for generating focus of plot (in km)

	# nam <- vect('C:/Ecology/Research Data/GADM/Version 4.1/High Res North America Level 1 sans Great Lakes SpatVector WGS84.gpkg')
	# nam <- as(nam, 'Spatial')
	# # usa <- raster::getData('GADM', country='USA', level=1, path='C:/!Scratch')
	# # mex <- raster::getData('GADM', country='MEX', level=1, path='C:/!Scratch')
	# # nam <- rbind(usa, mex)
	# nam <- sp::spTransform(nam, getCRS('climateNA', TRUE))

	# ### spatial data
	# load('./Study Region/GADM California, Nevada, Oregon States.rda')
	# load('./Study Region/GADM California, Nevada, Oregon Counties.rda')
	# load('./Study Region/GADM Plumas County.rda')
	# load(file='./Study Region/GADM GAP Counties.rda')
	# hs <- raster('./Data/Topography - SRTM/hillshade_srtm_ea.tif')

	# gapCounties <- sp::spTransform(gapCounties, getCRS('climateNA', TRUE))
	# plumas <- sp::spTransform(plumas, getCRS('climateNA', TRUE))
	# west1 <- sp::spTransform(west1, getCRS('climateNA', TRUE))
	# west2 <- sp::spTransform(west2, getCRS('climateNA', TRUE))

	# ### training occurrences
	# load('./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.rda')
	# surveys$origUsable <- as.logical(surveys$origUsable)
	# surveys <- surveys[surveys$origUsable, ]	
	# surveys <- surveys[surveys$origRecentPikaOrSignsObserved, ]
	# trainPresSp <- SpatialPointsDataFrame(surveys[ , ll], data=surveys, proj4=getCRS('wgs84', TRUE))
	# trainPresSpEa <- sp::spTransform(trainPresSp, getCRS('climateNA', TRUE))

	# ### test occurrences
	# testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# testSurveys <- as.data.frame(testSurveys)
	# testSurveysSp <- SpatialPointsDataFrame(testSurveys[ , ll], data=testSurveys, proj4=getCRS('nad83', TRUE))
	# testSurveysSpEa <- sp::spTransform(testSurveysSp, getCRS('climateNA', TRUE))
	
	# ### plot focus
	# focus <- gBuffer(testSurveysSpEa, width=1000 * focusBuff)

	# ### crop hillshade
	# plot(focus)
	# usr <- par('usr')
	# dev.off()

	# ext <- extent(usr)
	# ext <- as(ext, 'SpatialPolygons')
	# projection(ext) <- getCRS('climateNA')

	# hs <- crop(hs, ext)
	# hs <- hs - minValue(hs)
	# hs <- hs / maxValue(hs)
		
	# ### ENM prediction
	# pred <- raster('./ENMs/predictionEnsemble.tif')
	# predVals <- extract(pred, testSurveysSp)
	# beginCluster(4)
		# pred <- projectRaster(pred, crs=getCRS('climateNA'))
	# endCluster()

	# ### lakes
	# lakes <- shapefile(paste0(drive, '/Research Data/North America Rivers and Lakes USGS/hydrography_p_lakes_v2.shp'))
	# for (i in 1:20) {
		# maxs <- which.max(lakes$Shape_Leng)
		# lakes <- lakes[-maxs, ]
	# }
	# lakes <- sp::spTransform(lakes, getCRS('climateNA', TRUE))

	
	# ### colors for points
	# trainPresCol <- 'black'
	# testPresCol <- 'darkgreen'
	# testPresFill <- 'chartreuse3'
	# testShortTermAbsCol <- 'darkorange4'
	# testShortTermAbsFill <- 'darkgoldenrod3'
	# testLongTermAbsCol <- 'darkred'
	# incorrectPredCol <- 'cyan'
		
	# ### colors for rasters
	# predPresCol <- alpha('forestgreen', 0.6)
	# predShortTermAbsCol <- alpha('darkgoldenrod3', 0.6)
	# predLongTermAbsCol <- alpha('red', 0.4)
	
	# trainPresPch <- 3

	# testLongTermPch <- 25
	# testShortTermPch <- 23
	# testPresPch <- 21

	# predLongTerm <- predVals[testSurveys$status == '0 long absence'] / 1000
	# predShortTerm <- predVals[testSurveys$status == '1 recent absence'] / 1000
	# predPres <- predVals[testSurveys$status == '2 detected'] / 1000


	# tholdShortVsLong <- thresholdWeighted(predShortTerm, predLongTerm)[['mdss']]
	# tholsShortVsPres <- thresholdWeighted(predPres, predShortTerm)[['mdss']]

	# ### plot
	# png('./Figures & Tables/ENM Prediction.png', width=2 * 1200, height=2 * 1200, res=300)
		
		# par(mfrow=c(2, 2), mar=c(0.1, 0.1, 0.1, 0.1), oma=c(3, 3, 0.1, 0.1), mgp=c(3, 0.4, 0))

		# ### each class by itself
		
		# for (class in c(0:2)) {
		
			# plot(focus, border='white')
			# usr <- par('usr')
			
			# if (class == 0) {

				# predClass <- (pred / 1000) < tholdShortVsLong
				# longTermRast <- predClass
				# pointPreds <- sort(predLongTerm)

				# col <- predLongTermAbsCol
				# pointOutlineCol <- testLongTermAbsCol
				# pointFillCol <- rep(NA, length(pointPreds))
				# pointFillCol[pointPreds >= tholdShortVsLong & pointPreds < tholsShortVsPres] <- testShortTermAbsFill
				# pointFillCol[pointPreds >= tholsShortVsPres] <- testPresFill
				# rastLab <- 'No evidence predicted'
				# pointLab <- c('Correctly predicted', 'Former occ. predicted', 'Current occ. predicted')
				# pch <- testLongTermPch
				# axis(2, at=pretty(c(usr[3:4]), 3), tck=-0.01)
				
				# legColBg <- c(NA, NA, testShortTermAbsFill, testPresFill)
				
			# } else if (class == 1) {

				# predClass <- (pred / 1000) >= tholdShortVsLong & (pred / 1000) < tholsShortVsPres
				# shortTermRast <- predClass
				# pointPreds <- sort(predShortTerm)

				# col <- predShortTermAbsCol
				# pointOutlineCol <- testShortTermAbsCol
				# pointFillCol <- rep(NA, length(pointPreds))
				# pointFillCol[pointPreds < tholdShortVsLong] <- testLongTermAbsCol
				# pointFillCol[pointPreds < tholdShortVsLong & pointPreds < tholsShortVsPres] <- NA
				# pointFillCol[pointPreds >= tholsShortVsPres] <- testPresFill
				# rastLab <- 'Old evidence predicted'
				# pointLab <- c('No evidence predicted', 'Correctly predicted', 'Current occ. predicted')
				# pch <- testShortTermPch
				
				# legColBg <- c(NA, testLongTermAbsCol, NA, testPresFill)
				
			# } else if (class == 2) {
			
				# predClass <- (pred / 1000) >= tholsShortVsPres
				# presRast <- predClass
				# pointPreds <- sort(predPres, FALSE)

				# col <- predPresCol
				# pointOutlineCol <- testPresCol
				# pointFillCol <- rep(NA, length(pointPreds))
				# pointFillCol[pointPreds < tholdShortVsLong] <- testLongTermAbsCol
				# pointFillCol[pointPreds >= tholdShortVsLong & pointPreds < tholsShortVsPres] <- testShortTermAbsFill
				# pointFillCol[pointPreds >= tholsShortVsPres] <- NA
				# rastLab <- 'Presence predicted'
				# pointLab <- c('No evidence predicted', 'Former occ. predicted', 'Correctly predicted')
				# pch <- testPresPch
				# axis(1, at=pretty(c(usr[1:2]), 3)[1:3], tck=-0.01)
				# axis(2, at=pretty(c(usr[3:4]), 3), tck=-0.01)
				
				# legColBg <- c(NA, testLongTermAbsCol, testShortTermAbsFill, NA)
				
			# }
			
			# plot(hs, col=grays, legend=FALSE, add=TRUE)
			# plot(lakes, col = 'gray60', border = NA, add = TRUE)
			# plot(predClass, col=c(NA, col), legend=FALSE, add=TRUE)
			# plot(west2, border='gray40', add=TRUE)
			# plot(west1, lwd=2, add=TRUE)
			# points(trainPresSpEa, pch=3, cex=0.8, col=trainPresCol)

			# pchs <- rep(NA, nrow(testSurveys))
			
			# pchs[testSurveysSpEa$status == '0 long absence'] <- testLongTermPch
			# pchs[testSurveysSpEa$status == '1 recent absence'] <- testShortTermPch
			# pchs[testSurveysSpEa$status == '2 detected'] <- testPresPch



			# if (class == 0) {
				# points(testSurveysSpEa[testSurveysSpEa$status == '0 long absence', ], bg=pointFillCol, col=pointOutlineCol, pch=pch, cex=1.2)
			# } else if (class == 1) {
				# points(testSurveysSpEa[testSurveysSpEa$status == '1 recent absence', ], bg=pointFillCol, col=pointOutlineCol, pch=pch, cex=1.2)
			# } else if (class == 2) {
				# points(testSurveysSpEa[testSurveysSpEa$status == '2 detected', ], bg=pointFillCol, col=pointOutlineCol, pch=pch, cex=1.5)
			# }
				
			# box()
			
			# legend(
				# 'bottomleft',
				# inset=0.01,
				# legend=c(rastLab, pointLab),
				# pch=c(NA, pch, pch, pch),
				# fill=c(col, NA, NA, NA),
				# border=c('black', NA, NA, NA),
				# col=c(NA, pointOutlineCol, pointOutlineCol, pointOutlineCol),
				# pt.bg=legColBg,
				# bg=alpha('white', 0.9),
				# cex=0.85,
				# pt.cex=1.2
			# )
			
			# # scale bar
			# size <- 50000 # length of scale bar in meters
			# width <- usr[2] - usr[1]
			# height <- usr[4] - usr[3]
			# x <- usr[2] - size - 0.02 * width
			# x <- c(x, x + size)
			# y <- usr[3] + rep(0.02 * height, 2)
			
			# lines(x, y, lwd=4, xpd=NA, col='black', lend=1)
			
			# x <- mean(x)
			# y <- usr[3] + 0.05 * height
			# text(x, y[1], labels=paste(size / 1000, 'km'), cex=1)
			
		# }

		# ### all together

		# plot(focus, border='white')
		# axis(1, at=pretty(c(usr[1:2]), 3)[1:3], tck=-0.01)
		# plot(hs, col=alpha(grays, 0.8), legend=FALSE, add=TRUE)
		# plot(lakes, col = 'gray60', border = NA, add = TRUE)
		# plot(longTermRast, col=c(NA, predLongTermAbsCol), legend=FALSE, add=TRUE)
		# plot(shortTermRast, col=c(NA, predShortTermAbsCol), legend=FALSE, add=TRUE)
		# plot(presRast, col=c(NA, predPresCol), legend=FALSE, add=TRUE)
		# plot(west2, border='gray40', add=TRUE)
		# plot(west1, lwd=2, add=TRUE)
		# points(trainPresSpEa, pch=trainPresPch, cex=0.8)

		# testSurveysSpEa <- testSurveysSpEa[order(testSurveysSpEa$status), ]
		
		# cols <- rep(NA, nrow(testSurveys))
		# cols[testSurveysSpEa$status == '0 long absence'] <- testLongTermAbsCol
		# cols[testSurveysSpEa$status == '1 recent absence'] <- testShortTermAbsCol
		# cols[testSurveysSpEa$status == '2 detected'] <- testPresCol
		
		# cexs <- pchs <- rep(NA, nrow(testSurveys))
		
		# pchs[testSurveysSpEa$status == '0 long absence'] <- testLongTermPch
		# pchs[testSurveysSpEa$status == '1 recent absence'] <- testShortTermPch
		# pchs[testSurveysSpEa$status == '2 detected'] <- testPresPch
		
		# cexs[testSurveysSpEa$status == '0 long absence'] <- 1.2
		# cexs[testSurveysSpEa$status == '1 recent absence'] <- 1.2
		# cexs[testSurveysSpEa$status == '2 detected'] <- 1.5

		# points(testSurveysSpEa, col=cols, pch=pchs, cex=cexs)
			
		# box()
		
		# legend(
			# 'bottomleft',
			# inset=0.01,
			# legend=c('No evidence predicted', 'Former occ. predicted', 'Current occ. predicted', 'No evidence observed', 'Former occ. observed', 'Current occ. observed', 'Previously-known occ.'),
			# pch=c(NA, NA, NA, testLongTermPch, testShortTermPch, testPresPch, trainPresPch),
			# fill=c(predLongTermAbsCol, predShortTermAbsCol, predPresCol, NA, NA, NA, NA),
			# border=c('black', 'black', 'black', NA, NA, NA, NA),
			# col=c(NA, NA, NA, testLongTermAbsCol, testShortTermAbsCol, testPresCol, trainPresCol),
			# bg=alpha('white', 0.4),
			# cex=0.8,
			# pt.cex=1.2
		# )
		
			# # scale bar
			# size <- 50000 # length of scale bar in meters
			# width <- usr[2] - usr[1]
			# height <- usr[4] - usr[3]
			# x <- usr[2] - size - 0.02 * size
			# x <- c(x, x + size)
			# y <- usr[3] + rep(0.02 * height, 2)
			
			# lines(x, y, lwd=4, xpd=NA, col='black', lend=1)
			
			# x <- mean(x)
			# y <- usr[3] + 0.05 * height
			# text(x, y[1], labels=paste(size / 1000, 'km'), cex=1)

		# ### inset map
		
			# usr <- par('usr')
			# par(fig=c(0.22, 0.64, 0.18, 0.54), bg='white', new=TRUE)
			
			# # get extent of inset
			# largeFocus <- nam[nam@data$NAME_1 %in% c('California', 'Nevada', 'Oregon', 'Baja California'), ]
			# largeFocus <- gBuffer(largeFocus, width=50000)
			# insetExt <- extent(largeFocus)
			# insetExt <- as(insetExt, 'SpatialPolygons')
			# projection(insetExt) <- getCRS('climateNA')
			# insetExt <- vect(insetExt)
			# namVect <- vect(nam)
			# namCrop <- terra::crop(namVect, insetExt)
			
			# # get polygon to highlight focal area
			# focusExt <- extent(usr)
			# focusExt <- as(focusExt, 'SpatialPolygons')
			# projection(focusExt) <- getCRS('climateNA')

			# plot(insetExt, border='gray', col='white', axes=FALSE)
			# plot(namCrop, col='gray', border='gray40', ann=FALSE, add=TRUE)
			# plot(focusExt, lwd=2.2, border='black', add=TRUE)
			# plot(insetExt, lwd=1.2, axes=FALSE, add=TRUE)


		# title(sub=date(), cex.sub=0.3, outer=TRUE, line=2)
			
	# dev.off()
	
# say('#######################################')
# say('### confusion matrix of predictions ###')
# say('#######################################')

	# ### test occurrences
	# testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# testSurveys <- as.data.frame(testSurveys)
	# testSurveysSp <- SpatialPointsDataFrame(testSurveys[ , ll], data=testSurveys, proj4=getCRS('nad83', TRUE))
	# testSurveysSpEa <- sp::spTransform(testSurveysSp, getCRS('climateNA', TRUE))
	
	# ### ENM prediction
	# pred <- raster('./ENMs/predictionEnsemble.tif')
	# predVals <- extract(pred, testSurveysSp)
	
	# predPres <- predVals[testSurveys$status == '2 detected'] / 1000
	# predShortTerm <- predVals[testSurveys$status == '1 recent absence'] / 1000
	# predLongTerm <- predVals[testSurveys$status == '0 long absence'] / 1000

	# confuse <- data.frame(
		# truePres=rep(NA, 4),
		# trueRecentAbs=rep(NA, 4),
		# trueLongAbs=rep(NA, 4),
		# sum=rep(NA, 4),
		# row.names=c('predictedPres', 'predictedRecentAbs', 'predictedLongAbs', 'sum')
	# )
	
	# tholdPres_vs_shortTermAbs <- thresholdWeighted(predPres, predShortTerm)[['mdss']]
	# tholdPres_vs_longTermAbs <- thresholdWeighted(predPres, predLongTerm)[['mdss']]
	# tholdShortTerm_vs_longTermAbs <- thresholdWeighted(predShortTerm, predLongTerm)[['mdss']]
	
	# confuse$truePres <- c(
		# sum(predPres >= tholdPres_vs_shortTermAbs),
		# sum(predPres < tholdPres_vs_shortTermAbs & predPres >= tholdPres_vs_longTermAbs),
		# sum(predPres < tholdPres_vs_longTermAbs),
		# length(predPres)
	# )
	
	# confuse$trueRecentAbs <- c(
		# sum(predShortTerm >= tholdPres_vs_shortTermAbs),
		# sum(predShortTerm < tholdPres_vs_shortTermAbs & predShortTerm >= tholdShortTerm_vs_longTermAbs),
		# sum(predShortTerm < tholdShortTerm_vs_longTermAbs),
		# length(predShortTerm)
	# )
	
	# confuse$trueLongAbs <- c(
		# sum(predLongTerm >= tholdPres_vs_longTermAbs),
		# sum(predLongTerm < tholdPres_vs_longTermAbs & predLongTerm >= tholdShortTerm_vs_longTermAbs),
		# sum(predLongTerm < tholdShortTerm_vs_longTermAbs),
		# length(predLongTerm)
	# )
	
	# confuse$sum <- rowSums(confuse[ , 1:3])
	
	# # in/correct classification rate
	# props <- confuse[1:3, 1:3]
	# props <- props / length(predVals)
	
	# ccr <- round(sum(diag(as.matrix(props))), 3)
	# icr <- round(sum(props) - sum(diag(as.matrix(props))), 3)
	
	# # classification rates by class
	# propsByCol <- confuse[1:3, 1:3]
	# for (i in 1:3) {
		# propsByCol[i, ] <- propsByCol[i, ] / confuse[4, 1:3]
	# }
	
	# propsByCol <- round(propsByCol, 3)
	
	# sink('./Figures & Tables/Confusion Matrix - MDSS Threshold.txt', split=TRUE)
		# print(confuse)
		# say('')
		# print(propsByCol)
		# say('')
		# say('Total correct classification rate: ', ccr)
		# say('Total incorrect classification rate: ', icr)
		# say('')
		
		# say('Threshold minimizing difference between correct predictions of presences and short-term absences: ', tholdPres_vs_shortTermAbs)
		# say('Threshold minimizing difference between correct predictions of presences and long-term absences: ', tholdPres_vs_longTermAbs)
		# say('Threshold minimizing difference between correct predictions of short- and long-term absences: ', tholdShortTerm_vs_longTermAbs)
		
		# say(date())
	# sink()

# # # say('###############################################')
# # # say('### make KML of test sites with predictions ###')
# # # say('###############################################')
	
	# # # tests <- fread('./Data/Occurrences/Test Surveys 02 Cleaned.csv')
	# # # tests <- as.data.frame(tests)

	# # # tests <- SpatialPointsDataFrame(tests[ , c('longWgs84', 'latWgs84')], data=tests, proj4string=getCRS('wgs84', TRUE))
	# # # pred <- raster('./ENMs/predictionEnsemble.tif')
	
	# # # tests$prediction <- extract(pred, tests) / 1000
	
	# # # tests$predictedStatus <- NA
	# # # tests$predictedStatus[tests$prediction < 0.59] <- '0 long absence'
	# # # tests$predictedStatus[tests$prediction >= 0.59 & tests$prediction < 0.74] <- '1 recent absence'
	# # # tests$predictedStatus[tests$prediction >= 0.74] <- '2 detection'

	# # # classFx <- function(x) {
		# # # if (is.na(x)) {
			# # # NA
		# # # } else if (x < 0.59) {
			# # # 0
		# # # } else if (x >= 0.59 & x < 0.74) {
			# # # 1
		# # # } else if (x >= 0.74) {
			# # # 2
		# # # }
	# # # }
	
	# # # tholdPred <- calc(pred / 1000, classFx)
	
	# # # KML(tests, './Figures & Tables/California Gap Test Sites', overwrite=TRUE)
	# # # KML(pred, './Figures & Tables/California Gap Model Predictions', mappixel=10 * 100000, overwrite=TRUE)
	# # # KML(tholdPred, './Figures & Tables/California Gap Model Predictions Thresholded', mappixel=100 * 100000, overwrite=TRUE)
	
# say('############################################')
# say('### calculate climate means for gap area ###')
# say('############################################')

	# swe <- rast(listFiles('./Data/Climate - Monthly Means/swe', pattern='.tif'))
	# maxs <- rast(listFiles('./Data/Climate - Monthly Means/tmax', pattern='.tif'))
	# mins <- rast(paste0('./Data/Climate - Monthly Means/tmin/tmin_month', prefix(c(12, 1, 2), 2), '_mean2010to2019.tif'))

	# load('./Study Region/GADM Plumas County.rda')
	# plumas <- vect(plumas)
	
	# maxs <- mean(maxs)
	# mins <- mean(mins)
	# swe <- sum(swe)
	
	# maxs <- extract(maxs, plumas, ID=FALSE)
	# mins <- extract(mins, plumas, ID=FALSE)
	# swe <- extract(swe, plumas, ID=FALSE)

	# maxs <- colMeans(maxs)
	# mins <- colMeans(mins)
	# swe <- colMeans(swe)

# say('##################################################################')
# say('### graph of elevation of calibration and survey sites by type ###')
# say('##################################################################')

	# # Make a violin plot of occurrences/test sites where y-axis is elevation and x is category (train, test no evidence, test recent evidence, test occupied)

	# # prismElev <- rast('./Data/Topography - SRTM/elevation_srtm_m.tif')
	# prismElev <- rast(paste0(drive, '/Research Data/PRISM/PRISM_us_dem_800m.tif'))

	# # records
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.rda')
	# trainPres$origUsable <- as.logical(trainPres$origUsable)
	# trainPres <- trainPres[trainPres$origUsable, ]
	# test <- fread('./Data/Occurrences/Test Surveys 04 Extracted Environmental Values with PCA.csv')

	# # convert to spatial
	# test <- vect(test, geom = c('longWgs84', 'latWgs84'), crs = getCRS('wgs84'))
	# trainPres <- vect(trainPres, geom = c('longWgs84', 'latWgs84'), crs = getCRS('wgs84'))

	# # remove training presences far from test presenecs (to accomodate latitudinal shift in minimum pika elevation)
	# buff <- buffer(test, 100 * 1000)
	# extent <- ext(buff)
	# extent <- vect(extent, crs = getCRS('wgs84'))

	# inBox <- extract(extent, trainPres)
	# trainPres <- trainPres[!is.na(inBox[ , 'id.x']), ]

	# testElev <- extract(prismElev, test, ID = FALSE)
	# trainElev <- extract(prismElev, trainPres, ID = FALSE)

	# testElev$status <- test$status
	# testElev$status[testElev$status == '0 long absence'] <- 'No evidence'
	# testElev$status[testElev$status == '1 recent absence'] <- 'Formerly occupied'
	# testElev$status[testElev$status == '2 detected'] <- 'Currently occupied'

	# trainElev$status <- 'Previously-known presence'

	# elevs <- rbind(trainElev, testElev)
	
	# elevs$status <- factor(elevs$status, levels = c('Previously-known presence', 'No evidence', 'Formerly occupied', 'Currently occupied'), ordered = TRUE)
	
	# graph <- ggplot(elevs, aes(x = status, y = PRISM_us_dem_800m, fill = status)) +
	
		# geom_violin() +
		# scale_x_discrete(
			# labels = c('Previously-known presence' = 'Previously-known\npresence','No evidence' = 'No\nevidence', 'Formerly occupied' = 'Formerly\noccupied', 'Currently occupied' = 'Currently\noccupied')
		# ) +
		# scale_fill_manual(
			# values = c('Previously-known presence' = 'black', 'No evidence' = 'darkred', 'Formerly occupied' = 'yellow', 'Currently occupied' = 'darkgreen')
		# ) + 
		# xlab('') + ylab('Elevation (m)') +
		# theme(
			# legend.position = 'none',
			# axis.text = element_text(size = 18),
			# axis.title = element_text(size = 22)
		# )
	
	# ggsave(graph, file = './Figures & Tables/Elevation by Sites.png', dpi = 600, width = 8, height = 8)

# say('###################################')
# say('### ancillary analyses of sites ###')
# say('###################################')

	# ### inter-point distances between survey occurrences in the extreme NW corner of survey sites
	# #############################################################################################

	# sink('./Figures & Tables/Ancillary Analysis - Distances between Clusters of Test Detections.txt', split = TRUE)
	# say(date())

	# # northwestern test presences
	# testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# testPres <- testSurveys[status == '2 detected' & latWgs84 > 40]
	# testPres <- vect(testPres, geom = ll, crs = prismCrs)
	
	# dists <- distance(testPres)
	
	# say('Number of survey detections in northwest part of gap: ', nrow(testPres), pre = 2)
	# say('Pairwise distance(s) between them in meters:')
	# print(round(dists))

	# # center/north of southern edge of gap test presences
	# testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# testPres <- testSurveys[status == '2 detected' & latWgs84 > 39.6 & latWgs84 < 39.8]
	# testPres <- vect(testPres, geom = ll, crs = prismCrs)
	
	# dists <- distance(testPres)
	
	# say('Number of survey detections in northern part of southern gap edge region: ', nrow(testPres), pre = 2)
	# say('Pairwise distance(s) between them in meters:')
	# print(round(dists))

	# # southeastern-most test presences
	# testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# testPres <- testSurveys[status == '2 detected' & latWgs84 < 39.6]
	# testPres <- vect(testPres, geom = ll, crs = prismCrs)
	
	# dists <- distance(testPres)
	
	# say('Number of survey detections in southeastern-most part of gap: ', nrow(testPres), pre = 2)
	# say('Pairwise distance(s) between them in meters:')
	# print(round(dists))

	# sink()

	# # # ### shortest distance between training occurrences and recent test sites in almost southeastern-most locale
	# # # ###########################################################################################################
	
	# # # load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.rda')
	# # # trainPres <- vect(trainPres, geom = ll, crs = prismCrs)
	
	# # # testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# # # testRecentAbs <- testSurveys[status == '1 recent absence' & latWgs84 > 39.55 & latWgs84 < 39.62 & longWgs84 > -120.2]
	# # # testRecentAbs <- vect(testRecentAbs, geom = ll, crs = prismCrs)

	# # # plot(testRecentAbs)
	# # # plot(trainPres, pch = 1, col = 'red', add = TRUE)
		
	# # # dists_m <- distance(testRecentAbs, trainPres)
	# # # nearestDist_m <- apply(dists_m, 2, min)
	# # # nearby <- which(nearestDist_m <= 3000)
	
	# # # trainPresNearby <- trainPres[nearby]
	# # # dists_m <- distance(testRecentAbs, trainPresNearby)
	
	# # # sink('./Figures & Tables/Ancillary Analysis - Distances between Southwestern-most Recent Absences and Training Presences.txt', split = TRUE)
	# # # say(date(), post = 2)
	# # # say('Number of recent absences: ', nrow(testRecentAbs))
	# # # say('Number of nearby previously-known presences within 3000 m: ', nrow(trainPresNearby), post = 2)
	
	# # # say('Distance(s) between southwestern-most recent absence(s) and training occurrence(s) (in meters):', post = 2)
	# # # print(round(dists_m))
	# # # say('')
	# # # say('Training occurrence closest to this recent absence site:')
	# # # print(as.data.frame(trainPres[nearby]))
	
	# # # sink()
	
	# ### analysis of pellet elevations
	# #################################
	
	# buried <- c(7822, 7500, 7259, 7008, 6978, 6870, 6229, 7546)
	# nonburied <- c(9400, 8255, 7556, 7449, 7297, 7192, 7008, 6535, 6705, 7418, 7546)

	# pellets <- data.frame(elev_m = c(buried, nonburied), status = c(rep('buried', length(buried)), rep('non-buried', length(nonburied))))
	# boxplot(elev_m ~ status, data = pellets)

	# ### straight-line distance between only presences found in gap and nearest southern presence and nearest northern presence
	# ##########################################################################################################################

	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.rda')
	# trainPres <- vect(trainPres, geom = ll, crs = getCRS('WGS84'))

	# load('./Study Region/GADM GAP Counties.rda')
	# gapCounties <- vect(gapCounties)
	# crs(gapCounties) <- getCRS('WGS84')

	# testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# testPres <- testSurveys[status == '2 detected']
	# testPres <- vect(testPres, geom = ll, crs = getCRS('WGS84'))

	# plot(gapCounties)
	# plot(trainPres, pch = 16, cex = 0.6, col = 'red', add = TRUE)
	# plot(testPres, pch = 1, cex = 0.8, , col = 'blue', add = TRUE)

	# pts <- click(n = 3)
	# pts <- vect(pts, crs = getCRS('WGS84'))
	# dists_m <- distance(pts)
	
# say('########################')
# say('### measure gap area ###')
# say('########################')

	# # Measuring "gap" area using area of triangles of Delauney triangulation
	
	# sn <- vect('./Study Region/epa3_sierraNevadaSarrModified.shp')
	# sn <- project(sn, prismCrs)

	# ### area of gap prior to sampling (model training occurrences)
	# ##############################################################
	# load('./Study Region/GADM GAP Counties.rda')
	
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.rda')
	# trainPres <- vect(trainPres, geom = c('longWgs84', 'latWgs84'), crs = prismCrs)
	
	# del <- delaunay(trainPres)
	
	# plot(gapCounties, border = 'red')
	# plot(trainPres, add = TRUE)
	# plot(del, add = TRUE)
	# delArea_m2 <- expanse(del) / 1000^2
	
	# # manually select triangles from Delauney triangulation spanning gap
	# # pts <- click(n = 4) # how I got the points below
	# pts <- matrix(c(
		# -120.4289, 40.37693,
		# -120.9500, 40.33288,
		# -121.2298, 40.22275,
		# -121.2395, 40.36225
	# ), ncol = 2, byrow = TRUE)

	# gapDel <- extract(del, pts)
	# gapDel <- as.data.frame(gapDel)
	# gap <- del[gapDel$id.x, ]
	# plot(sn)
	# plot(gap, border=3, add=TRUE)
	
	# gap <- crds(gap)
	
	# # extend gap east and west
	# gapTop <- gap[gap[ , 2] > 40, ]
	# gapBottom <- gap[gap[ , 2] < 40, ]
	
	# extends <- matrix(c(
		# min(gapTop[ , 1]) - 1, gapTop[which.min(gapTop[ , 1]), 2],
		# max(gapTop[ , 1]) + 1, gapTop[which.max(gapTop[ , 1]), 2],
		# min(gapBottom[ , 1]) - 2, gapBottom[which.min(gapBottom[ , 1]), 2],
		# max(gapBottom[ , 1]) + 1, gapBottom[which.max(gapBottom[ , 1]), 2]
	# ), ncol = 2, byrow = TRUE)
	
	# gap <- rbind(gap, extends)
	
	# gapTop <- gap[gap[ , 2] > 40, ]
	# gapBottom <- gap[gap[ , 2] < 40, ]
	# gapTop <- gapTop[order(gapTop[ , 1]), ]
	# gapBottom <- gapBottom[order(gapBottom[ , 1], decreasing = TRUE), ]
	
	# gap <- rbind(gapTop, gapBottom)
	# gap <- vect(gap, type = 'polygons', crs = prismCrs)
	
	# plot(sn)
	# plot(gap, border=4, add=TRUE)

	# gap <- intersect(sn, gap)
	# gapArea_m2 <- expanse(gap) / 1000^2
	# snArea_m2 <- expanse(sn) / 1000^2
	
	# ### area of "recently-occupied" gap (area of "recent occupied" area predicted by model within gap Delauney triangles)
	# #####################################################################################################################
	
	# testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# testSurveys <- vect(testSurveys, geom = ll, crs = prismCrs)

	# ### ENM predictions in the "yellow" zone
	# pred <- rast('./ENMs/predictionEnsemble.tif')
	# pred <- pred / 1000
	# predVals <- extract(pred, testSurveys, ID = FALSE)
	# predVals <- unlist(predVals)
	
	# predPres <- predVals[testSurveys$status == '2 detected']
	# predShortTerm <- predVals[testSurveys$status == '1 recent absence']
	# predLongTerm <- predVals[testSurveys$status == '0 long absence']

	# tholdPres_vs_shortTermAbs <- enmSdmX::evalThreshold(predPres, predShortTerm)[['mdss']]
	# tholdPres_vs_longTermAbs <- enmSdmX::evalThreshold(predPres, predLongTerm)[['mdss']]
	# tholdShortTerm_vs_longTermAbs <- enmSdmX::evalThreshold(predShortTerm, predLongTerm)[['mdss']]
	
	# predRecentOcc <- pred > tholdShortTerm_vs_longTermAbs & pred < tholdPres_vs_shortTermAbs
		
	# cellAreas_m2 <- cellSize(pred) / 1000^2
	# cellAreas_m2 <- predRecentOcc * cellAreas_m2
	
	# plot(cellAreas_m2)
	# plot(sn, add = TRUE)
	# plot(gap, border = 'blue', add = TRUE)

	# recentAbsAreaInGap_m2 <- extract(cellAreas_m2, gap, ID = FALSE, weight = TRUE)
	# recentAbsAreaInGap_m2 <- sum(apply(recentAbsAreaInGap_m2, 1, prod))

	# sink('./Figures & Tables/Area of Gap.txt', split = TRUE)
		# say(date())
		# say('')
		# say('Area of pika gap estimated from previously-known occurrences using Delauney triangles extended to east/west extent of EPA III Sierra Nevada ecoregion is is ', round(gapArea_m2), ' km2.')
		# say('Area of Sierra Nevada ecoregion is is ', round(snArea_m2), ' km2.')
		# say('')
		# say('Area that was recently suitable in the gap but is now unoccupied (area predicted to be "recent absence" class) is ', round(recentAbsAreaInGap_m2), ' km2.')
		# say('')
		# say('Current gap : Sierra Nevada total: ', round(gapArea_m2 / snArea_m2, 3))
		# say('Suitable area lost / gap area: ', round(recentAbsAreaInGap_m2 / gapArea_m2, 3))
		# say('Gap was until "recently" ', round(100 * (gapArea_m2 - gapArea_m2 / (1 + recentAbsAreaInGap_m2 / gapArea_m2)) / gapArea_m2, 1), '% smaller.')
		
	# sink()

# say('###########################################################')
# say('### compare inter-point distances between sets of sites ###')
# say('###########################################################')
	
	# say('This section explores the distribution of inter-site distances.')

	# testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# testSurveys <- vect(testSurveys, geom = ll, crs = prismCrs)
	
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.rda')
	# trainPres <- vect(trainPres, geom = c('longWgs84', 'latWgs84'), crs = prismCrs)

	# testPres <- testSurveys[testSurveys$status == '2 detected']
	# recentAbs <- testSurveys[testSurveys$status == '1 recent absence']
	# oldAbs <- testSurveys[testSurveys$status == '0 long absence']

	# trainPres_m <- distance(trainPres)
	# testPres_m <- distance(testPres)
	# recentAbs_m <- distance(recentAbs)
	# oldAbs_m <- distance(oldAbs)
	
	# trainPres_m <- as.matrix(trainPres_m)
	# testPres_m <- as.matrix(testPres_m)
	# recentAbs_m <- as.matrix(recentAbs_m)
	# oldAbs_m <- as.matrix(oldAbs_m)
	
	# trainPres_m[upper.tri(trainPres_m)] <- NA
	# testPres_m[upper.tri(testPres_m)] <- NA
	# recentAbs_m[upper.tri(recentAbs_m)] <- NA
	# oldAbs_m[upper.tri(oldAbs_m)] <- NA
	
	# trainPres_m <- c(trainPres_m)
	# testPres_m <- c(testPres_m)
	# recentAbs_m <- c(recentAbs_m)
	# oldAbs_m <- c(oldAbs_m)
	
	# trainPres_m <- trainPres_m[!is.na(trainPres_m)]
	# testPres_m <- testPres_m[!is.na(testPres_m)]
	# recentAbs_m <- recentAbs_m[!is.na(recentAbs_m)]
	# oldAbs_m <- oldAbs_m[!is.na(oldAbs_m)]
	
	# trainPresVsTestPres_m <- distance(trainPres, testPres)
	# trainPresVsRecentAbs_m <- distance(trainPres, recentAbs)
	# trainPresVsOldAbs_m <- distance(trainPres, oldAbs)
	
	# testPresVsRecentAbs_m <- distance(testPres, recentAbs)
	# testPresVsOldAbs_m <- distance(testPres, oldAbs)
	
	# recentAbsVsOldAbs_m <- distance(recentAbs, oldAbs)

	# trainPresVsTestPres_m <- c(trainPresVsTestPres_m)
	# trainPresVsRecentAbs_m <- c(trainPresVsRecentAbs_m)
	# trainPresVsOldAbs_m <- c(trainPresVsOldAbs_m)
	# testPresVsRecentAbs_m <- c(testPresVsRecentAbs_m)
	# testPresVsOldAbs_m <- c(testPresVsOldAbs_m)
	# recentAbsVsOldAbs_m <- c(recentAbsVsOldAbs_m)
	
	
	# dists <- data.frame(
		# Group = c(
			# # rep('Previously known occurrences', length(trainPres_m)),
			# rep('Currently occupied', length(testPres_m)),
			# rep('Formerly occupied', length(recentAbs_m)),
			# rep('No evidence', length(oldAbs_m)),
			# rep('Previously known occs. vs. currently occupied', length(trainPresVsTestPres_m)),
			# rep('Previously known occs. vs. formerly occupied', length(trainPresVsRecentAbs_m)),
			# rep('Previously known occs. vs. no evidence', length(trainPresVsOldAbs_m)),
			# rep('Currently occupied vs. formerly occupied', length(testPresVsRecentAbs_m)),
			# rep('Currently occupied vs. no evidence', length(testPresVsOldAbs_m)),
			# rep('Formerly occupied vs. no evidence', length(recentAbsVsOldAbs_m))
		# ),
		# Set = c(
			# # rep('vs. Self', length(trainPres_m) + length(testPres_m) + length(recentAbs_m) + length(oldAbs_m)),
			# rep('vs. Self', length(testPres_m) + length(recentAbs_m) + length(oldAbs_m)),
			# rep('vs. Previously known', length(trainPresVsTestPres_m) + length(trainPresVsRecentAbs_m) + length(trainPresVsOldAbs_m)),
			# rep('Test vs. Test', length(testPresVsRecentAbs_m) + length(testPresVsOldAbs_m) + length(recentAbsVsOldAbs_m))
		# ),
		# dists = c(
			# # trainPres_m,
			# testPres_m,
			# recentAbs_m,
			# oldAbs_m,
			# trainPresVsTestPres_m,
			# trainPresVsRecentAbs_m,
			# trainPresVsOldAbs_m,
			# testPresVsRecentAbs_m,
			# testPresVsOldAbs_m,
			# recentAbsVsOldAbs_m
		# )
	# )
	
	# dists$dists <- dists$dists / 1000
	
	# vsSelf <- dists[dists$Set == 'vs. Self', ]
	# vsSelfPlot <- ggplot(vsSelf, aes(x = dists, color = Group, fill = Group)) +
		# geom_density(alpha = 0.4) +
		# scale_x_continuous(trans = 'sqrt') +
		# xlab('Distance (km)') + ylab('Density') +
		# ggtitle('Pairwise point distances\nbetween sites in same group')
		
	# testVsTest <- dists[dists$Set == 'Test vs. Test', ]
	# testVsTestPlot <- ggplot(testVsTest, aes(x = dists, color = Group, fill = Group)) +
		# geom_density(alpha = 0.4) +
		# scale_x_continuous(trans = 'sqrt') +
		# xlab('Distance (km)') + ylab('Density') +
		# ggtitle('Pairwise point distances\nbetween sites in different groups')
		
	# combo <- vsSelfPlot + testVsTestPlot
	# ggsave(combo, file = './Figures & Tables/Inter-point Distances.png', width = 12, height = 6, dpi = 600)
	
# say('######################################################################')
# say('### pairwise distances between surveyed test sites with detections ###')
# say('######################################################################')	

	# testSurveys <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# testSurveys <- vect(testSurveys, geom = ll, crs = prismCrs)
	# detectSurveys <- testSurveys[testSurveys$status == '2 detected']
	
	# dists <- distance(detectSurveys)
	# dists <- as.matrix(dists)
	# dists <- dists / 1000
	# dists <- round(dists, 3)
	# diag(dists) <- NA

	# minDist_km <- apply(dists, 2, min, na.rm = TRUE)
	
	# sink('./Figures & Tables/Inter-point Distances between Test Sites with Detections in km.txt', split = TRUE)
	# say('Minimum pairwise distances between each survey test site with a detection and any other survey site with a detection (km):')
	# say(minDist_km)
	# sink()

# say('############################')
# say('### response curve plots ###')
# say('############################')

	# # taking advantage of 1:1 correspondence between PC scores, predictions, and each environmental datum

	# # models
	# load('./ENMs/ENM BRT.rda')
	# brt <- model
	# load('./ENMs/ENM GAM.rda')
	# gam <- model
	# load('./ENMs/ENM GLM.rda')
	# glm <- model

	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.rda')
	# load('./ENMs/Background Sites/Random Background Sites from across Study Region Mask.rda')
	
	# testSites <- fread('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.rda')
	
	# bgScores <- predict(pca, bg)[ , 1]
	# trainScores <- predict(pca, trainPres)[ , 1]
	# testScores <- predict(pca, testSites)[ , 1]
	
	# # predictions
	# bgPreds <- predictEnsemble(bg, brt = brt, gam = gam, glm = glm)
	# trainPreds <- predictEnsemble(trainPres, brt = brt, gam = gam, glm = glm)
	# testPreds <- predictEnsemble(testSites, brt = brt, gam = gam, glm = glm)
	
	# resps <- list()
	# for (pred in predictors) {
		
		# # background sites
		# bgDf <- data.table(
			# var = bg[ , pred],
			# pc = bgScores,
			# prediction = bgPreds,
			# status = 'background site'
		# )
		
		# # training presences
		# trainDf <- data.table(
			# var = trainPres[ , pred],
			# pc = trainScores,
			# prediction = trainPreds,
			# status = 'previously-known presence'
		# )
		
		# # test sites
		# testDf <- data.table(
			# var = testSites[[pred]],
			# pc = testScores,
			# prediction = testPreds,
			# status = testSites$status
		# )
		# testDf <- testDf[order(status)]
		# testDf$status[testDf$status == '0 long absence'] <- 'long-term test absence'
		# testDf$status[testDf$status == '1 recent absence'] <- 'recent test absence'
		# testDf$status[testDf$status == '2 detected'] <- 'test presence'
		
		# df <- rbind(bgDf, trainDf, testDf)
		# bgDf <- bgDf[order(var)]
		
		# resps[[length(resps) + 1]] <- ggplot(df, aes(x = var, y = prediction, color = status, size = status, fill = status, shape = status)) +
			# geom_point() +
			# scale_color_manual(
				# values = c(
					# 'background site' = 'gray',
					# 'previously-known presence' = 'black',
					# 'long-term test absence' = 'black',
					# 'recent test absence' = 'black',
					# 'test presence' = 'black'
				# )
			# ) +
			# scale_fill_manual(
				# values = c(
					# 'background site' = 'gray',
					# 'previously-known presence' = 'black',
					# 'long-term test absence' = 'red',
					# 'recent test absence' = 'yellow',
					# 'test presence' = 'chartreuse'
				# )
			# ) +
			# scale_size_manual(
				# values = c(
					# 'background site' = 1,
					# 'previously-known presence' = 1.4,
					# 'long-term test absence' = 1.4,
					# 'recent test absence' = 1.4,
					# 'test presence' = 1.4
				# )
			# ) +
			# scale_shape_manual(
				# values = c(
					# 'background site' = 16,
					# 'previously-known presence' = 3,
					# 'long-term test absence' = 25,
					# 'recent test absence' = 23,
					# 'test presence' = 21
				# )
			# ) +
			# geom_smooth(data = bgDf, se = FALSE, aes(color = 'gray40')) +
			# xlab(pred) +
			# ylab('Predicted suitability') +
			# ylim(c(0, 1))

	# }
	
	# legend <- get_legend(
	  # # create some space to the left of the legend
	  # resps[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
	# )
	
	# respsSansLegend <- resps
	# for (i in seq_along(respsSansLegend)) respsSansLegend[[i]] <- respsSansLegend[[i]] + theme(legend.position = 'none')

	# responses <- plot_grid(plotlist = respsSansLegend, nrow = 3, labels = LETTERS[seq_along(resps)])
	# all <- plot_grid(responses, legend, ncol = 1, rel_heights = c(1, 0.15))
	# ggsave(all, filename = './Figures & Tables/Response Curves.png', width = 8, height = 11, dpi = 300)
	
say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!', level=1, pre=1)
