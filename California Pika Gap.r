### Ochotona princeps - Spatially-varying importance of variables
### Adam B. Smith | Missouri Botanical Garden | adam.smith@mobot.org | 2017-02
###
### source('C:/Ecology/Drive/Research/Pikas - California Gap (Erik Beever et al)/Code/California Pika Gap.r')
### source('E:/Ecology/Drive/Research/Pikas - California Gap (Erik Beever et al)/Code/California Pika Gap.r')
###
### CONTENTS ###
### libraries, variables, and functions ###
### define study region ###
### process predictor rasters ###
### collate TRAINING detections and non-detections ###
### collate TEST detections and non-detections ###
### construct KDE on TRAINING occurrences to identify gap ###

### map of gap sampling ###

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

######################
### generalization ###
######################


###########################################
### libraries, variables, and functions ###
###########################################

	rm(list=ls())
	memory.limit(memory.limit() * 2^30)
	gc()

	set.seed(1234567890)
	
	# working drive
	# drive <- 'C:'
	# drive <- 'D:'
	drive <- 'E:'

	# PRISM
	# prismDrive <- 'F:'
	prismDrive <- 'I:'
	
	# TerraClimate
	tcDrive <- 'F:'

	setwd(paste0(drive, '/ecology/Drive/Research/Pikas - California Gap (Erik Beever et al)'))

	prismCrs <- '+proj=longlat +datum=NAD83 +no_defs'

	###############################
	### libraries and functions ###
	###############################

		library(scales)
		library(sp)
		library(rgdal)
		library(mgcv)
		library(geosphere)
		library(cluster)
		library(rJava)
		library(fossil)
		library(rgeos)
		library(raster)
		library(dismo)
		library(spatialEco)
		library(gbm)
		library(party)
		library(terra)

		library(openxlsx)
		
		library(omnibus)
		library(enmSdm)
		library(legendary)
		library(statisfactory)
		
		rasterOptions(format='GTiff', overwrite=TRUE, tmpdir=paste0(drive ,'/ecology/!Scratch/_raster'))
		
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
		grassDir <- c('C:/OSGeo4W64/', 'grass-7.4.1', 'osgeo4W')
	
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
		rasterOptions(format='GTiff', overwrite=TRUE)

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
		
		# epa3 <- shapefile(paste0(drive, '/Ecology/Drive/Research/Iconic Species/Extents_Masks_Maps/EcoRegions/!OmernikEcoregions/us_eco_l3SarrREV'))

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
	
		# snBuffLarger <- gBuffer(snEa, width=2 * 1000 * studyRegionBuffer)
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
		
		# slope <- raster(slope)
		# aspect <- raster(aspect)
		
		# hs <- hillShade(slope, aspect, direction=45)
		# names(hs) <- 'hillshade_srtm'
		# writeRaster(hs, './Data/Topography - SRTM/hillshade_srtm', overwrite=TRUE)
		
		# beginCluster(4)
			# hsEa <- projectRaster(hs, crs=getCRS('climateNA'))
		# endCluster()
		# hs <- round(hsEa, 5)
		# names(hsEa) <- 'hillshade_srtm_ea'
		# writeRaster(hsEa, './Data/Topography - SRTM/hillshade_srtm_ea', datatype='FLT4S')

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

# say('#################################')
# say('### process predictor rasters ###')
# say('#################################')

	# ### calculate raster representing mean monthly values for each variable across a set time period
	# ### from these calculate specific predictor variable rasters

	# # calculate time period for climate layers
	# load('./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.rda')
	
	# png(paste0('./Figures & Tables/Year of Observations.png'), height=900, width=900, res=200)
		# hist(surveys$obsYear, breaks=(min(surveys$obsYear) - 1):(max(surveys$obsYear) + 1), main='Observations', xlab='Year')
	# dev.off()

	# mask <- terra::rast(paste0('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus', studyRegionBuffer, 'kmBuffer.tif'))

	# data('doyNonLeap', package='omnibus')

	# # calculate monthly means of max/min temperature and ppt from PRISM and TerraClimate
	# for (variable in c('tmin', 'tmax', 'srad', 'def', 'swe', 'vpdmax', )) { # ALL
	# # for (variable in c('tmin', 'tmax', 'vpdmax')) { # PRISM variables
	# # for (variable in c('srad', 'def', 'swe')) { # TerraClimate variables
		
		# dirCreate('./Data/Climate - Monthly Means/', variable)

		# months <- if (variable %in% c('tmax', 'srad', 'def', 'vpdmax')) {
			# 6:9
		# } else if (variable %in% c( 'tmin', 'swe')) {
			# 1:12
		# }

		# for (month in months) {
		
			# say(variable, ' | month ', month, ' | year', post=0)
		
			# for (year in years) {
		
				# say(year, post=ifelse(year == max(years), 1, 0))

				# # mid-day of this month
				# midMonthDoy <- doyNonLeap[15, paste0('month', month)]
			
				# if (variable %in% c('ppt', 'tmin', 'tmax', 'vpdmax')) {
					# thisMonthYear <- terra::rast(paste0(prismDrive, '/ecology/Climate/PRISM/working/an81/' , variable, '/monthly/', year, '/prism_', variable, '_us_30s_' , year, prefix(month, 2), '.tif'))
				# } else if (variable %in% c('pet', 'srad', 'swe')) {
					# thisMonthYear <- terra::rast(paste0(tcDrive, '/ecology/Climate/TerraClimate/ORIGINALS/', variable, '/TerraClimate_', variable, '_', year, '.nc'))
					# thisMonthYear <- thisMonthYear[[month]]
				# }
				
				# thisMonthYear <- terra::resample(thisMonthYear, mask)
				# thisMonthYear <- terra::crop(thisMonthYear, mask)
				
				# thisMonth <- if (!exists('thisMonth', inherits=FALSE)) {
					# thisMonthYear
				# } else {
					# c(thisMonth, thisMonthYear)
				# }
				
			# } # next year
			
			# # calculate mean for this month across years
			# meanForMonth <- mean(thisMonth)
			# names(meanForMonth) <- variable
			
			# terra::writeRaster(meanForMonth, paste0('./Data/Climate - Monthly Means/', variable, '/', variable, '_month', prefix(month, 2), '_mean', min(years), 'to', max(years), '.tif'), overwrite=TRUE, wopt=c(wopt, names=variable))
			
			# rm(thisMonth, meanForMonth); gc()
			
		# } # next month
		
	# } # next variable

	# dirCreate('./Data/Climate - Derived')
	
	# say('ASPECT')
	
		# elev <- terra::rast(paste0(drive, '/Ecology/Climate/PRISM/PRISM_us_dem_800m.tif'))
		# aspect <- terra::terrain(elev, v='aspect')
		# aspect <- terra::crop(aspect, mask)
		# aspect <- pi * aspect / 180
		# northness <- sin(aspect)
		# eastness <- cos(aspect)
		# names(northness) <- 'northness'
		# names(eastness) <- 'eastness'
		
		# terra::writeRaster(northness, './Data/Climate - Derived/northness.tif', overwrite=TRUE, wopt=c(wopt, names='northness'))
		# terra::writeRaster(eastness, './Data/Climate - Derived/eastness.tif', overwrite=TRUE, wopt=c(wopt, names='eastness'))
		
	# say('CHRONIC HEAT')
	
		# rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/tmax', pattern='.tif'))
		# chronicHeat <- mean(rasts)
		# terra::writeRaster(chronicHeat, './Data/Climate - Derived/chronicHeat_C.tif', overwrite=TRUE, wopt=c(wopt, names='chronicHeat_C'))
		
	# say('SUMMER NIGHTTIME HEAT')
	
		# rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/tmin', pattern='.tif'))
		# rasts <- rasts[[6:9]]
		# summerNightHeat <- mean(rasts)
		# terra::writeRaster(summerNightHeat, './Data/Climate - Derived/summerNightTemp_C.tif', overwrite=TRUE, wopt=c(wopt, names='summerNightTemp_C'))
		
	# say('ACUTE COLD')
	
		# rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/tmin', pattern='.tif'))
		# acuteCold <- min(rasts)
		# terra::writeRaster(acuteCold, './Data/Climate - Derived/acuteCold_C.tif', overwrite=TRUE, wopt=c(wopt, names='acuteCold_C'))
		
	# say('VPDMAX')
	
		# rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/vpdmax', pattern='.tif'))
		# vpdmax_haPa <- mean(rasts)
		# terra::writeRaster(vpdmax_haPa, './Data/Climate - Derived/vpdmax_haPa.tif', overwrite=TRUE, wopt=c(wopt, names='vpdmax_haPa'))
		
	# say('SWE')

		# rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/swe', pattern='.tif'))
		# swe <- sum(rasts)
		# terra::writeRaster(swe, './Data/Climate - Derived/swe_mm.tif', overwrite=TRUE, wopt=c(wopt, names='swe_mm'))
		
	# say('SUMMER SRAD')

		# rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/srad', pattern='.tif'))
		# srad <- sum(rasts)
		# terra::writeRaster(srad, './Data/Climate - Derived/summerSrad.tif', overwrite=TRUE, wopt=c(wopt, names='summerSrad'))

	# say('GSWD')

		# rasts <- terra::rast(listFiles('./Data/Climate - Monthly Means/def', pattern='.tif'))
		# gswb <- sum(rasts)
		# terra::writeRaster(gswb, './Data/Climate - Derived/growSeasonWaterDef_mm.tif', overwrite=TRUE, wopt=c(wopt, names='growSeasonWaterDef_mm'))
	
	# say('NDVI')
		
		# ndvi <- terra::rast('./Data/NDVI 1990-2010/Raw/NDVI.tif')
		# ndvi <- mean(ndvi)
		# ndvi <- terra::project(ndvi, mask)
		# ndvi <- crop(ndvi, mask)
		# terra::writeRaster(ndvi, paste0('./Data/NDVI 1990-2010/ndvi.tif'), overwrite=TRUE, wopt=c(wopt, names='ndvi'))

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
# say('### collate training detections and non-detections ###')
# say('######################################################')	
	
	# dirCreate('./Data/Occurrences')
	# surveys <- readRDS(paste0(drive, '/Ecology/Drive/Research/Iconic Species/Species Records - Pika/!Collated Data 2016-06-30 1256/00 Pika - Cleaned Using R - Usable Records for All Ochotona.rds'))

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
			# surveys$contact == 'Jason Brewer & Mary Flores' |
			# surveys$contact == 'Clint Epps, Jessica Castillo' |
			# surveys$contact == 'Jess Castillo, Clint Epps' |
			# surveys$contact == 'Connie Millar'
		# ), ]

	# ## subset to presences in study region
	# studyRegion <- raster(paste0('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus', studyRegionBuffer, 'kmBuffer.tif'))
	# inStudyRegion <- raster::extract(mask, surveys[ , ll])
	
	# surveys <- surveys[!is.na(inStudyRegion), ] # just in Sierra Nevada + buffer study region

	# save(surveys, file='./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.rda')

# say('##################################################')
# say('### collate test detections and non-detections ###')
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

	# # occurrences
	# load('./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.rda')
	# pres <- surveys[surveys$origRecentPikaOrSignsObserved, ]
	# presSp <- SpatialPointsDataFrame(pres[ , ll], data=pres, proj4=getCRS('wgs84', TRUE))
	# presSpEa <- sp::spTransform(presSp, getCRS('albersNA', TRUE))
	
	# # mask
	# mask <- raster(paste0('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus', studyRegionBuffer, 'kmBuffer.tif'))
	# maskEa <- projectRaster(mask, crs=getCRS('albersNA'))
	# maskEa <- round(maskEa)
	
	# ### KDE
	
	# # # Gaussian bandwidth
	# # bw <- 0.5 * sum(apply(coordinates(presSpEa), 2, sd)) * nrow(presSp)^(-1/6)
	
	# # Epanechnikov bandwidth
	# bw <- 1.77 * 0.5 * sum(apply(coordinates(presSpEa), 2, sd, na.rm=T)) * nrow(presSp)^(-1/6)
	
	# kde <- spatialEco::sp.kde(presSpEa, bw=bw, newdata=maskEa, standardize=TRUE)
	
	# writeRaster(kde, './KDE/kde')
	
# say('###########################')
# say('### map of gap sampling ###')
# say('###########################')

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
	# png('./Figures & Tables/Gap Sampling.png', width=1200, height=1200, res=300)
		
		# par(mar=c(2, 1, 1, 1), cex.axis=0.4, mgp=c(3, 0, 0.2))
		
		# plot(focus, border='white')
		
		# usr <- par('usr')
		# xs <- pretty(c(usr[1], usr[2]))
		# ys <- pretty(c(usr[3], usr[4]))
		# axis(1, at=xs, tck=0.01, labels=xs, col=NA, col.ticks='black')
		# axis(2, at=ys, tck=0.01, labels=ys, col=NA, col.ticks='black')
		
		# plot(hs, col=grays, legend=FALSE, add=TRUE)
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
	# testSites <- read.csv('./Data/Occurrences/Test Surveys 02 Cleaned.csv')
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
	# study <- shapefile('./Study Region/epa3_sierraNevadaSarrModified_200kmBuffer')
	# sn <- shapefile('./Study Region/epa3_sierraNevadaSarrModified')
	# load('./Study Region/GADM California, Nevada, Oregon States.rda')
	# load('./Study Region/GADM California, Nevada, Oregon Counties.rda')
	
	# # training occurrences
	# load('./Data/Occurrences/Training Surveys 02 Pika - Assigned Weights Based on Pairwise Distance Distribution.rda')
	# sacDist <- read.csv('./Figures & Tables/Spatial Autocorrelation between Survey Sites/!Characteristic Clustering Distance for Training Sites.csv')
	
	# trainPresSp <- SpatialPointsDataFrame(trainPres[ , ll], data=trainPres, proj4=CRS(prismCrs))

	# providers <- sacDist$provider
	
	# png('./Figures & Tables/Spatial Autocorrelation between Survey Sites/!Maps of Weighting by Provider.png', width=1600, height=800, res=200)
		
		# par(mfrow=c(2, 4), mar=c(1, 1, 4, 1))
		
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
	# study <- shapefile('./Study Region/epa3_sierraNevadaSarrModified_200kmBuffer')
	# sn <- shapefile('./Study Region/epa3_sierraNevadaSarrModified')
	# load('./Study Region/GADM California, Nevada, Oregon States.rda')
	# load('./Study Region/GADM California, Nevada, Oregon Counties.rda')
	
	# # test occurrences
	# sites <- read.csv('./Data/Occurrences/Test Surveys 03 Assigned Weights Based on Pairwise Spatial Autocorrelation.csv')
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
	# testSurveys <- read.csv('./Data/Occurrences/Test Surveys 03 Assigned Weights Based on Pairwise Spatial Autocorrelation.csv')
	
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
	# model <- trainBrt(train, resp='presBg', preds=pcs, w='weight', verbose=TRUE)
	# save(model, file='./ENMs/ENM BRT.rda')
	# rm(model); gc()
	
	# ### GLM
	# say('glm')
	# model <- trainGlm(train, resp='presBg', preds=pcs, w='weight', initialTerms=5, interaction=TRUE, verbose=TRUE)
	# save(model, file='./ENMs/ENM GLM.rda')
	# rm(model); gc()
	
	# ### GAM
	# say('gam')
	# model <- trainGam(train, resp='presBg', preds=pcs, w='weight', verbose=TRUE)
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
	
	# writeRaster(pred, './ENMs/predictionEnsemble', datatype='INT2U')
	# writeRaster(predStack, './ENMs/predictionByAlgorithm', datatype='INT2U')

# say('################################################')
# say('### extract predictions to TEST survey sites ###')
# say('################################################')	

	# ensPred <- raster('./ENMs/predictionEnsemble.tif')
	# algosPred <- stack('./ENMs/predictionByAlgorithm.tif')
	# names(algosPred) <- algos
	
	# testSurveys <- read.csv('./Data/Occurrences/Test Surveys 04 Extracted Environmental Values with PCA.csv')
	
	# predictEnsemble <- extract(ensPred, testSurveys[ , ll]) / 1000
	# predictAlgos <- extract(algosPred, testSurveys[ , ll]) / 1000
	
	# predicts <- cbind(predictEnsemble, predictAlgos)
	
	# testSurveys <- cbind(testSurveys, predicts)
	
	# write.csv(testSurveys, './Data/Occurrences/Test Surveys 05 Extracted Predictions.csv', row.names=FALSE)
	
# say('##################################################################')
# say('### calculate performance statistics against TEST survey sites ###')
# say('##################################################################')	

	# testSurveys <- read.csv('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
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

	# testSurveys <- read.csv('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	
	# testSurveys$weight <- pmin(testSurveys$weightTestVsTest, testSurveys$weightTestVsTrain)
	
	# png('./Figures & Tables/Predicted Suitability by Test Site Class - Unweighted.png', width=1200, height=1200, res=300)
		# par(mar=c(3, 4, 2, 1) + 0.1, cex.lab=0.8, cex.axis=0.8, cex.sub=0.4)
		# boxplot(predictEnsemble ~ status, data=testSurveys, names=c('Long-term\nabsence', 'Recent\nabsence', 'Presence'), ylim=c(0, 1), ylab='Predicted suitability')
		# title(sub=date(), line=2)
	# dev.off()
	
	# png('./Figures & Tables/Predicted Suitability by Test Site Class - Weighted.png', width=1200, height=1200, res=300)
		# par(mar=c(3, 4, 2, 1) + 0.1, cex.lab=0.8, cex.axis=0.8, cex.sub=0.4)
		# boxplot(predictEnsemble * weight ~ status, data=testSurveys, names=c('Long-term\nabsence', 'Recent\nabsence', 'Presence'), ylab='Predicted suitability (weighted)')
		# title(sub=date(), line=2)
	# dev.off()
	
# say('###############################################')
# say('### PCA of background and TEST site classes ###')
# say('###############################################')
	
	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.rda')
	
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.rda')
	# testSurveys <- read.csv('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	
	# pcs <- c('pc1', 'pc2')
	# cex <- 2.1

	# ### figure with background as background
	
	# png('./Figures & Tables/PCA with Test Classes.png', width=900, height=900)
	
		# par(cex.axis=1.8, cex.lab=2, c(5, 7, 4, 2) + 0.1)
		
		# smoothScatter(pca$scores[ , 1:2], pch=16, nrpoints=0, xlab='PC 1', ylab='PC 2')
		# points(trainPres[ , pcs], pch=3, col=alpha('black', 0.4), cex=1.6)
		
		# pres <- testSurveys[testSurveys$status == '2 detected', ]
		# recentAbs <- testSurveys[testSurveys$status == '1 recent absence', ]
		# longTermAbs <- testSurveys[testSurveys$status == '0 long absence', ]
		
		# points(longTermAbs[ , pcs], pch=25, bg=NA, col='red', cex=cex)
		# points(recentAbs[ , pcs], pch=23, bg=NA, col='yellow', cex=cex)
		# points(pres[ , pcs], pch=21, bg=NA, col='green', cex=cex)
		
		# legend('topleft', inset=0.01, legend=c('Background', 'Training presence', 'Test presence', 'Test recent absence', 'Test long-term absence'), pch=c(NA, 3, 1, 5, 6), fill=c(blues9[6], NA, NA, NA, NA), col=c(NA, 'black', 'green', 'yellow', 'red'), border=c('black', NA, NA, NA, NA), bg='white', cex=cex)
		
	# dev.off()
	
# say('##############################################################################')
# say('### visual comparison of background, TRAINING, and TEST sites by predictor ###')
# say('##############################################################################')	
	
	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.rda')
	# pcs <- paste0('pc', whichPcs(pca))

	# load('./ENMs/Background Sites/Random Background Sites from across Study Region Mask.rda')
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.rda')
	# testSurveys <- read.csv('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	
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

	# usa <- raster::getData('GADM', country='USA', level=1, path='C:/ecology/!Scratch')
	# mex <- raster::getData('GADM', country='MEX', level=1, path='C:/ecology/!Scratch')
	# nam <- rbind(usa, mex)
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
	# trainPres <- surveys[surveys$origRecentPikaOrSignsObserved, ]
	# trainPresSp <- SpatialPointsDataFrame(trainPres[ , ll], data=trainPres, proj4=getCRS('wgs84', TRUE))
	# trainPresSpEa <- sp::spTransform(trainPresSp, getCRS('climateNA', TRUE))

	# ### test occurrences
	# testSurveys <- read.csv('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
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
				# pointLab <- c('Correctly predicted', 'Old evidence predicted', 'Presence predicted')
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
				# pointLab <- c('No evidence predicted', 'Correctly predicted', 'Presence predicted')
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
				# pointLab <- c('No evidence predicted', 'Old evidence predicted', 'Correctly predicted')
				# pch <- testPresPch
				# axis(1, at=pretty(c(usr[1:2]), 3)[1:3], tck=-0.01)
				# axis(2, at=pretty(c(usr[3:4]), 3), tck=-0.01)
				
				# legColBg <- c(NA, testLongTermAbsCol, testShortTermAbsFill, NA)
				
			# }
			
			# plot(hs, col=grays, legend=FALSE, add=TRUE)
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
			# legend=c('No evidence predicted', 'Old evidence predicted', 'Occurrence predicted', 'No evidence observed', 'Old evidence observed', 'Occurrence observed', 'Training occurrence'),
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
	# testSurveys <- read.csv('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
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
	
say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!', level=1, pre=1)
