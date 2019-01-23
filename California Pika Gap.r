### Ochotona princeps - Spatially-varying importance of variables
### Adam B. Smith | 2017-02

# source('C:/ecology/Drive/Research/Iconic Species/Analysis - California Gap/Code/California Pika Gap.r')

	rm(list=ls())
	memory.limit(memory.limit() * 2^30)
	gc()

	set.seed(1234567890)
	
### CONTENTS ###
### libraries, variables, and functions ###
### define study region ###
### process predictor rasters ###
### collate training detections and non-detections ###
### collate test detections and non-detections ###
### construct KDE on training occurrences to identify gap ###
### map of gap sampling ###
### assign weights to training sites based on spatial autocorrelation ###
### calculate spatial autocorrelation within test and between test and training survey sites ###
### plot training sites scaled by SAC-based weight ###
### plot test sites scaled by SAC-based weight ###
### generate background sites and calculate PCA ###
### PCA biplot ###
### match predictors with training and test sites ###
### train ENMs ###
### assess predictor importance ###
### write ENM rasters ###
### extract predictions to test survey sites ###
### calculate performance statistics against test survey sites ###
### boxplot of distribution of predictions at test sites ###
### PCA of background and test site classes ###
### visual comparison of background, training, and test sites by predictor ###
### map of predictions ###

######################
### generalization ###
######################


###########################################
### libraries, variables, and functions ###
###########################################

	setwd('C:/ecology/Drive/Research/Iconic Species/Analysis - California Gap')

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

		library(openxlsx)
		
		library(omnibus)
		library(enmSdm)
		library(legendary)
		library(statisfactory)
		
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
		
			raster::stack(c(listFiles('./Data/Climate - Derived', pattern='.tif'), './Data/NDVI 1990-2010/ndvi.tif'))
		
		}

		# get vector of PCs to use for predictors
		whichPcs <- function(pca, thold=0.95) {
			
			# pca		object of class "princomp"
			# thold		cumulative variance explaned must be <= this for PC to be used
			
			varExplained <- pca$sdev^2 / sum(pca$sdev^2)
			cumVarExplained <- cumsum(varExplained)
			
			usePcs <- which(cumVarExplained <= 0.95)
			usePcs
		
		}

	# ensemble prediction across ENMs
	predictEnsemble <- function(data, brt=NULL, gam=NULL, glm=NULL) {
		
		# data		data frame
		# brt, gam, glm	model objects, declare to speed up
		
		preds <- matrix(NA, ncol=length(algos), nrow=nrow(data))
		
		if (is.null(brt)) load('./ENMs/ENM BRT.RData')
		brt <- model
		if (is.null(gam)) load('./ENMs/ENM GAM.RData')
		gam <- model
		if (is.null(glm)) load('./ENMs/ENM GLM.RData')
		glm <- model
	
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
		predictors <- c('eastness', 'northness', 'acuteCold', 'chronicHeat', 'gsaisr', 'gswb', 'swe', 'ndvi')
		
		# length of distance bins for calculating characteristic scale of clustering of survey sites (in meters)
		sacInterval <- 20000
		
		# gray color palette
		grays <- paste0('gray', 0:100)
	
		# ENM algorithms
		algos <- c('brt', 'glm', 'gam')
	
	###############
	### options ###
	###############
	
		options(stringsAsFactors=FALSE)
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

	# ## study region extent derived from modified EPA Level III "Sierra Nevada" ecoregion plus 200-km buffer
	# say('study region')
		
		# epa3 <- shapefile('C:/Ecology/Drive/Research/Iconic Species/Extents_Masks_Maps/EcoRegions/!OmernikEcoregions/us_eco_l3SarrREV')

		# snEa <- epa3[epa3$L3_KEY == '5  Sierra Nevada', ]
		# snBuff <- gBuffer(snEa, width=1000 * studyRegionBuffer)

		# snBuff <- sp::spTransform(snBuff, getCRS('prism', TRUE))
		# sn <- sp::spTransform(snEa, getCRS('prism', TRUE))

		# shapefile(snBuff, './Study Region/epa3_sierraNevadaSarrModified_200kmBuffer', overwrite=TRUE)
		# shapefile(sn, './Study Region/epa3_sierraNevadaSarrModified', overwrite=TRUE)
		
		# rm(epa3); gc()
		
	# ### PRISM elevation
	# say('PRISM elevation')
	
		# elevPrism <- raster('D:/Ecology/Climate/PRISM/30 arcsec/elevation.tif')
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
	
		# srtms <- listFiles('D:/Ecology/Topography/SRTM', pattern='.tif')
		
		# srtm <- raster(srtms[1])
		# for (i in 2:length(srtms)) {
			# add <- raster(srtms[i])
			# srtm <- mosaic(srtm, add, fun=max)
		# }
	
		# projection(srtm) <- getCRS('wgs84')
	
		# snBuffWgs84 <- sp::spTransform(snBuff, getCRS('wgs84', TRUE))
		# srtmSn <- crop(srtm, snBuff)
		# names(srtmSn) <- 'elevation_srtm_m'
	
		# dirCreate('./Data/Topography - SRTM')
		# writeRaster(srtmSn, './Data/Topography - SRTM/elevation_srtm_m', datatype='INT2S')
		
		# rm(srtmSn); gc()
		
	# ### SRTM hill shade
	# say('SRTM hill shade')
	
		# snBuffLarger <- gBuffer(snEa, width=2 * 1000 * studyRegionBuffer)
		# snBuffLargerWgs84 <- sp::spTransform(snBuffLarger, getCRS('wgs84', TRUE))
		# expandedSrtm <- crop(srtm, snBuffLargerWgs84)
	
		# slope <- terrain(expandedSrtm, 'slope')
		# aspect <- terrain(expandedSrtm, 'aspect')
		
		# names(slope) <- 'slope_srtm_rad'
		# names(aspect) <- 'aspect_srtm_rad'
		
		# hs <- hillShade(slope, aspect, direction=45)
		# hs <- round(hs, 5)
		# names(hs) <- 'hillshade_srtm'
		# writeRaster(hs, './Data/Topography - SRTM/hillshade_srtm', datatype='FLT4S')
		
		# beginCluster(4)
			# hsEa <- projectRaster(hs, crs=getCRS('climateNA'))
		# endCluster()
		# hs <- round(hsEa, 5)
		# names(hsEa) <- 'hillshade_srtm_ea'
		# writeRaster(hsEa, './Data/Topography - SRTM/hillshade_srtm_ea', datatype='FLT4S')

	# ## political entities
	# say('political geography')
	
		# # california state
		# gadm <- raster::getData('GADM', level=1, country='USA', path='C:/ecology/!Scratch')
		# west1 <- gadm[gadm$NAME_1 %in% c('California', 'Nevada', 'Oregon'), ]
		# save(west1, file='./Study Region/GADM California, Nevada, Oregon States.RData')
		
		# # counties
		# gadm <- raster::getData('GADM', level=2, country='USA', path='C:/ecology/!Scratch')
		# west2 <- gadm[gadm$NAME_1 %in% c('California', 'Nevada', 'Oregon'), ]
		# plumas <- gadm[gadm$NAME_1 == 'California' & gadm$NAME_2 == 'Plumas', ]
		# sierra <- gadm[gadm$NAME_1 == 'California' & gadm$NAME_2 == 'Sierra', ]
		# washoe <- gadm[gadm$NAME_1 == 'Nevada' & gadm$NAME_2 == 'Washoe', ]
		
		# nearGap <- gTouches(plumas, west2, byid=TRUE)
		# gapCounties <- west2[c(nearGap), ]
		# gapCounties <- rbind(plumas, sierra, washoe, gapCounties)
		
		# save(west2, file='./Study Region/GADM California, Nevada, Oregon Counties.RData')
		# save(plumas, file='./Study Region/GADM Plumas County.RData')
		# save(gapCounties, file='./Study Region/GADM GAP Counties.RData')

		# # plumas <- sp::spTransform(plumas, getCRS('climateNA', TRUE))
		# # plumasBuffer <- gBuffer(plumas, width=35000)
		# # plumasBuffer <- sp::spTransform(plumasBuffer, CRS(projection(gadm)))
		
		# # save(plumasBuffer, file='./Study Region/GADM Plumas County + 35-km Buffer.RData')

# say('#################################')
# say('### process predictor rasters ###')
# say('#################################')

	# # generalization 
	# prismDriveLetter <- 'F'

	# # # calculate time period for climate layers
	# # load('./Data/Occurrences for California Gap.RData')
	
	# # dirCreate('./Figures & Tables/Occurrences')
	# # png(paste0('./Figures & Tables/Occurrences/Year of Observations.png'), height=900, width=900, res=200)
		# # hist(pres$obsYear, breaks=(min(pres$obsYear) - 1):(max(pres$obsYear) + 1), main='Observations', xlab='Year')
	# # dev.off()

	# say('Using 20-yr window for climate layers (1996-2015).')
	
	# years <- 1996:2015
	
	# say('Using: MEAN MONTHLY PPT, TMIN, TMAX, ET0, AISR, and SWE')
	
	# mask <- raster(paste0('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus', studyRegionBuffer, 'kmBuffer.tif'))

	# # calculate monthly means of max/min temperature, et0, actual incoming solar radiation, and snow water equivalent
	# for (variable in c('ppt', 'tmin', 'tmax', 'et0', 'aisr', 'swe')) {
		
		# dirCreate('./Data/Climate - Monthly Means/', variable)

		# for (month in 1:12) {
		
			# say(variable, ' | month ', month, ' | year', post=0)
		
			# for (year in years) {
		
				# say(year, post=ifelse(year == 2015, 1, 0))

				# # mid-day of this month
				# midMonthDoy <- doyNonLeap[15, paste0('month', month)]
			
				# thisMonthYear <- if (variable %in% c('ppt', 'tmin', 'tmax')) {
					# raster(paste0(prismDriveLetter, ':/ecology/Climate/PRISM/AN81m 1981-2015/', variable, '/', year, '/prism_', variable, '_us_30s_', year, prefix(month, 2), '.tif'))
				# } else if (variable == 'et0') {
					# raster(paste0(prismDriveLetter, ':/ecology/Climate/PRISM/ET0 - Monthly - Dobrowski et al 2013/et0_year', year, '_month', prefix(month, 2), '.tif'))
				# } else if (variable == 'aisr') {
					# raster(paste0(prismDriveLetter, ':/ecology/Climate/PRISM/Actual Incoming Solar Radiation - Dobrowski et al 2013/globalActualDownwardShortwaveRadiation_year', year, '_month', prefix(month, 2), '_doy', prefix(midMonthDoy, 3), '_MJperM2perDay.tif'))
				# } else if (variable == 'swe') {
					# raster(paste0(prismDriveLetter, ':/ecology/Climate/PRISM/Snow Water Equivalent - Monthly - Dobrowski et al 2013/snowWaterEquivalent_year', year, '_month', prefix(month, 2), '_mm.tif'))
				# }
				
				# thisMonthYear <- crop(thisMonthYear, mask)
				
				# if (!exists('thisMonth', inherits=FALSE)) {
					# thisMonth <- stack(thisMonthYear)
				# } else {
					# thisMonth <- stack(thisMonth, thisMonthYear)
				# }
				
			# } # next year
			
			# # calculate mean for this month across years
			# meanForMonth <- mean(thisMonth)
			# meanForMonth <- setMinMax(meanForMonth)
			# names(meanForMonth) <- variable
			
			# writeRaster(meanForMonth, paste0('./Data/Climate - Monthly Means/', variable, '/', variable, '_month', prefix(month, 2), '_mean1996to2015'))
			
			# rm(thisMonth, meanForMonth); gc()
			
		# } # next month
		
	# } # next variable
	
	# dirCreate('./Data/Climate - Derived')
	
	# say('ASPECT')
	
		# aspect <- raster(paste0(prismDriveLetter, ':/ecology/Climate/PRISM/PRISM_us_dem_800m_aspect_degrees.tif'))
		# aspect <- crop(aspect, mask)
		# aspect <- pi * aspect / 180
		# northness <- sin(aspect)
		# eastness <- cos(aspect)
		# names(northness) <- 'northness'
		# names(eastness) <- 'eastness'
		
		# writeRaster(northness, './Data/Climate - Derived/northness')
		# writeRaster(eastness, './Data/Climate - Derived/eastness')
		
	# say('CHRONIC HEAT')
	
		# rasts <- stack(listFiles('./Data/Climate - Monthly Means/tmax', pattern='.tif'))
		# rasts <- subset(rasts, 6:9)
		# chronicHeat <- mean(rasts) / 100000
		# chronicHeat <- setMinMax(chronicHeat)
		# names(chronicHeat) <- 'chronicHeat'
		# writeRaster(chronicHeat, './Data/Climate - Derived/chronicHeat')
		
	# say('ACUTE COLD')
	
		# rasts <- stack(listFiles('./Data/Climate - Monthly Means/tmin', pattern='.tif'))
		# acuteCold <- min(rasts)
		# acuteCold <- acuteCold / 100000
		# acuteCold <- setMinMax(acuteCold)
		# names(acuteCold) <- 'acuteCold'
		# writeRaster(acuteCold, './Data/Climate - Derived/acuteCold')
		
	# say('SWE')

		# rasts <- stack(listFiles('./Data/Climate - Monthly Means/swe', pattern='.tif'))
		# swe <- log10(mean(rasts) + 1)
		# swe <- setMinMax(swe)
		# names(swe) <- 'swe'
		# writeRaster(swe, './Data/Climate - Derived/swe')
		
	# say('GSAISR')

		# rasts <- stack(listFiles('./Data/Climate - Monthly Means/aisr', pattern='.tif'))
		# rasts <- subset(rasts, 6:9)
		# gsaisr <- sum(rasts)
		# gsaisr <- setMinMax(gsaisr)
		# names(gsaisr) <- 'gsaisr'
		# writeRaster(gsaisr, './Data/Climate - Derived/gsaisr')
		
	# say('GSWB')

		# ppt <- stack(listFiles('./Data/Climate - Monthly Means/ppt', pattern='.tif'))
		# et0 <- stack(listFiles('./Data/Climate - Monthly Means/et0', pattern='.tif'))
		# ppt <- subset(ppt, 6:9)
		# et0 <- subset(et0, 6:9)
		
		# gswb <- ppt - et0
		
		# gswb <- sum(gswb)
		# gswb <- gswb / 100000
		# gswb <- setMinMax(gswb)
		# names(chronicHeat) <- 'gswb'
		# writeRaster(gswb, './Data/Climate - Derived/gswb')
	
	# say('NDVI')
		
		# ndvi <- raster('./Data/NDVI 1990-2010/Raw/NDVI.tif')
		# beginCluster(4)
			# ndvi <- projectRaster(ndvi, mask)
		# endCluster()
		# names(ndvi) <- 'ndvi'
		# writeRaster(ndvi, paste0('./Data/NDVI 1990-2010/ndvi'))

	# # say('PLOT ALL PREDICTORS')
	
	# # pres <- readRDS(paste0('./Presences for California Gap.rds'))
	# # elevPrism <- raster(paste0('./Predictors/elevPrism_studyRegionExtent.tif'))
	
	# # usa1 <- readOGR('C:/ecology/Political Geography/GADM/ver2pt8/WGS84', 'USA_adm1', verbose=FALSE)
	# # usa2 <- readOGR('C:/ecology/Political Geography/GADM/ver2pt8/WGS84', 'USA_adm2', verbose=FALSE)
	
	# # plumas <- usa2[usa2$NAME_2 == 'Plumas', ]
	
	# # studyRegion <- readOGR(paste0('./Study Region'), 'Study Region - EPA Level III Sierra Nevada Plus 200-km Buffer', verbose=FALSE)

	# # say('study region')
	# # png(paste0('./Maps/Predictor Maps - Study Region.png'), width=4 * 250, height= 2 * 350, res=200)
	
		# # par(mfrow=c(2, 4), mgp=c(1, 0.2, 0), mar=c(3, 2, 4, 1) + 0.1, cex.main=0.6, cex.lab=0.5, cex.axis=0.5, tck=-0.02)

		# # overlays <- function() {
		
			# # plot(usa2, border='gray40', lwd=0.7, add=TRUE)
			# # plot(usa1, border='gray20', lwd=0.7, add=TRUE)
			# # plot(plumas, border='blue', lwd=0.7, add=TRUE)
			# # plot(studyRegion, border='darkgreen', lwd=0.7, add=TRUE)
			# # points(pres$longWgs84, pres$latWgs84, pch=16, cex=0.2, col='red')
			
		# # }
		
		# # plot(elevPrism, main='Elevation'); overlays()
		# # plot(northness, main='Northness'); overlays()
		# # plot(chronicHeat, main='Chronic Heat'); overlays()
		# # plot(acuteCold, main='Acute Cold'); overlays()
		# # plot(swe, main='SWE'); overlays()
		# # plot(gswb, main='Growing Season\nAtmospheric Water Balance'); overlays()
		# # plot(gsaisr, main='Growing Season Actual\nIncoming Solar Radiation'); overlays()
		# # plot(ndvi, main='NDVI'); overlays()
		
	# # dev.off()
	
	# # say('gap')
	# # png(paste0('./Maps/Predictor Maps - Gap.png'), width=4 * 300, height= 2 * 300, res=200)
	
		# # par(mfrow=c(2, 4), mgp=c(1, 0.2, 0), mar=c(3, 2, 4, 1) + 0.1, cex.main=0.6, cex.lab=0.5, cex.axis=0.5, tck=-0.02)

		# # overlays <- function() {
		
			# # plot(usa2, border='gray40', lwd=0.7, add=TRUE)
			# # plot(usa1, border='gray20', lwd=0.7, add=TRUE)
			# # plot(plumas, border='blue', lwd=0.7, add=TRUE)
			# # plot(studyRegion, border='darkgreen', lwd=0.7, add=TRUE)
			# # points(pres$longWgs84, pres$latWgs84, pch=16, cex=0.2, col='red')
			
		# # }
		
		# # plot(plumas, main='Elevation'); plot(elevPrism, add=TRUE); overlays()
		# # plot(plumas, main='Northness'); plot(northness, add=TRUE); overlays()
		# # plot(plumas, main='Chronic Heat'); plot(chronicHeat, add=TRUE); overlays()
		# # plot(plumas, main='Acute Cold'); plot(acuteCold, add=TRUE); overlays()
		# # plot(plumas, main='SWE'); plot(swe, add=TRUE); overlays()
		# # plot(plumas, main='Growing Season\nAtmospheric Water Balance'); plot(gswb, add=TRUE); overlays()
		# # plot(plumas, main='Growing Season Actual\nIncoming Solar Radiation'); plot(gsaisr, add=TRUE); overlays()
		# # plot(plumas, main='NDVI'); plot(ndvi, add=TRUE); overlays()
		
	# # dev.off()

# say('######################################################')
# say('### collate training detections and non-detections ###')
# say('######################################################')	
	
	# dirCreate('./Data/Occurrences')
	# surveys <- readRDS('C:/Ecology/Drive/Research/Iconic Species/Species Records - Pika/!Collated Data 2016-06-30 1256/00 Pika - Cleaned Using R - Usable Records for All Ochotona.rds')
	# save(surveys, file='./Data/Occurrences/Training Surveys 00 Pika - Cleaned Using R - Usable Records for All Ochotona.RData')

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

	# save(surveys, file='./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.RData')

# say('##################################################')
# say('### collate test detections and non-detections ###')
# say('##################################################')
	
	# ### clean
	# file.copy('C:/Ecology/Drive/Research/Iconic Species/Analysis - California Gap/Data/Occurrences/Gap Surveys - Raw/TR_Edits_NC hole pika   locations__updated16Sept2016_KBO__Status....xlsx', './Data/Occurrences/Test Surveys 00 TR_Edits_NC hole pika   locations__updated16Sept2016_KBO__Status....xlsx')
	
	# testSurveys <- read.xlsx('./Data/Occurrences/Test Surveys 01 Edits in Excel to Make Machine-Readable.xlsx', sheet='All but unsuitable, Sorted Elev')
	
	# testSurveys <- testSurveys[!is.na(testSurveys$status), ]

	# # ### plot
	# # load('./Study Region/GADM California, Nevada, Oregon Counties.RData')
	# # load('./Study Region/GADM California, Nevada, Oregon States.RData')
	
	# # testSurveysSp <- SpatialPointsDataFrame(testSurveys[ , ll], data=testSurveys, proj4=getCRS('wgs84', TRUE))
	
	# # plot(testSurveysSp, col='white')
	# # plot(west2, border='gray', add=TRUE)
	# # plot(west1, add=TRUE)
	# # points(testSurveysSp, pch=1)
	
	# write.csv(testSurveys, './Data/Occurrences/Test Surveys 02 Cleaned.csv', row.names=FALSE)
	
# say('#############################################################')
# say('### construct KDE on training occurrences to identify gap ###')
# say('#############################################################')

	# # occurrences
	# load('./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.RData')
	# pres <- surveys[surveys$origRecentPikaOrSignsObserved, ]
	# presSp <- SpatialPointsDataFrame(pres[ , ll], data=pres, proj4=getCRS('wgs84', TRUE))
	# presSpEa <- sp::spTransform(presSp, getCRS('climateNA', TRUE))
	
	# # mask
	# mask <- raster(paste0('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus', studyRegionBuffer, 'kmBuffer.tif'))
	# maskEa <- projectRaster(mask, crs=getCRS('climateNA'))
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
	# load('./Study Region/GADM California, Nevada, Oregon States.RData')
	# load('./Study Region/GADM California, Nevada, Oregon Counties.RData')
	# load('./Study Region/GADM Plumas County.RData')
	# load(file='./Study Region/GADM GAP Counties.RData')
	# hs <- raster('./Data/Topography - SRTM/hillshade_srtm_ea.tif')
	
	# gapCounties <- sp::spTransform(gapCounties, getCRS('climateNA', TRUE))
	# plumas <- sp::spTransform(plumas, getCRS('climateNA', TRUE))
	# west1 <- sp::spTransform(west1, getCRS('climateNA', TRUE))
	# west2 <- sp::spTransform(west2, getCRS('climateNA', TRUE))

	# # training occurrences
	# load('./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.RData')
	# trainPres <- surveys[surveys$origRecentPikaOrSignsObserved, ]
	# trainPresSp <- SpatialPointsDataFrame(trainPres[ , ll], data=trainPres, proj4=getCRS('wgs84', TRUE))
	# trainPresSpEa <- sp::spTransform(trainPresSp, getCRS('climateNA', TRUE))

	# # test occurrences
	# testSurveys <- read.csv('./Data/Occurrences/Test Surveys 02 Cleaned.csv')
	# testSurveysSp <- SpatialPointsDataFrame(testSurveys[ , ll], data=testSurveys, proj4=getCRS('nad83', TRUE))
	# testSurveysSpEa <- sp::spTransform(testSurveysSp, getCRS('climateNA', TRUE))
	
	# # plot focus
	# focus <- gBuffer(testSurveysSpEa, width=1000 * focusBuff)
	
	# # kde
	# kde <- raster('./KDE/kde.tif')
	# kdeVals <- extract(kde, trainPresSpEa)
	# quants <- quantile(kdeVals, c(0.01, 0.05))
	# breaks <- c(0, min(kdeVals), quants, 1)
	# kdeClass <- cut(kde, breaks=breaks)

	# # colors
	# kdeCols <- c(NA, 'steelblue1', 'steelblue3', 'steelblue4')
	# for (i in seq_along(kdeCols)) kdeCols[i] <- alpha(kdeCols[i], 0.5)
	
	# # plot
	# dirCreate('./Figures & Tables/Gap Sampling')
	# png('./Figures & Tables/Gap Sampling/Gap Sampling.png', width=1200, height=1200, res=300)
		
		# trainPresCol <- '#1b9e77'
		# testPresCol <- '#7570b3'
		# testAbsCol <- '#d95f02'
		
		# par(mar=c(2, 2, 2, 2))
		
		# plot(focus, border='white')
		# plot(hs, col=grays, legend=FALSE, add=TRUE)
		# plot(kdeClass, col=kdeCols, legend=FALSE, add=TRUE)
		# plot(west2, border='gray40', add=TRUE)
		# plot(west1, lwd=2, add=TRUE)
		# points(trainPresSpEa, pch=1, cex=0.6, bg='white')

		# testSurveysSpEa <- testSurveysSpEa[order(testSurveysSpEa$status), ]
		
		# bgs <- rep(NA, nrow(testSurveys))
		# bgs[testSurveysSpEa$status == '0 long absence'] <- 'red'
		# bgs[testSurveysSpEa$status == '1 recent absence'] <- 'yellow'
		# bgs[testSurveysSpEa$status == '2 detected'] <- 'chartreuse'
		
		# pchs <- rep(NA, nrow(testSurveys))
		# pchs[testSurveysSpEa$status == '0 long absence'] <- 22
		# pchs[testSurveysSpEa$status == '1 recent absence'] <- 22
		# pchs[testSurveysSpEa$status == '2 detected'] <- 21
		
		# points(testSurveysSpEa, pch=pchs, bg=bgs, cex=0.6)
		# box()
		
		# # legend
		# legendBreaks(
			# 'bottomleft',
			# inset=0.01,
			# height=0.44,
			# width=0.31,
			# title='Training density',
			# titleAdj=c(0.5, 0.92),
			# col=kdeCols,
			# adjX=c(0.05, 0.225),
			# adjY=c(0.38, 0.85),
			# labels=c('\U2265Min presence', paste0('\U2265', '1st percentile'), paste0('\U2265', '5th percentile')),
			# labAdjX=0.18,
			# cex=0.45
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
			# pch=c(1, 22, 22, 21),
			# # col=c('black', 'black', testAbsCol),
			# pt.bg=c(NA, 'red', 'yellow', 'chartreuse'),
			# bty='n',
			# title='Surveys',
			# cex=0.45,
			# pt.cex=0.8
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
# say('### assign weights to training sites based on spatial autocorrelation ###')
# say('#########################################################################')

	# say('I want to use inverse-p weighting applied to survey sites to correct for spatial sampling bias. I will assume that two survey sites at the exact same location each have a weight of 1/2, three each have a weight of 1/3, and so on. At the other end of the spectrum survey sites that are far enough away should each have a weight of 1. I will define "far enough away" as the distance at which the proportion of pairwise observed distances falls below the upper 95th quantile of the distribution of pairwise distances from randomly located sites (one-tail). I will draw a number of randomly located sites in each iteration so it is equal to the total number of survey sites.', breaks=80)

	# ### generalization
	# minSurveys <- 0 # minimum number of survey sites by one provider to use in SAC assessment
	
	# ### mask
	# mask <- raster(paste0('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus', studyRegionBuffer, 'kmBuffer.tif'))

	# ### occurrences
	# load('./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.RData')
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
	
	# save(trainPres, file='./Data/Occurrences/Training Surveys 02 Pika - Assigned Weights Based on Pairwise Distance Distribution.RData')

# say('################################################################################################')
# say('### calculate spatial autocorrelation within test and between test and training survey sites ###')
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
	# load('./Data/Occurrences/Training Surveys 02 Pika - Assigned Weights Based on Pairwise Distance Distribution.RData')
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
# say('### plot training sites scaled by SAC-based weight ###')
# say('######################################################')

	# # spatial data
	# study <- shapefile('./Study Region/epa3_sierraNevadaSarrModified_200kmBuffer')
	# sn <- shapefile('./Study Region/epa3_sierraNevadaSarrModified')
	# load('./Study Region/GADM California, Nevada, Oregon States.RData')
	# load('./Study Region/GADM California, Nevada, Oregon Counties.RData')
	
	# # training occurrences
	# load('./Data/Occurrences/Training Surveys 02 Pika - Assigned Weights Based on Pairwise Distance Distribution.RData')
	# sacDist <- read.csv('./Figures & Tables/Spatial Autocorrelation between Survey Sites/!Characteristic Clustering Distance for Training Sites.csv')
	
	# trainPresSp <- SpatialPointsDataFrame(trainPres[ , ll], data=trainPres, proj4=getCRS('prism', TRUE))

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
# say('### plot test sites scaled by SAC-based weight ###')
# say('##################################################')

	# # spatial data
	# study <- shapefile('./Study Region/epa3_sierraNevadaSarrModified_200kmBuffer')
	# sn <- shapefile('./Study Region/epa3_sierraNevadaSarrModified')
	# load('./Study Region/GADM California, Nevada, Oregon States.RData')
	# load('./Study Region/GADM California, Nevada, Oregon Counties.RData')
	
	# # test occurrences
	# sites <- read.csv('./Data/Occurrences/Test Surveys 03 Assigned Weights Based on Pairwise Spatial Autocorrelation.csv')
	# sacDist <- read.csv('./Figures & Tables/Spatial Autocorrelation between Survey Sites/!Characteristic Clustering Distance for Test vs Test and Test vs Training Sites.csv')
	
	# sitesSp <- SpatialPointsDataFrame(sites[ , ll], data=sites, proj4=getCRS('prism', TRUE))

	# hs <- raster('./Data/Topography - SRTM/hillshade_srtm.tif')
	
	# # training occurrences
	# load('./Data/Occurrences/Training Surveys 02 Pika - Assigned Weights Based on Pairwise Distance Distribution.RData')
	
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
	
	# save(bg, file='./ENMs/Background Sites/Random Background Sites from across Study Region Mask.RData')
	# save(pca, file='./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.RData')

# say('##################')
# say('### PCA biplot ###')
# say('##################')

	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.RData')
	
	# totalVar <- sum(pca$sdev^2)
	# varPc1 <- pca$sdev[1]^2 / totalVar
	# varPc2 <- pca$sdev[2]^2 / totalVar
	
	# xlab <- paste0('PC1 (', sprintf('%.1f', 100 * varPc1), '%)')
	# ylab <- paste0('PC2 (', sprintf('%.1f', 100 * varPc2), '%)')
	
	# png('./Figures & Tables/Background Environment as PCA.png', width=1000, height=1000, res=200)
	
		# par(mar=c(4, 4, 4, 4), cex.lab=0.6, cex.axis=0.6)
		
		# plot(pca$scores[ , 1:2], type='p', pch=16, col=alpha('darkblue', 0.05), xlab=xlab, ylab=ylab)
		
		# stretch <- 3
		# for (i in 1:nrow(pca$loadings)) {
		
			# arrows(x0=0, y0=0, x1=stretch * pca$loadings[i, 1], y1=stretch * pca$loadings[i, 2], angle=15, length=0.07)
			# text(stretch * pca$loadings[i, 1], stretch * pca$loadings[i, 2], labels=rownames(pca$loadings)[i])
			
		# }
		
	# dev.off()
	
	
# say('#####################################################')
# say('### match predictors with training and test sites ###')
# say('#####################################################')

	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.RData')

	# load('./Data/Occurrences/Training Surveys 02 Pika - Assigned Weights Based on Pairwise Distance Distribution.RData')
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
	
	# save(trainPres, file='./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.RData')
	# write.csv(testSurveys, './Data/Occurrences/Test Surveys 04 Extracted Environmental Values with PCA.csv', row.names=FALSE)
	
# say('##################')
# say('### train ENMs ###')
# say('##################')	

	# # training presences and background
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.RData')
	# load('./ENMs/Background Sites/Random Background Sites from across Study Region Mask.RData')
	
	# # PCA
	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.RData')
	
	# bg$weight <- 1
	
	# pcs <- paste0('pc', whichPcs(pca))
	
	# ### collate training sites
	# train <- rbind(
		# trainPres[ , c(pcs, 'weight')],
		# bg[ , c(paste0(pcs, 'weight')]
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
	# model <- trainBrt(train, resp='presBg', preds=pcs, w='weight', verbose=TRUE)
	# save(model, file='./ENMs/ENM BRT.RData')
	
	# # ### CRF
	# # model <- trainCrf(train, resp='presBg', preds=pcs, w='weight', verbose=TRUE)
	# # save(model, file='./ENMs/ENM CRF.RData')
	
	# ### GLM
	# model <- trainGlm(train, resp='presBg', preds=pcs, w='weight', verbose=TRUE)
	# save(model, file='./ENMs/ENM GLM.RData')
	
	# ### GAM
	# model <- trainGam(train, resp='presBg', preds=pcs, w='weight', verbose=TRUE)
	# save(model, file='./ENMs/ENM GAM.RData')

# say('###################################')
# say('### assess predictor importance ###')
# say('###################################')

	# # PCA
	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.RData')
	# pcs <- paste0('pc', whichPcs(pca))
	
	# # training presences and background
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.RData')
	# load('./ENMs/Background Sites/Random Background Sites from across Study Region Mask.RData')
	
	# # PCA
	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.RData')
	
	# bg$weight <- 1
	
	# pcs <- paste0('pc', whichPcs(pca))
	
	# ### collate training sites
	# train <- rbind(
		# trainPres[ , pcs],
		# bg[ , pcs]
	# )

	# load(paste0('./ENMs/ENM BRT.RData'))
	# brt <- model
				
	# load(paste0('./ENMs/ENM GAM.RData'))
	# gam <- model
				
	# load(paste0('./ENMs/ENM GLM.RData'))
	# glm <- model
				
	# obsPred <- predictEnsemble(train, brt=brt, gam=gam, glm=glm)

	# varImp <- matrix(NA, ncol=length(pcs), nrow=100)
	# colnames(varImp) <- pcs
	
	# for (pc in seq_along(pcs)) {
	
		# say(paste0('pc', pc))
	
		# corTest <- rep(NA, 100)
		# for (iter in 1:100) {
		
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

	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.RData')

	# # stack of PC rasters
	# env <- stackEnv()
	
	# envPc <- predict(env, pca, index=1:nlayers(env))
	# names(envPc) <- paste0('pc', 1:nlayers(envPc))
	# envPc <- subset(envPc, whichPcs(pca))

	# # predict
	# if (exists('predStack')) rm(predStack)
	
	# for (algo in toupper(algos)) {
	
		# say(algo)
	
		# load(paste0('./ENMs/ENM ', algo, '.RData'))
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
# say('### extract predictions to test survey sites ###')
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
# say('### calculate performance statistics against test survey sites ###')
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
		# say('Presences over recent absences: ...... ', sprintf('%.3f', aucs[['pres_over_recentAbs']]))
		# say('Presences over long-term absences: ... ', sprintf('%.3f', aucs[['pres_over_longTermAbs']]))
		# say('Recent over long-term absences: ...... ', sprintf('%.3f', aucs[['recentAbs_over_longTermAbs']]), post=2)
		
		# say('Number of presences: ..................', numPres)
		# say('Number of recent absences: ............', numRecentAbs)
		# say('Number of long-term absences: .........', numLongTermAbs)
		
		# say(date())
	# sink()
	
# say('############################################################')
# say('### boxplot of distribution of predictions at test sites ###')
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
# say('### PCA of background and test site classes ###')
# say('###############################################')
	
	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.RData')
	
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.RData')
	# testSurveys <- read.csv('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	
	# pcs <- c('pc1', 'pc2')

	# ### figure with background as background
	
	# smoothScatter(pca$scores[ , 1:2], pch=16, nrpoints=0, xlab='PC 1', ylab='PC 2')
	# points(trainPres[ , pcs], pch=16, col=alpha('black', 0.1))
	
	# pres <- testSurveys[testSurveys$status == '2 detected', ]
	# recentAbs <- testSurveys[testSurveys$status == '1 recent absence', ]
	# longTermAbs <- testSurveys[testSurveys$status == '0 long absence', ]
	
	# points(longTermAbs[ , pcs], pch=25, bg='red')
	# points(recentAbs[ , pcs], pch=23, bg='yellow')
	# points(pres[ , pcs], pch=21, bg='green')
	
	# legend('topleft', inset=0.01, legend=c('Background', 'Training presence', 'Test presence', 'Test recent absence', 'Test long-term absence'), pch=c(NA, 16, 21, 23, 25), fill=c(blues9[6], NA, NA, NA, NA), pt.bg=c(NA, NA, 'green', 'yellow', 'red'), border=c('black', NA, NA, NA, NA), bg='white')

# say('##############################################################################')
# say('### visual comparison of background, training, and test sites by predictor ###')
# say('##############################################################################')	
	
	# dirCreate('./Figures & Tables/Comparison of Environment across Sites')

	# load('./ENMs/Background Sites/PCA on Random Background Sites from across Study Region Mask.RData')
	# pcs <- paste0('pc', whichPcs(pca))

	# load('./ENMs/Background Sites/Random Background Sites from across Study Region Mask.RData')
	# load('./Data/Occurrences/Training Surveys 03 Pika - Extracted Environmental Values with PCA.RData')
	# testSurveys <- read.csv('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	
	# for (pc in pcs) {
		
		# png(paste0('./Figures & Tables/Comparison of Environment across Sites/', toupper(pc), '.png'), width=1000, height=800, res=300)
			
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
			
			# plot(thisBg, bgCs, xlab=toupper(pc), ylab='Cumulative proportion', type='l', col='blue', lwd=2, main=toupper(pc))
			# lines(thisTrainPres, trainPresCs, col='black', lwd=2)
			# lines(longTermAbs, longTermAbsCs, col='red', lwd=2)
			# lines(recentAbs, recentAbsCs, col='yellow', lwd=2)
			# lines(detected, detectedCs, col='chartreuse', lwd=2)
			
			# legend('topleft', inset=0.01, legend=c('Background', 'Training presences', 'Test presences', 'Recent absences', 'Long-term absences'), col=c('blue', 'black', 'chartreuse', 'yellow', 'red'), lwd=2, bty='n', cex=0.4)
		
		# dev.off()
	
	# }

say('##########################')
say('### map of predictions ###')
say('##########################')

	# generalize
	focusBuff <- 25 # buffer size around test sites for generating focus of plot (in km)

	### spatial data
	load('./Study Region/GADM California, Nevada, Oregon States.RData')
	load('./Study Region/GADM California, Nevada, Oregon Counties.RData')
	load('./Study Region/GADM Plumas County.RData')
	load(file='./Study Region/GADM GAP Counties.RData')
	hs <- raster('./Data/Topography - SRTM/hillshade_srtm_ea.tif')
	
	gapCounties <- sp::spTransform(gapCounties, getCRS('climateNA', TRUE))
	plumas <- sp::spTransform(plumas, getCRS('climateNA', TRUE))
	west1 <- sp::spTransform(west1, getCRS('climateNA', TRUE))
	west2 <- sp::spTransform(west2, getCRS('climateNA', TRUE))

	### training occurrences
	load('./Data/Occurrences/Training Surveys 01 Pika - Selected Detections and Non-Detections from Data Providers in Study Region.RData')
	trainPres <- surveys[surveys$origRecentPikaOrSignsObserved, ]
	trainPresSp <- SpatialPointsDataFrame(trainPres[ , ll], data=trainPres, proj4=getCRS('wgs84', TRUE))
	trainPresSpEa <- sp::spTransform(trainPresSp, getCRS('climateNA', TRUE))

	### test occurrences
	testSurveys <- read.csv('./Data/Occurrences/Test Surveys 05 Extracted Predictions.csv')
	testSurveysSp <- SpatialPointsDataFrame(testSurveys[ , ll], data=testSurveys, proj4=getCRS('nad83', TRUE))
	testSurveysSpEa <- sp::spTransform(testSurveysSp, getCRS('climateNA', TRUE))
	
	### plot focus
	focus <- gBuffer(testSurveysSpEa, width=1000 * focusBuff)
	
	### ENM prediction
	pred <- raster('./ENMs/predictionEnsemble.tif')
	beginCluster(4)
		pred <- projectRaster(pred, crs=getCRS('climateNA', TRUE))
	endCluster()
	predVals <- extract(pred, trainPresSpEa)
	
	### colors
	trainPresCol <- '#1b9e77'
	testPresCol <- '#7570b3'
	testAbsCol <- '#d95f02'
		
	### plot
	png('./Figures & Tables/ENM Prediction.png', width=2 * 1200, height=2 * 1200, res=300)
		
		par(mfrow=c(2, 2), mar=c(1, 1, 1, 1))

		### each class by itself
		
		for (class in c(0:2)) {
		
			if (class == 0) {
				predClass <- pred >= min(predVals[testSurveys$status == '0 long absence']) &
					pred <= max(predVals[testSurveys$status == '0 long absence'])
				col <- alpha('red', 0.5)
				lab <- 'Long-term absence predicted'
			} else if (class == 1) {
				predClass <- pred >= min(predVals[testSurveys$status == '1 recent absence']) &
					pred <= max(predVals[testSurveys$status == '1 recent absence'])
				col <- alpha('darkgoldenrod1', 0.5)
				lab <- 'Recent absence predicted'
			} else if (class == 2) {
				predClass <- pred >= min(predVals[testSurveys$status == '2 detected'])
				col <- alpha('forestgreen', 0.5)
				lab <- 'Presence predicted'
			}
			
			plot(focus, border='white')
			plot(hs, col=grays, legend=FALSE, add=TRUE)
			plot(predClass, col=c(NA, col), legend=FALSE, add=TRUE)
			plot(west2, border='gray40', add=TRUE)
			plot(west1, lwd=2, add=TRUE)
			points(trainPresSpEa, pch=1.2, cex=1, bg='white')

			bgs <- rep(NA, nrow(testSurveys))
			bgs[testSurveysSpEa$status == '0 long absence'] <- 'red'
			bgs[testSurveysSpEa$status == '1 recent absence'] <- 'yellow'
			bgs[testSurveysSpEa$status == '2 detected'] <- 'chartreuse'
			
			pchs <- rep(NA, nrow(testSurveys))
			
			pchs[testSurveysSpEa$status == '0 long absence'] <- 24
			pchs[testSurveysSpEa$status == '1 recent absence'] <- 24
			pchs[testSurveysSpEa$status == '2 detected'] <- 21

			if (class == 0) {
				
				points(testSurveysSpEa[testSurveysSpEa$status == '2 detected', ], bg='chartreuse', pch=21, cex=1.5)
				points(testSurveysSpEa[testSurveysSpEa$status == '1 recent absence', ], bg='yellow', pch=24, cex=1.2)
				points(testSurveysSpEa[testSurveysSpEa$status == '0 long absence', ], bg='red', pch=24, cex=1.2)
				
			} else if (class == 1) {
			
				points(testSurveysSpEa[testSurveysSpEa$status == '0 long absence', ], bg='red', pch=24, cex=1.2)
				points(testSurveysSpEa[testSurveysSpEa$status == '2 detected', ], bg='chartreuse', pch=21, cex=1.5)
				points(testSurveysSpEa[testSurveysSpEa$status == '1 recent absence', ], bg='yellow', pch=24, cex=1.2)
				
			} else if (class == 2) {
			
				points(testSurveysSpEa[testSurveysSpEa$status == '0 long absence', ], bg='red', pch=24, cex=1.2)
				points(testSurveysSpEa[testSurveysSpEa$status == '1 recent absence', ], bg='yellow', pch=24, cex=1.2)
				points(testSurveysSpEa[testSurveysSpEa$status == '2 detected', ], bg='chartreuse', pch=21, cex=1.5)
				
			}
				
			box()
			
			legend(
				'bottomleft',
				inset=0.01,
				legend=c(lab, 'Long-term test absence', 'Recent test absence', 'Test presence', 'Training presence'),
				pch=c(NA, 24, 24, 21, 1),
				fill=c(col, NA, NA, NA, NA),
				border=c('black', NA, NA, NA, NA),
				pt.bg=c(NA, 'red', 'yellow', 'chartreuse', NA),
				bg=alpha('white', 0.3),
				cex=0.75,
				pt.cex=1.2
			)
			
			# scale bar
			size <- 50000 # length of scale bar in meters
			x <- usr[2] - size - 0.02 * width
			x <- c(x, x + size)
			y <- usr[3] + rep(0.02 * height, 2)
			
			lines(x, y, lwd=4, xpd=NA, col='black', lend=1)
			
			x <- mean(x)
			y <- usr[3] + 0.05 * height
			text(x, y[1], labels=paste(size / 1000, 'km'), cex=1)
			
		}

		### all together

		predLongTerm <- pred >= min(predVals[testSurveys$status == '0 long absence']) &
			pred <= max(predVals[testSurveys$status == '0 long absence'])
		predRecent <- pred >= min(predVals[testSurveys$status == '1 recent absence']) &
			pred <= max(predVals[testSurveys$status == '1 recent absence'])
		predPresent <- pred >= min(predVals[testSurveys$status == '2 detected'])
			
		predClass <- predLongTerm + predRecent + predPresent
		
		predCols <- c(NA, 'red', 'yellow', 'forestgreen')
		for (i in seq_along(predCols)) predCols[i] <- alpha(predCols[i], 0.5)
		
		plot(focus, border='white')
		plot(hs, col=grays, legend=FALSE, add=TRUE)
		plot(predClass, col=predCols, legend=FALSE, add=TRUE)
		plot(west2, border='gray40', add=TRUE)
		plot(west1, lwd=2, add=TRUE)
		points(trainPresSpEa, pch=1.2, cex=1, bg='white')

		testSurveysSpEa <- testSurveysSpEa[order(testSurveysSpEa$status), ]
		
		bgs <- rep(NA, nrow(testSurveys))
		bgs[testSurveysSpEa$status == '0 long absence'] <- 'red'
		bgs[testSurveysSpEa$status == '1 recent absence'] <- 'yellow'
		bgs[testSurveysSpEa$status == '2 detected'] <- 'chartreuse'
		
		cexs <- pchs <- rep(NA, nrow(testSurveys))
		
		pchs[testSurveysSpEa$status == '0 long absence'] <- 24
		pchs[testSurveysSpEa$status == '1 recent absence'] <- 24
		pchs[testSurveysSpEa$status == '2 detected'] <- 21
		
		cexs[testSurveysSpEa$status == '0 long absence'] <- 1.2
		cexs[testSurveysSpEa$status == '1 recent absence'] <- 1.2
		cexs[testSurveysSpEa$status == '2 detected'] <- 1.5

		points(testSurveysSpEa, bg=bgs, pch=pchs, cex=cexs)
			
		box()
		
		legend(
			'bottomleft',
			inset=0.01,
			legend=c('Long-term absence predicted', 'Recent absence predicted', 'Presence predicted', 'Long-term test absence', 'Recent test absence', 'Test presence', 'Training presence'),
			pch=c(NA, NA, NA, 24, 24, 21, 1),
			fill=c(alpha('red', 0.5), alpha('darkgoldenrod1', 0.5), alpha('forestgreen', 0.5), NA, NA, NA, NA),
			border=c('black', 'black', 'black', NA, NA, NA, NA),
			pt.bg=c(NA, NA, NA, 'red', 'yellow', 'chartreuse', NA),
			bg=alpha('white', 0.3),
			cex=0.75,
			pt.cex=1.2
		)
		
		# scale bar
		size <- 50000 # length of scale bar in meters
		x <- usr[2] - size - 0.02 * width
		x <- c(x, x + size)
		y <- usr[3] + rep(0.02 * height, 2)
		
		lines(x, y, lwd=4, xpd=NA, col='black', lend=1)
		
		x <- mean(x)
		y <- usr[3] + 0.05 * height
		text(x, y[1], labels=paste(size / 1000, 'km'), cex=1)

		title(sub=date(), cex.sub=0.3, outer=TRUE, line=-1)
		
	dev.off()
	
	
say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!', level=1, pre=1)
