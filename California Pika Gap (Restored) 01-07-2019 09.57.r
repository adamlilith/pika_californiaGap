### Ochotona princeps - Spatially-varying importance of variables
### Adam B. Smith | 2017-02

# source('C:/ecology/Drive/Research/Iconic Species/Analysis - California Gap/California Pika Gap.r')

	rm(list=ls())
	memory.limit(memory.limit() * 2^30)
	gc()

### CONTENTS ###
### libraries, variables, and functions ###
### define study region ###
### process predictor rasters ###
### collate presences for KDE ###
### train KDE ###
### make figure of KDE ###

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

	#################
	### variables ###
	#################
	
		# doyNonLeapYear <- read.csv('C:/ecology/Drive/R/Ancillary Files/daysOfYear_nonLeapYear.csv')
		# doyLeapYear <- read.csv('C:/ecology/Drive/R/Ancillary Files/daysOfYear_leapYear.csv')

		ll <- c('longWgs84', 'latWgs84')
		
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
		
		# epa3 <- shapefile('C:/Ecology/Drive/Research/Iconic Species/Extents_Masks_Maps/EcoRegions/!OmernikEcoregions/us_eco_l3SarrREV')

		# sn <- epa3[epa3$L3_KEY == '5  Sierra Nevada', ]
		# snBuff <- gBuffer(sn, width=200000)

		# snBuff <- sp::spTransform(snBuff, getCRS('prism', TRUE))
		# sn <- sp::spTransform(sn, getCRS('prism', TRUE))

		# shapefile(snBuff, './Study Region/epa3_sierraNevadaSarrModified_200kmBuffer', overwrite=TRUE)
		# shapefile(sn, './Study Region/epa3_sierraNevadaSarrModified', overwrite=TRUE)
		
	# ### PRISM elevation
	
		# elevPrism <- raster('D:/Ecology/Climate/PRISM/30 arcsec/elevation.tif')
		# elevPrism <- crop(elevPrism, snBuff)
		# names(elevPrism) <- 'elevation_prism_m'
		
		# dirCreate('./Data/Elevation - PRISM')
		# writeRaster(elevPrism, './Data/Elevation - PRISM/elevation_prism_m_studyRegion')

	# ### mask
	
		# maskPrism <- 1 + elevPrism * 0
		# maskPrism <- rasterize(snBuff, maskPrism)
		# names(maskPrism) <- 'mask'
		# writeRaster(maskPrism, './Study Region/mask_prism_sierraNevadaEpaLevel3Plus200kmBuffer', datatype='INT1U')

	# ### SRTM elevation
	
		# srtms <- listFiles('./Data/Elevation - SRTM/Raw', pattern='.tif')
		
		# srtm <- raster(srtms[1])
		# for (i in 2:length(srtms)) {
			# add <- raster(srtms[i])
			# srtm <- mosaic(srtm, add, fun=max)
		# }
	
		# srtm <- crop(srtm, snBuff)
		# names(srtm) <- 'elevation_srtm_m'
	
		# writeRaster(srtm, './Data/Elevation - SRTM/elevation_srtm_m', datatype='INT2S')

	# ### SRTM hillshade
	
		# slope <- terrain(srtm, 'slope')
		# aspect <- terrain(srtm, 'aspect')
		
		# names(slope) <- 'slope_srtm_rad'
		# names(aspect) <- 'aspect_srtm_rad'
		
		# hs <- hillShade(slope, aspect)
		# names(hs) <- 'hillshade_srtm'
		# writeRaster(hs, './Data/Elevation - SRTM/hillshade_srtm')
		
	# ## political entities
		
		# # california state
		# gadm <- raster::getData('GADM', level=1, country='USA', path='C:/ecology/!Scratch')
		# califNev1 <- gadm[gadm$NAME_1 %in% c('California', 'Nevada'), ]
		# save(califNev1, file='./Study Region/GADM California & Nevada States.RData')
		
		# # counties
		# gadm <- raster::getData('GADM', level=2, country='USA', path='C:/ecology/!Scratch')
		# califNev2 <- gadm[gadm$NAME_1 %in% c('California', 'Nevada'), ]
		# plumas <- gadm[gadm$NAME_1 == 'California' & gadm$NAME_2 == 'Plumas', ]
		# sierra <- gadm[gadm$NAME_1 == 'California' & gadm$NAME_2 == 'Sierra', ]
		# washoe <- gadm[gadm$NAME_1 == 'Nevada' & gadm$NAME_2 == 'Washoe', ]
		
		# nearGap <- gTouches(plumas, califNev2, byid=TRUE)
		# gapCounties <- califNev2[c(nearGap), ]
		# gapCounties <- rbind(plumas, sierra, washoe, gapCounties)
		
		# save(califNev2, file='./Study Region/GADM California & Nevada Counties.RData')
		# save(plumas, file='./Study Region/GADM Plumas County.RData')
		# save(gapCounties, file='./Study Region/GADM GAP Counties.RData')

		# plumas <- sp::spTransform(plumas, getCRS('climateNA', TRUE))
		# plumasBuffer <- gBuffer(plumas, width=35000)
		# plumasBuffer <- sp::spTransform(plumasBuffer, CRS(projection(gadm)))
		
		# save(plumasBuffer, file='./Study Region/GADM Plumas County + 35-km Buffer.RData')

say('#################################')
say('### process predictor rasters ###')
say('#################################')

	# generalization 
	prismDriveLetter <- 'F'

	# calculate time period for climate layers
	load('./Data/Occurrences for California Gap.RData')
	
	dirCreate('./Figures & Tables/Occurrences')
	png(paste0('./Figures & Tables/Occurrences/Year of Observations.png'), height=900, width=900, res=200)
		hist(pres$obsYear, breaks=(min(pres$obsYear) - 1):(max(pres$obsYear) + 1), main='Observations', xlab='Year')
	dev.off()

	say('Using 20-yr window for climate layers (1996-2015).')
	
	years <- 1996:2015
	
	say('Using: MEAN MONTHLY PPT, TMIN, TMAX, ET0, AISR, and SWE')
	
	mask <- raster('./Study Region/mask_prism_sierraNevadaEpaLevel3Plus200kmBuffer.tif')

	# calculate monthly means of max/min temperature, et0, actual incoming solar radiation, and snow water equivalent
	for (variable in c('ppt', 'tmin', 'tmax', 'et0', 'aisr', 'swe')) {
		
		dirCreate('./Data/Climate - Monthly Means/', variable)

		for (month in 1:12) {
		
			say(variable, ' | month ', month, ' | year', post=0)
		
			for (year in years) {
		
				say(year, post=ifelse(year == 2015, 1, 0))

				# get mid-say of month (for AISR)
				doyThisYear <- if (year %% 4 == 0) {
					doyLeapYear
				} else {
					doyNonLeapYear
				}
			
				# mid-day of this month
				midMonthDoy <- doyNonLeapYear[15, paste0('month', month)]
			
				thisMonthYear <- if (variable %in% c('ppt', 'tmin', 'tmax')) {
					raster(paste0(prismDriveLetter, ':/ecology/Climate/PRISM/AN81m 1981-2015/', variable, '/', year, '/prism_', variable, '_us_30s_', year, prefix(month, 2), '.tif'))
				} else if (variable == 'et0') {
					raster(paste0(prismDriveLetter, ':/ecology/Climate/PRISM/ET0 - Monthly - Dobrowski et al 2013/et0_year', year, '_month', prefix(month, 2), '.tif'))
				} else if (variable == 'aisr') {
					raster(paste0(prismDriveLetter, ':/ecology/Climate/PRISM/Actual Incoming Solar Radiation - Dobrowski et al 2013/globalActualDownwardShortwaveRadiation_year', year, '_month', prefix(month, 2), '_doy', prefix(midMonthDoy, 3), '_MJperM2perDay.tif'))
				} else if (variable == 'swe') {
					raster(paste0(prismDriveLetter, ':/ecology/Climate/PRISM/Snow Water Equivalent - Monthly - Dobrowski et al 2013/snowWaterEquivalent_year', year, '_month', prefix(month, 2), '_mm.tif'))
				}
				
				thisMonthYear <- crop(thisMonthYear, mask)
				
				if (!exists('thisMonth', inherits=FALSE)) {
					thisMonth <- stack(thisMonthYear)
				} else {
					thisMonth <- stack(thisMonth, thisMonthYear)
				}
				
			} # next year
			
			# calculate mean for this month across years
			meanForMonth <- mean(thisMonth)
			meanForMonth <- setMinMax(meanForMonth)
			names(meanForMonth) <- variable
			
			writeRaster(meanForMonth, paste0('./Data/Climate - Monthly Means/', variable, '/', variable, '_month', prefix(month, 2), '_mean1996to2015'))
			
			rm(thisMonth, meanForMonth); gc()
			
		} # next month
		
	} # next variable
	
	dirCreate('./Data/Climate - Derived')
	
	say('ASPECT')
	
		aspect <- raster('G:/ecology/Climate/PRISM/PRISM_us_dem_800m_aspect_degrees.tif')
		aspect <- crop(aspect, mask)
		aspect <- pi * aspect / 180
		northness <- sin(aspect)
		eastness <- cos(aspect)
		names(northness) <- 'northness'
		names(eastness) <- 'eastness'
		
		writeRaster(northness, './Data/Climate - Derived/northness')
		writeRaster(eastness, './Data/Climate - Derived/eastness')
		
	say('CHRONIC HEAT')
	
		rasts <- stack(listFiles('./Data/Climate - Monthly Means/tmax', pattern='.tif'))
		rasts <- subset(rasts, 6:9)
		chronicHeat <- mean(rasts) / 100000
		chronicHeat <- setMinMax(chronicHeat)
		names(chronicHeat) <- 'chronicHeat'
		writeRaster(chronicHeat, './Data/Climate - Derived/chronicHeat')
		
	say('ACUTE COLD')
	
		rasts <- stack(listFiles('./Climate - Monthly Means/tmin', pattern='.tif'))
		acuteCold <- min(rasts)
		acuteCold <- acuteCold / 100000
		acuteCold <- setMinMax(acuteCold)
		names(acuteCold) <- 'acuteCold'
		writeRaster(acuteCold, './Data/Climate - Derived/acuteCold')
		
	say('SWE')

		rasts <- stack(listFiles('./Climate - Monthly Means/swe', pattern='.tif'))
		swe <- log10(mean(rasts) + 1)
		swe <- setMinMax(swe)
		names(swe) <- 'swe'
		writeRaster(swe, './Data/Climate - Derived/swe')
		
	say('GSAISR')

		rasts <- stack(listFiles('./Climate - Monthly Means/aisr', pattern='.tif'))
		rasts <- subset(rasts, 6:9)
		gsaisr <- sum(rasts)
		gsaisr <- setMinMax(gsaisr)
		names(gsaisr) <- 'gsaisr'
		writeRaster(gsaisr, './Data/Climate - Derived/gsaisr')
		
	say('GSWB')

		ppt <- stack(listFiles('./Climate - Monthly Means/ppt', pattern='.tif'))
		et0 <- stack(listFiles('./Climate - Monthly Means/et0', pattern='.tif'))
		ppt <- subset(ppt, 6:9)
		et0 <- subset(et0, 6:9)
		
		gswb <- ppt - et0
		
		gswb <- sum(gswb)
		gswb <- gswb / 100000
		gswb <- setMinMax(gswb)
		names(chronicHeat) <- 'gswb'
		writeRaster(gswb, './Data/Climate - Derived/gswb')
	
	say('NDVI')
		
		ndvi <- raster('./Data/NDVI 1990-2010/Raw/NDVI.tif')
		beginCluster(4)
			ndvi <- projectRaster(ndvi, mask)
		endCluster()
		names(ndvi) <- 'ndvi'
		writeRaster(ndvi, paste0('./Data/NDVI/ndvi'))

	# say('PLOT ALL PREDICTORS')
	
	# pres <- readRDS(paste0('./Presences for California Gap.rds'))
	# elevPrism <- raster(paste0('./Predictors/elevPrism_studyRegionExtent.tif'))
	
	# usa1 <- readOGR('C:/ecology/Political Geography/GADM/ver2pt8/WGS84', 'USA_adm1', verbose=FALSE)
	# usa2 <- readOGR('C:/ecology/Political Geography/GADM/ver2pt8/WGS84', 'USA_adm2', verbose=FALSE)
	
	# plumas <- usa2[usa2$NAME_2 == 'Plumas', ]
	
	# studyRegion <- readOGR(paste0('./Study Region'), 'Study Region - EPA Level III Sierra Nevada Plus 200-km Buffer', verbose=FALSE)

	# say('study region')
	# png(paste0('./Maps/Predictor Maps - Study Region.png'), width=4 * 250, height= 2 * 350, res=200)
	
		# par(mfrow=c(2, 4), mgp=c(1, 0.2, 0), mar=c(3, 2, 4, 1) + 0.1, cex.main=0.6, cex.lab=0.5, cex.axis=0.5, tck=-0.02)

		# overlays <- function() {
		
			# plot(usa2, border='gray40', lwd=0.7, add=TRUE)
			# plot(usa1, border='gray20', lwd=0.7, add=TRUE)
			# plot(plumas, border='blue', lwd=0.7, add=TRUE)
			# plot(studyRegion, border='darkgreen', lwd=0.7, add=TRUE)
			# points(pres$longWgs84, pres$latWgs84, pch=16, cex=0.2, col='red')
			
		# }
		
		# plot(elevPrism, main='Elevation'); overlays()
		# plot(northness, main='Northness'); overlays()
		# plot(chronicHeat, main='Chronic Heat'); overlays()
		# plot(acuteCold, main='Acute Cold'); overlays()
		# plot(swe, main='SWE'); overlays()
		# plot(gswb, main='Growing Season\nAtmospheric Water Balance'); overlays()
		# plot(gsaisr, main='Growing Season Actual\nIncoming Solar Radiation'); overlays()
		# plot(ndvi, main='NDVI'); overlays()
		
	# dev.off()
	
	# say('gap')
	# png(paste0('./Maps/Predictor Maps - Gap.png'), width=4 * 300, height= 2 * 300, res=200)
	
		# par(mfrow=c(2, 4), mgp=c(1, 0.2, 0), mar=c(3, 2, 4, 1) + 0.1, cex.main=0.6, cex.lab=0.5, cex.axis=0.5, tck=-0.02)

		# overlays <- function() {
		
			# plot(usa2, border='gray40', lwd=0.7, add=TRUE)
			# plot(usa1, border='gray20', lwd=0.7, add=TRUE)
			# plot(plumas, border='blue', lwd=0.7, add=TRUE)
			# plot(studyRegion, border='darkgreen', lwd=0.7, add=TRUE)
			# points(pres$longWgs84, pres$latWgs84, pch=16, cex=0.2, col='red')
			
		# }
		
		# plot(plumas, main='Elevation'); plot(elevPrism, add=TRUE); overlays()
		# plot(plumas, main='Northness'); plot(northness, add=TRUE); overlays()
		# plot(plumas, main='Chronic Heat'); plot(chronicHeat, add=TRUE); overlays()
		# plot(plumas, main='Acute Cold'); plot(acuteCold, add=TRUE); overlays()
		# plot(plumas, main='SWE'); plot(swe, add=TRUE); overlays()
		# plot(plumas, main='Growing Season\nAtmospheric Water Balance'); plot(gswb, add=TRUE); overlays()
		# plot(plumas, main='Growing Season Actual\nIncoming Solar Radiation'); plot(gsaisr, add=TRUE); overlays()
		# plot(plumas, main='NDVI'); plot(ndvi, add=TRUE); overlays()
		
	# dev.off()
	

# say('#################################')
# say('### generate background sites ###')
# say('#################################')

	# dirCreate('./Background Sites')

	# kde <- raster(paste0('./ENM - KDE/Ochotona princeps/kde/ochPri_kde_multivariate_allSites.tif'))
	# mask <- raster(paste0('./Masks/maskPrism_sierraNevadaEpaLevel3Plus200kmBuffer.tif'))
	# kde <- kde * mask
	# kde <- stretch(kde, 0, 1000)

	# bgSitesKde <- as.data.frame(randomPoints(kde, 10000, prob=TRUE))
	# bgSitesRand <- as.data.frame(randomPoints(kde, 10000))

	# names(bgSitesKde) <- c('longWgs84', 'latWgs84')
	# names(bgSitesRand) <- c('longWgs84', 'latWgs84')

	# env <- loadPredRasters()

	# bgSitesKdeEnv <- as.data.frame(extract(env, bgSitesKde))
	# bgSitesKdeRand <- as.data.frame(extract(env, bgSitesRand))

	# bgKde <- cbind(bgSitesKde, bgSitesKdeEnv)
	# bgRand <- cbind(bgSitesRand, bgSitesKdeRand)

	# saveRDS(bgKde, paste0('./Background Sites/Background Sites - Drawn from KDE.rds'))
	# saveRDS(bgRand, paste0('./Background Sites/Background Sites - Random.rds'))

say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!', level=1, pre=1)
