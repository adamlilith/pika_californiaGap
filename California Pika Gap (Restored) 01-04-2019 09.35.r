### Ochotona princeps - Spatially-varying importance of variables
### Adam B. Smith | 2017-02

# source('C:/ecology/Drive/Research/Iconic Species/Analysis - California Gap/California Pika Gap.r')

	rm(list=ls())
	memory.limit(memory.limit() * 2^30)
	gc()

### CONTENTS ###
### libraries, variables, and functions ###
### define study region ###
### collate presences for KDE ###
### train KDE and write raster ###
### train KDE, evaluate, and write rasters ###
### map gap ###
### calculate predictor layers ###
### generate background sites ###


######################
### generalization ###
######################


###########################################
### libraries, variables, and functions ###
###########################################

	setwd('C:/ecology/Drive/Research/Iconic Species/')

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

		library(omnibus)
		library(enmSdm)
		library(legendary)
		library(statisfactory)
		
		### load predictor rasters
		loadPredRasters <- function() {
			
			env <- stack(c(
				listFiles('./Analysis - California Gap/Predictors/Derived', pattern='.tif'),
				'./Analysis - California Gap/Predictors/elevPrism_studyRegionExtent.tif'
			))
			
			names(env)[nlayers(env)] <- 'elevation'
			env
			
		}

	#################
	### variables ###
	#################
	
		# doyNonLeapYear <- read.csv('C:/ecology/Drive/R/Ancillary Files/daysOfYear_nonLeapYear.csv')
		# doyLeapYear <- read.csv('C:/ecology/Drive/R/Ancillary Files/daysOfYear_leapYear.csv')

	###############
	### options ###
	###############
	
		options(stringsAsFactors=FALSE)
		rasterOptions(format='GTiff', overwrite=TRUE)

say('###########################')
say('### define study region ###')
say('###########################')
	
	# study region extent derived from modified EPA Level III "Sierra Nevada" ecoregion plus 200-km buffer
	dirCreate('./Analysis - California Gap/Study Region')
	dirCreate('./Analysis - California Gap/Predictors')

	epa3 <- shapefile('C:/Ecology/Drive/Research/Iconic Species/Extents_Masks_Maps/EcoRegions/!OmernikEcoregions/!Pika_Omernik3REV_Join')

	sn <- epa3[epa3$L3_KEY == '5  Sierra Nevada', ]
	sn <- sp::spTransform(sn, getCRS('climateNA'))
	snBuff <- gBuffer(sn, width=200000)

	snBuff <- sp::spTransform(snBuff, CRS(crsNad83))
	sn <- sp::spTransform(sn, CRS(crsNad83))

	writeOGR(sn, paste0('./Analysis - California Gap/Study Region'), 'Study Region - EPA Level III Sierra Nevada Plus 200-km Buffer', driver='ESRI Shapefile', overwrite_layer=TRUE, verbose=FALSE)
	
	# mask PRISM elevation
	elevPrism <- raster('G:/ecology/Climate/PRISM/PRISM_us_dem_800m.tif')
	elevPrism <- crop(elevPrism, snBuff)
	
	writeRaster(elevPrism, paste0('./Analysis - California Gap/Predictors/elevPrism_studyRegionExtent'), format='GTiff', overwrite=TRUE)

	# mask from PRISM for extent of study region
	dirCreate('./Analysis - California Gap/Masks')
	
	maskPrism <- elevPrism * 0
	maskPrism <- rasterize(snBuff, maskPrism)
	writeRaster(maskPrism, paste0('./Analysis - California Gap/Masks/maskPrism_sierraNevadaEpaLevel3Plus200kmBuffer'), format='GTiff', overwrite=TRUE)

	# mask GMTED2010 elevation
	elevGmted2010 <- raster('G:/ecology/Topography/GMTED2010/7pt5_arcsec/conus_20101117_gmted_mea075.tif')
	elevGmted2010 <- crop(elevGmted2010, snBuff)
	
	writeRaster(elevGmted2010, paste0('./Analysis - California Gap/Predictors/elevGmted2010_studyRegionExtent'), format='GTiff', overwrite=TRUE)

	# mask from GMTED2010
	dirCreate('./Analysis - California Gap/Masks')
	
	maskGmted2010 <- elevGmted2010 * 0
	writeRaster(maskGmted2010, paste0('./Analysis - California Gap/Masks/maskGmted2010_studyRegionExtent'), format='GTiff', overwrite=TRUE)
	maskGmted2010 <- rasterize(snBuff, maskGmted2010)
	
	writeRaster(maskGmted2010, paste0('./Analysis - California Gap/Masks/maskGmted2010_sierraNevadaEpaLevel3'), format='GTiff', overwrite=TRUE)
	
# say('#################################')
# say('### collate presences for KDE ###')
# say('#################################')

	# pres <- readRDS(paste0('./Species Records - Pika/!Collated Data 2016-06-30 1256/02 Ochotona princeps - Usable - Presences & Non-detections 1895-2015 - PRISM & DayMet Climate Data Extracted.rds'))

	# pres <- pres[ , -which(grepl(names(pres), pattern='Daily'))]
	# pres <- pres[ , -which(grepl(names(pres), pattern='Monthly'))]
	# pres <- pres[ , -which(grepl(names(pres), pattern='Annual'))]
	
	# ## subset to presences AND years 1981-2015 AND providers who agreed to allow their data to be used and publicly-available databases
	# pres <- pres[
		# pres$origRecentPikaOrSignsObserved & pres$obsYear >= 1981 & pres$obsYear <= 2015 & (
			# pres$contact == 'NA' |
			# pres$contact == 'Adam Smith;' |
			# pres$contact == 'Chris Ray (Colorado University, Boulder)' |
			# pres$contact == 'Erik Beever, Chris Ray' |
			# pres$contact == 'Erik Beever' |
			# pres$contact == 'Erik Beever;' |
			# pres$contact == 'Ken Goehring;' |
			# pres$contact == 'Jason Brewer & Mary Flores' |
			# pres$contact == 'Clint Epps, Jessica Castillo' |
			# pres$contact == 'Jess Castillo, Clint Epps' |
			# pres$contact == 'Connie Millar'
		# ), ]

	# ## subset to presences in Sierra Nevada EPALevel III ecoregion
	# epa3 <- readOGR(paste0('./Extents_Masks_Maps/EcoRegions/OmernikEcoregions'), 'us_eco_l3SarrREV', verbose=FALSE)
	# epa3 <- sp::spTransform(epa3, CRS(crsNad83))
	
	# epa4 <- readRDS('C:/ecology/Ecoregionalizations/Ecoregions - Level IV - EPA/Eco_Level_IV_US_wgs84_dissolved.rds')
	# studyRegion <- readOGR(paste0('./Analysis - California Gap/Study Region'), 'Study Region - EPA Level III Sierra Nevada Plus 200-km Buffer', verbose=FALSE)
	
	# pres$epa3 <- raster::extract(epa3, cbind(pres$longWgs84, pres$latWgs84))$L3_KEY
	# pres$epa4 <- raster::extract(epa4, cbind(pres$longWgs84, pres$latWgs84))$US_L4NAME
	# pres$studyRegion <- extract(studyRegion, cbind(pres$longWgs84, pres$latWgs84))$L3_KEY
	
	# pres <- pres[!is.na(pres$studyRegion), ] # just in Sierra Nevada + buffer study region
	# pres$lassen <- pres$latWgs84 > 40 & pres$latWgs84 < 41 # flag for Lassen and surrounding area

	# ## save
	# saveRDS(pres, paste0('.//Analysis - California Gap/Presences for California Gap.rds'))

	# collateData(
		# speciesList='Ochotona princeps',
		# outputDir=paste0('./Analysis - California Gap/ENM - KDE'),
		# allSpeciesData=pres[ , 1:which(names(pres) %in% 'cellAreaPrism_km2')],
		# allTrainAbsenceEnvData=NULL,
		# allTestAbsenceEnvData=NULL,
		# maxNumTrainingAbsences=10000,
		# kFolds=10,
		# split='subsample',
		# userSplitField=NULL,
		# trainingProportion=0.80,
		# envThinPres=FALSE,
		# envThinAbs=FALSE,
		# envThinPerc=c(0.1, 0.2),
		# spatialThinPres=FALSE,
		# spatialThinTrainAbs=FALSE,
		# spatialThinTestAbs=FALSE,
		# minDist=100000,
		# distFunct=NULL,
		# minNumTrainSites=50,
		# minNumTestSites=10,
		# splitAbs=FALSE,
		# responseVar=NULL,
		# responseAsFactors=NULL,
		# predictorsToUse=NULL,
		# predictorsAsFactors=NULL,
		# redundRastPres=NULL,
		# redundRastAbsTrain=NULL,
		# redundRastAbsTest=NULL,
		# longLatFields=c('longWgs84', 'latWgs84'),
		# CRS=crsNad83,
		# speciesField='species',
		# speciesDataFileAppend=NULL,
		# verbose=1
	# )
	
# say('##############################################')
# say('### train KDE, evaluate, and write rasters ###')
# say('##############################################')

	# templateRast <- raster(paste0('./Analysis - California Gap/Masks/maskPrism_sierraNevadaEpaLevel3Plus200kmBuffer.tif'))

	# say('model...')
	# trainSdm(
		# speciesList='Ochotona princeps',
		# outputDir=paste0('./Analysis - California Gap/ENM - KDE'),
		# predictorsToUse=c('longWgs84', 'latWgs84'),
		# modelType='multivariate',
		# modelMethodsToUse='kde',
		# kFolds=TRUE,
		# outNames='short',
		# modelParams=list(
			# kde=list(
				# h='epa',
				# kern='epa',
				# templateRast=templateRast,
				# templateProj4=crsNad83,
				# eaProj4=crsClimateNA,
				# resMult=4,
				# epsilon=100
			# )
		# ),
		# tempDir='C:/ecology/!Scratch',
		# verbose=0
	# )

	# say('write rasters...')
	# writeModelRasters(
		# speciesList='Ochotona princeps',
		# baseModelDir=paste0('./Analysis - California Gap/ENM - KDE'),
		# rasterOutputDirName='kde',
		# predictorStack=longLatRasters(template=templateRast),
		# modelMethodsToUse='kde',
		# dataType=c('k', 'all sites'),
		# modelType='multivariate',
		# rescale=NULL,
		# outNames='short',
		# scratchDir='C:/ecology/!Scratch',
		# verbose=0,
		# format='GTiff',
		# overwrite=TRUE
	# )	
	
	# say('evaluate...')
	# testBg <- as.data.frame(randomPoints(templateRast, 10000))
	# names(testBg) <- c('longWgs84', 'latWgs84')
	
	# tHoldMaxSSS <- tHold10PercTest <- tHoldMinPred <- rep(NA, 10)
	
	# for (k in 1:10) {

		# say(k)
	
		# load(paste0('./Analysis - California Gap/ENM - KDE/Ochotona princeps/KDE multivariate k=', prefix(k, 2), '.Rdata'))
		
		# kde <- raster(paste0('./Analysis - California Gap/ENM - KDE/Ochotona princeps/kde/ochPri_kde_multivariate_k', prefix(k, 2), '.tif'))
		
		# kde <- stretch(kde, 0, 1)
		
		# predBg <- extract(kde, testBg)
		# predPres <- extract(kde, model$testPresences[ , c('longWgs84', 'latWgs84')])

		# thisEval <- evaluate(p=as.vector(predPres), a=as.vector(predBg), tr=seq(0, 1, by=10E-5))
		
		# tHoldMaxSSS[k] <- threshold(thisEval, stat='spec_sens')
		# tHold10PercTest[k] <- quantile(predPres, 0.1)
		# tHoldMinPred[k] <- min(predPres)
		
	# }
	
	# tHoldMaxSSS <- median(tHoldMaxSSS)
	# tHold10PercTest <- median(tHold10PercTest)
	# tHoldMinPred <- median(tHoldMinPred)
	
	# tHolds <- data.frame(tHoldMaxSSS=tHoldMaxSSS, tHold10PercTest=tHold10PercTest, tHoldMinPred=tHoldMinPred)
	# write.csv(tHolds, paste0('./Analysis - California Gap/KDE Thresholds.csv'), row.names=FALSE) 
	
# say('###############')
# say('### map gap ###')
# say('###############')

	# dirCreate('./Analysis - California Gap/Maps')

	# tHolds <- read.csv(paste0('./Analysis - California Gap/KDE Thresholds.csv'))

	# usa1 <- readOGR('C:/ecology/Political Geography/GADM/ver2pt8/WGS84', 'USA_adm1', verbose=FALSE)
	# usa2 <- readOGR('C:/ecology/Political Geography/GADM/ver2pt8/WGS84', 'USA_adm2', verbose=FALSE)
	
	# plumas <- usa2[usa2$NAME_2 == 'Plumas', ]
	
	# elevPrism <- raster('G:/ecology/Climate/PRISM/PRISM_us_dem_800m.tif')
	# elevGmted2010 <- raster(paste0('./Analysis - California Gap/Predictors/elevGmted2010_studyRegionExtent.tif'))

	# epa3 <- readOGR(paste0('./Extents_Masks_Maps/EcoRegions/OmernikEcoregions'), 'us_eco_l3SarrREV', verbose=FALSE)
	# epa3 <- sp::spTransform(epa3, CRS(crsNad83))
	
	# studyRegion <- readOGR(paste0('./Analysis - California Gap/Study Region'), 'Study Region - EPA Level III Sierra Nevada Plus 200-km Buffer', verbose=FALSE)

	# kde <- raster(paste0('./Analysis - California Gap/ENM - KDE/Ochotona princeps/kde/ochPri_kde_multivariate_allSites.tif'))
	# kde <- stretch(kde, 0, 1)
	
	# kdeMinPred <- kde < tHolds$tHoldMinPred
	# kdeMinPred <- calc(kdeMinPred, fun=function(x) ifelse(x == 1, 1, NA))
	
	# pres <- readRDS(paste0('./Analysis - California Gap/Presences for California Gap.rds'))
	
	# presOkToShowErik <- pres[pres$contact=='NA' | pres$contact=='Adam Smith;' | pres$contact=='Erik Beever, Chris Ray' | pres$contact=='Erik Beever' | pres$contact=='Erik Beever;' | pres$contact=='NA', ]
	
	# # plot gap
	# png(paste0('./Analysis - California Gap/Maps/KDE - Centered on Plumas County.png'), width=1600, height=1600, res=300)
	
		# plumasBuffSmall <- gBuffer(sp::spTransform(plumas, CRS(crsClimateNA)), width=20000)
		# plumasBuffSmall <- sp::spTransform(plumasBuffSmall, CRS(crsNad83))
		
		# plumasBuffLarge <- gBuffer(sp::spTransform(plumas, CRS(crsClimateNA)), width=30000)
		# plumasBuffLarge <- sp::spTransform(plumasBuffLarge, CRS(crsNad83))
		
		# elevGmted2010Crop <- crop(elevGmted2010, plumasBuffLarge)
		
		# plot(plumasBuffSmall, border=NA)
		# plot(elevGmted2010Crop, add=TRUE, breaks=seq(minValue(elevGmted2010Crop), maxValue(elevGmted2010Crop), length.out=100), col=paste0('gray', 0:100), legend=FALSE)
		# plot(epa3[epa3$L3_KEY == '5  Sierra Nevada', ], add=TRUE)
		# plot(kdeMinPred, add=TRUE, col=alpha('blue', 0.3), legend=FALSE)
		# plot(usa2[usa2$NAME_1 == 'California' | usa2$NAME_1 == 'Nevada', ], lwd=0.8, add=TRUE, border='orange')
		# # points(pres$longWgs84, pres$latWgs84, col='white', pch=1)
		# points(presOkToShowErik$longWgs84, presOkToShowErik$latWgs84, col='red', pch=1)
		
		# legend('topright', inset=0.01, legend=c('Sierra Nevada', 'Counties', 'Gap', 'Erik\'s/Public Presence'), col=c('black', 'orange', NA, 'red'), pch=c(NA, NA, NA, 1), lwd=c(1, 1, NA, NA), fill=c(NA, NA, alpha('blue', 0.5), NA), border=c(NA, NA, 'blue', NA), bg='white', cex=0.6, title='California Pika Gap')
		
		# title('Kernel Density Estimator of California Gap', sub=date(), cex.main=0.9, cex.sub=0.4)
		
	# dev.off()
	
# say('##################################')
# say('### calculate predictor layers ###')
# say('##################################')

	# dirCreate('./Analysis - California Gap/Predictors/Derived')
	
	# # calculate time period for climate layers
	# pres <- readRDS(paste0('./Analysis - California Gap/Presences for California Gap.rds'))
	
	# png(paste0('./Analysis - California Gap/Year of Observations.png'))
		# hist(pres$obsYear, breaks=(min(pres$obsYear) - 1):(max(pres$obsYear) + 1), main='Observations')
	# dev.off()

	# say('Using 20-yr window for climate layers (1996-2015).')
	
	# studyRegionMask <- raster(paste0('./Analysis - California Gap/Masks/maskPrism_sierraNevadaEpaLevel3Plus200kmBuffer.tif'))

	# say('MEAN MONTHLY PPT, TMIN, TMAX, ET0, AISR, and SWE')

	# # calculate monthly means of max/min temperature, et0, actual incoming solar radiation, and snow water equivalent
	# for (variable in c('ppt', 'tmin', 'tmax', 'et0', 'aisr', 'swe')) {
		
		# dirCreate('./Analysis - California Gap/Predictors/Monthly Means/', variable)

		# for (month in 1:12) {
		
			# say(variable, ' | month ', month, ' | year', post=0)
		
			# for (year in 1996:2015) {
		
				# say(year, post=ifelse(year == 2015, 1, 0))

				# # get mid-say of month (for AISR)
				# doyThisYear <- if (year %% 4 == 0) {
					# doyLeapYear
				# } else {
					# doyNonLeapYear
				# }
			
				# # mid-day of this month
				# midMonthDoy <- doyNonLeapYear[15, paste0('month', month)]
			
				# thisMonthYear <- if (variable %in% c('ppt', 'tmin', 'tmax')) {
					# raster(paste0('G:/ecology/Climate/PRISM/AN81m 1981-2015/', variable, '/', year, '/prism_', variable, '_us_30s_', year, prefix(month, 2), '.tif'))
				# } else if (variable == 'et0') {
					# raster(paste0('G:/ecology/Climate/PRISM/ET0 - Monthly - Dobrowski et al 2013/et0_year', year, '_month', prefix(month, 2), '.tif'))
				# } else if (variable == 'aisr') {
					# raster(paste0('G:/ecology/Climate/PRISM/Actual Incoming Solar Radiation - Dobrowski et al 2013/globalActualDownwardShortwaveRadiation_year', year, '_month', prefix(month, 2), '_doy', prefix(midMonthDoy, 3), '_MJperM2perDay.tif'))
				# } else if (variable == 'swe') {
					# raster(paste0('G:/ecology/Climate/PRISM/Snow Water Equivalent - Monthly - Dobrowski et al 2013/snowWaterEquivalent_year', year, '_month', prefix(month, 2), '_mm.tif'))
				# }
				
				# thisMonthYear <- crop(thisMonthYear, studyRegionMask)
				
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
			
			# writeRaster(meanForMonth, paste0('./Analysis - California Gap/Predictors/Monthly Means/', variable, '/', variable, '_month', prefix(month, 2), '_mean1996to2015'), format='GTiff', overwrite=TRUE)
			
			# rm(thisMonth, meanForMonth); gc()
			
		# } # next month
		
	# } # next variable
	
	# say('ASPECT')
	# aspect <- raster('G:/ecology/Climate/PRISM/PRISM_us_dem_800m_aspect_degrees.tif')
	# aspect <- crop(aspect, studyRegionMask)
	# aspect <- pi * aspect / 180
	# northness <- sin(aspect)
	# names(northness) <- 'northness'
	
	# writeRaster(northness, paste0('./Analysis - California Gap/Predictors/Derived/northness'), format='GTiff', overwrite=TRUE)
	
	# say('CHRONIC HEAT')
	# rasts <- stack(listFiles(paste0('./Analysis - California Gap/Predictors/Monthly Means/tmax'), pattern='.tif'))
	# rasts <- subset(rasts, 6:9)
	# chronicHeat <- mean(rasts) / 100000
	# chronicHeat <- setMinMax(chronicHeat)
	# names(chronicHeat) <- 'chronicHeat'
	# writeRaster(chronicHeat, paste0('./Analysis - California Gap/Predictors/Derived/chronicHeat'), format='GTiff', overwrite=TRUE)
	
	# say('ACUTE COLD')
	# rasts <- stack(listFiles(paste0('./Analysis - California Gap/Predictors/Monthly Means/tmin'), pattern='.tif'))
	# acuteCold <- min(rasts)
	# acuteCold <- acuteCold / 100000
	# acuteCold <- setMinMax(acuteCold)
	# names(acuteCold) <- 'acuteCold'
	# writeRaster(acuteCold, paste0('./Analysis - California Gap/Predictors/Derived/acuteCold'), format='GTiff', overwrite=TRUE)
	
	# say('SWE')
	# rasts <- stack(listFiles(paste0('./Analysis - California Gap/Predictors/Monthly Means/swe'), pattern='.tif'))
	# swe <- log10(mean(rasts) + 1)
	# swe <- setMinMax(swe)
	# names(swe) <- 'swe'
	# writeRaster(swe, paste0('./Analysis - California Gap/Predictors/Derived/swe'), format='GTiff', overwrite=TRUE)
	
	# say('GSAISR')
	# rasts <- stack(listFiles(paste0('./Analysis - California Gap/Predictors/Monthly Means/aisr'), pattern='.tif'))
	# rasts <- subset(rasts, 6:9)
	# gsaisr <- sum(rasts)
	# gsaisr <- setMinMax(gsaisr)
	# names(gsaisr) <- 'gsaisr'
	# writeRaster(gsaisr, paste0('./Analysis - California Gap/Predictors/Derived/gsaisr'), format='GTiff', overwrite=TRUE)
	
	# say('GSWB')
	# ppt <- stack(listFiles(paste0('./Analysis - California Gap/Predictors/Monthly Means/ppt'), pattern='.tif'))
	# et0 <- stack(listFiles(paste0('./Analysis - California Gap/Predictors/Monthly Means/et0'), pattern='.tif'))
	# ppt <- subset(ppt, 6:9)
	# et0 <- subset(et0, 6:9)
	
	# gswb <- ppt - et0
	
	# gswb <- sum(gswb)
	# gswb <- gswb / 100000
	# gswb <- setMinMax(gswb)
	# names(chronicHeat) <- 'gswb'
	# writeRaster(gswb, paste0('./Analysis - California Gap/Predictors/Derived/gswb'), format='GTiff', overwrite=TRUE)
	
	# say('NDVI')
	# ndvi <- raster(paste0('./Analysis - California Gap/Predictors/NDVI 1990-2010 Pre-processed/NDVI.tif'))
	# ndvi <- projectRaster(ndvi, studyRegionMask)
	# writeRaster(ndvi, paste0('./Analysis - California Gap/Predictors/Derived/ndvi'), format='GTiff', overwrite=TRUE)

	# say('PLOT ALL PREDICTORS')
	
	# pres <- readRDS(paste0('./Analysis - California Gap/Presences for California Gap.rds'))
	# elevPrism <- raster(paste0('./Analysis - California Gap/Predictors/elevPrism_studyRegionExtent.tif'))
	
	# usa1 <- readOGR('C:/ecology/Political Geography/GADM/ver2pt8/WGS84', 'USA_adm1', verbose=FALSE)
	# usa2 <- readOGR('C:/ecology/Political Geography/GADM/ver2pt8/WGS84', 'USA_adm2', verbose=FALSE)
	
	# plumas <- usa2[usa2$NAME_2 == 'Plumas', ]
	
	# studyRegion <- readOGR(paste0('./Analysis - California Gap/Study Region'), 'Study Region - EPA Level III Sierra Nevada Plus 200-km Buffer', verbose=FALSE)

	# say('study region')
	# png(paste0('./Analysis - California Gap/Maps/Predictor Maps - Study Region.png'), width=4 * 250, height= 2 * 350, res=200)
	
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
	# png(paste0('./Analysis - California Gap/Maps/Predictor Maps - Gap.png'), width=4 * 300, height= 2 * 300, res=200)
	
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

	# dirCreate('./Analysis - California Gap/Background Sites')

	# kde <- raster(paste0('./Analysis - California Gap/ENM - KDE/Ochotona princeps/kde/ochPri_kde_multivariate_allSites.tif'))
	# mask <- raster(paste0('./Analysis - California Gap/Masks/maskPrism_sierraNevadaEpaLevel3Plus200kmBuffer.tif'))
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

	# saveRDS(bgKde, paste0('./Analysis - California Gap/Background Sites/Background Sites - Drawn from KDE.rds'))
	# saveRDS(bgRand, paste0('./Analysis - California Gap/Background Sites/Background Sites - Random.rds'))

say('###############################################')
say('### collate detections and background sites ###')
say('###############################################')

	bg <- readRDS(paste0('./Analysis - California Gap/Background Sites/Background Sites - Drawn from KDE.rds'))

	pres <- readRDS(paste0('./Analysis - California Gap/Presences for California Gap.rds'))
	pres <- pres[!pres$lassen, ]
	env <- loadPredRasters()
	presEnv <- as.data.frame(extract(env, cbind(pres$longWgs84, pres$latWgs84)))
	pres <- cbind(pres, presEnv)

	collateData(
		speciesList='Ochotona princeps',
		outputDir=paste0('./Analysis - California Gap/ENM - Sans Gap & Lassen'),
		allSpeciesData=pres,
		allTrainAbsenceEnvData=bg,
		allTestAbsenceEnvData=bg,
		maxNumTrainingAbsences=10000,
		kFolds=5,
		split='geographic',
		userSplitField=NULL,
		trainingProportion=0.80,
		envThinPres=FALSE,
		envThinAbs=FALSE,
		envThinPerc=c(0.1, 0.2),
		spatialThinPres=FALSE,
		spatialThinTrainAbs=FALSE,
		spatialThinTestAbs=FALSE,
		minDist=100000,
		distFunct=NULL,
		minNumTrainSites=30,
		minNumTestSites=6,
		splitAbs=TRUE,
		responseAsFactors=NULL,
		predictorsToUse=names(env),
		predictorsAsFactors=NULL,
		redundRastPres=NULL,
		redundRastAbsTrain=NULL,
		redundRastAbsTest=NULL,
		longLatFields=c('longWgs84', 'latWgs84'),
		CRS=crsNad83,
		speciesField='species',
		speciesDataFileAppend=NULL,
		verbose=1
	)

say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!', pre=1)
say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
