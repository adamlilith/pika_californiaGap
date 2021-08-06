### Ochotona princeps - Calculate evapotranspiration rasters
### Adam B. Smith | 2016-08
###
### Calculate evapotranspiration rasters.

# source('F:/ecology/Climate/PRISM/Calculate Incoming Solar Radiation, SWE, and ET0 from PRISM, DayMet, and NLDAS.r')

	memory.limit(memory.limit() * 2^30)
	rm(list=ls())
	options(stringsAsFactors=FALSE)
	gc()

### CONTENTS ###
### libraries, variables, and functions ###
### calculate MONTHLY actual incoming solar radiation (see Appendix of Dobrowski et al 2013) ###
### calculate MONTHLY snow snow water equivalent (see Appendix of Dobrowski et al 2013) ###
### calculate mean annual SWE from monthly SWE (Dobrowski et al 2013) ###
### calculate DAILY ET0 (see Appendix of Dobrowski et al 2013) ###
### calculate MONTHLY ET0 (see Appendix of Dobrowski et al 2013) ###
### calculate ANNUAL ET0 from monthly ET0 (see Appendix of Dobrowski et al 2013) ###

###########################################
### libraries, variables, and functions ###
###########################################

	library(omnibus)
	library(enmSdm)
	library(rgdal)
	library(raster)
	library(rgrass7)
	# grassDir <- c('C:/OSGeo4W64/', 'grass-7.4.1', 'osgeo4W')
	grassDir <- c('C:/OSGeo4W64/', 'grass-78', 'osgeo4W')

	# days of year
	data(doyLeap, package='omnibus')
	data(doyNonLeap, package='omnibus')
	
	### functions from Dobrowski et al 2013 (modified) ###

	# This script contains 4 functions used to model ET0 and water balance:
	# 1. 'snowmod' estimates snowfall and snowpack and net moisture input as a function of temperature, precip, and existing
	#   snowpack.  It also outputs a vector of albedo values, generally 0.2 if there is no snow, or 0.8 if there is snow.
	# 2. 'monthlyET0' for calculating monthly reference evapotranspiration
	# 3. 'dailyET0' for calculating daily reference evapotranspiration
	# 4. 'aetmod' estimates actual et, deficit, soil moisture and runoff as a function of moisture input, existing
	#   soil moisture, and soil water capacity. 
	#
	# Author: Alan Swanson 2012
	###############################################################################

	snowmodModified <- function(tmean, ppt, radiation=NULL, snowpack_prev=NULL, albedo=0.23, albedo_snow=0.8){
		# This function computes monthly estimated snowfall and snowmelt. Output includes end-of-month snowpack,
		# water "input" (snowmelt plus rain), and albedo.
		# Arguments:
		#  tmean - vector of mean monthly temperatures
		#  radiation - vector of shortwave solar radiation in MJ/m^2/day.
		#  snowpack_prev - vector of snowpack at the beginning of the month.  If NULL this is 
		#    taken to be zero.
	    #  albedo - a single value for albedo in the absence of snow cover.  
	    #  albedo_snow - single value of albedo given snow cover
		#
		# Value:  dataframe with three columns for end-of-month snowpack, H2O input (rain plus snowmelt),
		#  and albedo.
		
		
	  N <- length(tmean)
		if(is.null(snowpack_prev)) snowpack_prev <- rep(0,N)
		snowpack <- rep(NA,N)
		input <- rep(NA,N)
	
	  # this is for radiation in MJ/m^2/day
		mf <- function(t,t0,t1) pmin(pmax((t-t0)/(t1-t0),0),1)
		linrmelt <- function(temp,radiation,b0,b1,b2) pmax((b0+temp*b1+radiation*b2),0)
		parvec <- c(-4.604,6.329,-398.4,81.75,25.05)
		mfsnow <- mf(tmean,parvec[1],parvec[2])
		mfmelt <- linrmelt(tmean,radiation,parvec[3],parvec[4],parvec[5])

	  # calculate values
		snow <- (1-mfsnow)*ppt
		rain <- mfsnow*ppt	
		melt <- pmin(mfmelt,snow+snowpack_prev) 
		snowpack <- snowpack_prev+snow-melt 
		input <-rain+melt
		
		# make vector of albedo values
		albedo <- rep(albedo,N)
		albedo[snowpack>0 | (snowpack_prev>0)] <- albedo_snow
		
		# return(data.frame(snowpack=snowpack,input=input,albedo=albedo))
		return(snowpack)
	}

	monthlyET0 <- function(radiation, tmax, tmin, wind, lat, elev, dpt, tmean_prev, albedo=0.23, month){
	# This function runs Reference ET estimates for monthly timesteps using methods based on
	# the Penman-Montieth equation as presented in Allen et al (1998).  It incorporates a 
	# modification which adjusts stomatal conductance downwards at temperatures below 5 C.
	#
	# Arguments:
	# radiation: vector of monthly average shortwave radiation in MJ/m^2/day
	# tmax, tmin: vectors of monthly average maximum and minimum temperatures in C, 
	# wind: vector of monthly average wind speed in m/s at 10m above ground, 
	# tmean_prev: vector of mean temp for the previous month, 
	# lat: vector of latitude in degrees 
	# elev: vector of elevation in meters, 
	# dpt: vector of dewpoint temperature in C.
	# tmean_prev: vector of mean temp of previous month in C
	# albedo: vector or scalar of albedo values [Appendix B says modified to 0.8 if snow present], 
	# month: scalar 1-12.

	#
	# Value: 
	# Returns a vector of ET0 values.

		t0 <- unclass(Sys.time())
		daysinmonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
		d2 <- c(31,59,90,120,151,181,212,243,273,304,334,365)
		d1 <- c(1,32,60,91,121,152,182,213,244,274,305,335)
		DoY <- (d1[month]+d2[month])/2 # use middle day of month to represent monthly average. 
		n_days <- daysinmonth[month]
		
		# calculate soil heat flux (total for the month) using change in temperature from previous month
		tmean <- (tmax+tmin)/2 
		G <- 0.14*(tmean-tmean_prev) # fixed from previous version
		
	  # convert to wind height at 2m
	  hw=10 # height of wind measurements 
	  wind <- wind*(4.87/log(67*hw-5.42))  # convert to wind height at 2m
		
	  # stomatal conductance adjustment for low temperatures
	  sr=100 # stomatal resistance sec/m
		ks_min=0.01 # minimum value for temps below T1
		Tl=-10       # minimum temp (sc goes to ks_min below this temp)
		T0=5		# optimal temp
		Th=100     # maximum temp (sc goes to zero above this)
		thresh=5   # temperature threshold below which to apply Jarvis equation (ks=1 above this temp)
		b4 <- (Th-T0)/(Th-Tl)
		b3 <- 1/((T0-Tl)*(Th-T0)^b4)
		ks <- pmax(pmin(b3*(tmean-Tl)*(Th-tmean)^b4,1),ks_min)
		ks[is.na(ks)] <- ks_min
		ks[tmean>=thresh] <- 1
			
		# convert to stomatal resistance.
		sr  <- sr/ks

		# ra is aerodynamic resistance, rs is bulk surface resistance
		# ra  <- 208/wind #(log((2-2/3*0.12)/(0.123*0.12))*log((2-2/3*0.12)/(0.1*0.123*0.12)))/(0.41^2*wind) # equal to 208/wind for hh=hw=2.
		ra <- (log((2-2/3*0.12)/(0.123*0.12))*log((2-2/3*0.12)/(0.1*0.123*0.12)))/(0.41^2*wind)
		rs <- sr/(0.5*24*0.12) # value of 70 when sr=100
		
		# Saturation vapor pressure
		es <- 0.6108*exp(tmin*17.27/(tmin+237.3))/2+0.6108*exp(tmax*17.27/(tmax+237.3))/2     
		ea <- 0.6108*exp((dpt)*17.27/((dpt)+237.3))
		vpd <- es - ea
		vpd[vpd<0] <- 0    # added because this can be negative if dewpoint temperature is greater than mean temp (implying vapor pressure greater than saturation).
		
		# delta - Slope of the saturation vapor pressure vs. air temperature curve at the average hourly air temperature 
		delta  <- (4098 * es)/(tmean + 237.3)^2  

		P <- 101.3*((293-0.0065*elev)/293)^5.26  # Barometric pressure in kPa
		lambda <- 2.501-2.361e-3*tmean # latent heat of vaporization    
		cp  <- 1.013*10^-3 # specific heat of air
		gamma <- cp*P/(0.622*lambda) # Psychrometer constant (kPa C-1)
		pa <- P/(1.01*(tmean+273)*0.287) # mean air density at constant pressure
		
		# Calculate potential max solar radiation or clear sky radiation	
		GSC=0.082      # MJ m -2 min-1 (solar constant)
		phi <- pi*lat/180 
		dr <- 1+0.033*cos(2*pi/365*DoY)      
		delt <- 0.409*sin(2*pi/365*DoY-1.39)     
		omegas <- acos(-tan(phi)*tan(delt)) 
		Ra <- 24*60/pi*GSC*dr*(omegas*sin(phi)*sin(delt) +cos(phi)*cos(delt)*sin(omegas))    # Daily extraterrestrial radiation
		Rso <- Ra*(0.75+2e-5*elev)     #For a cloudless day, Rs is roughly 75% of extraterrestrial radiation (Ra)
		
		
		# radfraction is a measure of relative shortwave radiation, or of
		# possible radiation (cloudy vs. clear-sky)
		radfraction <- radiation/Rso
		radfraction[radfraction>1] <- 1
		
		# longwave  and net radiation
		longw <- 4.903e-9*n_days*((tmax+273.15)^4+(tmin+273.15)^4)/2*(.34-.14*sqrt(ea))*(1.35*radfraction-.35)     
		netrad <- radiation*n_days*(1-albedo)-longw     
		
		# ET0
		et0 <- 0.408*((delta*(netrad-G))+(pa*cp*vpd/ra*3600*24*n_days))/(delta+gamma*(1+rs/ra))
		return(et0)
	} 

	dailyET0 <- function(radiation, tmax, tmin, wind, lat, elev, albedo=0.23, dpt, doy){
		# This is a ET0 function designed for daily inputs.  
	  
	  # Arguments:
	  # radiation: vector of monthly average shortwave radiation in MJ/m^2/day
	  # tmax, tmin: vectors of monthly average maximum and minimum temperatures in C, 
	  # wind: vector of monthly average wind speed in m/s, 
	  # tmean_prev: vector of mean temp for the previous month, 
	  # lat: vector of latitude in degrees 
	  # elev: vector of elevation in meters, 
	  # albedo: scaler or vector of albedo values, 
	  # doy: scalar day of year 1-365,
	  
	  #
	  # Value: 
	  # Returns a vector of ET0 values.
		tmean <- (tmin+tmax)/2
		n_days <- 1 
		G <- 0 # assume soil heat flux to be zero
		
		# wind adjustment to 2m from 10m output
		hw=10 # height at which wind is measured
        wind <- wind*(4.87/log(67*hw-5.42))  
	  
	    # stomatal conductance adjustment for low temperatures
	    sr=100 # stomatal resistance sec/m
		ks_min=.01 # minimum value for temps below T1
		Tl=-10       # minimum temp (sc goes to ks_min below this temp)
		T0=5		# optimal temp
		Th=100     # maximum temp (sc goes to zero above this)
		thresh=5   # temperature threshold below which to apply Jarvis equation (ks=1 above this temp)
		b4=(Th-T0)/(Th-Tl)  # from Jarvis 1978
		b3=1/((T0-Tl)*(Th-T0)^b4)
		ks=pmax(pmin(b3*(tmean-Tl)*(Th-tmean)^b4,1),ks_min)
		ks[is.na(ks)] <- ks_min
		ks[tmean>=thresh] <- 1
			
		# convert to stomatal resistance.
		sr <- sr/ks
			
		# ra is aerodynamic resistance, rs is bulk surface resistance
		ra <- 208/wind # 
		rs <- sr/(0.5*24*0.12) # value of 70 when sr=100
		
		# Saturation vapor pressure , 
		es <- 0.6108*exp(tmin*17.27/(tmin+237.3))/2+0.6108*exp(tmax*17.27/(tmax+237.3))/2  # saturation vapor pressure
		ea <- 0.6108*exp((dpt)*17.27/((dpt)+237.3))                                        # actual vapor pressure
		vpd <- es - ea
		vpd[vpd<0] <- 0    
	  
		# delta - Slope of the saturation vapor pressure vs. air temperature curve
		delta <- (4098 * es)/(tmean + 237.3)^2  
		P <- 101.3*((293-0.0065*elev)/293)^5.26  # Barometric pressure in kPa
		lambda <- 2.501-2.361e-3*tmean # latent heat of vaporization    
		cp <- 1.013*10^-3 # specific heat of air
		gamma <- cp*P/(0.622*lambda) # Psychrometer constant (kPa C-1)
		pa <- P/(1.01*(tmean+273)*0.287) # mean air density at constant pressure
		
		# Calculate potential max solar radiation or clear sky radiation.
		GSC = 0.082      # MJ m -2 min-1 (solar constant)
		phi <- pi*lat/180 
		dr <- 1+0.033*cos(2*pi/365*doy)      
		delt <- 0.409*sin(2*pi/365*doy-1.39)     
		omegas <- acos(-tan(phi)*tan(delt))     
		Ra <- 24*60/pi*GSC*dr*(omegas*sin(phi)*sin(delt) +cos(phi)*cos(delt)*sin(omegas)) # daily extraterrestrial radiation
		Rso <- Ra*(0.75+2e-5*elev)     #For a cloudless day, Rs is roughly 75% of extraterrestrial radiation (Ra)

		# radfraction is a measure of relative shortwave radiation, or of
		# possible radiation (cloudy vs. clear-sky), needs to be less than 1
		radfraction <- radiation/Rso
		radfraction[radfraction>1] <- 1
		
		# longwave radiation
		longw <- 4.903e-9*n_days*((tmax+273.15)^4+(tmin+273.15)^4)/2*(.34-.14*sqrt(ea))*(1.35*radfraction-.35)     
		
	    # net radiation
		netrad <- radiation*n_days*(1-albedo)-longw     
		
		# ET from long-form P-M eqn.
		et0 <- 0.408*((delta*(netrad-G))+(pa*cp*vpd/ra*3600*24*n_days))/(delta+gamma*(1+rs/ra))
		return(et0)
	} 


	aetmod <- function(et0,input,awc,soil_prev=NULL){
	  # This function computes AET given ET0, H2O input, soil water capacity, and beginning-of-month soil moisture
	  # Arguments:
	  # et0: vector of monthly reference evapotranspiration in mm
	  # input: vector of monthly water input to soil in mm
	  # awc: vector of soil water capacity in mm
	  # soil_prev: vector of soil water content for the previous month (mm).  If left NULL this is assigned to be zero.
	  #
	  # Value:
	  # returns a data frame with columns for AET, deficit, end-of-month soil moisture, and runoff.
	  
	  N <- length(et0)
	  runoff <- def <-  aet <- soil <- rep(NA,N) # 
	  if(is.null(soil_prev)) soil_prev <- rep(0,N)
	  
	  deltasoil <- input-et0 # positive=excess H2O, negative=H2O deficit
	  
	  # Case when there is a moisture surplus:
	  Case <- deltasoil>=0
	  if(sum(Case)>0){
		aet[Case] <- et0[Case]
		def[Case] <- 0
		soil[Case] <- pmin(soil_prev[Case]+deltasoil[Case],awc[Case])	# increment soil moisture, but not above water holding capacity
		runoff[Case] <- pmax(soil_prev[Case]+deltasoil[Case]-awc[Case],0) # when awc is exceeded, send the rest to runoff
	  }
	  
	  # Case where there is a moisture deficit:  soil moisture is reduced
	  Case <- deltasoil<0
	  if(sum(Case)>0){
		soildrawdown <- soil_prev[Case]*(1-exp(-(et0-input)[Case]/awc[Case]))	# this is the net change in soil moisture (neg)
		aet[Case] <- pmin(input[Case] + soildrawdown,et0[Case])
		def[Case] <- et0[Case] - aet[Case]
		soil[Case] <- soil_prev[Case]-soildrawdown
		runoff[Case] <- 0
	  }
	  
	  return(data.frame(aet=aet,def=def,soil=soil,runoff=runoff))
	  
	}

say('################################################################################################')
say('### calculate MONTHLY actual incoming solar radiation (see Appendix of Dobrowski et al 2013) ###')
say('################################################################################################')
	
	say('This script calculates actual incoming solar radiation from beam and diffuse radiation combined based on the scripts in the Appendix of Dobrowski et al. 2013. They show how to calculate coefficients of the fraction of beam and diffuse radiation striking the surface given semi-empirical estimates of cloudiness from the NDLAS forcing variables.  I have modified the procedure to use albedo estimated from average of snow cover across each month from the DayMet data set (Dobrowski calculate albedo *after* snow, whereas I introduce albedo in the calculation of incoming actual solar radiation, which is then used to calculate Dobrowski snow--yes, this is somewhat circular, but it seems better than assuming constant albedo). Almost all intermediate files are saved in GRASS format.', breaks=80)
	
	say('In GRASS create GRASS project "GRASS", location "nad83", and mapset "Adam" in "F:\\Climate\\PRISM\\GRASS" using PROJ4 from PRISM mask.', breaks=80)
	
	## set working region/projection/resolution
	maskPrism <- raster('F:/ecology/Climate/PRISM/PRISM_us_dem_800m.tif')
	maskPrism <- maskPrism * 0L + 1L
	maskPrismGrid <- as(maskPrism, 'SpatialGrid')

	# initialize GRASS
	rgrass7::use_sp()
	grass <- initGRASS(
		# gisBase='C:/Program Files/GRASS GIS 7.0.4',
		gisBase='C:/Program Files/GRASS GIS 7.8',
		gisDbase='F:/ecology/Climate/PRISM/GRASS',
		SG=maskPrismGrid,
		location='nad83',
		mapset='Adam',
		home='F:/ecology/Climate/PRISM/GRASS',
		remove_GISRC=TRUE,
		override=TRUE
	)

	# get elevation raster
	execGRASS(
		cmd='r.in.gdal',
		flags=c('o', 'overwrite'),
		parameters=list(
			input=paste0('F:/ecology/Climate/PRISM/PRISM_us_dem_800m.tif'),
			output='demPrism'
		)
	)

say('   ### create horizon rasters ###')

	# r.horizon -d --overwrite elevation=gmted2010_croppedToOchotonaPrincepsIucnRangePlus400km_prism@Adam step=5 maxdistance=1 output=horizonAngle_gmted2010_prism_radians_fromDegree
	execGRASS(
		cmd='r.horizon',
		flags=c('overwrite'),
		parameters=list(
			elevation='demPrism',
			output='horizonAngle_fromDegree',
			step=5,
			maxdistance=1
		)
	)

say('   ### calculate slope/aspect ###')

	# r.slope.aspect --overwrite elevation=gmted2010_croppedToOchotonaPrincepsIucnRangePlus400km_prism@Adam slope=gmted2010_slope_prism aspect=gmted2010_aspect_prism
	execGRASS(
		cmd='r.slope.aspect',
		flags=c('overwrite'),
		parameters=list(
			elevation='demPrism',
			slope='slope',
			aspect='aspect'
		)
	)

	# import mask raster to use as faux coefdh and coefbh values
	writeRAST(as(maskPrism, 'SpatialGridDataFrame'), 'coeff_bh_valEquals1', overwrite=TRUE)
	writeRAST(as(maskPrism, 'SpatialGridDataFrame'), 'coeff_dh_valEquals1', overwrite=TRUE)
	
say('   ### calculate potential incoming shortwave radiation ###')

	for (month in 1:12) {
	
		midDoyThisMonth <- doyNonLeap[15, paste0('month', month)]

		say('month ', month, ' doy ', midDoyThisMonth)

		# calculate potential incoming solar radiation
		execGRASS(
			cmd='r.sun',
			flags=c('overwrite', 'verbose'),
			parameters=list(
				elevation='demPrism',
				aspect='aspect',
				slope='slope',
				linke_value=1,
				albedo_value=0.2,
				coeff_bh='coeff_bh_valEquals1',
				coeff_dh='coeff_dh_valEquals1',
				horizon_basename='horizonAngle_fromDegree',
				horizon_step=5,
				glob_rad=paste0('globalClearSkyRadiation_month', prefix(month, 2), '_doy', prefix(midDoyThisMonth, 3)),
				day=midDoyThisMonth,
				step=0.05,
				distance_step=1
			)
		)
			
	} # next month
	
say('   ### calculate actual incoming solar radiation ###')

	dirCreate('F:/ecology/Climate/PRISM/Actual Incoming Solar Radiation - Dobrowski et al 2013')

	source('C:/ecology/Drive/R/Climate/Calculate Extrasolar Radiation Ra.r')
	
	# calculate long/lat rasters
	# source('C:/ecology/Drive/R/Geography/Rasters - Calculate Latitude and Longitude Rasters.r')
	lat <- longLatRasters(template=maskPrism, mask=maskPrism)[[2]]
	lat <- c(as.matrix(lat))

	for (year in 1980:2015) {
	
		for (month in 1:12) {
			
			# calculate radiation for 15th day of day each month
			doyThisMonth <- doyNonLeap[15, paste0('month', month)]
				
			say('   ', year, ' ', month, ' ', doyThisMonth, ' ', date(), ' ------------------------------------')
			
			say('   ### process Link turbidity coefficient raster ###')
			say('   ### projection/cropping to PRISM done outside this script--see "./Climate/Linke Turbidity Factors - Redmund et al 2003"')

				# # linke <- raster(paste0('C:/ecology/Climate/Linke Turbidity Factors - Redmund et al 2003/ORIGINAL (divide by 20 to get real values)/month', prefix(month, 2), '.tif'))
				
				# # projection(linke) <- crsWgs84
				# # extent(linke) <- extent(raster())
				# # linke <- linke / 20
				# # linke <- projectRaster(linke, maskPrism)
				
				# linke <- raster(paste0('C:/ecology/Climate/Linke Turbidity Factors - Redmund et al 2003/PRISM NAD83/linke_month', prefix(month, 2), '.tif'))
				# writeRAST(x=as(linke, 'SpatialGridDataFrame'), vname=paste0('linke_month', prefix(month, 2)), overwrite=TRUE)
				
				execGRASS(
					cmd='r.in.gdal',
					flags=c('o', 'overwrite'),
					parameters=list(
						input=paste0('F:/ecology/Climate/Linke Turbidity Factors - Redmund et al 2003/PRISM NAD83/linke_month', prefix(month, 2), '.tif'),
						output=paste0('linke_month', prefix(month, 2))
					)
				)

			say('   ### extrasolar radiation ###')
				
				SWtotal <- extraSol(J=doyThisMonth, phi=lat)
				SWtotal <- raster(matrix(SWtotal, ncol=ncol(maskPrism), byrow=F))
				projection(SWtotal) <- projection(maskPrism)
				extent(SWtotal) <- extent(maskPrism)
				names(SWtotal) <- paste0('solarRad_doy', prefix(doyThisMonth, 3))
				SWtotal <- 277.777777778 * SWtotal # convert MJ/m2/day to Wh/m2/day
		
			say('   ### import NLDAS downward incoming shortwave radiation ###')

				# get this year/month's NLDAS layers
				nldas <- readGDAL(paste0('C:/ecology/Integrated/NASA Land Data Assimilation Systems/NLDAS_FORA0125_M NLDAS Primary Forcing Data L4 Monthly 0.125 x 0.125 degree V002/NLDAS_FORA0125_M.A', year, prefix(month, 2), '.002.grb'), silent=TRUE) # read GRIB file

				nldas <- brick(nldas)
				
				# get incoming direct shortwave downward radiation layer
				SWsurface <- nldas[[11]]
				rm(nldas); gc()

				# convert NLDAS shortwave radiation from average W/m2 each hour to Wh/m2/day
				SWsurface <- SWsurface * 24
				
				# project/crop/resample to PRISM
				SWsurface <- projectRaster(SWsurface, maskPrism)
				gc()

				# writeRAST(x=as(SWsurface, 'SpatialGridDataFrame'), vname=paste0('SWsurface_year', year, '_month', prefix(month, 2)), overwrite=TRUE)
				# SWsurface <- readRAST(paste0('SWsurface_year', year, '_month', prefix(month, 2)))
				# SWsurface <- raster(SWsurface)
				
			say('   ### import clear-sky downward shortwave radiation ###')
			
				# port clear-sky radiation from GRASS to R
				execGRASS(
					cmd='r.out.gdal',
					flags=c('overwrite', 'quiet'),
					parameters=list(
						input=paste0('globalClearSkyRadiation_month', prefix(month, 2), '_doy', prefix(doyThisMonth, 3)),
						output=paste0('C:/ecology/!Scratch/globalClearSkyRadiation_month', prefix(month, 2), '_doy', prefix(doyThisMonth, 3), '.tif'),
						format='GTiff',
						type='Float32'
					)
				)
				
				SWclear <- raster(paste0('C:/ecology/!Scratch/globalClearSkyRadiation_month', prefix(month, 2), '_doy', prefix(doyThisMonth, 3), '.tif'))
				SWclear <- projectRaster(SWclear, maskPrism)
				gc()
				SWclear <- SWclear * maskPrism

			say('   ### calculate direct/diffuse coefficients ###')

			# calculate coefficients
				coefbh <- min(maskPrism, (SWsurface / SWclear) * (1 - (-1.43 * (SWsurface / SWtotal) + 1.16)))
				coefdh <- min(maskPrism, (SWsurface / SWclear) * (-1.43 * (SWsurface / SWtotal) + 1.16))
			
				# port coefficient rasters to GRASS
				writeRAST(x=as(coefbh, 'SpatialGridDataFrame'), vname=paste0('coeff_bh_year', year, '_month', prefix(month, 2)), overwrite=TRUE)
				writeRAST(x=as(coefdh, 'SpatialGridDataFrame'), , vname=paste0('coeff_dh_year', year, '_month', prefix(month, 2)), overwrite=TRUE)
			
				rm(coefbh, coefdh); gc()

			say('   ### get albedo raster ###')

				albedo <- raster(paste0('C:/ecology/Climate/DayMet/V3/Albedo/meanAlbedo_year', year, '_month', prefix(month, 2), '.tif'))
				albedo <- projectRaster(albedo, maskPrism)
				albedo <- albedo * maskPrism
			
				writeRAST(x=as(albedo, 'SpatialGridDataFrame'), vname=paste0('meanAlbedo_year', year, '_month', prefix(month, 2)), overwrite=TRUE)
			
			say('   ### calculate actual incoming solar radiation ###')
			
				# calculate ACTUAL incoming solar radiation
				execGRASS(
					cmd='r.sun',
					flags=c('overwrite', 'quiet'),
					parameters=list(
						elevation='demPrism',
						aspect='aspect',
						slope='slope',
						linke=paste0('linke_month', prefix(month, 2)),
						albedo=paste0('meanAlbedo_year', year, '_month', prefix(month, 2)),
						coeff_bh=paste0('coeff_bh_year', year, '_month', prefix(month, 2)),
						coeff_dh=paste0('coeff_dh_year', year, '_month', prefix(month, 2)),
						horizon_basename='horizonAngle_fromDegree',
						horizon_step=5,
						glob_rad=paste0('globalActualRadiation_year', year, '_month', prefix(month, 2), '_doy', prefix(doyThisMonth, 3)),
						day=doyThisMonth,
						step=0.05,
						distance_step=1
					)
				)
				
				actualRadiation <- readRAST(paste0('globalActualRadiation_year', year, '_month', prefix(month, 2), '_doy', prefix(doyThisMonth, 3)))
				actualRadiation <- raster(actualRadiation)
				actualRadiation <- actualRadiation / 277.777777778 # convert to MJ/m2/day
				actualRadiation <- projectRaster(actualRadiation, maskPrism)
				gc()
				
				writeRaster(actualRadiation, paste0('F:/ecology/Climate/PRISM/Actual Incoming Solar Radiation - Dobrowski et al 2013/globalActualDownwardShortwaveRadiation_year', year, '_month', prefix(month, 2), '_doy', prefix(doyThisMonth, 3), '_MJperM2perDay'), format='GTiff', overwrite=TRUE)
				
		} # next month
		
	} # next year

# say('###########################################################################################')
# say('### calculate MONTHLY snow snow water equivalent (see Appendix of Dobrowski et al 2013) ###')
# say('###########################################################################################')

	# say('For January 1980 assuming previous month snowpack was 0.')

	# dirCreate('F:/ecology/Climate/PRISM/Snow Water Equivalent - Monthly - Dobrowski et al 2013')
	
	# maskPrism <- raster('F:/ecology/Climate/PRISM/PRISM_us_dem_800m.tif')
	# maskPrism <- maskPrism * 0
	
	# for (year in 1980:2015) {
	
		# for (month in 1:12) {
		
			# say('Monthly SWE for ', year, ' ', month, ' ', date())
	
			# # mid-day of this month
			# doyThisMonth <- doyNonLeap[15, paste0('month', month)]
	
			# # mean precip
			# tmin <- raster(paste0('F:/ecology/Climate/PRISM/AN81m 1981-2015/tmin/', year, '/prism_tmin_us_30s_', year, prefix(month, 2), '.tif'))
			# tmax <- raster(paste0('F:/ecology/Climate/PRISM/AN81m 1981-2015/tmax/', year, '/prism_tmax_us_30s_', year, prefix(month, 2), '.tif'))
			# tmean <- (tmin + tmax) / 2
			# tmean <- tmean / 100000
			# names(tmean) <- 'tmean'
			# rm(tmin, tmax); gc()
	
			# # previous month's snowpack
			# snowpack_prev <- if (year == 1980 & month == 1) {
				# maskPrism
			# } else if (month == 1 ) {
				# raster(paste0('F:/ecology/Climate/PRISM/Snow Water Equivalent - Monthly - Dobrowski et al 2013/snowWaterEquivalent_year', year - 1, '_month', prefix(12, 2), '_mm.tif'))
			# } else {
				# raster(paste0('F:/ecology/Climate/PRISM/Snow Water Equivalent - Monthly - Dobrowski et al 2013/snowWaterEquivalent_year', year, '_month', prefix(month - 1, 2), '_mm.tif'))
			# }
			# names(snowpack_prev) <- 'snowpack_prev'
	
			# # precip
			# ppt <- raster(paste0('F:/ecology/Climate/PRISM/AN81m 1981-2015/ppt/', year, '/prism_ppt_us_30s_', year, prefix(month, 2), '.tif'))
			# ppt <- ppt / 100000
			# names(ppt) <- 'ppt'
			
			# # actual incoming solar radiation
			# radiation <- raster(paste0('F:/ecology/Climate/PRISM/Actual Incoming Solar Radiation - Dobrowski et al 2013/globalActualDownwardShortwaveRadiation_year', year, '_month', prefix(month, 2), '_doy', prefix(doyThisMonth, 3), '_MJperM2perDay.tif'))
			# names(radiation) <- 'radiation'
	
			# sweVec <- snowmodModified(
				# tmean=c(as.matrix(tmean)),
				# ppt=c(as.matrix(ppt)),
				# radiation=c(as.matrix(radiation)),
				# snowpack_prev=c(as.matrix(snowpack_prev)),
				# albedo=0.23,
				# albedo_snow=0.8
			# )
		
			# rm(ppt, snowpack_prev); gc()
			
			# swe <- matrix(sweVec, nrow=nrow(tmean), ncol=ncol(tmean), byrow=FALSE)
			# rm(sweVec); gc()
			# swe <- raster(swe)
			# extent(swe) <- extent(tmean)
			# projection(swe) <- projection(tmean)
			# names(swe) <- paste0('snowWaterEquivalent_year', year, '_month', prefix(month, 2), '_mm')

			# writeRaster(swe, paste0('F:/ecology/Climate/PRISM/Snow Water Equivalent - Monthly - Dobrowski et al 2013/snowWaterEquivalent_year', year, '_month', prefix(month, 2), '_mm'), format='GTiff', overwrite=TRUE)
			
			# rm(swe); gc()
			
		# } # next month
		
	# } # next year

# say('#########################################################################')
# say('### calculate mean annual SWE from monthly SWE (Dobrowski et al 2013) ###')
# say('#########################################################################')

# dirCreate('F:/ecology/Climate/PRISM/Snow Water Equivalent - Annual - Mean - Dobrowski et al 2013')

# say('Skipping 1980 because no available snow pack from 1979!')

# for (year in 1981:2015) {

	# say(year)

	# monthlys <- stack(listFiles('F:/ecology/Climate/PRISM/Snow Water Equivalent - Monthly - Dobrowski et al 2013', pattern=paste0('snowWaterEquivalent_year', year)))
	# meanAnnual <- mean(monthlys)
	# names(meanAnnual) <- paste0('meanAnnualSwe_year', year)
	# writeRaster(meanAnnual, paste0('F:/ecology/Climate/PRISM/Snow Water Equivalent - Annual - Mean - Dobrowski et al 2013/snowWaterEquivalent_year', year, '_mm'), format='GTiff', overwrite=TRUE)
	# rm(monthlys, meanAnnual); gc()

# }
	
# say('##################################################################')
# say('### calculate DAILY ET0 (see Appendix of Dobrowski et al 2013) ###')
# say('##################################################################')

	# dirCreate('F:/ecology/Climate/PRISM/ET0 - Daily - Dobrowski et al 2013')
	# sink('F:/ecology/Climate/PRISM/ET0 - Daily - Dobrowski et al 2013/!Divide all values by 100,000 to get mm per day.txt')
	# sink()
	
	# maskPrism <- raster('F:/ecology/Climate/PRISM/PRISM_us_dem_800m.tif')
	# maskPrism <- maskPrism * 0 + 1

	# # latitude
	# lat <- longLatRasters(template=maskPrism)[[2]]
	# lat <- c(as.matrix(lat))
	
	# # elevation
	# elev <- raster('F:/ecology/Climate/PRISM/PRISM_us_dem_800m.tif')
	# elev <- c(as.matrix(elev))
	
	# for (year in 1981:2015) {

		# for (month in 1:12) {
		
			# doyThisYear <- if (year %% 4 == 0) {
				# doyLeapYear
			# } else {
				# doyNonLeap
			# }
		
			# # mid-day of this month
			# midMonthDoy <- doyNonLeap[15, paste0('month', month)]

			# # get monthly input rasters
			# radiation <- raster(paste0('F:/ecology/Climate/PRISM/Actual Incoming Solar Radiation - Dobrowski et al 2013/globalActualDownwardShortwaveRadiation_year', year, '_month', prefix(month, 2), '_doy', prefix(midMonthDoy, 3), '_MJperM2perDay.tif'))
			# radiation <- c(as.matrix(radiation))
			# gc()
			
			# wind <- raster(paste0('F:/ecology/Integrated/LDAS/NASA LDAS/NLDAS_FORA0125_M NLDAS Primary Forcing Data L4 Monthly 0.125 x 0.125 degree V002/Mean Wind Speed/meanWindSpeed_', year, '_', prefix(month, 2), '.tif'))
			# wind <- projectRaster(wind, maskPrism)
			# wind <- c(as.matrix(wind))
			# gc()
		
			# dpt <- raster(paste0('F:/ecology/Climate/PRISM/AN81m 1981-2015/tdmean/', year, '/prism_tdmean_us_30s_', year, prefix(month, 2), '.tif'))
			# dpt <- c(as.matrix(dpt)) / 100000
			# gc()
			
			# for (doy in na.omit(doyThisYear[ , paste0('month', month)])) {
			
				# say('Daily ET0 for ', year, ' ', month, ' ', doy, ' ', date())
			
				# snow <- raster(paste0('F:/ecology/Climate/DayMet/V3/GeoTIFF/swe_year', year, '_doy', prefix(doy, 3), '.tif'))
				# snow <- projectRaster(snow, maskPrism)
				# albedo <- calc(snow, fun=function(x) ifelse(x > 0, 0.8, 0.23))
				# albedo <- c(as.matrix(albedo))
				# gc()

				# tmin <- raster(paste0('F:/ecology/Climate/PRISM/AN81d 1981-2015/tmin/', year, '/prism_tmin_us_30s_', year, prefix(month, 2), prefix(which(doyThisYear[ , paste0('month', month)] == doy), 2), '.tif'))
				# tmin <- c(as.matrix(tmin)) / 100000
				# gc()
				
				# tmax <- raster(paste0('F:/ecology/Climate/PRISM/AN81d 1981-2015/tmax/', year, '/prism_tmax_us_30s_', year, prefix(month, 2), prefix(which(doyThisYear[ , paste0('month', month)] == doy), 2), '.tif'))
				# tmax <- c(as.matrix(tmax)) / 100000
				# gc()
				
				# # get daily input rasters
				# et0 <- dailyET0(
					# radiation=radiation,
					# tmax=tmax,
					# tmin=tmin,
					# wind=wind,
					# lat=lat,
					# elev=elev,
					# albedo=albedo,
					# dpt=dpt,
					# doy=doy
				# )

				# # handle very large/small values by forcing to NA then (later) through interpolation
				# et0 <- ifelse(et0 < 0, 0, et0)
				# et0 <- ifelse(et0 > 1000, NA, et0)
				
				# et0 <- matrix(et0, nrow=nrow(maskPrism), ncol=ncol(maskPrism), byrow=FALSE)
				# et0 <- raster(et0)
				# extent(et0) <- extent(maskPrism)
				# projection(et0) <- projection(maskPrism)

				# # interpolate
				# et0 <- focal(et0, w=matrix(rep(1/9, 9), nrow=3, ncol=3), na.rm=TRUE, NAonly=TRUE)
				# et0 <- et0 * maskPrism

				# et0 <- et0 * 100000
				# et0 <- round(et0)
				
				# names(et0) <- paste0('et0_year', year, '_month', prefix(month, 2), '_mm')

				# writeRaster(et0, paste0('F:/ecology/Climate/PRISM/ET0 - Daily - Dobrowski et al 2013/et0_year', year, '_doy', prefix(doy, 3)), format='GTiff', datatype='INT4S', overwrite=TRUE)
				
				# rm(et0); gc()
				
			# }

		# } # next month

		# rm(radiation, wind, albedo, tmin, tmax, dpt); gc()
		
	# } # next year

# say('####################################################################')
# say('### calculate MONTHLY ET0 (see Appendix of Dobrowski et al 2013) ###')
# say('####################################################################')

	# dirCreate('F:/ecology/Climate/PRISM/ET0 - Monthly - Dobrowski et al 2013')
	# sink('F:/ecology/Climate/PRISM/ET0 - Monthly - Dobrowski et al 2013/!Divide all values by 100,000 to get mm per month.txt')
	# sink()
	
	# maskPrism <- raster('F:/ecology/Climate/PRISM/PRISM_us_dem_800m.tif')
	# maskPrism0 <- maskPrism * 0
	# maskPrism1 <- maskPrism + 1

	# # latitude
	# lat <- longLatRasters(template=maskPrism)[[2]]
	# lat <- c(as.matrix(lat))
	
	# # elevation
	# elev <- raster('F:/ecology/Climate/PRISM/PRISM_us_dem_800m.tif')
	# elev <- c(as.matrix(elev))
	
	# for (year in 1981:2015) {
	# # for (year in 2015) {
	
		# for (month in 1:12) {
		# # for (month in 1:6) {
		# # for (month in 7:12) {
		# # for (month in 12) {
		
			# say('Monthly ET0 for ', year, ' ', month, ' ', date())
			
			# # mid-day of this month
			# midMonthDoy <- doyNonLeap[15, paste0('month', month)]

			# # get monthly input rasters
			# radiation <- raster(paste0('F:/ecology/Climate/PRISM/Actual Incoming Solar Radiation - Dobrowski et al 2013/globalActualDownwardShortwaveRadiation_year', year, '_month', prefix(month, 2), '_doy', prefix(midMonthDoy, 3), '_MJperM2perDay.tif'))
			# radiation <- c(as.matrix(radiation))
			# gc()
			
			# wind <- raster(paste0('F:/ecology/Integrated/LDAS/NASA LDAS/NLDAS_FORA0125_M NLDAS Primary Forcing Data L4 Monthly 0.125 x 0.125 degree V002/Mean Wind Speed/meanWindSpeed_', year, '_', prefix(month, 2), '.tif'))
			# wind <- projectRaster(wind, maskPrism)
			# wind <- c(as.matrix(wind))
			# gc()
		
			# # albedo <- raster(paste0('F:/ecology/Climate/DayMet/V3/Albedo/meanAlbedo_year', year, '_month', prefix(month, 2), '.tif'))
			# # albedo <- projectRaster(albedo, maskPrism)
			# # albedo <- c(as.matrix(albedo))
			# albedo <- 0.23
			# gc()

			# tmin <- raster(paste0('F:/ecology/Climate/PRISM/AN81m 1981-2015/tmin/', year, '/prism_tmin_us_30s_', year, prefix(month, 2), '.tif'))
			# tmin <- c(as.matrix(tmin)) / 100000
			# gc()
			
			# tmax <- raster(paste0('F:/ecology/Climate/PRISM/AN81m 1981-2015/tmax/', year, '/prism_tmax_us_30s_', year, prefix(month, 2), '.tif'))
			# tmax <- c(as.matrix(tmax)) / 100000
			# gc()
			
			# tmin_prev <- if (year == 1981 & month == 1) {
				# raster('F:/ecology/Climate/PRISM/AN81m 1895-2015/tmin/1981/prism_tmin_us_30s_198112.tif')
			# } else if (month == 1) {
				# raster(paste0('F:/ecology/Climate/PRISM/AN81m 1981-2015/tmin/', year - 1, '/prism_tmin_us_30s_', year - 1, '12.tif'))
			# } else {
				# raster(paste0('F:/ecology/Climate/PRISM/AN81m 1981-2015/tmin/', year, '/prism_tmin_us_30s_', year, prefix(month, 2), '.tif'))
			# }
			
			# tmax_prev <- if (year == 1981 & month == 1) {
				# raster('F:/ecology/Climate/PRISM/AN81m 1895-2015/tmax/1981/prism_tmax_us_30s_198112.tif')
			# } else if (month == 1) {
				# raster(paste0('F:/ecology/Climate/PRISM/AN81m 1981-2015/tmax/', year - 1, '/prism_tmax_us_30s_', year - 1, '12.tif'))
			# } else {
				# raster(paste0('F:/ecology/Climate/PRISM/AN81m 1981-2015/tmax/', year, '/prism_tmax_us_30s_', year, prefix(month, 2), '.tif'))
			# }
			# tmean_prev <- tmin_prev + tmax_prev / 2
			# rm(tmin_prev, tmax_prev)
			# tmean_prev <- c(as.matrix(tmean_prev)) / 100000
			# gc()
			
			# dpt <- raster(paste0('F:/ecology/Climate/PRISM/AN81m 1981-2015/tdmean/', year, '/prism_tdmean_us_30s_', year, prefix(month, 2), '.tif'))
			# dpt <- c(as.matrix(dpt)) / 100000
			# gc()

			# # PET
			# vals <- monthlyET0(
				# radiation=radiation,
				# tmax=tmax,
				# tmin=tmin,
				# wind=wind,
				# lat=lat,
				# elev=elev,
				# dpt=dpt,
				# tmean_prev=tmean_prev,
				# albedo=albedo,
				# month=month
			# )

			# rm(radiation, wind, albedo, tmean_prev, dpt); gc()

			# # convert to raster
			# vals[vals < 0] <- 0
			
			# et0 <- matrix(vals, nrow=nrow(maskPrism), ncol=ncol(maskPrism), byrow=FALSE)
			# et0 <- raster(et0)
			# extent(et0) <- extent(maskPrism)
			# projection(et0) <- projection(maskPrism)

			# # interpolate and force NAs to 0
			# # et0 <- focal(et0, w=matrix(rep(1/9, 9), nrow=3, ncol=3), na.rm=TRUE, NAonly=TRUE)
			# # et0 <- calc(et0, fun=function(x) ifelse(is.na(x), 0, x))
			# # et0 <- et0 * maskPrism1

			# et0 <- et0 * 100000
			# et0 <- round(et0)
			# et0 <- setMinMax(et0)
			
			# names(et0) <- paste0('et0_year', year, '_month', prefix(month, 2), '_mm')

			# writeRaster(et0, paste0('F:/ecology/Climate/PRISM/ET0 - Monthly - Dobrowski et al 2013/et0_year', year, '_month', prefix(month, 2)), format='GTiff', datatype='INT4S', overwrite=TRUE)
			
			# rm(et0); gc()

		# } # next month
		
	# } # next year

# say('####################################################################################')
# say('### calculate ANNUAL ET0 from monthly ET0 (see Appendix of Dobrowski et al 2013) ###')
# say('####################################################################################')

	# dirCreate('F:/ecology/Climate/PRISM/ET0 - Annual - Sum - Dobrowski et al 2013')
	# sink('F:/ecology/Climate/PRISM/ET0 - Annual - Sum - Dobrowski et al 2013/!Divide all values by 100,000 to get mm per month.txt')
	# sink()
	
	# # for (year in 1981:2015) {
	
		# say('Annual ET0 for ', year, ' ', date())

		# monthly <- stack(list.files('F:/ecology/Climate/PRISM/ET0 - Monthly - Dobrowski et al 2013', pattern=paste0('et0_year', year), full.names=TRUE))
		# annual <- sum(monthly)
		
		# names(annual) <- paste0('et0_year', year, '_mm')

		# writeRaster(annual, paste0('F:/ecology/Climate/PRISM/ET0 - Annual - Sum - Dobrowski et al 2013/et0_year', year), format='GTiff', datatype='INT4S', overwrite=TRUE)
		
		# rm(annual, monthly); gc()

	# } # next year

	
say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
say('DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')