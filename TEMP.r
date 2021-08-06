# source('E:/Ecology/Drive/Research/Pikas - California Gap (Erik Beever et al)/Code/TEMP.r')


	### colors for points
	trainPresCol <- 'black'
	testPresCol <- 'darkgreen'
	testPresFill <- 'chartreuse3'
	testShortTermAbsCol <- 'darkorange4'
	testShortTermAbsFill <- 'darkgoldenrod3'
	testLongTermAbsCol <- 'darkred'
	incorrectPredCol <- 'cyan'
		
	### colors for rasters
	predPresCol <- alpha('forestgreen', 0.6)
	predShortTermAbsCol <- alpha('darkgoldenrod3', 0.6)
	predLongTermAbsCol <- alpha('red', 0.4)
	
	trainPresPch <- 3

	testLongTermPch <- 25
	testShortTermPch <- 23
	testPresPch <- 21

	predLongTerm <- predVals[testSurveys$status == '0 long absence'] / 1000
	predShortTerm <- predVals[testSurveys$status == '1 recent absence'] / 1000
	predPres <- predVals[testSurveys$status == '2 detected'] / 1000


	tholdShortVsLong <- thresholdWeighted(predShortTerm, predLongTerm)[['mdss']]
	tholsShortVsPres <- thresholdWeighted(predPres, predShortTerm)[['mdss']]

	### plot
	png('./Figures & Tables/ENM Prediction.png', width=2 * 1200, height=2 * 1200, res=300)
		
		par(mfrow=c(2, 2), mar=c(0.1, 0.1, 0.1, 0.1), oma=c(3, 3, 0.1, 0.1), mgp=c(3, 0.4, 0))

		### each class by itself
		
		for (class in c(0:2)) {
		
			plot(focus, border='white')
			usr <- par('usr')
			
			if (class == 0) {

				predClass <- (pred / 1000) < tholdShortVsLong
				longTermRast <- predClass
				pointPreds <- sort(predLongTerm)

				col <- predLongTermAbsCol
				pointOutlineCol <- testLongTermAbsCol
				pointFillCol <- rep(NA, length(pointPreds))
				pointFillCol[pointPreds >= tholdShortVsLong & pointPreds < tholsShortVsPres] <- testShortTermAbsFill
				pointFillCol[pointPreds >= tholsShortVsPres] <- testPresFill
				rastLab <- 'No evidence predicted'
				pointLab <- c('Correctly predicted', 'Old evidence predicted', 'Presence predicted')
				pch <- testLongTermPch
				axis(2, at=pretty(c(usr[3:4]), 3), tck=-0.01)
				
				legColBg <- c(NA, NA, testShortTermAbsFill, testPresFill)
				
			} else if (class == 1) {

				predClass <- (pred / 1000) >= tholdShortVsLong & (pred / 1000) < tholsShortVsPres
				shortTermRast <- predClass
				pointPreds <- sort(predShortTerm)

				col <- predShortTermAbsCol
				pointOutlineCol <- testShortTermAbsCol
				pointFillCol <- rep(NA, length(pointPreds))
				pointFillCol[pointPreds < tholdShortVsLong] <- testLongTermAbsCol
				pointFillCol[pointPreds < tholdShortVsLong & pointPreds < tholsShortVsPres] <- NA
				pointFillCol[pointPreds >= tholsShortVsPres] <- testPresFill
				rastLab <- 'Old evidence predicted'
				pointLab <- c('No evidence predicted', 'Correctly predicted', 'Presence predicted')
				pch <- testShortTermPch
				
				legColBg <- c(NA, testLongTermAbsCol, NA, testPresFill)
				
			} else if (class == 2) {
			
				predClass <- (pred / 1000) >= tholsShortVsPres
				presRast <- predClass
				pointPreds <- sort(predPres, FALSE)

				col <- predPresCol
				pointOutlineCol <- testPresCol
				pointFillCol <- rep(NA, length(pointPreds))
				pointFillCol[pointPreds < tholdShortVsLong] <- testLongTermAbsCol
				pointFillCol[pointPreds >= tholdShortVsLong & pointPreds < tholsShortVsPres] <- testShortTermAbsFill
				pointFillCol[pointPreds >= tholsShortVsPres] <- NA
				rastLab <- 'Presence predicted'
				pointLab <- c('No evidence predicted', 'Old evidence predicted', 'Correctly predicted')
				pch <- testPresPch
				axis(1, at=pretty(c(usr[1:2]), 3)[1:3], tck=-0.01)
				axis(2, at=pretty(c(usr[3:4]), 3), tck=-0.01)
				
				legColBg <- c(NA, testLongTermAbsCol, testShortTermAbsFill, NA)
				
			}
			
			plot(hs, col=alpha(grays, 0.8), legend=FALSE, add=TRUE)
			plot(predClass, col=c(NA, col), legend=FALSE, add=TRUE)
			plot(west2, border='gray40', add=TRUE)
			plot(west1, lwd=2, add=TRUE)
			points(trainPresSpEa, pch=3, cex=0.8, col=trainPresCol)

			pchs <- rep(NA, nrow(testSurveys))
			
			pchs[testSurveysSpEa$status == '0 long absence'] <- testLongTermPch
			pchs[testSurveysSpEa$status == '1 recent absence'] <- testShortTermPch
			pchs[testSurveysSpEa$status == '2 detected'] <- testPresPch



			if (class == 0) {
				points(testSurveysSpEa[testSurveysSpEa$status == '0 long absence', ], bg=pointFillCol, col=pointOutlineCol, pch=pch, cex=1.2)
			} else if (class == 1) {
				points(testSurveysSpEa[testSurveysSpEa$status == '1 recent absence', ], bg=pointFillCol, col=pointOutlineCol, pch=pch, cex=1.2)
			} else if (class == 2) {
				points(testSurveysSpEa[testSurveysSpEa$status == '2 detected', ], bg=pointFillCol, col=pointOutlineCol, pch=pch, cex=1.5)
			}
				
			box()
			
			legend(
				'bottomleft',
				inset=0.01,
				legend=c(rastLab, pointLab),
				pch=c(NA, pch, pch, pch),
				fill=c(col, NA, NA, NA),
				border=c('black', NA, NA, NA),
				col=c(NA, pointOutlineCol, pointOutlineCol, pointOutlineCol),
				pt.bg=legColBg,
				bg=alpha('white', 0.9),
				cex=0.8,
				pt.cex=1.2
			)
			
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

			
		}

		### all together

		plot(focus, border='white')
		axis(1, at=pretty(c(usr[1:2]), 3)[1:3], tck=-0.01)
		plot(hs, col=alpha(grays, 0.8), legend=FALSE, add=TRUE)
		plot(longTermRast, col=c(NA, predLongTermAbsCol), legend=FALSE, add=TRUE)
		plot(shortTermRast, col=c(NA, predShortTermAbsCol), legend=FALSE, add=TRUE)
		plot(presRast, col=c(NA, predPresCol), legend=FALSE, add=TRUE)
		plot(west2, border='gray40', add=TRUE)
		plot(west1, lwd=2, add=TRUE)
		points(trainPresSpEa, pch=trainPresPch, cex=0.8)

		testSurveysSpEa <- testSurveysSpEa[order(testSurveysSpEa$status), ]
		
		cols <- rep(NA, nrow(testSurveys))
		cols[testSurveysSpEa$status == '0 long absence'] <- testLongTermAbsCol
		cols[testSurveysSpEa$status == '1 recent absence'] <- testShortTermAbsCol
		cols[testSurveysSpEa$status == '2 detected'] <- testPresCol
		
		cexs <- pchs <- rep(NA, nrow(testSurveys))
		
		pchs[testSurveysSpEa$status == '0 long absence'] <- testLongTermPch
		pchs[testSurveysSpEa$status == '1 recent absence'] <- testShortTermPch
		pchs[testSurveysSpEa$status == '2 detected'] <- testPresPch
		
		cexs[testSurveysSpEa$status == '0 long absence'] <- 1.2
		cexs[testSurveysSpEa$status == '1 recent absence'] <- 1.2
		cexs[testSurveysSpEa$status == '2 detected'] <- 1.5

		points(testSurveysSpEa, col=cols, pch=pchs, cex=cexs)
			
		box()
		
		legend(
			'bottomleft',
			inset=0.01,
			legend=c('No evidence predicted', 'Old evidence predicted', 'Occurrence predicted', 'No evidence observed', 'Old evidence observed', 'Occurrence observed', 'Training occurrence'),
			pch=c(NA, NA, NA, testLongTermPch, testShortTermPch, testPresPch, trainPresPch),
			fill=c(predLongTermAbsCol, predShortTermAbsCol, predPresCol, NA, NA, NA, NA),
			border=c('black', 'black', 'black', NA, NA, NA, NA),
			col=c(NA, NA, NA, testLongTermAbsCol, testShortTermAbsCol, testPresCol, trainPresCol),
			bg=alpha('white', 0.4),
			cex=0.8,
			pt.cex=1.2
		)
		
			# scale bar
			size <- 50000 # length of scale bar in meters
			x <- usr[2] - size - 0.02 * size
			x <- c(x, x + size)
			y <- usr[3] + rep(0.02 * height, 2)
			
			lines(x, y, lwd=4, xpd=NA, col='black', lend=1)
			
			x <- mean(x)
			y <- usr[3] + 0.05 * height
			text(x, y[1], labels=paste(size / 1000, 'km'), cex=1)

		title(sub=date(), cex.sub=0.3, outer=TRUE, line=2)
			
		### inset map
		
			usr <- par('usr')
			par(fig=c(0.22, 0.64, 0.18, 0.54), bg='white', new=TRUE)
			
			# get extent of inset
			largeFocus <- nam[nam@data$NAME_1 %in% c('California', 'Nevada', 'Oregon', 'Baja California'), ]
			largeFocus <- gBuffer(largeFocus, width=50000)
			insetExt <- extent(largeFocus)
			insetExt <- as(insetExt, 'SpatialPolygons')
			projection(insetExt) <- getCRS('climateNA')
			insetExt <- vect(insetExt)
			namVect <- vect(nam)
			namCrop <- terra::crop(namVect, insetExt)
			
			# get polygon to highlight focal area
			focusExt <- extent(usr)
			focusExt <- as(focusExt, 'SpatialPolygons')
			projection(focusExt) <- getCRS('climateNA')

			plot(insetExt, border='gray', col='white', axes=FALSE)
			plot(namCrop, col='gray', border='gray40', ann=FALSE, add=TRUE)
			plot(focusExt, lwd=2.2, border='black', add=TRUE)
			plot(insetExt, lwd=1.2, axes=FALSE, add=TRUE)
				
	dev.off()

