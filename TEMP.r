# source('C:/Ecology/Drive/Research/Iconic Species/Analysis - California Gap/Code/TEMP.r')

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
