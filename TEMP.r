# source('C:/Ecology/Drive/Research/Pikas - California Gap (Erik Beever et al)/pika_californiaGap/TEMP.r')

	# plot
	png('./Figures & Tables/Gap Sampling.png', width=1200, height=1200, res=300)
		
		par(mar=c(2, 1, 1, 1), cex.axis=0.4, mgp=c(3, 0, 0.2))
		
		plot(focus, border='white')
		
		usr <- par('usr')
		xs <- pretty(c(usr[1], usr[2]))
		ys <- pretty(c(usr[3], usr[4]))
		axis(1, at=xs, tck=0.01, labels=xs, col=NA, col.ticks='black')
		axis(2, at=ys, tck=0.01, labels=ys, col=NA, col.ticks='black')
		
		plot(hs, col=grays, legend=FALSE, add=TRUE)
		plot(lakes, col = 'gray60', border = NA, add = TRUE)
		plot(kdeClass, col=kdeCols, legend=FALSE, add=TRUE)
		plot(west2, border='gray40', add=TRUE)
		plot(west1, border='gray40', lwd=2, add=TRUE)

		testSurveysSpEa <- testSurveysSpEa[order(testSurveysSpEa$status), ]
		
		cols <- rep(NA, nrow(testSurveys))
		cols[testSurveysSpEa$status == '0 long absence'] <- testLongTermAbsCol
		cols[testSurveysSpEa$status == '1 recent absence'] <- testShortTermAbsCol
		cols[testSurveysSpEa$status == '2 detected'] <- testPresCol
		
		pchs <- rep(NA, nrow(testSurveys))
		pchs[testSurveysSpEa$status == '0 long absence'] <- testLongTermPch
		pchs[testSurveysSpEa$status == '1 recent absence'] <- testShortTermPch
		pchs[testSurveysSpEa$status == '2 detected'] <- testPresPch
		
		points(testSurveysSpEa, pch=pchs, col=cols, cex=0.6)
		points(trainPresSpEa, pch=trainPresPch, cex=0.5, col=trainPresCol)
		box()
		
		# legend
		legendBreaks(
			'bottomleft',
			inset=0.01,
			height=0.44,
			width=0.27,
			title='Previously-known\noccurrences',
			titleAdj=c(0.5, 0.92),
			col=kdeCols,
			adjX=c(0.05, 0.225),
			adjY=c(0.38, 0.85),
			labels=c('\U2265min presence', paste0('\U2265', '5th percentile'), paste0('  \U2265', '10th percentile')),
			labAdjX=0.58,
			cex=0.54,
			boxBg=alpha('white', 0.8)
		)
		
		usr <- par('usr')
		width <- usr[2] - usr[1]
		height <- usr[4] - usr[3]
		x <- usr[1] + 0.01 * width
		y <- usr[3] + 0.17 * height
		
		legend(
			x,
			y,
			legend=c('Previously-known occ.', 'Long-term test absence', 'Recent test absence', 'Test presence'),
			pch=c(trainPresPch, testLongTermPch, testShortTermPch, testPresPch),
			col=c(trainPresCol, testLongTermAbsCol, testShortTermAbsCol, testPresCol),
			bty='n',
			title='Surveys',
			cex=0.45,
			pt.cex=0.7
		)
		
		# scale bar
		size <- 50000 # length of scale bar in meters
		x <- usr[2] - size - 0.02 * width
		x <- c(x, x + size)
		y <- usr[3] + rep(0.02 * height, 2)
		
		lines(x, y, lwd=4, xpd=NA, col='black', lend=1)
		
		x <- mean(x)
		y <- usr[3] + 0.05 * height
		text(x, y[1], labels=paste(size / 1000, 'km'), cex=0.5)

		title(sub=date(), cex.sub=0.3, outer=TRUE, line=-1)
		
	dev.off()
