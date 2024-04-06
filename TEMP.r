# source('C:/Ecology/Research/Pikas - California Gap (Erik Beever et al)/pika_californiaGap/TEMP.r')

	# plot
	png('./Figures & Tables/Gap Sampling Plumas Emphasized with Range Maps & No Survey Sites.png', width=1200, height=1200, res=300)
		
		par(mar=c(2, 1, 1, 1), cex.axis=0.4, mgp=c(3, 0, 0.2), fig=c(0, 1, 0, 1))
		
		plot(focus, border='white')
		
		usr <- par('usr')
		xs <- pretty(c(usr[1], usr[2]))
		ys <- pretty(c(usr[3], usr[4]))
		axis(1, at=xs, tck=0.01, labels=xs, col=NA, col.ticks='black')
		axis(2, at=ys, tck=0.01, labels=ys, col=NA, col.ticks='black')
		
		plot(hs, col=grays, legend=FALSE, add=TRUE)
		plot(lakes, col = 'gray60', border = NA, add = TRUE)

		iucnSp <- as(iucn, 'Spatial')
		iucnSp <- sp::spTransform(iucnSp, CRS(proj4string(west1)))
		iucnSp <- gBuffer(iucnSp, width=20*1000)
		iucnSp <- gBuffer(iucnSp, width=-20*1000)
		plot(iucnSp, col=alpha('forestgreen', 0.5), add=TRUE)

		testSurveysSpEa <- testSurveysSpEa[order(testSurveysSpEa$status), ]
		
		plot(kdeClass, col=kdeCols, legend=FALSE, add=TRUE)
		plot(west2, border='gray40', add=TRUE)
		plot(west1, border='gray40', lwd=2, add=TRUE)
		plot(west2[west2$NAME_2 == 'Plumas', ], border='gray40', lwd=3, add=TRUE)

		# # cols <- rep(NA, nrow(testSurveys))
		# # cols[testSurveysSpEa$status == '0 long absence'] <- testLongTermAbsCol
		# # cols[testSurveysSpEa$status == '1 recent absence'] <- testShortTermAbsCol
		# # cols[testSurveysSpEa$status == '2 detected'] <- testPresCol
		
		# # pchs <- rep(NA, nrow(testSurveys))
		# # pchs[testSurveysSpEa$status == '0 long absence'] <- testLongTermPch
		# # pchs[testSurveysSpEa$status == '1 recent absence'] <- testShortTermPch
		# # pchs[testSurveysSpEa$status == '2 detected'] <- testPresPch
		
		points(trainPresSpEa, pch=trainPresPch, cex=0.5, col=trainPresCol)
		# points(testSurveysSpEa, pch=pchs, col=cols, cex=0.6)
		box()
		
		# legend
		legendBreaks(
			'bottomleft',
			inset=0.01,
			height=0.39,
			width=0.30,
			title='Density of previously-\nknown occurrences',
			titleAdj=c(0.5, 0.9),
			col=kdeCols,
			adjX=c(0.08, 0.225),
			adjY=c(0.35, 0.77),
			labels=c('\U2265min presence', paste0('\U2265', '5th percentile'), paste0('  \U2265', '10th percentile')),
			labAdjX=0.58,
			cex=0.57,
			boxBg=alpha('white', 0.8)
		)
		
		usr <- par('usr')
		width <- usr[2] - usr[1]
		height <- usr[4] - usr[3]
		x <- usr[1] + 0.01 * width
		y <- usr[3] + 0.15 * height
		
		legend(
			x,
			y,
			legend=c('Previously-known\noccurrence', 'IUCN range'),
			pch=c(trainPresPch, NA),
			col=c(trainPresCol, 'black'),
			fill=c(NA, 'darkolivegreen3'),
			border=c(NA, 'black'),
			bty='n',
			cex=0.57,
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
		text(x, y[1], labels=paste(size / 1000, 'km'), cex=0.7)

		# inset
		insetNam <- nam1[nam1$NAME_1 %in% c('California', 'Oregon', 'Baja California', 'Baja California Sur', 'Washington'), ]
		insetNam <- ext(insetNam)
		insetNam <- vect(insetNam, crs=crs(nam1))
		insetNam <- buffer(insetNam, 100 * 1000)
		insetNam <- ext(insetNam)
		insetNam <- vect(insetNam, crs=crs(nam1))

		par(fig = c(0.6, 1, 0.0, 0.65), bg='white', new=TRUE)
		plot(insetNam, col='white', axes=FALSE, bty='o')
		plot(crop(nam1, insetNam), col='gray80', lwd=0.1, add=TRUE)
		plot(crop(iucn, insetNam), col='forestgreen', lwd=0.2, add=TRUE)
		
		counties <- nam2[nam2$NAME_2 %in% gapCounties$NAME_2, ]
		counties <- ext(counties)
		foc <- vect(counties, crs=crs(nam2))
		plot(foc, lwd=2, add=TRUE)

		# north arrow
		x0 <- -2350000
		x1 <- x0 + 7500
		y0 <- 6005000 - 20000
		y1 <- y + 115000
		arrows(x0 = x0, y0 = y0, x1 = x1, y1 = y1, angle = 20, length = 0.1, lwd = 2, xpd = NA)
		text(x1 + 2000, y1 - 5000, labels = 'N', cex = 1, pos = 3, srt = -23, xpd = NA)

		title(sub=date(), cex.sub=0.3, outer=TRUE, line=-1)		
				
	dev.off()
