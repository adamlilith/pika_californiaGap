# source('E:/Ecology/Drive/Research/Pikas - California Gap (Erik Beever et al) V2/Code/TEMP.r')


	xlab <- paste0('PC1 (', round(100 * pca$sdev[1]^2 / sum(pca$sdev^2), 1), '%)')
	ylab <- paste0('PC2 (', round(100 * pca$sdev[2]^2 / sum(pca$sdev^2), 1), '%)')
	png('./Figures & Tables/PCA with Test Classes.png', width=1800, height=1800, res=600)

		par(cex.axis=0.5, cex.lab=0.6, mar=rep(1.8, 4), oma=rep(0, 4), lwd=0.6, bty='n', mgp=c(0.7, -0, 0), tck=-0.01)
		
		smoothScatter(pca$scores[ , 1:2], pch=16, nrpoints=0, xlab=xlab, ylab=ylab)
		points(trainPres[ , pcs], pch=3, col=alpha('black', 0.4), cex=0.5 * cex)
		
		pres <- testSurveys[testSurveys$status == '2 detected', ]
		recentAbs <- testSurveys[testSurveys$status == '1 recent absence', ]
		longTermAbs <- testSurveys[testSurveys$status == '0 long absence', ]
		
		points(longTermAbs[ , pcs], pch=25, bg='red', cex=cex)
		points(recentAbs[ , pcs], pch=23, bg='yellow', cex=cex)
		points(pres[ , pcs], pch=21, bg='green', cex=cex)

		mult <- 4
		x0 <- 10.3
		y0 <- -3
		for (i in 1:9) {
		
			arrows(x0, y0, x0 + loads[i, 1] * mult, y0 + loads[i, 2] * mult, angle=20, length=0.07, xpd=NA)
			
			label <- rownames(loads)[i]
			if (label == 'summer respite') {
				ydelta <- + 0.2
			} else {
				ydelta <- 0
			}
			
			text(x0 + loads[i, 1] * mult, y0 + loads[i, 2] * mult + ydelta, labels = label, cex=0.5, xpd=NA)
		}

		legend('topleft', inset=c(0.01, 0.1), legend=c('Background', 'Previously-known occurrence', 'Occupied', 'Formerly occupied', 'No evidence'), pch=c(NA, 3, 21, 23, 25), fill=c(blues9[6], NA, NA, NA, NA), col=c(NA, 'black', 'black', 'black', 'black'), pt.bg = c(NA, NA, 'green', 'yellow', 'red'), border=c('black', NA, NA, NA, NA), bty='n', bg=NA, cex=cex, xpd=NA)


	dev.off()