"plot.cjs" <- 
function( x, type="n", ci=TRUE, smooth=TRUE, occasions=-1, animals=-1,
	smubass=5, ... ){

	if( type == "s" ){
		#  Plot survival
		y <- x$s.hat
		occasion <- 1:ncol(y)
		if( any(occasions <= 0) ){ 
			occasions <- 1:ncol(y)
		}
		occasions <- occasions[ occasions != ncol(y) ]
		if( is.null(dimnames(y)[[2]]) ){
			nms <- 1:ncol(y)
		} else {
			nms <- dimnames(y)[[2]]
		}
		y <- y[,occasions]
		occasion <- occasion[occasions]
		nms   <- nms[occasions]
		if( (length(animals) != 1) | (animals[1] > 0 ) ){
			y <- y[animals,]
		}


		survival <- y
		if( is.matrix(y) ){
			if( nrow(y) > 1 ){ 
				survival <- t( y )
			}
		}


		matplot(occasion, survival, type="l", xaxt="n", xlab="Period", 
		ylab="Survival", ...)
		axis(1, at=occasion, label=nms )
	
		ans <- survival
	} else {
		# Plot N hat
		y <- x$n.hat
		occasion <- seq( along=y )
		ns <- length( y )
		if( length(names(y)) != ns ) names(y) <- 1:ns
		nms <- names(y) 
		if( any(occasions <= 0) ){ 
			occasions <- 1:length(y)
		}
		occasions <- occasions[ occasions != 1 ]

		y <- y[occasions]
		occasion <- occasion[occasions]
		nms   <- nms[occasions]

		n.hat <- y
		if( ci ){
			lower.ci <- x$n.hat.lower[occasions]
			upper.ci <- x$n.hat.upper[occasions]
	
			#if( ylim[1] == -1 ){
			#	myylim <- range(c(upper.ci, lower.ci), na.rm=T)
			#} else {
			#	myylim <- ylim
			#}

			plot( range(occasion), range(n.hat,lower.ci,upper.ci), 
				type="n", xaxt="n", xlab="Occasion", ylab="N estimate", ...)
			axis( side=1, at=occasion, labels=nms )
			lines(occasion, lower.ci, type="l", lty=2)
			lines(occasion, upper.ci, type="l", lty=2)
			lines(occasion, n.hat, type="b")
		} else {

			#if( ylim[1] == -1 ){
			#	myylim <- range(n.hat)
			#} else {
			#	myylim <- ylim
			#}
			#
			#n.hat[ n.hat < min(myylim) ] <- min(myylim)
			#n.hat[ n.hat > max(myylim) ] <- max(myylim)
	
			plot(occasion, n.hat, type="b", xaxt="n",  ...)
			axis(1, at=occasion, label=nms )

		}


		if( smooth & exists("supsmu") ){
			sm <- supsmu( occasion, n.hat, bass=smubass )
			lines( sm, lwd=3, col=3 )
		} else {
			sm <- NA
		}

		ans <- sm
	}


	invisible( ans )
}

