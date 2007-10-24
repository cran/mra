residuals.cjs <- function( object, type="pearson", ... ){

if( !is.null(object$residuals) ){
	return( object$residuals )
}

#	Compute expected cell means
if( !is.null(object$fitted) ){
	cellmeans <- object$fitted
} else {
	cellmeans <- predict.cr( cjsobj )
}

#	Compute residuals
hists <- object$histories
hists[ hists >=2 ] <- 1  # Don't need '=' of '>=' because only 0,1, and 2 present, but included it just to be sure we only have 0's and 1's
if( type == "deviance" ){
	resids <- -sqrt(2*abs( log(1-cellmeans) ))   # these are for hists == 0
	resids[ hists == 1 ] <- sqrt(2*abs(log(cellmeans[ hists == 1 ])))  # these are for hists == 1
} else {
	resids <- (hists - cellmeans) / sqrt(cellmeans*(1-cellmeans))
}

resids
}

