F.fit.table <- function( fits=ls(pattern="^fit"), rank.by= "qaicc" ){
#
#	Compile a table of fit statistics for estimated models. 
#
#   fits = vector of fit object names
#   rank.by = the statistic to use to order the table by 
#

if(length(fits) == 0) stop("Vector of fitted objects must be specified.")

fits <- sort(fits)

n <- length(fits)

ans <- data.frame( model.num=rep(NA,n), 
    model.name=rep("", n), 
	converged = rep(NA, n),
	n.est.parameters=rep(NA,n), 
	n.coefficients=rep(NA,n), aicc=rep(NA,n), delta.aicc=rep(NA,n),
	aicc.wgt=rep(NA,n), qaicc = rep(NA,n), delta.qaicc=rep(NA,n), qaicc.wgt=rep(NA,n) )

ans$model.name <- as.character(ans$model.name)

for( i in 1:n ){
	fitname <- fits[i]
	fit <- get( fitname, pos=.GlobalEnv )
	ans$model.num[i] = i
	ans$model.name[i] = fitname
	if( (fit$exit.code == 1) & (fit$cov.code == 0) & (fit$df > 0) ){
		ans$converged[i] <- TRUE
		ans$aicc[i] = fit$aicc
		ans$qaicc[i] = fit$qaicc
		ans$n.est.parameters[i] = fit$df 
		ans$n.coefficients[i] = length( fit$parameters )
	} else {
		ans$converged[i] <- FALSE
		ans$aicc[i] = Inf
		ans$qaicc[i] = Inf
		ans$n.est.parameters[i] = NA 
		ans$n.coefficients[i] = length( fit$parameters )
	}		
}

ans$delta.aicc  <- ans$aicc - min(ans$aicc, na.rm=TRUE)
ans$delta.qaicc <- ans$qaicc - min(ans$qaicc, na.rm=TRUE)

ans$aicc.wgt <- exp(-.5*ans$delta.aicc)
ans$aicc.wgt <- ans$aicc.wgt / sum(ans$aicc.wgt, na.rm=TRUE)

ans$qaicc.wgt <- exp(-.5*ans$delta.qaicc)
ans$qaicc.wgt <- ans$qaicc.wgt / sum(ans$qaicc.wgt, na.rm=TRUE)

ord <- order( ans[,rank.by] )
ans <- ans[ord,]

ans
}


