F.cr.model.avg <- function( fits=ls(pat="^fit"), what="survival", fit.stat="qaicc" ){
#
#	Perform model averaging on a list of capture-recapture models.
#
#	Input:
#	fits = vector listing the names of fitted objects to be included in the averaging.  E.g., 
#		c("fit1","fit2","fit3"), or ls(pat="fit"). All these fitted objects, of type "cr", 
#		are assumed to reside in the 1st search position.
#	what = the statistic to model average.  Choices are "survival", "capture", "n.hat", and "beta"
#	fit.stat = fit statistic that should be used to compute weights.  I.e., if fit.stat = "aic", 
#		then all models in list are ordered by AIC, and weights are based on delta AIC values.
#		Valid values are any fit statistic available in the fitted CR object (i.e., aic, aicc, 
#		qaic, or qaicc )
#
#	Output: 
#	A "cr" object with model averaged estimates, conditional (on model selection) variances, and 
#	unconditional variances.  This object is suitable for plotting with plot.cr.
#	
#	Details: 
#	Original routine by Eric Regehr. 
#
 
 
#	Get all the statistics from the fitted objects
stats <- se.stats <- all.fit.stat <- good.fits <- NULL

for( f in fits ){
	fit <- get( f, pos=.GlobalEnv )	
	if( (fit$exit.code != 1) | (fit$cov.code != 0) | (fit$df == 0) ) next

	if( substring(what,1,1) == "s" ){
		nan      <- nrow(fit$s.hat)
		ns       <- ncol(fit$s.hat)
		stats    <- rbind( stats, c(fit$s.hat) )
		se.stats <- rbind( se.stats, c(fit$se.s.hat) )
	} else if( substring(what, 1,1) == "c" ){
		nan      <- nrow(fit$p.hat)
		ns       <- ncol(fit$p.hat)
		stats    <- rbind( stats, c(fit$p.hat) )
		se.stats <- rbind( se.stats, c(fit$se.p.hat) )
	} else if( substring(what, 1,1) == "n" ){
		stats    <- rbind( stats, c(fit$n.hat) )
		se.stats <- rbind( se.stats, c(fit$se.n.hat) )
	} else {
        stop(paste("Invalid option. Cannot model average '", what, "'", sep=""))
	}



	all.fit.stat <- rbind( all.fit.stat, fit[[fit.stat]] )
	good.fits <- c(good.fits, f)

}

# 	At this point, all.fit.stat is a matrix size length(list) X length(statistics vector) containing 
#	the statistics to average.  Average over rows, using weights. 

#	Calculate AIC weights (Burnam and Anderson 2002 page XX.):	
delta.AIC <- all.fit.stat - min( all.fit.stat, na.rm = TRUE )
wi.array <- exp( -0.5 * delta.AIC ) / sum( exp( -0.5 * delta.AIC ), na.rm = TRUE )


#	Calculate the model averaged real parameters and standard errors (Burnham and Anderson 2002 pages 150 and 162):
#	Note the automatic replication of w.array across columns.  stats is n x m, w.array is  n x 1
wi.array <- matrix( rep(wi.array, ncol(stats)), nrow(stats))
a1 <- stats * wi.array
theta.average <- apply( a1, 2, sum, na.rm = TRUE )

var.theta <- se.stats^2
a1 <- matrix( theta.average, nrow=nrow(stats), ncol=ncol(stats), byrow=TRUE )


a2 <- wi.array * ( var.theta + ( stats - a1 )^2 )^0.5 
se.theta.average <- apply( a2, 2, sum, na.rm = TRUE ) 




#This is the average conditional standard error among models (i.e.,
#does not include a variance component for model selection uncertainty). This information
#is useful in estimating how much of the overall unconditional standard error was 
#due to variation among models:

a2 <- wi.array * se.stats  
se.conditional.theta.average <- apply( a2, 2, sum, na.rm = TRUE ) 

AIC.table <- matrix( c( all.fit.stat, delta.AIC, wi.array[,1] ), nrow = length(delta.AIC), ncol = 3 )
dimnames(AIC.table) <- list( good.fits, c( fit.stat, paste("delta.", fit.stat, sep=""), paste( fit.stat, ".weight", sep="")))
AIC.table <- AIC.table[ order(AIC.table[,1]), ]


if( substring(what,1,1) == "s" ){
	a1 <- list( fit.table = AIC.table, 
		s.hat = matrix( theta.average, nan, ns ), 
		se.s.hat = matrix( se.theta.average, nan, ns), 
		se.s.hat.conditional = matrix( se.conditional.theta.average, nan, ns ),
		mod.selection.proportion = matrix( (se.theta.average - se.conditional.theta.average) / se.theta.average, nan, ns )
		)
    a1$s.hat <- a1$s.hat[ , -ns ]
    a1$se.s.hat <- a1$se.s.hat[ , -ns ]
    a1$se.s.hat.conditional <- a1$se.s.hat.conditional[ , -ns ]
    a1$mod.selection.proportion <- a1$mod.selection.proportion[ , -ns ]
} else if( substring(what, 1,1) == "c" ){
	a1 <- list( fit.table = AIC.table, 
		p.hat = matrix( theta.average, nan, ns ), 
		se.p.hat = matrix( se.theta.average, nan, ns), 
		se.p.hat.conditional = matrix( se.conditional.theta.average, nan, ns ),
		mod.selection.proportion = matrix( (se.theta.average - se.conditional.theta.average) / se.theta.average, nan, ns )
		)
    a1$p.hat <- a1$p.hat[ , -1 ]
    a1$se.p.hat <- a1$se.p.hat[ , -1 ]
    a1$se.p.hat.conditional <- a1$se.p.hat.conditional[ , -1 ]
    a1$mod.selection.proportion <- a1$mod.selection.proportion[ , -1 ]
} else if( substring(what, 1,1) == "n" ){
	a1 <- list( fit.table = AIC.table, 
		n.hat =  theta.average, 
		se.n.hat = se.theta.average, 
		se.n.hat.conditional = se.conditional.theta.average, 
		mod.selection.proportion = (se.theta.average - se.conditional.theta.average) / se.theta.average
		)
    a1$n.hat[1] <- NA
    a1$se.n.hat[1] <- NA
    a1$se.n.hat.conditional[1] <- NA
    a1$mod.selection.proportion[1] <- NA
    a1$n.hat.lower <- a1$n.hat - 1.96*a1$se.n.hat
    a1$n.hat.upper <- a1$n.hat + 1.96*a1$se.n.hat

} 

	
a1 
}












