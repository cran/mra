F.cr.model.avg <- function( fits=ls(pattern="^fit"), what="survival", fit.stat="qaicc" ){
#
#	Perform model averaging on a list of capture-recapture models.
#	Details: 
#	Original routine by Eric Regehr. 
#
 
 
#	Get all the statistics from the fitted objects
stats <- se.stats <- all.fit.stat <- good.fits <- NULL

#   Convert fits to a list if necessary
if( !is.list(fits) ){
    #   Assume fits is a character vector of fit names. 
    #   Pull the objects into a list
    fits.names <- fits
    fits <- NULL
    for( f in fits.names ){
    	fit <- get( f, pos=.GlobalEnv )	
        if( "cr" %in% class(fit) ){
    	   if( (fit$exit.code == 1) & (fit$cov.code == 0) & (fit$df > 0) ){ 
              fits <- c(fits, list(fit))
              names(fits)[ length(fits) ] <- f
           }
        } else {
            warning(paste("Object", f, "in fits is not a CR object and has been ignored."))
        }
    }        
} 


#   Find the dimensions
nan <- lapply( fits, function(x) x$aux$nan )
nan <- nan[[1]]
ns  <- lapply( fits, function(x) x$aux$ns )
ns <- ns[[1]]


#   Extract statistics
f.pull.stats <- function( fit, stat ){

    if( is.atomic( fit ) ){
        return( NA )  # fit is probably NA here
    } else {
        if( (fit$exit.code != 1) | (fit$cov.code != 0) | (fit$df == 0) ){
            return(NA)  # fit did not converge
        } else {
            return( fit[[stat]] )
        }
    }       
}

if( substring(what,1,1) == "s" ){
	x.name <- "s.hat"
    x.se   <- "se.s.hat"
} else if( substring(what, 1,1) == "c" ){
	x.name <- "p.hat"
    x.se   <- "se.p.hat"
} else if( substring(what, 1,1) == "n" ){
	x.name <- "n.hat"
    x.se   <- "se.n.hat"
} else {
    stop(paste("Invalid option. Cannot model average '", what, "'", sep=""))
}

stats    <- lapply( fits, f.pull.stats, stat=x.name )
se.stats <- lapply( fits, f.pull.stats, stat=x.se )

one.stat <- stats[[1]]

#   Put statistics into arrays
n.stats <- length(stats[[1]])
stats <- matrix( unlist(stats), length(stats), n.stats, byrow=T )
se.stats <- matrix( unlist(se.stats), length(se.stats), n.stats, byrow=T )



#   retreive the fit statistics
all.fit.stat <- lapply( fits, f.pull.stats, stat=fit.stat )
all.fit.stat <- unlist( all.fit.stat )
 
 

# 	At this point, all.fit.stat is a vector containing the fit statistics. 
#   stats is an array of size length(fits) X length(statistics vector) containing 
#	the statistics to average.  Average over rows, using weights. 

#	Calculate AIC weights (Burnam and Anderson 2002 page XX.):	
delta.AIC <- all.fit.stat - min( all.fit.stat, na.rm = TRUE )
wi.array <- exp( -0.5 * delta.AIC ) / sum( exp( -0.5 * delta.AIC ), na.rm = TRUE )


#	Calculate the model averaged real parameters and standard errors (Burnham and Anderson 2002 pages 150 and 162):
#	stats is n x m, w.array is  n x 1
wi.array <- matrix( wi.array, nrow(stats), ncol(stats) )
a1 <- stats * wi.array
theta.average <- apply( a1, 2, sum )

var.theta <- se.stats^2
a1 <- matrix( theta.average, nrow=nrow(stats), ncol=ncol(stats), byrow=TRUE )


a2 <- wi.array * sqrt( var.theta + ( stats - a1 )^2 ) 
se.theta.average <- apply( a2, 2, sum ) 




#This is the average conditional standard error among models (i.e.,
#does not include a variance component for model selection uncertainty). This information
#is useful in estimating how much of the overall unconditional standard error was 
#due to variation among models:

a2 <- wi.array * se.stats  
se.conditional.theta.average <- apply( a2, 2, sum ) 


AIC.table <- data.frame( names(fits), all.fit.stat, delta.AIC, wi.array[,1], stringsAsFactors=F )
names(AIC.table) <- c( "model", fit.stat, paste("delta.", fit.stat, sep=""), paste( fit.stat, ".weight", sep=""))
AIC.table <- AIC.table[ order(AIC.table[,fit.stat]), ]

if( substring(what,1,1) == "s" ){

	hat = matrix( theta.average, nan, ns )
	se.hat = matrix( se.theta.average, nan, ns) 
	se.hat.conditional = matrix( se.conditional.theta.average, nan, ns )
	mod.selection.proportion = matrix( (se.theta.average - se.conditional.theta.average) / se.theta.average, nan, ns )

    dim.nms <- dimnames(one.stat)
    dimnames(hat) <- dim.nms
    dimnames(se.hat) <- dim.nms
    dimnames(se.hat.conditional) <- dim.nms
    dimnames(mod.selection.proportion) <- dim.nms

	a1 <- list( fit.table = AIC.table, 
		s.hat = hat, 
		se.s.hat = se.hat, 
		se.s.hat.conditional = se.hat.conditional,
		mod.selection.proportion = mod.selection.proportion        
		)

} else if( substring(what, 1,1) == "c" ){

	hat = matrix( theta.average, nan, ns ) 
	se.hat = matrix( se.theta.average, nan, ns) 
	se.hat.conditional = matrix( se.conditional.theta.average, nan, ns )
	mod.selection.proportion = matrix( (se.theta.average - se.conditional.theta.average) / se.theta.average, nan, ns )

    dim.nms <- dimnames(one.stat)
    dimnames(hat) <- dim.nms
    dimnames(se.hat) <- dim.nms
    dimnames(se.hat.conditional) <- dim.nms
    dimnames(mod.selection.proportion) <- dim.nms

	a1 <- list( fit.table = AIC.table, 
		p.hat = hat, 
		se.p.hat = se.hat, 
		se.p.hat.conditional = se.hat.conditional,
		mod.selection.proportion = mod.selection.proportion
		)

} else if( substring(what, 1,1) == "n" ){

	mod.selection.proportion = (se.theta.average - se.conditional.theta.average) / se.theta.average

    dim.nms <- names(one.stat)
    names(theta.average) <- dim.nms
    names(se.theta.average) <- dim.nms
    names(se.conditional.theta.average) <- dim.nms
    names(mod.selection.proportion) <- dim.nms

	a1 <- list( fit.table = AIC.table, 
		n.hat =  theta.average, 
		se.n.hat = se.theta.average, 
		se.n.hat.conditional = se.conditional.theta.average, 
		mod.selection.proportion = mod.selection.proportion
		)

    a1$n.hat.lower <- a1$n.hat - 1.96*a1$se.n.hat
    a1$n.hat.upper <- a1$n.hat + 1.96*a1$se.n.hat
    a1$nhat.v.meth <- fits[[1]]$nhat.v.meth + 3
    a1$n.hat.conf <- 0.95
    a1$intervals <- fits[[1]]$intervals
    
    class(a1) <- c("nhat", "cr")

} 

	
a1 
}












