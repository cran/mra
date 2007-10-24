"F.cr.model.matrix" <-
function( capture, survival ){
    call <- match.call()
    contrasts <- NULL
    mf <- match.call(expand.dots = FALSE)
    mf$family <- mf$start <- mf$control <- mf$maxit <- NULL
    mf$model <- mf$method <- mf$x <- mf$y <- mf$contrasts <- NULL
    mf$... <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")


    sf <- mf
    sf$formula <- survival
    sf$survival <- NULL
    sf$capture <- NULL
    form <- as.character(formula(sf))[2]
    if( nchar(form)==1 & form == "1" ){
	# this is a model with intercept only
	survX <- matrix(1,1,1)
	surv.intercept <- TRUE
	ny <- 0
	sur.names <- character(0)
    } else {
    	#sf <- eval(sf, parent.frame())
	sf <- eval( sf )
    	mt <- attr(sf, "terms")
	if( exists("R.version") ){
		xvars <- as.character(attr(mt, "variables"))[-1]
	} else {
		xvars <- as.character(attr(mt, "variables"))
	}
	if ((yvar <- attr(mt, "response")) > 0) 
        	xvars <- xvars[-yvar]
	xlev <- if (length(xvars) > 0) {
        	xlev <- lapply(sf[xvars], levels)
	        xlev[!sapply(xlev, is.null)]
    	}

    	survX <- NA
    	survX <- if ( length(attr( mt, "order")) != 0 ) {
        	model.matrix(mt, sf, contrasts)
	} else {
		stop("Empty models not allowed")
	}
        #assign.col <- F.get.assign( survX )
	assign.col <- attr( survX, "assign" )

	if( sum( assign.col == 0 ) == 1 ){
		surv.intercept <- TRUE
	    	ny <- sum( unique(assign.col) > 0 )
	} else {
		surv.intercept <- FALSE
		ny <- length( unique(assign.col) )
	}	
        sur.names <- xvars
    }

    cf <- mf
    cf$formula <- capture
    cf$survival <- NULL
    cf$capture <- NULL
    form <- as.character(formula(cf))[2]
    if( nchar(form)==1 & form == "1" ){
	# this is a model with intercept only
	capX <- matrix(1,1,1)
	cap.intercept <- TRUE
	nx <- 0
	cap.names <- character(0)
    } else {
    	#cf <- eval(cf, parent.frame())
	cf <- eval(cf)
	mt <- attr(cf, "terms")
	if( exists("R.version") ){
		xvars <- as.character(attr(mt, "variables"))[-1]
	} else {
		xvars <- as.character(attr(mt, "variables"))
	}
    	if ((yvar <- attr(mt, "response")) > 0) 
        	xvars <- xvars[-yvar]
    	xlev <- if (length(xvars) > 0) {
        	xlev <- lapply(cf[xvars], levels)
        	xlev[!sapply(xlev, is.null)]
    	}


    	capX <- NA
    	capX <- if (length(attr( mt, "order")) != 0) {
        	model.matrix(mt, cf, contrasts)
	} else {
		stop("Empty models not allowed")
	}


    	#assign.col <- F.get.assign( capX )
	assign.col <- attr( capX, "assign" )
	if( sum( assign.col == 0 ) == 1 ){
		cap.intercept <- TRUE
	    	nx <- sum( unique(assign.col) > 0 )
	} else {
		cap.intercept <- FALSE
		nx <- length( unique(assign.col) )
	}	
    	cap.names <- xvars
    }


    ans <- list( 
	capX=capX, 
	survX=survX, 
	n.cap.covars=nx, 
	n.sur.covars=ny, 
	cap.intercept= cap.intercept, 
	sur.intercept = surv.intercept, 
	cap.vars = cap.names, 
	sur.vars = sur.names)
    ans
}

