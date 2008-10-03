F.huggins.estim <- function(capture, recapture=NULL, histories, remove=FALSE, cap.init, recap.init,
    algorithm=1, cov.meth=1, nhat.v.meth=1, df=NA){

start.sec <- proc.time()[3]

run.date <- Sys.time()

if( missing(histories) ){
    stop("Capture histories matrix must be specified")
}
if( missing(capture) ){
    stop("Capture covariates must be specified")
}
if( missing(recapture) ){
    stop("Re-capture covariates must be specified")
}

if( length(union( unique(histories), c(0,1))) > 2 ) stop("Capture histories must consist of 0's and 1's only.")

# Remove rows of all zeros. These are errors, and we could stop, but I'll just remove.
zero.ind <- apply( histories, 1, sum ) == 0
if( any(zero.ind) ){
    stop(paste("Rows of all zeros not allowed in history matrix.", sum(zero.ind), "rows found."))
}


algorithm <- 1
hist.name <- deparse(substitute(histories))
cr.call <- match.call()
nan <- nrow( histories )
ns <- ncol( histories )

# return the X and Y matrix.  Model for initial captures is in 'capture'. 
# Model for recaptures is in 'recapture'. 
# After this, both matrices are huge NAN by (ns*nx) matrices.

capX <- F.3d.model.matrix( as.formula(capture), nan, ns )
cap.intercept <- attr(capX, "intercept") == 1
cap.names <- attr(capX, "variables")
nx <- length(cap.names)

if( !is.null(recapture) ){
    recapX <- F.3d.model.matrix( as.formula(recapture), nan, ns )
    recap.intercept <- attr(recapX, "intercept") == 1
    recap.names <- attr(recapX, "variables")
    ny <- length(recap.names)
} else {
    recapX <- rep(1,nan)
    recap.intercept <- FALSE
    recap.names <- "NULL"
    ny <- 0
}

# Rep out the remove vector, and convert to integer
if(length(remove) < nx){
    remove <- rep(remove, nx)
} else if (length(remove) > nx ){
    remove <- remove[1:nx]
}
remove.vec <- as.numeric( remove )

#   Set initial values if missing or short
if( missing(cap.init) ){
    cap.init <- rep(0,nx)
} else if(length(cap.init) < nx ){
    cap.init <- c(cap.init, rep(0, nx-length(cap.init)))
} 

if( missing(recap.init) ){
    recap.init <- rep(0,ny)
} else if(length(recap.init) < ny ){
    recap.init <- c(recap.init, rep(0, ny-length(recap.init)))
} 



#   Determine covariance method.  1 = numeric 2nd derivatives, 2 = inverse of optimization Hessian
if( !any( cov.meth == c(1,2) ) ){
    cat("Covariance method must be either 1 = numeric 2nd derivatives (default), or 2 = Hessian of optimazation\n")
    cat("Using 2nd derivative method (cov.meth = 1).\n")
    cov.meth <- 1 
}

#   Do the estimation, but first allocate room for answers
loglik <- deviance <- aic <- qaic  <- lower.ci <- upper.ci <- 0
parameters <- se.param  <- rep(0, nx+ny )
covariance <- matrix( 0, nx+ny, nx+ny )
p.hat <- se.p.hat <- c.hat <- se.c.hat <- matrix( 0, nan, ns )
exit.code <- cov.code <- 0
n.hat <- se.n.hat <- 0
if( is.na(df) ){
    df.estimated <- 1  # Have MRA estimate rank of var-covar matrix
} else {
    df.estimated <- 0  # Don't bother, df either set by user or will use nx+ny
}

cat( "Calling MRA DLL to maximize likelihood.  Please wait...\n")

ans <- .Fortran( "hugginsmodel", 
        as.integer(nan), 
        as.integer(ns), 
        as.integer(nx), 
        as.integer(ny),
        as.integer(histories),  
        as.integer(remove.vec),
        as.integer(algorithm), 
        as.integer(cov.meth), 
        as.integer(nhat.v.meth), 
        as.double(capX),
        as.double(recapX), 
        as.double(cap.init),
        as.double(recap.init), 
        as.double(loglik), 
        as.double(deviance), 
        as.double(aic),   
        as.double(parameters),
        as.double(se.param), 
        as.double(covariance), 
        as.double(p.hat), 
        as.double(se.p.hat), 
        as.double(c.hat), 
        as.double(se.c.hat),
        as.double(n.hat), 
        as.double(se.n.hat), 
        as.double(lower.ci),
        as.double(upper.ci),
        as.integer(exit.code), 
        as.integer(cov.code), 
        as.integer(df.estimated), 
        PACKAGE="mra" )

cat(paste("Returned from MRA. Details in MRA.LOG.\n", sep=""))

#   Set to place in ans of first output.  This makes it easier to add or subtract 
#   parameters to .Fortran call.  All outputs come after this.
out.ind <- 14

loglik     <- ans[[out.ind + 0]]
deviance   <- ans[[out.ind + 1]] 
aic        <- ans[[out.ind + 2]] 
parameters <- ans[[out.ind + 3]]
se.param   <- ans[[out.ind + 4]] 
covariance <- ans[[out.ind + 5]] 
p.hat      <- ans[[out.ind + 6]] 
se.p.hat   <- ans[[out.ind + 7]] 
c.hat      <- ans[[out.ind + 8]]
se.c.hat   <- ans[[out.ind + 9]] 
n.hat      <- ans[[out.ind + 10]]
se.n.hat   <- ans[[out.ind + 11]]
lower.ci   <- ans[[out.ind + 12]]
upper.ci   <- ans[[out.ind + 13]]
exit.code  <- ans[[out.ind + 14]]
cov.code   <- ans[[out.ind + 15]]
df.estimated <- ans[[out.ind + 16]]

# ----- Fortran sets missing standard errors < 0. Reset missing standard errors to NA.
se.param[ se.param < 0 ] <- NA

# ----- R does not preserve the matrix structures in .Fortran call.  Put matricies, 
#   which are now vectors, back to matricies.
covariance <- matrix( covariance, nrow=nx+ny ) 
p.hat      <- matrix( p.hat, nrow=nan )
se.p.hat   <- matrix( se.p.hat, nrow=nan )
c.hat      <- matrix( c.hat, nrow=nan )
se.c.hat   <- matrix( se.c.hat, nrow=nan )



# ----- Work out exit codes.  Code is returned by VA09AD.
if( exit.code==0 ){
    exit.mess = "FAILURE: Initial Hessian not positive definite"
} else if( exit.code == 1 ){
    exit.mess = "SUCCESS: Convergence criterion met"
} else if( exit.code == 2 ){
    exit.mess = "FAILURE: G'dX > 0, rounding error"
} else if( exit.code == 3 ){
    exit.mess = "FAILURE: Likelihood evaluated too many times"
} else if( exit.code == -1 ){
    exit.mess = "FAILURE: Algorithm 2 not implimented yet.  Contact Trent McDonald."
} else {
    exit.mess = "Unknown exit code"
}

if(algorithm == 1){
    message <- "Optimization by VA09AD."
} else {
    message <- "Unknown optimization routine."
}

if(cov.meth == 1){
    cov.mess = "Covariance from numeric derivatives."
} else if (cov.meth == 2){
    cov.mess = "Covariance from optimization Hessian."
} else {
    cov.mess = "Unkown covariance method."
}

if(cov.code == 0){
    cov.mess = paste( cov.mess,  "SUCCESS: Non-singular covariance matrix.")
} else if( cov.code == 1) {
    cov.mess = paste( cov.mess,  "ERROR: COVARIANCE MATRIX IS SINGULAR.")
} else {
    cov.mess = paste( cov.mess,  "ERROR: NON-NEGATIVE DEFINITE COVARIANCE MATRIX.")
}

# ----- Remove trailing blanks from message
message <- c( paste( message, exit.mess), cov.mess)
cat(paste( "\t(", message, ")\n", sep="" ))

# ----- Wipe out the first recapture probability.  Can't have recapture in first period.
c.hat[,1] <- NA
se.c.hat[,1] <- NA

# ----- Transfer over coefficients
capcoef <- parameters[1:nx]
se.capcoef <- se.param[1:nx]
names(capcoef) <- cap.names
names(se.capcoef) <- cap.names
nms <- paste("cap.",names( capcoef ),sep="")

if( ny >= 1 ){
    recapcoef <- parameters[(nx+1):(nx+ny)]
    se.recapcoef <- se.param[(nx+1):(nx+ny)]
    names(recapcoef) <- recap.names
    names(se.recapcoef) <- recap.names
    nms <- c( nms, paste("recap.",names( recapcoef ), sep="") )
} else {
    recapcoef <- NULL
    se.recapcoef <- NULL
}    

dimnames(covariance) <- list( nms, nms )

# ----- Fix up number of parameters.
if( is.na(df) ){
    df <- df.estimated   # use the rank estimated by MRAWIN
} else if( df <= 0 ){
    df <- nx+ny  # assume full rank
} # else use the values supplied by user (unchanged from input)

# ----- Now that df is computed, recompute fit statistics
aic <- -2*loglik + 2*df
n.eff <- nan * ns
aicc <- aic + (2*df*(df+1))/(n.eff - df - 1)



# ----- Put all the "auxillary" info together
aux <- list( call=cr.call, nan=nan, ns=ns, nx=length(capcoef), ny=length(recapcoef),
        cov.name=nms, 
        ic.name=hist.name, 
        mra.version=packageDescription("mra")$Version, 
        R.version=R.version.string,
        run.date=run.date )


# ----- Estimates of N are computed in Fortran.  No need to modify.

# ----- Done. Put into a 'hug' object.
ans <- list( histories=histories, 
    aux=aux, 
    loglike=loglik, 
    deviance=deviance, 
    aic=aic, 
    aicc=aicc, 
    capcoef=capcoef, 
    se.capcoef=se.capcoef, 
    recapcoef=recapcoef, 
    se.recapcoef=se.recapcoef,
    remove=remove,
    covariance=covariance,
    p.hat=p.hat, 
    se.p.hat=se.p.hat, 
    c.hat=c.hat, 
    se.c.hat=se.c.hat,
    df=df, 
    message=message, 
    exit.code=exit.code, 
    cov.code=cov.code, 
    cov.meth=cov.meth, 
    n.hat = n.hat, 
    se.n.hat = se.n.hat, 
    n.hat.lower = lower.ci, 
    n.hat.upper = upper.ci,
    n.hat.conf = 0.95, 
    nhat.v.meth = nhat.v.meth, 
    num.caught=nan,
    n.effective=n.eff)
class(ans) <- c("hug", "cr")




ex.time <- (proc.time()[3] - start.sec) / 60
if( ex.time < 0.01666667 ){
    cat("\t(Execution time < 1 second)\n")
} else {
    cat(paste("\t(Execution time =", round(ex.time,2), "minutes)\n"))
}


ans

}
