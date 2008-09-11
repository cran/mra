F.cjs.estim <- function(capture, survival, histories, cap.init, sur.init, 
    group, algorithm=1, cov.meth=1, nhat.v.meth=1, 
    c.hat=-1.0, df=NA, intervals=rep(1,ncol(histories)-1), conf=0.95){

start.sec <- proc.time()[3]

run.date <- Sys.time()

if( missing(histories) ){
    stop("Capture histories matrix must be specified")
}
if( missing(capture) ){
    stop("Capture covariates must be specified")
}
if( missing(survival) ){
    stop("Survival covariates must be specified")
}

if( length(union( unique(histories), c(0,1,2))) > 3 ) stop("Capture histories must consist of 0's, 1's, and 2's only.")

if( !any( algorithm == c(1,2,3) ) ){
    warning(paste("Algorithm", algorithm, "not recognized. Algorithm 1 used"))
    algorithm <- 1
}

hist.name <- deparse(substitute(histories))
cr.call <- match.call()

nan <- nrow( histories )
ns <- ncol( histories )

if( length(intervals) < (ns-1)){
    stop(paste("Too few time intervals specified. INTERVALS vector should have length", ns-1))
} else if(length(intervals) >= (ns-1)){
    intervals <- intervals[1:(ns-1)]
    intervals <- c(intervals, 0)  # Make this vector length ns, but we never use the last element in estimation.
}


# ---- Get the X and Y matricies.  After this, the x and Y matricies are in huge NAN by (ns*nx) matricies.
covars <- F.cr.model.matrix( capture, survival, nan, ns )  

nx <- covars$n.cap.covars  #nx and ny include intercept.  This is total number of parameters
ny <- covars$n.sur.covars


#   Set initial values if missing or short
if( missing(cap.init) ){
    cap.init <- rep(0,nx)
} else if(length(cap.init) < (nx) ){
    cap.init <- c(cap.init, rep(0, nx-length(cap.init)))
} 

if( missing(sur.init) ){
    sur.init <- rep(0,ny)
} else if(length(sur.init) < (ny) ){
    sur.init <- c(sur.init, rep(0, ny-length(sur.init)))
} 


#   Fix up the group variable
if( missing( group )){
    group <- rep(1, nan)
    ng <- 1
} else {
    ng <- length(unique(group))
}



#   Determine covariance method.  1 = numeric 2nd derivatives, 2 = inverse of optimization Hessian
if( !any( cov.meth == c(1,2) ) ){
    cat("Covariance method must be either 1 = numeric 2nd derivatives (default), or 2 = Hessian of optimazation\n")
    cat("Using 2nd derivative method (cov.meth = 1).\n")
    cov.meth <- 1 
}


#   Transfer over c.hat
vif <- c.hat

#   Do the estimation, but first allocate room for answers
loglik <- deviance <- aic <- qaic <- chisq.vif <- df.vif <- 0
parameters <- se.param <- rep(0, nx + ny )
covariance <- matrix( 0, nx+ny, nx+ny )
p.hat <- se.p.hat <- s.hat <- se.s.hat <- matrix( 0, nan, ns )
exit.code <- cov.code <- 0
n.hat <- se.n.hat <- rep(0, ns)
if( is.na(df) ){
    df.estimated <- 1  # Have MRAWIN estimate rank of var-covar matrix
} else {
    df.estimated <- 0  # Don't bother, df either set by user or will use nx+ny
}

cat( "Calling MRA DLL to maximize likelihood.  Please wait...\n")

ans <- .Fortran( "cjsmod", 
        as.integer(nan), 
        as.integer(ns), 
        as.integer(nx), 
        as.integer(ny), 
        as.integer(ng), 
        as.integer(histories), 
        as.integer(group), 
        as.integer(algorithm), 
        as.integer(cov.meth), 
        as.integer(nhat.v.meth), 
        as.double(covars$capX), 
        as.double(covars$survX), 
        as.double(cap.init), 
        as.double(sur.init), 
        as.double(loglik), 
        as.double(deviance), 
        as.double(aic), 
        as.double(qaic), 
        as.double(vif), 
        as.double(chisq.vif), 
        as.double(df.vif), 
        as.double(parameters),
        as.double(se.param), 
        as.double(covariance), 
        as.double(p.hat), 
        as.double(se.p.hat), 
        as.double(s.hat), 
        as.double(se.s.hat), 
        as.double(n.hat), 
        as.double(se.n.hat), 
        as.integer(exit.code), 
        as.integer(cov.code), 
        as.integer(df.estimated), 
        as.double(intervals), 
        PACKAGE="mra" 
        ) 

#   Remember to add PACKAGE="mra" back into the .Fortran call above.

cat(paste("Returned from MRAWIN. Details in MRA.LOG.\n", sep=""))

loglik     <- ans[[15]]
deviance   <- ans[[16]] 
aic        <- ans[[17]] 
qaic       <- ans[[18]] 
vif        <- ans[[19]] 
chisq.vif  <- ans[[20]] 
df.vif     <- ans[[21]] 
parameters <- ans[[22]]
se.param   <- ans[[23]] 
covariance <- ans[[24]] 
p.hat      <- ans[[25]] 
se.p.hat   <- ans[[26]] 
s.hat      <- ans[[27]] 
se.s.hat   <- ans[[28]]
n.hat      <- ans[[29]]
se.n.hat   <- ans[[30]]
exit.code  <- ans[[31]]
cov.code   <- ans[[32]]
df.estimated <- ans[[33]]

# ----- Reset missing standard errors to NA
se.param[ se.param < 0 ] <- NA

# ----- R does not preserve the matrix structures in .Fortran call.  Put matricies, 
#   which are now vectors, back to matricies.
covariance <- matrix( covariance, nrow=nx+ny ) 
p.hat      <- matrix( p.hat, nrow=nan )
se.p.hat   <- matrix( se.p.hat, nrow=nan )
s.hat      <- matrix( s.hat, nrow=nan )
se.s.hat   <- matrix( se.s.hat, nrow=nan )



# ----- Work out exit codes
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

# Remove trailing blanks from message
message <- c( paste( message, exit.mess),
        cov.mess)
cat(paste( "\t(", message, ")\n", sep="" ))


# ----- Fix up capture and survival estimates
#   Wipe out the first capture probability and the last survival probability. These are computed 
#   by mrawin5 if there are covariates out there, but they are not used in the likelihood, and 
#   we can't really estimate them.
p.hat[,1] <- NA
s.hat[,ns] <- NA
se.s.hat[,ns] <- NA
se.p.hat[,1] <- NA

capcoef <- parameters[1:nx]
se.capcoef <- se.param[1:nx]
surcoef <- parameters[(nx+1):(nx+ny)]
se.surcoef <- se.param[(nx+1):(nx+ny)]

names(capcoef) <- covars$cap.vars
names(se.capcoef) <- covars$cap.vars
names(surcoef) <- covars$sur.vars
names(se.surcoef) <- covars$sur.vars


dimnames(covariance) <- list( c( paste("cap.",names( capcoef ),sep=""), paste("sur.",names( surcoef ), sep="")),
    c( paste("cap.",names( capcoef ),sep=""), paste("sur.",names( surcoef ), sep="")))

# ----- Fix up number of parameters.
if( is.na(df) ){
    df <- df.estimated   # use the rank estimated by MRAWIN
} else if( df <= 0 ){
    df <- nx+ny  # assume full rank
} # else use the values supplied by user (unchanged from input)

# ----- Now that df is computed, recompute fit statistics
aic <- -2*loglik + 2*df
qaic <- -2*(loglik/vif) + 2*df
aicc <- aic + (2*df*(df+1))/(nan - df - 1)
qaicc <- qaic + (2*df*(df+1))/(nan - df - 1)


# ----- LEAVE THIS CODE HERE FOR FUTURE REFERENCE 
# ----- Compute EBC of Peterson (1986), Stat & Prob letters, p227
#qf <- t(parameters) %*% covariance %*% parameters
#q1 <- qf/((df)*nan)
#q2 <- qf/(df)
#ebc <- -2*loglik + (df)*log(nan) + 
#    (df)*log( max(c( 1/nan, q1 ))) +
#    (df)*min( c( 1, q2 ) )
    
# ----- Put all the "auxillary" info together
aux <- list( call=cr.call, 
        nan=nan, 
        ns=ns, 
        nx=nx, 
        ny=ny, 
        cov.name=c(names(capcoef), names(surcoef)), 
        ic.name=hist.name, 
        mra.version=packageDescription("mra")$Version, 
        R.version=R.version.string,
        run.date=run.date )


# ----- Fix up the estimates of N.
names(n.hat) <- dimnames(histories)[[2]]
num.observed <- c( t( rep(1, nrow(histories)) ) %*% (histories >= 1) )
crit.value <- qnorm( 1-((1-conf)/2) )
lower.ci <- n.hat - crit.value * se.n.hat
lower.ci <- ifelse(lower.ci < num.observed, num.observed, lower.ci)
upper.ci <- n.hat + crit.value * se.n.hat

n.hat[1] <- NA  # can't estimate the first N
se.n.hat[1] <- NA
lower.ci[1] <- NA
upper.ci[1] <- NA


# ----- Done. Put into a 'CR' object.

ans <- list( histories=histories, 
    aux=aux, 
    intervals=intervals[1:(ns-1)], 
    loglike=loglik, 
    deviance=deviance, 
    aic=aic, 
    qaic=qaic, 
    aicc=aicc, 
    qaicc=qaicc,
    vif=vif, 
    chisq.vif=chisq.vif, 
    vif.df=df.vif, 
    parameters=parameters, 
    se.param=se.param, 
    capcoef=capcoef, 
    se.capcoef=se.capcoef, 
    surcoef=surcoef, 
    se.surcoef=se.surcoef, 
    covariance=covariance,
    p.hat=p.hat, 
    se.p.hat=se.p.hat, 
    s.hat=s.hat, 
    se.s.hat=se.s.hat, 
    df=df, 
    message=message, 
    exit.code=exit.code, 
    cov.code=cov.code, 
    cov.meth=cov.meth, 
    n.hat = n.hat, 
    se.n.hat=se.n.hat, 
    n.hat.lower = lower.ci, 
    n.hat.upper = ceiling(upper.ci),
    n.hat.conf = conf, 
    nhat.v.meth = nhat.v.meth, 
    num.caught=num.observed)
    
class(ans) <- c("cjs", "cr")



#   Compute fitted and residual components
ans$fitted <- predict( ans )
ans$residuals <- residuals( ans, type="pearson" )  
ans$resid.type <- "pearson"


ex.time <- (proc.time()[3] - start.sec) / 60
if( ex.time < 0.01666667 ){
    cat("\t(Execution time < 1 second)\n")
} else {
    cat(paste("\t(Execution time =", round(ex.time,2), "minutes)\n"))
}


ans

}
