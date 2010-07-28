\name{F.cjs.simulate}
\alias{F.cjs.simulate}

\title{
F.cjs.simulate - Generation of capture histories that follow a CJS model.
}

\description{
This function will generate a large number of capture history matrices that 
exactly follow an open Cormack-Jolly-Seber (CJS) model.  A super-population 
approach is taken wherein individuals with unique capture and survival 
probabilites that follow a specified model are randomly 'born' into the 
realized population and captured. Any CJS model, including those with 
heterogeneity, can be 
simulated using this approach.  
}

\usage{
F.cjs.simulate(super.p, super.s, fit, N1 = 1000, birth.rate = "constant.popln", R = 100)
}

\arguments{
  \item{super.p}{A matrix or vector of true capture probabilities in the super-population of individuals.
         If \code{super.p} is a vector, all individuals in the realized population will have the 
         same true capture probabilities, specified in \code{super.p}.  If \code{super.p} is a 
         vector, it is assumed that \code{length(super.p)} capture occasions will 
         be simulated.  If \code{super.p} is a matrix, the rows of \code{super.p} will be 
         randomly selected and used for the capture probabilities of individuals in the 
         realized population. That is, when animals are 'born' into the population, 
         a random row from \code{super.p} is picked and associated with the newly 
         'born' animal as it's vector of capture probabilities, when it survives. 
         If \code{super.p} is a matrix, it is assumed that \code{ncol(super.p)} 
         capture occasions will be simulated.  Number of rows in \code{super.p} can be anything
         greater than or equal to 1, and does not need to match number of rows in \code{super.s}. }

  \item{super.s}{A matrix or vector of true survival probabilities in the super-population of individuals.
         If \code{super.s} is a vector, all individuals in the realized population will have the 
         same true survival probabilities, specified in \code{super.s}.  If \code{super.p} is a matrix, 
         the rows of \code{super.p} will be 
         randomly selected and used as survival probabilities for individuals in the 
         realized population. That is, when animals are 'born' into the population, 
         a random row of \code{super.s} will be picked and associated with the newly 
         'born' animal as it's vector of survival probabilities.  The number of 
         capture occasions to simulate is set by the length or number of columns of \code{super.p}
         (see above).  Assuming number of occasions to simulate is NS, \code{super.s} must either 
         be of length \eqn{NS - 1} (if a vector) or have \eqn{NS - 1} columns (if a matrix).
         Number of columns in \code{super.s} is one less than NS because 
         survival probabilities apply between capture occasions. The vector
         \code{super.s[,j]} is the set of survival probabilities for animals alive just after 
         occasion \code{j} until just before occasion \code{j+1}. Number of rows in \code{super.s} can be anything
         greater than or equal to 1, and does not need to match number of rows in \code{super.p}.}

  \item{fit}{A previously estimated CJS object.  Instead of specifying \code{super.p} and \code{super.s} above,
        a fitted CJS model can be specified.  If either one of \code{super.p} or \code{super.s} is missing, 
        the estimated capture probabilities (\code{fit$p.hat}) and survival probabilities (\code{fit$s.hat}) will 
        be used as \code{super.p} and \code{super.s}, respectively.  Because capture probabilities for the 
        first occasion are not estimable by CJS models (unless you extrapolate using covariates), capture 
        probabilities for the first occasion are set equal to 1.0.  All members of the realized population 
        will be observed on the first occasion in this case. }
        
  \item{N1}{A scalar specifying the initial population size.  I.e., \code{N1} individuals will be 
        'born' into the realized population just before the first sampling occasion.}
  
  \item{birth.rate}{Either a vector of birth rates, or the string "constant.popln" (the default).  If 
        \code{birth.rate} = "constant.popln", the number of births between each capture occasion 
        will equal the number of deaths between each occasion.  In this way,  
        true population size will be (exactly) constant through time.  If \code{birth.rate} is 
        a vector of length \eqn{NS - 1}, then round( \eqn{N_j}{N(j)}*\code{birth.rate[,j]} ) births 
        will occur between occasions \code{j} and \code{j+1}, where \eqn{N_j}{N(j)} is the true number of 
        individuals in the population at occasion \code{j}.  Values in \code{birth.rate} can 
        be any non-negative number.  All animals in the population have one offspring  
        between occasions if \code{birth.rate = rep(1,NS)}}.  
  
  \item{R}{A scalar specifying the number of replications for the simulation. A total of \code{R} 
        independent capture history matricies will be generated. }
}

\details{Some examples: A two-group heterogeneity model has one group of individuals that have one common
        set of capture probability, and another group of individuals that have another set of common 
        capture probabilities.  A closed (no births or deaths) two-group heterogeneity model with 10 
        occasions can be simulated 
        using 
        \itemize{ \item F.cjs.simulate( rbind( rep(0.25,10),rep(0.75,10) ), rep(1,9) ). }
        One group of individuals will have probability of capture equal to 0.25, the other group will 
        have probability of capture equal to 0.75, and there will be approximately the same number of 
        individuals in both groups. Unequal groups can be simulated using        
        \itemize{ \item F.cjs.simulate( rbind( matrix(0.25,1,10),matrix(0.75,9,10) ), rep(1,9) ). }
        Using this call, approximatley 10\% of the individuals will have capture probabilities 
        equal to 0.25, while 90\% will have capture probabilities equal to 0.75. Additional 
        groups can be included by including more rows with distinct probabilities for \code{super.p}.
        
        A population with capture heterogeneity proportional to a vector \code{w} can be simulated using 
        \itemize{ \item F.cjs.simulate( matrix( w/sum(x), 100, 10), rep(1,9) ) }.  
    
        A stocastic population size with constant expected size can be simulated with a statement like 
        \itemize{ \item F.cjs.simulate( rep(0.25,10), rep(.5,9), birth.rate=rep(1,9) ). }
        Or, in general, setting survival to Phi and birth rates to (1-Phi)/Phi. 

        Because of the super-population approach taken here, it is not possible to specify which individuals 
        have which survival or capture probabilities, nor to guarentee that a certain number of individuals 
        in the realized population have capture probabilites equal to any particular value.  

}


\value{A list of length \code{R}.  Each component of this list is a list of length 2. 
        Each of the \code{R} sublists contains the following components:  

    \item{hists}{The simulated capture histories for a particular iteration. This is 
        a matrix with a random number of rows (due to stocastic nature of captures and 
        removal of un-observed individuals) and NS columns.}
    \item{popln.n}{A vector of length NS containing the true population sizes at each 
        sampling occasion.}}

\author{Trent McDonald, WEST Inc. (tmcdonald@west-inc.com)}



\seealso{
\code{\link{F.cjs.estim}}
}

\examples{

\donttest{  # dont test because it takes too long.  Code works, however.

## Simulate constant model, and analyze

ns <- 10
N <- 100
sim.list <- F.cjs.simulate( rep(0.3,ns), rep(0.9,ns-1), N1=N, R=100 )  # should do at least R = 500, but takes a while...

f.analyze <- function(x){
    fit <- F.cjs.estim( ~1, ~1, x$hists, control=mra.control(maxfn=200, cov.meth=2) )
    if( fit$exit.code == 1 ){
        return( fit$n.hat )
    } else {
        return( rep(NA,ncol(x$hists)) )
    }
}
results <- t(sapply(sim.list, f.analyze))
plot( 1:10, colMeans(results, na.rm=TRUE), xlab="Occasion", ylab="Mean population estimate", col="red", type="b")
abline( h=N )

## Plot RMSE by occasion
std <- apply(results, 2, sd, na.rm=TRUE)
bias <- apply(results - N, 2, mean, na.rm=TRUE)
plot( std, bias, type="n" )
text( std, bias, 2:10 )
abline(h=0)
title(main="RMSE by Sample Occasion")

## Show bias when heterogeneity is present
sim.list <- F.cjs.simulate( matrix(c(0.3,.7,.7,.7),4,ns), rep(0.9,ns-1), N1=N, R=100 ) 
results <- t(sapply(sim.list, f.analyze))
mean.N <- colMeans(results, na.rm=TRUE)
plot( 1:length(mean.N), mean.N, ylim=range(c(mean.N,N),na.rm=TRUE), xlab="Occasion", ylab="Mean population estimate", col="red", type="b")
abline( h=N )
abline( h=mean(mean.N), col="red", lty=2)
title(main="Heterogeniety causes negative bias")


## Simulate CJS model, first estimate one
data(dipper.histories)
ct <- as.factor( paste("T",1:ncol(dipper.histories), sep=""))
attr(ct,"nan")<-nrow(dipper.histories)
dipper.cjs <- F.cjs.estim( ~tvar(ct,drop=c(1,2)), ~tvar(ct,drop=c(1,6,7)), dipper.histories )

## Now generate histories from it.
sim.list <- F.cjs.simulate( fit=dipper.cjs, N1=100, birth.rate=rep(1,6), R=100 )

## Now analyze generated histories using lapply or sapply.  Can fit any model. 
## Here we fit the correct model.
f.analyze <- function(x){
    #   write a counter to console, this is not necessary
    i <- get("i", env=.GlobalEnv) + 1
    cat(paste("Iteration", i, "\n"))
    assign("i",i,env=.GlobalEnv)
    
    ct <- as.factor( 1:ncol(x$hists) )
    fit <- F.cjs.estim( ~tvar(ct,nan=nrow(x$hists),drop=c(1,2)), ~tvar(ct,nan=nrow(x$hists),drop=c(1,6,7)), 
        x$hists, control=mra.control(maxfn=200, cov.meth=2) )
    if( fit$exit.code == 1 ){
        return( fit$n.hat )
    } else {
        return( rep(NA,ncol(x$hists)) )
    }
}
i <- 0
results <- t(sapply(sim.list, f.analyze))
mean.N <- colMeans(results, na.rm=TRUE)
plot( 1:length(mean.N), mean.N, ylim=range(c(mean.N,N),na.rm=TRUE), xlab="Occasion", ylab="Mean population estimate", col="red", type="b")
abline( h=N )
title(main="Time varying CJS model")
}

}

\keyword{ survival }
\keyword{ models }% __ONLY ONE__ keyword per line
