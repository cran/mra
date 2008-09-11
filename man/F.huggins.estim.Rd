\name{F.huggins.estim}
\alias{F.huggins.estim}
\alias{huggins.estim}
\alias{huggins}
\alias{Huggins}

\title{ F.huggins.estim - Estimation of Huggins closed 
    population capture-recapture model. }

\description{
   Estimates Huggin's closed population capture-recapture models with individual, time, and
    individual-time varying covariates using the "regression" parameterization of 
    Amstrup et al (2006, Ch 9). For live recaptures only. 
    A logistic link function is used to relate probability of capture 
    to external covariates.  
}

\usage{
F.huggins.estim(capture, recapture, histories, cap.init, recap.init,
    algorithm=1, cov.meth=1, nhat.v.meth=1, df=NA)
}


\arguments{
  \item{capture}{ Formula specifying covariates to included in the initial 
    capture probability model.  Must be a formula object without 
    a response. Specify ~, followed by the names of 2-D arrays of covariates to relate to 
    initial capture probability. 
    For example: 'capture = ~ age + sex', where age and sex are matrices of size NAN X NS 
    containing the age and sex covariate values.  NAN = number of individuals = number of rows in 
    \code{histories} matrix (see below).  NS = number of samples = number of columns in histories
    matrix (see below). Number of matrices specified in the initial capture model is 
    assumed to be NX. Specify 
    intercept only model as 'capture = ~ 1'.  Specify model without an intercept 
    using 'capture = ~ -1 + x'.  }
  \item{recapture}{ Formula specifying covariates to included in the 
    recapture probability model.  Should be specified the same way as the 
    \code{capture} model.  
    For example: 'recapture = ~ behave + sex'. The number of covariates specified in
    the recapture model is NY.  Thus, total number of parameters this routine 
    attempts to estimate is NX+NY. See \code{df} argument. }
  \item{histories}{ A NAN X NS = (number of individuals) X (number of capture occasions) matrix 
    containing capture histories.  Capture histories are comprised of 0's and 1's only. 
    0 in cell (i,j) of this matrix means individual i was not captured on occasion j, 
    1 in cell (i,j) means individual 
    i was captured on occasion j and released live back into the population.  Because 
    the population 
    being sampled is assumed closed, deaths on capture (known removals) are not allowed. 
    If deaths on capture occured and an estimate of N at the beginning of the 
    study is sought, remove the entire history, estimate N using this routine from
    the remaining histories, and add back the number of deleted histories.}
  \item{cap.init}{ (optional) Vector of initial values for coefficients in the initial 
    capture model.  One element 
    per covariate in \code{capture}. This parameter does not usually need to be specified. }
  \item{recap.init}{ (optional) Vector of initial values for coefficients in the 
    recapture model.  One element 
    per covariate in \code{recapture}. This parameter does not usually need to be specified. }
  \item{algorithm}{ Integer specifying the maximization algorithm to use.  If \code{algorithm} = 1, 
    the VA09AD algorithm from the HSL library is used. The VA09AD algorithm is very reliable, 
    has been tested extensively (same algorithm as Program MARK), and will almost always find the 
    maximum likelihood estimates. This parameter was added to allow easy addition of other 
    algorithms in the future, but no other options are currently implemented because VA09AD 
    works so well. Anything other than 1 is ignored. }
  \item{cov.meth}{ Integer specifying the covariance estimation method. \code{cov.meth} = 1 
    takes numerical 2nd partial derivatives.  \code{cov.meth} = 2 inverts the Hessian 
    of the maximization.  Method 1 (numeric 2nd derivatives) is the preferred method, but
        is computationally expensive and therefore slower than method 2.  The difference in 
        speed is minimal when number of parameters < 15. Method 2 variances are 
        generally very close to method 1 variances, and could be used during initial model 
        fitting stages when exact estimation of variance may not be necessary.    }
  \item{nhat.v.meth}{ Integer specifying method for computing variance estimate 
    for the population size estimate. Currently, only \code{nhat.v.meth} = 1 is implemented. 
    Plans are for \code{nhat.v.meth} = 2 to be a boot strap estimate of variance.
    \code{nhat.v.meth} = 1 is a delta method estimator utilizing the derivative of 
    P(ever captured) w.r.t. the capture parameters.  This is the same estimator as 
    used in program MARK. }
  \item{df}{ External (override) model degrees of freedom to use during estimation.  
    If \code{df} == NA, the number of parameters is estimated from the rank of the 
    matrix of 2nd derivatives or Hessian, depending on \code{cov.meth} parameter.  
    If \code{df} <= 0, the number of parameters will be 
    set to NX+NY = the number of estimated coefficients. Otherwise, if \code{df} > 0, 
    the supplied value is used.  Only AIC, QAIC, AICc, and QAICc are dependent on 
    this value (in their penalty terms). }
}

\details{
    This routine compiles all the 
    covariate matrices, then calls a Fortran routine to maximize the 
    Huggins closed population likelihood. So-called heterogenious models 
    that utilize mixture distributions for probability of capture cannot be
    fitted via this routine.  
  
    See Details of \code{help(F.cjs.estim)} for ways that 2-d matrices, 1-d 
    vectors, and 1-d factors can be specified in the capture and recapture
    models.  
    
    A log file, named \code{mra.log}, is written to the current directory.  
    This file contains additional details in a semi-friendly 
    format.  CAREFUL: This file is overwritten each run. 
}


\value{
An object (list) of class c("hug","cr") with many components.
    Use \code{print.hug} to print
    it nicely. Use \code{names(fit)}, where the call was \code{fit <- F.huggins.estim(...)},
    to see names of all returned components. To see values of individual components,
    issue commands like fit\$n.hat, fit\$se.n.hat, etc.

    Components of the returned object are as follows:
  \item{histories}{The input capture history matrix. Size NAN x NS}
  \item{aux}{Auxiliary information, mostly stored input values.
    This is a list containing:
    \$call,
    \$nan=number of individuals,
    \$ns=number of samples,
    \$nx=number of coefficients in the initial capture model,
    \$ny=number of coefficients in recapture model,
    \$cov.name=names of the covariates in both models (initial capture covariates first, then recapture
    covariates),
    \$ic.name=name of capture history matrix,
    \$mra.version=version number of this package,
    \$R.version=R version used,
    \$run.date=date the model was estimated.}
  \item{loglik}{Value of the Huggins log likelihood at it's maximum.}
  \item{deviance}{Model deviance = -2*\code{loglik}.  This is relative deviance
    because all covariates are individual and time varying.  When individual
    covariates are present, a saturated likelihood cannot be computed.  Use this
    to compute deviance differences only.}
  \item{aic}{AIC for the model = \code{deviance} + 2*(df).  df is either the estimated
    number of independent parameters (by default), or NX+NY, or a specified value,
    depending on the input value of \code{df} parameter.}
  \item{aicc}{AIC with small sample correction = AIC + (2*\code{df}*(\code{df}+1))
    / (\code{NAN} - \code{df} - 1)}
  \item{capcoef}{Vector of estimated coefficients in the initial capture model.
    Length NX.}
  \item{se.capcoef}{Vector of estimated standard errors for coefficients in
    initial capture model. Length NX.}
  \item{recapcoef}{Vector of estimated coefficients in the recapture model.
    Length NY.}
  \item{se.surcoef}{Vector of standard errors for coefficients in recapture model.
    Length NY.}
  \item{covariance}{Variance-covariance matrix for the estimated model coefficients.
    Size (NX+NY) X (NX+NY).}
  \item{p.hat}{Matrix of estimated initial capture probabilities
    computed from the model.  Size of this matrix is NAN x NS.
    Cell (i,j) is estimated probability of first capture for individual i
    during capture occasion j.}
  \item{se.p.hat}{Matrix of standard errors for estimated initial capture probabilities.
    Size NAN x NS.}
  \item{c.hat}{Matrix of estimated recapture probabilities computed from the model.
    Size NAN x NS.  Cell (i,j) is estimated probability of capturing individual i
    during occasion j given that it was initially captured prior to j. }
  \item{se.c.hat}{Matrix of standard errors for estimated recapture probabilities.
    Size NAN X NS.  }
  \item{df}{Number of estimable parameters in the model. \code{df} is either
    the estimated number of independent parameters (by default) based on rank
    of the variance matrix,
    or NX+NY, or a specified value, depending on the input value of \code{df} parameter.}
  \item{message}{A string indicating whether the maximization routine converged.}
  \item{exit.code}{Exit code from the maximization routine.
    Interpretation for \code{exit.code} is in \code{message}.}
  \item{cov.code}{A code indicating the method used to compute the covariance matrix.}
  \item{cov.meth}{String indicating method used to compute covariance matrix.
    Interprets \code{cov.code}.}
  \item{n.hat}{The Huggins estimate of population size.  This estimate is
    sum( 1/ pstar(i) ), where pstar(i) is probability of observing individual i,
    which equals 1 - p.hat[i,1]*p.hat[i,2]* ... *p.hat[i,NS], where p.hat is the
    returned value of \code{p.hat}.}
  \item{se.n.hat}{Estimated standard error of \code{n.hat}.  Computed using method
    specified in \code{nhat.v.meth}.}
  \item{n.hat.lower}{Lower limit of log based confidence interval for
    \code{n.hat}.  }
  \item{n.hat.upper}{Upper limit of log based confidence interval for
    \code{n.hat}.}
  \item{n.hat.conf}{Confidence level for the interval on \code{n.hat}.
    Currently, confidence level cannot be changed
    from 95\%. }
  \item{nhat.v.meth}{Code for method used to compute variance of \code{n.hat}.  Currently,
    this is 1 only.}
  \item{num.caught}{Number of individuals ever captured = number of
    rows in the \code{histories} matrix.}
  \item{n.effective}{Effective sample size = number of observed individuals times
    number of occasions = NAN * NS}

}    


\references{ 
Huggins, R. M. 1989. On the statistical analysis of capture experiments. Biometrika 76:133-140.

Amstrup, S. C., T. L. McDonald, and B. F. J. Manly (editors). 2005. Handbook of 
    Capture-Recapture Analysis. Princeton University Press.

 }
\author{ Trent McDonald, WEST-INC, tmcdonald@west-inc.com }

\seealso{ \code{\link{print.hug}}, 
  \code{\link{F.cjs.estim}} }

\examples{

\donttest{
# Fit time-constant model (i.e., model M(0)).
ch.mat <- matrix( round(runif(30*5)), nrow=30, ncol=5)
ch.mat <- ch.mat[ apply(ch.mat,1,sum) > 0, ]  # no zero rows allowed
hug.m0 <- F.huggins.estim( ~1, ~1, ch.mat )

# Fit time varying model (i.e., model M(t))
ct <- as.factor(1:ncol(ch.mat))
attr(ct,"nan") <- nrow(ch.mat)
hug.mt <- F.huggins.estim( ~tvar(ct), ~tvar(ct), ch.mat)

}

}
\keyword{ survival }% at least one, from doc/KEYWORDS
\keyword{ models }% __ONLY ONE__ keyword per line
