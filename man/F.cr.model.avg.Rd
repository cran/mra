\name{F.cr.model.avg}
\alias{F.cr.model.avg}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{ F.cr.model.avg - Model averaging of mark-recapture parameters. }

\description{
    Computes model averaged estimates of survival, capture probability, or population size
    estimates from a set of previously fitted MRA objects.
}

\usage{
F.cr.model.avg( fits=ls(pat="^fit"), what="survival", fit.stat="qaicc" )
}

%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{fits}{A character vector of MRA fitted object names to include in the
    model averaged estimates. The default value
    will use any object whose name starts with "fit" in the working directory (.GlobalEnv).
    An an example, if fitted objects are named
    "fit1.01", "fit1.02", and "fit1.03", \code{fits} should equal \code{c("fit1.01", "fit1.02", "fit1.03")},
    or \code{ls(pat="^fit1")}, assuming no other objects in the working directory start with "fit1".}
  \item{what}{A text string naming the parameter of model fits to average.  Choices are "survival" (the 
    default), "capture", and "n.hat".  Only the first character is inspected (e.g., "c" is 
    equivalent to "capture"). }  
  \item{fit.stat}{A string (scalar) naming the model fit statistic to use when computing 
    model weights.  Possible values
    are: "qaicc" (the default), "qaic", "aicc", and "aic".}
}

\details{
    A rudimentary check for convergence is done on each fitted model.  If the routine determines
    that a model did not converge, the model is excluded from the model averaged estimates. 
    The test for whether a model converged is
    \code{(fit$exit.code == 1) & (fit$cov.code == 0) & (fit$df > 0)}, where \code{fit} is
    the fitted object.

    Conditional and unconditional variance estimates are computed following Burnham and 
    Anderson 2002 (pages 150 and 162 and surrounding). 
    
    The returned object is suitable for plotting using \code{plot.cr}. 


}

\value{A data frame, sorted by \code{fit.stat.to.use}, with the following columns

  \item{model.num}{Model number assigned by this routine, equal to the position of the
    model in the input list of fits.}
  \item{model.name}{Name of the fitted object.}
  \item{converged}{Logical values indicating whether this routine thinks the model
    converged or not.  Value is TRUE if the this routine thinks the model converged,
    FALSE otherwise. }
  \item{n.est.parameters}{Number of estimable parameters in the model.  This is MRA's guess at the number
    of estimable parameters in the model, not length of the coefficient vector.}
  \item{n.coefficients}{Number of coefficients in the model.  This is length of the coefficient
    vector without regard to estimability.  If \code{n.coefficients} > \code{n.est.parameters},
    the model is not full rank, and at least one coefficient is probably not estimable.}
  \item{aicc}{AIC of the model including the small sample correction =
    AIC + (2*\code{df}*(\code{df}+1)) / (\code{nan} - \code{df} - 1)}
  \item{delta.aicc}{Difference between AICc for the model and the minimum AICc in the table.}
  \item{aicc.wgt}{AICc model weights.  These weights equal exp(-.5*(delta.aicc)), scaled to sum to 1.0,}
  \item{qaicc}{QAIC of the model including the small sample correction =
    QAIC + (2*\code{df}*(\code{df}+1))/(\code{nan} - \code{df} - 1)}
  \item{delta.qaicc}{Difference between QAICc for the model and the minimum QAICc in the table.}
  \item{qaicc.wgt}{QAICc model weights.  These weights equal exp(-.5*(delta.qaicc)), scaled to sum to 1.0,}

}



\author{ Original routine by Eric Regehr, US Fish and Wildlife.  
Modified for MRA by Trent McDonald, WEST-INC, tmcdonald@west-inc.com }


\seealso{ \code{\link{F.cjs.estim}}, \code{\link{F.huggins.estim}}, \code{\link{F.fit.table}}, 
    \code{\link{plot.cjs}} }

\examples{

\donttest{
## Fit several CJS model to dipper data. Model average survival

## Time varying survival and capture (true CJS model)
data(dipper.histories)
ct <- as.factor( paste("T",1:ncol(dipper.histories), sep=""))
attr(ct,"nan")<-nrow(dipper.histories)
dipper.01 <- F.cjs.estim( ~tvar(ct,drop=c(1,2)), ~tvar(ct,drop=c(1,6,7)), dipper.histories )

## Linear trend in survival
cT <- 1:ncol(dipper.histories)
dipper.02 <- F.cjs.estim( ~tvar(ct,drop=c(1,2)), ~tvar(cT, nan=nrow(dipper.histories)), dipper.histories )

## No trend in survival
dipper.03 <- F.cjs.estim( ~tvar(ct,drop=c(1,2)), ~1, dipper.histories )

## Model average
mod.avg.surv <- F.cr.model.avg( ls(pat="^dipper"), what="s", fit.stat="aicc" )

## Plot
plot.cjs( mod.avg.surv )

}

}

\keyword{ survival }% at least one, from doc/KEYWORDS

\keyword{ models }% __ONLY ONE__ keyword per line
