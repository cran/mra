\name{F.cjs.estim}
\alias{F.cjs.estim}
\alias{cr.estim}
\alias{CJS}
\alias{cjs}
\alias{Cormack}
\alias{Jolly}
\alias{Seber}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ F.cjs.estim - Cormack-Jolly-Seber estimation }
\description{
   Estimates Cormack-Jolly-Seber (CJS) capture-recapture models with individual, time, and
	individual-time varying covariates using the "regression" parameterization of 
	Amstrup et al (2006, Ch 9). For live recaptures only. Losses on capture allowed.
	Uses a logistic link function to relate probability of capture and survival 
	to external covariates.  
}
\usage{
F.cjs.estim(capture, survival, histories, cap.init, sur.init, group,  
	algorithm = 1, cov.meth = 1, nhat.v.meth = 1, c.hat = -1, df = NA,
	intervals=rep(1,ncol(histories)-1))
}

%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{capture}{ Formula specifying the capture probability model.  Must be a formula object with 
	no response, then ~, followed by the names of 2-D arrays of covariates to fit in the capture model. 
	For example: 'capture = ~ age + sex', where age and sex are matrices of size NAN X NS 
	containing the age and sex covariate values.  NAN = number of animals = number of rows in 
	\code{histories} matrix (see below).  NS = number of samples = number of columns in histories
	matrix (see below). Number of matrices specified in the capture model is assumed to be NX. Specify 
	intercept only model as 'capture = ~ 1'.  Specify model without an intercept using 'capture = ~ -1 + x'.  }
  \item{survival}{ Formula specifying the survival probability model.  Must be a formula object with 
	no response, then ~, followed by the names of 2-D arrays of covariates to fit in the survival model. 
	For example: survival = ~ year + ageclass where year and ageclass are matrices of size NAN X NS 
	containing year and ageclass covariate values.  NAN = number of animals = number of rows in 
	\code{histories} matrix (see below).  NS = number of samples = number of columns in histories
	matrix (see below). Number of matrices specified in the survival model is assumed to be NY. Specify 
	intercept only model as 'survival = ~ 1'.  Specify model without an intercept using 'survival = -1 + y'.}
  \item{histories}{ A NAN X NS = (number of animals) X (number of capture occasions) matrix 
	containing capture histories.  Capture histories are comprised of 0's, 1',s and 2's. 
	0 in cell (i,j) means animal i was not captured on occasion j, 1 in cell (i,j) means animal 
	i was captured on occasion j and released live back into the population, 2 in cell (i,j) means 
	animal i was captured on occasion j and was not released back into the population (e.g., it died).
	Animals with '2' as the last non-zero in thier histories are 'censored'. I.e., Animals with 
	2's are 'removed' from the likelihood after the occasion with the 2. }
  \item{cap.init}{ (optional) Vector of initial values for coefficients in the capture model.  One element 
	per covariate in \code{capture}. The default value usually works. }
  \item{sur.init}{ (optional) Vector or initial values for coefficients in the survival model.  One element 
	per covariate in \code{survival}. The default value usually works. }
  \item{group}{ (optional) A vector of length NAN giving the (non-changing) group membership of every 
	captured animal (e.g., sex). Group is used only for computing TEST 2 and TEST 3. 
	TEST 2 and TEST 3 are computed separately for each group. E.g., if group=sex, 
	TEST 2 and TEST 3 are computed for each sex. TEST 2 and TEST3 are used only
	to estimate C-hat. Results of TEST 2 and TEST 3 are output only 
	to the log file ('mrawin.log'). See \code{c.hat} for pooling rules for these 
	test components to estimate C-hat.}
  \item{algorithm}{ Integer specifying the maximization algorithm to use.  If \code{algorithm} = 1, 
	the VA09AD algorithm from the HSL library is used. The VA09AD algorithm is very reliable, 
	has been tested extensively (same algorithm as Program MARK), and will almost always find the 
	maximum likelihood estimates. This parameter was added to allow easy addition of other 
	algorithms in the future, but no other options are currently implemented because VA09AD 
	works so well. Anything other than 1 is ignored. }
  \item{cov.meth}{ Integer specifying the covariance estimation method. \code{cov.meth} = 1 
	takes numerical 2nd partial derivatives.  \code{cov.meth} = 2 inverts the Hessian 
	of the maximization.  Method 1 (numeric 2nd derivatives) is the preferred method, but
        is very computationally expensive and therefore much slower than method 2.  Difference in 
        speed of methods should be minimal when number of parameters < 15. Method 2 variances are 
     	generally very close to method 1 variances, and could be used during initial model fitting 
	stages when exact estimation of variance may not be necessary.    }
  \item{nhat.v.meth}{ Integer specifying method for computing variance estimates 
	of population size estimates. \code{nhat.v.meth} = 1 uses the variance estimator of 
	Taylor et al. 2002, Ursus, p. 188 which is the so-called Huggins variance estimator, and incorporates 
	covariances.  \code{nhat.v.meth} = 2 uses the variance estimator of Amstrup et al. 
	2005 (p. 244, Eqn. 9.10), which is the same variance estimator as \code{nhat.v.meth} = 1 
	with more 2nd order approximation terms included.  Method 2 should provide better variances 
	than method 1, especially if the coefficient of variation of capture probabilities 
	are >1.0, but method 2 has not been studied as much as method 1. 
	\code{nhat.v.meth} = 3 uses the variance estimator of McDonald and Amstrup, 1999, JABES, which is a 1st order approximation
	that does not incorporate covarainces. Method 3 is much faster than methods 1 and 2 
	and could be easily calculated by hand, but 
	should only be used when there is little capture heterogeneity. }
  \item{c.hat}{ External (override) estimate of variance inflation factor (\code{c.hat}) to use
	during estimation.  If input value of \code{c.hat} is <= 0, MRAWIN computes 
	an estimate of variance inflation based on TEST 2 and TEST 3
	applied to groups (if called for, see \code{group} above) 
	using Manly, McDonald, and McDonald, 1993, rules for pooling. I.e., 
	all cells in each TEST 2 or TEST 3 Chi-square component table must be 
	>= 5 before that component contributes to the estimate of C-hat.   This rules is  
	slightly different than program MARK's pooling rules, so MRA's and MARK's
	estimates of \code{c.hat} will generally be different.  If the input \code{c.hat} > 0, 
	MRAWIN does not estimate C.hat, and uses the supplied value. }
  \item{df}{ External (override) model degrees of freedom to use during estimation.  
	If \code{df} == NA, the number of parameters is estimated from the rank of the 
	matrix of 2nd derivatives or Hessian, depending on \code{cov.meth} parameter.  
	If \code{df} <= 0, the number of parameters will be 
	set to NX+NY = the number of estimated coefficients. Otherwise, if \code{df} > 0, 
	the supplied value is used.  Only AIC, QAIC, AICc, and QAICc are dependent on 
	this value (in their penalty terms). }
  \item{intervals}{Time intervals.  This is a vector of length \code{ncol(histories)-1} (i.e., 
	number of capture occasions minus 1) specifying relative time intervals between occasions. 
	For example, if capture occasions occurred in 1999, 2000, 2005, and 2007 \code{intervals} would 
	be set to \code{c(1,5,2)}.  Estimates of survival are adjusted for time intervals between 
	occasions assuming an exponential lifetime model, i.e., probability of surviving 
	from occasion \code{j} to occasion \code{j+1} is \code{Phi(j)^(jth interval length)}, and 
	it is the \code{Phi(j)}'s that are related to covariates through the survival model.  In 
	other words, all survival estimates are for an interval of length 1.  If an interval 
	of 1 is one year, then all survival estimates will be annual survival, with probability 
	of surviving 2 years equal to annual survival squared, probability of surviving 3 years 
	equal to annual survival cubed, etc.    
	 }.
}

\details{
	This is the work-horse routine for estimating CJS models. It compiles all the 
	covariate matrices, then calls a Fortran routine to maximize the CJS likelihood 
	and perform goodness-of-fit tests.  Horvitz-Thompson-type population size estimates 
	are also computed by default. 
  
	A log file, named \code{mra.log}, is written to the current directory.  This file contains
     	additional details, such as individual Test 2 and Test 3 components, in a semi-friendly 
     	format.  This file is overwritten each run. 
}

\value{An object (list) of class c("cjs","cr") with many components.  Use \code{print.cr} to print 
	it nicely. Use \code{names(fit)}, where the call was \code{fit <- F.cr.estim(...)}, 
	to see names of all returned components. To see values of individual components, 
	issue commands like fit\$s.hat, fit\$se.s.hat, fit\$n.hat, etc. 

 	Components of the returned object are as follows:
  \item{aux}{Auxiliary information, mostly stored input values. This is a list containing:
	\$call, \$nan=number of animals, \$ns=number of samples, \$nx=number of coefficients 
	in capture model, \$ny=number of coefficients in survival model, \$cov.name=names 
	of coefficients, \$ic.name=name of capture history matrix, \$mrawin.version=
	version number for this package, \$package=software (S-Plus or R) used to 
	estimate the model, \$package.version=S-Plus or R version used for estimation,
	\$run.date=date the model was estimated.}
  \item{loglik}{Maximized CJS likelihood value for the model}
  \item{deviance}{Model deviance = -2*\code{loglik}.  This is relative deviance, see help for 
	\code{F.sat.lik}.}
  \item{aic}{AIC for the model = \code{deviance} + 2*(df).  df is either the estimated number of independent 
	parameters (by default), or NX+NY, or a specified value, depending on the input value of DF parameter.}
  \item{qaid}{QAIC (quasi-AIC) = (\code{deviance} / \code{vif}) + 2(df)}
  \item{aicc}{AIC with small sample correction = AIC + (2*\code{df}*(\code{df}+1)) / (\code{nan} - \code{df} - 1)}
  \item{qaicc}{QAIC with small sample correction = QAIC + (2*\code{df}*(\code{df}+1))/(\code{nan} - \code{df} - 1)}
  \item{ebc}{Empirical Bayes Criterion of Peterson (1986), Stat and Prob letters, p227.}
  \item{vif}{Variance inflation factor used = estimate of c.hat = \code{chisq.vif} / \code{chisq.df} }
  \item{chisq.vif}{Composite Chi-square statistic from Test 2 and Test 3 used to compute \code{vif}, based 
	on pooling rules.}
  \item{vif.df}{Degrees of freedom for composite chi-square statistic from Test 2 and Test 3, based on pooling rules.}
  \item{parameters}{Vector of all coefficient estimates, NX capture probability coefficients first, 
	then NY survival coefficients.  This vector is length NX+NY regardless of estimated DF.}   
  \item{se.param}{Standard error estimates for all coefficients.  Length NX+NY.}
  \item{se1.param}{Standard error estimates from 'alternate' method.}
  \item{capcoef}{Vector of coefficients in the capture model.  Length NX.}
  \item{se.capcoef}{Vector of standard errors for coefficients in capture model. Length NX.}
  \item{surcoef}{Vector of coefficients in the survival model. Length NY.}
  \item{se.surcoef}{Vector of standard errors for coefficients in survival model. Length NY.}
  \item{covariance}{Variance-covariance matrix for the estimated model coefficients.  
	Size (NX+NY) X (NX+NY).}
  \item{p.hat}{Matrix of estimated capture probabilities computed from the model.  
	One for each animal each occasion.  Cell (i,j) is estimated capture probability 
	for animal i during capture occasion j.
	Size NAN X NS.  First column corresponding to first capture probability is NA because
	cannot estimate P1 in a CJS model.}
  \item{se.p.hat}{Matrix of standard errors for estimated capture probabilities. One for each animal 
	each occasion.	Size NAN X NS.  First column is NA.}
  \item{s.hat}{Matrix of estimated survival probabilities computed from the model.  One for each animal each occasion.
	Size NAN X NS.  Cell (i,j) is estimated probability animal i survives from 
	occasion j to j+1.  There are only NS-1 intervals between occasions. Last 
	column corresponding to survival between occasion NS and NS+1 is NA.}
  \item{se.s.hat}{Matrix of standard errors for estimated survival probabilities. 
	Size NAN X NS.  Last column is NA.}
  \item{df}{Number of estimable parameters in the model. df is either the estimated number of independent 
	parameters (by default) based on rank of the variance matrix, 
	or NX+NY, or a specified value, depending on the input value of DF parameter.} 
  \item{message}{A string indicating whether the maximization routine converged.}
  \item{exit.code}{Exit code from the maximization routine.  Interpretation for \code{exit.code} is in \code{message}.}
  \item{cov.code}{A code indicating the method used to compute the covariance matrix.}
  \item{cov.meth}{String indicating method used to compute covariance matrix.  Interprets \code{cov.code}.}
  \item{n.hat}{Vector of Horvitz-Thompson estimates of population size.  Length NS.  No estimate for 
	first occasion.}
  \item{se.n.hat}{Estimated standard errors for \code{n.hat} estimates.  Computed using method 
	specified in \code{nhat.v.meth}.}
  \item{n.hat.lower}{Lower limit of \code{n.hat.conf} percent on \code{n.hat}.  Length NS.} 
  \item{n.hat.upper}{Upper limit of \code{n.hat.conf} percent on \code{n.hat}.  Length NS.} 
  \item{n.hat.conf}{Confidence level of intervals on \code{n.hat}}
  \item{nhat.v.meth}{Code for method used to compute variance of \code{n.hat}}
  \item{num.caught}{Vector of observed number of animals captured each occasion.  Length NS.}
  \item{fitted}{Matrix of fitted values for the capture histories.  Size NAN X NS.  Cell (i,j)
	is expected value of capture indicator in cell (i,j) of \code{histories} matrix.}
  \item{residuals}{Matrix of Pearson residuals = ($h_{ij} - \psi_{ij})^2 / \psi_{ij}$, where 
	$\psi_{ij}$ is fitted value for cell (i,j) and $h_{ij}$ is the capture indicator for 
	animal i at occasion j.  Size NAN X NS. See help for "overall test" in 
	documentation for \code{F.cjs.gof} for a description of $\psi_{ij}$. }
  \item{resid.type}{String describing the type of residuals computed.  Currently, only Pearson
	residuals are returned.}
}



\references{ 
Taylor, M. K., J. Laake, H. D. Cluff, M. Ramsay, and F. Messier. 2002. Managing the risk from hunting for 
	the Viscount Melville Sound polar bear population. Ursus 13:185-202.

Manly, B. F. J., L. L. McDonald, and T. L. McDonald. 1999. The robustness of mark-recapture methods: a 
	case study for the northern spotted owl. Journal of Agricultural, Biological, and Environmental 
	Statistics 4:78-101.

Huggins, R. M. 1989. On the statistical analysis of capture experiments. Biometrika 76:133-140.

Amstrup, S. C., T. L. McDonald, and B. F. J. Manly (editors). 2005. Handbook of 
	Capture-Recapture Analysis. Princeton University Press.

Peterson. 1986. Statistics and Probability Letters. p.227.

McDonald, T. L., and S. C. Amstrup. 2001. Estimation of population size using 
	open capture-recapture models. Journal of Agricultural, Biological, and Environmental 
	Statistics 6:206-220.
 }
\author{ Trent McDonald, WEST-INC, tmcdonald@west-inc.com }

\seealso{ \code{\link{print.cjs}}, \code{\link{residuals.cjs}}, \code{\link{plot.cjs}}, 
  \code{\link{F.cjs.covars}}, \code{\link{F.cjs.gof}} }

\examples{

# Fit CJS model to dipper data, time-varying capture and survivals.
data(dipper.histories)
xy <- F.cjs.covars( nrow(dipper.histories), ncol(dipper.histories) )
for(j in 1:ncol(dipper.histories)){ assign(paste("x",j,sep=""), xy$x[,,j]) } # Extract 2-D matricies of 0s and 1s
dipper.cjs <- F.cjs.estim( ~x2+x3+x4+x5+x6, ~x1+x2+x3+x4+x5, dipper.histories )
print(dipper.cjs)
}

\keyword{ survival }% at least one, from doc/KEYWORDS

\keyword{ models }% __ONLY ONE__ keyword per line
