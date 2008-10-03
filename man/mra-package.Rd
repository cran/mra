\name{mra-package}
\alias{mra-package}
\alias{mra}
\docType{package}
\title{
MRA - Mark Recapture Analysis
}
\description{
\bold{Description} - 
This package contains analysis functions, and associated routines, to conduct 
analyses of mark-recapture (capture-recapture) data using individual, 
time, and individual-time varying covariates. In general, these routines 
relate vectors of capture histories to vectors of covariates using 
a regression approach (Amstrup et al. 2005, Ch 9).  All capture, survival, 
transition, etc. parameters are functions of individual and time 
specific covariates, and the estimated parameters 
are coefficients in logistic-linear equations. 

\bold{Relationship to MARK} - 
For the most part, these routines perform a subset of the analyses available in 
program MARK or via the MARK front-end package, RMark.  
However, there are differences.  The most significant difference between this package
and MARK is parameterization.  The parameterization used here
does not utilize triangular 
"parameter information matrices" (PIMs) as MARK (and RMark) does.  
Because of this, the "design" matrix utilized by 
this package is not parallel to the "design" matrix of program MARK.  For those new 
to mark-recapture analysis, this parameterization difference will be inconsequential.  
The approach taken here provides equivalent modeling flexibility, yet is 
easier to grasp and visualize, in our opinion. 
For those already familiar with the PIMs used 
by program MARK, it is helpful to view the "PIMs" of this package as 
rectangular matrices of the real parameters. I.e., the "PIMs" of this package are 
rectangular matrices where cell (i,j) contains the real parameter (capture or survival) 
for individual i at capture occasion j. 

Analyses available here that are \emph{not} included in program MARK include:
\itemize{
\item Estimation of population size from open population CJS models via 
the Horvitz-Thompson estimator. 
\item Residuals, goodness of fit tests, and associated plots for assessing
model fit in open CJS models.
}

\bold{History} - 
These routines grew from consulting work on multiple mark-recapture projects. 
The original Fortran code, upon which the package is based, was written
by Dr. Bryan Manly for a northern spotted owl similation project in 1991.  
Dr. Manly is the one 
who originally envisioned and programed the non-PIM (or rectangular PIM, if you 
prefer) approach.  However, Dr. Manly did not realize what he had done. 
In 1997, Dr. Trent McDonald almost completely revised the original 
Fortran routines for use on a polar bear mark-recapture project.  At that time, 
the routines were stand-alone Fortran executables.  Subsequent revisions 
required by other projects included addition of closed-form variance estimates (originally, 
variances were estimated by the bootstrap), the Horvitz-Thompson size estimates, 
and goodness of fit testing.   In 2003, it finally dawned on Dr. McDonald how 
to call a Fortran DLL from S-Plus and R, thus 
eliminating the executable version and 
allowing S-Plus or R to do front-end data manipulation and plotting.  
S-Plus was abandoned in favor of R in 2004. 
After publication of Amstrup et al. (2005), 
Dr. McDonald realized that an official R package with documentation was 
needed, and learned how to make a package (not an easy process for him). 
Version 1.X of MRA contained routines for open population CJS models only.
Version 2.X of MRA added closed population Huggins model estimation routines. In future, 
addition of a major model type will increment the primary version number, revisions 
to routines for existing models will increment the secondary version number.  Thus, 
the number of major model types available in MRA will be the primary version number.
Throughout the process, several statisticians, including Dr. Manly, 
Dr. Jeff Laake and Dr. Gary White, have 
provided comments that helped shape the approach.

\bold{Ways You Can Help} - 
It is a no-brainer that others have R routines that perform various mark-recapture analyzes. 
The author of MRA views this effort as collaborative and 
welcomes routines, comments, and assistance in developing MRA. The overall goal 
is to develop MRA into a comprehensive repository for 
mark-recapture analyzes in R.  Indeed, collaboration is the strength of R and open source
code. 

Along these lines, the author of MRA would especially welcome routines that 
perform the following analyzes:
\itemize{
\item Continuous time models.  Especially those that allow inclusion of covariates.
\item Band recovery models.
\item Baysian models.
\item Joint live-dead recovery models.
\item MCMC methods or routines that can be applied to exiting models. 
\item Plotting methods for exiting models.
\item Model selection methods for existing models. 
\item Simulation methods and routines.
}

The above is only a partial list.  All routines are welcome and will be considered.  
Original authors 
will of course be acknowledged in the routine's documentation.  Assistance in formating 
documentation files is available.  I.e., either the author of MRA will write the 
documentation or the routine's original author will be sent a text file to 
edit that contains all the necessary mark-up 
and sectioning. 

If you have a routine that would be useful, email MRA's maintainer 
at the address below.  Supporting papers and data sets can be distributed 
with MRA, and are encouraged. 

}

\details{
\tabular{ll}{
Package: \tab mra\cr
Type: \tab Package\cr
License: \tab GNU General Public License\cr
}

List of routines:
\preformatted{
F.cjs.covars            Returns matricies that can be used to fit a CJS model
F.cjs.estim             Estimation of Cormack-Jolly-Seber (CJS) open population model
F.cjs.gof               Goodness-of-fit tests for CJS models
F.cr.model.matrix       Function that returns two 3-D design matricies.
F.3d.model.matrix       Function that expands a formula into the type of 3-D array
                        needed for estimation.                        
F.huggins.estim         Estimation of Huggin's closed population model
F.sat.lik               Returns the saturated likelihood for a CJS model that does 
                        not contain individual covariates
dipper.data             European Dipper data set
lines.cjs               Lines method for cjs objects
plot.cjs                Plot method for cjs objects
predict.cjs             Predicted values for active cells of a CJS model. 
print.cjs               Print method for cjs objects
print.hug               Print method for hug (Huggin's model) objects
print.nhat              Pring method for size estimates from a CJS model
residuals.cjs           Residuals for CJS models
ivar                    Function for specifying individual varying effects
tvar                    Function for specifying time varying effects
}
}

\author{
Trent McDonald

Maintainer: Trent McDonald <tmcdonald@west-inc.com>
}

\references{
Amstrup, S.C., T.L. McDonald, and B.F.J. Manly. 2005. \emph{Handbook of 
Capture-Recapture Analysis}, Princeton: Princeton University Press.
}

\keyword{ package }
