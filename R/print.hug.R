print.hug <- function( x, ... ){

nx <-   x$aux$nx 
ny <-   x$aux$ny

cap.coef <- round( x$capcoef, 5 )
recap.coef <- round( x$recapcoef, 5 )
se.cap <- round( x$se.capcoef, 5 )
se.recap <- round( x$se.recapcoef, 5 )

if(nx > ny){
    recap.coef<- c(recap.coef, rep("", nx-ny))
    se.recap  <- c(se.recap, rep("", nx-ny))
} else {
    cap.coef<- c(cap.coef, rep("", ny-nx))
    se.cap  <- c(se.cap, rep("", ny-nx))    
}

cat("Call:\n")
print(x$aux$call)
cat("\n")

cat(paste( format( c(" Capture var", names(cap.coef))), 
    format( c(" Est", cap.coef) ),  
    format( c(" SE", se.cap) ),
    "  ",
    format( c(" Recapture var", names(recap.coef))), 
    format( c(" Est", recap.coef) ), 
    format( c(" SE", se.recap) ),
    "\n", sep="  "))

cat("\nPopulation Size Estimate (se): ")
cat(paste( round(x$n.hat, 4), " (", round(x$se.n.hat,4), ")\n", sep=""))
cat(paste( x$n.hat.conf*100, "% confidence interval for population size: ", 
    round(x$n.hat.lower,2), " to ", round(x$n.hat.upper,2), "\n", sep=""))
cat(paste("Individuals observed: ", round(x$num.caught), "\n", sep=""))
cat(paste("Effective sample size: ", round(x$n.effective), "\n", sep=""))


cat(paste("\nMessage =", x$message ))
cat(paste("\nNumber of estimable coefficients (estimated) = ", x$df))
cat(paste("\nLog likelihood = ", x$loglike))
cat(paste("\nDeviance = ", x$dev))
cat(paste("\nAIC = ", x$aic))
cat(paste("\nAICc = ", x$aicc))
cat("\n")

invisible()

}
