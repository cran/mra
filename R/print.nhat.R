"print.nhat" <- 
function( x, ... ){

if( !is.null(names(x$n.hat)) ){
	nms <- names(x$n.hat) 
} else {
	nms <- 1:length(x$n.hat)
}

if( x$nhat.v.meth == 1 ){
	mess <- "(Huggins estimator)"
} else if (x$nhat.v.meth == 2){
	mess <- "(Huggins estimator + higher terms)"
} else if (x$nhat.v.meth == 3){
	mess <- "(McDonald and Amststrup estimator)"
} else {
	mess <- "(Unknown)"
}
	
ind <- !is.na( x$n.hat )

cat(paste( x$n.hat.conf*100, "% Confidence Intervals\n", sep="" ))
cat(paste( "Standard errors computed by method=", x$nhat.v.meth, " ", mess, "\n", sep=""))

se.n.hat <- x$se.n.hat
z.lower.ci <- x$n.hat.lower
z.upper.ci <- x$n.hat.upper
n.hat <- x$n.hat
se.n.hat <- x$se.n.hat

z.lower.ci[ind] <- round(x$n.hat.lower[ind])
z.upper.ci[ind] <- round(x$n.hat.upper[ind])
n.hat[ind] <- round(n.hat[ind], 1)
se.n.hat[ind] <- round(se.n.hat[ind], 1)


if( is.null( se.n.hat ) ) se.n.hat<-rep(NA, length(x$n.hat)) 
if( is.null( z.lower.ci ) ) z.lower.ci<-rep(NA, length(x$n.hat))
if( is.null( z.upper.ci ) ) z.upper.ci<-rep(NA, length(x$n.hat))

out <- paste( format( c(" ", "Occasion", nms)), 
	format( c(" ", "N est", n.hat) ),  
	format( c(" ", "Se(N)", se.n.hat ) ),
	format( c("Norm", "Low", z.lower.ci)), 
	format( c("Norm", "High", z.upper.ci)), 
	sep="\t")

if( exists( "supsmu" ) ){
	smu.n <- supsmu( seq( along=nms[ind] ), x$n.hat[ind], bass=.5 )
	out <- paste( out, 
		format( c("Smooth","N est", c(NA, round(smu.n$y,1)) )),
		"\n", sep="\t")
} else {
	out <- paste( out, "\n")
}


cat( out )

invisible()

}

