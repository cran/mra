.First.lib<-function(libname, pkgname){

	library.dynam("mra", pkgname)

	cat(paste( "Mark-Recapture Analysis.", "\n\nTrent McDonald, WEST Inc.\n(tmcdonald@west-inc.com, www.west-inc.com)\n") )
}
