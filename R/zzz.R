.First.lib<-function(libname, pkgname){

	library.dynam("mra", pkgname)

	cat(paste( "This is (M)ark (R)ecapture (A)nalysis.", "\n\nTrent McDonald, WEST Inc.\n(tmcdonald@west-inc.com, www.west-inc.com)\n") )
}
