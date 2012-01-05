.onAttach<-function(libname, pkgname){

	#library.dynam("mra", pkgname, libname)  # deprecated
    #v <- packageDescription("mra", fields="Version")  # this requires utils package, and I don't want to make mra dependent on utils, just for this.

    packageStartupMessage( "Mark-Recapture Analysis (vers 2.10)" )  # You have to change this every version
	packageStartupMessage("\nTrent McDonald, WEST Inc.\n(tmcdonald@west-inc.com, www.west-inc.com)") 
}


