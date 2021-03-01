# Miscellaneous
.onUnload <- function (libpath) {
	library.dynam.unload("divDyn", libpath)
}

.onAttach<-function(libpath, pkgname){
	packageStartupMessage(
	"v0.8.1. See latest changes with \'news(package=\"divDyn\")\'. Default: \nTime flows forward with \'bin\', and backward with \'age\' and \'slice()\'.")
}
