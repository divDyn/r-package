# Miscellaneous
.onUnload <- function (libpath) {
	library.dynam.unload("divDyn", libpath)
}

.onAttach<-function(libpath, pkgname){
	packageStartupMessage(
	"v0.7.1. See latest changes with 'news(package=\"divDyn\")'.\nReminder: Time flows from smaller to larger numbers in \'bin\'.\nYou can reverse this with 'ages=TRUE', where available.")
}
