#' @importFrom Rcpp evalCpp
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics abline plot polygon rect text clip par points segments
#' @importFrom stats binom.test dbinom loess median pbinom predict quantile
#' @importFrom utils flush.console
#' @importFrom methods new
#' @useDynLib divDyn, .registration=TRUE

# Miscellaneous
.onUnload <- function (libpath) {
	library.dynam.unload("divDyn", libpath)
}

.onAttach<-function(libpath, pkgname){
	packageStartupMessage(
	"v0.8.2. See latest changes with \'news(package=\"divDyn\")\'. Default: \nTime flows forward with \'bin\', and backward with \'age\' and \'slice()\'.")
}
