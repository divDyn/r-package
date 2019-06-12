
#' Discretization of continuous time dimension - slicing
#' 
#' The function will slices time with a given set of boundaries and produce a time scale object if desired.
#' 
#' Due to stratigraphic constraints, we can only process deep time data, when it is sliced to discrete bins. It is suggested that you do this separately for most of your analyses. This function is also used by the code{\link{divDyn}} function when \code{age} entries are provided.
#' 
#' @param x (\code{numeric}) Vector of continouos age/time estimates.
#' @param breaks (\code{numeric}) Vector of boundaries, the breaks argument of the \code{\link[base]{cut}} function
#' @param offset (\code{numeric}) Single  value. If desired the resulting integer bin numbers can be offset by some amount.
#' @param ts (\code{logical}) Should a time scale object be also produced when the function is run? 
#' @param revtime (\code{logical}) Should the time dimension be reversed? This argument is set to \code{TRUE} by default, meaning that the function will reverse the order of time: smaller values of \code{x} will be translated to higher values (\code{slc}) in the function output.
#' @examples
#' y<- runif(200, 0,100)
#' au <- slice(y, breaks=seq(0, 100, 10))
#' withTs <- slice(y, breaks=seq(0, 100, 10), ts=TRUE)
#' @export
slice<-function(x, breaks, offset=0, ts=TRUE, revtime=TRUE){
	# only allow to pass through if it is numeric
	if(!is.numeric(x)) stop("The provided time vector is not numeric.")
	if(revtime){
		x<- -x
		breaks <- -breaks
	}

	if(!is.numeric(offset) | length(offset)>1) stop("Invalid 'offset' parameter.")
	# first, sort the breaks
	breaks<-sort(breaks)
	# cut the vector
	y <- cut(x, breaks=breaks, labels=FALSE)
	
	timeMat <- cbind(
			breaks[1:(length(breaks)-1)],
			breaks[2:length(breaks)])
	lev<-apply(timeMat,1, mean)

	y<- y+offset
	if(ts==FALSE){
		if(revtime) lev <- -lev
		return(list(slc=y, lev=lev))
	}else{
		timedim<-1:nrow(timeMat)
		timedim<- timedim+offset

		# depending on whether you use ages or other
		if(revtime) {
			tMat <- cbind(-timeMat, -lev,timedim)
		}else{
			tMat <- cbind(timeMat, lev,timedim)
		}

		tMat <- as.data.frame(tMat)
		colnames(tMat) <- c("bottom", "top", "mid", "slc")
		return(list(slc=y, ts=tMat[,c("bottom","mid","top", "slc")]))
	}
}

