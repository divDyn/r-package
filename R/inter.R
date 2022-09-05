#' Match the dates of a time-dependent variable with a predefined vector
#' 
#' The function takes a variable \code{x} (e.g. a vector or a list object), and reorders it to best match the dates provided in a vector \code{y}.
#' 
#' @param x Object to be reordered to match \code{y}.
#' 
#' @param y (\code{numeric}) The vector of dates (numeric values) to order to.
#' 
#' @param index (\code{logical}) If this argument is \code{TRUE}, only the indices will be returned that refer to the new order, rather than the reordered \code{x} variable.
#' 
#' @param ... Additional arguments passed to class-specific methods.
#' @rdname matchtime
#' @return An object of the class as \code{x} or a \code{numeric} vector.
#' @examples
#' # original vector
#' orig <- 1:10
#' # target values
#' targ <- c(5.1,4.2, 3.4, 2.7, 2.3)
#' # how do the two series match the best?
#' matchtime(orig, targ)
#' @export matchtime
setGeneric("matchtime", function(x,y,...) standardGeneric("matchtime"))

#' @rdname matchtime
setMethod(
	"matchtime", 
	signature="numeric",
	function(x, y, index=FALSE, ...){

		newIndex <- rep(NA, length(y))

		for(i in 1:length(y)){
			absDiff <- abs(y[i]-x)
			# which is closest
			newIndex[i] <- which(min(absDiff)==absDiff)[1]
		}

		if(!index){
			return(x[newIndex])
		}else{
			return(newIndex)
		}
	}
)

#' @rdname matchtime
setMethod(
	"matchtime", 
	signature="character",
	function(x, y, index=FALSE, ...){
		a<-as.numeric(x)
		newIndex<- matchtime(x=a,y=y,index=TRUE)

	if(!index){
		return(x[newIndex])
	}else{
		return(newIndex)
	}
})

#' @rdname matchtime
setMethod(
	"matchtime", 
	signature="list",
	function(x, y, index=FALSE, ...){
		newIndex <- matchtime(names(x), y=y, index=TRUE)
		if(!index){
			return(x[newIndex])
		}else{
			return(newIndex)
		}
	}
)


