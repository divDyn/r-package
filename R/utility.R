# utility functions
# @param unused logical - should unused arguments be returned - \&code{"0"}
argdist <- function(x, xNames, argList, unused=FALSE){
	if(length(x)!=length(xNames)) stop("Number of functions does not match function names.")
	# get the arguments of every function
	argNames <- lapply(x, function(d){
		if(!is.null(d)) names(formals(d))
	})

	# look up which functions arguments are there
	distributed <- lapply(argNames, function(z){
	#	z<- argNames[[1]]
		argList[names(argList)%in%z]
	})
	names(distributed) <- xNames
	
	# check the unused arguments
	if(unused){
		distributed<-c(distributed, list(argList[!names(argList)%in%unlist(argNames)]))
		names(distributed)[length(distributed)] <- "0"
	}
	
	return(distributed)
}


# function to calculate the geometric mean of a vector
geom <- function(x, na.rm=TRUE){
	suppressWarnings(res<- exp(mean(log(x), na.rm=na.rm)))

	return(res)
}
