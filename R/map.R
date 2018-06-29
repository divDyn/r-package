#' Mapping multiple entries to categories
#' 
#' This basic function replaces groups of values in a vector with single values with the help of a key object.
#' 
#' Online datasets usually contain overly detailed information, as enterers intend to conserve as much data in the entry process, as possible. However, in analyses some values are treated to represent the same, less-detailed information, which is then used in further procedures. The \code{map} function allows users to do this type of multiple replacement using a specific object called a \code{'key'}. 
#'
#' A \code{key} is an informal class and is essentially a \code{list} of \code{vectors}. In the case of \code{character} vectors as \code{x}, each vector element in the \code{list} corresponds to a set of entries in \code{x}. These will be replaced by the name of the \code{vector} in the \code{list}, to indicate their assumed identity. 
#'
#' In the case of \code{numeric} \code{x} vectors, if the \code{list} elements of the \code{key}
#' are \code{numeric} vectors with 2 values, then this vector will be treated as an interval. The same value will be assigned to the entries that are in this interval (Example 2). If \code{x} contains values that form the boundary of an interval, than either only the one of the two boundary values can be considered to be in the interval (see the \code{incbound} argument to set which of the two).
#' The elements of \code{key} are looped through in sequence. If values of \code{x} occur in multiple elements of \code{key}, than the last one will be used (Example 3).
#' 
#' Examples of this data type have been included (\code{\link{keys}}) to help process Paleobiology Database occurrences.
#' @param x \code{(vector)} Object containing the values to be replaced. 
#' 
#' @param key \code{(list)} A list of vectors. Each \code{vector} includes the possible elements that will be replaced in a group, the \code{names} of the \code{vector}s will be the replacement values. Also has to include an element named 'default' with a single value. (see examples)
#' 
#' @param incbound \code{(character)} Either \code{"lower"} or \code{"higher"}. Interval identifiers will be treated with different interval rules. \code{"lower"} will treat the lowest entry as included, \code{"higher"} works the opposite. The argument will be renamed to 'include.lowest' to make the interface easier to remember.
#' @examples
#' # Example 1
#' # x, as character
#'    set.seed(1000)
#'    toReplace <- sample(letters[1:6], 15, replace=T)
#' # a and b should mean 'first', c and d 'second' others: NA
#'    key<-list(first=c("a", "b"), second=c("c", "d"), default=NA)
#' # do the replacement
#'   map(toReplace, key)
#' 
#' # Example 2 - numeric entries and mixed types
#' # basic vector to be grouped
#'   toReplace2<-1:16
#' 
#' # replacement rules: 5,6,7,8,9 should be "more", 11 should be "eleven" the rest: "other"
#'   key2<-list(default="other", more=c(5,10),eleven=11)
#'   map(toReplace2, key2)
#' 
#' # Example 3 - multiple occurrences of same values
#' # a and b should mean first, a and should mean 'second' others: NA
#'   key3<-list(first=c("a", "b"), second=c("a", "d"), default=NA)
#' # do the replacement (all "a" entries will be replaced with "second")
#'   map(toReplace, key3)
#'    
#' @export
map<- function(x, key, incbound="lower"){
	# assign the defaults
	default<-key$default
	if(is.null(default)) default <- NA
	
	y<- rep(default, length(x))
	
	# delete the default from the set
	key$default <- NULL
	
	for(i in 1:length(key)){
		tempVect <-key[[i]]
		
		# NA is present here
		if(sum(is.na(tempVect)!=0)){
			y[is.na(x)] <- names(key)[i]
			tempVect <- tempVect[!is.na(tempVect)]
		}
		
		if(is.numeric(tempVect)){
			if(length(tempVect)==2){
				if(incbound=="lower"){
					y[min(tempVect)<= x & x<max(tempVect)] <- names(key)[i]
				}
				if(incbound=="higher"){
					y[min(tempVect)< x & x<=max(tempVect)] <- names(key)[i]
				}
			
			}else{
				tempVect <- as.character(tempVect)
			}
		
		}
		if(is.character(tempVect)){
			# look up
			y[x%in%tempVect] <- names(key)[i]
		}
		
	}
	
	return(y)
	
}
