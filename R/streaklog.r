#' Utility functions for slicing gappy time series
#' 
#' The function returns where the continuous streaks start and how long they are, which can be used for efficient and flexible subsetting.
#'
#' The output list of \code{streaklog} contains the following elements:
#'
#' \code{starts}: the indices where the streaks start.
#'
#' \code{streaks}: the lengths of the individual streaks (number of values).
#'
#' \code{runs}: the number of streaks.
#'
#' The function whichmaxstreak() will return the indices of those values that are in the longest continuous streak.
#' 
#' @param x (\code{vector}) A vector with missing values.
#' 
#' @examples
#'	# generate a sequence of values
#'	  b<-40:1
#'	# add some gaps
#'	  b[c(1:4, 15, 19, 23:27)]  <- NA
#'	# the functions
#'	  streaklog(b)
#'	  whichmaxstreak(b)
#' @export
streaklog<-function(x)
{
	if(!is.logical(x)) x<-!is.na(x)
	
	starts<-numeric()
	streaks<-numeric()
	runs<-0
		
	for (i in 1:length(x))
	{
		
		#first (no i-1)
		if(i==1)
		{
			#if (x[i]==0)
			if (x[i]==1)
			{
				streaks<-c(streaks,1)
				runs<-1
				starts<-c(starts,i)
			}
		} 
		else
		{
			
			if (x[i]==1)
			{
				if (x[i-1]==1)
				{
					streaks[runs]<-streaks[runs]+1
					
				}
				if(x[i-1]==0)
				{
					streaks<-c(streaks,1)
					runs<-runs+1
					starts<-c(starts,i)
				
				}	
			}
			
			
		}
	
	}
		
	return(list(starts=starts, streaks=streaks,runs=runs))
}

#' @param which \code{(integer)} In case multiple streaks of the same length are found, which of them should be returned by the vector (integer).
#' @rdname streaklog
#' @export
whichmaxstreak<-function(x, which=-1)
{
	#x<-bZ
	lVec<-streaklog(x)

	start<-lVec$starts[max(lVec$streaks)==lVec$streaks]
	end<-max(lVec$streaks)+start-1
	
	
	if(which==-1 & length(start)>1)
	{
		warning(paste(length(start), "maximum length streaks were found.\nSet which= argument to select one, the first is returned now.\n"))
	}
	if (which==-1)
	{
		result<-start[1]:end[1]
	}else{
		if(which%in%1:length(start))
		{
			
			result<-start[which]:end[which]
			
		}else{
			result<-NULL
			cat(paste("Only ", length(start), " streaks were found, returning NULL.\n", sep=""))
		}
		
	}
	
	
	return(result)
}

vectorFromLog<-function(streak, len=T){
	#how long the final vector should be 
	if(len==T) len<-streak$starts[length(streak$starts)]+streak$streaks[length(streak$starts)]-1
	
	#empty values
		vec<-c()
		
	for(i in 1:length(streak$starts)){
		
		#add the false values before the start
		if(i==1)falses<-rep(F, streak$starts[i]-1)
		if(i>1) falses<-rep(F, streak$starts[i]-1-streak$starts[i-1])
		
		vec<-c(vec, falses)
		
		trues<-rep(T, streak$streaks[i])
		vec<-c(vec, trues)
	}
	#add the final falses (if any)
		falses<-rep(F,len-length(vec))
		vec<-c(vec, falses)
	
	return(vec)
}


#' Determination and omission of consecutive duplicates in a vector.
#'
#' \code{seqduplicated()} The function determines which elements of a vector are duplicates (similarly to \code{\link[base]{duplicated}}) in consecutive rows.
#' 
#' These functions are essentially about checking whether a value in a vector at index is the same as the value at the previous index. This seamingly primitive task had to be rewritten with Rcpp for speed and the appropriate handling of \code{NA} values. 
#' @param x (\code{vector}): input object.
#' @param na.rm (\code{logical}): Are \code{NA} entries to be treated as duplicates (\code{TRUE}) or just like a normal value (\code{FALSE})?
#' @param na.breaks (\code{logical}): If \code{na.rm=TRUE} and the \code{NA} values are surrounded by the same values, should the streak be treated as broken? Running \code{seqduplicated(, na.rm=TRUE)} on \code{(2, 1,NA, 1)} while setting \code{na.breaks} to \code{TRUE} will return \code{(FALSE, FALSE, TRUE, FALSE)}, and with \code{TRUE} it will return \code{(FALSE, FALSE, TRUE, TRUE)}. The results with the same argumentation of \code{collapse()} will be \code{(2,1)} and \code{(2,1,1)}.
#' @rdname collapse
#' @examples
#'   
#' # example vector
#'   examp <- c(4,3,3,3,2,2,1,NA,3,3,1,NA,NA,5, NA, 5)
#' 
#' # seqduplicated()
#'   seqduplicated(examp)
#' 
#'   # contrast with 
#'   duplicated(examp)
#' 
#'   # with NA removal
#'   seqduplicated(examp, na.rm=TRUE)
#' @export
seqduplicated <- function(x, na.rm=FALSE, na.breaks=TRUE){
	if(!is.vector(x)) stop("x has to be a vector.")
	# cast to numeric
	y <- as.integer(factor(x))
	# has to treat NAs the same

	if(na.rm & !na.breaks){
		boolNA <- !is.na(y)
		y<-y[boolNA]
	}
	logic<- .Call('_divDyn_seqduplicated', PACKAGE = 'divDyn', y)
	
	if(na.rm & !na.breaks){
		oldLogic <- logic
		logic <- rep(TRUE, length(x))
		logic[boolNA] <- oldLogic
	}
	# treat NAs as duplicates
	if(na.rm & na.breaks) logic[is.na(x)] <- TRUE

	return(logic)
}

#' Determination and omission of consecutive duplicates in a vector.
#'
#' \code{collapse()} Omits duplicates similarly to \code{\link[base]{unique}}, but only in consecutive rows, so the sequence of state changes remains, but without duplicates.
#'
#' @rdname collapse
#' @examples
#'  
#' # the same with collapse()
#'   collapse(examp)
#' 
#'   # contrast with 
#'   unique(examp)
#' 
#'   # with NA removal
#'   collapse(examp, na.rm=TRUE)
#'
#'   # with NA removal, no breaking
#'   collapse(examp, na.rm=TRUE, na.breaks=FALSE)
#' 
#' 
#' @export
collapse<- function(x, na.rm=FALSE, na.breaks=TRUE){
	if(!is.vector(x)) stop("x has to be a vector.")

	if(na.rm & !na.breaks) x<-x[!is.na(x)]
	
	# cast to numeric
	y <- as.integer(factor(x))
	# has to treat NAs the same

	logic<- .Call('_divDyn_seqduplicated', PACKAGE = 'divDyn', y)

	if(na.rm & na.breaks) logic[is.na(x)] <- TRUE

	x[!logic]

}


#' Filling of missing values in a vector, based on the marginal values of the gaps
#'
#' The function will loop through a vector and will substitute \code{NA} values with the value it last encountered or replaced.
#'
#' \code{NA}s won't be substituted when they are the first values the loop encounters. 
#'
#' @param x \code{(vector)}  Vector to be filled.
#' 
#' @param inc \code{(numeric)} Only if \code{x} is \code{numeric}, the function will increase the substituted value by this amount (useful for filling in sequences).
#' 
#' @param forward \code{(logical)} Should the loop go forward or backward?
#' @examples
#' # forward, replace with previous
#' dummy<- c(TRUE, FALSE, NA, TRUE, FALSE, NA)
#' fill(dummy)
#' 
#' # forward, replace with previous+1
#' dummy2 <- c(1,NA, 3, 1, 2, NA, NA, 9, NA,3)
#' fill(dummy2, inc=1)
#' 
#' # backward, replace with previous in loop direction
#' fill(dummy2, inc=0, forward=FALSE)
#' @export
fill <- function(x,forward=TRUE, inc=0){
	if(is.numeric(x)) x2 <- .Call('_divDyn_fillNumeric', PACKAGE = 'divDyn',x, dir=forward, inc)

	if(is.character(x)) x2 <- .Call('_divDyn_fillCharacter', PACKAGE = 'divDyn',x, dir=forward)

	if(is.logical(x)) x2 <- .Call('_divDyn_fillLogical', PACKAGE = 'divDyn',x, dir=forward)

	if(is.factor(x)){
		# convert to characte
		x2<-as.character(x)

		# the do recursion
		x2 <- fill(x2)

		# then covert back to factor
		x2<-factor(x2)
	}
	return(x2)
}


