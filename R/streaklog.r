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


