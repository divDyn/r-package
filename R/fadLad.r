#' FAD - LAD matrix from occurrence data
#' 
#' Function to generate range data from an occurrence dataset.
#' 
#' The following variables are produced:
#'
#' row names: the name of the taxon.
#'
#' FAD: first appearance date in the given time slices.
#'
#' LAD: last appearance date in the given time slices.
#'
#' range: the range of the taxo in the given time slices
#'
#' FAD_DATE: optionally, first appearance dates in numeric ages.
#'
#' LAD_DATE: optionally, last appearance dates in numeric ages.
#'
#' DUR: optionally, the duration of taxa in numeric ages. It is given for single-interval taxa as well assuming that the taxa's range span over the entire time slice.
#'
#' nThrough: the number of through-ranging taxa
#'
#' @param dat (data.frame): Occurrence data.
#' @param tax (char): the column name of taxon names.
#' @param bin (char): the column name of bin names.
#' @param maxdate (char): optional, the column names of maximum radiometric dates of occurrences.
#' @param mindate (char): optional, the column names of minimum radiometric dates of occurrences.
#' @examples
#'	data(corals)
#'	fl <- fadLad(corals, tax="genus", bin="slc", maxdate="max_ma", mindate="min_ma")
#' @export
fadLad<-function(dat, tax, bin, maxdate=NULL, mindate=NULL){
	# for the prototype
#	dat <- corals
#	tax<- "genus"
#	bin <- "slc"
#	maxdate <- "max_ma"
#	mindate <- "min_ma"

	# defense
		if(!is.matrix(dat) & !is.data.frame(dat) ) stop("Invalid 'dat' argument.")
		if(!tax%in%colnames(dat))stop("Invalid 'tax' argument.")
		if(!bin%in%colnames(dat)) stop("Invalid 'bin' argument.")
		if(!is.numeric(dat[,bin])) stop("'bin' must be a numeric variable.")
	
	# in case there are numerical dates
	if(!is.null(maxdate) & !is.null(mindate)){
		# check whether the columns are there
		if(sum(colnames(dat)%in%c(maxdate, mindate))!=2) stop("Invalid 'mindate' or 'maxdate' arguments.")
		if(!is.numeric(dat[,maxdate]) | !is.numeric(dat[,mindate])) stop("Invalid 'mindate' and/or 'maxdate' arguments.")
		
		# single iteration
		tempFL<-tapply(INDEX=dat[,tax], X=1:nrow(dat), FUN=function(x){
			suppressWarnings(c(
				min(dat[x,bin], na.rm=T), 
				max(dat[x,bin], na.rm=T),
				max(dat[x,maxdate], na.rm=T), 
				min(dat[x,mindate], na.rm=T)))
		})
		fl<-matrix(unlist(tempFL), ncol=4, byrow=T)
		fl[!is.finite(fl)]<-NA
		fl<-as.data.frame(fl)
		colnames(fl) <- c("FAD", "LAD", "myFAD", "myLAD")
		
		# calculate durations
		fl$duration <- fl[,"myFAD"]-fl[,"myLAD"]
		fl$range <- fl[,"LAD"]-fl[,"FAD"]
		
		# order properly
		fl<-fl[,c("FAD", "LAD", "range", "myFAD", "myLAD", "duration")]
		
	# in case there aren't
	}else{
		
		# single iteration
		tempFL<-tapply(INDEX=dat[,tax], X=dat[,bin], FUN=function(x){
			c(min(x, na.rm=T), max(x, na.rm=T))
		})
		fl<-matrix(unlist(tempFL), ncol=2, byrow=T)
		fl<-data.frame(FAD=fl[,1], LAD=fl[,2])
		fl$range <- fl[,"LAD"]-fl[,"FAD"]
	}
	
	# add the species names
	rownames(fl)<- names(tempFL)
	return(fl)
	

}



#' Estimation of survivorship probabilities
#' 
#' This function will calculate both forward and backward survivorship probabilities from a given occurrence dataset or FAD-LAD matrix.
#' 
#' @param dat (data.frame): the data frame containing PBDB occurrences.
#' 
#' @param tax (char): variable  name of the occurring taxa (variable type: factor) - such as "occurrence.genus_name"
#' 
#' @param bin (char): variable name of the time slice numbers of the particular occurrences (variable type: int)- such as "slc" or whatever. Bin numbers should be in ascending order,can contain NA's, it can start from a number other than 1 and must not start with 0.
#' @param noNAStart (bool): useful when the dataset does not start from bin No. 1. Then noNAStart=TRUE will cut the first part of the resulting table, 
#' 						so the first row will contain the estimates for the lowest bin number.
#' 
#' @param fl (matrix or data.frame). If so desired, the function can be run on an FAD-LAD dataset. (fadLad)
#' 
#' @param method (character value): either "forward" or "backward".
#' 
#' @examples
#' data(corals)
#' surv<-survivors(corals, tax="genus", bin="slc", method="forward")
#' 
#' # plot
#' data(stages)
#' plotTS(stages, shading="series", boxes="per", xlim=c(260,0), 
#'   ylab="proportion of survivors present", ylim=c(0.01,1),plot.args=list(log="y"))
#'   
#' for(i in 1:ncol(surv)) lines(stages$mid, surv[,i])
#' 
#' @export
survivors<-function(dat, tax="genus", bin="slc", method="forward", noNAStart=FALSE,fl=NULL){
	
	# calculate FAD-LAD matrix first, if not provided
	if(is.null(fl)){
		fl <- fadLad(dat, tax, bin)
	}else{
		if(!"FAD"%in%colnames(fl) | !"LAD"%in%colnames(fl)) stop("Invalid FAD-LAD matrix.")
	}
	
	if(noNAStart){
		inter<-min(fl$FAD):max(fl$LAD)
	}else{
		inter <- 1:max(fl$LAD)
	}
	
	template<-rep(NA, max(inter))
	
	if(method=="forward"){
		
		sProbList<-sapply(inter, FUN=function(x){
			# who are there?
			bOrig <- fl$FAD<=x & x<=fl$LAD
			original<-fl[bOrig,]
			subVector <- x:max(inter)
			curSurv<-sapply(subVector, FUN=function(y){
				return(sum(original$FAD<=y & y<=original$LAD))
				
			})
			
			retVal <- template
			retVal[subVector]<-curSurv/sum(bOrig)
			return(retVal)
		})
	}
	
	if(method=="backward"){
		
		sProbList<-sapply(inter, FUN=function(x){
			# who are there?
			bOrig <- fl$FAD<=x & x<=fl$LAD
			original<-fl[bOrig,]
			subVector <- min(inter):x
			curSurv<-sapply(subVector, FUN=function(y){
				return(sum(original$FAD<=y & y<=original$LAD))
				
			})
			
			retVal <- template
			retVal[subVector]<-curSurv/sum(bOrig)
			return(retVal)
		})
	
	}
	
	return(sProbList)
}



