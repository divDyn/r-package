#' FAD - LAD matrix from occurrence data
#' 
#' Function to generate range data from an occurrence dataset.
#' 
#' The function will output First and Last Appearance Data of the taxa in the dataset. Keep in mind that incomplete sampling will influence these data and will make the ranges appear shrunken.
#'
#' The following variables are produced:
#'
#' attribute \code{row.names}: The names of the taxa.
#'
#' \code{FAD}: First appearance data, given in time slice numbers.
#'
#' \code{LAD}: Last appearance data, given in time slice numbers.
#'
#' \code{range}: the range of the taxo in the given time slices
#'
#' \code{myFAD}: optionally, first appearance dates in numeric ages.
#'
#' \code{myLAD}: optionally, last appearance dates in numeric ages.
#'
#' \code{duration}: optionally, the duration of taxa in numeric ages. It is given for single-interval taxa as well assuming that the taxa's range span over the entire time slice.
#'
#' @param dat \code{(data.frame)} Occurrence data.
#' @param tax \code{(character)} The column name of taxon names.
#' @param bin \code{(character)} The column name of bin names, has to be an integer 
#' @param maxdate \code{(character)} Optional, the column names of maximum estimated ages of occurrences.
#' @param mindate \code{(character)} Optional, the column names of minimum estimated ages of occurrences.
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



#' Proportions of survivorship
#' 
#' This function will calculate both forward and backward survivorship proportions from a given occurrence dataset or FAD-LAD matrix.
#' 
#' Proportions of survivorship are great tools to visualize changes in the composition of a group over time (Raup, 1978). The curves show how a once coexisting set of taxa, called a cohort, looses its participants (forward survivorship) as time progress, or gains its elements as time is analyzed backwards.
#' Each value corresponds to a cohort in a bin (\emph{a}) and one other bin (\emph{b}). The value expresses what proportion of the analyzed cohort (present together in bin \emph{a}) is present in bin \emph{b}.
#' 
#' References:
#'
#' Raup, D. M. (1978). Cohort analysis of generic survivorship. Paleobiology, 4(1), 1-15.
#' 
#' @param dat \code{(data.frame)} The data frame containing fossil occurrences.
#' 
#' @param tax \code{(character)} The variable  name of the occurring taxa (variable type: \code{factor} or \code{character}).
#' 
#' @param bin \code{(character)} The variable name of the time slice numbers of the particular occurrences (variable type: \code{numeric}). Bin numbers should be in ascending order,can contain \code{NA}s, it can start from a number other than 1 and must not start with 0.
#' @param noNAStart \code{(logical)} Useful when the dataset does not start from bin \code{1}. Then \code{noNAStart=TRUE} will cut the first part of the resulting table, 
#' 						so the first row will contain the estimates for the lowest bin number.
#' 
#' @param fl \code{(matrix} or \code{data.frame}). If so desired, the function can be run on an FAD-LAD dataset, output by the \code{\link{fadLad}} function. 
#' 
#' @param method \code{(character)} Either \code{"forward"} or \code{"backward"}.
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



