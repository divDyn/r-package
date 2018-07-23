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
#' \code{duration}: optionally, the duration of taxa in numeric ages. It is given for single-interval taxa as well assuming that the taxa's range span over the entire time slice.
#' 
#' @param dat \code{(data.frame)}: Occurrence data.
#' @param tax \code{(character)}: The column name of taxon names.
#' @param bin \code{(character)}: The column name of bin variable. If two column names are entered, then they will interpreted as the minimum and maximum age uncertainty. (see examples) 
#' @param ages \code{(logical)}: Are the bin entries ages (reversed time axis)? Setting ages to TRUE will replace the entries in 'bin' column(s) with their additive inverses.
#' @examples 
#' data(corals)
#' 
#' # binned data
#'   flBinned <- fadLad(corals, tax="genus", bin="slc")
#' 
#' # using basic bin lengths
#'   flDual <- fadLad(corals, tax="genus", bin=c("max_ma", "min_ma"), ages=TRUE)
#'
#' # single age esimate 
#'   data(stages)
#'   corals$mid <- stages$mid[corals$slc]
#'   flSingle <- fadLad(corals, tax="genus", bin="mid", ages=TRUE)
#' 
#' 
#' @export
fadLad<-function(dat, tax, bin, ages=FALSE){
	# for the prototype
#	dat <- corals
#	tax<- "genus"
#	bin <- "slc"

	# defense
		if(!is.matrix(dat) & !is.data.frame(dat) ) stop("Invalid 'dat' argument.")
		if(!tax%in%colnames(dat))stop("Invalid 'tax' argument.")
		if(sum(bin%in%colnames(dat))!=length(bin)) stop("Invalid 'bin' argument(s).")
		for(i in 1:length(bin)) if(!is.numeric(dat[,bin[i]])) stop("'bin' must be a numeric variable.")
	
	# demote the tax column to character
		taxVar <- as.character(dat[,tax])
	
	# coerce bin to matrix
		binVar <- as.matrix(dat[,bin])

	# are ages?
	if(ages){
		for(i in 1:length(bin)) binVar[,i] <- -binVar[,i]
	}
		
	# single iteration
	tempFL<-tapply(INDEX=taxVar, X=1:nrow(dat), FUN=function(x){
		# taxon specific part
		taxDat <- binVar[x,]
		
		suppressWarnings(c(
			min(taxDat, na.rm=T), 
			max(taxDat, na.rm=T)
			)
		)
	})
	
	# expand
	fl<-matrix(unlist(tempFL), ncol=2, byrow=T)
	fl[!is.finite(fl)]<-NA
	fl<-as.data.frame(fl)
	
	# column names
	colnames(fl) <- c("FAD", "LAD")
		
	# calculate durations
		fl$duration <- abs(fl[,"FAD"]-fl[,"LAD"])
		
	# add the species names
	rownames(fl)<- names(tempFL)
	return(fl)
	
}



#' Proportions of survivorship
#' 
#' This function will calculate both forward and backward survivorship proportions from a given occurrence dataset or FAD-LAD matrix.
#' 
#' Proportions of survivorship are great tools to visualize changes in the composition of a group over time (Raup, 1978). The curves show how a once coexisting set of taxa, called a cohort, loses its participants (forward survivorship) as time progress, or gains its elements as time is analyzed backwards.
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



