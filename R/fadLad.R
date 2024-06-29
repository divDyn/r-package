#' FAD - LAD matrix from occurrence data
#' 
#' Function to generate range data from an occurrence dataset.
#' 
#' The function will output First and Last Appearance Dates of the taxa in the dataset. Keep in mind that incomplete sampling will influence these data and will make the ranges appear shrunken.
#'
#' The following variables are produced:
#'
#' \code{row.names} attribute: The names of the taxa.
#'
#' \code{FAD}: First appearance dates in time bin nmbers or ages.
#'
#' \code{LAD}: Last appearance dates in time bin numbers or ages.
#'
#' \code{duration}: The durations of taxa in bin numbers or ages.
#' 
#' @param x \code{(data.frame)}: Occurrence data.
#' @param tax \code{(character)}: Column name of taxon names.
#' @param bin \code{(character)}: Column name(s) of the discreet bin variable(s). If two column names are entered, then they will interpreted as minimum and maximum uncertainty (see examples) By default, time flows from lower to higher numbers. You can change this behavior by setting \code{revtime=FALSE}. Either a \code{bin} or \code{age} argument is mandatory. 
#' @param age \code{(character)}: Column name(s) of the continuous age variable(s). If two column names are entered, then they will interpreted as the minimum and maximum age uncertainty. (see examples) By default, time flows from higher to lower numbers. You can change this behavior by setting \code{revtime=FALSE}. Either a \code{bin} or \code{age} argument is mandatory. 
#' @param revtime \code{(logical)}: Should time be reversed? 
#' @param na.rm \code{(logical)}: Should taxa that have no valid FADs or LADs (due to \code{NA} entries) be removed from the output?
#' @param diffbin \code{(logical)}: Difference-based duration for discreet time (only applicabble to cases when \code{bin} is provided). If set to \code{TRUE}, single-interval taxa will \code{0} durations. Setting this argument to \code{FALSE} will add \code{1} to the durations of all taxa.
#' @examples 
#' data(corals)
#' 
#' # binned data
#'   flBinned <- fadlad(corals, tax="genus", bin="stg")
#' 
#' # using basic bin lengths
#'   flDual <- fadlad(corals, tax="genus", age=c("max_ma", "min_ma"))
#'
#' # single age esimate 
#'   data(stages)
#'   corals$mid <- stages$mid[corals$stg]
#'   flSingle <- fadlad(corals, tax="genus", age="mid")
#' 
#' 
#' @rdname fadlad
#' @return A data.frame, with rows corresponding to \code{tax} entries.
#' @export
fadlad<-function(x, tax, bin=NULL,age=NULL, revtime=FALSE, na.rm=TRUE, diffbin=TRUE){
	# for the prototype
#	x <- corals
#	tax<- "genus"
#	bin <- "stg"
	# correct data?
	if(!is.matrix(x) & !is.data.frame(x) ) stop("Invalid 'x' argument.")
	
	if(is.null(bin) & is.null(age)) stop("You need to provide either a 'bin' or an 'age' variable.")
	if(!is.null(bin) & !is.null(age)){
		warning("You provided both a 'bin' and an 'age' variable. Only 'bin' is used.")
		age <- NULL
	}
	if(!tax%in%colnames(x))stop("Invalid 'tax' argument.")


	# binned data
	if(!is.null(bin)){
		if(sum(bin%in%colnames(x))!=length(bin)) stop("Invalid 'bin' argument(s).")
		for(i in 1:length(bin)) if(!is.numeric(x[,bin[i], drop=TRUE])) stop("'bin' must be a numeric variable.")
		binned <- TRUE
		time <- bin
	# age data
	}else{
		if(sum(age%in%colnames(x))!=length(age)) stop("Invalid 'age' argument(s).")
		for(i in 1:length(age)) if(!is.numeric(x[,age[i], drop=TRUE])) stop("'age' must be a numeric variable.")
		binned <- FALSE
		if(!diffbin) warning("The 'diffbin=FALSE' option was discarded, because 'age' was provided.")
		time <- age
		x[, time] <- -x[,time]
	}

	# demote the tax column to character
		taxVar <- as.character(x[,tax, drop=TRUE])
	
	# coerce bin to matrix
		timeVar <- as.matrix(x[,time, drop=TRUE])

	# are ages?
	if(revtime){
		for(i in 1:length(time)) timeVar[,i] <- -timeVar[,i, drop=TRUE]
	}
	
	# single iteration
	tempFL<-tapply(INDEX=taxVar, X=1:nrow(x), FUN=function(w){
		# taxon specific part
		taxDat <- timeVar[w,]
		
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
		fl$duration <- abs(fl[,"FAD", drop=TRUE]-fl[,"LAD", drop=TRUE])

		# only for binned data, and when instructed!
		if(!diffbin & binned){
			fl$duration <- fl$duration+1
		}

	# if ages are used, the default is to revert time
	if(!binned){
		fl$FAD <- -fl$FAD
		fl$LAD <- -fl$LAD
	}

	# add the species names
	rownames(fl)<- names(tempFL)
		
	if(na.rm){
		fl<-fl[!is.na(fl[,"FAD", drop=TRUE]),]
	}
	
	if(any(""==rownames(fl))){
		warning("Taxon name <\"\"> (empty quotes) is detected!")
	}

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
#' @param x \code{(data.frame)} The data frame containing fossil occurrences.
#' 
#' @param tax \code{(character)} The variable  name of the occurring taxa (variable type: \code{factor} or \code{character}).
#' 
#' @param bin \code{(character)} The variable name of the time slice numbers of the particular occurrences (variable type: \code{numeric}). Bin numbers should be in ascending order,can contain \code{NA}s, it can start from a number other than 1 and must not start with 0.
#' @param noNAStart \code{(logical)} Useful when the dataset does not start from bin \code{1}. Then \code{noNAStart=TRUE} will cut the first part of the resulting table, 
#' 						so the first row will contain the estimates for the lowest bin number.
#' 
#' @param fl \code{(matrix} or \code{data.frame}). If so desired, the function can be run on an FAD-LAD dataset, output by the \code{\link{fadlad}} function. 
#' 
#' @param method \code{(character)} Either \code{"forward"} or \code{"backward"}.
#' 
#' @return A numeric matrix of survivorship probabilities.
#' @examples
#' data(corals)
#' surv<-survivors(corals, tax="genus", bin="stg", method="forward")
#' 
#' # plot
#' data(stages)
#' tsplot(stages, shading="series", boxes="sys", xlim=c(260,0), 
#'   ylab="proportion of survivors present", ylim=c(0.01,1),plot.args=list(log="y"))
#'   
#' for(i in 1:ncol(surv)) lines(stages$mid, surv[,i])
#' 
#' @export
survivors<-function(x, tax="genus", bin="stg", method="forward", noNAStart=FALSE,fl=NULL){
	
	# calculate FAD-LAD matrix first, if not provided
	if(is.null(fl)){
		fl <- fadlad(x, tax, bin)
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
		
		sProbList<-sapply(inter, FUN=function(w){
			# who are there?
			bOrig <- fl$FAD<=w & w<=fl$LAD
			original<-fl[bOrig,]
			subVector <- w:max(inter)
			curSurv<-sapply(subVector, FUN=function(y){
				return(sum(original$FAD<=y & y<=original$LAD))
				
			})
			
			retVal <- template
			retVal[subVector]<-curSurv/sum(bOrig)
			return(retVal)
		})
	}
	
	if(method=="backward"){
		
		sProbList<-sapply(inter, FUN=function(w){
			# who are there?
			bOrig <- fl$FAD<=w & w<=fl$LAD
			original<-fl[bOrig,]
			subVector <- min(inter):w
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

