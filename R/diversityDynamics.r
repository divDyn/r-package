#' Time series from metrics of diversity dynamics 
#' 
#' This function calculates various metrics from occurrence datasets in the form of time series.
#' 
#' The following variables are produced:
#'
#' \code{bin}: Bin number, or the numeric identifier of the bin.
#'
#' \code{tThrough}: Number of through-ranging taxa, taxa that have first occurrences before, and last occurrences after the focal bin.
#'
#' \code{tOri}: Number of originating taxa, taxa that have first occurrences in the focal bin, and last occurrences after it.
#'
#' \code{tExt}: Number of taxa getting extinct. These are taxa that have first occurrences before the focal bin, and last occurrences in it.
#'
#' \code{tSing}: Number of stratigraphic singleton (single-interval) taxa, taxa that only occur in the focal bin.
#'
#' \code{t2d}: Number of lower two timers (Alroy, 2008; 2014), taxa that are present in the \emph{i}-1th and the ith bin (focal bin). 
#'
#' \code{t2u}: Number of upper two timers (Alroy, 2008; 2014), taxa that are present in the \emph{i}th (focal) and the \emph{i}+1th bin. (Alroy, 2008; 2014)
#'
#' \code{tGFu}: Number of upper gap-fillers (Alroy, 2014), taxa that occurr in bin \emph{i}+2 and \emph{i}-1, but were not found in \emph{i}+1. (Alroy, 2014)
#'
#' \code{tGFd}: Number of lower gap-fillers (Alroy, 2014), taxa that occurr in bin \emph{i}-2 and \emph{i}+1, but were not found in \emph{i}-1. (Alroy, 2014)
#'
#' \code{t3}: Number of three timer taxa (Alroy, 2008; 2014), present in bin \emph{i}-1, \emph{i}, and \emph{i}+1. (Alroy, 2008; 2014)
#'
#' \code{tPart}: Part timer taxa (Alroy, 2008; 2014), present in bin \emph{i}-1,and \emph{i}+1, but not in bin \emph{i}. 
#'
#' \code{extProp}: Proportional extinctions including single-interval taxa: \emph{(tExt + tSing) / (tThrough + tOri + tExt + tSing)}.
#'
#' \code{oriProp}: Proportional originations including single-interval taxa:  \emph{(tOri + tSing) / (tThrough + tOri + tExt + tSing)}.
#' 
#' \code{extPC}: Per capita extinction rates of Foote (1999). \emph{-log(tExt/(tExt + tThrough))}.  Values are not normalized with bin lengths. Similar equations were used by Alroy (1996) but without taking the logarithm.
#'
#' \code{oriPC}: Per capita origination rates of Foote (1999). \emph{-log(tOri/(tOri + tThrough))}. Values are not normalized with bin lengths. Similar equations were used by Alroy (1996) but without taking the logarithm.
#'
#' \code{ext3t}: Three-timer extinction rates of Alroy (2008). \emph{log(t2d/t3)}.
#'
#' \code{ori3t}: Three-timer origination rates of Alroy (2008). \emph{log(t2u/t3)}.
#'
#' \code{extC3t}: Corrected three-timer extinction rates of Alroy (2008). \emph{ext3t[\emph{i}] + log(samp3t[\emph{i}+1])}.
#'
#' \code{oriC3t}: Corrected three-timer origination rates of Alroy (2008). \emph{ori3t[\emph{i}] + log(samp3t[\emph{i}-1])}.
#'
#' \code{divSIB}: Sampled-in-bin diversity (richness), the number of genera sampled in the focal bin.
#'
#' \code{divCSIB}: Corrected sampled-in-bin diversity (richness). \emph{divSIB/samp3t*totSamp3t}, where \emph{totSamp3t} is total three-timer sampling completeness of the dataset (Alroy, 2008). 
#'
#' \code{divBC}: Boundary-crosser diversity (richness), the number of taxa with ranges crossing the boundaries of the interval. \emph{tExt + tOri + tThrough}.
#'
#' \code{divRT}: Range-through diversity (richness), all taxa in the interval, based on the range-through assumption. \emph{(tSing + tOri + tExt + tThrough)}.
#'
#' \code{sampRange}: Range-based sampling probability, without observed range end-points (Foote), \emph{(divSIB - tExt - tOri- t-Sing)/tThrough}
#'
#' \code{samp3t}: Three-timer sampling completeness of Alroy (2008). \emph{t3/(t3+tPart)}
#'
#' \code{extGF}: Gap-filler extinction rates of Alroy(2014). \emph{log((t2d + tPart)/(t3+tPart+tGFd))}
#'
#' \code{oriGF}: Gap-filler origination rates of Alroy(2014). \emph{log((t2u + tPart)/(t3+tPart+tGFd))}
#'
#' \code{E2f3}: Second-for-third extinction propotions of Alroy (2015). As these metrics are based on an algorithmic approach, for the equations please refer to the Alroy (2015, p. 634, right column and Eq. 4)). See source code (\url{https://github.com/divDyn/r_package}) for the exact implementation, found in the \code{Metrics} function in the diversityDynamics.R file.
#'
#' \code{O2f3}: Second-for-third origination propotions of Alroy (2015). Please see \code{E2f3}.
#'
#' \code{ext2f3}: Second-for-third extinction rates (based on Alroy, 2015). Transformed to the usual rate form with \emph{log(1/(1-E2f3))}.
#'
#' \code{ori2f3}: Second-for-third origination rates (based on Alroy, 2015). Transformed to the usual rate form with \emph{log(1/(1-O2f3))}.
#' 
#' \strong{References:}
#'
#' Foote, M. (1999) Morphological Diversity In The Evolutionary Radiation Of Paleozoic and Post-Paleozoic Crinoids. Paleobiology 25, 1–115. doi:10.1017/S0094837300020236.
#'
#' Alroy, J. (2008) Dynamics of origination and extinction in the marine fossil record. Proceedings of the National Academy of Science 105, 11536-11542. doi: 10.1073/pnas.0802597105
#'
#' Alroy, J. (2014) Accurate and precise estimates of origination and extinction rates. Paleobiology 40, 374-397. doi: 10.1666/13036
#'
#' Alroy, J. (2015) A more precise speciation and extinction rate estimator. Paleobiology 41, 633-639. doi: 10.1017/pab.2015.26
#'
#' @param x \code{(data.frame)} Fossil occurrence table.
#' 
#' @param tax \code{(character)} Variable name of the occurring taxa (variable type: \code{factor} or \code{character} - such as \code{"genus"}
#' 
#' @param bin \code{(character)} Variable name of the discrete bin numbers of the occurrences. This variable should be \code{numeric}. Time flows from lower to higher values by default. Use \code{revtime} to reverse this order.
#' @param age \code{(character)} Variable name of the ages of the occurrences that will be sliced with the \code{slice} function using the intervals provided in \code{breaks}. This variable should be \code{numeric}. Time flows from higher to lower values by default. Use \code{revtime} to reverse this order.
#' @param revtime \code{(logical)} Argument for reversing the default direction of time. Setting this argument to \code{TRUE} will make time flow from higher to lower values when \code{bin} is used, and from lower to higher values, when \code{age} is given. CAUTION: Failing to set this argument properly can make originations become extinctions and vice versa! 
#' @param breaks \code{(numeric)} If \code{NULL} (default) the used values in the \code{bin} variable will designate independent time slices that follow each other in succession. If a vector is provided, than the numeric entries in \code{bin} will be binned similarly to the \code{\link[graphics]{hist}} or \code{\link[base]{cut}} function. The order of elements in this vector is arbitrary.
#' @param noNAStart (logical) Useful when the entries in the \code{bin} variable do not start from bin no. 1, but positive integer bin numbers are provided. Then \code{noNAStart=TRUE} will cut the first part of the resulting table, so the first row will contain the estimates for the lowest bin number. In case of positive integer bin identifiers, and if \code{noNAStart=FALSE}, the index of the row will be the bin number. 
#' 
#' @param data.frame \code{(logical)} Should the output be a \code{data.frame} or a \code{matrix}?
#' 
#' @param om \code{(character)} The \code{om} argument of the \code{omit()} function. If set to \code{NULL} (default), then no occurrences will be omitted before the execution of the function.
#' @param filterNA \code{(logical)} The \code{filterNA} parameter of the \code{\link{omit}} function.
#' @param coll \code{(character)} The variable name of the collection identifiers. (optional, only for use with the internal \code{\link{omit}} function)
#' @param ref \code{(character)} The variable name of the reference identifiers. (optional, only for use with the internal \code{\link{omit}} function)
#' @examples
#'	# import data
#'	  data(corals)
#'	  data(stages)
#'
#'	# calculate metrics of diversity dynamics
#'    dd <- divDyn(corals, tax="genus", bin="stg")
#'
#'	# plotting
#'	  tsplot(stages, shading="series", boxes="sys", xlim=c(260,0), 
#'	    ylab="range-through diversity (genera)", ylim=c(0,230))
#'	  lines(stages$mid, dd$divRT, lwd=2)
#' 
#'  # with omission of single reference taxa  
#'    ddNoSing <- divDyn(corals, tax="genus", bin="stg", om="ref", ref="reference_no")
#'    lines(stages$mid, ddNoSing$divRT, lwd=2, col="red")
#'
#'  # using the estimated ages (less robust) - 10 million years
#'    # mean ages
#'    corals$me_ma <- apply(corals[, c("max_ma", "min_ma")], 1, mean)
#'    # ages reverse the direction of time! set ages to TRUE in this case
#'    ddRadio10 <- divDyn(corals, tax="genus", age="me_ma", 
#'		breaks=seq(250,0,-10))
#'    lines(ddRadio10$me_ma, ddRadio10$divRT, lwd=2, col="green")
#'       
#'  # legend
#'    legend("topleft", legend=c("all", "no single-ref. taxa", "all, estimated ages"), 
#'      col=c("black", "red", "green"), lwd=c(2,2,2), bg="white")
#'    
#'
#' @export
#' @return A data.frame object, with every row corresponding to a time bin. 
divDyn <- function(x, tax, bin=NULL, age=NULL, revtime=FALSE, breaks=NULL, coll=NULL, ref=NULL, om=NULL,noNAStart=FALSE, data.frame=TRUE, filterNA=FALSE)
{
	
#	# tral run
#	x<-corals
#	bin<-"stg"
#	tax<- "genus"
#	ages <- FALSE
#	breaks<-NULL

	# 1. Binning organization
	if(is.null(bin) & is.null(age)) stop("You have to provide either a 'bin' or an 'age' argument.")

	if(!is.null(bin) & !is.null(age)) stop("The 'bin' and 'age' arguents cannot be used together.")

	# calculation offset. Bins are shifted for the time of the calculations. 
	# Used to make sure that the Rcpp counting algorithm works accurately.
	offset <- 4

	# by default do not use factors:
	useFactor <-FALSE

	# 1a. the bin variable is used
	if(!is.null(bin)){
		# the breaks argument cannot be used
		if(!is.null(breaks)) stop("The 'breaks' argument cannot be used when discreet bins are provided. ")

		# for testing, checking and binning calculations
		binVarFull <- x[,bin, drop=TRUE]
	
		# is the variable the right type?
		if(!is.numeric(binVarFull) & !is.factor(binVarFull)){
			stop("The 'bin' variable is not numeric or a factor.")
		}
	
		#if bin is a factor, convert it to numeric here
		if(is.factor(binVarFull)){
			# conserve the levels
			factLevs <- levels(binVarFull)

			#save the names
			binVarFull<-as.numeric(binVarFull)

			useFactor<- TRUE
		}
		
		# reverse time at this stage
		if(revtime){
			x[, bin] <- -binVarFull
			binVarFull <- -binVarFull
		}

		# get levels
		lev <- sort(unique(binVarFull)) # can NAs get past this? - yes, but the code below handles them

		# check whether there are enough entries?
		if(length(lev)<3) stop("You have less than 3 bins.")

		# if there are integer entries coerce continuous numbering for the output
		# and use only shifts of values, rather than replacement (much faster)
		if(sum(lev%%1==0)==length(lev)) {
			minLev<-min(lev)
			maxLev<-max(lev)
			lev<-minLev:maxLev

			# shift the ranks so 5 should be the first index
			y<- binVarFull-minLev+1+offset
		
		# used for specific bin entries e.g. -300, -295
		# this might have to be forced by factorization! 
		}else{

			ranks <- rank(lev)+offset
			names(ranks) <- lev
	
			# the NAs
			notNA <- !is.na(binVarFull)
			y<-rep(NA, length(binVarFull))
	
			# create the new bin vector
			y[notNA]<-ranks[as.character(binVarFull[notNA])]
			# y became bin
			# and lev became lev
		}
		
		# copy in the new variable
		x[,bin] <- y

		# use a unified variable name from now on
		time <- bin

	}

	# 1b. age data
	if(!is.null(age)){
		# the breaks argument must be used
		if(is.null(breaks)) stop("If you use ages, you have to use provide a 'breaks' argument.")
		if(noNAStart) warning("You provided ages, 'noNAStart' will be ignored.")
		
		# for testing, checking and binning calculations
		ageVarFull <- x[,age, drop=TRUE]
	
		# is the variable the right type?
		if(!is.numeric(ageVarFull)){
			stop("The 'age' variable is not numeric.")
		}

		# reverse time at this stage by default, opposite the way as for bins
		if(!revtime){
			x[, age] <- -x[,age, drop=TRUE]
			ageVarFull <- -ageVarFull
			breaks <- -breaks
		}

		# call the automatic binning function
		aubi <- slice(x[,age, drop=TRUE], breaks=breaks, offset=offset, ts=FALSE, revtime=FALSE)
		x[,age] <- aubi$slc
		lev <- aubi$lev

		# use a unified variable name from now on
		time <- age
	}
	

	# in case the scale defined by breaks has younger parts than the data
	# the number of de facto bins for which metrics will be calculated
	maxVal<- max(x[,time], na.rm=TRUE)-offset

	# the omission phase - does this work?
	if(!is.null(om)){
		x<-x[!omit(x, tax=tax, bin=time,om=om, ref=ref, coll=coll, filterNA=filterNA),]
	}
	
	# sub dataset
	subDat <- unique(x[,c(tax, time)])
	
	# omit NAs
	bNeed<- !(is.na(subDat[,tax, drop=TRUE]) | is.na(subDat[,time, drop=TRUE]))

	# taxon vars
	taxVar<-subDat[bNeed, tax, drop=TRUE]
	binVar<-subDat[bNeed, time, drop=TRUE]

	if(any(""==taxVar)){
		warning("Taxon name <\"\"> (empty quotes) is detected!")
	}

	# factorize
	taxVar <- as.numeric(factor(taxVar))
	
	# shift to 0 indices
	taxVar<-taxVar-1
	binVar<-binVar-1
	
	# here come the counts variables
	counts <- Counts(taxVar,binVar)

	# the metrics
	metrics <- Metrics(counts)
	
	# replace infinites with NAs
	metrics[!is.finite(metrics)]<-NA

	# concatenate the columns
	dCountsAndMetrics<-cbind(counts, metrics)
		
	# delete the safety offset
	dCountsAndMetrics<- dCountsAndMetrics[(offset+1):nrow(dCountsAndMetrics), ]

	# cbind all together
	dCountsAndMetrics<-cbind(time=lev[1:maxVal], dCountsAndMetrics)
	colnames(dCountsAndMetrics)[colnames(dCountsAndMetrics)=="time"] <- time
	
	#create the returned table
	columns <- c(
		time,
		"t2d",
		"t2u",
		"t3",
		"tPart",
		"tGFd",
		"tGFu",
		"tSing",
		"tOri",
		"tExt",
		"tThrough",
		"divSIB",
		"divCSIB",
		"divRT",
		"divBC",
		"extProp",
		"oriProp",
		"extPC",
		"oriPC",
		"ext3t",
		"ori3t",
		"extC3t",
		"oriC3t",
		"extGF",
		"oriGF",
		"E2f3",
		"O2f3",
		"ext2f3",
		"ori2f3",
		"samp3t",
		"sampRange"
	)
	dCountsAndMetrics<-dCountsAndMetrics[,columns]

	# coerce NAs in variables and intervals where they were not supposed to be
		# first two rows
		dCountsAndMetrics[1,c("t2d","t3","tPart","tGFd","tGFu", "tThrough", "tExt", "divBC", "O2f3")] <- NA
		dCountsAndMetrics[2,c("tGFd","oriGF", "ori2f3", "O2f3")] <- NA
	
		# last two rows
		nEnd<-nrow(dCountsAndMetrics)
		dCountsAndMetrics[nEnd,c("t2u","t3","tPart","tGFd","tGFu", "ext2f3", "tOri", "tThrough", "E2f3")] <- NA
		dCountsAndMetrics[nEnd-1,c("tGFu","extGF", "ext2f3", "E2f3")] <- NA

	# bin/age variable in the output table
	# reverse time for bins if asked so
	if(!is.null(bin)){
		if(revtime){
			dCountsAndMetrics[,bin]<- -dCountsAndMetrics[,bin]
		}
	# reverse time by default for ages
	}else{
		if(!revtime){
			dCountsAndMetrics[,age]<- -dCountsAndMetrics[,age]
		}
	}

	# by default,use an output data.frame 
	if(data.frame){				
		dCountsAndMetrics<-as.data.frame(dCountsAndMetrics, stringsAsFactors=F)
		if(useFactor) dCountsAndMetrics[,bin] <- factLevs
	}
	
	#want to see the NA's at the beginnning? 
	# only use noNAStart when the simple indices are used	
	if(!is.null(bin) & !useFactor){
		if(!noNAStart){
			# only if it is applicable - e.g. 
			firstBin <- min(dCountsAndMetrics[,bin]) 
			
			# (when the time series does not start with bin 1 and bins are positive integers)
			if(sum((dCountsAndMetrics[,bin]%%1)==0)==nEnd & firstBin>1){
				empty <- matrix(NA, nrow=firstBin-1, ncol=ncol(dCountsAndMetrics))
				colnames(empty) <- columns
				
				# insert the NAs
				dCountsAndMetrics<-rbind(empty,dCountsAndMetrics)
				
				# make the first column complete so it can be used as rownames
				dCountsAndMetrics[1:(firstBin-1),bin] <- 1:(firstBin-1)
			
			}
			
		}
	}
	
	

	# after this the rownames should be alright, copy bin to that, so it can be used in subsampling
	rnames <- as.character(dCountsAndMetrics[,time])

	# weird problem for ages: sometimes during the character conversion and rounding 
	# the bin IDs become identical...
	bDuplum <- duplicated(rnames)
	if(sum(bDuplum)>0){
		# append some junk to make them different
		rnames[bDuplum] <- paste(rnames[bDuplum], 1:sum(bDuplum), sep="")
	}
	
	# and copy the rownames
	rownames(dCountsAndMetrics) <- rnames

	# use factor levels as rownames when the data.frame option is FALSE
	if(!data.frame & useFactor) rownames(dCountsAndMetrics)<- factLevs
		
	#return the table					
	return(dCountsAndMetrics)
}

# function version 2.0
Counts <- function(tax, bin){
	counts<- .Call('_divDyn_Counts', PACKAGE = 'divDyn', tax, bin)
	colnames(counts)<- c(
		"t1",
		"t2d",
		"t2u",
		"t3",
		"tPart",
		"tGFd",
		"tGFu",
		"s1d",
		"s2d",
		"s3d",
		"s1u",
		"s2u",
		"s3u",
		"tSing",
		"tOri",
		"tExt",
		"tThrough",
		"divSIB"
	)
	return(counts)

}

Metrics<- function(counts){
##########################################
	# the metrics
	metNames<-c("divCSIB",
		"divRT",
		"divBC",
		"extProp",
		"oriProp",
		"extPC",
		"oriPC",
		"ext3t",
		"ori3t",
		"extC3t",
		"oriC3t",
		"extGF",
		"oriGF",
		"E2f3",
		"O2f3",
		"ext2f3",
		"ori2f3",
		"samp3t",
		"sampRange")
	
	#container
	metrics<-matrix(NA, ncol=length(metNames), nrow=nrow(counts))
	colnames(metrics)<-metNames
	
	#BC diversity
	metrics[,"divBC"] <- counts[,"tThrough", drop=TRUE]+counts[,"tExt", drop=TRUE]

	# total diversity
	sib. <- counts[,"divSIB", drop=TRUE]-counts[,"tSing", drop=TRUE]
	div <- metrics[,"divBC", drop=TRUE]+counts[,"tOri", drop=TRUE]                  
	metrics[,"divRT"]<-div+counts[,"tSing", drop=TRUE]
	
	#Sampling parameters
	#sampling probability (Foote, 2000)
	obs <- sib.-counts[,"tExt", drop=TRUE]-counts[,"tOri", drop=TRUE]
	metrics[,"sampRange"] <- obs/counts[,"tThrough", drop=TRUE] 

	#Three-timer sampling completeness (Alroy, 2008)	
	metrics[,"samp3t"] <- counts[,"t3", drop=TRUE]/(counts[,"t3", drop=TRUE]+counts[,"tPart", drop=TRUE])
			
			#Sampling completeness of the entire time series
			nTot3tSampComp<- sum(counts[,"t3", drop=TRUE], na.rm=T)/(sum(counts[,"t3", drop=TRUE], na.rm=T)+sum(counts[,"tPart", drop=TRUE], na.rm=T))
			
	# proportional extinctions
		metrics[,"extProp"]<-(counts[,"tExt", drop=TRUE]+counts[,"tSing", drop=TRUE])/metrics[,"divRT", drop=TRUE]
		metrics[,"oriProp"]<-(counts[,"tOri", drop=TRUE]+counts[,"tSing", drop=TRUE])/metrics[,"divRT", drop=TRUE]
	
	#Foote (2000) rates
		metrics[,"extPC"]<- -log(counts[,"tThrough", drop=TRUE]/(counts[,"tThrough", drop=TRUE]+counts[,"tExt", drop=TRUE]))
		metrics[,"oriPC"]<- -log(counts[,"tThrough", drop=TRUE]/(counts[,"tThrough", drop=TRUE]+counts[,"tOri", drop=TRUE]))
		
	#Three-timer rates by Alroy (2008)
		#uncorrected:
		metrics[,"ext3t"] <- log(counts[,"t2d", drop=TRUE]/counts[,"t3", drop=TRUE]) # two-timer/three-timer ratio (bottom)
		metrics[,"ori3t"] <- log(counts[,"t2u", drop=TRUE]/counts[,"t3", drop=TRUE]) # two-timer/three-timer ration (top)
		
		#corrected:
			#extinction rates:
				thtSampCompNext <- c(metrics[2:nrow(counts),"samp3t"],NA) # Sampling probability in subsequent bin (BIN5)
				metrics[,"extC3t"] <- metrics[,"ext3t", drop=TRUE] + log(thtSampCompNext)
				#nC3tExt[nC3tExt<0] <- 0 # omit negative values
			
			#origination rates:
				thtSampCompPrev <- c(NA, metrics[1:(nrow(counts)-1),"samp3t", drop=TRUE]) # Sampling probability in previous bin (BIN5)
				metrics[,"oriC3t"] <- metrics[,"ori3t", drop=TRUE] + log(thtSampCompPrev)
				#nC3tOri[nC3tOri<0] <- 0 #omit negative values
	 
	
	#Gap filler rates by Alroy(2014)
		metrics[,"extGF"]<-log((counts[,"t2d", drop=TRUE]+counts[,"tPart", drop=TRUE])/(counts[,"t3", drop=TRUE]+counts[,"tPart", drop=TRUE]+counts[,"tGFu", drop=TRUE]))
		metrics[,"oriGF"]<-log((counts[,"t2u", drop=TRUE]+counts[,"tPart", drop=TRUE])/(counts[,"t3", drop=TRUE]+counts[,"tPart", drop=TRUE]+counts[,"tGFd", drop=TRUE]))
	
	# second- for third (Alroy, 2015)
	
	#	substituteF<-function(y){
	#		# in such cases ... 
	#		# reversed ordering
	#		first <- y[1]<y[2] & y[2]<y[3]
	#		
	#		# s2 is minimum
	#		second <- y[2]<y[1] & y[2]<y[3]
	#		
	#		# s2 is maximum
	#		third <- y[2]>y[1] & y[2]>y[3]
	#		
	#		# substituting the second lowest of the three counts for s3
	#		if(first | second | third){
	#			sort(y)[2]
	#		}else{
	#			y[3]
	#		}
	#	}
		# this is not enough, because the proportion can be still very negative,
		# he correctly says at the end of the paragraph, it should be a second in the order
	
		# extinction proportions (Eq. 4 in Alroy, 2015)
		sSubD<-apply(counts[,c("s1d","s2d","s3d")],1,function(y) sort(y)[2])
		metrics[,"E2f3"] <- (counts[,"s1d", drop=TRUE]-sSubD)/(counts[,"t2d", drop=TRUE]+counts[,"tPart", drop=TRUE])
		
		# origination proportions
		sSubU<-apply(counts[,c("s1u","s2u","s3u")],1,function(y) sort(y)[2])
		metrics[,"O2f3"] <- (counts[,"s1u", drop=TRUE]-sSubU)/(counts[,"t2u", drop=TRUE]+counts[,"tPart", drop=TRUE])
	
		# transform to classical rate form
		metrics[,"ext2f3"] <-log(1/(1-metrics[,"E2f3", drop=TRUE]))
		metrics[,"ori2f3"] <-log(1/(1-metrics[,"O2f3", drop=TRUE]))
		
	
	#corrected sampled-in-bin diversity
		metrics[,"divCSIB"]<-counts[,"divSIB", drop=TRUE]*nTot3tSampComp/metrics[,"samp3t", drop=TRUE]
	

	return(metrics)
}


#' Omission of taxa that have a poor occurrence record
#' 
#' Function to quickly omit single-collection and single-reference taxa.
#' 
#' The function returns a \code{logical} vector, with a value for each row. \code{TRUE} values indicate rows to be omitted, \code{FALSE} values indicate rows to be kept. The function is embedded in the \code{\link{divDyn}} function, but can be called independently.
#' 
#' @param bin \code{(character)} The name of the bin variable (has to be \code{numeric} for the function to run). For time series, this is the time slice variable.
#' 
#' @param tax \code{(character)} The name of the taxon variable.
#' 
#' @param x \code{(data.frame)} Occurrence dataset, with \code{bin}, \code{tax} and \code{coll} as column names.
#' @param coll \code{(character)} The variable name of the collection identifiers. 
#' @param ref \code{(character)} The variable name of the reference identifiers (optional). 
#' @param om \code{(character)} The type of omission. \code{"coll"} omits occurrences of taxa that occurr only in one collection. \code{"ref"} omits occurrences of taxa that were described only in one reference. \code{"binref"} will omit the set of single reference taxa that were described by more than one references, but appear in only one reference in a time bin.
#' @param filterNA \code{(logical)} Additional entries can be added to influence the dataset that might not have reference or collection information (\code{NA} entries). These occurrences are treated as single-collection or single-reference taxa if the \code{na.rm} argument is set to \code{FALSE} (default). Setting this argument to \code{TRUE} will keep these entries. (see example)
#' 
#' @return A logical vector. 
#' @examples
#' # omit single-reference taxa
#'   data(corals)
#'   data(stages)
#'   toOmit <- omit(corals, bin="stg", tax="genus", om="ref", ref="reference_no")
#'   x <- corals[!toOmit,]
#' 
#' # within divDyn
#'   # plotting
#'	  tsplot(stages, shading="series", boxes="sys", xlim=c(260,0), 
#'	    ylab="range-through diversity (genera)", ylim=c(0,230))
#'   # multiple ref/slice required
#'   ddNoSing <- divDyn(corals, tax="genus", bin="stg", om="binref", ref="reference_no")
#'   lines(stages$mid, ddNoSing$divRT, lwd=2, col="red")
#'
#'   # with the recent included (NA reference value)
#'   ddNoSingRec <- divDyn(corals, tax="genus", bin="stg",
#'     om="binref", filterNA=TRUE,ref="reference_no")
#'   lines(stages$mid, ddNoSingRec$divRT, lwd=2, col="blue")
#'   
#'   # legend
#'   legend("topleft", legend=c("no single-ref. taxa", 
#'     "no single-ref. taxa,\n with recent"), 
#'     col=c("red", "blue"), lwd=c(2,2))
#' @export
omit <- function(x, om="ref", tax="genus", bin="bin", coll=NULL, ref=NULL, filterNA=FALSE){


	if(!om%in%c("coll", "ref","binref")) stop("Invalid om argument.")
	
	if(om=="coll"){
		if(is.null(coll)) stop("You have to provide a 'coll' argument to omit by collection.")
		if(is.null(x[[coll]])) stop("Invalid collection variable ('coll'). ")
		# omit multiple occ rows of same tax (genus) and same coll
		nonDupl <- !duplicated(x[,c(tax, coll)])
		
		# which taxa come from just one collection?
		tabSing <- table(x[nonDupl,tax, drop=TRUE])
		
		# single collection taxa
		singTax<- names(tabSing)[tabSing==1]
		
		# omit the single collection taxa
		boolEnd<-x[,tax, drop=TRUE]%in%singTax
		
		# if na.rm TRUE than do not omit NA collection stuff
		if(filterNA){
			boolEnd <- boolEnd & !is.na(x[,coll, drop=TRUE])
		}
	}
	
	if(om=="ref"){
		if(is.null(ref)) stop("You have to provide a 'ref' argument to omit by reference")
		if(is.null(x[[ref]])) stop("Invalid reference variable ('ref'). ")
		
		# omit multiple occ rows of same tax (genus) and same coll
		nonDupl <- !duplicated(x[,c(tax, ref)])
		
		# which taxa come from just one collection?
		tabSing <- table(x[nonDupl,tax, drop=TRUE])
		
		# single collection taxa
		singTax<- names(tabSing)[tabSing==1]
		
		# omit the single collection taxa
		boolEnd<-x[,tax, drop=TRUE]%in%singTax
		
		# if na.rm TRUE than do not omit NA collection stuff
		if(filterNA){
			boolEnd <- boolEnd & !is.na(x[,ref, drop=TRUE])
		}
	}
	
	
	
	if(om=="binref"){
		if(is.null(ref)) stop("You have to provide a 'ref' argument to omit by reference.")
		if(is.null(x[[ref]])) stop("Invalid reference variable ('ref'). ")
		if(is.null(bin)) stop("You have to provide a 'bin' argument to omit by bin-reference.")
		if(is.null(x[[bin]])) stop("Invalid bin variable ('bin'). ")
		

		nonDupl <- !duplicated(x[,c(tax, ref, bin)])
		
		activeDat <- x[nonDupl, ]
		
		tap<-tapply(INDEX=activeDat[,bin, drop=TRUE], X=activeDat[,tax, drop=TRUE], function(y){
					
			tabSing <- table(y)
		
			nonSingTax<- names(tabSing)[tabSing!=1]
			return(nonSingTax)
			
		
		})
		# list of taxa that do not occur in in only reference/slice
		taxaMoreThanOne <- unique(unlist(tap))
		
		tap2<-(1:nrow(x))[!x[, tax, drop=TRUE]%in%taxaMoreThanOne]
		
	
		boolEnd<-rep(FALSE, nrow(x))
		boolEnd[tap2]<-TRUE
		
		# if na.rm TRUE than do not omit NA reference stuff
		if(filterNA){
			boolEnd <- boolEnd & !is.na(x[,ref, drop=TRUE])
		}
		
	}
	
	return(boolEnd)

}



