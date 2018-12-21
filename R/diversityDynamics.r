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
#' \code{tOri}: Number of originating taxa, taxa that have first occurrences in the focal bin, and last occurrences after.
#'
#' \code{tExt}: Number of taxa getting extinct. These are taxa that have first occurrences before the focal bin, and last occurrences after it.
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
#' \code{sampRange}: Range-based sampling probability (Foote), \emph{(divSIB - tExt - tOri- t-Sing)/tThrough}
#'
#' \code{samp3t}: Three-timer sampling completeness of Alroy (2008). \emph{t3/(t3+tPart)}
#'
#' \code{extGF}: Gap-filler extinction rates of Alroy(2014). \emph{log((t2u + tPart)/(t3+tPart+tGFd))}
#'
#' \code{oriGF}: Gap-filler origination rates of Alroy(2014). \emph{log((t2u + tPart)/(t3+tPart+tGFd))}
#'
#' \code{E2f3}: Second-for-third extinction propotions of Alroy (2015). As these metrics are based on an algorithmic approach, for the equations please refer to the Alroy (2015, p. 634, right column and Eq. 4)). See source code (\url{http://www.github.com/adamkocsis/divDyn}) for the exact implementation, found in the \code{Metrics} function in the diversityDynamics.R file.
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
#' @param dat \code{(data.frame)} Fossil occurrence table.
#' 
#' @param tax \code{(character)} Variable name of the occurring taxa (variable type: \code{factor} or \code{character} - such as \code{"genus"}
#' 
#' @param bin \code{(character)} Variable name of the bin numbers of the occurrences. This variable should be \code{numeric}. 
#'
#' @param ages \code{(logical)} Argument for setting the direction of time in \code{bin}. The default setting, \code{FALSE} will execute the function with higher numbers corresponding to later intervals. \code{TRUE} will reverse the direction of values in the \code{bins} variable, making originations extinctions and vice versa (use this for age input!, e.g 400 Ma, see examples). 
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
#'	  tsplot(stages, shading="series", boxes="per", xlim=c(260,0), 
#'	    ylab="range-through diversity (genera)", ylim=c(0,230))
#'	  lines(stages$mid, dd$divRT, lwd=2)
#' 
#'  # with omission of single reference taxa  
#'    ddNoSing <- divDyn(corals, tax="genus", bin="stg", om="ref")
#'    lines(stages$mid, ddNoSing$divRT, lwd=2, col="red")
#'
#'  # using the estimated ages (less robust) - 10 million years
#'    # mean ages
#'    corals$me_ma <- apply(corals[, c("max_ma", "min_ma")], 1, mean)
#'    # ages reverse the direction of time! set ages to TRUE in this case
#'    ddRadio10 <- divDyn(corals, tax="genus", bin="me_ma", 
#'		breaks=seq(250,0,-10), ages=TRUE)
#'    lines(ddRadio10$bin, ddRadio10$divRT, lwd=2, col="green")
#'       
#'  # legend
#'    legend("topleft", legend=c("all", "no single-ref. taxa", "all, estimated ages"), 
#'      col=c("black", "red", "green"), lwd=c(2,2,2), bg="white")
#'    
#'
#' @export
divDyn <- function(dat, tax, bin, ages=FALSE, breaks=NULL, coll="collection_no", ref="reference_no", om=NULL,noNAStart=FALSE, data.frame=TRUE, filterNA=FALSE)
{
	
	# checking the binning argument
	# is numeric
	if(!is.numeric(dat[,bin])){
		stop("The bin variable is not numeric.")
	}
	
	if(ages){
		dat[, bin] <- -dat[,bin]
	}
	
	# what do you want to do with the bin numbers?
	# nothing, individual time slices
	if(is.null(breaks)){
		# smallest bin value
		minBin<-min(dat[,bin], na.rm=T)
			
		# if bin values are positive integers
		if(sum(dat[,bin]%%1, na.rm=T)==0 & minBin>0){
			
			binID<-NULL
			
		# non-integers: factorization
		}else{

			#### redo
			# potential bin names in order
			levs <- sort(unique(dat[,bin]))

			# their index in the processed script
			vals <- 1:length(levs)
			names(vals)<-levs

			# NA bin entries
			charBins <- as.character(dat[,bin])
			naBins <- is.na(charBins)
			
			# replace bin numbers with positive integers
			dat[!naBins,bin] <- vals[charBins[!naBins]]+4 # add some offset 
		
			# use later
			binID <-c(rep(NA, 4), levs)


		}
	
	# use a predefined binning of numeric values
	}else{	
		if(!is.numeric(breaks)) stop("The breaks argument has to be a numeric vector. ")
		
		# revert the breaks too, in case ages are provided
		if(ages) breaks <- -breaks
		# order the breaking vector 
		breaks<-sort(breaks)
		
		# calculate the bin averages
		both<-cbind(c(breaks,NA),c(NA, breaks))
		means<-apply(both, 1, mean)
		means<-means[!is.na(means)]
		
		# bin the variable
		fact<-cut(dat[,bin], breaks)
		
		# the function will use this to output the data
		newBin<-as.numeric(fact)
		
		if(length(unique(newBin))<4) stop("At least 4 time slices are necessary.")
		
		# replace bin numbers with positive integers
		dat[,bin] <- newBin
		
		# save the relevant means, use the binID later to identify the position
	
		to<-max(newBin,na.rm=T)
		binID<-c(rep(NA,4),means[1:to])
		
		 # add 4 as an offset
		 dat[,bin]<-dat[,bin]+4
		
	}
	
	
	# the omission phase
	if(!is.null(om)){
		dat<-dat[!omit(dat, tax=tax, bin=bin,om=om, ref=ref, coll=coll, filterNA=filterNA),]
	}
	
	# sub dataset
	subDat <- unique(dat[,c(tax, bin)])
	
	# omit NAs
	bNeed<- !(is.na(subDat[,tax]) | is.na(subDat[,bin]))

	# taxon vars
	taxVar<-subDat[bNeed, tax]
	binVar<-subDat[bNeed, bin]

	if(any(""==taxVar)){
		warning("Taxon name <\"\"> (empty quotes) is detected!")
	}

	#the maximum number of time slices
	nVectorLength<-max(binVar, na.rm=T)
	
	#starting interval
	nStart<-min(binVar, na.rm=T)
	
	#ending interval
	nEnd<-max(binVar, na.rm=T)
	
	#the vector of time slice numbers
	nTimeSlice<-c(rep(NA, nStart-1), nStart:nEnd)
	
	
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
		
	# cbind all together
	dCountsAndMetrics<-cbind(bin=nTimeSlice,  counts, metrics)
	
	if(!is.null(binID)){
		dCountsAndMetrics[,"bin"]<-binID
		
		# and get rid of the NAs
		dCountsAndMetrics<-dCountsAndMetrics[5:nrow(dCountsAndMetrics),]

		nStart <- 1
		nEnd <- nEnd-4
	}
	if(ages){
		dCountsAndMetrics[,"bin"]<- -dCountsAndMetrics[,"bin"]
	}
		
	#create the returning table
	dCountsAndMetrics<-dCountsAndMetrics[,c(
		"bin",
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
	)]
	
	if(data.frame){				
		dCountsAndMetrics<-as.data.frame(dCountsAndMetrics, stringsAsFactors=F)
	}
	
	
	#!!!nTot3tSampComp

	# coerce NAs in variables and intervals where they were not supposed to be
		# starting part of the table, only if there is no pseudofactorization involved, i.e. index conservation
		if(nStart>1 & is.null(binID)) dCountsAndMetrics[1:nStart-1,] <- NA
	
		# first row
		dCountsAndMetrics[nStart,c("t2d","t3","tPart","tGFd","tGFu", "tThrough", "tExt", "divBC", "O2f3")] <- NA
		dCountsAndMetrics[nStart+1,c("tGFd","oriGF", "ori2f3", "O2f3")] <- NA
	
		# last row
		dCountsAndMetrics[nEnd,c("t2u","t3","tPart","tGFd","tGFu", "ext2f3", "tOri", "tThrough", "E2f3")] <- NA
		dCountsAndMetrics[nEnd-1,c("tGFu","extGF", "ext2f3", "E2f3")] <- NA

	#want to see the NA's at the beginnning? (when the time series does not start with bin 1)
		
	# only use noNAStart when the simple indices are used	
	if (noNAStart==TRUE & is.null(binID))
	{
		dCountsAndMetrics<-dCountsAndMetrics[nStart:nEnd,]
		rownames(dCountsAndMetrics)<-NULL
	}
			
	
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
	metrics[,"divBC"] <- counts[,"tThrough"]+counts[,"tExt"]

	# total diversity
	sib. <- counts[,"divSIB"]-counts[,"tSing"]
	div <- metrics[,"divBC"]+counts[,"tOri"]                  
	metrics[,"divRT"]<-div+counts[,"tSing"]
	
	#Sampling parameters
	#sampling probability (Foote, 2000)
	obs <- sib.-counts[,"tExt"]-counts[,"tOri"]
	metrics[,"sampRange"] <- obs/counts[,"tThrough"] 

	#Three-timer sampling completeness (Alroy, 2008)	
	metrics[,"samp3t"] <- counts[,"t3"]/(counts[,"t3"]+counts[,"tPart"])
			
			#Sampling completeness of the entire time series
			nTot3tSampComp<- sum(counts[,"t3"], na.rm=T)/(sum(counts[,"t3"], na.rm=T)+sum(counts[,"tPart"], na.rm=T))
			
	# proportional extinctions
		metrics[,"extProp"]<-(counts[,"tExt"]+counts[,"tSing"])/metrics[,"divRT"]
		metrics[,"oriProp"]<-(counts[,"tOri"]+counts[,"tSing"])/metrics[,"divRT"]
	
	#Foote (2000) rates
		metrics[,"extPC"]<- -log(counts[,"tThrough"]/(counts[,"tThrough"]+counts[,"tExt"]))
		metrics[,"oriPC"]<- -log(counts[,"tThrough"]/(counts[,"tThrough"]+counts[,"tOri"]))
		
	#Three-timer rates by Alroy (2008)
		#uncorrected:
		metrics[,"ext3t"] <- log(counts[,"t2d"]/counts[,"t3"]) # two-timer/three-timer ratio (bottom)
		metrics[,"ori3t"] <- log(counts[,"t2u"]/counts[,"t3"]) # two-timer/three-timer ration (top)
		
		#corrected:
			#extinction rates:
				thtSampCompNext <- c(metrics[2:nrow(counts),"samp3t"],NA) # Sampling probability in subsequent bin (BIN5)
				metrics[,"extC3t"] <- metrics[,"ext3t"] + log(thtSampCompNext)
				#nC3tExt[nC3tExt<0] <- 0 # omit negative values
			
			#origination rates:
				thtSampCompPrev <- c(NA, metrics[1:(nrow(counts)-1),"samp3t"]) # Sampling probability in previous bin (BIN5)
				metrics[,"oriC3t"] <- metrics[,"ori3t"] + log(thtSampCompPrev)
				#nC3tOri[nC3tOri<0] <- 0 #omit negative values
	 
	
	#Gap filler rates by Alroy(2014)
		metrics[,"extGF"]<-log((counts[,"t2d"]+counts[,"tPart"])/(counts[,"t3"]+counts[,"tPart"]+counts[,"tGFu"]))
		metrics[,"oriGF"]<-log((counts[,"t2u"]+counts[,"tPart"])/(counts[,"t3"]+counts[,"tPart"]+counts[,"tGFd"]))
	
	# second- for third (Alroy, 2015)
	
	#	substituteF<-function(x){
	#		# in such cases ... 
	#		# reversed ordering
	#		first <- x[1]<x[2] & x[2]<x[3]
	#		
	#		# s2 is minimum
	#		second <- x[2]<x[1] & x[2]<x[3]
	#		
	#		# s2 is maximum
	#		third <- x[2]>x[1] & x[2]>x[3]
	#		
	#		# substituting the second lowest of the three counts for s3
	#		if(first | second | third){
	#			sort(x)[2]
	#		}else{
	#			x[3]
	#		}
	#	}
		# this is not enough, because the proportion can be still very negative,
		# he correctly says at the end of the paragraph, it should be a second in the order
	
		# extinction proportions (Eq. 4 in Alroy, 2015)
		sSubD<-apply(counts[,c("s1d","s2d","s3d")],1,function(x) sort(x)[2])
		metrics[,"E2f3"] <- (counts[,"s1d"]-sSubD)/(counts[,"t2d"]+counts[,"tPart"])
		
		# origination proportions
		sSubU<-apply(counts[,c("s1u","s2u","s3u")],1,function(x) sort(x)[2])
		metrics[,"O2f3"] <- (counts[,"s1u"]-sSubU)/(counts[,"t2u"]+counts[,"tPart"])
	
		# transform to classical rate form
		metrics[,"ext2f3"] <-log(1/(1-metrics[,"E2f3"]))
		metrics[,"ori2f3"] <-log(1/(1-metrics[,"O2f3"]))
		
	
	#corrected sampled-in-bin diversity
		metrics[,"divCSIB"]<-counts[,"divSIB"]*nTot3tSampComp/metrics[,"samp3t"]
	

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
#' @param dat \code{(data.frame)} Occurrence dataset, with \code{bin}, \code{tax} and \code{coll} as column names.
#' @param coll \code{(character)} The variable name of the collection identifiers. 
#' @param ref \code{(character)} The variable name of the reference identifiers (optional). 
#' @param om \code{(character)} The type of omission. \code{"coll"} omits occurrences of taxa that occurr only in one collection. \code{"ref"} omits occurrences of taxa that were described only in one reference. \code{"binref"} will omit the set of single reference taxa that were described by more than one references, but appear in only one reference in a time bin.
#' @param filterNA \code{(logical)} Additional entries can be added to influence the dataset that might not have reference or collection information (\code{NA} entries). These occurrences are treated as single-collection or single-reference taxa if the \code{na.rm} argument is set to \code{FALSE} (default). Setting this argument to \code{TRUE} will keep these entries. (see example)
#' 
#' @examples
#' # omit single-reference taxa
#'   data(corals)
#'   data(stages)
#'   toOmit <- omit(corals, bin="stg", tax="genus", om="ref")
#'   dat <- corals[!toOmit,]
#' 
#' # within divDyn
#'   # plotting
#'	  tsplot(stages, shading="series", boxes="per", xlim=c(260,0), 
#'	    ylab="range-through diversity (genera)", ylim=c(0,230))
#'   # multiple ref/slice required
#'   ddNoSing <- divDyn(corals, tax="genus", bin="stg", om="binref")
#'   lines(stages$mid, ddNoSing$divRT, lwd=2, col="red")
#'
#'   # with the recent included (NA reference value)
#'   ddNoSingRec <- divDyn(corals, tax="genus", bin="stg",
#'     om="binref", filterNA=TRUE)
#'   lines(stages$mid, ddNoSingRec$divRT, lwd=2, col="blue")
#'   
#'   # legend
#'   legend("topleft", legend=c("no single-ref. taxa", 
#'     "no single-ref. taxa,\n with recent"), 
#'     col=c("red", "blue"), lwd=c(2,2))
#' @export
omit <- function(dat, tax="genus", bin="bin", coll="collection_no", ref="reference_no", om="ref", filterNA=FALSE){

	if(!om%in%c("coll", "ref","binref")) stop("Invalid om argument.")
	
	if(om=="coll"){
		# omit multiple occ rows of same tax (genus) and same coll
		nonDupl <- !duplicated(dat[,c(tax, coll)])
		
		# which taxa come from just one collection?
		tabSing <- table(dat[nonDupl,tax])
		
		# single collection taxa
		singTax<- names(tabSing)[tabSing==1]
		
		# omit the single collection taxa
		boolEnd<-dat[,tax]%in%singTax
		
		# if na.rm TRUE than do not omit NA collection stuff
		if(filterNA){
			boolEnd <- boolEnd & !is.na(dat[,coll])
		}
	}
	
	if(om=="ref"){
		# omit multiple occ rows of same tax (genus) and same coll
		nonDupl <- !duplicated(dat[,c(tax, ref)])
		
		# which taxa come from just one collection?
		tabSing <- table(dat[nonDupl,tax])
		
		# single collection taxa
		singTax<- names(tabSing)[tabSing==1]
		
		# omit the single collection taxa
		boolEnd<-dat[,tax]%in%singTax
		
		# if na.rm TRUE than do not omit NA collection stuff
		if(filterNA){
			boolEnd <- boolEnd & !is.na(dat[,ref])
		}
	}
	
	
	
	if(om=="binref"){
		nonDupl <- !duplicated(dat[,c(tax, ref, bin)])
		
		
		activeDat <- dat[nonDupl, ]
		
		tap<-tapply(INDEX=activeDat[,bin], X=activeDat[,tax], function(x){
					
			tabSing <- table(x)
		
			nonSingTax<- names(tabSing)[tabSing!=1]
			return(nonSingTax)
			
		
		})
		# list of taxa that do not occur in in only reference/slice
		taxaMoreThanOne <- unique(unlist(tap))
		
		tap2<-(1:nrow(dat))[!dat[, tax]%in%taxaMoreThanOne]
		
	
		boolEnd<-rep(FALSE, nrow(dat))
		boolEnd[tap2]<-TRUE
		
		# if na.rm TRUE than do not omit NA reference stuff
		if(filterNA){
			boolEnd <- boolEnd & !is.na(dat[,ref])
		}
		
	}
	
	return(boolEnd)

}


.onUnload <- function (libpath) {
	library.dynam.unload("divDyn", libpath)
}
