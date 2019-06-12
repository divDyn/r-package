#' Subsampling trial functions
#' 
#' These functions create one subsampling trial dataset with a desired subsampling method
#' 
#' The essence of these functions are present within the subsampling wrapper function \code{\link{subsample}}. Each function implements a certain subsampling type.
#' The return value of the funcfions by default is a \code{logical} vector indicating which rows of the original dataset should be present in the subsample. 
#' The inexact method for SQS is implemented here as it is computationally less demanding. 
#'
#' \strong{References:}
#'
#' Alroy, J., Marshall, C. R., Bambach, R. K., Bezusko, K., Foote, M., Fürsich, F. T., … Webber, A. (2001). Effects of sampling standardization on estimates of Phanerozoic marine diversification. Proceedings of the National Academy of Science, 98(11), 6261-6266.
#' 
#' Alroy, J. (2010). The Shifting Balance of Diversity Among Major Marine Animal Groups. Science, 329, 1191-1194. https://doi.org/10.1126/science.1189910
#'
#' Raup, D. M. (1975). Taxonomic Diversity Estimation Using Rarefaction. Paleobiology, 1, 333-342. https: //doi.org/10.2307/2400135
#' 
#' @param q (\code{numeric)}: Subsampling level argument (mandatory). Depends on the subsampling function, it is the number of occurrences for \code{"cr"}, and the number of desired occurrences to the power of \code{xexp} for O^x^W. It is also the quorum of the SQS method.
#' 
#' @param x (\code{data.frame}): Occurrence dataset, with \code{bin}, \code{tax} and \code{coll} as column names.
#' 
#' @param bin (\code{character}): The name of the subsetting variable (has to be integer). For time series, this is the time-slice variable. Rows with \code{NA} entries in this column will be omitted.
#' @param useFailed (\code{logical}): If the bin does not reach the subsampling quota, should the bin be used? If \code{bin!=NULL} and \code{useFailed=TRUE} then only \code{TRUE} values will be output (indicating the use of the full dataset). 
#' @param keep (\code{numeric}): The bins, which will not be subsampled but will be added to the subsampling trials. NIf the number of occurrences does not reach the subsampling quota, by default it will not be represented in the subsampling trials. You can force their inclusion with the \code{keep} argument separetely (for all, see the \code{useFailed} argument). Only applicable when \code{bin!=NULL}. 
#' @param unit (\code{character}): Argument of the CR subsampling type. The name of the variable that designates the subsampling units. In every bin, CR selects a certain number (quota) of entries from the dataset. By default (\code{unit=NULL}), the units will be the rows, and the \code{q} number of rows will be selected in each bin.
#' However, this can be a higher level category that has multiple entries in the each bin. If \code{unit} is a valid column of the dataset \code{x}, then CR will select \code{q} number entries in this variable, and will return all the corresponding rows.
#' @param showFailed (\code{logical}): Toggles the output of the function. If set to \code{TRUE} the output will be a list, including both the default output (logical vector of rows) and the \code{numeric} vector of bins that did not have enough entries to reach the quota \code{q}. Only applicable when \code{bin!=NULL}. 
#' @rdname subtrial
#' @examples
#' #one classical rarefaction trial
#'   data(corals)
#' # return 5 references for each stage
#'   bRows<-subtrialCR(corals, bin="stg", unit="reference_no", q=5)
#'   # control
#'   unCor<-unique(corals[bRows,c("stg", "reference_no")])
#'   table(unCor$stg)
#'
#' @export
subtrialCR<-function(x, q, bin=NULL, unit=NULL, keep=NULL, useFailed=FALSE, showFailed=FALSE){
	quoVar <- q
	if(!is.null(bin)){
		# omit NA bin entries
		x <- x[!is.na(x[,bin]),]
	
		if(is.null(unit)){
			# order the rows 
			oBin <- order(x[,bin, drop=TRUE])
			
			# the only important variable
			binVar <- x[oBin,bin, drop=TRUE]
			
			# intact rows
			keepRows<-binVar%in%keep
			
			# output of cpp function
			binCRres<-.Call('_divDyn_CRbinwise', PACKAGE = 'divDyn', binVar, quoVar)
			
			# increase indexing (R)
			binCRres[,1]<-binCRres[,1, drop=TRUE]+1
			
			# rows where quota not reached
			failOutput <- binCRres[,1, drop=TRUE]==-8
			
			# bins where quota not reached
			fail<-binCRres[failOutput,2, drop=TRUE]
			
			# the rows where 
			usedRows<-rep(FALSE, nrow(x))
			
			# output by subsampling
			usedRows[binCRres[!failOutput,1, drop=TRUE]] <- TRUE
			
			usedRows[keepRows] <- TRUE
			# should be kept intact
			res<-list(rows=usedRows, fail=fail)
		
					
		}else{
			if(!unit%in%colnames(x)) stop("variable unit not found in x.")
			
			binVar<-x[,bin, drop=TRUE]
			unitVar <- x[,unit, drop=TRUE]
			
			rows <- 1:nrow(x)
			
			tap <- tapply(INDEX=binVar, X=rows, function(w){
				unitSlice<-unitVar[w]
				
				# random sample
				uniqueUnit<-sample(unique(unitVar[w]))
				
				if(quoVar<=length(uniqueUnit)){
					# rows corresponding t
					return(w[unitSlice%in%uniqueUnit[1:quoVar]])
				}else{
					return(NULL)
				}
			
			
			})
			
			# what was not enough for the subsampling
			failed<-unlist(lapply(tap, is.null))
			failed<-as.numeric(names(failed[failed]))
			failed<-failed[!failed%in%keep]
			
			# the rows that should be passed
			trialRows <- rows%in%unlist(tap)
			
			# these rows should be present regardless of the subsampling
			intactRows<-binVar%in%keep
			
			# combine the two
			subsampleRows<-trialRows | intactRows
			
			# result
			res<-list(rows=subsampleRows, fail=failed)
		
		}
		if(useFailed){
			res$rows[x[, bin, drop=TRUE]%in%res$fail] <- TRUE
		}
		if(showFailed){
			return(res)
		}else{
			return(res$rows)
		}
	}else{
		if(is.null(unit)){
			# use only the first column - like a vector
			if(is.data.frame(x)) x <- x[,1, drop=TRUE]

			# then write for a vector
			if(is.vector(x)){
				# quoVar too big?
				# if useFailed=FALSE, check whether the quota is low enoug
				if(!useFailed & length(x)<quoVar) stop("The entered quota is too low. Set 'useFailed=TRUE' to get a result with the unsubsampled dataset.")
				if(length(x)< quoVar) quoVar <- length(x)
		
				# do the actual trial - just some random rows
				index<-sample(1:length(x), quoVar)
				res <- rep(FALSE, length(x))
				res[index] <- TRUE		
				return(res)
			}

		}else{
			if(is.data.frame(x)) x <- x[,unit, drop=TRUE]
			if(is.vector(x)){
				tUn <- table(x)
				# if useFailed=FALSE, check whether the quota is low enoug
				if(!useFailed & length(tUn)<quoVar) stop("The entered quota is too low. Set 'useFailed=TRUE' to get a result with the unsubsampled dataset.")
				if(length(tUn)< quoVar) quoVar <- length(tUn)
				
				# randomize the unit table
				randTab <- sample(tUn)
				passed <- randTab[1:quoVar]

				# boolean index
				res <- as.character(unitVar)%in%names(passed)
				return(res)

			}
		}

	}
}


#' @param coll (\code{character}): The variable name of the collection identifiers. 
#' @param xexp (\code{numeric}): Argument of the OxW type. The exponent of by-list subsampling, by default it is 1.
#' @rdname subtrial
#' @export
subtrialOXW<-function(x, q, bin=NULL, coll=NULL,xexp=1, keep=NULL, useFailed=FALSE, showFailed=FALSE){
	quoVar <- q
	
	if(!is.null(bin)){
		# omit NA bin entries
		x <- x[!is.na(x[,bin, drop=TRUE]),]
	
		binVar <- x[,bin, drop=TRUE]
		collVar <- x[,coll, drop=TRUE]
		
		# future argument: appr represents how the subsampling quota is approached
		# appr can be "over", "under", "optimize"
		appr="over"
		
		# the chosen quota
		if(xexp==0){
			appr <- "other"
		}
		rows<- 1:length(binVar)
		
		#keep track of not enough
		ne<-NULL
		tap<-tapply(X=rows, INDEX=binVar, FUN=function(y){
		#	y<-rows[binVar==i]
			# because indices are subsetted
			collection<-collVar[y]
			
			# get the weight of each collection
			expTab<-table(collection)^xexp
			
			# shuffle
			expTab<-sample(expTab)
			
			# add up
			cumulative<-cumsum(expTab)
	
	
			# there are enough 
			if(length(cumulative)>0){
				if(cumulative[length(cumulative)]>quoVar){
					bSelect<-cumulative<=quoVar
					# potential forking!!! 
					if(appr=="over"){
						if(sum(bSelect)!=length(bSelect)){
							bSelect[sum(bSelect)+1] <- TRUE
						}
					}
					
					# the collections to keep
					keepCollections<-as.numeric(names(cumulative)[bSelect])
					return(y[collection%in%keepCollections])
					
				}else{
					return(NULL)
				}
			}else{
				return(NULL)
				
			}
		
		})
		
		# what was not enough for the subsampling
		failed<-unlist(lapply(tap, is.null))
		failed<-as.numeric(names(failed[failed]))
		failed<-failed[!failed%in%keep]
		
		# the rows that should be passed
		trialRows <- rows%in%unlist(tap)
		
		# these rows should be present regardless of the subsampling
		intactRows<-binVar%in%keep
		
		# combine the two
		subsampleRows<-trialRows | intactRows
		
		# result
		res<-list(rows=subsampleRows, fail=failed)
	
		if(useFailed){
			res$rows[x[, bin, drop=TRUE]%in%res$fail] <- TRUE
		}
		
		if(showFailed){
			return(res)
		}else{
			return(res$rows)
		}
	}else{
		# you only need the collection vector
		if(is.data.frame(x)){
			if(is.null(coll)) stop("If you provide a data.frame, you also have to provide a collection column.")
 			x <- x[,coll, drop=TRUE]
		}

		# then write for a vector
		if(is.vector(x)){

			if(any(is.na(x))) stop("The 'coll' variable includes NAs. ")

			# the table
			tColl <- table(x)
			tColl<- tColl^xexp

			# quoVar too big?
			# if useFailed=FALSE, check whether the quota is low enoug
			if(!useFailed & sum(tColl) < quoVar) stop("The entered quota is too high. Set 'useFailed=TRUE' to get a result with the unsubsampled dataset.")
			if(sum(tColl)< quoVar) quoVar <- sum(tColl)

			# randomly shuffle
			tCollRand <- sample(tColl)

			# cumulate them
			tCollCu<- cumsum(tCollRand)

			# select those collections that make the cut
			passed <- tCollCu<=quoVar

			# rows that belong to this set
			res<- as.character(x)%in%names(tCollCu)[passed]
		}
		return(res)
	}
}


#' @param tax (\code{character}): The name of the taxon variable.
#' @param ref (\code{character}): The name of the reference variable, optional - depending on the subsampling method.
#' @param excludeDominant \code{(logical)}: Argument of SQS. This parameter sets whether the dominant taxon should 
#'	be excluded from all calculations involving frequencies (this is the second correction of Alroy, 2010).
#' 
#' @param largestColl \code{(logical)}: Parameter of SQS. This parameter sets whether the occurrences of taxa only ever
#'  found in the most diverse collection should be excluded from the count of 
#'	single-publication occurrences. (this is the third correction of Alroy, 2010) Note that \code{largestColl=TRUE} is dependent on \code{excludeDominant=TRUE}. Setting \code{excludeDominant} to \code{FALSE} will turn this correction off.
#' 
#' @param fcorr \code{(character)}: Parameter for the inexact method of SQS. either "good" or "alroy". This argument changes the frequency correction procedure of the 
#'  'inexact' version of SQS (Alroy 2010). As not all taxa are present in the samples, 
#'  the sampled frequencies of taxa tend overestimate their frequencies in the sampling pool. 
#'  In Alroy (2010) these are corrected using Good's u ("good", default), in the later versions 
#'  of SQS this metric is changed to a different method using single occurrence and double occurrence taxa ("alroy"). 
#' @param byList (\code{character}): A parameter of the \code{"inexact"} method of SQS. Sets whether occurrences should be subsampled with (\code{FALSE}) or without (\code{TRUE}) breaking the collection integrity. 
#'	
#' @param appr (\code{character}): A parameter of the inexact method of SQS. Either "over" (default) or ("under"). The current 
#' version is not concerned with small fluctuations around the drawn subsampling quorum. 
#' Therefore, in the inexact algorithm, sampling is finished when the subset 
#' either is immediately below the quorum (\code{"under"}) or above it (\code{"over"}).
#' @param singleton \code{(character)}: A parameter of SQS. Either \code{"ref"}, \code{"occ"} or \code{FALSE}. If set to \code{"occ"}, the coverage estimator (e.g. Good's u) will be calculated based on the number of single-occurrence taxa. 
#' If set to "ref" the number of occurrences belonging to single-reference taxa will be used instead. In case of the inexact algorithm, if set to \code{FALSE} then coverage corrections of frequencies will not be applied.
#' @rdname subtrial
#' @export
subtrialSQS <- function(x, tax, q, bin=NULL, coll=NULL, ref=NULL, singleton="occ", excludeDominant=FALSE, 
	largestColl=FALSE, fcorr="good",byList=FALSE, keep=NULL, useFailed=FALSE, showFailed=FALSE, appr="under"){
	if(byList & is.vector(x)) stop("You cannot supply a vector for 'x' with by-list SQS.")
	if(byList & is.null(coll)) stop("You cannot do by-list subsampling without a collection variable!")
	if(byList & !is.null(coll)){
		if(!coll%in%colnames(x)) stop("Please provide a valid collection variable name.")
	} 
	# the quota variable
		quoVar <- q
		
	# is there a bin?
	if(!is.null(bin)){
		# omit NA bin entries
		x <- x[!is.na(x[,bin, drop=TRUE]),]
	
		# the frequencies
		freqVar <- frequencies(x=x,bin=bin, tax=tax, coll=coll, ref=ref, 
			singleton=singleton, excludeDominant=excludeDominant, largestColl=largestColl, fcorr=fcorr)
	
		# the bins
		binVar=x[,bin, drop=TRUE]
		
		# the collections (if present)
		if(is.null(coll)){
			collVar <- NULL
		}else{
			collVar=x[, coll, drop=TRUE]
		}
		
		res <- subsampleSQSinexact(binVar=binVar, freqVar=freqVar, quoVar=quoVar, collVar=collVar, byList=byList, intact=keep, appr=appr, trialRet="occ")
		
		if(useFailed){
			res$rows[x[, bin, drop=TRUE]%in%res$fail] <- TRUE
		}
		if(showFailed){
			return(res)
		}else{
			return(res$rows)
		}
	# no bins just a single trial
	}else{
		# vectors are possible but only when byList=FALSE
		if(is.vector(x)){
			if(byList | !is.null(coll)) stop("You cannot supply a vector for 'x' with by-list SQS.")
			x<-as.data.frame(x, stringsAsFactors=FALSE)
			colnames(x) <- "taxon"
			tax<-"taxon"
		}

		# frequency calculation , but without the bin=NULL
		freqVar<- frequenciesOne(x, tax=tax, coll=coll, ref=ref, singleton=singleton, 
			excludeDominant=excludeDominant, largestColl=largestColl, fcorr=fcorr)
		
		# by-list is by default: FALSE - to make coll necessary, you have to include it in the added arguments
		collVar <- NULL

		if(byList){
			# and reassing the collection variable
			collVar <- x[, coll, drop=TRUE]
		}

		# do the actual subsampling
		res <- subsampleSQSinexactOne(freqVar=quoVar, quoVar=quoVar, collVar=collVar, byList=byList, appr=appr, trialRet="occ")

		# if failed
		if(is.null(res)){
			if(useFailed){
				res <- rep(TRUE, nrow(x))
			}else{
				stop("The entered quota is too high. Set 'useFailed=TRUE' to get a result with the unsubsampled dataset.")
			}
		}

		return(res)

	}
}
