#' Origination/extinction response table for statistical modelling.
#' 
#' This function takes an occurrence dataset and reformats it to a table that can be used as input for logistic models.
#' 
#' Every entry in the output table corresponds to one cell in the \code{bin}/\code{tax} matrix. This function omits duplicates and concatenates two \code{logical} vectors (response variables) to the occurrence dataset:  
#' The \code{ori} vector is \code{TRUE} in the interval when the taxon first appeared, and \code{FALSE} in all others. The \code{ext} vector is \code{TRUE} in the interval the taxon appeared for the last time, and \code{FALSE} in the rest.
#'
#' The true date of extinction and origination is unknown, therefore these events can only be expressed as probabilities. The argument \code{probs} allows the replacement of a binary response with two probability values, which are based on the apparent sampling patterns. For extinctions, when \code{probs} is set to \code{"samp3t"}, the response parameter for extinctions in the last bin of appearance is set to the three-timer sampling compelteness of the following bin. Assuming that the taxon'as range offset is not larger than a whole bin, if the taxon did not go extinct in the bin in which it appeared the last time, it is assumed to be going extinct in the following bin, and the remainder (1 -  sampling completeness) is assigned to that bin. The pattern is reversed for originations. For \code{probs="sampRange"}, the range-based completeness measures are applied in a similar fashion. For Phanerozoic-scale analyses, a whole bin difference between apparent event and the actual event is reasonable. See more in Reddin et al. 2021. Note that the response probabilities are set to missing values (\code{NA}s) when the probabilities cannot be calculated. The variable \code{ext} is also set to \code{NaN} for the early virtual extension of the range, and \code{ori} is treated the same for the late-extension. 
#'
#' \strong{References:}
#'
#' Reddin, C. J., Kocsis, Á. T., Aberhan, M., & Kiessling, W. (2021). Victims of ancient hyperthermal events herald the fates of marine clades and traits under global warming. Global Change Biology, 27(4), 868–878. https://doi.org/10.1111/gcb.15434
#' 
#' @param x \code{(data.frame)} Fossil occurrence data.frame.
#' 
#' @param tax \code{(character)} Variable name of the occurring taxa (variable type: \code{factor} or \code{character} - such as \code{"genus"})
#' 
#' @param bin \code{(character)} Variable name of the bin numbers of the occurrences. This variable should be \code{numeric} and should increase as time passes by (use negative values for age estimates). The current version only supports discreet, non-negative integer interval numbers.
#' 
#' @param taxvars \code{(character)} Taxon-specific column names of the variables that should be included in the output table. Only one entry/taxon is used, make sure that the data are clean.
#' 
#' @param rt \code{(logical)} Should the range-through assumption be applied within the function? If set to \code{TRUE} then missing occurrences will be interpolated with \code{FALSE} values in both the \code{ext} and \code{ori} variables. .
#' 
#' @param singletons \code{(logical)} Should single-interval taxa be included from the final table? This is not recommended, as it is impossible to get a \code{FALSE} response for these taxa. 
#' @param probs \code{(logical)} When set to \code{NULL}, the response variable will be binary. When set to \code{"samp3t"} or \code{"sampRange"} the response results will be probabilities, based on the three-timer sampling completeness, or the range-based sampling completeness, respectively. 
#'
#' 
#' @rdname modeltab
#' @examples
#' # load necessary data
#' data(corals)
#' # simple table
#' modTab<-modeltab(corals, bin="stg", tax="genus", taxvars=c("ecology", "family"))
#' # probabilities for extinction modeling
#' modTab2 <- modeltab(corals, bin="stg", tax="genus", probs="samp3t")
#' # only extinction response (omit virtual origination extensions)
#' extTab <- modTab2[!is.nan(modTab2$ext), ]
#' # only extinction response (omit virtual extinction extensions)
#' oriTab <- modTab2[!is.nan(modTab2$ori), ]
#'
#' @export
#' @return A data.frame with binary response variables.
modeltab<-function(x, tax, bin, taxvars=NULL,  rt=FALSE, singletons=FALSE, probs=NULL){
#	x <- corals[corals$stg!=95,]
#	bin<- "stg"
#	tax <- "genus"
#	taxvars <- "ecology"

	# sub dataset
	# the cell-variable
	unVar<- paste(x[, tax, drop=TRUE], x[, bin, drop=TRUE], sep="")
	subDat <- x[,c(tax, bin, taxvars)]
	subDat <- subDat[!duplicated(unVar),]

	# omit NAs
	bNeed<- !(is.na(subDat[,tax, drop=TRUE]) | is.na(subDat[,bin, drop=TRUE]))
	# taxon vars
	subDat<- subDat[bNeed,]

	# the data are complied - reorder
	subDat <- subDat[order(subDat[,tax, drop=TRUE], subDat[, bin, drop=TRUE]),]
	

	# filter singleton taxa
	if(!singletons){
		taxInd <- as.numeric(factor(subDat[, tax, drop=TRUE]))

		#if two 1 follow each other that indicates that a singleton is present
		diffInd <- c(1,diff(taxInd))
		diffShift<-c(diffInd[2:length(diffInd)], 1)

		# subset
		subDat<- subDat[!(diffInd & diffShift),]

	}

	# for extinctions
	rows <- 1:nrow(subDat)

	# for every taxon
	respInd<-tapply(INDEX=subDat[, tax], X=rows, function(w){
		# sampled bins
		bins<- subDat[w, bin, drop=TRUE]
		# return the FAD,LAD row numbers
		c(w[bins==min(bins)], w[bins==max (bins)])
	})

	# FAD-LAD row numbers 
	numResp <- matrix(unlist(respInd), ncol=2, byrow=T)

	# construct response variables
	ori <- rep(FALSE, length(rows))
	ext <- rep(FALSE, length(rows))
	ori[numResp[,1]] <- TRUE
	ext[numResp[,2]] <- TRUE
		

	# the model table to be returned
	datMod <- cbind( ori,ext, subDat)

	# if range interpolation is allowed
	if(rt){
		occurrence <- rep(TRUE, length(ori))
		datMod <- cbind(occurrence, datMod)

		# the rows of original data, gaps are NAs
		indInNew <- tapply(subDat[, tax], X=rows, function(w){
			bins <- subDat[w, bin, drop=TRUE]
			ran <- (range(bins)[1]:range(bins)[2])

			one <- rep(NA, length(ran))

			one[ran%in%bins] <- w
			return(one)
		})

		# in a single vector
		indInNew <- unlist(indInNew)

		# use this to subset the sampled part
		datMod<-datMod[indInNew,]
		datMod$occurrence[is.na(datMod$occurrence)] <- FALSE
		# fill in the obvious gaps
		# taxon name
		datMod[, tax] <- fill(datMod[, tax, drop=TRUE])
		# bin - increment by one
		datMod[, bin] <- fill(datMod[, bin, drop=TRUE], inc=1)

		# extinction/origination response
		datMod[, "ext"] <- fill(datMod[, "ext", drop=TRUE])
		datMod[, "ori"] <- fill(datMod[, "ori", drop=TRUE], forward=FALSE)

		if(!is.null(taxvars)){
			# the taxon variables
			for(i in 1:length(taxvars)){
				datMod[, taxvars[i]]<- fill(datMod[, taxvars[i]])
			}
		}
	}

	# final object
	rownames(datMod) <- NULL

	# replacement with probabilities?
	if(!is.null(probs)){
		dd<- divDyn(x, tax, bin)
		if(probs=="samp3t"){
			samp <- dd$samp3t
		}
		if(probs=="sampRange"){
			samp <- dd$sampRange
		}

		# new response
		datMod$newOri <- datMod$ori
		datMod$newExt <- datMod$ext

		# add occurrence column to separate from virtual entries
		if(!rt) datMod$occurrence <- TRUE

		# the direct completeness
		datMod$newOri[which(datMod$ori)] <- samp[datMod$stg[datMod$ori]-1]
		datMod$newExt[which(datMod$ext)] <- samp[datMod$stg[datMod$ext]+1]

		# the remainders
		# origination
		preOri <- datMod[datMod$ori,]
		preOri[,bin] <- preOri[,bin]-1
		preOri$newOri <- 1-preOri$newOri
		preOri$newExt<- NaN

		preOri$occurrence <- FALSE
		
		# extinction
		postExt<- datMod[datMod$ext,]
		postExt[,bin] <- postExt[,bin]+1
		postExt$newExt <- 1-postExt$newExt
		postExt$newOri<- NaN

		postExt$occurrence <- FALSE
		

		# bind it to the rest
		datMod <- rbind(preOri, datMod, postExt)
		datMod<- datMod[order(datMod[,tax], datMod[,bin]),]

		# replace old with new
		datMod$ext <- datMod$newExt
		datMod$ori <- datMod$newOri

		# get rid of
		datMod$newExt <- NULL
		datMod$newOri <- NULL
		# limit scope of data
		datMod <- datMod[datMod[,bin] >= min(x[,bin], na.rm=T) & datMod[,bin] <= max(x[,bin], na.rm=T) ,]
		rownames(datMod) <- NULL
		
	}
	
	return(datMod)


}

# It is very easy to append additional variables to this table. Such as the output of affinity()
