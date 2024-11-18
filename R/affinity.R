#' Environmental affinities of taxa
#'
#' This function will return the preferred environment of the taxa, given the distribution of occurrences.
#'
#' Sampling patterns have an overprinting effect on the frequency of taxon occurrences in different environments. The environmental affinity (Foote, 2006; Kiessling and Aberhan, 2007; Kiessling and Kocsis, 2015) expresses whether the taxa are more likely to occur in an environment, given the sampling patterns of the dataset at hand. The function returns the likely preferred environment for each taxon as a vector. \code{NA} outputs indicate that the environmental affinity is equivocal based on the selected method.
#'
#' \strong{The following methods are implemented:}
#'
#' \code{'majority'}: Environmental affinity will be assigned based on the number of occurrences of the taxon in the different environments, without taking sampling of the entire dataset into account. If the taxon has more occurrences in \emph{environment 1}, the function will return \emph{environment 1} as the preferred habitat. 
#'
#' \code{'binom'}: The proportion of occurrences of a taxon in \emph{environment 1} and \emph{environment 2} will be compared to a null model, which is based on the distribution of all occurrences from the stratigraphic range of the taxon (in \code{x} or if provided, in \code{reldat}). Then a binomial test is run on with the numbers of the most likely preference (against all else). The \code{alpha} value indicates the significance of the binomial tests, setting \code{alpha} to \code{1} will effectively switch the testing off: if the ratio of occurrences for the taxon is different from the ratio observed in the dataset, an affinity will be assigned. This is the default method. If an environment is not sampled at all in the dataset to which the taxon's occurrences are compared to, the binomial method returns \code{NA} for the taxon's affinity. 
#' 
#' \strong{References}
#'
#' Foote, M. (2006). Substrate affinity and diversity dynamics of Paleozoic marine animals. Paleobiology, 32(3), 345-366.
#'
#' Kiessling, W., & Aberhan, M. (2007). Environmental determinants of marine benthic biodiversity dynamics through Triassic–Jurassic time. Paleobiology, 33(3), 414-434.
#'
#' Kiessling, W., & Kocsis, Á. T. (2015). Biodiversity dynamics and environmental occupancy of fossil azooxanthellate and zooxanthellate scleractinian corals. Paleobiology, 41(3), 402-414.
#' 
#' @param x \code{(data.frame)} The occurrence dataset containing the taxa with unknown environmental affinities.
#' @param env \code{(character)} The environmental variable of the occurrences.
#' @param method \code{(character)} The method used for affinity calculations. Can be either \code{"binom"} or \code{"majority"}.
#' @param tax \code{(character)} The column name of taxon names.
#' @param bin \code{(character)} The column name of bin names.
#' @param coll \code{(character)} The column name of collection identifiers (optional). If this is provided, then then the multiple entries of a taxon within the collections will be treated as 1.
#' @param alpha \code{(numeric)} The alpha value of the binomial tests. By default binomial testing is off (\code{alpha=1}) and the methods returns that environment as the preferred one, which has the highest likelihood (odds ratio). 
#' @param reldat \code{(data.frame)} Database with the same structure as \code{x}.  \code{x} is typically a subset of \code{reldat}. If given, the occurrence distribution of \code{reldat} is used 
#' as the null model of sampling. Defaults to \code{NULL}, which means that \code{x} itself will be used as \code{reldat}.
#' @param na.rm \code{(logical)} Should the \code{NA} entries in the relevant columns of \code{x} be omitted automatically?
#' @param bycoll \code{(logical)} If set to \code{TRUE}, the number of collections (or samples, in \code{coll}) will be used rather than the number of occurrences.
#'
#' @examples
#'	data(corals)
#'	# omit values where no occurrence environment entry is present, or where unknown
#'	  fossils<-subset(corals, stg!=95)
#'	  fossilEnv<-subset(fossils, bath!="uk")
#'	# calculate affinities
#'	  aff<-affinity(fossilEnv, env="bath", tax="genus", bin="stg", alpha=1)
#'	
#' @export
#' @return A named vector, values corresponding to affinities.
affinity<-function(x, tax, bin, env, coll=NULL, method="binom", alpha=1,reldat=NULL, na.rm=FALSE, bycoll=FALSE){
# version 2.0
#	x <- fossilEnv
#	env <- "bath"
#	tax <- "genus"
#	bin <- "stg"
#	method <- "binom"
#	alpha <- 1
#	reldat <- NULL
#	coll <- NULL
#	na.rm<-TRUE

	if(method=="majority" & !is.null(reldat)) warning("Majority rule selected, reldat will be ignored.")
	
	if(bycoll) if(any(!coll%in%colnames(x))) stop("The \'coll\' argument has to be a column name of \'x\' if \'bycoll=TRUE\'.")

	# omit everything from x that is not necessary
		x<-unique(x[,c(coll, tax,bin, env)]) # this can be faster!!
	
	match.arg(method, c("binom", "majority"))
	
	# omit NA bins
	naBin <- is.na(x[,bin, drop=TRUE])
	naTax <- is.na(x[,tax, drop=TRUE])
	naEnv <- is.na(x[,env, drop=TRUE])
	if(!is.null(coll)) naColl <- is.na(x[,coll, drop=TRUE])

	if(!na.rm){
		# stop execution
		if(any(naBin)) stop("The \'bin\' variable contains NAs. Omit these or set \'na.rm=TRUE\'")
		if(any(naTax)) stop("The \'tax\' variable contains NAs. Omit these or set \'na.rm=TRUE\'")
		if(any(naEnv)) stop("The \'env\' variable contains NAs. Omit these or set \'na.rm=TRUE\'")
		if(!is.null(coll)) if(any(naColl)) stop("The \'coll\' variable contains NAs. Omit these or set \'na.rm=TRUE\'")
	}else{
		# go forward
		if(!is.null(coll)){
			x <- x[!naBin & !naTax & !naEnv & !naColl, ]
		}else{
			x <- x[!naBin & !naTax & !naEnv, ]
		}
	}
	
	# the affinity variable
		affLevels<-levels(factor(x[,env, drop=TRUE]))
	
	# create an FAD-LAD matrix first
		dFL<-fadlad(x, tax, bin)
	
	# by this time fadLad will probably give a warning, but process this nevertheless
		if(any(""==x[,tax, drop=TRUE])){
			rownames(dFL)[rownames(dFL)==""] <- "emptyQuotes"
			x[x[,tax]=="",tax] <- "emptyQuotes"
		}
	
	# add the names to the matrix so apply can process it
		dFL$taxon<-rownames(dFL)

	# the bins
		allBins<-sort(unique(x[,bin, drop=TRUE]))
		
	# make a 3D array - bin, taxon, environment
		occArr <- array(0, dim=c(length(allBins), nrow(dFL), length(affLevels)))
		dimnames(occArr) <- list(allBins, rownames(dFL), affLevels)

	# NOTE: nothing needs to be done if bycoll=TRUE - if coll is provided, every occurrence will mean one collection, anyway
	# fill it with values
	for(i in 1:length(affLevels)){
		# environment-specific subset
		firstDat<-x[x[,env, drop=TRUE]==affLevels[i],]
	
		# tabulate
		fTab<-table(firstDat[,bin, drop=TRUE], firstDat[,tax, drop=TRUE])
		class(fTab)<-"matrix"

		# add to the array
		occArr[rownames(fTab), colnames(fTab),i]<-fTab
	}

	# reference/relative dataset
		# as x is not used from now on, use that to save time and memory
		# in case a relative dataset is added
		if(!is.null(reldat)){
			x <- reldat
		}
		
		# by-collection
		if(bycoll){
			# omit occurrence-level data, use only collections
			x <- unique(x[,c(bin, env, coll)])
		}

		# tabulate
		relTab <- table(x[,bin, drop=TRUE], x[,env, drop=TRUE])
		class(relTab) <-"matrix"
		
		# check if reldat was a different dataset
		if(!is.null(reldat)){
			# POTENTIAL FORKING POINT!
			# select only those columns that are present in the analyzed dataset
			if(any(!affLevels%in%colnames(relTab))) stop("The provided \'reldat\' does not contain all \'env\' entries of \'x\'")
			if(any(!allBins%in%rownames(relTab))) stop("The provided \'reldat\' does not contain all \'bin\' entries of \'x\'")
		}

		# pass only the relevant entries
		relTab<-relTab[as.character(allBins),affLevels]

		
	# calculate the affinity of every taxon
	affVarTaxon<-apply(dFL, 1, FUN=function(w){

#	affVarTaxon<-character(length(dFL$taxon))
#	for (i in 1:length(affVarTaxon)){
#		w<- dFL[i,]

	#subset of the taxons ranges
		# bins where the taxon is present
		vectRange<-as.numeric(w[1]):as.numeric(w[2])

		################ emergency bug fix to cover the gaps
		# subset it to those that are present in the rownames
		vectRange <- vectRange[vectRange%in%rownames(relTab)]
		################
		
		# what you relate to
			thisRel <- relTab[as.character(vectRange), , drop=FALSE]
			
			# the occurrence probabilities
			relProbs<-apply(thisRel, 2, sum)/sum(thisRel)
			
		# taxon occurrences in the different evnironment	
			taxOccBin <- occArr[as.character(vectRange), as.character(w[length(w)]),, drop=FALSE]
			taxOcc <- apply(taxOccBin, c(2,3), sum)

		# all occurrences of the taxon
			all <- sum(taxOcc)

		#no occurrences from known habitats
		if (sum(taxOcc)==0)		{
			return(NA)
#			affVarTaxon[i] <- NA
		}
		
		# otherwise, continue
		# majority rule
		if(method=="majority"){
			# highest number of occurrences
 			maxVal<-max(taxOcc, na.rm=T)

 			# which
			mostOcc <- which(maxVal==taxOcc)

			# if undecided, return NA
			if(length(mostOcc)>1){
				return(NA)
#				affVarTaxon[i] <- NA
			}else{
				return(affLevels[mostOcc])
			}
		} # end majority
		
		# binomial basic
		if(method=="binom"){

			# an important imperative: if the environment is not sampled, then
			# what do we know about the taxon's affinity?
			if(any(relProbs==0)) return(NA)

			# 1. which is the most likely (odds ratios) 
			# taxon probabilities
			taxProbs <- taxOcc/all

			# the odds ratios
			OR <-taxProbs/relProbs

			# which environment has the highest odds ratio?
			# replace infinties with NAs
			OR[!is.finite(OR)] <- NA

			# which is the highest?
			maxVal <- max(OR, na.rm=T)
			highestOR <- which(OR==maxVal)

			# in case it is undecided
			if(length(highestOR)>1){
				return(NA)
#				affVarTaxon[i] <- NA
			# if there is one largest
			# do a binomial test for this in constrast with the others
			}else{
				focus <-taxOcc[highestOR] 
				pFocus <- relProbs[highestOR]

				binTest <- stats::binom.test(focus, all, pFocus, "greater")$p.val

				if(!is.null(alpha)){
					if(binTest<=alpha)
					{
						return(affLevels[highestOR])
	#					affVarTaxon[i] <- affLevels[highestOR]
					}else{
						return(NA)
	#					affVarTaxon[i] <- NA
					}
				}else{
					names(binTest) <-  affLevels[highestOR]
					return(binTest)
				}
			}
			
		} # end binom
	}

	)
	
	return(affVarTaxon)
}


