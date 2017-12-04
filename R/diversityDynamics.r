#' Cleanse Species Vector
#' 
#' This function will get rid of 
#'
#' This version will not keep subgenera, and will assign species to the base genus. The following qualifiers will be omitted:
#' "n.", "sp.", "?", "gen.", "aff.", "cf.", "ex", "gr.", "subgen.", "spp" and informal species designated with letters. Entries with "informal" and "indet." in them will also be invalidated. 
#' 
#' @param vec Character vector, the vector containing species names with qualifiers of open taxonomy.
#'
#' @param mode Character value, either "simple" or "debug". "simple" will return the cleaned species name vector, "debug" returns a data table that allows one by one checking.
#' @param collapse Character value, this argument will be passed to the paste function's argument of the same name. The character value to be inserted between the genus and species names.
#'
#' @examples
#' examp <- c("Genus cf. species", "Genus spp.", "Family indet.", "Mygenus yourspecies", "Okgenus ? questionsp")
#' spCleanse(examp) 
#' @export
# function to cleanse a noisy species name vector
spCleanse <- function(vec, mode="simple", collapse="_"){
	
	# keep the original
	vecOrig<-vec
	# presplit for parenthesis error
	vec<-gsub("\\( ", "\\(", vec)
	
	# presplit the the double space error
	vec<-gsub("  ", " ", vec)
	
	# split the string
	split<-strsplit(vec, " ")
	
	# excluded parts
		# these entries will be omitted without further issues	
		exclude <- c("n.", "sp.", "?", "gen.", "aff.", "cf.", "ex", "gr.", "subgen.", paste(LETTERS, ".", sep=""), "spp.")
		
		# if elements of this list are found in pairs, those will be excluded (but only if both are found, so species name "lato" can be left)
		jointExclude<-list(c("sensu", "lato"), c("sensu", "stricto"))
	
		# if these entries are round, the species name will be invalid
		special <- c("sp.1","sp.2", "informal", "indet.", letters)
	
	
	dual<-lapply(split, function(x){
	# missing entries
		if(sum(is.na(x))==length(x)){
			return(NA, NA)
		}
	
	#is a name starting with quotes - remove quotes
		
		quotes<-sapply(x, function(y){
			substr(y, 1,1)%in%c("\"")
		})
		if(sum(quotes)>0){
			tem<-unlist(strsplit(x[quotes], "\""))[2]
			x[quotes]<-tem
		}
		
	# is there a subgenus name - omit it
		# first character is parenthesis- omit
		parenth1<-sapply(x, function(y){
			substr(y, 1,1)=="("
		})
		
		if(sum(parenth1)>0){
			x<-x[!parenth1]
		}
		
		#last character is parenthesis- omit
		parenth2<-sapply(x, function(y){
			substr(y, nchar(y),nchar(y))==")"
		})
		
		if(sum(parenth2)>0){
			x<-x[!parenth2]
		}
		
	# omit the prefixes and suffixes
		x<-x[!x%in%exclude]
		
	# omit the jointly occurring notes (e.g. 'sensu lato')
		jointOcc<-unlist(lapply(jointExclude, function(y){
			sum(y%in%x)==length(y)
		
		}))
		if(sum(jointOcc)>0){
			je<-jointExclude[jointOcc]
			for(i in 1:length(je)){
				x<-x[!x%in%je[[i]]]
			}
			
		}
		
	# if there is a non-valid species name indicator - remove the entry
		if(sum(x%in%special)>0){
			return(c(NA, NA))
		}
		numConvert<-suppressWarnings(as.numeric(x))
		if(sum(!is.na(numConvert))>0){
			return(c(NA, NA))
		}
		
		if(length(x)==1){
			return(c(NA, NA))
		}
	
		
	# return genus and species name - potentially subspecies and crap at the end
		return(x)
	
	})
	
#	# not two
#	len<-unlist(lapply(dual, length))!=2
#	
#	
#	View(cbind(gen,sp)[len,])
#	
#	prob<-vec[len]
#	View(prob)
	
	# merge the genus and species in one column
	singular<-unlist(lapply(dual, function(x){
		if(is.na(x[1])){
			return(NA)
		}else{
			paste(x[1:2], collapse=collapse)
		}
	}))
	
	# if names start with " omit those as well
	# omit parentheses
	if(mode=="debug"){
		
		gen<-unlist(lapply(dual, function(x){
			x[1]
		}))
		
		sp<-unlist(lapply(dual, function(x){
			x[2]
		
		}))
	
		dat<-data.frame(original=vecOrig,genus=gen,species=sp, omitted=rep(FALSE, length(gen)))
		dat$omitted[is.na(singular)] <- TRUE
		dat$binomen <- singular
		return(dat)
	}
	if(mode=="simple"){
		return(singular)
	}
}


#' Time series from metrics of diversity dynamics 
#' 
#' This function calculates 
#' various metrics in the form of time series from resolved PBDB occurrences.
#' 
#' The following variables are produced:
#'
#' bin: the time slice number
#'
#' tThrough: the number of through-ranging taxa
#'
#' tOri: the number of originating taxa
#'
#' tExt: the number of taxa getting extinct
#'
#' t1: the number of stratigraphic singleton (single-interval) taxa
#'
#' t2d: the number of taxa that are present in the i-1th and the ith interval (lower two timers)
#'
#' t2u: the number of taxa that are present in the ith and the i+1th interval (upper two timers)
#'
#' tGFu: upper gap-fillers ($later)
#'
#' tGFd: lower gap-fillers ($later)
#'
#' t3: the number of three timer taxa, present in time slice i-1, i, and i+1
#'
#' tPart: the part timer taxa, present in time slice i-1,and i+1, but not in i
#'
#' extPC: the per capita extinction rates of foote (not normalized with interval length!)
#'
#' oriPC: the per capita origination rates of foote (not normalized with interval length!)
#'
#' ext3t: the three-timer extinction rates 
#'
#' ori3t: the three-timer origination rates 
#'
#' extC3t: the corrected three-timer extinction rates 
#'
#' oriC3t: the corrected three-timer origination rates
#'
#' divSIB: sampled-in-bin diversity (richness)
#'
#' divCSIB: corrected sampled-in-bin diversity (richness)
#'
#' divBC: boundary-crosser diversity (richness)
#'
#' divRT: range-through diversity (richness)
#'
#' sampRange: sampling probability (Foote)
#'
#' samp3t: three-timer sampling completeness (used as a correcting value in nC3tExt, nC3tOri and nCorrSIB)
#'
#' extGF: gap-filler extinction rates
#'
#' oriGF: gap-filler origination rates
#'
#' @param dat (data.frame): the data frame containing PBDB occurrences.
#' 
#' @param tax (char): variable  name of the occurring taxa (variable type: factor) - such as "occurrence.genus_name"
#' 
#' @param bin (char): variable name of the time slice numbers of the particular occurrences (variable type: int)- such as "slc" or whatever. Bin numbers should be in ascending order,can contain NA's, it can start from a number other than 1 and must not start with 0.
#' @param noNAStart (bool): useful when the dataset does not start from bin No. 1. Then noNAStart=TRUE will cut the first part of the resulting table, 
#' 						so the first row will contain the estimates for the lowest bin number.
#' 
#' @param Inf (logical): should infinites be converted to NAs?
#' 
#' @examples
#'	# import data
#'	  data(scleractinia)
#'	  data(stages)
#'
#'	# calculate metrics of diversity dynamics
#'    dd <- divDyn(scleractinia, tax="genus", bin="slc")
#'
#'	# plotting
#'	  plotTS(stages, shading="series", boxes="period", xlim=c(260,0), 
#'	  ylab="range-through diversity (genera)", ylim=c(0,230))
#'	  lines(stages$mid, dd$divRT, lwd=2)
#'
#'
#' @export
divDyn <- function(dat, tax="occurrence.genus_name", bin="bin", noNAStart=F, inf=F)
{
	#bin<-""
	#tax<-"occurrence.genus_name"
	
	#independent vector for the taxa names
	fTaxa<-dat[, tax]
	#independent vector for the time slice numbers
	nBins<-dat[, bin]
	
	#the maximum number of time slices
	nVectorLength<-max(nBins, na.rm=T)
	#starting interval
	nStart<-min(nBins, na.rm=T)
	
	#ending interval
	nEnd<-max(nBins, na.rm=T)
	
	#the vector of time slice numbers
	nTimeSlice<-c(rep(NA, nStart-1), nStart:nEnd)
	
	
	# Modified to two-timer metrics backward
	# Data vectors - definitions
	ex <- numeric(nVectorLength) #vector of extinct taxa
	or <- numeric(nVectorLength) #vector of originating taxa
	su <- numeric(nVectorLength) #vector of total survivors
	nD2tTaxa <- numeric(nVectorLength) #vector of two timers (bottom)
	nU2tTaxa <- numeric(nVectorLength) #vector of two timers (top)
	n3tTaxa <- numeric(nVectorLength) #vector of three timers 
	nPtTaxa <- numeric(nVectorLength) #vector for part timers
	n1tTaxa <- numeric(nVectorLength) #vector for singletons
	nThrough <- numeric(nVectorLength) #vector for through-ranging taxa 
	nSIB <- numeric(nVectorLength) #vector for SIB diversities
	nUGfTaxa<-numeric(nVectorLength) #vector of gap fillers (top for extinctions)
	nDGfTaxa<-numeric(nVectorLength) #vector of gap fillers (bottom, for originations)

	
	# for every time slice
	for (i in nStart:nEnd) 
	{
		rfl<-levels(factor(fTaxa[nBins==i])) #All genera sampled in bin
		rfupl<-levels(factor(fTaxa[nBins>i])) # All genera in younger bins
		rfdol<-levels(factor(fTaxa[nBins<i])) # All genera in older bins
		rfotl<-levels(factor(fTaxa[nBins!=i])) # All genera outside bin
		rf1d<-levels(factor(fTaxa[nBins==(i-1)])) # All genera in previous bin
		rf1u<-levels(factor(fTaxa[nBins==(i+1)])) # All genera in subsequent bin
		rf2d<-levels(factor(fTaxa[nBins==(i-2)])) # All genera in bin i-2
		rf2u<-levels(factor(fTaxa[nBins==(i+2)])) # All genera in bin i+2
		
	
	
		s1 <- length(rfl [rfl %in% rfupl]) #Survivors from bin, top crosser
		s2 <- length(rfl [rfl %in% rfdol]) #Survivors to bin, bottom crosser
		ns <- length(rfl [rfl %in% rfotl]) #non-single interval taxa
		nD2tTaxa[i] <- length(rfl [rfl %in% rf1d]) # Two-timers at bottom
		nU2tTaxa[i]<- length(rfl [rfl %in% rf1u]) # Two-timers at top
		n3tTaxa[i] <- length(rfl [rfl %in% rf1d & rfl %in% rf1u]) # Three-timers
		nPtTaxa[i] <- length(rf1d [rf1d %in% rf1u])-n3tTaxa[i] # Part timers
		nUGfTaxa[i] <-length(rf2u[rf2u %in% rf1d[!(rf1d %in% rf1u)]]) # gap fillers at the top (not present in i+1, but present in i-1 and i+2)
		nDGfTaxa[i] <-length(rf2d[rf2d %in% rf1u[!(rf1u %in% rf1d)]]) # gap fillers at the bottom (not present in i-1, but present in i+1 and i-2)
		nThrough[i] <- length(rfdol [rfdol %in% rfupl]) # Through-rangers


		n1tTaxa[i]<-length(rfl)-ns # single interval taxa
		nSIB[i]<-length(rfl) # sampled in bin taxa
		ex[i]<- length(rfl)-s1 # extinct taxa, sib-s1
		or[i]<- length(rfl)-s2 # originating taxa, sib-s2
		
		
	}
	
	
## when the lowest number of time slice is not one!
#	ttd<-ttd[start:end]
#	ttu<-ttu[start:end]
#	tht<-tht[start:end]
#	ptm<-ptm[start:end]
#	nThrough<-nThrough[start:end]
#	sg<-sg[start:end]
#	sib<-sib[start:end]
#	ex<-ex[start:end]
#	or<-or[start:end]
  
	
  
	nOri <- or-n1tTaxa					#top crosser, originating
	nExt <- ex-n1tTaxa					#bottom crosser, extinguishing
	nBoundaryCrosser <- nThrough+nExt                   #Boundary crossers
	div <- nBoundaryCrosser+nOri                   # total diversity
	sib. <- nSIB-n1tTaxa
	nRangeThrough<-div+n1tTaxa
	
	sib.[sib.==0]<-NA

	#Sampling parameters
		#sampling probability (Foote, 2000)
		ss <- nThrough
		obs <- sib.-nExt-nOri
		nobs <- ss-obs
		nFooteSampProb <- obs/ss 

		#Three-timer sampling completeness (Alroy, 2008)	
		n3tSampComp <- n3tTaxa/(n3tTaxa+nPtTaxa)
			
			#Sampling completeness of the entire time series
			nTot3tSampComp<- sum(n3tTaxa, na.rm=T)/(sum(n3tTaxa, na.rm=T)+sum(nPtTaxa, na.rm=T))
			
	
	#Foote (2000) rates
		nFooteExt<- -log(nThrough/(nThrough+nExt))
		nFooteOri<- -log(nThrough/(nThrough+nOri))
		
	#Three-timer rates by Alroy (2008)
		#uncorrected:
		n3tExt <- log(nD2tTaxa/n3tTaxa) # two-timer/three-timer ratio (bottom)
		n3tOri <- log(nU2tTaxa/n3tTaxa) # two-timer/three-timer ration (top)
		
		#corrected:
			#extinction rates:
				thtSampCompNext <- c(n3tSampComp[2:nVectorLength],NA) # Sampling probability in subsequent bin (BIN5)
				nC3tExt <- n3tExt + log(thtSampCompNext)
				#nC3tExt[nC3tExt<0] <- 0 # omit negative values
			
			#origination rates:
				thtSampCompPrev <- c(NA, n3tSampComp[1:nVectorLength-1]) # Sampling probability in previous bin (BIN5)
				nC3tOri <- n3tOri + log(thtSampCompPrev)
				#nC3tOri[nC3tOri<0] <- 0 #omit negative values
	 
	
	#Gap filler rates by Alroy(2013)
		nGfExt<-log((nD2tTaxa+nPtTaxa)/(n3tTaxa+nPtTaxa+nUGfTaxa))
		nGfOri<-log((nU2tTaxa+nPtTaxa)/(n3tTaxa+nPtTaxa+nDGfTaxa))
	
	#corrected sampled-in-bin diversity
		nCorrSIB<-nSIB*nTot3tSampComp/n3tSampComp
	
		
	#create the returning table
	dRatesAndMetrics<-cbind(bin=nTimeSlice, tThrough=nThrough, tOri=nOri, tExt=nExt, t1=n1tTaxa, t2d=nD2tTaxa, t2u=nU2tTaxa, tGFu=nUGfTaxa, tGFd=nDGfTaxa, t3=n3tTaxa, tPart=nPtTaxa, extPC=nFooteExt, oriPC=nFooteOri, ext3t=n3tExt, ori3t=n3tOri, extC3t=nC3tExt, oriC3t=nC3tOri, divSIB=nSIB, 
							divCSIB=nCorrSIB, divBC=nBoundaryCrosser,divRT=nRangeThrough, sampRange=nFooteSampProb, samp3t=n3tSampComp, extGF=nGfExt, oriGF=nGfOri)
	if(Inf!=T)
	{
		dRatesAndMetrics[is.infinite(dRatesAndMetrics)]<-NA
		
	}
							
	dRatesAndMetrics<-as.data.frame(dRatesAndMetrics, stringsAsFactors=F)
	
	
				#!!!nTot3tSampComp

	#want to see the NA's at the beginnning? (when the time series does not start with bin 1)
		if (missing(noNAStart)) {}
		else
		{
			if (noNAStart==TRUE)
			{
				dRatesAndMetrics<-dRatesAndMetrics[nStart:nEnd,]
			}
			
			if (noNAStart!=TRUE & noNAStart!=FALSE)
			{
				print("You have entered an invalid argument for noNAStart, no cropping will occurr in the final table.")
			}
		}	
	
	
	#return the table					
	return(dRatesAndMetrics)
}
