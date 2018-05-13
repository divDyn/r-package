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
#' tSing: the number of stratigraphic singleton (single-interval) taxa
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
#' extProp: proportional extinctions, including single-interval taxa
#'
#' oriProp: proportional originations, including single-interval taxa
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
#' extGF: gap-filler extinction rates (Alroy, 2014)
#'
#' oriGF: gap-filler origination rates (Alroy, 2014)
#'
#' E2f3: second-for-third extinction propotions (Alroy, 2015)
#'
#' O2f3: second-for-third origination propotions (Alroy, 2015)
#
#' ext2f3: second-for-third extinction rates (based on Alroy, 2015)
#'
#' ori2f3: second-for-third origination rates (based on Alroy, 2015)
#' 
#' References:
#' Foote, M. 2000.
#' Alroy, J. 2008
#' Alroy, J. 2014.
#' Alroy, J. 2015.
#'
#' @param dat (data.frame): the data frame of fossil occurrences.
#' 
#' @param tax (character): variable  name of the occurring taxa (variable type: factor) - such as "occurrence.genus_name"
#' 
#' @param bin (character): variable name of the time slice numbers of the particular occurrences (variable type: int)- such as "slc" or whatever. Bin numbers should be in ascending order,can contain NA's, it can start from a number other than 1 and must not start with 0.
#' @param noNAStart (logical): useful when the dataset does not start from bin No. 1. Then noNAStart=TRUE will cut the first part of the resulting table, 
#' 						so the first row will contain the estimates for the lowest bin number.
#' 
#' @param inf (logical): should infinites be converted to NAs?
#' @param data.frame (logical): should be a data frame or a matrix?
#' 
#' @param om (character): the om parameter of the omit() function. If set to NULL (default), then no occurrences will be omitted before the execution of the function.
#' @param filterNA (logical): the filterNA() parameter of the omit function.
#' @param coll (character): the variable name of the collection identifiers. (optional, only for filtering!)
#' @param ref (character): the variable name of the reference identifiers. (optional, only for filtering!)
#' @examples
#'	# import data
#'	  data(scleractinia)
#'	  data(stages)
#'
#'	# calculate metrics of diversity dynamics
#'    dd <- divDyn(scleractinia, tax="genus", bin="slc")
#'
#'	# plotting
#'	  plotTS(stages, shading="series", boxes="per", xlim=c(260,0), 
#'	    ylab="range-through diversity (genera)", ylim=c(0,230))
#'	  lines(stages$mid, dd$divRT, lwd=2)
#' 
#'  # with omission of single reference taxa  
#'    ddNoSing <- divDyn(scleractinia, tax="genus", bin="slc", om="ref")
#'    lines(stages$mid, ddNoSing$divRT, lwd=2, col="red")
#'    
#'    # legend
#'    legend("topleft", legend=c("all", "no single-ref. taxa"), 
#'      col=c("black", "red"), lwd=c(2,2))
#'    
#'
#' @export
divDyn <- function(dat, tax="occurrence.genus_name", bin="bin", coll="collection_no", ref="occurrence.reference_no", om=NULL,noNAStart=F, inf=F, data.frame=T, filterNA=FALSE)
{
	
	# is binVar numeric
	if(!is.numeric(dat[,bin])){
		stop("The bin variable is not numeric.")
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
	
	# here comes the counts variables
	counts <- Counts(taxVar,binVar)

	# the metrics
	metrics <- Metrics(counts)
							 
	if(inf!=T)
	{
		metrics[is.infinite(metrics)]<-NA
		
	}
	
	# cbind all together
	dCountsAndMetrics<-cbind(bin=nTimeSlice,  counts, metrics)
	
		
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

	#want to see the NA's at the beginnning? (when the time series does not start with bin 1)
		if (missing(noNAStart)) {}
		else
		{
			if (noNAStart==TRUE)
			{
				dCountsAndMetrics<-dCountsAndMetrics[nStart:nEnd,]
			}
			
			if (noNAStart!=TRUE & noNAStart!=FALSE)
			{
				print("You have entered an invalid argument for noNAStart, no cropping will occurr in the final table.")
			}
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


#' Omission of occurrences that belong to poorly sampled taxa
#' 
#' Function to quickly omit single-collection and single-reference taxa.
#' 
#' The function returns a logical vector of the rows to be omitted. The function is embedded in the divDyn() function, but can be called independently.
#' 
#' @param bin (character value): The name of the subsetting variable (has to be integer). For time series, this is the time-slice variable. If set to NULL, the function performs unbinned subsampling.
#' 
#' @param tax (character value): The name of the taxon variable.
#' 
#' @param dat (data.frame): Occurrence dataset, with bin, tax and coll as column names.
#' @param coll (character value): the variable name of the collection identifiers. 
#' @param ref (character value): the variable name of the reference identifiers. 
#' @param om (character value): the type of omission. "coll" omits occurrences of taxa that occurr only in one collection. "ref" omits occurrences of taxa that were described only in one reference. "binref"  will omit the set of single reference taxa that were described by more than one references, but appear in only one reference in a time bin.
#' @param filterNA (logical value): additional entries can be added to influence the dataset that might not have reference or collection information (NA entries). These occurrences are treated as single-collection or single-reference taxa if the na.rm argument is set to FALSE (default). Setting this argument to TRUE will keep these entries. (see example)
#' 
#' @examples
#' # omit single-reference taxa
#'   data(scleractinia)
#'   toOmit <- omit(scleractinia, bin="slc", tax="genus", om="ref")
#'   dat <- scleractinia[!toOmit,]
#' 
#' # within divDyn
#'   # plotting
#'	  plotTS(stages, shading="series", boxes="per", xlim=c(260,0), 
#'	    ylab="range-through diversity (genera)", ylim=c(0,230))
#'   # multiple ref/slice required
#'   ddNoSing <- divDyn(scleractinia, tax="genus", bin="slc", om="binref")
#'   lines(stages$mid, ddNoSing$divRT, lwd=2, col="red")
#'
#'   # with the recent included (NA reference value)
#'   ddNoSingRec <- divDyn(scleractinia, tax="genus", bin="slc",
#'     om="binref", filterNA=TRUE)
#'   lines(stages$mid, ddNoSingRec$divRT, lwd=2, col="blue")
#'   
#'   # legend
#'   legend("topleft", legend=c("no single-ref. taxa", 
#'     "no single-ref. taxa,\n with recent"), 
#'     col=c("red", "blue"), lwd=c(2,2))

#' @export
omit <- function(dat, tax="occurrence.genus_name", bin="bin", coll="collection_no", ref="occurrence.reference_no", om="ref", filterNA=FALSE){

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
		
		rows<-1:nrow(dat)
	
		tap<-tapply(INDEX=dat[nonDupl,bin], X=rows[nonDupl], function(x){
			sliceTax<-dat[x,tax]
			
			tabSing <- table(dat[x,tax])
		
			singTax<- names(tabSing)[tabSing==1]
			
			x[sliceTax%in%singTax]
		
		})
		
		boolEnd<-rep(FALSE, length(rows))
		boolEnd[unlist(tap)]<-TRUE
		
		# if na.rm TRUE than do not omit NA reference stuff
		if(filterNA){
			boolEnd <- boolEnd & !is.na(dat[,ref])
		}
		
	}
	
	return(boolEnd)

}
