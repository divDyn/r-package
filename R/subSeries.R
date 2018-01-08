#' Serial subsampling wrapper function
#' 
#' The function will take a desired function that takes an occurrence dataset as an argument, and reruns it iteratively on its subsets. 
#' 
#' The procedure of 'subsampling' or 'sampling standarization' has to be applied in cases 
#' In paleontological questions... 
#' 
#' The function is used to test whether a specific statement holds when sampling standardization is applied. The procedure calculates the variables in question given a certain sampling standardization and routine and level.
#' §more
#' 
#' 
#' @param q (numeric value): Subsampling level argument. Depends on the subsampling function. This is mandatory argument.
#' 
#' @param dat (data.frame): Occurrence dataset.
#' 
#' @param iter (numeric value): The number of iterations to be executed.
#' 
#' @param bin (character value): The name of the subsetting variable (has to be numeric). For time series, this is the time-slice variable.
#' 
#' @param tax (character value): The name of the taxon variable.
#' 
#' @param method (character value): The type of subsampling to be implemented. By default this is classical rarefaction ("cr").
#' 
#' @param FUN The function to be iteratively executed during the subsampling trials. IF set to NULL, no function will be executed, and the subsampled datasets will be returned as a list. If set to the default "divDyn", the divDyn(), function will be run with default settings. The function should be defined with a single argument, the database subset resulting from a trial. 
#' 
#' @param output (character value): If the function output are vectors or matrices, the 'arit' and 'geom' values will trigger simple averaging with arithmetic or geometric means. If the function output of a single trial is again a vector or a matrix, setting the output to 'dist' will return the calculated results of every trial, organized in a list of independent variables (e.g. if the function output is value, the return will contain a list vectors, if it is a vector, the output will be a list of vectors, if the function output is a data.frame, the output will be a list of matrices). If output="list", the structure of the original function output will be retained, and the results of the individual trials will be concatenated to a list.
#' 
#' @param intact (numeric vector): The bins, which will not be subsampled but will be added to the subsampling trials. Negative values will be treated as indications on which bins to omit. If the number of occurrences does not reach the subsampling quota, by default it will not be represented in the subsampling trials. You can force their inclusion with the intact argument.
#' 
#' @param implement (character value): Either "lapply", "for" or "foreach". The iterator function of the subsampling trials. Use "lapply" for speed, but this may result in memory issues. Currently only the "for" implementation works.
#' 
#' @param duplicates (logical value): Toggles whether multiple entries from the same taxon ("tax") and collection ("coll") variables should be omitted. Useful for omitting occurrences of multiple species-level occurrences of the same genus.
#' @param coll (character value): the variable name of the collection identifiers. 
#' 
#' @param ... arguments passed to the method-specific subsampling functions.
#' 
#' @examples
#' 
#' data(scleractinia)
#' data(stages)
#' # Example 1-calculate metrics of diversity dynamics
#'   dd <- divDyn(scleractinia, tax="genus", bin="slc")
#'   rarefDD<-subseries(scleractinia,iter=50, q=50,
#'   tax="genus", bin="slc", output="dist", intact=95)
#' 	
#' # plotting
#'   plotTS(stages, shading="series", boxes="per", xlim=c(260,0), 
#'   ylab="range-through diversity (genera)", ylim=c(0,230))
#'   lines(stages$mid, dd$divRT, lwd=2)
#'   shades(stages$mid, rarefDD$divRT, col="blue")
#'   legend("topleft", legend=c("raw","rarefaction"),
#'     col=c("black", "blue"), lwd=c(2,2), bg="white")
#'   
#'   # compare with previous function (will be obsolete)
#'   rarefDD <- crDD(scleractinia, quota=50, iter=100, intactBins=95)
#'   lines(stages$mid,rarefDD$divRT, col="red", lwd=2)
#'   
#'   
#' # Example 2-SIB diversity
#' # draft a simple function to calculate SIB diversity
#' sib<-function(x){
#'   tapply(INDEX=x$slc, X=x$genus, function(y){
#'     length(levels(factor(y)))
#'   })
#' }
#' sibDiv<-sib(scleractinia)
#' 
#' # calculate it with subsampling
#' rarefSIB<-subseries(scleractinia,iter=50, q=50,
#'   tax="genus", bin="slc", output="arit", intact=95, FUN=sib)
#'   
#' @export
subseries <- function(dat, FUN="divDyn", q,iter=10,  bin="SLC",intact=NULL, tax="genus", duplicates=FALSE, coll="collection_no", output="arit", implement="for", method="cr",...){
	
#	bin <- "slc"
#	tax<- "genus"
#	dat <- scleractinia
#	q<-40
#	output<- "arit"
#	intact<-c(94,95)
	# defense
	if(!output%in%c("arit", "geom", "list", "dist")) stop("Invalid output argument.")

	# check the presences of vectors
	if(!duplicates){
		bDupl<-duplicated(dat[, c(tax, coll)])
		dat<-dat[!bDupl,]
	}
	
	# what to remove and include (intact)
	negative<-NULL
	keep<-NULL
	if(!is.null(intact)){
		negative<-sign(intact)==-1
		level<-unique(dat[,bin])
		level<-level[!is.na(level)]
	
		rem<--intact[negative]
		keep<-intact[!negative]
		
		if(sum(!rem%in%level)>0 | sum(!keep%in%level)>0) stop("Invalid intact argument.")
		
		if(length(negative)>0){
			dat<-dat[!dat[, bin]%in%rem,]
		}
			
	}
	
	
	#a. given that the final output is not a list
	if (output!="list"){
		#a. perform the function on the total dataset to assess final structure
		if(!is.null(FUN)){
			if(is.function(FUN)){
					FUN<-match.fun(FUN)
					wholeRes<-FUN(dat)
			}else{
				if(FUN=="divDyn"){
					wholeRes<-divDyn(dat, tax=tax, bin=bin)
				}else{
					stop("Invalid FUN argument.")
				}
			}
		}else{
			output<- "list"
			wholeRes<-NULL
		}
		
			}else{
		wholeRes<-NULL
	}
	# default
	out<-"list"

	if(output!="list"){
	
		# if the output is a vector
		if(is.vector(wholeRes) | is.table(wholeRes)| is.array(wholeRes)){
			if(length(dim(wholeRes))==1){
				holder <- matrix(NA, ncol=iter, nrow=length(wholeRes))
				rownames(holder)<-names(wholeRes)
				out<-"vector"
			}
		}
		# if the output is matrix
		if(is.matrix(wholeRes)| is.array(wholeRes)){
			if(length(dim(wholeRes))==2){
				holder<-array(NA, dim=c(dim(wholeRes),iter))
				colnames(holder) <- colnames(wholeRes)
				rownames(holder) <- rownames(wholeRes)
				out<-"matrix"
			}
		
		}
		#if the output is a data.frame
		if(is.data.frame(wholeRes)){
			nVar<-ncol(wholeRes)
			# repeat original structure
			holder <- matrix(NA, ncol=nVar*iter, nrow=nrow(wholeRes))
			holder<-as.data.frame(holder)
			colnames(holder)<-rep(1:iter, each=nVar)
			out<-"data.frame"
		}
		
		if(out=="list"){
			message("U have to do the averaging yourself")
		}

	}
	
	#checkpoint 1
#	return(holder)
	

	# if it is a list, the averaging has to be done by the user
	if(output=="list" | out=="list"){
		
		# force list output
		out <-"list"
		
		# create a final container
		containList<-list()
	}
	
	#b. iteration
	if(implement=="for"){
	
		for(k in 1:iter){
			#1. produce one subsample
			if(method=="cr"){
				trialDat<-dat[subsampleCR(binVar=dat[,bin], q=q, intact=keep),]
			}
			
			#2. run the function on subsample
			if(is.null(FUN)){
				trialRes<-trialDat
			}else{
				if(is.function(FUN)){
					trialRes<-FUN(trialDat)
				}else{
					if(FUN=="divDyn"){
						trialRes<- divDyn(trialDat, tax=tax, bin=bin)
					} # else has been defended before
				}
			}
				
			
			#3. save the results depending on the output type
			if(out=="vector"){
				if(sum(names(trialRes)%in%rownames(holder))==length(trialRes)){
					holder[names(trialRes),k]<-trialRes
				}else{
					if(nrow(holder)==length(trialRes)){
						holder[,k]<-trialRes
					}else{
						message("Subsampling trial result structure does not match original. Returning a list.")
						out<-"list"
					}
					
				}
				
			}
			if(out=="matrix"){
				if(sum(dim(trialRes)==dim(holder[,,k]))==2){
					holder[,,k] <- trialRes
				}else{
					message("Subsampling trial result structure does not match original. Returning a list.")
					out<-"list"	
				}
			}
			if(out=="data.frame"){
				posVec<-(1:nVar)+((k-1)*nVar)
				if(sum(dim(holder[,posVec])%in%dim(trialRes))==2){
					holder[,posVec] <- trialRes
				}else{
					message("Subsampling trial result structure does not match original. Returning a list.")
					out<-"list"	
				}
			}
			if(out=="list"){
				containList<-c(containList, list(trialRes))
			}
				
		}
		
	}
	# Checkpoint2
#	return(holder)

	#c. averaging and return
	if(out=="vector"){
		if(output=="arit"){
			totalResult<-apply(holder,1,mean, na.rm=T)
		}
		if(output=="geom"){
			loggedVars<-log(holder)
			loggedVars[is.infinite(loggedVars)]<-NA
			totalResult<- apply(loggedVars, 1, mean, na.rm=T) 
			totalResult<- exp(totalResult)
		}
		if(output=="dist"){
			totalResult<-holder
		}
	}
	
	if(out=="matrix"){
	
		if(output=="arit"){
			totalResult<-apply(holder,c(1,2),mean, na.rm=T)
		}
		
		if(output=="geom"){
			loggedVars<-log(holder)
			loggedVars[is.infinite(loggedVars)]<-NA
			totalResult<- apply(loggedVars, c(1,2), mean, na.rm=T) 
			totalResult<- exp(totalResult)
		}
		
		if(output=="dist"){
			totalResult<-holder
		}
	}
	
	if(out=="data.frame"){
	
		# if the original format is a data frame
		totalResult <- matrix(NA, ncol=ncol(wholeRes), nrow=nrow(wholeRes))
		totalResult<-as.data.frame(totalResult)
		colnames(totalResult) <- colnames(wholeRes)
		rownames(totalResult) <- rownames(wholeRes)
		
		# in case the output format should be a list of distributions
		if(output=="dist"){
			totalResult<-list()
		}
		
		for(i in 1:nVar){
			
			# variable specific columns
			varCol<-seq(i,ncol(holder), by=nVar)
			varDat<-holder[,varCol]
				
			if(output=="dist"){
				totalResult<-c(totalResult, list(varDat))
				names(totalResult)[length(totalResult)] <- colnames(wholeRes)[i]
				
			}else{
				# if it is numeric at all!!
				if(!is.character(wholeRes[,i])){
					
					# which meaning should be used
					if(output=="arit"){
						totalResult[,i]<-apply(varDat, 1, mean, na.rm=T)
					}
					
					if(output=="geom"){
						loggedVars<-log(as.matrix(varDat))
						loggedVars[is.infinite(loggedVars)]<-NA
						totalResult[,i]<- apply(loggedVars, 1, mean, na.rm=T) 
						totalResult[,i]<- exp(totalResult[,i])
					}
					
					# some hacking required for the divDyn for the 'outer rates'
					
				}else{
					# some functionality can be added
					totalResult[,i] <- wholeRes[,i]
				}
			}
			
	
		}
	
		
	
	}
	if(out=="list"){
		totalResult<- containList
		names(totalResult)<-paste("trial", 1:iter, sep="")
	}
	
	return(totalResult)
}


subsampleCR <- function(binVar,q,intact=NULL){
#	binVar <- dat[,bin]
	rows<- 1:length(binVar)
	
	tap<-tapply(INDEX=binVar, X=rows, FUN=function(x, q){
		if(length(x)>=q){
			return(sample(x, q, replace=F))
		}else{
			return(NULL)
		}
	
	}, q=q)
	
	# the rows that should be passed
	trialRows <- rows%in%unlist(tap)
	
	# these rows should be present regardless of the subsampling
	intactRows<-binVar%in%intact
	
	# combine the two
	subsampleRows<-trialRows | intactRows

	# return
	return(subsampleRows)
}
