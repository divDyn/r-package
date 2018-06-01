#' Serial subsampling wrapper function
#' 
#' The function will take a desired function that takes an occurrence dataset as an argument, and reruns it iteratively on its subsets. 
#' 
#' The procedure of 'subsampling' or 'sampling standarization' has to be applied in cases 
#' In paleontological questions... 
#' 
#' The function is used to test whether a specific statement holds when sampling standardization is applied. The procedure calculates the variables in question given a certain sampling standardization and routine and level.
#' The implementation of SQS has more variables than necessary, as I intended to reproduce the work of Alroy, which was a continuous process. This method will be simplified in the future after explicit simulation trials.
#' 
#' 
#' @param q (numeric value): Subsampling level argument (mandatory). Depends on the subsampling function, it is the number of occurrences for "cr"
#' 
#' @param dat (data.frame): Occurrence dataset, with bin, tax and coll as column names.
#' 
#' @param iter (numeric value): The number of iterations to be executed.
#' 
#' @param bin (character value): The name of the subsetting variable (has to be integer). For time series, this is the time-slice variable. If set to NULL, the function performs unbinned subsampling.
#' 
#' @param tax (character value): The name of the taxon variable.
#' 
#' @param useFailed (logical): If the interval does not reach the subsampling quota, should the data be used? 
#' @param method (character value): The type of subsampling to be implemented. By default this is classical rarefaction ("cr"). "oxw" stands for occurrence weighted by list subsampling. If set to "sqs", the program will execute the shareholder quorum subsampling algorithm as it was suggested by Alroy (2010).
#' 
#' @param FUN The function to be iteratively executed on the results of the subsampling trials. If set to NULL, no function will be executed, and the subsampled datasets will be returned as a list. By default set to the divDyn() function. The function must have an argument called 'dat', that represents the dataset resulting from a subsampling trial (or the entire dataset). Arguments of the subsample() function call will be searched for potential arguments of this function, which means that already provided variables (e.g. 'bin' and 'tax') will also be used. You can also provide additional arguments (similarly to the apply iterator). Functions that allow arguments to pass through (that have argument '...') are not allowed.
#' 
#' @param output (character value): If the function output are vectors or matrices, the 'arit' and 'geom' values will trigger simple averaging with arithmetic or geometric means. If the function output of a single trial is again a vector or a matrix, setting the output to 'dist' will return the calculated results of every trial, organized in a list of independent variables (e.g. if the function output is value, the return will contain a list vectors, if it is a vector, the output will be a list of vectors, if the function output is a data.frame, the output will be a list of matrices). If output="list", the structure of the original function output will be retained, and the results of the individual trials will be concatenated to a list.
#' 
#' @param intact (numeric vector): The bins, which will not be subsampled but will be added to the subsampling trials. Negative values will be treated as indications on which bins to omit. If the number of occurrences does not reach the subsampling quota, by default it will not be represented in the subsampling trials. You can force their inclusion with the intact argument.
#' 
#' @param implement (character value): Either "lapply", "for" or "foreach". The iterator function of the subsampling trials. Use "lapply" for speed, but this may result in memory issues. Currently only the "for" implementation is finished.
#' 
#' @param duplicates (logical value): Toggles whether multiple entries from the same taxon ("tax") and collection ("coll") variables should be omitted. Useful for omitting occurrences of multiple species-level occurrences of the same genus.
#' @param coll (character value): the variable name of the collection identifiers. 
#' 
#' @param ... arguments passed to FUN and the method-specific subsampling functions:
#' 
#' @examples
#' 
#' data(corals)
#' data(stages)
#' # Example 1-calculate metrics of diversity dynamics
#'   dd <- divDyn(corals, tax="genus", bin="slc")
#'   rarefDD<-subsample(corals,iter=50, q=50,
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
#'   # compare with previous function (is obsolete)
#'   rarefDD <- crDD(corals, quota=50, iter=100, intactBins=95)
#'   lines(stages$mid,rarefDD$divRT, col="red", lwd=2)
#'   
#'   
#' 
#' # Example 2-SIB diversity (§correct the indexing!)
#' # draft a simple function to calculate SIB diversity
#' sib<-function(dat, bin, tax){
#'   calc<-tapply(INDEX=dat[,bin], X=dat[,tax], function(y){
#'     length(levels(factor(y)))
#'   })
#'   return(calc[as.character(stages$num)])
#' }
#' sibDiv<-sib(corals, bin="slc", tax="genus")
#' 
#' # calculate it with subsampling
#' rarefSIB<-subsample(corals,iter=50, q=50,
#'   tax="genus", bin="slc", output="arit", intact=95, FUN=sib)
#' rarefDD<-subsample(corals,iter=50, q=50,
#'   tax="genus", bin="slc", output="arit", intact=95)
#' 
#' # plot
#' plotTS(stages, shading="series", boxes="per", xlim=c(260,0), 
#'   ylab="SIB diversity (genera)", ylim=c(0,230))
#' 
#' lines(stages$mid, rarefDD$divSIB, lwd=2, col="black")
#' lines(stages$mid, rarefSIB, lwd=2, col="blue")
#'     
#' 
#' # Example 3 - different subsampling methods with default function (divDyn)
#' # compare different subsampling types
#'   # classical rarefaction
#'   cr<-subsample(corals,iter=50, q=20,tax="genus", bin="slc", output="dist", intact=95)
#'   # by-list subsampling (unweighted) - 3 collections
#'   UW<-subsample(corals,iter=50, q=3,tax="genus", bin="slc", output="dist", intact=95, method="oxw", x=0)
#'   # occurrence weighted by list subsampling
#'   OW<-subsample(corals,iter=50, q=20,tax="genus", bin="slc", output="dist", intact=95, method="oxw", x=1)
#'  
#'   SQS<-subsample(corals,iter=50, q=0.4,tax="genus", bin="slc", output="dist", intact=95, method="sqs", ref="reference_no")
#'
#' # plot
#'   plotTS(stages, shading="series", boxes="per", xlim=c(260,0), 
#'   ylab="range-through diversity (genera)", ylim=c(0,100))
#'   shades(stages$mid, cr$divRT, col="red")
#'   shades(stages$mid, UW$divRT, col="blue")
#'   shades(stages$mid, OW$divRT, col="green")
#'   shades(stages$mid, SQS$divRT, col="cyan")
#'   
#'   legend("topleft", bg="white", legend=c("CR (20)", "UW (3)", "OW (20)", "SQS (0.4)"), 
#'     col=c("red", "blue", "green", "cyan"), lty=c(1,1,1,1), lwd=c(2,2,2,2))
#' 
#' @rdname subsample
#' @export
subsample <- function(dat, q, tax="genus", coll="collection_no", bin="SLC", FUN=divDyn,  iter=50,  method="cr", intact=NULL, duplicates=FALSE,  output="arit", implement="for", useFailed=FALSE, ...){
	
#	bin <- "slc"
#	tax<- "genus"
#	dat <- corals
#	q<-40
#	output<- "arit"
#	intact<-c(94,95)
	# defense
	
	quoVar<-q
	
	if(!output%in%c("arit", "geom", "list", "dist")) stop("Invalid output argument.")
	
	# prepare for the fact that method should be a function
	# distribute method -> 

	# method specific defence
	if(method=="sqs"){
		if(quoVar>1 | quoVar<=0) stop("Shareholder Quorum has to be in the 0-1 interval.")
	}
	if(method%in%c("cr", "oxw")){
		if(quoVar<=1) stop("The quota has to be a natural number larger than 1.")
	}
	
	# should be null or positive integers
	if(!is.null(bin)){
		bVar<-dat[,bin]
		if(sum(is.na(bVar))>0) stop("The bin column contains missing (NA) values.")
		if(sum(bVar<=0)>0 | sum(bVar%%1)>0) stop("The bin column may only contain positive integers")
	
	}
	
	# check the presences of vectors
	if(!duplicates){
		bDupl<-duplicated(dat[, c(tax, coll)])
		dat<-dat[!bDupl,]
	}
	
	# before anything, remove the unneeded rows (no taxon/bin)
	dat<-dat[!is.na(dat[,tax]) & !is.na(dat[,bin]),]
	
	# what to remove and include (intact argument)
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
				# match function
				FUN<-match.fun(FUN)
				
				# match function arguments to the big call
					# all function arguments
					allArgs<-as.list(match.call()) 
					appliedArgs<-formals(FUN)
					
					# arguments already set
					overlapArgs<-allArgs[names(allArgs)%in%names(appliedArgs)]
					
					appliedArgs[names(overlapArgs)]<-overlapArgs
					
					appliedArgs$dat<-dat
					
					wholeRes<-do.call(FUN, appliedArgs)
			}else{
			#	if(FUN=="divDyn"){
			#		wholeRes<-divDyn(dat, tax=tax, bin=bin)
			#	}else{
					stop("Invalid FUN argument.")
			#	}
			}
		}else{
			output<- "list"
			# empty list to save the subsampling output
			appliedArgs<-list()
			
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
	
	
	
	# subsampling specific preparatory phase
	# CR
	if(method=="cr"& !is.null(bin)){
		# order the rows 
		oBin <- order(dat[,bin])
		
		# reorder dataset
		dat <- dat[oBin,]
		
		# the only important variable
		binVar <- dat[,bin]
		
		# intact rows
		keepRows<-binVar%in%keep
		
		# just organize
		oneResult<-list()
	
	}
	
	
	# SQS
	if(method=="sqs" & !is.null(bin)){
		# distribute arguments
			addArgs<-list(...)
			
		if(!"vers"%in%names(addArgs)){
			vers<-"inexact"
			
		}else{
			# extract argument if added
			vers<-addArgs$vers
			addArgs<-addArgs[names(addArgs)!="vers"]
		}
		if(vers=="inexact"){	
			
			# present in the prepSQS
			prepSQSargs<-c("ref", "singleton","excludeDominant", "largestColl", "fcorr")
			
			# get the arguments for prepSQS
			prepArgs<-addArgs[names(addArgs)%in%prepSQSargs]
			prepArgs<-c(prepArgs, list("dat"=dat, "bin"=bin, "tax"=tax, "coll"=coll))
			
			# frequencies
			freqVar<-do.call(frequencies, prepArgs)
			
			# argument distribution for the subsampler (sqsinexact)
			sqsArgs<-addArgs[names(addArgs)%in%c("byList", "appr", "trialRet")]
			
			sqsArgs<-c(sqsArgs,list(
				"binVar"=dat[,bin],
				"freqVar"=freqVar,
				"collVar"=dat[, coll],
				"q"=quoVar,
				"intact"=keep
			))
		}
		
		if(vers=="exact"){
			# distribute arguments
			addArgs<-list(...)
			
			# the reference column argument
			if("ref"%in%names(addArgs)){
				ref<-addArgs$ref
			}else{
				ref<-"occurrence.reference_no"
			}
			
			sqsArgs<-addArgs[names(addArgs)%in%c("byList", "intact", "singleton")]

			sqsArgs<-c(sqsArgs, list(
				"binVar"=dat[,bin],
				"taxVar"=dat[,tax],
				"collVar"=dat[,coll],
				"refVar"=dat[,ref],
				"q"=quoVar
			
			))

		}
		
	}	
	
	# later: replace all character entries with integers. will allow faster implementation and later C++ conversion
	
	# matrix to keep track failure
	failMat<-matrix(FALSE, ncol=iter, nrow=max(dat[, bin]))
	
	
	#b. iteration
	if(implement=="for"){
	
		for(k in 1:iter){
			#1. produce one subsample
	#		if(method=="cr" & !is.null(bin)){
	#			oneResult<-subsampleCR(binVar=dat[,bin], q=quoVar, intact=keep)
	#			appliedArgs$dat<-dat[oneResult$rows,]
	#		}
			
			if(method=="cr" & !is.null(bin)){
			
				# output of cpp function
				binCRres<-.Call('_divDyn_CRbinwise', PACKAGE = 'divDyn', binVar, quoVar)
			
				# increase indexing (R)
				binCRres[,1]<-binCRres[,1]+1
			
				# rows where quota not reached
				failOutput <- binCRres[,1]==-8
				
				# bins where quota not reached
				oneResult$fail<-binCRres[failOutput,2]
				
				# the rows where 
				usedRows<-rep(FALSE, nrow(dat))
				
				# output by subsampling
				usedRows[binCRres[!failOutput,1]] <- TRUE
				
				usedRows[keepRows] <- TRUE
				# should be kept intact
				appliedArgs$dat<-dat[usedRows,]
				
			}
			
			if(method=="oxw"& !is.null(bin)){
				oneResult<-subsampleOXW(binVar=dat[,bin], collVar=dat[, coll], q=quoVar, intact=keep,...)
				appliedArgs$dat<-dat[oneResult$rows,]
			}
			
			if(method=="sqs"& !is.null(bin)){
			#	oneResult<-subsampleSQS2010(binVar=dat[,bin], freqVar=freqVar, collVar=dat[, coll], q=quoVar, intact=keep)
				if(vers=="inexact"){
					oneResult<-do.call(subsampleSQSinexact, args=sqsArgs)
					appliedArgs$dat<-dat[oneResult$rows,]
				}
				if(vers=="exact"){
					oneResult<-do.call(subsampleSQSexact, args=sqsArgs)
					appliedArgs$dat<-dat[oneResult$rows,]
				}
			
			}
			
			# add the failed stuff is so required
			if(useFailed){
				appliedArgs$dat<-rbind(appliedArgs$dat, dat[dat[,bin]%in%oneResult$fail,])
			}else{
			# keep track
				failMat[oneResult$fail,k] <-TRUE
			}
			
			#2. run the function on subsample
			if(is.null(FUN)){
				trialRes<-appliedArgs$dat
			}else{
				if(is.function(FUN)){
					trialRes<-do.call(FUN,appliedArgs)
				}else{
				#	if(FUN=="divDyn"){
				#		trialRes<- divDyn(trialDat, tax=tax, bin=bin)
				#	} # else has been defended before
				}
			}
				
			
			#3. save the results depending on the output type
			if(out=="vector"){
				if(sum(names(trialRes)%in%rownames(holder))==length(trialRes)){
				#	holder[names(trialRes),k]<-trialRes
					holder[,k]<-trialRes
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
		cat(k, "\r")
		flush.console()
		}
		
	}
	
#	if(implement=="lapply" & method=="cr"){
#		listSource<-rep(list(list(binVar,quoVar, keepRows)), iter)
#		allResults<-lapply(listSource, function(x){
#			# get one trial dataset
#			oneResult<-do.call(subsampleCRpp,x)
#			
#			# add it to function args
#			appliedArgs$dat<-dat[oneResult$row,]
#					
#			# add the failed stuff is so required
#			if(useFailed){
#				appliedArgs$dat<-rbind(appliedArgs$dat, dat[dat[,bin]%in%oneResult$fail,])
#			}
#			
#			#2. run the function on subsample
#			if(is.null(FUN)){
#				trialRes<-appliedArgs$dat
#			}else{
#				if(is.function(FUN)){
#					trialRes<-do.call(FUN,appliedArgs)
#				}
#			}
#			
#			return(list(res=trialRes, fail=oneResult$fail))
#		})
#		
#	}
	
	# Checkpoint2
#	return(holder)

	if(!useFailed){
		failedBins<-which(apply(failMat,1, sum)>0)
	}

	#c. averaging and return
	if(out=="vector"){
		if(output=="arit"){
			totalResult<-apply(holder,1,mean, na.rm=T)
			if(!useFailed){
				totalResult[failedBins]<-NA
			}
		}
		if(output=="geom"){
			loggedVars<-log(holder)
			loggedVars[is.infinite(loggedVars)]<-NA
			totalResult<- apply(loggedVars, 1, mean, na.rm=T) 
			totalResult<- exp(totalResult)
			
			if(!useFailed){
				totalResult[failedBins]<-NA
			}
		}
		if(output=="dist"){
			totalResult<-holder
			if(!useFailed){
				totalResult[failedBins,]<-NA
			}
		}
	}
	
	if(out=="matrix"){
	
		if(output=="arit"){
			totalResult<-apply(holder,c(1,2),mean, na.rm=T)
			
			if(!useFailed){
				totalResult[failedBins,]<-NA
			}
		}
		
		if(output=="geom"){
			loggedVars<-log(holder)
			loggedVars[is.infinite(loggedVars)]<-NA
			totalResult<- apply(loggedVars, c(1,2), mean, na.rm=T) 
			totalResult<- exp(totalResult)
			if(!useFailed){
				totalResult[failedBins,]<-NA
			}
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
						if(!useFailed){
							totalResult[failedBins,]<-NA
						}
					}
					
					if(output=="geom"){
						loggedVars<-log(as.matrix(varDat))
						loggedVars[is.infinite(loggedVars)]<-NA
						totalResult[,i]<- apply(loggedVars, 1, mean, na.rm=T) 
						totalResult[,i]<- exp(totalResult[,i])
						
						if(!useFailed){
							totalResult[failedBins,]<-NA
						}
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


#' @rdname subsample
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
	
	# what was not enough for the subsampling
	failed<-unlist(lapply(tap, is.null))
	failed<-as.numeric(names(failed[failed]))
	failed<-failed[!failed%in%intact]
	
	# the rows that should be passed
	trialRows <- rows%in%unlist(tap)
	
	# these rows should be present regardless of the subsampling
	intactRows<-binVar%in%intact
	
	# combine the two
	subsampleRows<-trialRows | intactRows

	
	# result
	res<-list(rows=subsampleRows, fail=failed)

	# return
	return(res)
}



#' @param x the exponent of by-list subsampling (rarefaction) method.
#' @rdname subsample
subsampleOXW<-function(binVar, collVar, q, intact=NULL,x=1){
	
	# future argument: appr represents how the subsampling quota is approached
	# appr can be "over", "under", "optimize"
	appr="over"
	
	# the chosen quota
	if(x!=0){
		newQ<-q^x
	}else{
		newQ<-q
	}
	rows<- 1:length(binVar)
	
	#keep track of not enough
	ne<-NULL
	tap<-tapply(X=rows, INDEX=binVar, FUN=function(y){
	#	y<-rows[binVar==i]
		# because indices are subsetted
		collection<-collVar[y]
		
		# get the weight of each collection
		expTab<-table(collection)^x
		
		# shuffle
		expTab<-sample(expTab)
		
		# add up
		cumulative<-cumsum(expTab)


		# there are enough 
		if(length(cumulative)>0){
			if(cumulative[length(cumulative)]>newQ){
				bSelect<-cumulative<=newQ
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
	failed<-failed[!failed%in%intact]
	
	# the rows that should be passed
	trialRows <- rows%in%unlist(tap)
	
	# these rows should be present regardless of the subsampling
	intactRows<-binVar%in%intact
	
	# combine the two
	subsampleRows<-trialRows | intactRows
	
	# result
	res<-list(rows=subsampleRows, fail=failed)

	# return
	return(res)

}




#' @param excludeDominant (logical) This parameter sets whether the dominant taxon should 
#'	be excluded from all calculations involving frequencies (this is the second correction of Alroy, 2010).
#' @param largestColl (logical) This parameter sets whether the occurrences of taxa only ever
#'  found in the most diverse collection should be excluded from the count of 
#'	single-publication occurrences. (this is the third correction of Alroy, 2010)
#' @param fcorr (character value) either "good" or "alroy". This argument changes the frequency correction procedure of the 
#'  'inexact' version of SQS (Alroy 2010). As not all taxa are present in the samples, 
#'  the sampled frequencies of taxa tend overestimate their frequencies in the sampling pool. 
#'  In Alroy (2010) these are corrected using Good's u ("good", default), in the later versions 
#'  of SQS this metric is changed to a different method using single occurrence and double occurrence taxa ("alroy"). 
#'@rdname subsample
# excludeDominant sets the Alroy's 2nd correction, sets the frequency of the dominant taxon's to 0
frequencies<-function(dat,bin, tax, coll, ref="reference_no", singleton="ref", excludeDominant=TRUE, largestColl=TRUE, fcorr="good"){

	# the number of occurrences
		O<-table(dat[,bin])
		rowIndex<-1:nrow(dat)
	
	#calculate the frequencies
	if(!excludeDominant){
		# calculate the taxon frequencies in the bin, and assign the value to each occurrence
		freqVarList <-tapply(INDEX=dat[,bin], X=rowIndex, function(x){
		#	x<-rowIndex[unique(dat[,bin])[1]==dat[,bin]]
			
			sliceTaxa<-dat[x,tax]
			taxTab<-table(sliceTaxa)
			freq<-taxTab/length(sliceTaxa)
			
			occWise <- freq[sliceTaxa]
			names(occWise)<-x
			return(occWise)
		
		})
		
		# frequencies
		freqVar<-unlist(freqVarList, use.names=F)
		
		# row identifiers
		nameVec<-unlist(lapply(freqVarList,function(x){
			names(x)
	
		}))
		names(freqVar)<-nameVec
	
		# reorder to reflect the order of occurrences
		freqVar<-freqVar[as.character(rowIndex)]
		
	# adjust dominant to 0
	}else{
		freqVarList <-tapply(INDEX=dat[,bin], X=rowIndex, function(x){
			sliceTaxa<-dat[x,tax]
			taxTab<-table(sliceTaxa)
			freq<-taxTab/length(sliceTaxa)
			domIndex<-which((max(freq)==freq))
			freq[domIndex[1]]<-0
			
			occWise <- freq[sliceTaxa]
			names(occWise)<-x
			return(occWise)
		
		})
		
		freqVar<-unlist(freqVarList, use.names=F)
		nameVec<-unlist(lapply(freqVarList,function(x){
			names(x)
	
		}))
		names(freqVar)<-nameVec
	
		# reorder to reflect the order of occurrences
		freqVar<-freqVar[as.character(rowIndex)]
		
	}
	
	# Alroy 1st correction
	if(singleton%in%c("ref", "occ")){
		
		# number of sampled taxa
		vectS<-tapply(INDEX=dat[,bin], X=dat[,tax], function(x){
			return(length(unique(x)))	
		})
			
		# use single reference taxa
		if(singleton=="ref"){
			# count single reference taxa
			refTax<-unique(dat[,c(ref,tax)])
			tTax<-table(refTax[,tax])
			
			# single reference taxa
			singRefTaxa<-names(tTax)[tTax==1]
			
			# the number of occurrences of single reference taxa
			p1<-tapply(INDEX=dat[,bin], X=rowIndex, function(x){
				sliceTax<-dat[x,tax]
				return(sum(sliceTax%in%singRefTaxa))
			})
			
			# two reference taxa
			twoRefTaxa<-names(tTax)[tTax==2]
			
			# the number of occurrences of single reference taxa
			p2<-tapply(INDEX=dat[,bin], X=rowIndex, function(x){
				sliceTax<-dat[x,tax]
				return(sum(sliceTax%in%twoRefTaxa))
			})
			
			
		}
		
		# use single occurrence taxa
		if(singleton=="occ"){
			occTax<-unique(dat[,c(tax, bin, coll)])
			p1<-tapply(INDEX=occTax[,bin], X=occTax[,tax], function(x){
				tTax<-table(x)
				return(sum(tTax==1))
			})
			
			p2<-tapply(INDEX=occTax[,bin], X=occTax[,tax], function(x){
				tTax<-table(x)
				return(sum(tTax==2))
			})
			
		
		}
	
		# Alroy 2nd correction
		if(excludeDominant){
			# count the dominant taxon 
			n1 <-tapply(INDEX=dat[,bin], X=rowIndex, function(x){
				sliceTaxa<-dat[x,tax]
				taxTab<-table(sliceTaxa)
				domIndex<-which((max(taxTab)==taxTab))
				ret<-taxTab[domIndex]
				return(ret[1])
			})
			
			# then set 
			# for the original solution
			uP <-(O-p1-n1)/(O-n1)
			
			# for Alroy's solution
			correction<-(((p1+p2)/2)/vectS)/(O-n1)
			
			# Alroy's 3rd correction
			if(largestColl){
			
				#number of collections a taxon belong to (total dataset)
				taxTab<-table(dat[,tax])
				
				# taxa in the most diverse collection - that are single-reference and do not occur anywhere else
				tMax<-tapply(INDEX=dat[,bin], X=rowIndex, function(x){
				#	x<-rowIndex[dat[,bin]==61]
				
					# unique collection/taxon table
					collTax<-unique(dat[x,c(coll,tax)])
					# how many taxon/collection
					tColl<-sort(table(collTax[,coll]), decreasing=T)
					
					# which are the taxa that are there
					taxInMost<-collTax[collTax[,coll]==names(tColl[1]), tax]
					
					# which are single reference taxa
					taxInMostAndSingRef<-taxInMost[taxInMost%in%singRefTaxa]
					bTaxInMostAndSingRef<-names(taxTab)%in%taxInMostAndSingRef
					
					# taxa only ever found in the most diverse collection
					return(sum(taxTab==1 & bTaxInMostAndSingRef))
					
					
				})
				
				# calculate u'
				uP <- (O-p1-n1+tMax)/(O-n1)
				
				# for Alroy's solution (exclude taxa found in the largest collection from the single-reference taxa)
				correction<-((((p1-tMax)+p2)/2)/vectS)/(O-n1)
			}
		
		}else{
			uP = 1-p1/O
			correction<-(((p1+p2)/2)/vectS)/O
		}
		# map uP to freqVar and adjust
		if(fcorr=="good"){
			freqVar<-freqVar*uP[as.character(dat[,bin])]
		}
		if(fcorr=="alroy"){
			
			freqVar <- freqVar-correction[as.character(dat[,bin])]
		
		}
	}
	names(freqVar)<-dat[,tax]
	return(freqVar)
}

#' @param singleton (character value): Either "ref" or "occ". If set to "occ", the coverage estimator (e.g. Good's u) will be calculated based on the number of single-occurrence taxa. 
#' If set to "ref" the number of occurrences belonging to single-reference taxa will be used instead. In case of the inexact algorithm, if set to FALSE, then coverage corrections of frequencies will not be applied (not advised).
#' @rdname subsample
subsampleSQSexact<-function(binVar, q, taxVar, collVar, refVar, byList=FALSE, intact=NULL, singleton="ref"){
	rows<- 1:length(binVar)
	#keep track of not enough
	ne<-NULL

	# by occurences approach
	if(!byList){
		# row indices of all time slices
		tap<-tapply(X=rows, INDEX=binVar, FUN=function(y){
		
			shuffledRows<-sample(y)
			taxa<-factor(taxVar[shuffledRows])
			refs<-factor(refVar[shuffledRows])
			colls<-factor(collVar[shuffledRows])
			
			# the number of occurrences in the slice
			sliceOcc<-length(y)
		
			if(singleton=="ref"){
				#refTax<-cbind(refs,taxa)
				refTax<-paste(refs, taxa)
				uVect<-sapply(1:sliceOcc, function(O){
					firstTax<-taxa[1:O]
					
					# Alroy's u'
					1-sum(table(firstTax[!duplicated(refTax[1:O])])==1)/O
					
				})
			}
			if(singleton=="occ"){
				colTax<-paste(colls, taxa)
				uVect<-sapply(1:sliceOcc, function(O){
					firstTax<-taxa[1:O]
					
					# Good's u
					1-sum(table(firstTax[!duplicated(colTax[1:O])])==1)/O
					
				})
			}
			
			# where does uVect cross the quorum? 
			zeroDiff<-uVect-q
			
			# failed trial (the target quorum is never reached)
			if(sum(zeroDiff>=0)==0){
				return(NULL)
			}else{
				
				# exactly the quorum sampled
				indices<-which(zeroDiff==0)
				
				# more likely going above and below
				signs<-sign(zeroDiff)
				
			
				signDiff<-c(NA,diff(signs))
				
				# all indices that are likely
				indices<-unique(c(indices,which(signDiff!=0), which(signDiff!=0)-1))
				
				# which is the middle?
				middle<-median(indices)
				
				#more
				if(!middle%in%indices){
					closeness<-abs(indices-middle)
					take<-which(min(closeness)==closeness)
				
					# single closest
					if(length(take)==1){
						middle<-indices[take]
					}else{
						# select a random one, it will probably notmatter
						middle<-sample(indices[take],1)
					}
				
				}
			
			
				# return indices
				return(shuffledRows[1:middle])
			}
			
#		}
		
		})
	}
	
	# what was not enough for the subsampling
	failed<-unlist(lapply(tap, is.null))
	failed<-as.numeric(names(failed[failed]))
	failed<-failed[!failed%in%intact]
	
	# the rows that should be passed
	trialRows <- rows%in%unlist(tap)
	
	# these rows should be present regardless of the subsampling
	intactRows<-binVar%in%intact
	
	# combine the two
	subsampleRows<-trialRows | intactRows

	
	# result
	res<-list(rows=subsampleRows, fail=failed)

	# return
	return(res)

}


#' @param byList (character value): A parameter of SQS. Sets whether occurrences should be subsampled with (FALSE) or without (TRUE) breaking the collection integrity. (not yet for the exact algorithm.)
#'	
#' @param appr (character value): A parameter of SQS, either "over" (default) or ("under"). The current 
#' version is not concerned with small fluctuations around the drawn subsampling quorum. 
#' Therefore, in the inexact algorithm, sampling is finished when the subset 
#' either is immediately below the quorum ("under") or above it ("over").
#' @rdname subsample
subsampleSQSinexact<-function(binVar, freqVar, q, collVar=NULL, byList=FALSE, intact=NULL, appr="over", trialRet="occ"){
	
	rows<- 1:length(binVar)
	#keep track of not enough
	ne<-NULL
	
	# non by-list SQS
	if(!byList){
		tap<-tapply(X=rows, INDEX=binVar, FUN=function(y){
			# the frequencies
			frequencies<-freqVar[y]
			
			# the taxa
			taxa<-names(freqVar)[y]
			
			# shuffle, we want to go by occurrences (random order),
			shuffledOccs<-sample(1:length(taxa))
			shuffledSp<-taxa[shuffledOccs]
			shuffledFreq<-frequencies[shuffledOccs]
			
			# but should only count new species' coverage
			speciesFirstEncounterOrder<-shuffledSp[!duplicated(shuffledSp)]
			coverageOrder<-shuffledFreq[!duplicated(shuffledSp)]
			
			# add up
			cumulative<-cumsum(coverageOrder)
			
			if(length(cumulative)>0){
				if(cumulative[length(cumulative)]>q){
					bSelect<-cumulative<=q
					# potential forking!!! 
					if(appr=="over"){
						if(sum(bSelect)!=length(bSelect)){
							bSelect[sum(bSelect)+1] <- TRUE
						}
					}
					
					# the species to keep
					keepTaxa<-speciesFirstEncounterOrder[bSelect]
					
					# what should the function return? all occurrences of the 
					# species, or just those that were skimmed by the algorithm?
					if(trialRet=="occ"){
						# first encounter of the last:
						lastTaken<-keepTaxa[max(which(bSelect))]
						# keep occurrences before that
						lastIndex<-min(which(shuffledSp==lastTaken))
						rowShuffle<-y[shuffledOccs]
						return(rowShuffle[1:lastIndex])
					}
					
					if(trialRet=="taxon"){
						return(y[taxa%in%keepTaxa])
					}
					
				}else{
					return(NULL)
				}
			}else{
				return(NULL)
			}
				
			
		})
	# by-list SQS
	}else{
		tap<-tapply(X=rows, INDEX=binVar, FUN=function(y){
		#	y<-rows[binVar==tempVar[zui]]
			# the frequencies
			frequencies<-freqVar[y]
			
			# the taxa
			taxa<-names(freqVar)[y]
			
			# collections
			colls<-collVar[y]
			
			colTab<-table(taxa, colls)
			class(colTab) <- "matrix"
			
			# ensure that it will run if there is no collection info
			if(ncol(colTab)>1){
				
				# shuffle collections
				newOrd<-sample(1:ncol(colTab))
				colTab<-colTab[,newOrd]
				
				firstOccur<-apply(colTab, 1, function(b){
					min(which(as.logical(b)))
				})
				# empty
				colTab[,]<-0
				
				# frequencies should be calculated when species occur for the first time
				colTab[cbind(1:nrow(colTab),firstOccur)]<-frequencies[names(firstOccur)]
				
				# frequencies to be cumulated
				coverageOrder<-apply(colTab, 2, sum)
				
				# add up
				cumulative<-cumsum(coverageOrder)
				
				if(length(cumulative)>0){
					if(cumulative[length(cumulative)]>q){
						bSelect<-cumulative<=q
						# potential forking!!! 
						if(appr=="over"){
							if(sum(bSelect)!=length(bSelect)){
								bSelect[sum(bSelect)+1] <- TRUE
							}
						}
						
						# the collections to keep
						keepColls<-colnames(colTab)[bSelect]
						
						return(y[colls%in%keepColls])
						
						
					}else{
						return(NULL)
					}
				}else{
					return(NULL)
				}
			}
			
		}) # end tapply
		
	} # end bylist
	
	# what was not enough for the subsampling
	failed<-unlist(lapply(tap, is.null))
	failed<-as.numeric(names(failed[failed]))
	failed<-failed[!failed%in%intact]
	
	# the rows that should be passed
	trialRows <- rows%in%unlist(tap)
	
	# these rows should be present regardless of the subsampling
	intactRows<-binVar%in%intact
	
	# combine the two
	subsampleRows<-trialRows | intactRows

	
	# result
	res<-list(rows=subsampleRows, fail=failed)

	# return
	return(res)

}




