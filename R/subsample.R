#' Serial subsampling wrapper function
#' 
#' The function will take a desired function that has an occurrence dataset as an argument, and reruns it iteratively on the subsets of the dataset. 
#' 
#'
#'
#' The \code{subsample} function implements the iterative framework of the sampling standardization procedure. 
#' The function 1. takes the dataset \code{dat}, 2. runs function \code{FUN} on the dataset and creates a container for results of trials
#' 3. runs one of the subsampling trial functions (e.g. \code{\link{subtrialCR}}) to get a subsampled 'trial dataset'
#' 4. runs \code{FUN} on the trial dataset and
#' 5. averages the results of the trials. For a detailed treatment on what the function does, please see the vignette ('Handout to the R package 'divDyn' v0.5.0 for diversity dynamics from fossil occurrence data'). Currently the Classical Rarefaction (\code{"cr"}, Raup, 1975), the occurrence weighted by-list subsampling (\code{"oxw"}, Alroy et al., 2001) and the Shareholder Quorum Subsampling methods are implemented (\code{"sqs"}, Alroy, 2010).
#'
#' \strong{References:}
#'
#' Alroy, J., Marshall, C. R., Bambach, R. K., Bezusko, K., Foote, M., Fürsich, F. T., … Webber, A. (2001). Effects of sampling standardization on estimates of Phanerozoic marine diversification. Proceedings of the National Academy of Science, 98(11), 6261-6266.
#' 
#' Alroy, J. (2010). The Shifting Balance of Diversity Among Major Marine Animal Groups. Science, 329, 1191-1194. https://doi.org/10.1126/science.1189910
#'
#' Raup, D. M. (1975). Taxonomic Diversity Estimation Using Rarefaction. Paleobiology, 1, 333-342. https: //doi.org/10.2307/2400135
#' 
#' @param q (\code{numeric)}: Subsampling level argument (mandatory). Depends on the subsampling function, it is the number of occurrences for \code{"cr"}, and the number of desired occurrences to the power of \code{x} for O^x^W. It is also the quorum of the SQS method.
#' 
#' @param dat (\code{data.frame}): Occurrence dataset, with \code{bin}, \code{tax} and \code{coll} as column names.
#' 
#' @param iter (\code{numeric}): The number of iterations to be executed.
#' 
#' @param bin (\code{character}): The name of the subsetting variable (has to be integer). For time series, this is the time-slice variable. Rows with \code{NA} entries in this column will be omitted.
#' 
#' @param tax (\code{character}): The name of the taxon variable.
#' @param ref (\code{character}): The name of the reference variable, optional - depending on the subsampling method.
#' @param useFailed (\code{logical}): If the bin does not reach the subsampling quota, should the bin be used? 
#' @param type (\code{character}): The type of subsampling to be implemented. By default this is classical rarefaction (\code{"cr"}). (\code{"oxw"}) stands for occurrence weighted by-list subsampling. If set to (\code{"sqs"}), the program will execute the shareholder quorum subsampling algorithm as it was suggested by Alroy (2010). Setting the argument to \code{"none"} will invoke no subsamling, but the applied function will be iterated on the trials, nevertheless.
#' 
#' @param FUN (\code{function}): The function to be iteratively executed on the results of the subsampling trials. If set to \code{NULL}, no function will be executed, and the subsampled datasets will be returned as a \code{list}. By default set to the \code{\link{divDyn}} function. The function must have an argument called \code{dat}, that represents the dataset resulting from a subsampling trial (or the entire dataset). Arguments of the \code{subsample} function call will be searched for potential arguments of this function, which means that already provided variables (e.g. \code{bin} and \code{tax}) will also be used. You can also provide additional arguments (similarly to the \code{\link[base]{apply}} iterator). Functions that allow arguments to pass through (that have argument '...') are not allowed, as well as functions that have the same arguments as \code{subsample} but would require different values.
#' 
#' @param output (\code{character}): If the function output are vectors or matrices, the \code{"arit"} and \code{"geom"} values will trigger simple averaging with arithmetic or geometric means. If the function output of a single trial is again a \code{vector} or a \code{matrix}, setting the output to \code{"dist"} will return the calculated results of every trial, organized in a \code{list} of independent variables (e.g. if the function output is value, the return will contain a single \code{vector}, if it is a \code{vector}, the output will be a list of \code{vector}s, if the function output is a \code{data.frame}, the output will be a \code{list} of \code{matrix} class objects). If \code{output="list"}, the structure of the original function output will be retained, and the results of the individual trials will be concatenated to a \code{list}.
#' 
#' @param keep (\code{numeric}): The bins, which will not be subsampled but will be added to the subsampling trials. NIf the number of occurrences does not reach the subsampling quota, by default it will not be represented in the subsampling trials. You can force their inclusion with the \code{keep} argument separetely (for all, see the \code{useFailed} argument).
#' @param rem  (\code{numeric}): The bins, which will be removed from the dataset before the subsampling trials.
#' @param duplicates (\code{logical} ): Toggles whether multiple entries from the same taxon (\code{"tax"}) and collection (\code{"coll"}) variables should be omitted. Useful for omitting occurrences of multiple species-level occurrences of the same genus. By default these are allowed through analyses (\code{duplicates=TRUE}), setting this to \code{FALSE} will require you to provide a collection variable. (\code{coll})
#' @param coll (\code{character}): The variable name of the collection identifiers. 
#' 
#' @param ... arguments passed to \code{FUN} and the type-specific subsampling functions: \code{\link{subtrialCR}}, \code{\link{subtrialOXW}}, \code{\link{subtrialSQS}} 
#' 
#' @examples
#' 
#' data(corals)
#' data(stages)
#' # Example 1-calculate metrics of diversity dynamics
#'   dd <- divDyn(corals, tax="genus", bin="stg")
#'   rarefDD<-subsample(corals,iter=50, q=50,
#'   tax="genus", bin="stg", output="dist", keep=95)
#' 	
#' # plotting
#'   tsplot(stages, shading="series", boxes="sys", xlim=c(260,0), 
#'   ylab="range-through diversity (genera)", ylim=c(0,230))
#'   lines(stages$mid, dd$divRT, lwd=2)
#'   shades(stages$mid, rarefDD$divRT, col="blue")
#'   legend("topleft", legend=c("raw","rarefaction"),
#'     col=c("black", "blue"), lwd=c(2,2), bg="white")
#'   
#'   
#'   
#' \donttest{
#' # Example 2-SIB diversity 
#' # draft a simple function to calculate SIB diversity
#' sib<-function(dat, bin, tax){
#'   calc<-tapply(INDEX=dat[,bin], X=dat[,tax], function(y){
#'     length(levels(factor(y)))
#'   })
#'   return(calc[as.character(stages$stg)])
#' }
#' sibDiv<-sib(corals, bin="stg", tax="genus")
#' 
#' # calculate it with subsampling
#' rarefSIB<-subsample(corals,iter=50, q=50,
#'   tax="genus", bin="stg", output="arit", keep=95, FUN=sib)
#' rarefDD<-subsample(corals,iter=50, q=50,
#'   tax="genus", bin="stg", output="arit", keep=95)
#' 
#' # plot
#' tsplot(stages, shading="series", boxes="sys", xlim=c(260,0), 
#'   ylab="SIB diversity (genera)", ylim=c(0,230))
#' 
#' lines(stages$mid, rarefDD$divSIB, lwd=2, col="black")
#' lines(stages$mid, rarefSIB, lwd=2, col="blue")
#' 
#' 
#' # Example 3 - different subsampling types with default function (divDyn)
#' # compare different subsampling types
#'   # classical rarefaction
#'   cr<-subsample(corals,iter=50, q=20,tax="genus", bin="stg", output="dist", keep=95)
#'   # by-list subsampling (unweighted) - 3 collections
#'   UW<-subsample(corals,iter=50, q=3,tax="genus", bin="stg", coll="collection_no", 
#'     output="dist", keep=95, type="oxw", x=0)
#'   # occurrence weighted by list subsampling
#'   OW<-subsample(corals,iter=50, q=20,tax="genus", bin="stg", coll="collection_no", 
#'     output="dist", keep=95, type="oxw", x=1)
#'  
#'   SQS<-subsample(corals,iter=50, q=0.4,tax="genus", bin="stg", output="dist", keep=95, type="sqs")
#'
#' # plot
#'   tsplot(stages, shading="series", boxes="sys", xlim=c(260,0), 
#'   ylab="range-through diversity (genera)", ylim=c(0,100))
#'   shades(stages$mid, cr$divRT, col="red")
#'   shades(stages$mid, UW$divRT, col="blue")
#'   shades(stages$mid, OW$divRT, col="green")
#'   shades(stages$mid, SQS$divRT, col="cyan")
#'   
#'   legend("topleft", bg="white", legend=c("CR (20)", "UW (3)", "OW (20)", "SQS (0.4)"), 
#'     col=c("red", "blue", "green", "cyan"), lty=c(1,1,1,1), lwd=c(2,2,2,2))
#' }
#'
#' @rdname subsample
#' @export
subsample <- function(dat, q, tax="genus", bin="stg",  FUN=divDyn, coll=NULL, ref=NULL, iter=50,  type="cr", keep=NULL, rem=NULL, duplicates=TRUE,  output="arit",  useFailed=FALSE, ...){
	
#	bin <- "stg"
#	tax<- "genus"
#	dat <- corals
#	q<-40
#	output<- "arit"
	# defense
	
	# temporary omission
	implement<- "for"
	quoVar<-q
	
	# function defense
		# iteration
		if(length(iter)!=1) stop("Only a single number of iterations is allowed.")
		if(iter<1 | iter%%1!=0) stop("Only a positive integers are allowed.")
				
		# output
		if(!output%in%c("arit", "geom", "list", "dist")) stop("Invalid output argument.")
		
		#type
		if(!type%in%c("sqs", "cr", "oxw", "none")) stop("Invalid subsampling type.")
	
		# type specific defense
		if(type=="sqs"){
			if(quoVar>1 | quoVar<=0) stop("Shareholder Quorum has to be in the 0-1 interval.")
		}
		if(type%in%c("cr", "oxw")){
			if(quoVar<=1) stop("The quota has to be a natural number larger than 1.")
		}
		

		# the bin identifier - omit BINS
		if(!is.null(bin)){
			bVar<-is.na(dat[,bin])
		
			# omit NA values
			if(sum(bVar)>0) dat <- dat[!bVar,]
		#	if(sum(bVar<=0)>0 | sum(bVar%%1)>0) stop("The bin column may only contain positive integers")
		
		}else{
			stop("You must provide a 'bin' variable.")
		}
		
		
	# list additional arguments for defense of invalid combinations
		addArgs<-list(...)

		if(!is.null(addArgs$largestColl)){
			if(is.null(coll) & addArgs$largestColl) stop("You cannot correct with the largest collection, if you do not provide collection variable.")
		}
		
		if(!is.null(addArgs$singleton)){
			if(!addArgs$singleton%in%c("ref", "occ", FALSE)) stop("Invalid singleton value.")
			if(is.null(ref) & addArgs$singleton=="ref") stop("You cannot calculate coverage based on references, if you do not provide a references variable.")
		}
		
		if(!is.null(addArgs$byList)){
			if(addArgs$byList & is.null(coll)) stop("You cannot do by-list subsampling without providing a collection variable.")
		}
		if(type=="oxw"){
			if(is.null(coll))stop("You cannot do oxw subsampling without providing a collection 'coll' variable.")
			if(sum(is.na(dat[!dat[,bin]%in%keep,coll]))>0) stop("Collection variable 'coll' includes NAs, please omit these rows manually,\n or indicate the bins where subsampling is not desired in 'keep'")
		} 
	
	# cr and unit column
		if("unit"%in%names(addArgs) & type=="cr"){
			if(!addArgs$unit%in%colnames(dat)) stop("Invalid unit column.")
			if(sum(is.na(dat[,addArgs$unit])>0)) stop("Variable 'unit' includes NAs, please omit these rows manually.")
		} 
		

	# the bin levels
		level<-unique(dat[,bin])
		level<-level[!is.na(level)]
				
		
		# coerce keeping intervals
		if(!is.null(keep)){
			if(sum(!keep%in%level)>0) stop("Invalid keep argument.")
		}
		
		if(!is.null(rem)){
			if(sum(!rem%in%level)>0) stop("Invalid rem argument.")
			if(length(rem)>0){
				dat<-dat[!dat[, bin]%in%rem,]
			}
		
		}
	
	
	
	# omit multiple instances of the same genus, if desired!
	if(!duplicates){
		if(length(dat[[coll]])==0) stop("Invalid 'coll' argument. If you want to omit duplicates, provide a valid 'coll' value.")
		bDupl<-duplicated(dat[, c(tax, coll)])
		dat<-dat[!bDupl,]
	}
	
	# before anything, remove the unneeded rows (no taxon/bin)
	dat<-dat[!is.na(dat[,tax]) & !is.na(dat[,bin]),]

	
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
		
				if(output!="list"){				
					wholeRes<-do.call(FUN, appliedArgs)
				}else{
					wholeRes<-NULL
				}
		}else{
			stop("Invalid FUN argument.")
		
		}
	}else{
		output<- "list"
		# empty list to save the subsampling output
		appliedArgs<-list()
		
		# result template won't be necessary
		wholeRes<-NULL
	}
	
	
	# the output object will be 
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
	if(type=="cr"& !is.null(bin)){
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
		
		if("unit"%in%names(addArgs)){
			byUnit <- TRUE
			unitVar <- dat[,addArgs$unit]
			rows <- 1:nrow(dat)
					
		}else{
			byUnit <- FALSE
		}
	}
	
	if(type=="oxw"){

		oxwArgs<-list(
			binVar=dat[,bin],
			collVar=dat[, coll],
			q=quoVar,
			intact=keep
		)
		
		if("x"%in%names(addArgs)){
			oxwArgs<-c(oxwArgs, list(x=addArgs$x))
		}
	}
	
	
	# SQS
	if(type=="sqs" & !is.null(bin)){
		
		# present in the prepSQS
		prepSQSargs<-c("singleton","excludeDominant", "largestColl", "fcorr")
		
		# get the arguments for prepSQS
		prepArgs<-addArgs[names(addArgs)%in%prepSQSargs]
		prepArgs<-c(prepArgs, list("dat"=dat, "bin"=bin, "tax"=tax, "coll"=coll, "ref"=ref))
		
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

	if(type=="none"){
		oneResult<- list()
	}
	
	# later: replace all character entries with integers. will allow faster implementation and later C++ conversion
	
	# matrix to keep track failure
	failMat<-matrix(FALSE, ncol=iter, nrow=length(level))
	rownames(failMat) <- level
	
	#b. iteration
	if(implement=="for"){
	
		for(k in 1:iter){
			#1. produce one subsample
			if(type=="cr" & !is.null(bin)){
			
				# switch, which subsampling method should be used - by-unit or by-row
				#by-row method
				if(!byUnit){
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
				# by-unit
				}else{
								
					tap <- tapply(INDEX=binVar, X=rows, function(x){
						unitSlice<-unitVar[x]
						
						# random sample
						uniqueUnit<-sample(unique(unitVar[x]))
						
						if(quoVar<=length(uniqueUnit)){
							# rows corresponding t
							return(x[unitSlice%in%uniqueUnit[1:quoVar]])
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
					oneResult<-list(rows=subsampleRows, fail=failed)
					
					#the final output
					appliedArgs$dat<-dat[oneResult$rows,]
				
				}
				
			}
			
			if(type=="oxw"& !is.null(bin)){
					
				oneResult<-do.call(subsampleOXW, args=oxwArgs)
				appliedArgs$dat<-dat[oneResult$rows,]
			}
			
			if(type=="sqs"& !is.null(bin)){
			
				oneResult<-do.call(subsampleSQSinexact, args=sqsArgs)
				appliedArgs$dat<-dat[oneResult$rows,]
				
			}

			if(type=="none"){
				appliedArgs$dat <- dat

			}
			
			# add the failed stuff if so required
			if(useFailed){
				appliedArgs$dat<-rbind(appliedArgs$dat, dat[dat[,bin]%in%oneResult$fail,])
			}else{
			# keep track
				failMat[as.character(oneResult$fail),k] <-TRUE
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
		utils::flush.console()
		}
		
	}
	
#	if(implement=="lapply" & type=="cr"){
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
		failedBins<-rownames(failMat)[which(apply(failMat,1, sum)>0)]
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
						# switch warnings off
						warnVal<-options()$warn
						options(warn=-1)
						
						loggedVars<-log(as.matrix(varDat))
						loggedVars[is.infinite(loggedVars)]<-NA
						totalResult[,i]<- apply(loggedVars, 1, mean, na.rm=T) 
						totalResult[,i]<- exp(totalResult[,i])
						
						#warnings as before
						options(warn=warnVal)
						
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

##########################################################################################
# Subsampling trial functions follow:


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
#' @param q (\code{numeric)}: Subsampling level argument (mandatory). Depends on the subsampling function, it is the number of occurrences for \code{"cr"}, and the number of desired occurrences to the power of \code{x} for O^x^W. It is also the quorum of the SQS method.
#' 
#' @param dat (\code{data.frame}): Occurrence dataset, with \code{bin}, \code{tax} and \code{coll} as column names.
#' 
#' @param bin (\code{character}): The name of the subsetting variable (has to be integer). For time series, this is the time-slice variable. Rows with \code{NA} entries in this column will be omitted.
#' @param useFailed (\code{logical}): If the bin does not reach the subsampling quota, should the bin be used? 
#' @param keep (\code{numeric}): The bins, which will not be subsampled but will be added to the subsampling trials. NIf the number of occurrences does not reach the subsampling quota, by default it will not be represented in the subsampling trials. You can force their inclusion with the \code{keep} argument separetely (for all, see the \code{useFailed} argument).
#' @param unit (\code{character}): Argument of the CR subsampling type. The name of the variable that designates the subsampling units. In every bin, CR selects a certain number (quota) of entries from the dataset. By default (\code{unit=NULL}), the units will be the rows, and the \code{q} number of rows will be selected in each bin.
#' However, this can be a higher level category that has multiple entries in the each bin. If \code{unit} is a valid column of the dataset \code{dat}, then CR will select \code{q} number entries in this variable, and will return all the corresponding rows.
#' @param showFailed (\code{logical}): Toggles the output of the function. If set to \code{TRUE} the output will be a list, including both the default output (logical vector of rows) and the \code{numeric} vector of bins that did not have enough entries to reach the quota \code{q}.
#' @rdname subtrial
#' @examples
#' #one classical rarefaction trial
#'   data(corals)
#' # return 5 references for each stage
#'   bRows<-subtrialCR(corals, bin="stg", unit="reference_no", q=5)
#'   # control
#'   unCor<-unique(corals[bRows,c("stg", "reference_no")])
#'   table(unCor$stg)
#' # classical rarefaction can also be used to do by reference subsampling of collections
#'   fossils <- corals[corals$stg!=95, ]
#'   # rows in the sampling standardization
#'   threeCollPerRef <- subtrialCR(dat=fossils, unit="collection_no", 
#'     bin="reference_no", q=3, useFailed=TRUE)
#'   # check: table of collections/references in the new subset (with failed!!)
#'   crosstab<-unique(fossils[threeCollPerRef,c("collection_no", "reference_no")])
#'   table(crosstab$reference_no)
#' @export
subtrialCR<-function(dat, bin, q, unit=NULL, keep=NULL, useFailed=FALSE, showFailed=FALSE){
	quoVar <- q
	
	# omit NA bin entries
	dat <- dat[!is.na(dat[,bin]),]

	if(is.null(unit)){
		# order the rows 
		oBin <- order(dat[,bin])
		
		# the only important variable
		binVar <- dat[oBin,bin]
		
		# intact rows
		keepRows<-binVar%in%keep
		
		# output of cpp function
		binCRres<-.Call('_divDyn_CRbinwise', PACKAGE = 'divDyn', binVar, quoVar)
		
		# increase indexing (R)
		binCRres[,1]<-binCRres[,1]+1
		
		# rows where quota not reached
		failOutput <- binCRres[,1]==-8
		
		# bins where quota not reached
		fail<-binCRres[failOutput,2]
		
		# the rows where 
		usedRows<-rep(FALSE, nrow(dat))
		
		# output by subsampling
		usedRows[binCRres[!failOutput,1]] <- TRUE
		
		usedRows[keepRows] <- TRUE
		# should be kept intact
		res<-list(rows=usedRows, fail=fail)
	
				
	}else{
		if(!unit%in%colnames(dat)) stop("variable unit not found in dat.")
		
		binVar<-dat[,bin]
		unitVar <- dat[,unit]
		
		rows <- 1:nrow(dat)
		
		tap <- tapply(INDEX=binVar, X=rows, function(x){
			unitSlice<-unitVar[x]
			
			# random sample
			uniqueUnit<-sample(unique(unitVar[x]))
			
			if(quoVar<=length(uniqueUnit)){
				# rows corresponding t
				return(x[unitSlice%in%uniqueUnit[1:quoVar]])
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
		res$rows[dat[, bin]%in%res$fail] <- TRUE
	}
	if(showFailed){
		return(res)
	}else{
		return(res$rows)
	}
}


#' @param coll (\code{character}): The variable name of the collection identifiers. 
#' @param x (\code{numeric}): Argument of the OxW type. The exponent of by-list subsampling, by default it is 1.
#' @rdname subtrial
#' @export
subtrialOXW<-function(dat, bin, q, coll,x=1, keep=NULL, useFailed=FALSE, showFailed=FALSE){
	
	# omit NA bin entries
	dat <- dat[!is.na(dat[,bin]),]

	binVar <- dat[,bin]
	collVar <- dat[,coll]
	
	# future argument: appr represents how the subsampling quota is approached
	# appr can be "over", "under", "optimize"
	appr="over"
	
	# the chosen quota
	if(x==0){
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
		expTab<-table(collection)^x
		
		# shuffle
		expTab<-sample(expTab)
		
		# add up
		cumulative<-cumsum(expTab)


		# there are enough 
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
		res$rows[dat[, bin]%in%res$fail] <- TRUE
	}
	
	if(showFailed){
		return(res)
	}else{
		return(res$rows)
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
subtrialSQS <- function(dat,bin, tax, q, coll=NULL, ref=NULL, singleton="occ", excludeDominant=FALSE, 
	largestColl=FALSE, fcorr="good",byList=FALSE, keep=NULL, useFailed=FALSE, showFailed=FALSE, appr="under"){
	
	if(byList & is.null(coll)) stop("You cannot do by-list subsampling without a collection variable!")
	if(byList & !is.null(coll)){
		if(!coll%in%colnames(dat)) stop("Please provide a valid collection variable name.")
	} 
	
	# omit NA bin entries
	dat <- dat[!is.na(dat[,bin]),]

	# the quota variable
	quoVar <- q
	
	
	# the frequencies
	freqVar <- frequencies(dat=dat,bin=bin, tax=tax, coll=coll, ref=ref, 
		singleton=singleton, excludeDominant=excludeDominant, largestColl=largestColl, fcorr=fcorr)

	# the bins
	binVar=dat[,bin]

	# the collections (if present)
	collVar=dat[, coll]

	res <- subsampleSQSinexact(binVar=binVar, freqVar=freqVar, q=quoVar, collVar=collVar, byList=byList, intact=keep, appr=appr, trialRet="occ")
	
	if(useFailed){
		res$rows[dat[, bin]%in%res$fail] <- TRUE
	}
	if(showFailed){
		return(res)
	}else{
		return(res$rows)
	}
}

##########################################################################################
# Subsampling utility functions:
subsampleOXW<-function(binVar, collVar, q, intact=NULL,x=1){
	
	# future argument: appr represents how the subsampling quota is approached
	# appr can be "over", "under", "optimize"
	appr="over"
	
	# the chosen quota
	if(x==0){
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
		expTab<-table(collection)^x
		
		# shuffle
		expTab<-sample(expTab)
		
		# add up
		cumulative<-cumsum(expTab)


		# there are enough 
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


# excludeDominant sets the Alroy's 2nd correction, sets the frequency of the dominant taxon's to 0
frequencies<-function(dat,bin, tax, coll=NULL, ref=NULL, singleton="occ", excludeDominant=FALSE, largestColl=FALSE, fcorr="good"){

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
		
			# depending on whether there is a collection entry
			if(!is.null(coll)){
				occTax<-unique(dat[,c(tax, bin, coll)])
			}else{
				occTax<-dat[, c(tax, bin)]
			}
			
			
			p1<-tapply(INDEX=occTax[,bin], X=occTax[,tax], function(x){
				tTax<-table(x)
				return(sum(tTax==1))
			})
			
			p2<-tapply(INDEX=occTax[,bin], X=occTax[,tax], function(x){
				tTax<-table(x)
				return(sum(tTax==2))
			})
			
			# 'singRefTaxa' taxa that occurr only in single collections
			singRefTaxa<-tapply(INDEX=occTax[,bin], X=occTax[,tax], function(x){
				return(names(table(x)))
				
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


subsampleSQSinexact<-function(binVar, freqVar, q, collVar=NULL, byList=FALSE, intact=NULL, appr="under", trialRet="occ"){
	
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

					# if first selected species have larger then quorum frequency 
					if(appr=="under" & sum(bSelect)==0){
						bSelect[1]<-TRUE	
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
