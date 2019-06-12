#' Subsampling wrapper function
#' 
#' The function will take a function that has an occurrence dataset as an argument, and reruns it iteratively on the subsets of the dataset. 
#' 
#' The \code{subsample} function implements the iterative framework of the sampling standardization procedure. 
#' The function 1. takes the dataset \code{x}, 2. runs function \code{FUN} on the dataset and creates a container for results of trials
#' 3. runs one of the subsampling trial functions (e.g. \code{\link{subtrialCR}}) to get a subsampled 'trial dataset'
#' 4. runs \code{FUN} on the trial dataset and
#' 5. averages the results of the trials for a simple output of step 4. such as \code{vector}s, \code{matrices} and \code{data.frames}. For averaging, the \code{vectors} and \code{matrices} have to have the same output dimensions in the subsampling, as in the original object. For \code{data.frames}, the bin-specific information have to be in rows and the \code{bin} numbers have to be given in a variable \code{bin} in the output of \code{FUN}.
#' For a detailed treatment on what the function does, please see the vignette ('Handout to the R package 'divDyn' v0.5.0 for diversity dynamics from fossil occurrence data'). Currently the Classical Rarefaction (\code{"cr"}, Raup, 1975), the occurrence weighted by-list subsampling (\code{"oxw"}, Alroy et al., 2001) and the Shareholder Quorum Subsampling methods are implemented (\code{"sqs"}, Alroy, 2010).
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
#' @param x (\code{data.frame}): Occurrence dataset, with \code{bin}, \code{tax} and \code{coll} as column names.
#' @param iter (\code{numeric}): The number of iterations to be executed. 
#' @param bin (\code{character}): The name of the subsetting variable (has to be integer). For time series, this is the time-slice variable. Rows with \code{NA} entries in this column will be omitted.
#' @param tax (\code{character}): The name of the taxon variable.
#' @param useFailed (\code{logical}): If the bin does not reach the subsampling quota, should the bin be used? 
#' @param type (\code{character}): The type of subsampling to be implemented. By default this is classical rarefaction (\code{"cr"}). (\code{"oxw"}) stands for occurrence weighted by-list subsampling. If set to (\code{"sqs"}), the program will execute the shareholder quorum subsampling algorithm as it was suggested by Alroy (2010). Setting the argument to \code{"none"} will invoke no subsamling, but the applied function will be iterated on the trials, nevertheless.
#' @param FUN (\code{function}): The function to be iteratively executed on the results of the subsampling trials. If set to \code{NULL}, no function will be executed, and the subsampled datasets will be returned as a \code{list}. By default set to the \code{\link{divDyn}} function. The function must have an argument called \code{x}, that represents the dataset resulting from a subsampling trial (or the entire dataset). Arguments of the \code{subsample} function call will be searched for potential arguments of this function, which means that already provided variables (e.g. \code{bin} and \code{tax}) will also be used. You can also provide additional arguments (similarly to the \code{\link[base]{apply}} iterator). Functions that allow arguments to pass through (that have argument '...') are not allowed, as well as functions that have the same arguments as \code{subsample} but would require different values.
#' @param output (\code{character}): If the function output are vectors or matrices, the \code{"arit"} and \code{"geom"} values will trigger simple averaging with arithmetic or geometric means. If the function output of a single trial is again a \code{vector} or a \code{matrix}, setting the output to \code{"dist"} will return the calculated results of every trial, organized in a \code{list} of independent variables (e.g. if the function output is value, the return will contain a single \code{vector}, if it is a \code{vector}, the output will be a list of \code{vector}s, if the function output is a \code{data.frame}, the output will be a \code{list} of \code{matrix} class objects). If \code{output="list"}, the structure of the original function output will be retained, and the results of the individual trials will be concatenated to a \code{list}.
#' @param keep (\code{numeric}): The bins, which will not be subsampled but will be added to the subsampling trials. NIf the number of occurrences does not reach the subsampling quota, by default it will not be represented in the subsampling trials. You can force their inclusion with the \code{keep} argument separetely (for all, see the \code{useFailed} argument). 
#' @param rem  (\code{numeric}): The bins, which will be removed from the dataset before the subsampling trials.
#' @param duplicates (\code{logical} ): Toggles whether multiple entries from the same taxon (\code{"tax"}) and collection (\code{"coll"}) variables should be omitted. Useful for omitting occurrences of multiple species-level occurrences of the same genus. By default these are allowed through analyses (\code{duplicates=TRUE}), setting this to \code{FALSE} will require you to provide a collection variable. (\code{coll})
#' @param coll (\code{character}): The variable name of the collection identifiers. 
#' @param na.rm (\code{logical}): The function call includes more column names that might contain missing values. If this flag is set to \code{TRUE}, all rows will be dropped that have missig values in the specificed columns. This might lead to the exclusion of some data you do not want to exclude. 
#' @param FUN.args (\code{list}): Arguments passed to the applied function \code{FUN} but not used by the subsampling wrapper. Normally, the arguments of \code{FUN} can be added to the call of \code{subsample}, but in case you want to use different values for the same argument, then the arguments added here will be used for the call of \code{FUN}. For instance, if you want to call \code{subsample} with \code{bin=NULL}, but want to run \code{FUN=divDyn} with a valid \code{bin} column then you can add the column name here, e.g. \code{FUN.args=list(bin="stg")}.
#' @param counter (\code{logical}): Should the loop counting be visible?
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
#' \donttest{
#' # Example 2-SIB diversity 
#' # draft a simple function to calculate SIB diversity
#' sib<-function(x, bin, tax){
#'   calc<-tapply(INDEX=x[,bin], X=x[,tax], function(y){
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
#'     output="dist", keep=95, type="oxw", xexp=0)
#'   # occurrence weighted by list subsampling
#'   OW<-subsample(corals,iter=50, q=20,tax="genus", bin="stg", coll="collection_no", 
#'     output="dist", keep=95, type="oxw", xexp=1)
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
subsample<- function(x, q, tax=NULL, bin=NULL,  FUN=divDyn, coll=NULL, iter=50,  type="cr", keep=NULL, rem=NULL, duplicates=TRUE,  output="arit",  useFailed=FALSE, FUN.args=NULL, na.rm=FALSE, counter=TRUE, ...){
	
	# 0. rename arguments to make sure it is not interpreted as a closure anywhere
	quoVar <- q

	# arguments to be available later
	loop <- "for"

	# 1. SYSTEMATIC DEFENSE of main arguments
		# x
		datVect <- FALSE
		if(is.null(tax) & is.null(bin) & is.vector(x)){
			x <- data.frame(x, stringsAsFactors=FALSE)
			# keep in mind that the user will expect FUN to be applicable to a vector!
			datVect <- TRUE
			if(identical(FUN, divDyn)) stop("The default 'divDyn()' function cannot be run on vectors.")
		}

		if(!is.data.frame(x)) stop("'x' has to be a data.frame class object.")
		
		# get the call as list
		ca <- as.list(match.call())
			
		# tax
		if(any("tax"==names(ca))){
			# is it even a column?
			if(is.null(x[[tax]])) stop("You have to provide a valid taxon column.")
			if(length(tax) > 1) stop("Only one taxon column is allowed.")
			
			# logical
			bTaxNA <- is.na(x[,tax, drop=TRUE])
			if(na.rm){
				if(sum(bTaxNA)>0) x <- x[!bTaxNA,]
			}else{
				if(sum(bTaxNA)>0) stop("The 'tax' variable includes NAs.")
			}
		}

		# bin
		if(any("bin"==names(ca))){
			if(!is.null(bin)){
				if(is.null(x[[bin]])) stop("The argument 'bin' is neither a column name in 'x' nor NULL.")
				if(length(bin) > 1) stop("Only one 'bin' column is allowed.")	
	

				# logical
				bVar<- x[,bin, drop=TRUE]
				bBinNA <- is.na(bVar)
				if(na.rm){
					if(sum(bBinNA)>0) x <- x[!bBinNA,]
				}else{
					if(sum(bBinNA)>0) stop("The 'bin' variable includes NAs.")
				}

				if(sum(x[,bin, drop=TRUE]%%1)>0) stop("The bin column may only contain integers.")
			}
		}

		# q=quoVar
		if(!is.numeric(quoVar) | length(quoVar)>1) stop("Only a single numeric quota is allowed.")

		# FUN
		# check validity
		if(!is.null(FUN) & !is.function(FUN)) stop("The argument 'FUN' has to be a function or 'NULL'.")
		
		# match function
		if(!is.null(FUN)) FUN<-match.fun(FUN)
	
		# iter - complete
		if(length(iter)!=1) stop("Only a single number of iterations is allowed.")
		if(!is.numeric(iter)) stop("The entered number of iterations is not numeric.")
		if(iter<1 | iter%%1!=0) stop("Only a positive integers are allowed.")
		
		# type
		if(!type%in%c("sqs", "cr", "oxw", "none") | length(type)!=1) stop("Invalid subsampling type.")
		
		# keep- below

		# rem -below
	
		# duplicates
		if(!is.logical(duplicates)) stop("'duplicates' has to be a logical value.")
	
		# output
		if(length(output)!=1) stop("A single 'output' type is needed/allowed.")
		if(!output%in%c("arit", "geom", "list", "dist")) stop("Invalid 'output' argument.")
		
		# useFailed
		if(!is.logical(useFailed) | length(useFailed)>1) stop("'useFailed' has to be a logical value.")
		

	# 2. DATAFRAME MODIFICATION, FILTERING and further checks
		
		# duplicates - omit multiple instances of the same genus, if desired!
		if(!duplicates){
			
			# requires the coll argument
			if(!is.null(coll)){
				if(is.null(x[[coll]])) stop("The argument 'coll' has to be either 'NULL' or a column name of 'x'. ")
				if(length(coll)>1) stop("Only a single collection variable is allowed.")
			} 
			# no need to check for NAs

			# execute the function
			bDupl<-duplicated(x[, c(tax, coll)])
			x<-x[!bDupl,]
		}

		# additional columns might need to be emptied for subsampling, but the prototype should run with these. 
		# as the columns will be omitted below, save the original dataset for prototyping
		origDat <- x
	
		# interval-specific checking
		if(!is.null(bin)){
			
			# for testing and later reuse
			level<-unique(x[,bin, drop=TRUE])
			level<-level[!is.na(level)]
		
			# keep
			if(!is.null(keep)){
				if(sum(!keep%in%level)>0) stop("Invalid keep argument.")
			}
			
			# rem
			if(!is.null(rem)){
				if(sum(!rem%in%level)>0) stop("Invalid rem argument.")
				if(length(rem)>0){
					x<-x[!x[, bin, drop=TRUE]%in%rem,]
				}
			}
		}

	# 3. ARGUMENT DISTRIBUTION and SUBFUNCTION DEFENSE

		# the list of all function arguments
		passedArgs <- list(...)
		
		# additional columns that need to be checked
		check <- NULL

		if(type=="none"){
			funList <- list(FUN)
			fuNames <- "FUN"
		} 
		if(type=="cr"){
			funList <- list(FUN, checkCR)
			fuNames <- c("FUN", "checkCR")
			
			# if present, add the 'unit' column to the columns to be checked
		
			# is the unit variable NA or erroneous?
			if(any("unit"==names(ca))){

				if(!is.null(passedArgs$unit)){
					if(is.null(x[[passedArgs$unit]])) stop("The argument 'unit' has to be either 'NULL' or a column name of 'x'. ")
					if(length(passedArgs$unit)>1) stop("Only a single 'unit' variable is allowed.")

					check<- c(check, (passedArgs$unit))
				} 
			}
		
		} 
	
		if(type=="oxw"){
			funList <- list(FUN, checkOXW)
			fuNames <- c("FUN", "checkOXW")
			
			# collection is required!
			if(is.null(coll)) stop("You cannot do oxw subsampling without providing a collection 'coll' variable.")
			if(is.null(x[[coll]])) stop("The argument 'coll' has to be either 'NULL' or a column name of 'x'. ")
			if(length(coll)>1) stop("Only a single collection variable is allowed.")

			check <- c(check,coll)
			
		} 
		if(type=="sqs"){
			funList <- list(FUN, checkSQSinexact, checkFreq)
			fuNames <- c("FUN", "checkSQSinexact", "checkFreq")

			if(quoVar>1 | quoVar<=0) stop("Shareholder Quorum has to be in the 0-1 interval.")

			# by-list argument
			if(any("byList"==names(ca))){
				if(!is.logical(passedArgs$byList) | length(passedArgs$byList)!=1) stop("Invalid 'byList' argument.")
				
				# this is incomplete!
				if(passedArgs$byList){
					if(is.null(coll)) stop("The byList method requires a collection variable.")
					if(is.null(x[[coll]])) stop("The argument 'coll' has to be a column name of 'x'. ")
					if(length(coll)>1) stop("Only a single collection variable is allowed.")
					
					check <- c(check,coll)
				}
			}

			# largestColl argument
			if(any("largestColl"==names(ca))){
				if(!is.logical(passedArgs$largestColl) | length(passedArgs$largestColl)!=1) stop("Invalid 'largestColl' argument.")
				
				# this is incomplete!
				if(passedArgs$largestColl){
					if(is.null(coll)) stop("You cannot correct with the largest collection, if you do not provide collection variable.")
					if(is.null(x[[coll]])) stop("The argument 'coll' has to be a column name of 'x'. ")
					if(length(coll)>1) stop("Only a single collection variable is allowed.")
					check <- c(check,coll)
				}
			}

			if(any("singleton"==names(ca))){
				if(!is.null(passedArgs$singleton)){
					if(passedArgs$singleton=="ref"){
						if(is.null(passedArgs$ref)) stop("The 'singleton=\"ref\"' option requires a reference variable.")
						if(is.null(x[[passedArgs$ref]])) stop("The argument 'ref' has to be a column name of 'x'. ")
						if(length(passedArgs$ref)>1) stop("Only a single reference variable is allowed.")
						check <- c(check,passedArgs$ref)
					}

				}

			}

		}
	
		# check these columnns
		if(length(check)>0){
			# check the entered argument
			for(i in 1:length(check)){
				# check for NAs
				bColNA <- is.na(x[,check[i], drop=TRUE])

				# remove every NA that can be found...
				if(na.rm){
					if(sum(bColNA)>0) x <- x[!bColNA,]
				# or check whether NAs pose a threat
				}else{
					# do you have to run subsampling on the interval where the NAs are?
					# separate check when not all intervals are included
					if(!is.null(bin) & length(keep)>0){
						if(sum(is.na(x[!x[,bin, drop=TRUE]%in%keep,check[i]]))>0) 
							stop(paste("The variable '",check[i],"' includes NAs in at least one interval-subset that you want to subsample.\nOmit these, or include the interval in 'keep'. "), sep="")
					}else{
						if(sum(bColNA)>0) stop(paste("The '", check[i],"'' variable includes NAs.", sep=""))
					}
				}
			}

		}
		
		allArgs<- c(list(
			x=x, 
			quoVar=q, 
			tax=tax,
			bin=bin,  
			FUN=FUN, 
			coll=coll, 
			iter=iter, 
			type=type, 
			keep=keep, 
			rem=rem, 
			duplicates=duplicates,  
			output=output,  
			useFailed=useFailed, 
			FUN.args=FUN.args), passedArgs)

		# the distributed output
		distributedArgs <- argdist(funList, fuNames, allArgs, unused=TRUE)

		# evaluate these
		appliedArgs <- distributedArgs$FUN

		# the prototype should use the original data structure
		# in case was a vector
		if(datVect){
			appliedArgs$x <- origDat[,1, drop=TRUE]
		# or a data.frame
		}else{
			appliedArgs$x <- origDat
		}
		
		# arguments that for the applied function (FUN)
		if(is.list(FUN.args)){
			# remove more specific arguments
			appliedArgs<- appliedArgs[!names(appliedArgs)%in%names(FUN.args)]
			appliedArgs<- c(appliedArgs, FUN.args)
		}

		# check the formals of the function itself, and what is not used just state as a warning
		unUsedArgs <- distributedArgs[["0"]]
		unUsed <- names(unUsedArgs)
		remain <- unUsed[!unUsed%in%names(formals(subsample))]

		# there is an empty quote '' if you do not do this. This is something I do not understand...
	#	remain<- remain[-1]

		if(length(remain)>0) warning(paste("The argument(s) '", paste(remain, collapse="', '"), "' is/are not used with current configuration.", sep=""))

		# defend the subfunctions
			if(type=="cr") do.call("checkCR", distributedArgs$checkCR)
		#	if(type=="bs") do.call("checkBS", distributedArgs$checkBS)
			if(type=="oxw") do.call("checkOXW", distributedArgs$checkOXW)
			if(type=="sqs"){
				do.call("checkFreq", distributedArgs$checkFreq)
				do.call("checkSQSinexact", distributedArgs$checkSQSinexact)
			} 
	
	# 4. PROTOTYPING - create a replicate based on the total dataset
	# only if there is a function and if there is a matching phase at the end
	if(!is.null(FUN) & output!="list"){
		prototype<-do.call(FUN, appliedArgs)

		# criteria for FUN!!! now is the time to restrain people from using bullshit
		# matching directive, and failreplacement settings
	

		# should the bin-specific results be replaced? default:NO
		replaceFailed <- FALSE
			
		# vectors
		if(is.vector.like(prototype)){

			# try dimensions for matching if no names are present
			if(is.null(names(prototype))){
				matchdir <- "dim"
			}else{
				matchdir <- "name"
				
				
				if(!is.null(bin)){
					# if all the bins are in the names, or if they are omitted they should be omitted if failed
					if(sum(!level%in%names(prototype))==0){
						replaceFailed <- TRUE
					}
				}
			}
		}
		
		# matrix and data.frames
		if(is.data.frame(prototype) | is.matrix(prototype)){
			
			# use dimensions for matching if at least one name is missing
			if(is.null(colnames(prototype)) | is.null(rownames(prototype))){
				matchdir <- "dim"
			}else{
				matchdir <- "name"

				if(!is.null(bin)){
					# if bins represent rows
					if(sum(!level%in%rownames(prototype))==0){
							replaceFailed <- TRUE
					}
				}
			}
		}

		# special prototyping should be used for divDyn()
		if(identical(FUN, divDyn)) a<-NA

	}


	# 5. PREPARATORY CALCULATIONS and SETTINGS

	# 5A. binned subsampling
	if(!is.null(bin)){
		
		# 1. Classical Rarefaction - CR (already implemented in Rcpp)
		if(type=="cr"){
			# order the rows 
			oBin <- order(x[,bin, drop=TRUE])
			
			# reorder dataset
			x <- x[oBin,]
			
			# by-unit subsampling?
			if("unit"%in%names(distributedArgs$checkCR)){
				byUnit <- TRUE
				unitVar <- x[,distributedArgs$checkCR$unit, drop=TRUE]
				rows <- 1:nrow(x)
						
			}else{
				byUnit <- FALSE
			}

			# the only important variable
			binVar <- x[,bin, drop=TRUE]
			
			# intact rows
			keepRows<-binVar%in%keep
			
			# just organize
			oneResult<-NULL
			

		}

		# 2. Occurrence-weighted by list subsampling - OXW (not yet implemented in Rcpp)
		if(type=="oxw"){
			
			# add the arguments that were not available previously
			oxwArgs <- c(
				distributedArgs$checkOXW, 
				list(intact=keep, binVar=x[,bin, drop=TRUE],collVar=x[, coll, drop=TRUE]))
		}		

		# 3. Shareholder Quorum Subsamping - SQS (not yet implemented in Rcpp)
		if(type=="sqs"){
			# frequency calculation
			freqVar<-do.call(frequencies, distributedArgs$checkFreq)

			# by-list is by default: FALSE - to make coll necessary, you have to include it in the added arguments
			collVar <- NULL

			if(!is.null(distributedArgs$checkSQSinexact$byList)){
				if(distributedArgs$checkSQSinexact$byList){
					# and reassing the collection variable
					collVar <- x[, coll, drop=TRUE]
				}
			}
			
			# and create argument basis for SQS
			sqsArgs<-c(
				distributedArgs$checkSQSinexact,
				list(
					"binVar"=x[,bin, drop=TRUE],
					"freqVar"=freqVar,
					"collVar"=collVar,
					"intact"=keep
				)
			)

		}

		# 4. a fake subsampling trial dataset index container
		if(type=="none"){
			oneResult<- NULL
		}

		# x4. Failure-tracking setup
		failMat<-matrix(FALSE, ncol=iter, nrow=length(level))
		rownames(failMat) <- level


	# 5B. not bin-specific subsampling	
	}else{
		if(type=="cr"){
			# by-unit subsampling?
			if("unit"%in%names(distributedArgs$checkCR)){
				byUnit <- TRUE

				unitVar <- x[,distributedArgs$checkCR$unit, drop=TRUE]
				tUn <- table(unitVar)
				

				# if useFailed=FALSE, check whether the quota is low enough
				if(!useFailed & length(tUn)<quoVar) stop("The entered quota is too high. Set 'useFailed=TRUE' to get a result with the unsubsampled dataset.")
				if(length(tUn)< quoVar) quoVar <- length(tUn)


			}else{
				# if useFailed=FALSE, check whether the quota is low enoug
				if(!useFailed & nrow(x)<quoVar) stop("The entered quota is too high. Set 'useFailed=TRUE' to get a result with the unsubsampled dataset.")
				if(nrow(x)< quoVar) quoVar <- nrow(x)
				byUnit <- FALSE
				rows <- 1:nrow(x)
			}
		}
		if(type=="oxw"){
			# collecction - number present?
			if(is.null(coll))stop("You cannot do oxw subsampling without providing a collection 'coll' variable.")
			collVar <- x[[coll]]
			if(is.null(collVar)) stop("The 'coll' argument is not a variable in 'x'.")
			if(sum(is.na(collVar))>0) stop("The 'coll' variable includes NAs, please omit these manually.")

			# number of items in a collection + raise to the exponent
			collTab <- table(collVar)^distributedArgs$checkOXW$xexp
			quoVar <- quoVar
			
			# defense against high quotas
			if(!useFailed & sum(collTab)<quoVar) stop("The entered quota is too high. Set 'useFailed=TRUE' to get a result with the unsubsampled dataset.")
			if(sum(collTab)< quoVar) quoVar <- length(collTab)

		}
		if(type=="sqs"){

			# frequency calculation , but without the bin=NULL
			freqVar<-do.call(frequenciesOne, distributedArgs$checkFreq[names(distributedArgs$checkFreq)!="bin"])

			# by-list is by default: FALSE - to make coll necessary, you have to include it in the added arguments
			collVar <- NULL

			if(!is.null(distributedArgs$checkSQSinexact$byList)){
				if(distributedArgs$checkSQSinexact$byList){
					# and reassing the collection variable
					collVar <- x[, coll, drop=TRUE]
				}
			}
			
			# and create argument basis for SQS
			sqsArgs<-c(
				distributedArgs$checkSQSinexact,
				list(
					"freqVar"=freqVar,
					"collVar"=collVar
				)
			)



		}
		

	}

	# 5D. divDyn-specific optimization
	# override noNAStart, to make sure it that plotting-related empty rows do not affect performance
	if(identical(FUN, divDyn) & output!="list"){
		appliedArgs$noNAStart <- TRUE
	}

	
	# 6. ITERATION - rerun the subsampling function in loops, and run the applied function

	# 6A.  FOR-type iteration
	if(loop=="for"){
	
		# container for iteration results
		containList<-list()

		# a. bin-specific subsampling
		if(!is.null(bin)){
			
			# for every loop 
			for(k in 1:iter){

				# 1. Subsampling subphase
				# classical rarefaction
				if(type=="cr"){
					# switch, which subsampling method should be used - by-unit or by-row
					
					#by-row method
					if(!byUnit){
						# output of cpp function
						binCRres<-.Call('_divDyn_CRbinwise', PACKAGE = 'divDyn', binVar, quoVar)
					
						# increase indexing (R)
						binCRres[,1]<-binCRres[,1, drop=TRUE]+1
					
						# rows where quota not reached
						failOutput <- binCRres[,1, drop=TRUE]==-8
						
						# bins where quota not reached
						oneResult$fail<-binCRres[failOutput,2]
						
						# the rows where 
						oneResult$rows <-rep(FALSE, nrow(x))
						
						# output by subsampling
						oneResult$rows[binCRres[!failOutput,1, drop=TRUE]] <- TRUE
						
						oneResult$rows[keepRows] <- TRUE

					
					# by-unit
					}else{
						
						# subsampling a tapply() loop			
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
						oneResult<-list(rows=subsampleRows, fail=failed)
						
					} # end by-unit CR

				} # end CR

				# Occurrence-weighted by list subsampling
				if(type=="oxw"){
					
					# trial dataset rows + failed
					oneResult<-do.call(subsampleOXW, args=oxwArgs)

				} # end OXW
				
				# Shareholder-quorum subsampling
				if(type=="sqs"){
					
					# trial dataset rows + failed
					oneResult<-do.call(subsampleSQSinexact, args=sqsArgs)
					
				} # end SQS


				# 2. make the trial dataset from the subsampling trial outputs
					# subset the function-namespace x and to copy to appliedArgs
					appliedArgs$x<-x[oneResult$rows,]


				# 3. Failure registration subphase
				# add the failed stuff if so required
				if(useFailed){

					# add the failed bins' original data to the trial dataset
					appliedArgs$x<-rbind(appliedArgs$x, x[x[,bin, drop=TRUE]%in%oneResult$fail,])

				# or keep track
				}else{
					# save in the matrix
					failMat[as.character(oneResult$fail),k] <-TRUE
				}

				# 4. Application subphase

				# nothing is done
				if(is.null(FUN)){
					trialRes<-appliedArgs$x

				# apply FUN
				}else{
					trialRes<-do.call(FUN,appliedArgs)
				}

				# 5. Storage subphase
				containList[[k]]<-trialRes

				# 6. loop display
				if(counter){
					cat(k, "\r")
					utils::flush.console()
				}

			} # end for iteration

		# b. no bin subsampling
		}else{

			# for every iteration
			for(k in 1:iter){
				
				# A. CR method
				if(type=="cr"){
					# by-row CR subsampling
					if(!byUnit){
						index <- sample(rows,quoVar, replace=FALSE) 
	
					# by-unit CR subsampling
					}else{
						# randomize the unit table
						randTab <- sample(tUn)
						passed <- randTab[1:quoVar]

						# boolean index
						index <- as.character(unitVar)%in%names(passed)
					}
					

				} # end CR method

				if(type=="oxw"){
					randTab <- sample(collTab)

					# raise to exponent
					summed <- cumsum(randTab)
					bNeed<- summed <= quoVar

					index <- collVar%in%(names(randTab)[bNeed])
				}

				# Shareholder-quorum subsampling
				if(type=="sqs"){
					
					# trial dataset rows or FAILED
					index<-do.call(subsampleSQSinexactOne, args=sqsArgs)

					# NULL indicates that the trial failed
					if(is.null(index)){
						if(!useFailed){
							stop("The entered quorum is too high. Set 'useFailed=TRUE' to get a result with the unsubsampled dataset.")
						}else{
							index<- rep(TRUE, length(freqVar))
						}
					}
					
				} # end SQS

				# 2. the result of subsampling - trial dataset
				# as a vector, if it was so originally
				if(datVect){
					appliedArgs$x <- x[index,1, drop=TRUE]
				# or a data.frame
				}else{
					appliedArgs$x <- x[index,]
				}
			
				
				# 3. Application subphase
				# nothing is done
				if(is.null(FUN)){
					trialRes<-appliedArgs$x

				# apply FUN
				}else{
					trialRes<-do.call(FUN,appliedArgs)
				}

				# 4. Storage subphase
				containList[[k]]<-trialRes

				# 5. loop display
				if(counter){
					cat(k, "\r")
					utils::flush.console()
				}

			} # end iteration for loop

		}# end no bins

	} # end of for-type looping

	# 7. POST-PROCESS
	# which bins failed the trials?
	if(!is.null(bin)){
		failedBins<-sort(rownames(failMat)[which(apply(failMat,1, sum)>0)])
	}

	# 8. MATCHING subphase

	# flag, will the match succeed after going through everything?
		matchFailed <- FALSE

	# a. if output should not be a list
	# originally specified or function to be applied
	if(output!="list" & !is.null(FUN)){
		
		# choose the appropriate averagin function
		if(output=="arit") AFU <- match.fun(mean)
		if(output=="geom") AFU <- match.fun(geom)
		if(output=="dist") AFU <- NULL

		# try to match replicates to prototype, and apply the averaging function
		finalResult<- tryCatch(repmatch(containList, FUN=AFU, proto=prototype, direct=matchdir, na.rm=TRUE), error=function(a){
			message("Replicate matching failed, returning a list.")
			NULL
		})

		# if it succeeds
		if(!is.null(finalResult)){ 
			if(!is.null(bin)){
				# make the failed replacements.
				if(replaceFailed & !useFailed){
				
					# vector-specific 
					if(is.vector.like(prototype)){
						if(output!="dist"){
							finalResult[failedBins] <- NA
						}else{
							finalResult[failedBins,] <- NA
						}
					}

					# data.frame-specific things
					if(is.data.frame(prototype) | is.matrix(prototype)){
						if(output!="dist"){
							finalResult[failedBins,] <- NA
						}else{
							finalResult <- lapply(finalResult, function(w){
								w[failedBins,] <- NA
								w
							})
						}
					}
				}

				# divDyn-specific optimization
				# add the complete bin numbers back, except when distribution list!
				if(identical(divDyn, FUN) & output!="dist"){
					finalResult[,bin] <- prototype[,bin]
				}
			}

		}else {
			matchFailed <- TRUE
		}

	}

	# b. if 1.list output is asked, matching fails, or no applied FUN
	if(output=="list" | matchFailed | is.null(FUN)){
		if(!is.null(bin)){
			finalResult <- list(results=containList, failed=failedBins)
		}else{
			finalResult<-containList
		}
	}

	# 9. RETURN
	return(finalResult)
}





########################################################
# CHECKING FUNCTIONS


# These functions help the argdist() function to select the appropriate arguments
# after the arguments are separated they are called to check whether they are appropriate or not.

# only for checking!
checkCR<-function(x, quoVar, unit=NULL){
	# quota
	if(quoVar<=1 | quoVar%%1!=0) stop("The quota has to be a natural number larger than 1.")
}

#checkBS<-function(x, unit=NULL){
#	
#	if(!is.null(unit)){
#		# check whether this is present in the broader function scope 'x'
#		if(is.null(x[[unit]])) stop("The 'x' object has no column 'unit'. ")
#
#		# check if the unit column contains NAs
#		if(sum(is.na(x[,unit, drop=TRUE])>0)) stop("Variable 'unit' includes NAs, please omit these rows manually.")
#	}
#}


# for checking and distribution, additional arguments to subsampleOXW will be added within the subsample() function.
checkOXW<- function(quoVar, intact=NULL,xexp=1){
	# the quota variable
	if(quoVar<=1 | quoVar%%1!=0) stop("The quota has to be a natural number larger than 1.")
			
	# quota check
	# the xexp 
		if(!is.numeric(xexp) | length(xexp)!=1) stop("Invalid 'xexp' argument.")
		if(xexp<0) stop("Invalid 'xexp' argument.")

	# intact is an internal argument
	intact <- NULL

}

checkFreq <- function(x,bin, tax, coll=NULL, ref=NULL, singleton="occ", excludeDominant=FALSE, largestColl=FALSE, fcorr="good"){

	# ref is only checked here
	if(!is.null(ref)) if(is.null(x[[ref]])) stop("The argument 'ref' has to be either 'NULL' or a column name of 'x'. ")
	if(length(ref)>1) stop("Only a single reference variable is allowed.")

	# singleton
	if(!singleton%in%c("ref", "occ", FALSE) | length(singleton)!=1) stop("Invalid 'singleton' value.")
	if(singleton=="ref" & is.null(ref)) stop("You cannot calculate coverage based on references, if you do not provide a references variable.")

	# excludeDominant
	if(!is.logical(excludeDominant) | length(excludeDominant)!=1) stop("Invalid 'excludeDominant' argument.")
	
	if(!fcorr%in%c("good", "alroy") | length(fcorr)!=1) stop("Invalid 'fcorr' argument.")

	# things that are present but do not have to be checked - checked before
	x <- NULL
	bin <- NULL
	tax <- NULL
	
}

checkSQSinexact <- function(quoVar, byList=FALSE, intact=NULL, appr="under", trialRet="occ"){
	
	# check whether coll is existing in the broader namespace
	
	
	intact <- NULL
	appr <- NULL
	trialRet<-NULL
}


########################################################
# Type-specific functions

# excludeDominant sets the Alroy's 2nd correction, sets the frequency of the dominant taxon's to 0
frequencies<-function(x,bin, tax, coll=NULL, ref=NULL, singleton="occ", excludeDominant=FALSE, largestColl=FALSE, fcorr="good"){

	# the number of occurrences
		O<-table(x[,bin, drop=TRUE])
		rowIndex<-1:nrow(x)
	
	#calculate the frequencies
	if(!excludeDominant){
		# calculate the taxon frequencies in the bin, and assign the value to each occurrence
		freqVarList <-tapply(INDEX=x[,bin, drop=TRUE], X=rowIndex, function(w){
		#	w<-rowIndex[unique(x[,bin])[1]==x[,bin]]
			
			sliceTaxa<-x[w,tax, drop=TRUE]
			taxTab<-table(sliceTaxa)
			freq<-taxTab/length(sliceTaxa)
			
			occWise <- freq[sliceTaxa]
			names(occWise)<-w
			return(occWise)
		
		})
		
		# frequencies
		freqVar<-unlist(freqVarList, use.names=F)
		
		# row identifiers
		nameVec<-unlist(lapply(freqVarList,function(w){
			names(w)
	
		}))
		names(freqVar)<-nameVec
	
		# reorder to reflect the order of occurrences
		freqVar<-freqVar[as.character(rowIndex)]
		
	# adjust dominant to 0
	}else{
		freqVarList <-tapply(INDEX=x[,bin], X=rowIndex, function(w){
			sliceTaxa<-x[w,tax, drop=TRUE]
			taxTab<-table(sliceTaxa)
			freq<-taxTab/length(sliceTaxa)
			domIndex<-which((max(freq)==freq))
			freq[domIndex[1]]<-0
			
			occWise <- freq[sliceTaxa]
			names(occWise)<-w
			return(occWise)
		
		})
		
		freqVar<-unlist(freqVarList, use.names=F)
		nameVec<-unlist(lapply(freqVarList,function(w){
			names(w)
	
		}))
		names(freqVar)<-nameVec
	
		# reorder to reflect the order of occurrences
		freqVar<-freqVar[as.character(rowIndex)]
		
	}
	
	# Alroy 1st correction
	if(singleton%in%c("ref", "occ")){
		
		# number of sampled taxa
		vectS<-tapply(INDEX=x[,bin, drop=TRUE], X=x[,tax, drop=TRUE], function(w){
			return(length(unique(w)))	
		})
			
		# use single reference taxa
		if(singleton=="ref"){
			# count single reference taxa
			refTax<-unique(x[,c(ref,tax)])
			tTax<-table(refTax[,tax, drop=TRUE])
			
			# single reference taxa
			singRefTaxa<-names(tTax)[tTax==1]
			
			# the number of occurrences of single reference taxa
			p1<-tapply(INDEX=x[,bin, drop=TRUE], X=rowIndex, function(w){
				sliceTax<-x[w,tax, drop=TRUE]
				return(sum(sliceTax%in%singRefTaxa))
			})
			
			# two reference taxa
			twoRefTaxa<-names(tTax)[tTax==2]
			
			# the number of occurrences of single reference taxa
			p2<-tapply(INDEX=x[,bin, drop=TRUE], X=rowIndex, function(w){
				sliceTax<-x[w,tax, drop=TRUE]
				return(sum(sliceTax%in%twoRefTaxa))
			})
			
			
		}
		
		# use single occurrence taxa
		if(singleton=="occ"){
		
			# depending on whether there is a collection entry
			if(!is.null(coll)){
				occTax<-unique(x[,c(tax, bin, coll)])
			}else{
				occTax<-x[, c(tax, bin)]
			}
			
			
			p1<-tapply(INDEX=occTax[,bin, drop=TRUE], X=occTax[,tax, drop=TRUE], function(w){
				tTax<-table(w)
				return(sum(tTax==1))
			})
			
			p2<-tapply(INDEX=occTax[,bin, drop=TRUE], X=occTax[,tax, drop=TRUE], function(w){
				tTax<-table(w)
				return(sum(tTax==2))
			})
			
			# 'singRefTaxa' taxa that occurr only in single collections
			singRefTaxa<-tapply(INDEX=occTax[,bin, drop=TRUE], X=occTax[,tax, drop=TRUE], function(w){
				return(names(table(w)))
				
			})
		
		}
	
		# Alroy 2nd correction
		if(excludeDominant){
			# count the dominant taxon 
			n1 <-tapply(INDEX=x[,bin, drop=TRUE], X=rowIndex, function(w){
				sliceTaxa<-x[w,tax, drop=TRUE]
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
				taxTab<-table(x[,tax, drop=TRUE])
				
				# taxa in the most diverse collection - that are single-reference and do not occur anywhere else
				tMax<-tapply(INDEX=x[,bin, drop=TRUE], X=rowIndex, function(w){
				#	w<-rowIndex[x[,bin]==61]
				
					# unique collection/taxon table
					collTax<-unique(x[w,c(coll,tax)])
					# how many taxon/collection
					tColl<-sort(table(collTax[,coll, drop=TRUE]), decreasing=T)
					
					# which are the taxa that are there
					taxInMost<-collTax[collTax[,coll, drop=TRUE]==names(tColl[1]), tax]
					
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
			freqVar<-freqVar*uP[as.character(x[,bin, drop=TRUE])]
		}
		if(fcorr=="alroy"){
			
			freqVar <- freqVar-correction[as.character(x[,bin, drop=TRUE])]
		
		}
	}
	names(freqVar)<-x[,tax, drop=TRUE]
	return(freqVar)
}



# excludeDominant sets the Alroy's 2nd correction, sets the frequency of the dominant taxon's to 0
frequenciesOne<-function(x, tax=NULL, coll=NULL, ref=NULL, singleton="occ", excludeDominant=FALSE, largestColl=FALSE, fcorr="good"){
	# just a single variable passed
	if(ncol(x)==1) tax<-1

	# the number of occurrences
		O<-nrow(x)
		bin<-NULL
	#calculate the frequencies
	if(!excludeDominant){
		# calculate the taxon frequencies in the bin, and assign the value to each occurrence
		sliceTaxa<-x[,tax, drop=TRUE]
		taxTab<-table(sliceTaxa)
		freq<-taxTab/length(sliceTaxa)
		
		freqVar <- freq[sliceTaxa]
		
	# adjust dominant to 0
	}else{
		domIndex<-which((max(freqVar)==freqVar))
		freqVar[domIndex[1]]<-0
	}
	
	# Alroy 1st correction
	if(singleton%in%c("ref", "occ")){
		# number of sampled taxa
		vectS <- length(taxTab)
			
		# use single reference taxa
		if(singleton=="ref"){
			# count single reference taxa
			refTax<-unique(x[,c(ref,tax)])
			tTax<-table(refTax[,tax, drop=TRUE])
			
			# single reference taxa
			singRefTaxa<-names(tTax)[tTax==1]
			
			# the number of occurrences of single reference taxa
			sliceTax<-x[,tax, drop=TRUE]
			p1 <- sum(sliceTax%in%singRefTaxa)

			# two reference taxa
			twoRefTaxa<-names(tTax)[tTax==2]
			p2 <- sum(sliceTax%in%twoRefTaxa)
			
		}
		
		# use single occurrence taxa
		if(singleton=="occ"){
		
			# depending on whether there is a collection entry
			if(!is.null(coll)){
				occTax<-unique(x[,c(tax, coll)])
			}else{
				occTax<-x[, c(tax), drop=FALSE]
			}
			# table
			tTax<-table(occTax[,tax, drop=TRUE])
			p1 <- sum(tTax==1)
			p2 <- sum(tTax==2)

			# 'singRefTaxa' taxa that occurr only in single collections
			singRefTaxa<-names(table(occTax[,tax, drop=TRUE]))
		
		}
	
		# Alroy 2nd correction
		if(excludeDominant){
			# count the dominant taxon 
			sliceTaxa<-x[,tax, drop=TRUE]
			taxTab<-table(sliceTaxa)
			domIndex<-which((max(taxTab)==taxTab))
			ret<-taxTab[domIndex]
			n1<- ret[1]
			names(n1) <-NULL
			
			
			# then set 
			# for the original solution
			uP <-(O-p1-n1)/(O-n1)
			
			# for Alroy's solution
			correction<-(((p1+p2)/2)/vectS)/(O-n1)
			
			# Alroy's 3rd correction
			if(largestColl){
			
				#number of collections a taxon belong to (total dataset)
				taxTab<-table(x[,tax, drop=TRUE])
				
				# taxa in the most diverse collection - that are single-reference and do not occur anywhere else
			
					# unique collection/taxon table
					collTax<-unique(x[,c(coll,tax)])
					# how many taxon/collection
					tColl<-sort(table(collTax[,coll, drop=TRUE]), decreasing=T)
					
					# which are the taxa that are there
					taxInMost<-collTax[collTax[,coll, drop=TRUE]==names(tColl[1]), tax]
					
					# which are single reference taxa
					taxInMostAndSingRef<-taxInMost[taxInMost%in%singRefTaxa]
					bTaxInMostAndSingRef<-names(taxTab)%in%taxInMostAndSingRef
					
					# taxa only ever found in the most diverse collection
					tMax<- sum(taxTab==1 & bTaxInMostAndSingRef)
					
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
			freqVar<-freqVar*uP
		}
		if(fcorr=="alroy"){
			
			freqVar <- freqVar-correction
		
		}
	}
	names(freqVar)<-x[,tax, drop=TRUE]
	return(freqVar)
}


# Subsampling utility functions:
subsampleOXW<-function(binVar, collVar, quoVar, intact=NULL,xexp=1){
	
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

# single-bin sqs
subsampleSQSinexactOne<-function(freqVar, quoVar, collVar=NULL, byList=FALSE, intact=NULL, appr="under", trialRet="occ"){
	# non by-list SQS
	if(!byList){

		# the frequencies
		frequencies<-freqVar
		
		# the taxa
		taxa<-names(freqVar)
		
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
			if(cumulative[length(cumulative)]>quoVar){
				bSelect<-cumulative<=quoVar

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
					res<- rep(FALSE, length(freqVar))
					res[shuffledOccs[1:lastIndex]] <- TRUE
				}
				
				if(trialRet=="taxon"){
					res<- (taxa%in%keepTaxa)
				}
				
			}else{
				res<-NULL
			}
		}else{
			res<-NULL
		}
			
	# by-list SQS
	}else{
		# the frequencies
		frequencies<-freqVar
		
		# the taxa
		taxa<-names(freqVar)
		
		# collections
		colls<-collVar
		
		colTab<-table(taxa, colls)
		class(colTab) <- "matrix"
		
		# ensure that it will run if there is no collection info
		if(ncol(colTab)>1){
			
			# shuffle collections
			newOrd<-sample(1:ncol(colTab))
			colTab<-colTab[,newOrd, drop=FALSE]
			
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
				if(cumulative[length(cumulative)]>quoVar){
					bSelect<-cumulative<=quoVar
					# potential forking!!! 
					if(appr=="over"){
						if(sum(bSelect)!=length(bSelect)){
							bSelect[sum(bSelect)+1] <- TRUE
						}
					}
					
					# the collections to keep
					keepColls<-colnames(colTab)[bSelect]
					
					res<-(colls%in%keepColls)
					
					
				}else{
					res<-NULL
				}
			}else{
				res<-NULL
			}
		}
		
		
	} # end bylist
	
	# return
	return(res)

}


#- multi-bin inexact sqs
subsampleSQSinexact<-function(binVar, freqVar, quoVar, collVar=NULL, byList=FALSE, intact=NULL, appr="under", trialRet="occ"){
	
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
				if(cumulative[length(cumulative)]>quoVar){
					bSelect<-cumulative<=quoVar

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
				colTab<-colTab[,newOrd, drop=FALSE]
				
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
					if(cumulative[length(cumulative)]>quoVar){
						bSelect<-cumulative<=quoVar
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

