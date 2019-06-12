#' Apply function to TAxon/BIN subset of occurrences and iterATE
#'
#' The function takes another function and reruns it on every taxon- and/or bin-specific subsets of an occurrence dataset.
#'
#' The main \code{tabinate} function acts as a wrapper for any type of function that requires a subset of the occurrence dataset that represents either one \code{bin} or one \code{tax} entry or both.
#' For example, the iterator can be used to calculate geographic ranges from occurrence coordinates (\code{georange}).
#'
#' The output structure of FUN should be independent from the input subset, or the function will return an error. 
#' Setting both \code{bin} If \code{bin=NULL} and code{tax=NULL}, will run \code{FUN} on the entire dataset (no effect). Providing either \code{bin} or \code{tax} and keeping the other \code{NULL} will iterate \code{FUN} for every \code{bin} or \code{tax} entry (whichever is presented).
#' The function returns a vector of values if the return value of \code{FUN} is a single value. In case it is a vector, the final output will be a matrix. 
#' When both \code{bin} and \code{tax} is presented, the function output will be a matrix (one output value for a taxon/bin subset) or an array (3d, when \code{FUN} returns a vector).  Setting \code{FUN} to \code{NULL} will return the occurrence dataset as \code{list}s. 
#' 
#' @param x \code{(data.frame)} Fossil occurrence table.
#' 
#' @param FUN (\code{function}) The function applied to the subset of occurrences. The subset of occurence data will be passed to this function as \code{x}. 
#' 
#' @param bin \code{(character)} Variable name of the bin numbers of the occurrences. This variable should be \code{numeric}. 
#' 
#' @param tax \code{(character)} Variable name of the occurring taxa (variable type: \code{factor} or \code{character} - such as \code{"genus"}
#' 
#' @param ... arguments passed to \code{FUN} 
#' 
#' @rdname tabinate
#'  
#' @examples
#'  data(corals)
#'
#' # the number of different coordinate pairs in every time slice
#'   tabinate(corals, bin="stg", FUN=georange, lat="paleolat", 
#'     lng="paleolng", method="co")
#' # geographic range (site occupancy) of every taxon in every bin
#'   tabinate(corals, bin="stg", tax="genus", FUN=georange, 
#'     lat="paleolat", lng="paleolng", method="co")
#'
#' @export
tabinate <- function(x,bin=NULL, tax=NULL, FUN=NULL, ...){
#	x <- testData
#	bin <- binName
#	tax <- taxName
#	FUN <- function(x) c(nrow(x), nrow(x)+1)
#	FUN <- function(x) {one <- c(nrow(x), nrow(x)+1); names(one)<-c("a", "b"); one}
#	addArgs<-list()

	# additional arguments
	addArgs <- list(...)

	# all the different occurrences
	if(is.null(bin) & is.null(tax)){
	
		# call function based on this
		callArgs <-list(
			x=x
			)
		
		# append user-supplied arguments to those defined by this function
		callArgs <- c(callArgs, addArgs)
		
		if(!is.null(FUN)){
			# call applied function
			oneResult<- do.call(FUN, callArgs)
		}else{
			oneResult <- x
		}

		# the final result
		res <- oneResult

	# should be iterated across multiple bins or taxa
	}else{
		# single-binned analysis with multiple taxa
		if(is.null(bin)){
			# omit those entries where no taxon is given
			x<-x[!is.na(x[, tax, drop=TRUE]),]

 			# if factor convert to character
 			if(is.factor(x[, tax, drop=TRUE])) x[, tax, drop=TRUE] <- as.character(x[, tax, drop=TRUE])

			rows <- 1:nrow(x)
			iteratorListOutput<-tapply(INDEX=x[,tax, drop=TRUE], X=rows, function(w){
				callArgs <-list(
					x=x[w,]
				)

				callArgs <- c(callArgs, addArgs)
				if(!is.null(FUN)){
					# call the applied function
					oneResult<- do.call(FUN, callArgs)
				}else{
					oneResult <- x[w,]
				}
				
				return(oneResult)
			})

			res <- flattenList(iteratorListOutput)
	
		# multi-bin application
		}else{

			# multi-bin, one taxon
			if(is.null(tax)){
				# recursion, setting tax, to bin!
				callArgs <- list(
					x=x,
					bin=NULL, 
					tax=bin,
					FUN=FUN
				)
				callArgs<-c(callArgs, addArgs)
				res<- do.call(tabinate, callArgs)
				
			# multi-bin, multi taxon
			}else{
				
				# omit values where bin or x is missing
				x<-x[!is.na(x[,bin, drop=TRUE]) & !is.na(x[,tax, drop=TRUE]),]

				rows <- 1:nrow(x)
				# on every bin
				tabin <- paste(x[, tax, drop=TRUE], x[,bin, drop=TRUE], sep="_")

				iterRes<- tapply(INDEX=tabin, X=rows, FUN=function(w){
					
					if(!is.null(FUN)){
						callArgs <-list(
							x=x[w,]
						)
						callArgs <- c(callArgs, addArgs)
					
						# call the applied function
						return(do.call(FUN, callArgs))
					}else{
						return(x[w,])
					}
					
				})

				if(is.null(FUN)){
					res<-iterRes
				}


				# simple value output
				if(is.numeric(iterRes) | is.logical(iterRes) | is.character(iterRes)){
					allTax <- sort(unique(x[, tax, drop=TRUE]))
					allBin <- sort(unique(x[, bin, drop=TRUE]))
	
					taxpart <- rep(allTax, each=length(allBin))
					binpart <-rep(allBin, length(allTax))
					allNames<-paste(taxpart, binpart, sep="_")
	
					theMat <- rep(NA, length(allNames))
					names(theMat) <-allNames
					theMat[names(iterRes)]<-iterRes
	
					# proper dimensions
					theMat <- matrix(theMat, ncol=length(allTax),nrow=length(allBin),  byrow=FALSE)
					colnames(theMat) <- allTax
					rownames(theMat) <- allBin

					res<- theMat
				}

				# if it is a list (more complex output)
				if(is.list(iterRes)){
					# if the element is a list
					if(is.numeric(iterRes[[1]]) | is.logical(iterRes[[1]]) | is.character(iterRes[[1]])){
						vecLen <- unlist(lapply(iterRes, length))
	
						# simple array - all vector lengths are the same
						diffLen<-unique(vecLen)
						if(length(diffLen)==1){
							# length of the vectors
							theLen <- diffLen[1]
							theNames<-names(iterRes[[1]])
						}else{
							theLen<- max(vecLen)
							theNames <- names(iterRes[[which(vecLen==max(vecLen))[1]]])
						}
							
						if(is.null(theNames)) theNames<-1:theLen

						# all possible entries are combined from these
						allTax <- sort(unique(x[, tax, drop=TRUE]))
						allBin <- sort(unique(x[, bin, drop=TRUE]))
	
						# the original results (flat)
						names(iterRes)<-paste(names(iterRes), ".", sep="")
						flattened <- unlist(iterRes)

						# in case there are original names of the elements (delete this)
						namesSplit <- strsplit(names(flattened),"\\.")
						simplifiedNames <- unlist(lapply(namesSplit, function(w) w[[1]]))

						tempor<- unlist(sapply(vecLen, function(w) 1:w))

						names(flattened) <- paste(simplifiedNames, tempor, sep="_")

						# the name of the final container (flat)
						taxpart <- rep(allTax, each=length(allBin))
						binpart <-rep(allBin, length(allTax))

						# first two dimensions
						allNames<-paste(taxpart, "_",binpart,"_", sep="")
						
						# extend to three dimensions
							nameIndex<-rep(allNames, theLen)
							second<-rep(1:theLen, each=length(allNames))

						# the actual 3d flat container
							theArray<-rep(NA, length(nameIndex))
							names(theArray)<-paste(nameIndex, second, sep="")

						# fill in the container
						theArray[names(flattened)]<-flattened

						# make it 3d dimensional
						resArray<-array(theArray, dim=c(length(allBin), length(allTax),theLen))
						dimnames(resArray) <- list(allBin, allTax, theNames)
						# the final output
						res<-resArray
												
					}
				}

			} # end of multi-bin method
		}
	}
	return(res)
}


# the output of the tapply() loop is a list
flattenList <- function(iteratorListOutput){
	if(is.list(iteratorListOutput)){
		
		# first element is a vector
		singleElement <- iteratorListOutput[[1]]
		listOut <- TRUE
		# matrix final output
		if(is.numeric(singleElement) | is.character(singleElement) | is.logical(singleElement)){
			# are all the vectors of the same length?
			varlength<-lapply(iteratorListOutput, FUN=length)

			# all of them are of the same length
			if(length(unique(varlength))==1){
				listOut <- FALSE
				cMethods <- varlength[[1]]
				res <- matrix(NA, nrow=length(iteratorListOutput), ncol=cMethods)
				for(i in 1:cMethods){
					res[,i]<-unlist(lapply(iteratorListOutput, function(w) w[i]))
				}
				rownames(res) <- names(iteratorListOutput)
				colnames(res) <- names(iteratorListOutput[[1]])
			}
		}
		if(is.matrix(singleElement)){
			listOut <- FALSE
			res <- "NOT YET!"

		}

		# list final output
		if(listOut){
			res<- iteratorListOutput
		}
	# vector final output
	}else{
		res <- iteratorListOutput
		res<-as.numeric(iteratorListOutput)
		names(res)<- names(iteratorListOutput)
	}
	return(res)
}

# oldmethod
# outerRes <- tapply(INDEX=x[,bin], X=rows, function(w){
# 				#	w<-rows[x[,bin]==80]
# 
# 					binTax <- x[w,tax]
# 				
# 					# result for every taxon in a current bin - returns a list
# 					innerRes <- tapply(INDEX=binTax, X=w, function(y){
# 					#	y<- w[binTax==binTax[1]]
# 						callArgs <-list(
# 							x=x[y,, drop=FALSE]
# 						)
# 
# 						callArgs <- c(callArgs, addArgs)
# 						
# 						if(!is.null(FUN)){
# 							# call range calculation function
# 							oneResult<- do.call(FUN, callArgs)
# 						}else{
# 							
# 							oneResult <- x[y,]
# 						}
# 
# 						return(oneResult)
# 						
# 					})
# 
# 					# what is the output of a single iteration?
# 					if(is.numeric(innerRes) | is.character(innerRes) | is.logical(innerRes)){
# 						# do not do anything
# 						return(innerRes)
# 					}
# 					
# 					if(is.list(innerRes)){
# 						# flatten the output
# 						tempInner<- flattenList(innerRes)
# 						return(tempInner)
# 					}
# 				})# output is a list
