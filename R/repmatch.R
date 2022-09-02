#' Replicate matching and merging 
#' 
#' This pseudo-generic function iterates a function on the subelements of a list of objects that have the same class and matching dimensions/names and reorganizes the result to match the structure of the replicates or a prototype template.
#' 
#' The function is designed to unify/merge objects that result from the same function applied to different source data (e.g. the results of \code{subsample()}). In its current form, the function supports \code{vectors} (including one-dimensional \code{tables} and \code{arrays}), \code{matrix} and \code{data.frame} objects. 
#' @return If \code{FUN} is a \code{function}, the output is \code{vector} for \code{vector}-like replicates, \code{matrix} when \code{x} is a \code{list} of \code{matrix} objects, and \code{data.frame}s for \code{data.frame} replicates. In case \code{FUN=NULL}: if \code{x} is a list of \code{vectors}, the function will return a \code{matrix}; an \code{array} is returned, if \code{x} is a \code{list} of \code{matrix} class obejcts; if \code{x} is a list of \code{data.frame} objects, the function returns a \code{data.frame}.
#' @param x (\code{list}): A \code{list} of replicates. 
#' @param FUN (\code{function}): A function to merge with and to be applied to the values of identical positions in different replicates. This function must have a single output value, \code{vectors} are not allowed. The default \code{NULL} option returns an element-wise reorganization of the data.
#' @param proto (\code{same as x[[1]]}): The prototype for matching/merging. The prototype is used as a check (\code{"dim"}) or a template (\code{"name"}) during the matching process, depending on the used directive (\code{direct} argument).  It is an object with the same class as the replicates, and have the same dimensions and/or overlapping names. If the \code{"name"} directive is used and a \code{prototype} is provided, the funtion will force the output to have the same structure as the prototype, by omitting unnecessary information and inserting missing values (\code{NAs}). The prototype is expected to be an object that has more or equal elements than the replicates, otherwise the call will result in a warning.
#' @param direct (\code{character}): Matching directive(s). Can either be dimension-based (\code{"dim"}) and/or name-based (\code{name}). Dimension-based directive matches the replicates if they have the same dimesions. The \code{"name"} directive requires named input (for \code{matrices} and \code{data.frames} \code{colnames} and \code{rownames} attributes). Replicates will be matched if the values have the same names. In case both directives are specified (default), dimension-based directive takes higher priority, if matching is unsuccessful with dimensions, names will be tried after.
#' @param ... arguments passed to \code{FUN}. 
#, 
#' @examples
#' # basic example
#' vect <- rnorm(100)
#' # make 50 replicates
#' repl <- rep(list(vect), 50)
#' repmatch(repl, FUN=mean, direct="dim")
#'
#' # named input
#'   # two vectors
#'     # a
#'     a<- 1:10
#'     names(a) <- letters[1:length(a)]
#'     a[c(3,5,8)] <- NA
#'     a <- a[!is.na(a)]
#'   
#'     #b
#'     b<- 10:1
#'     names(b) <- letters[length(b):1]
#'     b[c(1, 3,6, length(b))]<- NA
#'     b <- b[!is.na(b)]
#' 
#'   # list
#'   x2 <- rep(c(list(a),list(b)), 3)
#' 
#' # simple match - falling through "dim" to "name" directive
#' repmatch(x2, FUN=NULL)
#' 
#' # prototyped
#' prot <- 1:10
#' names(prot) <-letters[1:10]
#' 
#' repmatch(x2, FUN=mean, proto=prot, na.rm=TRUE)
#' @export 
repmatch <- function(x, FUN=NULL, proto=NULL,direct=c("dim", "name"), ...){
	# x should be a list
	if(!is.list(x)) stop("'x' must be a list class object. ")
	allClass<- sapply(x, function(y){
		paste(class(y), collapse="_")
	})
	if(length(unique(allClass))!=1) stop("The replicates do not have the same class.")

	# the class of the first element
	cl<- class(x[[1]])

	# vector-like output
	res <- NULL
	if(is.vector.like(x[[1]])) res <- vrm(x=x, FUN=FUN, proto=proto, direct=direct,...)

	# matrix 
	if("matrix"%in%cl) res <- mrm(x=x, FUN=FUN, proto=proto, direct=direct,...)
	
	# i.e.multivariate 
	if("data.frame"%in%cl) res <- dfrm(x=x, FUN=FUN, proto=proto, direct=direct, ...)

	return(res)
}
 

# for vectors
vrm<-function(x,  FUN=NULL, proto=NULL,direct=c("dim", "name"), ...){
	# make sure that proto's class is the same as that of the other
	if(!is.null(proto)){
		if(!is.vector.like(proto)) stop("The prototype is not a vector.")
	}

	# assume no match will happen
	match <- FALSE

	# 1. first matching directive: the proto dimensions match those of all the replicates
	# replicate lengths are the same
	if("dim"%in%direct){
		# get the length of all replicates
		allLength <- unique(unlist(lapply(x, length)))

		# only if there is only one length type
		if(length(allLength)==1){

			# simply calculate it
			unrolled<-unlist(x)
			mat <-matrix(unrolled, ncol=length(x))
			if(!is.null(FUN)){
				res <-apply(mat, 1, FUN,...)
			}else{
				res <- mat

			}
			
			# the proto is only used for checking 
			if(!is.null(proto)){
				if(allLength==length(proto)){	
					match <- TRUE
				}else{
					match <- FALSE
				}
			}else{
				match <- TRUE
			}
		}
	}

	# 2. directive: use the names of the vectors for matching - lengths do not need to be the same
	if("name"%in%direct & !match){
		# all replicates should have names - otherwise this won't work
		
		repNames <- lapply(x, function(i){
			is.null(names(i))
		})

		# check whether all replicates have names
		if(sum(unlist(repNames))==0){
			if(!is.null(FUN)){

				# unfold the list
				all <- unlist(x)

				# some funcitons expect NAs and treat them differently
				# append NAs to the set, if the names are incomplete
				nameTab <- table(names(all))

				# where this is the maximum potential amount,that means no NA should be added
				bComplete <- nameTab==length(x)

				append <- rep(NA, sum(!bComplete))
				names(append) <- names(nameTab)[!bComplete]

				# add to the rest
				all<-c(all, append)

				# calculate FUN for all levels
				res <- tapply(INDEX=names(all), X=all,FUN=FUN,...)

				# get rid of the dimesnsion crap of the tapply()
				nameRes <- names(res)
				dim(res) <- NULL
				names(res) <-nameRes


				# if there is a proto
				if(!is.null(proto)){

					# only apply if proto has names and 
					if(is.null(names(proto))) stop("Prototype provided, but it has no names for the name-directive!")

					# compare the proto names with the replicates
						nameMissing <- !names(res)%in%names(proto)
						if(sum(nameMissing)==length(nameMissing)) stop("Prototype provided, but it has no shared names with the replicates!")

						# warning if the names of the replicates are among the names of the proto
						if(any(nameMissing)) warning("Replicates have names that the prototype doesn't have!")


					# then force proto structure
					res <- res[names(proto)]

					# cover the ugly, potential NAs in the names
					names(res) <- names(proto)
					match <- TRUE
					
				}else{
					match <- TRUE
				}
			}else{
				nameVect <- unlist(lapply(x, function(y){
					names(y)
				}))

				# unique identities
				newNames <- sort(unique(nameVect))

				# the number of entrien in every element
				noEntries <- sapply(x, length)

				# the unfolded list
				mess <- unlist(x)

				# prototype check for output
				if(!is.null(proto)){
					
					# compare the proto names with the replicates
					nameMissing <- !newNames%in%names(proto)
					if(sum(nameMissing)==length(nameMissing)) stop("Prototype provided, but it has no shared names with the replicates!")

					# warning if some of the names of the replicates are not among the names of the prototype
					if(any(nameMissing)) warning("Replicates have names that the prototype doesn't have!")

					# new names to be used
					newNames<- names(proto)

					# short mess - changed according to the prototype
					bMessKeep <- names(mess)%in%newNames
					mess <- mess[bMessKeep]


				# no prototype - all entries should be there
				}else{

					bMessKeep<-rep(T, length(mess))
				}

				# final, desired structure
				resIndex <- paste(rep(newNames, length(x)), rep(1:length(x),each=length(newNames)), sep="_")

				# index
				ind <- unlist(mapply(function(a, b){
					rep(a, b)
				}, 1:length(noEntries), noEntries))

				# rewrite the names of the object
				names(mess) <- paste(names(mess), ind[bMessKeep], sep="_")

				# the final matrix
				res<- matrix(mess[resIndex], ncol=length(x), byrow=FALSE)
				rownames(res) <- newNames

				match <- TRUE
			}
		}else{
			stop("The 'name' directive is used, but at least one replicate has no names.")
		}
	}

	# if still not found, report
	if(!match & !is.null(proto)) stop("Could not match vector replicates with the prototype.")
	if(!match & is.null(proto)) stop("Dimensions of vector replicates do not match.")

	return(res)

}



# for matrices  
mrm<-function(x,  FUN=NULL, proto=NULL,direct=c("dim", "name"), ...){
	# make sure that proto's class is the same as that of the other
	if(!is.null(proto)){
		if(!is.matrix(proto)) stop("The prototype is not a matrix.")
	}

	# assume no match will happen
	match <- FALSE

	# 1. first matching directive: the proto dimensions match those of all the replicates
	# replicate lengths are the same
	if("dim"%in%direct){
		# get the length of all replicates
		allDim <- t(sapply(x, dim))
		dimGuide <- unique(unique(allDim))

		# only if there is only one dimensiony type
		if(nrow(dimGuide)==1){

			# simply calculate it
			unrolled<-unlist(x)
			arr <-array(unrolled, dim=c(dimGuide[1,,drop=TRUE], length(x)))
			
			if(!is.null(FUN)){
				res <-apply(arr, c(1,2), FUN,...)
			}else{
				res <- arr

			}
			
			# the proto is only used for checking 
			if(!is.null(proto)){
				if(sum(dimGuide==dim(proto))==length(dimGuide)){	
					match <- TRUE
				}else{
					match <- FALSE
				}
			}else{
				match <- TRUE
			}
		}
	}

	# 2. directive: use the names of the vectors for matching - lengths do not need to be the same
	if("name"%in%direct & !match){
		# all replicates should have names - otherwise this won't work
		
		repNames <- lapply(x, function(i){
			is.null(colnames(i)) | is.null(rownames(i))
		})

		# check whether all replicates have names
		if(sum(unlist(repNames))==0){
			
			# get the rownames
			allRowNames<- unique(unlist(lapply(x, function(i){
				rownames(i)
			})))

			# get the colnames
			allColNames<- unique(unlist(lapply(x, function(j){
				colnames(j)
			})))

			# is there a prototype?
			if(!is.null(proto)){
				# is it useable?
				if(is.null(rownames(proto)) | is.null(colnames(proto))) stop("Prototype provided, but it has no 'colnames' and 'rownames' to be used with the name directive.")

				# overlap with the names of replicates?
				if(sum(rownames(proto)%in%allRowNames)==0) stop("The prototype and the replicates do not share any element in their 'rownames' attribute.")
				if(sum(colnames(proto)%in%allColNames)==0) stop("The prototype and the replicates do not share any element in their 'colnames' attribute.")

				if(sum(allColNames%in%colnames(proto))!=length(allColNames)) warning("The replicates have 'colnames' entries that the prototype doiesn't have.")
				if(sum(allRowNames%in%rownames(proto))!=length(allRowNames)) warning("The replicates have 'rownames' entries that the prototype doiesn't have.")

				# the use the prototype's names for the output
				allRowNames <- rownames(proto)
				allColNames <- colnames(proto)
			}
		
			# create an array for holding the data
			arr <- array(NA, dim=c(length(allRowNames), length(allColNames), length(x)))
			dimnames(arr) <- list(allRowNames,allColNames,1:length(x))

			# for starters, use a simple for loop to fill in the matrix
			for(i in 1:length(x)){
				cur <- x[[i]]
				matchRow <-  rownames(cur)[rownames(cur)%in%allRowNames]
				matchCol <- colnames(cur)[colnames(cur)%in%allColNames]

				arr[matchRow, matchCol,i] <- cur[matchRow, matchCol]
			}

			# then appply function if necessary
			if(!is.null(FUN)){
				res <-apply(arr, c(1,2), FUN,...)
			}else{
				res <- arr
			}
			
			match <- TRUE
		
			
		}else{
			stop("The 'name' directive is used, but at least one replicate has no names.")
		}
	}

	# if still not found, report
	if(!match & !is.null(proto)) stop("Could not match vector replicates with the prototype.")
	if(!match & is.null(proto)) stop("Dimensions of vector replicates do not match.")

	return(res)

}







# @param to (\code{character}) Argument of the \code{data.frame} matching method. Either \code{"num"}, \code{"log"} or \code{"char"}. This argument specifices which variable types the \code{FUN} will be applied to (to avoid problems with non-applicables, such \code{mean} to \code{character} \code{vectors}).
dfrm <- function(x, FUN=NULL, proto=NULL, direct=c("dim", "name"), to=c("num", "log"),...){
	if(!is.null(proto)){
		if(!class(proto)%in%"data.frame") stop("Prototype class is not data.frame!")
	}
	
	# assume a non-successful match
	match <- FALSE
	
	# 1.st directive: use the dimensions and types
	if("dim"%in%direct){
		# immediately check whether it is applicable or not - <make sure every rep has two dimensions!>
		allDim <- t(sapply(x, dim))
		dimGuide <- unique(unique(allDim))

		# if all have just one dimension vector
		if(nrow(dimGuide)==1){
			dim(dimGuide) <- NULL
			
			if(!is.null(proto)){
				if(sum(dim(proto)==dimGuide)!=length(dimGuide)) stop("The dimensions of prototype doesn't match that of the replicates.")
			}

			# check the vector types in each column of every replicate
			types<- t(sapply(x, function(i){
				sapply(i,function(j){
					if(is.numeric(j)) ret <- "num"
					if(is.character(j)) ret <- "char"
					if(is.logical(j)) ret <- "log"
					if(is.factor(j)) ret <- "fac"
					return(ret)
				})
				
			}))

			# are they consistent?
			typeGuide <- unique(types)

			origColNames <- colnames(typeGuide)


			# if the types also match 
			if(nrow(typeGuide)==1){
				dim(typeGuide)<- NULL
				if("fac"%in%typeGuide) stop("Factor variables are not supported, use 'character' instead.")

				# at this point assume that the match succeeded
				match <- TRUE

				# is there a prototype -  check whether it matches the replicates
				if(!is.null(proto)){

					# make sure that the prototype has the same types as the replicates
					typeProto<- sapply(proto,function(j){
						ret <- NA
						if(is.numeric(j)) ret <- "num"
						if(is.character(j)) ret <- "char"
						if(is.logical(j)) ret <- "log"
						if(is.factor(j)) ret <- "fac"
						return(ret)
					})
					# and check whether they match
					if(any(typeProto[origColNames]!=typeGuide)){
						# match is unsuccesful, you have to use the name directive
						match <- FALSE
					}

				#if there isn't, then use replicates' , it's cool 
				}

				# continue if not directed otherwise
				if(match){

					# unroll the whole thing:
					mess <- unlist(x)

					# three index vectors to separate the numeric, character and logical entries
					numInd <- which(typeGuide=="num")
					charInd <- which(typeGuide=="char")
					logInd <- which(typeGuide=="log")

					# boolean index values for these types  
					bNum <-selcolUnrolled(colInd=numInd, nRow=dimGuide[1], nCol=dimGuide[2], nPlane=length(x))
					bChar <-selcolUnrolled(colInd=charInd, nRow=dimGuide[1], nCol=dimGuide[2], nPlane=length(x))
					bLog <-selcolUnrolled(colInd=logInd, nRow=dimGuide[1], nCol=dimGuide[2], nPlane=length(x))
					
					#select the values and cast
					numPart <- as.numeric(mess[bNum])
					charPart <- as.character(mess[bChar])
					logPart <- as.logical(mess[bLog])

					# structure these to arrays
					arNum <- array(numPart, dim=c(dimGuide[1],length(numInd), length(x)))
					arChar <- array(charPart, dim=c(dimGuide[1],length(charInd), length(x)))
					arLog <- array(logPart, dim=c(dimGuide[1],length(logInd), length(x)))
					
					# regular function applied
					if(!is.null(FUN)){
						# and apply if asked for - default is not for characters!!!
						if("num"%in%to){
							resNum <- as.data.frame(apply(arNum, c(1,2), FUN,...))
						}else{
							resNum <- as.data.frame(matrix(as.numeric(NA), nrow=dimGuide[1], ncol=length(numInd)))
						}
						if("char"%in%to){
							resChar <- as.data.frame(apply(arChar, c(1,2), FUN,...), stringsAsFactors=FALSE)
						}else{
							resChar <- as.data.frame(matrix(as.character(NA), nrow=dimGuide[1], ncol=length(charInd)), stringsAsFactors=FALSE)
						}

						if("log"%in%to){
							resLog <- as.data.frame(apply(arLog, c(1,2), FUN,...))
						}else{
							resLog <- as.data.frame(matrix(as.logical(NA), nrow=dimGuide[1], ncol=length(logInd)))
						}
						
						# create an output dataframe
						res<- cbind(resChar, resLog, resNum)

						# since the dimensions did not change, you can just copy the row names
						rownames(res) <- rownames(x[[1]])

						
					# return variable-list construct!
					}else{
						# create a list for all Variables (columns)
						# numeric first
						numList <- list()
						for(i in 1:length(numInd)){
							numList[[i]]<- arNum[,i,]
							rownames(numList[[i]]) <- rownames(x[[1]])
						}
						charList <- list()
						for(i in 1:length(charInd)){
							charList[[i]]<- arChar[,i,]
							rownames(charList[[i]]) <- rownames(x[[1]])
						}
						logList <- list()
						for(i in 1:length(logInd)){
							logList[[i]]<- arLog[,i,]
							rownames(logList[[i]]) <- rownames(x[[1]])
						}

						# concatenate
						res <- c(charList, logList, numList)

					}
					# use previous indices to 
					origInd <- c(charInd, logInd, numInd)
					resInd <- rep(NA, length(origInd))
					resInd[origInd] <- 1:length(resInd)

					# order the types correctly
					res <- res[resInd]

					# use the names of the first names
					names(res) <- origColNames
				}

			}else{
				# error:dimensions match but types do not -use names!!!
			}
		}
	}

	# 2.nd directive: use the names (rownames) to match
	if("name"%in%direct & !match){
		# get the dimensions
		allDim <- t(sapply(x, dim)) # they can be different

		# if there is a prototype
		if(!is.null(proto)){
			# check whether prototype has names
			if(is.null(colnames(proto)) | is.null(rownames(proto))) stop("Prototype has no names, cannot be processed with names directive.")
		}

		# check whether everyhing has names
		repNames <- lapply(x, function(i){
			is.null(colnames(i)) | is.null(rownames(i))
		})
		
		# only if every replicate has both colnames and rownames
		if(sum(unlist(repNames))==0){
			
			# check the vector types in each column of every replicate
			types<- lapply(x, function(i){
				sapply(i,function(j){
					if(is.numeric(j)) ret <- "num"
					if(is.character(j)) ret <- "char"
					if(is.logical(j)) ret <- "log"
					if(is.factor(j)) ret <- "fac"
					return(ret)
				})
				
			})

			# unfold the types, just a vector of tpyes
			names(types) <- NULL
			messType <- unlist(types)

			# check whether the same variable names have just one type!
			typeMatch <- tapply(INDEX=names(messType), X=messType, FUN=function(i){
				length(unique(i))==1
			})

			# all types match
			if(!any(!typeMatch)){
				
				# unfold the data
				mess <- unlist(x)

				# rownames of every data frame
				na<-lapply(x, function(d){
					rownames(d)
				})

				# create a unique ID for every entry entry
				# if there is a function
				if(!is.null(FUN)){
					
					idList <- mapply(FUN=function(i,j){
					#	i<-na[[1]]
					#	j<-types[[1]]
						rowID <-rep(i, length(j))
						bothID<- paste(rowID, rep(names(j), each=length(i)), sep="_")
					}, na, types)
				
				# if there is NULL
				}else{

					idList <- mapply(FUN=function(i,j,k){
					#	i<-na[[1]]
					#	j<-types[[1]]
						rowID <-rep(i, length(j))
						planeID <- rep(k, length(rowID))
						bothID<- paste(rowID, rep(names(j), each=length(i)),planeID, sep="_")
					}, na, types, 1:length(x))
				
				}

				# the unique name of every entry
				names(mess) <- unlist(idList)

				# distribute data to numeric, logical and 
				numInd <- which(messType=="num")
				charInd <- which(messType=="char")
				logInd <- which(messType=="log")
				# factor and others?

				#how many entries are there in column
				colGuide <- apply(allDim, 1, function(d){
					rep(d[1], d[2])
				})
				if(is.list(colGuide)) colGuide <- unlist(colGuide)
				if(is.matrix(colGuide)) dim(colGuide)<-NULL
				
				bLog<-selcolDiffdat(colGuide, logInd)
				bNum<-selcolDiffdat(colGuide, numInd)
				bChar<-selcolDiffdat(colGuide, charInd)

				# get the value specific parts and cast
				messLog <-as.logical(mess[bLog])
				messNum <- as.numeric(mess[bNum])
				messChar <- as.character(mess[bChar])

				# the averaging can be done with tapply so we need the names
				nameLog <- names(mess)[bLog]
				nameNum <- names(mess)[bNum]
				nameChar <- names(mess)[bChar]

				# there is a function to process the data
				if(!is.null(FUN)){
					# and apply if asked for - default is not for characters!!!
					if("num"%in%to){
						resNum <- tapply(INDEX=nameNum, X=messNum, FUN=FUN,...)
						# omit stupid NaN
						resNum[is.nan(resNum)] <- NA
					}else{
						uNumName <- sort(unique(nameNum))
						resNum<- rep(as.numeric(NA), length(uNumName))
						names(resNum) <- uNumName
					}
					if("char"%in%to){
						resChar <- tapply(INDEX=nameChar, X=messChar, FUN=FUN,...)
					}else{
						uCharName <- sort(unique(nameChar))
						resChar<- rep(as.character(NA), length(uCharName))
						names(resChar) <- uCharName
					}
					if("log"%in%to){
						resLog <- tapply(INDEX=nameLog, X=messLog, FUN=FUN,...)
						resLog[is.nan(resLog)] <- NA
					}else{
						uLogName <- sort(unique(nameLog))
						resLog<- rep(as.logical(NA), length(uLogName))
						names(resLog) <- uLogName
					}
				
				# distribution output
				}else{
					names(messNum) <- nameNum
					names(messLog) <- nameLog
					names(messChar) <- nameChar
				}

				# final object structure
				if(is.null(proto)){
					# the final potential rownames
					totRows <- sort(unique(unlist(na)))

					# final potential variables
					totNumCols <- unique(names(messType)[numInd])
					totCharCols <- unique(names(messType)[charInd])
					totLogCols <- unique(names(messType)[logInd])

				}else{
					# column and rownames of the prototype
					totRows<- rownames(proto)

					# types of the prototype
					typeProto <- sapply(proto,function(j){
						if(is.numeric(j)) ret <- "num"
						if(is.character(j)) ret <- "char"
						if(is.logical(j)) ret <- "log"
						if(is.factor(j)) ret <- "fac"
						return(ret)
					})
					totNumCols <- names(typeProto)[typeProto=="num"]
					totCharCols <- names(typeProto)[typeProto=="char"]
					totLogCols <- names(typeProto)[typeProto=="log"]

				}

				# for function output
				if(!is.null(FUN)){
					# potential identifiers
					potNumNames <- paste(rep(totRows, length(totNumCols)), rep(totNumCols, each=length(totRows)),sep="_")
					potCharNames <- paste(rep(totRows, length(totCharCols)), rep(totCharCols, each=length(totRows)),sep="_")
					potLogNames <- paste(rep(totRows, length(totLogCols)), rep(totLogCols, each=length(totRows)),sep="_")

					# reorder the result output for folding
					resNumFold <- resNum[potNumNames]
					resCharFold <- resChar[potCharNames]
					resLogFold <- resLog[potLogNames]

					# fold and make a data.frame
					numMat <- as.data.frame(matrix(resNumFold, nrow=length(totRows), byrow=FALSE))
					charMat <- as.data.frame(matrix(resCharFold, nrow=length(totRows), byrow=FALSE), stringsAsFactors=FALSE)
					logMat <- as.data.frame(matrix(resLogFold, nrow=length(totRows), byrow=FALSE))

					# assign the columns
					colnames(numMat) <- totNumCols
					colnames(charMat) <- totCharCols
					colnames(logMat) <- totLogCols

					# output and rownames
					res <- cbind(numMat, charMat, logMat)
					rownames(res) <- totRows

				# dist output
				}else{
					# potential identifiers
					potNumNamesPlane <- paste(rep(totRows, length(totNumCols)), rep(totNumCols, each=length(totRows)),sep="_")
					potCharNamesPlane <- paste(rep(totRows, length(totCharCols)), rep(totCharCols, each=length(totRows)),sep="_")
					potLogNamesPlane <- paste(rep(totRows, length(totLogCols)), rep(totLogCols, each=length(totRows)),sep="_")

					# with the plane identifiers
					potNumNames <- paste(rep(potNumNamesPlane, length(x)), rep(1:length(x), each=length(potNumNamesPlane)), sep="_")
					potCharNames <- paste(rep(potCharNamesPlane, length(x)), rep(1:length(x), each=length(potCharNamesPlane)), sep="_")
					potLogNames <- paste(rep(potLogNamesPlane, length(x)), rep(1:length(x), each=length(potLogNamesPlane)), sep="_")

					# reorder input
					resNumFold <- messNum[potNumNames]
					resCharFold <- messChar[potCharNames]
					resLogFold <- messLog[potLogNames]

					# array representation
					numArr<- array(resNumFold, dim=c(length(totRows), length(totNumCols), length(x)))
					logArr<- array(resLogFold, dim=c(length(totRows), length(totLogCols), length(x)))
					charArr<- array(resCharFold, dim=c(length(totRows), length(totCharCols), length(x)))

					# make the variable specific data in a list
					# slowest solution ever...
					numList<-list()
					if(length(totNumCols)>0){
						for(i in 1:length(totNumCols)){
							temp <- numArr[,i,]
							rownames(temp) <- totRows
							numList[[i]] <- temp
						}
					}

					logList<-list()
					if(length(totLogCols)>0){
						for(i in 1:length(totLogCols)){
							temp <- logArr[,i,]
							rownames(temp) <- totRows
							logList[[i]] <- temp
						}
					}

					charList<-list()
					if(length(totCharCols)>0){
						for(i in 1:length(totCharCols)){
							temp <- charArr[,i,]
							rownames(temp) <- totRows
							charList[[i]] <- temp
						}
					}

					res<-c(numList,logList,charList)
					names(res) <- c(totNumCols,totLogCols,totCharCols)

					# if there is a prototype
					if(!is.null(proto)){
						res <- res[names(typeProto)]
						names(res) <- names(typeProto)
					}
					
				}
			
			

			} # else fail!
		
		}else{
			stop("matching failed.")

		}

	}
	
	return(res)
}


###### utility functions

# select columns from an unrolled 3d array
selcolUnrolled <- function(colInd,nRow, nCol, nPlane){
	
	#colInd <- numInd
	#nRow <- dimGuide[1]
	#nCol<- dimGuide[2]
	#nPlane<- length(x)

	# start with the numeric indices
	totNum <- rep(FALSE, nRow*nCol*nPlane)
	
	if(length(colInd)>0){
		# for all iterations
		for(i in 1:nPlane){
			# for of the indices
			for(j in 1:length(colInd)){
				# all entries from the column
				entries<- 1:nRow
				# offset the columns
				colOffset <- nRow*(colInd[j]-1)
				
				# offset the tables: shift columns by this amount time the number entries in a column
				iterOffset <- (i-1)*nCol*nRow

				totNum[entries+colOffset+iterOffset] <- TRUE
			}
		}
	}
	return(totNum)
}


# select column from unrolled data.frames with 
selcolDiffdat <- function(colGuide, colInd){
	# logical vector for subsetting
	bSubset <- rep(FALSE, sum(colGuide))
	if(length(colInd)>0){
		# for every relevant column
		for(i in 1:length(colInd)){
			if(colInd[i]==1){
				offset<- 0
			}else{
				offset<- sum(colGuide[1:(colInd[i]-1)])
			}
			bSubset[(1:colGuide[colInd[i]])+offset] <- TRUE
		}
	}
	return(bSubset)
}

is.vector.like <- function(x){
	cl <- class(x)
	if(
		"numeric"%in%cl | 
		"integer"%in%cl | 
		"double"%in%cl | 
		"character"%in%cl | 
		"logical"%in%cl | 
		("table"%in%cl & length(dim(x))==1) |
		("array"%in%cl & length(dim(x))==1)
	) state <- TRUE else state <- FALSE
	return(state)
}