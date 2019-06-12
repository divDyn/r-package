#' List of singleton taxa
#' 
#' The function returns lists of taxa that occurr with only one particular entry in a given variable. 
#' 
#' Singletons are defined in number of ways in the literature. True singletons are species that are represented by only one specimen, but one can talk about
#' single-occurrence, single-interval, single-reference or single collection taxa as well. These can be returned with this function. 
#' 
#' As the time bin has particular importance, it is possible to filter singleton taxa in the context of a single bin. These can be returned with the \code{bybin} argument, that constrains and iterates the filtering to every bin.
#' If this argument is set to \code{TRUE} and the variable in question is a references, than single-reference taxa will be taxa that occurred in only one reference within each bin - it does not necessarily mean that only one reference describes the taxon in the total database!
#' 
#' 
#' @param dat (\code{data.frame}): Occurrence dataset, with \code{bin}, \code{tax} and \code{coll} as column names.
#' @param bin (\code{character}): Lists of taxa can be tabulated in every bin. Rows with \code{NA} entries in this column will be omitted.
#' @param tax (\code{character}): The name of the taxon variable.
#' @param var (\code{character}): The variable that is used to define singletons. Use the reference variable for single-reference taxa, and the collection variable for single-collection taxa, the bin identifier for single-interval taxa and so forth. If you set this to the default \code{NULL}, the function will return single-occurrence taxa. 
#' @param bybin (\code{logical}): The type of the filtering process. Was it supposed to be applied to bin-specific subsets (\code{TRUE}), or the whole data (\code{FALSE})? Setting this argument to \code{TRUE} will return a \code{list} class object, where every element of the list is a bin-specific \code{character} vector. This settig also removes all \code{NA} entries form \code{bin} variable.
#' @param na.rm (\code{logical}): If \code{var} is not \code{NULL}, setting this argument to \code{TRUE} removes all rows where var is \code{NA}. Otherwise these will be returned as singletons.
#' @examples
#' # load example dataset
#'   data(corals)
#'
#' # Example 1. single-occurrence taxa
#'   singOcc <- singletons(corals, tax="genus", bin="stg")
#' 
#' # Example 2. output for every bin
#'   singOccBin <- singletons(corals, tax="genus", bin="stg", bybin=TRUE)
#' 
#' # Example 3. single-interval taxa (all)
#'   singInt <- singletons(corals, tax="genus", var="stg")
#' 
#' # Example 4. single interval taxa (for every bin)
#'   singIntBin <- singletons(corals, tax="genus", var="stg", bin="stg", bybin=TRUE)
#' 
#' # Example 5. single reference taxa (total dataset)
#'   singRef <- singletons(corals, tax="genus", var="reference_no")
#' 
#' #  Example 6. single reference taxa (see description for differences )
#'   singRefBin <- singletons(corals, tax="genus", var="reference_no", bin="stg", bybin=TRUE)
#' 
#' @export
singletons <- function(dat, tax="clgen", var=NULL, bin=NULL, bybin=FALSE, na.rm=TRUE){
	
	if(!is.logical(na.rm)) stop("The argument 'na.rm' has to be logical.")
	if(!is.logical(bybin)) stop("The argument 'bybin' has to be logical.")

	if(na.rm){
		if(!is.null(var)){
			dat <-dat[!is.na(dat[, var, drop=TRUE]),]
		}
	}

	if(bybin){
		if(is.null(bin)) stop("You cannot do by-bin filtering, if you do not provide a bin variable.")
		if(!bin%in%colnames(dat)) stop("Invalid bin argument.")

		dat <-dat[!is.na(dat[, bin, drop=TRUE]),]
		binVar <- dat[,bin, drop=TRUE]
	}

	# taxon variable
	if(!tax%in%colnames(dat)) stop("The argument 'tax' has to be a column name of 'dat'. ")
	taxVar <- dat[, tax, drop=TRUE]

	# single-occurrence taxa
	if(is.null(var)){
		if(bybin){
			# contingency table
			occTab <- table(binVar, taxVar)
			class(occTab) <- "matrix"
			
			# number of occurrences of every taxon
			occTaxa<-apply(occTab, 2, sum)
	
			# which are singletons
			singTab <- occTab[,occTaxa==1, drop=FALSE]
	
			# list of single-occurrence taxa in every bin
			taxList <- apply(singTab, 1, function(w){
				colnames(singTab)[as.logical(w)]
			})

		}else{
			tTax <- table(taxVar)
			taxList <- names(tTax)[tTax==1]
		}

	}else{
		if(!var%in%colnames(dat)) stop("The argument 'var' has to be a column name of 'dat'. ")

		# is it bybin - binVar is required
		if(bybin){
			# is the var bin itself?
			if(bin==var){
				# the same as below
				crossTab <- table(taxVar, binVar)
				class(crossTab)  <- "matrix"
		
				# make it presence absence
				crossTab[crossTab>0] <- 1
		
				# in how many intervals does a taxon occurr?
				taxInterval<- apply(crossTab, 1, sum)
		
				# which are single interval? 
				singTab <- crossTab[taxInterval==1,, drop=FALSE]
		
				# single interval taxa in every interval
				taxList<-apply(singTab, 2, function(w){
					rownames(singTab)[as.logical(w)]
		
				})
			# no, it is a different variable
			}else{
				varVar <- dat[,var, drop=TRUE]
				# define core variables
				rows <- 1:nrow(dat)
	
				taxList <- tapply(INDEX=binVar, X=rows, function(w){

					combTab <- unique(cbind(taxVar[w],varVar[w]))
					
					# how many references desribe a taxon - var==NA is also there!
					tTax <- table(combTab[,1, drop=TRUE])
	
					singVarTax <- names(tTax)[tTax==1]
					return(singVarTax)
				})
			
			}

		# no binVar is required
		}else{
			varVar <- dat[,var, drop=TRUE]
				
			crossTab <- table(taxVar, varVar)
			class(crossTab)  <- "matrix"
	
			# make it presence absence
			crossTab[crossTab>0] <- 1
	
			# in how many entries does a taxon occurr?
			taxVarone<- apply(crossTab, 1, sum)
	
			# which are single-unit, or have NA entries in var?
			singTab <- crossTab[taxVarone==1 | taxVarone==0,, drop=FALSE]
	
			# by every item
			taxList<-apply(singTab, 2, function(w){
				rownames(singTab)[as.logical(w)]
	
			})
	
			taxList<- sort(unlist(taxList))
			names(taxList) <- NULL
	
		}

	}
	return(taxList)
}
