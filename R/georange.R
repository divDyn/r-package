#' Estimation of geographic ranges from occurrence data
#'
#' Geographic range as a function of a set of coordinates or sample/site/cell membeships.
#' 
#' Multiple estimators of geographic ranges are implemented based on coordinates or cell identifiers. The function outputs a vector of the results based on the calculation methods specified in \code{methods}.
#'
#' @param x \code{(data.frame)} Occurrence table containing the coordinates/locality memberships as variables.
#'
#' @param lng (\code{character}) The variable name of the longitudes, required for the \code{"co"}, \code{"mst"} and \code{"mgcd"} methods. 
#' 
#' @param lat (\code{character}) The variable name of the latitudes, required for the \code{"co"}, \code{"mst"} and \code{"mgcd"} methods. 
#' 
#' @param loc (\code{character}) The variable name of the locality entries: cells, site or samples, required for the \code{"lo"} method.
#' @param method (\code{character}) Geographic range estimator method. Can take multiple entries (concatenating the results in a vector). The following methods are implemented: \code{"co"}: coordinate-based occupancy, the number of different coordinate pairs; \code{"lo"}: locality-based occupancy for sites, samples or geographic cells, number of different entries in a variable. \code{"mst"}, the total length of a minimum spanning tree of the point cloud, based on the great circle distances between points (requires the 'vegan' and 'icosa' packages). \code{"mgcd"}, maximum great-circle distance that can be measured in the point cloud (this version is limited to half the circumference of the equator, requires the 'icosa' package). 
#' @examples
#' data(corals)
#' # select a  taxon from a certain time slice
#'   bitax <- corals[corals$stg==69 & corals$genus=="Microsolena",]
#'   georange(bitax, lng="paleolng", lat="paleolat", method="co")
#' 
#' @export
georange <- function(x, lng=NULL, lat=NULL, loc=NULL, method="co"){
	
	if(length(method)>0){
		# the final output 
		allRanges <- rep(NA, length(method))
		names(allRanges) <- method
	
		#do not use dupliate entries
		unDat <- unique(x[,c(lng, lat, loc), drop=FALSE])
		
		# "lo"
		if(any(method=="lo")){
			if(!is.null(loc)){
				if(length(dim(unDat))>1){
					vec <- unique(unDat[,loc, drop=TRUE])
					allRanges["lo"]<-length(vec[!is.na(vec)])
					if(allRanges["lo"]==0) allRanges["lo"]<-NA
				}else{
					vec <- unique(unDat)
					allRanges["lo"]<-length(vec[!is.na(vec)])
					if(allRanges["lo"]==0) allRanges["lo"]<-NA
				}
			}else{
				warning("No locality variable (loc) is provided the \"lo\" method will return NA")
			}
		}
		
		#coordinates are provided
		notRun <- NULL

		# coorsdinate entries required
	#	if(any(method=="mgcd") | any(method=="mst") | any(method=="sch") | any(method=="co")){
		if(any(method=="mgcd") | any(method=="mst") | any(method=="co")){
			if( !is.null(lng) & !is.null(lat)){
				# the basic coordiantes
				coordat <- unDat[, c(lng, lat)]
				# omit missing entries
				coordat <- coordat[!is.na(coordat[,lng, drop=TRUE]) & !is.na(coordat[,lat, drop=TRUE]), ]
				occNo <- nrow(coordat)

			}
		}

		# co method
		if(any(method=="co")){
			if( !is.null(lng) & !is.null(lat)){
				if(occNo==0){
					allRanges["co"] <- NA
				}else{
					allRanges["co"] <- nrow(unique(coordat))
				}
			}else{
				notRun <- c(notRun, "co")
			}
		}
	
		if(any(method=="mgcd") | any(method=="mst")){
			if( !is.null(lng) & !is.null(lat)){
				# if there are actually occurrences passed to the function
				if(occNo>0){
					if(requireNamespace("icosa", quietly=TRUE)){
						# transform them
						transCoords <- icosa::PolToCar(coordat)
						# calculate distance matrix
						distanceMatrix<- icosa::arcdistmat(transCoords)
					}else{
						distanceMatrix<-NA
						warning("The icosa package is required to use methods 'mgcd' and 'mst'.")
					}
				}
			}
		}
	
		if(any(method=="mgcd")){
			if( !is.null(lng) & !is.null(lat)){
				if(occNo==0){
					allRanges["mgcd"]<-NA
				}
				if(occNo>0){
					allRanges["mgcd"]<-max(distanceMatrix)
				}

			}else{
				notRun <- c(notRun, "mgcd")
			}
		}
		if(any(method=="mst")){
			if( !is.null(lng) & !is.null(lat)){
				if(occNo==0) {
					allRanges["mst"]<-NA
				}
				if(occNo==1){
						allRanges["mst"]<-0
				}
				if(occNo>1){
					if(requireNamespace("vegan", quietly=TRUE)){
						stree<-vegan::spantree(distanceMatrix)
						allRanges["mst"]<-sum(stree$dist)
					}else{
						allRanges["mst"]<-NA
						warning("The vegan package is required to use method 'mst'.")
					}
				}
			}else{
				notRun <- c(notRun, "mst")
			}
		}
		if(!is.null(notRun)){
			warning(paste("No lat/lng variables are provided, the \"", paste(notRun, collapse="\", \""), "\" methods will return NAs.", sep=""))
		}
	#	if(any(method=="sch")){
	#		if( !is.null(lng) & !is.null(lat)){
	#			if(occNo==0) allRanges["sch"] <- NA
	#			if(occNo==1 | occNo==2){
	#				allRanges["sch"]<- 0
	#			}
	#			if(occNo>2){
	#				if(requireNamespace("icosa", quietly=TRUE)){
	#					allRanges["sch"]<- icosa::surfacechullsphere(unDat[, c(lng, lat)])
	#				}else{
	#					allRanges["sch"] <-NA
	#				}
	#			}
	#		}
	#	}
		return(allRanges)
	}else{
		NULL
	}
}




