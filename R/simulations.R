#' Tools useful for simulation dodels
#'
#' Helpful collection of functions.
#' 
#' The function \code{wienerts} generate a time series of integer values with a random walk (Wiener process). \code{expOri} and \code{expExt} are
#' used for scaling origination and extinction probabilities based on a carrying capacity to reach equilibrial dynamics.
#'
#' @param S Starting value
#' @param steps Number of steps to be generated
#' @param sd Standard deviation of a Gaussian distribution to be taken in every time step, will be rounded to integer values.
#' @param rounding Logical - should integer series be produced?
#' @return A numeric vector (time series) with steps as the the length, starting with S.
#' @rdname simultools
#' @examples
#' K<- wienerts(S=500, sd=10)
#' plot(K)
#' @export
wienerts <- function(S=1000, steps=22400, sd=0.5, rounding=TRUE){
	theDiff <- rnorm(steps-1,0, sd )
	if(rounding) theDiff <- round(theDiff)
	c(S, S+cumsum(theDiff))
}


#' Functions to scale originations
#'
#' @param p The base probabilities.
#' @param k Carrying capacity.
#' @param s Richness in previous slice
#' @rdname simultools
#' @export
expOri<- function(p, k, s){
	exp(1)^(k/s)*p / exp(1)
}

#' Functions to scale extinctions
#'
#' @rdname simultools
expExt<- function(p, k, s){
	exp(1)^(s/k)*p / exp(1)
}


#' Equilibrial Origination-Extinction (Birth-Death) dynamics with Mass Extinctions
#'
#' Extinction and origination cannot happen in the same timestep. The returned data.frame
#' records lineages as rows, and it has four columns:
#' - id: the unique identifier of the lineage
#' - fad: the first appearance date of the lineage (time step number).
#' - lad: the last appearance date of the lineage (time step number). The extant lineages
#' have missing values as their LAD
#' - parent: the parent (ancestor) linage id.
#' The model was originally described in Kocsis, 2015. Analysis of global diversity patterns and dynamics of selected Mesozoic marine invertebrate groups (ELTE, PhD Thesis)
#'
#' @param S The number of starting shites
#' @param K Carrying capacity. A single value or a vector of values.
#' If this is a single value, then it is replicated 'steps' times.
#' @param pExtBase Baseline per timestep, per lineage probability of extinction.
#' @param pOriBase Baseline per timestep, per lineage probability of origination.
#' @param fExt Function to scale extinction probabilities to achieve equilibrial dynamics.
#' The function has to have the arguments: the baseline probability (p), the equilibrial richness (k),
#' and standing richness (s) in the previous time step.
#' @param fExt Function to scale extinction probabilities to achieve equilibrial dynamics.
#' The function has to have the arguments: the baseline probability (p), the equilibrial richness (k),
#' and standing richness (s) in the previous time step.
#' @param steps The number of time steps.
#' @param meP The extinction probability during a mass extinctions.
#' @param me Mass extinction intervals (time step numbers).
#' @param plot Should a simple plot of the carrying capacity and the simulated richness be plotted?
#' @param extNA Should the starting and extant taxa be denoted with the LAD of NA (TRUE), or with the last time step?
#' @return A data.frame with four columns: id, fad, lad, and parent. See description for
#' for the details.
#' @examples
#' set.seed(64)
#' K<- wienerts(S=1000, sd=20)
#' plot(K)
#' one <- equilibrialBD(S=1000, K=K, pExt=0.001, pOri=0.001)
#' @rdname simulations
#' @export
equilibrialBD <- function(S=200, K=200, pExtBase=0.001, fExt=expExt,
	pOriBase=0.001, fOri=expOri, steps=22400, meP = 0.05, me = 9979:9999,
	plot=TRUE, counter=TRUE, extNA=TRUE){

	# some defense
	if(length(steps)>1) stop("Please provie a single positive integer as the time steps.")
	if(steps<1) stop("Please provie a single positive integer as the time steps.")

	#  carrying capcaity
	if(length(K)==1) K <- rep(K, steps)
	if(length(K)!=steps) stop("Wrong carrying capacity length.")

	# final data.frame components
	# record of extinct lineages
	id <- integer()
	fad <- integer()
	lad <- integer()
	parent <- integer()

	# live variables
	fadAlive <- rep(1L, S)
	idAlive <- 1:S
	parentAlive <- rep(NA, S)

	# standing richnesss
	richness <- rep(NA, steps)
	richness[1] <- S

	# total number of lineages
	n <- S

	# markov chain - non parallel
	for(i in 2:steps){
		# current carrying capacity
		k <- K[i]
		# previous richenss
		s <- richness[i-1]

		# 1. Probability scaling
		# origination probs
		pOri <- fOri(p=pOriBase, k=k, s=s)

		# extinctions probs
		# is this ME time
		if(any(i==me)){
			pExt <- meP
		# normal
		}else{
			pExt <- fExt(p=pExtBase, k=k, s=s)
		}

		# 2. Origination
		originations <- runif(s, 0, 1)<=pOri
		nNew <- sum(originations, na.rm=TRUE)


		# 3. Extinctions
		extinctions <- runif(s, 0,1) <=pExt
		indExtinctions <- which(extinctions)
		indExtinctions

		# if there are extinctions
		if(length(indExtinctions)>0){

			# write out the dead ones
			fad <- c(fad, fadAlive[indExtinctions])
			lad <- c(lad, rep(i, length(indExtinctions)))
			id <- c(id, idAlive[indExtinctions])
			parent <- c(parent, parentAlive[indExtinctions])

			# update lives ones
			fadAlive <- fadAlive[-indExtinctions]
			idAlive <- idAlive[-indExtinctions]
			parentAlive <- parentAlive[-indExtinctions]

		}

		# aggregate new lineages
		if(nNew>0){
			# new first appearance
			fadNew <- rep(i, nNew)
			parentNew <- idAlive[originations]

			# fad
			fadAlive <- c(fadAlive, fadNew)
			parentAlive <- c(parentAlive, parentNew)
			idAlive <- c(idAlive, (n+1):(n+nNew))

			# update the number of lineages
			n <- n+ nNew
		}

		# update richness
		richness[i] <- length(idAlive)

		if(counter){
			cat(i, "\r")
			flush.console()
		}

	}

	# add the survivng live ones to the record
	fad <- c(fad, fadAlive)
	if(extNA){
		lad <- c(lad, rep(NA, length(fadAlive)))
	}else{
		lad <- c(lad, rep(i, length(fadAlive)))
	}

	id <- c(id, idAlive)
	parent <- c(parent, parentAlive)

	if(plot){
		plot(NULL, NULL, ylim=c(0, max(c(K, richness))), xlim=c(0, steps),
			xlab="timesteps", ylab="Number of taxa")
		lines(1:steps, richness)
		lines(1:steps, K, col="blue")
		legend("bottomright", legend=c("Simulated richness", "Equilibrium richness"),
			bty="n", col=c("black", "blue"), lwd=1, lty=1)
	}

	result <- data.frame(
		id, fad,lad, parent
	)

	return(result)

}

#' Create pseudo-records from an FAD-LAD matrix
#'
#' Create a fully-sampled occurrence record data frame from an FAD-LAD matrix. 
#' Defaulted to be used with the fadlad function. Can also be used to work with \code{equilibrialBD}.
#' @param x FAD-LAD matrix or data.frame.
#' @param fad The column name of the first appearances.
#' @param lad The column name of the last appearances.
#' @param binned Logical flag indicating whether the FAD and LAD dates are given as integer bins or not (only binned).
#' @return An occurrence record data.frame
#' @export
pseudorecords <- function(x, tax="row.names", fad="FAD", lad="LAD", binned=TRUE){
	# separate taxon column from data.frame
	if(tax!="row.names") {
		taxon <- x[, tax]
	}else{
		taxon <- rownames(x)
		tax <- "taxon"
	}

	# defend fad and lad against no column

	# numeric vectors
	fadVec <- x[, fad]
	ladVec <- x[, lad]

	# test for numeric!

	# binned?
	if(binned){

		# calculate durations (min==1)
		duration <- ladVec - fadVec + 1

		# number of pseudorecords
		rows <- sum(duration)

		# allocate space for pseudorecords
		taxonResult <- rep(NA, rows)
		binResult <- rep(NA, rows)

		# position tracker
		pos <- 1

		for(i in 1:length(taxon)){

			# the number of bins
			offset <- duration[i]

			# position where these need to be saved
			indices <- pos:(pos+offset-1)

			# the bin information
			binResult[indices] <- fadVec[i]:ladVec[i]

			# taxon information
			taxonResult[indices] <- rep(taxon[i], offset)

			# update position tracker
			pos <- pos + offset

		}

		res <- data.frame(taxonResult, binResult)
		colnames(res) <- c(tax, "bin")


	}else{
		stop("Not yet!.")
	}

	return(res)

}


## #' Omitting records from an occurrence data frame.
## #'
## #' This method is uniform across lineages and time.
## #' @param x The occurrence data frame
## #' @param p Overall preseravation probability
## preserveUniform <-function(x, p){
## 	x[runif(1:nrow(x), 0,1)<=p,]
## }
