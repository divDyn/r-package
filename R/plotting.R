#' Time series plotting using a custom time scale 
#'
#' This function allows the user to quickly plot a time scale data table
#'
#' As most analysis use an individually compiled time scale object, in order to ensure compatibility between the analyzed and plotted values, the time scale table used for the analysis could be plotted rather than a standardized table. Two example tables have been included in the package (\code{\link{stages}} and \code{\link{tens}}) that can serve as templates.
#' @param tsdat \code{(data frame)}: The time scale data frame.
#' @param boxes \code{(character)}: Column name indicating the names that should be plotted as boxes of the timescale.
#' @param ylim \code{(numeric)}: The vertical extent of the plot, analogous to the same argument of \code{\link[graphics]{plot}}. By default it is set to the \code{[0,1]} interval.
#' @param xlim \code{(numeric)}: The horizontal extent of the plot, analogous to the same argument of \code{\link[graphics]{plot}}. By default it is set to plot the entire table. If a numeric vector with two values is supplied, it will be interpreted as the standard \code{xlim} argument and the plot will be displayed based on numerically constrained ages. If it is an integer vector with more than two values are plotted, the interval corresponding to the row indices of the table will be plotted.
#' @param prop \code{(numeric)}: Proportion of the vertical extent of the plot to display the the time scale at the bottom.
#' @param xlab \code{(character)}: The label of the time axis.
#' @param ylab \code{(character)}: The label of the data axis.
#' @param gap \code{(numeric)}: Proportion of the vertical extent of the plot that should be a gap betwen the time scale and the plot.
#' @param bottom \code{(character)}: Column name of the table for the variable that contains the older ages of intervals.
#' @param top \code{(character)}: Column name of the table for the variable that contains the earliest ages of intervals.
#' @param boxes.col \code{(character)}: Column name of the colour codes for the boxes. Each entry in this column has to correspond to an entry in the \code{boxes} column. It also overrides the \code{col} entries in the \code{boxes.args} argument.
#' @param shading \code{(character)}: Column name used for the shading. By default, no shading will be drawn (\code{shading = NULL}).
#' @param shading.col \code{(character)}: Colors that will be used for the shading, if shading is set. It is either a single column of the \code{tsdat} object with color codes, or multiple color entries. The provided colors will be repeated as many times as necessary.
#' @param plot.args \code{(list)}: Arguments that will be passed to the main \code{\link[graphics]{plot}} function.  Can be useful for the suppression of axes, font change etc.
#' @param labels.args \code{(list)}: Arguments that will be passed to the \code{\link[graphics]{text}} function that draws the labels. Can be \code{list}s of \code{list}s if multiple series of ``boxes`` are used.
#' @param labels \code{(logical)}: Should the labels within the boxes be drawn? Setting this argumnet to \code{FALSE} will not call the \code{\link[graphics]{text}} function that draws the labels. 
#' @param boxes.args \code{(list)}: Arguments that will be passed to the \code{\link[graphics]{rect}} function that draws the rectangles of time intervals.
#' @param rplab \code{logical}: When the right boundary of the plot does not match with any of the boundaries of the time scale boxes (and \code{labels=TRUE}), should the label of the partially drawn box be plotted?
#' @param lplab \code{logical}: When the left boundary of the plot does not match with any of the boundaries of the time scale boxes (and \code{labels=TRUE}), should the label of the partially drawn box be plotted?
#' @return The function has no return value.
#' @examples
#'	data(stages) 
#'	  tsplot(stages, boxes="sys", shading="series")
#'  # same with colours
#'	  tsplot(stages, boxes="sys", shading="series", boxes.col="systemCol") 
#' 
#'	# only the Mesozoic, custom axes
#'	  tsplot(stages, boxes="system", shading="stage", xlim=52:81, 
#'	    plot.args=list(axes=FALSE, main="Mesozoic"))
#'	  axis(1, at=seq(250, 75, -25), labels=seq(250, 75, -25))
#'	  axis(2)
#'	
#'	# only the Triassic, use the supplied abbreviations
#'	  tsplot(stages, boxes="short", shading="stage", xlim=c(250,199), 
#'	    ylab="variable", labels.args=list(cex=1.5, col="blue"), 
#'	    boxes.args=list(col="gray95"))
#'
#'  # colourful plot with two levels of hierarchy
#'    tsplot(stages, boxes=c("short", "system"), shading="series",
#'      boxes.col=c("col", "systemCol"), xlim=c(52:69))
#' @export
tsplot<-function(tsdat,  ylim=c(0,1), xlim=NULL, prop=0.05, gap=0,
	bottom="bottom", top="top",
	xlab="Age (Ma)", ylab="",
	boxes=NULL, boxes.col=NULL,
	shading=NULL,shading.col=c("white", "gray80"),
	plot.args=NULL,
	boxes.args=NULL,
	labels=TRUE, 
	labels.args=NULL,
	lplab=TRUE, rplab=TRUE){
	
#	tsdat<-stages
#	boxes<-"sys"
#	bottom="bottom"
#	top="top"
#	ylim=c(0,1)
#	xlab="Age (Ma)"
#	ylab=""
#	shading="series"
#	shading.col=c("white", "gray80")
#	xlim <- c(100,50)
#	boxes.col<-NULL
#	prop <- 0.05
#	gap<-0
#	
#	boxes.args<-NULL
#	plot.args<-NULL
#	labels<-TRUE
#	labels.args<-NULL
#	rplab=TRUE
#	lplab=TRUE

	tsCols<-colnames(tsdat)
	if(!bottom%in%tsCols) stop("The 'bottom' column is not found.")
	if(!top%in%tsCols) stop("The 'top' column is not found.")
	
	
	# length(ylim) should be 2, and numeric, default
	if(!is.numeric(ylim) | length(ylim)!=2){
		stop("ylim should be a numeric vector of length 2")
	}
	
	if(!labels & !is.null(labels.args)){
		warning("Labels are set not to plotted, yet there are plotting arguments provided!")
	}

	
	if(sum(boxes%in%colnames(tsdat))!=length(boxes)){
		stop("The referenced 'boxes' column does not exist.")
	}
	
	# no boxes
	if(length(boxes)==0) prop <- 0
	
	if(!is.null(boxes.col)){
		if(length(boxes)!=length(boxes.col)) stop("'boxes.col' should have the same length as 'boxes'")
		if(sum(boxes.col%in%colnames(tsdat))!=length(boxes.col)) stop("One of the columns in 'boxes.col' is not found.")
	}


	# multiple layer of boxes
	if(length(prop)==1) prop <-rep(prop, length(boxes))
	
	if(!is.null(shading)){
		if(length(shading)>1){
			stop("Only one column can be used to plot the shades.")
			
		}else{
			if(!shading%in%colnames(tsdat)){
				stop("The referenced 'shading' column does not exist.")
			}
		}
	}
	
	# if the log axis is planned
	yLog<-FALSE
	if(!is.null(plot.args)){
		if(!is.null(plot.args$log)){
			if(is.na(plot.args$log)) stop("Invalid log argument in plot.args.")
			if(plot.args$log=="y"){
				yLog<-TRUE
			
			}
		}
	
	}
	
	# xlim defense
	if(is.null(xlim)){
		xlim <- c(max(tsdat[,bottom, drop=TRUE],na.rm=T),min(tsdat[,top, drop=TRUE],na.rm=T))
	} else{
		if(is.numeric(xlim)){
			if(length(xlim)>2){
				if(sum(xlim%%1)==0 & sum(xlim<=0)==0 & sum(xlim>nrow(tsdat))==0){
					xlim <- c(max(tsdat[xlim,bottom, drop=TRUE],na.rm=T),min(tsdat[xlim,top, drop=TRUE],na.rm=T))
				}
			}
			if(length(xlim)==1){
				stop("Invalid xlim values.")
			}
		}
	}
	
	# the ylim adjustment
	# conditional for axis reversal - change the offset direction
		if(min(ylim[1])<max(ylim[2])){
			 signChanger <-1
		}else{
			signChanger <- -1
		}
	
	if(!yLog){
		boxesTop<-ylim[1]-(diff(range(ylim))*gap)*signChanger
		plotBottom<-ylim[1]
	
		# in case there are boxes entries
		if(length(boxes)>0){
			vertVector<-boxesTop
			
			# for every box level, subtract the necessary amount
			for(i in 1:length(boxes)){
				vertVector<-c(vertVector, vertVector[i]-(diff(range(ylim))*prop[i]*signChanger))
			}
			
			ylim[1] <-vertVector[length(vertVector)]
		}else{
			ylim[1] <- boxesTop
		}
		
	}else{
		boxesTop<-exp(log(ylim[1])-(diff(range(log(ylim)))*gap*signChanger))
		plotBottom<-ylim[1]
	
		# in case there are boxes entries
		if(length(boxes)>0){
			vertVector<-boxesTop
			
			# for every box level, subtract the necessary amount
			for(i in 1:length(boxes)){
				vertVector<-c(vertVector, exp(log(vertVector[i])-(diff(log(range(ylim)))*prop[i]*signChanger)))
			}
			
			ylim[1] <-vertVector[length(vertVector)]
		}else{
			ylim[1] <- boxesTop
		}
	
	}
	
	# check whether the additional arguments work ok or not
	if(!is.null(boxes.args)){
		if(!is.list(boxes.args)){
			stop("Invalid additional arguments supplied for the boxes.args/rect().")
		}
	}
	
	if(!is.null(plot.args)){
		if(!is.list(plot.args)){
			stop("Invalid additional arguments supplied for plot.args/plot().")
		}
	}
	
	if(!is.null(labels.args)){
		if(!is.list(labels.args)){
			stop("Invalid additional arguments supplied for labels.args/text() .")
		}
	}
	
		

	# the empty plot
		# invoke plot
		plotArgs<-list(
			x=1,
			y=1,
			xlim=xlim, 
			ylim=ylim,
			xlab=xlab, 
			ylab=ylab,
			xaxs="i",
			yaxs="i",
			type="n"
		)
		plotArgs<-c(plotArgs, plot.args)
		do.call(graphics::plot, plotArgs)
	
	if(!is.null(boxes)){
		# for every level of boxes
		for(j in 1:length(boxes)){
		
			# the box drawing
			boxLev<-collapse(tsdat[,boxes[j], drop=TRUE])
			xLeft<-rep(NA, length(boxLev))
			xRight<-rep(NA, length(boxLev))
			yTop<-rep(NA, length(boxLev))
			yBottom<-rep(NA, length(boxLev))
			labMid<-rep(NA, length(boxLev))
			for(i in 1:length(boxLev)){
				# the rectangles
					boolSel<-boxLev[i]==tsdat[,boxes[j]]
					curBottom<-max(tsdat[boolSel,bottom, drop=TRUE], na.rm=T)
					curTop<-min(tsdat[boolSel,top, drop=TRUE], na.rm=T)
					
					# export the positions
					xLeft[i]<-curBottom
					xRight[i]<-curTop
					yTop[i]<-vertVector[j]
					yBottom[i]<-vertVector[j+1]
					labMid[i]<-mean(c(curBottom, curTop))
				
			}
			
			# many cases the plot boundary does not coincide with the bounding boxes, and the labels are not positioned right
			# do the plot boundaries coincide with the boxes?
			#left edge
			# if not, then 
			  
			if(!any(xlim[1]==xLeft)){
				# only if the plot boundary is within the timescale itslef
				if(xlim[1]<max(xLeft)){
					# index of value that needs to be adjusted
					leftAdjust<-max(which(xlim[1]<xLeft))
					# put the value in the middle between the plot boundary and the box boundary
					if(lplab) labMid[leftAdjust] <-  (xlim[1]+xRight[leftAdjust])/2

					# if this is the case, the box left boundary has to be adjusted too
					xLeft[leftAdjust] <- xlim[1]
				}
			}
			# same with right edge
			
			if(!any(xlim[2]==xRight)){
				# only if the plot boundary is within the timescale itslef
				if(xlim[2]>min(xRight)){
					# index of value that needs to be adjusted
					rightAdjust<-min(which(xlim[2]>xRight))
					# put the value in the middle between the plot boundary and the box boundary
					if(rplab)  labMid[rightAdjust] <-  (xlim[2]+xLeft[rightAdjust])/2
				
					# if this is the case, the box right boundary has to be adjusted too
						xRight[rightAdjust] <- xlim[2]
				}
			}


			#boxes
			# inner arguments
			boxArgs<-list(
				xleft=xLeft, 
				xright=xRight, 
				ytop= yTop,
				ybottom=yBottom
			)

			# hierarchical input 
			if(class(boxes.args[[1]])=="list"){
				if(length(boxes.args)!=length(boxes)) stop("Number of lists in 'boxes.args' is not the same as the length of 'boxes'")
				boxes.argsIn <- boxes.args[[j]]
			}else{
				boxes.argsIn <- boxes.args	
			}
					
			
			# process the supplied arguments, in case they are in the table
			boolUnique<-boxes.argsIn%in%colnames(tsdat)
			if(sum(boolUnique)>0){
				ind<-which(boolUnique)
				for(i in ind){
					boxes.argsIn[[i]]<-unique(tsdat[,boxes.argsIn[[i]], drop=TRUE])
				}
			}

			if(!is.null(boxes.col)){
				boolRequire <- !seqduplicated(tsdat[,boxes[j], drop=TRUE])

				# the color values of the box
				boxes.argsIn$col <- tsdat[boolRequire, boxes.col[j], drop=TRUE]

			}
			
			# combine with outer arguments
			boxArgs<-c(boxArgs, boxes.argsIn)
			do.call(graphics::rect, boxArgs)
						
			
			
			# the labels
			# should the labels be plotted
			if(labels){

				# inner arguments
				labArgs<-list(
					label=boxLev, 
					y=(boxArgs$ytop+boxArgs$ybottom)/2, 
					x= labMid
				)

				# hierarchical input 
				if(class(labels.args[[1]])=="list"){
					if(length(labels.args)!=length(boxes)) stop("Number of lists in 'labels.args' is not the same as the length of 'boxes'")
					labels.argsIn <- labels.args[[j]]
				}else{
					labels.argsIn <- labels.args	
				}
					
				# combine inner and outer arguments
				# process the supplied arguments, in case they are in the table
				boolUnique<-labels.argsIn%in%colnames(tsdat)
				if(sum(boolUnique)>0){
					ind<-which(boolUnique)
					for(i in ind){
						labels.argsIn[[i]]<-unique(tsdat[,labels.argsIn[[i]], drop=TRUE])
					}
				}
				
				labArgs<-c(labArgs, labels.argsIn)
				do.call(text,labArgs)
			}
		}
	} 
	
	# shading
	if(!is.null(shading)){
		bDupl<-seqduplicated(tsdat[, shading, drop=TRUE])
		shadeLev<-tsdat[!bDupl,shading, drop=TRUE]

		# predefined color set
		shadeLevCol <- NULL

		#use a collumn to get this
		if(length(shading.col)==1){
			if(shading.col%in%colnames(tsdat)){
				shadeLevCol<-tsdat[!bDupl,shading.col, drop=TRUE]

			}
		}
			
		for(i in 1:length(shadeLev)){
			
			# the rectangles
			boolSel<-shadeLev[i]==tsdat[,shading, drop=TRUE]
			curBottom<-max(tsdat[boolSel,bottom, drop=TRUE], na.rm=T)
			curTop<-min(tsdat[boolSel,top, drop=TRUE], na.rm=T)
			# no predefined colors
			if(is.null(shadeLevCol)){
				if(i%%length(shading.col)){
					graphics::rect(xleft=curBottom, xright=curTop, ybottom=plotBottom, ytop=ylim[2], border=NA, col=shading.col[i%%length(shading.col)+1])
				}else{
					graphics::rect(xleft=curBottom, xright=curTop, ybottom=plotBottom, ytop=ylim[2], border=NA, col=shading.col[i%%length(shading.col)+1])
				}
			# predefined colors
			}else{
				graphics::rect(xleft=curBottom, xright=curTop, ybottom=plotBottom, ytop=ylim[2], border=NA, col=shadeLevCol[i])
			}
		
		}
	
		graphics::rect(xright=xlim[2], xleft=xlim[1], ytop=ylim[2], ybottom=plotBottom)
	}
	# force box()	
	if(gap>0){
		graphics::rect(xleft=xlim[1], xright=xlim[2], ytop=plotBottom, ybottom=boxesTop,border="white")
		graphics::abline(h=c(boxesTop, plotBottom))
	}
	
	clip(xlim[1], xlim[2], boxesTop, ylim[2])
	
}


#' Quantile plot of time series 
#'
#' This intermediate-level function will plot a time series with the quantiles shown with transparency values.
#'
#' @param interpolate \code{(logical)}: In case the symmetric method is chosen, the series of quantile values can be interpolated with a LOESS function. 
#' @param x \code{(numeric)}: The x coordinates.
#' @param y \code{(numeric matrix)}: The series of distributions to be plotted. Every row represents a distribution of values. The number of rows must equal to the length of \code{x}. 
#' @param res \code{(numeric)}: If a single value is entered, than this argument represents the number of quantiles to be shown (coerced to 150, if higher is entered). If it is vector of values, it will be interpreted as the vector of quantiles to be shown. If \code{method="symmetric"}, only an odd number of quantiles are plotted. 
#' @param border \code{(character)}: The color of the quantile lines. A single value, by default, no lines are drawn (\code{border=NA}).
#' @param col \code{(character)}: The color of the quantiles, currently just a single color is allowed.
#' @param method \code{(character)}: The default \code{"symmetric"} method will plot the mid quantile range with highest opacity and the shades will be more translucent at the tails of the distributions. The \code{"decrease"} method will decrease the opacity with higher quantiles, which can make more sense for bottom-bounded distributions (e.g. exponential).
#' @param na.rm \code{(logical)}: If set to \code{FALSE}, than rows that are missing from the dataset will be plotted as gaps in the shading. If set to \code{TRUE}, than these gaps will be skipped. 
#' @examples
#' # some random values accross the Phanerozoic
#'	data(stages)
#'	tsplot(stages, boxes="sys", shading="series", ylim=c(-5,5), ylab=c("normal distributions"))
#'	  randVar <- t(sapply(1:95, FUN=function(x){rnorm(150, 0,1)}))
#'	  shades(stages$mid, randVar, col="blue", res=10,method="symmetric")
#'	  
#' # a bottom-bounded distribution (log normal)
#'	tsplot(stages, boxes="sys", shading="series", ylim=c(0,30), ylab="log-normal distributions")
#'	  randVar <- t(sapply(1:95, FUN=function(x){rlnorm(150, 0,1)}))
#'	  shades(stages$mid, randVar, col="blue", res=c(0,0.33, 0.66, 1),method="decrease")	 
#' @export
#' @return The function has no return value.
shades <- function(x, y, col="black", res=10, border=NA,interpolate=FALSE, method="symmetric",na.rm=FALSE){
	if(nrow(y)!=length(x)) stop("length of x and y don't match")
	
	# omit missing?
	if(na.rm){
		logVector<- apply(y, 1, function(a){
			if(sum(is.na(a))==length(a)) FALSE else TRUE
		})
		
		x<-x[logVector]
		y<-y[logVector,]
	}
	
	if(length(res)==1){
		if(res>150) res<-150
	}else{
		# coerce interpolation
		interpolate <- FALSE
	}
	
	if(length(col)>1 | length(border)>1) stop("Please enter only one color as the col and border arguments.")
	
	
	if(interpolate){
		
		colRes<-100*2
		colQ<-seq(0,1,length.out=colRes)
		
		
		resolved<-apply(y, 1, function(x){
			x<-x[!is.na(x)]
			if(length(x)>3){
				quant<-stats::quantile(x,colQ)
				lmod<-stats::loess(quant ~ colQ, span=0.4)
				tap<-stats::predict(lmod, colQ)
			return(tap)
			}else{
				return(rep(NA, colRes))
			}
		
		})
		
		
		
		for(i in 1:ncol(resolved)){
			if(!is.na(sum(resolved[,i]))){
				ran<-range(y[i,], na.rm=T)
				resolved[resolved[,i]<ran[1],i]<-ran[1]
				resolved[resolved[,i]>ran[2],i]<-ran[2]
			}
		}
			
		
	
		hex<-t(grDevices::col2rgb(col))/255
		hex<-grDevices::rgb(red=hex[1],green=hex[2],blue=hex[3],alpha=1/(colRes/2))

	}else{
		
		
		if(length(res)>1){
			qS<-res
			res<-length(qS)
		}else{
			if(res%%2==1){
				res<-res-1
			}
			qS<-seq(0,1,length.out=res)
			
			
		}
		
		
		
		resolved<-apply(y, 1, function(x){
			x<-x[!is.na(x)]
			if(length(x)>3){
				quant<-stats::quantile(x,qS)
				return(quant)
			}else{
				return(rep(NA, res))
			}
		
		})
		
		hex<-t(grDevices::col2rgb(col))/255
		if(method=="symmetric"){
	
			hex<-grDevices::rgb(red=hex[1],green=hex[2],blue=hex[3],alpha=1/(res/2))
		}
		
		if(method=="decrease"){
			hex<-grDevices::rgb(red=hex[1],green=hex[2],blue=hex[3],alpha=1/(res))
		
		}
	
	}
	
	
	colSums<-apply(resolved, 2, sum)
	streaks<-!is.na(colSums)
	sL<-streaklog(streaks)
	
	for(j in 1:sL$runs){
		bSelect<-rep(FALSE, length(x))
		bSelect[sL$starts[j]:(sL$starts[j]+sL$streaks[j]-1)] <- TRUE
		
		xSub<-x[bSelect]
		resSub<-resolved[,bSelect, drop=FALSE]
		
		# in case a single row is subsetted
		if(ncol(resSub)==1){
			# the x coordinate
			actInd<-which(bSelect)
			xPrev<-x[actInd-1]
			xNext<-x[actInd+1]
			
			xPrev<-xSub-0.2*(xSub-xPrev)
			xNext<-xSub-0.2*(xSub-xNext)
			
			xSub<-c(xPrev, xNext)
			
			# double
			resSub<-cbind(resSub, resSub)
		
		
		}
		
		if(method=="symmetric"){
			for(i in 1:(nrow(resSub)/2)){
				xNew<-c(xSub,rev(xSub))
				yNew<-c(resSub[i,],rev(resSub[nrow(resSub)-i+1,])	)
			#	graphics::polygon(xNew,yNew, col=allCols[i+1], border=F)
				graphics::polygon(xNew,yNew, col=hex, border=border)
			}
		}
		
		if(method=="decrease"){
		for(i in nrow(resSub):1){
				xNew<-c(xSub,rev(xSub))
				yNew<-c(resSub[1,],rev(resSub[nrow(resSub)-i+1,])	)
				graphics::polygon(xNew,yNew, col=hex, border=border)
			}
		
		}
	}
	
}
	


#' Plot time series counts or proportions as polygons
#' 
#' This function plots the changing shares of categories in association with an independent variable. 
#' 
#' This function is useful for displaying the changing proportions of a category as time progresses. Check out the examples for the most frequent implementations.
#' 
#' To be added: missing portions are omitted in this version, but should be represented as gaps in the polygons. 
#' 
#' @param x \code{(numeric)}: The independent variable through which the proportion is tracked. Identical entries are used to assess which values belong together to a set. Their values represent the x coordinate over the plot.
#' 
#' @param b (\code{character} or \code{factor}): A single vector with the category designations. This vector will be segmented using the entries of \code{x}.
#' 
#' @param ord \code{(character)}: The parameter of the variable order. Either \code{"up"} (increasing alphabetical order), \code{"down"} (decreasing alphabetical order) or the vector of categories in the desired order.
#' 
#' @param col \code{(character)}: The color of polygons, has to a be a vector with as many entries as there are categories in \code{b}. By default \code{(col=NULL)} this is grayscale.
#' 
#' @param border \code{(character)}: The a single color of the polygon borders. By default (\code{border=NA}), no borders are drawn. 
#' 
#' @param prop \code{(logical)}: Should the diagram show proportions (\code{TRUE}) or counts (\code{FALSE})?
#' 
#' @param xlim \code{(numeric)}: Two values, analogous to the \code{xlim} argument of \code{\link[graphics]{plot}}, and has to exceed the range of \code{x}. The polygons that represent non-zero values with the lowest and highest values of \code{x} will be extended to these \code{x} coordinates. 
#' 
#' @param ylim \code{(numeric)}: If \code{prop=TRUE}, then the argument controls the position of the proportions in the plotting area (useful to show proportions as a sub plot in a plot). If \code{prop=FALSE}, then the entire plotting area will be shifted by a single \code{ylim} value.
#' 
#' @param labs \code{(logical)}: Should the category names be plotted?
#' 
#' @param na.valid \code{(logical)}: If \code{TRUE}, than the missing values will be treated as an independent category. Entries where \code{x} is \code{NA} will be omitted either way.
#' 
#' @param labs.args \code{(list)}: Arguments for the \code{\link[graphics]{text}} function. If one entry for each argument is provided, then it will be applied to all labels. If the number of elements in an argument equals the number of categories to be plotted, then one to one assignment will be used. For example, for 4 categories in total, if the \code{labs.args} \code{list} contains a \code{col} vector element with 4 values, see examples).
#' 
#' @param plot \code{(logical)}: If set to \code{TRUE}, then the function will plot the output. If set to \code{FALSE}, then a matrix with the relevant values will be returned. This output is similar to the output of \code{\link[base]{table}}, but handles proportions instantly.
#' 
#' @param vertical \code{(logical)}: Horizontal or vertical plotting? If \code{FALSE}, the independent variable will be horizontal, if \code{TRUE}, the count/proportion variable will be horizontal. In the latter case \code{xlim} and \code{ylim} has reversed roles.
#' @return The function has no return value.
#' @examples
#' 
#' # dummy examples 
#'   # independent variable
#'   slc<-c(rep(1, 5), rep(2,7), rep(3,6))
#' 
#'   # the categories as they change
#'   v1<-c("a", "a", "b", "c", "c") # 1
#'   v2<-c("a", "b", "b", "b", "c", "d", "d") # 2
#'   v3<-c("a", "a", "a", "c", "c", "d") #3
#'   va<-c(v1, v2,v3)
#' 
#'   # basic function
#'     plot(NULL, NULL, ylim=c(0,1), xlim=c(0.5, 3.5))
#'     parts(slc, va, prop=TRUE)
#'  
#'   # vertical plot
#'     plot(NULL, NULL, xlim=c(0,1), ylim=c(0.5, 3.5))
#'     parts(slc, va, col=c("red" ,"blue", "green", "orange"), xlim=c(0.5,3.5), 
#'       labs=TRUE, prop=TRUE, vertical=TRUE)
#' 
#'   # intensive argumentation
#'     plot(NULL, NULL, ylim=c(0,10), xlim=c(0.5, 3.5))
#'     parts(slc, va, ord=c("b", "c", "d", "a"), col=c("red" ,"blue", "green", "orange"), 
#' 	  xlim=c(0.5,3.5), labs=TRUE, prop=FALSE, 
#' 	  labs.args=list(cex=1.3, col=c("black", "orange", "red", "blue")))
#' 
#'   # just the values
#'     parts(slc, va, prop=TRUE,plot=FALSE)
#' 	
#' # real example
#'   # the proportion of coral occurrences through time in terms of bathymetry
#'   data(corals)
#'   data(stages)
#' 
#'   # time scale plot
#'   tsplot(stages, shading="series", boxes="sys", xlim=c(250,0), 
#'     ylab="proportion of occurrences", ylim=c(0,1))
#'   
#'   # plot of proportions	
#'   cols <- c("#55555588","#88888888", "#BBBBBB88")
#'   types <- c("uk", "shal", "deep")
#'   
#'   parts(x=stages$mid[corals$stg], b=corals$bath, 
#'    ord=types, col=cols, prop=TRUE,border=NA, labs=FALSE)
#'    
#'   # legend
#'   legend("left", inset=c(0.1,0), legend=c("unknown", "shallow", "deep"), fill=cols, 
#'     bg="white", cex=1.4) 
#' 
#' @export
parts<-function(x, b=NULL, ord="up", prop=FALSE, plot=TRUE,  col=NULL, xlim=NULL, border=NULL, ylim=c(0,1), na.valid=FALSE, labs=TRUE, labs.args=NULL, vertical=FALSE){
	
	if(is.factor(b)) b<-as.character(b)
	
	# starting arguments
	if(is.matrix(x)){
		if(ncol(x==2)){
			b<-x[,2]
			x<-x[,1]
		}
	}
	if(length(x)!=length(b)) stop("x and b have to have the same length.")
	
	#filter NAs
	if(!na.valid){
		bothNA<-is.na(b) | is.na(x)
		b<-b[!bothNA]
		x<-x[!bothNA]
	}else{
		# remove NAs in x
		b<-b[!is.na(x)]
		x<-x[!is.na(x)]
		# an assign another
		b[is.na(b)]<-"N/A"
		
	}
		
	# the unique entries
	bLevs<-unique(b)
	
	# the x values
	# increasing order
	if(!is.null(xlim)){
		if(xlim[1]<xlim[2]){
			x2<-sort(x)
		}
		if(xlim[1]>xlim[2]){
			x2<-sort(x, decreasing=T)
		}
	}else{
		x2<-sort(x)
		xlim<-range(x2)
		xR<-diff(xlim)
		xlim[1]<-xlim[1]-0.1*xR
		xlim[2]<-xlim[2]+0.1*xR
	
	}
	# the appropriate order of the bins
	xLev<-unique(x2)
	
	
	# the order of plotting of categoreies
	if(sum(bLevs%in%ord)==length(bLevs)){
		bLevs<-ord
	}else{
		if(length(ord)==1){
			if(ord=="up"){
				bLevs<-sort(bLevs)
			}
			if(ord=="down"){
				bLevs<-sort(bLevs, decreasing=T)
			}
		}else{
			stop("Invalid ord argument.")
		}
	}
	
	# arguments of the labels
	if(labs){
		if(!is.null(labs.args)){
			labs.args<-lapply(labs.args, function(x){
				if(length(x)!=length(bLevs) & length(x)!=1) stop("Provide each category lab an individual entry in arguments or one.")
				if(length(x)==1){
					return(rep(x, length(bLevs)))
				}
				if(length(x)==length(bLevs)){
					return(x)
				}
			
			})
		}
	}
	
	# default colors
	if(is.null(col)){
		dif<-floor(1/(length(bLevs)-1)*100)
		grayVal<-seq(0, 100, by=dif)
		col<-paste("gray", grayVal, sep="")
		names(col)<-bLevs
	}else{
		if(!length(col)==length(bLevs)) stop("The number of colour in 'col' doesn't match the number of categories.")
		names(col)<-bLevs
	
	}
	if(!is.null(border)){
		if(length(border)==1){
			border<-rep(border, length(bLevs))
			names(border)<-bLevs
		}else{
			if(!length(border)==length(ord)) stop("The number of colours in 'border' doesn't match the number of categories.")
		}
	}
	
	allProps<-tapply(INDEX=x, X=b, function(y){
		oneres<-sapply(bLevs, function(z){
			sum(z==y)
		})
		
		#if plotted, return the cumulative
		if(plot){
			if(prop){
				cumsum(oneres/length(y))
			}else{
				cumsum(oneres)
			}
		# if not plotted, then return the values
		}else{
			if(prop){
				oneres/length(y)
			}else{
				oneres
			}
		}
	})
	
	# reorder the by x
	allProps<-allProps[as.character(xLev)]
	
	# the labels
	if(labs){
		theProps<-tapply(INDEX=x, X=b, function(y){
			oneres<-sapply(bLevs, function(z){
				sum(z==y)
			})
			
			if(prop){
				oneres/length(y)
			}else{
				oneres			
			}
		})
		
		theProps<-theProps[as.character(xLev)]
		
		labPos <-sapply(bLevs, function(y){
			lap<-lapply(theProps, function(z){
				z[y]
			})
			res<-unlist(lap)
			xLev[which(res==max(res))[1]]
		})
	}
	
	if(prop){
		yRange<-range(ylim)
		if(length(diff(yRange))>1) stop('ylim argument must have two entries')
		if(diff(yRange)<0) stop("ylim values have to be ascending")
	}
	
	# get the values out
	allCum <-sapply(bLevs, function(y){
		lap<-lapply(allProps, function(z){
			z[y]
		})
		res<-unlist(lap)
		names(res)<-names(lap)
		res
	})
	
	# return values if plot is FALSE
	if(!plot){
		return(allCum)
	}
	
	# rescale if necessary	
	if(prop){
		allCum<-allCum*diff(yRange)
	}
	
	# and shift
	allCum <- allCum+ylim[1]
		
	for(i in length(bLevs):1){
		if(i>1){
			plotY<-c(
				allCum[1,bLevs[i]], 
				allCum[,bLevs[i]], 
				allCum[nrow(allCum),bLevs[i]],
				allCum[nrow(allCum),bLevs[i-1]],
				rev(allCum[,bLevs[i-1]]), 
				allCum[1,bLevs[i-1]]
			)
		
		}else{
			plotY<-c(
				allCum[1,bLevs[i]], 
				allCum[,bLevs[i]], 
				allCum[nrow(allCum),bLevs[i]],
				rep(ylim[1], length(xLev)+2))
		}
		# the xlims
		plotX<-c(xlim[1], xLev, xlim[2], xlim[2], rev(xLev), xlim[1])
	
		if(vertical){
			graphics::polygon(plotY, plotX, col=col[bLevs[i]], border=border[bLevs[i]])
		}else{
			graphics::polygon(plotX, plotY, col=col[bLevs[i]], border=border[bLevs[i]])
		}
	}
	
	if(labs){
		
		for(i in 1:length(labPos)){
			xCoord<-labPos[i]
			
			yVect<-allCum[as.character(labPos[i]), ]
			
			yPos<-which(names(yVect)==names(labPos)[i])
			if(yPos==1){
				yCoord<-mean(c(yVect[yPos],ylim[1]))
			}else{
				yCoord<-mean(c(yVect[yPos],yVect[yPos-1]))
			}
			if(vertical){
				newArgs<-list(
					label=names(labPos)[i],
					y=xCoord,
					x=yCoord
				)	
			}else{
				newArgs<-list(
					label=names(labPos)[i],
					x=xCoord,
					y=yCoord
				)	
			}
			
			plusArgs<-lapply(labs.args, function(u){
				u[i]
			})
			
			# plot the labels
			do.call(graphics::text, c(newArgs, plusArgs))
		}
	
	}
}

	
	


#' Plotting ranges and occurrence distributions through time
#' 
#' Visualization of occurrence data
#' 
#' This function will draw a visual representation of the occurrence dataset. The interpolated ranges will be drawn, as well as the occurrence points.
#' 
#' @param dat \code{(data.frame)}: The occurrence dataset or the FAD-LAD dataset that is to be plotted. The FAD dataset must have numeric variables named \code{"FAD"} and \code{"LAD"}. Taxon ranges will be searched for in the \code{row.names} attribute of the table. 
#' 
#' @param bin (\code{character}): The column(s) containing the entries of the time dimension. Use one column name if you have one estimate for the occurrences and use two if you have a minimum and a maximum estimate. Reveresed axis (ages) are supported too.
#' 
#' @param tax (\code{character}): The column containing the taxon entries.
#' 
#' @param group (\code{character}): By default, all ranges in the plot are treated as parts of the same group. However, one subsetting variable can be named, by which the ranges will be grouped. This has to be a column name in the dataset (see examples).
#' 
#' @param xlim (\code{numeric}) :This argument is used for the subsetting of the taxa. Only those taxa are shown that have ranges within the interval (but ranges are displayed outside of it, if you do not want to plot anything within an interval, use the \code{clip} function)
#' 
#' @param ylim (\code{numeric}): Ranges will be distributed equally between the assigned ylim values. If set to NULL, than it will be based on the plotting area of the open device.
#' 
#' @param occs (\code{logical}): Should the occurrence data be plotted? If you entered two bin column names, than occurrences will be plotted with the ranges of the estimates (segments).
#' 
#' @param gap (\code{numeric}): Evaluated only when the \code{group} argument points to a valid column. The amount of space between the group-specific range charts, expressed as the proportion of the entire plotting area. 
#' 
#' @param total (\code{character}): The name of the range group to be plotted. When multiple groups are used (see \code{group} argument), this is set by the \code{character} values in the column.
#' 
#' @param filt (\code{character}): When xlim filters the taxa, how should they be filtered. \code{"include"} (default) will show all ranges that have parts within the \code{xlim} interval. \code{"orig"} will show only those taxa that originate within the interval.
#'
#' @param labs (\code{character}): Should the taxon labels be plotted?
#'
#' 
#' @param total.args (\code{list}): Arguments that will be passed to the \code{\link[graphics]{text}} function that draws the \code{total} label. If valid grouping is present (see argument \code{group}), then vector entries will be distributed across the groups (see examples.)
#' 
#' @param ranges.args (\code{list}): Arguments that will be passed to the \code{\link[graphics]{segments}} function that draws ranges. If valid grouping is present (see argument \code{group}), then vector entries will be distributed across the groups (see examples.)
#' 
#' @param occs.args (\code{list}): Arguments that will be passed to the \code{\link[graphics]{points}} or \code{\link[graphics]{segments}} functions that draw the occurence points/lines. If you provided two \code{bin} columns, occurrence lines will be drawn instead of points. If valid grouping is present (see argument \code{group}), then vector entries will be distributed across the groups (see examples.)
#' 
#' @param labels.args (\code{list}): Arguments that will be passed to the \code{\link[graphics]{text}} function that draws the labels of taxa. If valid grouping is present (see argument \code{group}), then vector entries will be distributed across the groups (see examples.)
#'
#' @param decreasing (\code{logical}): This parameter sets whether the series of ranges should start from the top \code{decreasing=TRUE} or bottom of the plot \code{decreasing=FALSE}. 
#' 
#' @return The function has no return value.
#' @examples
#'  # import
#'  data(stages)
#'  data(corals)
#'  
#'  # all ranges - using the age uncertainties of the occurrences
#'  tsplot(stages, boxes="sys", xlim=c(250,0))
#'  ranges(corals, bin=c("max_ma", "min_ma"), tax="genus", occs=FALSE)
#'
#'  # or use single estimates: assign age esimates to the occurrences
#'  corals$est<-stages$mid[corals$stg]
#'  
#'  # all ranges (including the recent!!)
#'  tsplot(stages, boxes="sys", xlim=c(250,0))
#'  ranges(corals, bin="est", tax="genus", occs=FALSE)
#'  
#'  # closing on the Cretaceous, with occurrences
#'  tsplot(stages, boxes="series", xlim=c(145,65), shading="short")
#'  ranges(corals, bin="est", tax="genus", occs=TRUE, ranges.args=list(lwd=0.1))
#'  
#'  # z and az separately
#'  tsplot(stages, boxes="series", xlim=c(145,65), shading="short")
#'  ranges(corals, bin="est", tax="genus", occs=FALSE, group="ecology", 
#'    ranges.args=list(lwd=0.1))
#'  	
#'  # same, show only taxa that originate within the interval
#'  tsplot(stages, boxes="series", xlim=c(105,60), shading="short")
#'  ranges(corals, bin="est", tax="genus", occs=TRUE, group="ecology", filt="orig" ,
#'    labs=TRUE, labels.args=list(cex=0.5))
#'  
#' # same using the age uncertainties of the occurrence age estimates
#' tsplot(stages, boxes="series", xlim=c(105,60), shading="short")
#' ranges(corals, bin=c("max_ma", "min_ma"), tax="genus", occs=TRUE, group="ecology", filt="orig" , 
#'    labs=TRUE, labels.args=list(cex=0.5))
#'    
#' # fully customized/ annotated
#' tsplot(stages, boxes="series", xlim=c(105,60), shading="short")
#' ranges(
#'   corals, # dataset
#'   bin="est", # bin column
#'   tax="genus", # taxon column
#'   occs=TRUE, # occurrence points will be plotted
#'   group="growth", # separate ranges based on growth types
#'   filt="orig" , # show only taxa that originate in the interval
#'   ranges.args=list(
#'     lwd=1, # set range width to 1
#' 	   col=c("darkgreen", "darkred") # set color of the ranges (by groups)
#'   ), 
#'   total.args=list(
#'     cex=2, # set the size of the group identifier lablels
#'     col=c("darkgreen", "darkred") # set the color of the group identifier labels
#'   ),
#'   occs.args=list(
#'	   col=c("darkgreen", "darkred"),
#'	   pch=3
#'	 ),
#'   labs=TRUE, # taxon labels will be plotted
#'   labels.args=list(
#'     cex=0.4, # the sizes of the taxon labels
#' 	col=c("darkgreen", "darkred") # set the color of the taxon labels by group
#'   )
#' ) 
#'    
#' 
#' @export
ranges <- function(dat, bin=NULL, tax=NULL, xlim=NULL, ylim=c(0,1), total="", filt="include", occs=FALSE, labs=FALSE, decreasing=TRUE, group=NULL, gap=0,  labels.args=NULL, ranges.args=NULL, occs.args=NULL, total.args=NULL){
	
#	dat<-corals
#	bin<-"est"
#	tax<-"genus"
#	xlim<-NULL
#	group<-NULL
#	ylim<-c(0,1)
#	occs<-TRUE
#	labs<-FALSE
#	decreasing<-TRUE
#	gap <- 0
#	labs.args<- NULL
#	ranges.args<-NULL
#	occs.args <- NULL
#	filt<-"orig"
	
	# is input a FAD-LAD matrix, or an occurrence database
	if(!is.null(dat$FAD) & !is.null(dat$LAD)){
		# copy over
		fl<- dat
		
		# suppress occurrence plotting
		if(occs) message("Range-data were provided, occurrence plotting is not possible.")
		occs <- FALSE
	}else{
		
		if(is.null(bin)) stop("Argument bin is missing.")
		if(is.null(tax)) stop("Argument tax is missing.")
		if(sum(bin%in%colnames(dat))!=length(bin)) stop("Column bin not found.")
		if(!tax%in%colnames(dat)) stop("Column tax not found.")
		
		# calculate the FAD LAD 
		fl<-fadlad(dat, tax=tax, bin=bin, na.rm=T)
	}
	
	
	# get the 
	if(is.null(xlim)){
		rng <- par("usr")
		xlim<-rng[1:2]
	}
	if(is.null(ylim)){
		rng <- par("usr")
		ylim<-rng[3:4]
	}
	
	# determine orientation (1. ages)
	ladGreater <- sum(fl$FAD<fl$LAD)>0
	
	# which taxa to keep
	bKeep<-rep(TRUE, nrow(fl))
	
	if(filt=="include"){	
		if(!ladGreater){
			bKeep <- fl$FAD>= xlim[2] & fl$FAD <= xlim[1] |
			fl$LAD>=xlim[2] & fl$LAD <=xlim[1] |
			fl$FAD>xlim[1] & fl$LAD <=xlim[2]
		}
		if(ladGreater){
			bKeep <- fl$LAD>= xlim[2] & fl$LAD <= xlim[1] |
			fl$FAD>=xlim[2] & fl$FAD <=xlim[1] |
			fl$LAD>xlim[1] & fl$FAD <=xlim[2]
		
		}
	}
	if(filt=="orig"){
		if(!ladGreater){
			bKeep <- fl$FAD >=xlim[2] & fl$FAD<=xlim[1]
		}else{
			bKeep <- fl$LAD >=xlim[2] & fl$LAD<=xlim[1]
		}
	
	}
		
	# filter based on the limits of the plot
	newFL<- fl[bKeep,]
	
	# recursive case
	if(!is.null(group)){
		
		# create group sepcific subsets
		groupLevs<-levels(factor(dat[, group, drop=TRUE]))
		
		# if no color is provided, use default
		if(is.null(ranges.args$col) & is.null(labels.args$col) & is.null(occs.args$col) & length(groupLevs)<11){
			ranges.args$col <- mainHex[1:length(groupLevs)]
			labels.args$col <- mainHex[1:length(groupLevs)]
			occs.args$col <- mainHex[1:length(groupLevs)]
			total.args$col <- mainHex[1:length(groupLevs)]
		} 
		
		# rearrange
		ranges.argsDist<- distribute(ranges.args, length(groupLevs))
		labels.argsDist<- distribute(labels.args, length(groupLevs))
		occs.argsDist<- distribute(occs.args, length(groupLevs))
		total.argsDist<- distribute(total.args, length(groupLevs))
		
		# add the group variable to FAD-LAD 
		#only for occs!
		groupVar<-dat[, group, drop=TRUE]
		names(groupVar)<-dat[, tax, drop=TRUE] 
		newFL$group <- groupVar[row.names(newFL)]
		
		# where to plot
		yRange <- ylim[2]-ylim[1]
		gapBetween <- gap*yRange
		yAvailable <- yRange-(gapBetween*(length(groupLevs)-1))
		
		# the number of groups in the different taxa
		taxNoGroup<-table(newFL$group)
		
		# how much of the plot is available for the group
		yPartGroup<-yAvailable/sum(taxNoGroup)*taxNoGroup
		
		# calculate where the lots should be
		ylimMat <- matrix(c(ylim[1], yPartGroup[1]), ncol=2, nrow=1)
		
		subDat<-dat[dat[, group, drop=TRUE]==groupLevs[1] & dat[, tax, drop=TRUE]%in%row.names(newFL), ]
		
		# first draw everything for the first group		
		ranges(subDat, bin=bin, tax=tax, xlim=xlim, ylim=ylimMat[1,], occs=occs, labs=labs, total=groupLevs[1], decreasing=decreasing, group=NULL, 
			labels.args=labels.argsDist[[1]], ranges.args=ranges.argsDist[[1]], occs.args=occs.argsDist[[1]], total.args=total.argsDist[[1]])
		
		for(i in 2:length(groupLevs)){
			# where should the plot be
			first<-ylimMat[i-1,2]+gapBetween
			second<-first+yPartGroup[i]
			ylimMat<-rbind(ylimMat, c(first, second))
			
			subDat<-dat[dat[, group, drop=TRUE]==groupLevs[i] & dat[, tax, drop=TRUE]%in%row.names(newFL), ]
		
			ranges(subDat, bin=bin, tax=tax, xlim=xlim, ylim=ylimMat[i,],  occs=occs, labs=labs,  total=groupLevs[i],decreasing=decreasing, group=NULL,
				labels.args=labels.argsDist[[i]], ranges.args=ranges.argsDist[[i]], occs.args=occs.argsDist[[i]],total.args=total.argsDist[[i]])
		
		
		}
		
	# base case (one group)
	}else{
		# sort taxa by fad/lad
		if(ladGreater) newOrd <- order(newFL$LAD, newFL$FAD)
		if(!ladGreater) newOrd <- order(newFL$FAD, newFL$LAD)
	
	
		# the reordered ranges
		plotFL <- newFL[newOrd, ]
		
		# vertical positions of taxa
		taxWhere<-seq(ylim[1], ylim[2], length.out=nrow(plotFL)+2)
		
		taxWhere<-taxWhere[c(-1, -length(taxWhere))]
		
		# reverse order for different look
		if(!decreasing){
			taxWhere <-rev(taxWhere)
		}
		
		names(taxWhere)<-rownames(plotFL)
		
		# Drawing the ranges
			# use the supported arguments
			ranArgs<-ranges.args
			
			# overwrite with internals
			ranArgs$x0 <-plotFL[,1]
			ranArgs$x1 <-plotFL[,2]
			ranArgs$y0 <-taxWhere
			ranArgs$y1 <-taxWhere
			
			# function call
			do.call(segments, ranArgs)
			
		# draw the labels
		if(labs){
			# use the supported arguments
			textArgs<-labels.args
			
			# overwrite
			textArgs$pos<-2
			textArgs$label<-names(taxWhere)
			textArgs$y<-taxWhere
			
			if(ladGreater){
				textArgs$x<-plotFL$LAD
			}else{
				textArgs$x<-plotFL$FAD
			}
			# function call
			do.call(text, textArgs)
			
		}
		
		# drawing the total name
			totalArgs<-total.args
			# add the label itself
			totalArgs$label<-total
			
			# default settings
			if(is.null(totalArgs$x)) totalArgs$x<-0.2
			if(is.null(totalArgs$y)) totalArgs$y<-0.2
			
			totalArgs$x<-max(xlim)+diff(xlim)*totalArgs$x
			if(decreasing) totalArgs$y<-min(ylim)+diff(ylim)*totalArgs$y
			if(!decreasing) totalArgs$y<-max(ylim)-diff(ylim)*totalArgs$y
			
			do.call(text, totalArgs)
		
			
		# the occurrences
		if(occs){
			# do the subsetting of the occurence dataset (for plotting)
			plotDat <- dat[dat[,tax]%in%names(taxWhere), ]
			plotCoords<-unique(cbind(plotDat[,bin, drop=TRUE], taxWhere[plotDat[,tax, drop=TRUE]]))
			
			# use the supplied arguments
				occArgs<-occs.args
				
			# draw occurrence points
			if(length(bin)==1){			
		
				
				# overwrite with internals
				occArgs$x <-plotCoords[,1]
				occArgs$y <-plotCoords[,2]
				
				# defaults
				if(is.null(occArgs$pch)) occArgs$pch<-16
				
				# function call
				do.call(points, occArgs)
			}
			
			# draw occurrence ranges
			if(length(bin)==2){
			
				occArgs$x0 <-plotCoords[,1]
				occArgs$x1 <-plotCoords[,2]
				occArgs$y0 <-plotCoords[,3]
				occArgs$y1 <-plotCoords[,3]
				
				# if line width for occurrence not given
				if(is.null(occArgs$lwd)){
					# and if the width of the range is given
					if(!is.null(ranArgs$lwd)){
						occArgs$lwd<-ranArgs$lwd*4
					} else{
						occArgs$lwd<-4
					}
				}
				
				# function call
				do.call(segments, occArgs)
			}
			
		}
		
	}
}



distribute<-function(oneList, target){
	
	smallList<-lapply(oneList, FUN=function(x){
		if(length(x)==1) fin <- rep(x, target)
		if(length(x)==target) fin <- x
		if(length(x)!= 1 & length(x)!=target){
			fin <- x[1:target]
		}
		return(fin)
	
	})
	
	datF<-as.data.frame(smallList, stringsAsFactors=F)
	
	if(sum(is.na(datF))) warning("The number of supplied arguments doesn't match the number of groups")
	
	newList<-list()
	for(i in 1:target){
		one<-datF[i,, drop=FALSE]
		one<-as.list(one)
		newList<-c(newList, list(one))	
	}
	
	if(length(newList)==0) newList<-NULL
	
	
	return(newList)
}

#' Function to plot a series a values with bars that have variable widths
#' 
#' Function to use bars for time series.
#' 
#' People often present time series with connected points, although the visual depiction implies a certain process that describes how the values change between the points.
#' Instead of using simple scatter plots, Barplots can be used to describe series where a single value is the most descriptive of a discreet time bin. The \code{tsbars()} function
#' draws rectangles of different widths with the \code{\link[graphics]{rect}} function, to plot series in such a way.
#' 
#' @param x \code{(numeric)} Vector specifying where the centers of the bars should be on the x axis. 
#' @param y \code{(numeric)} Vector containing the heights of the bars.
#' @param width \code{(numeric)} Vector containing the widths of the bars. Recycling is not supported, has to be either a single numeric value, or a numeric vector with the same length as \code{x} and \code{y}. Automatic width calculation is possible, the default \code{"max"} option sets the bar width even and equal to the the maximum width that can be used evenly witout causing overlaps. The option \code{"half"}, places the boundaries of the bars halfway between the points. This will make the bars' width asymmetrical around the \code{x} coordinates.
#' @param yref \code{(numeric)} Single numeric value in the y dimension indicating common base for the bars.
#' @param gap \code{(numeric)} The amount of gap there should be between the bars (in the unit of the plotting). Defaults to no gaps. 
#' @param vertical \code{(logical)} Switching this option to \code{FALSE} will reverse the x and y dimensions of the plot.
#' @param ... Arguments passed to \code{\link[graphics]{rect}}.
#' @examples
#' # an occurrence-based example
#' # needed data
#'   data(stages)
#'   data(corals)
#' # calculate diversites
#'   dd <-divDyn(corals, tax="genus", bin="stg")
#' # plot range-through diversities
#'   tsplot(stages, xlim=51:94, ylim=c(0,250), boxes="sys")
#'   tsbars(x=stages$mid, y=dd$divRT, width=stages$dur, gap=1, col=stages$col)
#' 
#' @export
#' @return The function has no return value.
tsbars <- function(x, y, width="max", yref=0, gap=0, vertical=TRUE, ...){
#	x<- stages$mid
#	y<- dd$divRT
#	width <- stages$dur
#	yref<-0

	# check if if x has the right length
	if(length(x)!=length(y)) stop("'x' has to have the same length as 'y'")
	if(!is.numeric(x) | !is.numeric(y))  stop("'x' and 'y' have to be numeric")

	if(length(width)!=1){
		changeLeft<- width/2
		changeRight<- width/2
	}else{
		if(width=="max"){
			changeLeft <- min(abs(diff(x))/2, na.rm=TRUE)
			changeRight <- changeLeft
		}
		if(width=="half"){
			ser <- abs(diff(x))/2
			changeLeft<-c(ser[1], ser)
			changeRight<-c(ser, ser[length(ser)])
		}
		if(is.numeric(width)){
			changeLeft<- width/2
			changeRight<- width/2
		}	
	}
	
	# gap between bars
	if(gap!=0){
		changeLeft<-changeLeft-gap/2
		changeRight<-changeRight-gap/2
	}

	# core
	if(length(yref)==1){
		yBot <- rep(yref, length(x))
	}else{
		if(length(yref)==length(x)){
			yBot <- yref
		}
	}

	# left and right borders
	xLeft <- x+changeLeft
	xRight <- x-changeRight

	# function call
	if(vertical){
		rect(ytop=y, ybottom=yBot, xleft=xLeft, xright=xRight, ...)
	}else{
		rect(ytop=xLeft, ybottom=xRight, xleft=yBot, xright=y, ...)
	}
	
}


# default high contrast color palette
mainHex<-c(
	
	# red
	"#F21414",
	
	# green
	"#07F80C",
	
	# blue
	"#0712F8",
	
	# yellow
	"#F8F500",
	
	# purple
	"#9E00F8",
	
	# orange
	"#FF9C00",
	
	# cyan
	"#00F0FF",
	
	# pink
	"#F00FE3",
	
	# brown
	"#744932",
	
	# grass
	"#A5EA21")
