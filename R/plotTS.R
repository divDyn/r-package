#' Time series plotting using a custom timescale object
#'
#' This function allows the user to quickly plot a time scale data table
#'
#' As most analysis use an individually compiled time scale object, in order to ensure compatibility between the analyzed and plotted values, the time scale table used for the analysis should be plotted rather than a standardized table.
#' @param tsdat (data frame): time scale data frame
#' @param boxes (character): column name indicating the names that should be plotted as a timescale
#' @param ylim (numeric): the vertical extent of the plot, analogous to the same argument of plot(). By default it is set to the [0,1] interval.
#' @param xlim (numeric): the horizontal extent of the plot, analogous to the same argument of plot(). By default it is set to  plot the entire table. If a numeric vector of length=2 is supplied, it will be interpreted as the standard xlim argument and the plot will be displayed based on numerically constrained ages. If it an integer vector with length>2, the interval corresponding to the row indices of the table will be plotted.
#' @param prop (numeric value): proportion of the vertical extent of the plot for the timescale
#' @param gap (numeric value): proportion of the vertical extent of the plot that should be a gap betwen the timescale and the plot.
#' @param bottom (character): column name of the table for the variable that contains the older ages of intervals.
#' @param mid (character): column name of the table for the variable that contains the mean ages of intervals.
#' @param top (character): column name of the table for the variable that contains the earliest ages of intervals.
#' @param shading (character): column name used for the shading.
#' @param shading.col (character): name of colors that will be used for the shading, if shading is set.
#' @param plot.args (list): list of arguments that will be passed to the main plot() function.  Can be useful for the suppression of axes, font change etc.
#' @param labels.args (list): list of arguments that will be passed to the text() function that draws the labels. 
#' @param boxes.args (list): list of arguments that will be passed to the rect() function that draws the rectangles of time intervals.
#' @examples
#'	data(stages)
#'	  plotTS(stages, boxes="per", shading="series")
#'
#'	# only the Mesozoic, custom axes
#'	  plotTS(stages, boxes="period", shading="stage", xlim=52:81, 
#'	    plot.args=list(axes=F, main="Mesozoic"))
#'	  axis(1, at=seq(250, 75, -25), labels=seq(250, 75, -25))
#'	  axis(2)
#'	
#'	# only the Triassic, use the supplied abbreviations
#'	  plotTS(stages, boxes="short", shading="stage", xlim=c(250,199), 
#'	    ylab="variable", labels.args=list(cex=1.5, col="blue"), 
#'	    boxes.args=list(col="gray95"))
#' @export
plotTS<-function(tsdat,  boxes, ylim=c(0,1), xlim=NULL, prop=0.05, gap=0,
	bottom="bottom", mid="mid",top="top",
	xlab="age (Ma)", ylab="",
	shading=NULL,shading.col=c("white", "gray80"),
	plot.args=NULL,
	boxes.args=NULL,
	labels.args=NULL){
	
#	tsdat<-stages
#	boxes<-"per"
#	bottom="bottom"
#	mid="mid"
#	top="top"
#	ylim=c(0,1)
#	xlab="age (Ma)"
#	ylab=""
#	shading="series"
#	shading.col=c("white", "gray80")
#	
#	boxes.args<-NULL
#	plot.args<-NULL
#	labels.args<-NULL
	
	# length(ylim) should be 2, and numeric, default
	if(!is.numeric(ylim)){
		stop("ylim should be a numeric vector of length 2")
	}
	
	if(length(boxes)>1){
		stop("Only one column can be used to plot the boxes.")
		
	}else{
		if(!boxes%in%colnames(tsdat)){
			stop("The referenced 'boxes' column does not exist.")
		}
	}
	
	if(length(shading)>1){
		stop("Only one column can be used to plot the shades.")
		
	}else{
		if(!shading%in%colnames(tsdat)){
			stop("The referenced 'shading' column does not exist.")
		}
	}
	
	
	# xlim defense
	if(is.null(xlim)){
		xlim <- c(max(tsdat[,bottom],na.rm=T),min(tsdat[,top],na.rm=T))
	} else{
		if(is.numeric(xlim)){
			if(length(xlim)>2){
				if(sum(xlim%%1)==0 & sum(xlim<=0)==0 & sum(xlim>nrow(tsdat))==0){
					xlim <- c(max(tsdat[xlim,bottom],na.rm=T),min(tsdat[xlim,top],na.rm=T))
				}
			}
			if(length(xlim)==1){
				stop("Invalid xlim values.")
			}
		}
	}
	
	# the ylim adjustment
	boxesTop<-ylim[1]-diff(range(ylim))*gap
	plotBottom<-ylim[1]
	
	ylim[1]<-ylim[1]-diff(range(ylim))*prop-diff(range(ylim))*gap
	
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
			x=NULL,
			y=NULL,
			xlim=xlim, 
			ylim=ylim,
			xlab=xlab, 
			ylab=ylab,
			xaxs="i",
			yaxs="i"
		)
		plotArgs<-c(plotArgs, plot.args)
		do.call(plot, plotArgs)
		
	# the box drawing
	boxLev<-unique(tsdat[,boxes])
	xLeft<-rep(NA, length(boxLev))
	xRight<-rep(NA, length(boxLev))
	yTop<-rep(NA, length(boxLev))
	yBottom<-rep(NA, length(boxLev))
	labMid<-rep(NA, length(boxLev))
	for(i in 1:length(boxLev)){
		# the rectangles
			boolSel<-boxLev[i]==tsdat[,boxes]
			curBottom<-max(tsdat[boolSel,bottom], na.rm=T)
			curTop<-min(tsdat[boolSel,top], na.rm=T)
			
			# export the positions
			xLeft[i]<-curBottom
			xRight[i]<-curTop
			yTop[i]<-boxesTop
			yBottom[i]<-ylim[1]
			labMid[i]<-mean(c(curBottom, curTop))
	}
	
	#boxes
	# inner arguments
	boxArgs<-list(
		xleft=xLeft, 
		xright=xRight, 
		ytop= yTop,
		ybottom=yBottom
	)
	
	# process the supplied arguments, in case they are in the table
	boolUnique<-boxes.args%in%colnames(tsdat)
	if(sum(boolUnique)>0){
		ind<-which(boolUnique)
		for(i in ind){
			boxes.args[[i]]<-unique(tsdat[,boxes.args[[i]]])
		}
	}
	
	# combine with outer arguments
	boxArgs<-c(boxArgs, boxes.args)
	do.call(rect, boxArgs)
				
		
		
	# the labels
	# inner arguments
	labArgs<-list(
		label=boxLev, 
		y=(boxArgs$ytop+boxArgs$ybottom)/2, 
		x= labMid
	)
		
	# combine inner and outer arguments
	# process the supplied arguments, in case they are in the table
	boolUnique<-labels.args%in%colnames(tsdat)
	if(sum(boolUnique)>0){
		ind<-which(boolUnique)
		for(i in ind){
			labels.args[[i]]<-unique(tsdat[,labels.args[[i]]])
		}
	}
	
	labArgs<-c(labArgs, labels.args)
	do.call(text,labArgs)
	
	# shading
	if(!is.null(shading)){
		shadeLev<-unique(tsdat[,shading])
		for(i in 1:length(shadeLev)){
			
			# the rectangles
			boolSel<-shadeLev[i]==tsdat[,shading]
			curBottom<-max(tsdat[boolSel,bottom], na.rm=T)
			curTop<-min(tsdat[boolSel,top], na.rm=T)
			if(i%%length(shading.col)){
				rect(xleft=curBottom, xright=curTop, ybottom=plotBottom, ytop=ylim[2], border=NA, col=shading.col[i%%length(shading.col)+1])
			}else{
				rect(xleft=curBottom, xright=curTop, ybottom=plotBottom, ytop=ylim[2], border=NA, col=shading.col[i%%length(shading.col)+1])
			}
		
		}
		rect(xright=xlim[2], xleft=xlim[1], ytop=ylim[2], ybottom=plotBottom)
	}
	# force box()	
	if(gap>0){
		rect(xleft=xlim[1], xright=xlim[2], ytop=plotBottom, ybottom=boxesTop,border="white")
		abline(h=c(boxesTop, plotBottom))
	}
	
}

