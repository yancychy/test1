kmplot <- function(km, mark=3, simple=FALSE,
    xaxis.at=pretty(km$time), xaxis.lab=xaxis.at, 
    lty.surv=1, lwd.surv=1, col.surv=1,
    lty.ci=0, lwd.ci=1, col.ci=col.surv, #By default (lty.ci=0), confidence intervals are not plotted.
    group.names=NULL, group.order=seq(length(km$n)), extra.left.margin=4, 
    label.n.at.risk=FALSE, draw.lines=TRUE, cex.axis=1,
    xlab='', ylab='', main='', xlim=c(0,max(km$time)), ylim=c(0,1),
    grid=TRUE, lty.grid=1, lwd.grid=1, col.grid=grey(.9),
    legend=!is.null(km$strata), loc.legend='topright',
    ... # ... is passed to par()
    ){
	# Version 2.5.3: 2012/2/29

	# km is the output from survfit() function in survival package.
	
	# xaxis.at specifies where 'n at risk' will be computed and printed.
	# xaxis.lab specifies what will be printed at xaixs.at.  (see example)

	# If group names are long, add extra left margin by setting extra.left.margin to something greater than 0.

	# line specifications (lty.surv, lwd.surv, col.surv) will be recycled.
	# Set lty.ci to 1 if confidence intervals are needed.
	# group.names will overwrite whatever is specified in survfit() output.
	# group.order specifies the order of groups from top in 'n at risk'.  1 is at top, 2 next, and so on.

	op <- par(no.readonly = TRUE) 

	ng0 <- length( km$strata ) ; ng <- max(ng0,1) 
		# When only one group...
		if(ng0==0){ km$strata <- length(km$time) ; names(km$strata) <- 'All' ; legend <- draw.lines <- FALSE } 

	lty.surv <- rep(lty.surv, ng) ; lwd.surv <- rep(lwd.surv, ng) ; col.surv <- rep(col.surv, ng)
	lty.ci <- rep(lty.ci, ng) ;     lwd.ci <- rep(lwd.ci, ng)     ; col.ci <- rep(col.ci, ng)

	## group names and error checking	
	gr <- c(km$strata)
	if( is.null(group.names) ){ group.names <- names(km$strata) }
	if( length(unique(group.names)) != ng ){ stop('\n','length(unique(group.names)) != number of groups.') }
	if( suppressWarnings(any( sort(group.order) != 1:ng)) )
		{ stop('\n', 'Something wrong with group.order.','\n','sort(group.order) must equal 1:', ng, '.') }
	group.names <- gsub(' *$', '', group.names)  #to remove unwanted white spaces in group.names.
	if(ng==1 & (group.names[1]=='group.names') ){ group.names <- 'N at risk' ; label.n.at.risk = FALSE }

	## graphic parameters
	par(list(oma=c(1,1,1,1), mar=c(4+ng,4+extra.left.margin,4,2)+.1))
		if(simple) par( mar=c(3,4,2,1)+.1 )
	par( list(...) )

	## reformat survival estimates
	dat <- data.frame(time=km$time, n.risk=km$n.risk, n.event=km$n.event, survival=km$surv, std.err=km$std.err, 
						lower=km$lower, upper=km$upper, group=rep( group.names, gr) )
	dat.list <- split(dat, f=dat$group)

	## plot (but not survival curves) 
	plot(0,type='n', xlim=xlim, ylim=c(0,1), xaxt='n', yaxt='n', xlab='', ylab='' )
		if(grid){
			par('xpd'=FALSE)
			abline(v=xaxis.at,       lty=lty.grid, lwd=lwd.grid, col=col.grid )
			abline(h=pretty(c(0,1)), lty=lty.grid, lwd=lwd.grid, col=col.grid )
		}
	axis( side=2, at=pretty(c(0,1)), cex.axis=cex.axis )	
	axis( side=1, at=xaxis.at, label=xaxis.lab, line=-0.5, tick=FALSE, cex.axis=cex.axis )
	axis( side=1, at=xaxis.at, label=rep('',length(xaxis.at)), line=0, tick=TRUE )
	title(xlab=xlab, line=1.5, adj=.5, ...) ; title(ylab=ylab, ... )

	if(!simple){
	## write group names
	group.name.pos <- (par()$usr[2]-par()$usr[1]) / -8 ; padding <- abs( group.name.pos / 8 )
		line.pos <- (1:ng)[order(group.order)] + 2
	mtext( group.names, side=1, line=line.pos, at=group.name.pos, adj=1, col=1, las=1, cex=cex.axis )

	## draw matching lines for n at risk.		
	if(draw.lines){
		par('xpd'=TRUE)
		for(i in 1:ng){
		axis(side=1, at=c(group.name.pos+padding,0-2*padding), labels=FALSE, line=line.pos[i]+0.6, lwd.ticks=0,
					 col=col.surv[i], lty=lty.surv[i], lwd=lwd.surv[i] ) }
				}

	## numbers at risk
	kms <- summary(km, times=xaxis.at) ; if(is.null(kms$strata)) kms$strata <- rep(1,length(kms$time) )
	d1 <- data.frame(time = kms$time, n.risk = kms$n.risk, strata = c(kms$strata))
	d2 <- split(d1, f=d1$strata)

	## Right-justifying the numbers 
	ndigits <- lapply(d2, function(x) nchar(x[,2]) )
	max.len <- max( sapply(ndigits, length) )
	L <- do.call('rbind', lapply(ndigits, function(z){ length(z) <- max.len ; z} ) )
	nd <- apply( L, 2, max, na.rm=T )
	for( i in seq(ng) ){
		this <- d2[[i]] 
		w.adj <- strwidth('0', cex=cex.axis, font=par('font')) / 2 * nd[1:nrow(this)]
		mtext( side=1, at=this$time+w.adj, text=this$n.risk, line=line.pos[i], cex=cex.axis, adj=1, col=1, las=1)
				   }
	if(label.n.at.risk) mtext( side=1, text='N at risk', at=group.name.pos, line=1.5, adj=1, col=1, las=1, cex=cex.axis )
	} ## End of if(!simple)	

	# Legend
	rlp <- group.order
	if(legend){
		bgc <- ifelse( par('bg')=='transparent', 'white', par('bg') )
	legend(x=loc.legend, legend=group.names[rlp], col=col.surv[rlp], lty=lty.surv[rlp], lwd=lwd.surv[rlp],
		bty='o', cex=cex.axis, bg=bgc, box.col='transparent', inset=.01 )
			}

	## draw confidence intervals
	for(i in 1:ng){
		this <- dat.list[[i]] 
			x <- this$time ; L <- this$lower ; U <- this$upper ; S <- this$survival
				naL <- which( is.na(L) ) ; L[naL] <- L[naL-1] ; U[naL] <- U[naL-1]
		lines( x, L, type='s', col=col.ci[i], lty=lty.ci[i], lwd=lwd.ci[i] )
		lines( x, U, type='s', col=col.ci[i], lty=lty.ci[i], lwd=lwd.ci[i] )
				}
	# draw curves
	lines(km, conf.int=FALSE, col=col.surv, lty=lty.surv, lwd=lwd.surv, mark=mark )

	box(bty=par('bty'))
	sdiff <- survdiff(eval(km$call$formula), data = eval(km$call$data))
    pval <- pchisq(sdiff$chisq,length(sdiff$n) - 1,lower.tail = FALSE)
    pvaltxt <- ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))
    mtext(pvaltxt)
	# par(op) 
}
