#--------------------------------------------------------------------------------------------------------
#
# PURPOSE: produce a plot comparing two strains of v.cholerae
#
# VERSION: 4
#
# COPYRIGHT: Nicolas Guex & Christian Iseli, SIB - Swiss Institute of Bioinformatics
#
#
# DEPENDENCES: library('circlize')
# INPUT:       10kb chunks alignment results
# OUTPUT:      pdf file
#
#
# REMARK: this script is not generalized
#         (e.g. heavily specific for the figure, with several hardcoded values and heavy use of global variables)
#         It is provided solely with the intent to be able to reproduce the figure
#         To do this, simply type the following commands
#
#           R --vanilla --quiet < plot_A1552_vs_CP028827-v4.R
#
#--------------------------------------------------------------------------------------------------------

plotAnnot <- function(idx,offset)
{
	col.intensity <- 0.0
	circos.rect(outer.factor*(annot$startPos[idx]),rep(1.6+offset,length(idx)),outer.factor*(annot$endPos[idx]),rep(1.7+offset ,length(idx)),col=hsv(outer.col,1,1.0-col.intensity),border=hsv(outer.col,1,0.5))
	for (i in idx) { circos.text(outer.factor*(annot$startPos[i]+0.5*(annot$endPos[i]-annot$startPos[i])),1.9+offset,annot$NAME[i],cex=0.5*textsize) }
}
#--------------------------------------------------------------------------------------------------------

LoadFile <- function(fn,fn2,purge,filter)
{


	match <- read.table(fn,sep=",",header=F)
	colnames(match) <- c('log','orientation','chr','srcStart','srcEnd','dstStart','dstEnd','SNPcnt','delCnt','insCnt')

	idx <- which(match$insCnt >= 100000)
	if (length(idx) > 0)
	{
		match$SNPcnt[idx] <-   -1
	}
	

	match.v2 <- read.table(fn2,sep=",",header=F)
	colnames(match.v2) <- c('log','orientation','chr','srcStart','srcEnd','dstStart','dstEnd','SNPcnt','delCnt','insCnt')

	idx <- which(match.v2$insCnt >= 100000)
	if (length(idx) > 0)
	{
		match.v2$SNPcnt[idx] <-   -1
	}

	# keep best of two trials

	idx <- which((match.v2$SNPcnt < match$SNPcnt & match.v2$SNPcnt >= 0) |  (match.v2$SNPcnt != -1 & match$SNPcnt == -1))
	if (length(idx) > 0)
	{
		match[idx,] <- match.v2[idx,]
	}

	if (length(filter) > 0)
	{
		print('Filtering')
		print(match[filter,])
		match$SNPcnt[filter] <- -1
	}
	
	match <- match[match$SNPcnt != -1,]

	check <- which(diff(match$srcStart) == 0)
	if (purge==TRUE && length(check) > 0)
	{
		for (i in check)
		{
			if (match$SNPcnt[i] >= match$SNPcnt[i+1] ) { match$SNPcnt[i] = -1 } else { match$SNPcnt[i+1] = -1 }
		}
	}

	# split origin span

	match <- cbind(match,origin=0)
	idx <- which(match$dstEnd > (inner.len+512)) # do not split if this concern just a few nucleotides
	if (length(idx) == 1)
	{
		print('Splitting origin overlap')
		print(match[idx,])
		spillover <- -(inner.len-match$dstEnd[idx]+1)
		match <- rbind(match,match[idx,])
		match$dstEnd[idx] <- match$dstEnd[idx] - spillover
		match$srcEnd[idx] <- match$srcEnd[idx] - spillover
		match$origin[idx] <- 1

		match$srcStart[nrow(match)] <- match$srcEnd[idx] + 1
		match$dstStart[nrow(match)] <- 1
		match$dstEnd[nrow(match)] <- spillover
		match$origin[nrow(match)] <- 2
		print('split')
		print(match[idx,])
		print(match[nrow(match),])
	}

	return(match)

}
#--------------------------------------------------------------------------------------------------------
PlotLegend <- function()
{

	text(0.0,  -0.05, "diff. per 10Kb", adj = c(0.5,0.5), col='black' , cex=textsize/2)
	for (i in 0:20)
	{
		bottom <- -0.30+i*0.01
		rect(-0.05,bottom,0.00,bottom+0.01,col=hsv(outer.col,1,1.0-(0.05*i)),border=NA)
		rect(-0.00,bottom,0.05,bottom+0.01,col=hsv(inner.col,1,1.0-(0.05*i)),border=NA)
		if (i%%4 == 0 && i < 20) { text(0.09,  bottom+0.01,  (maxSNPcnt/20) * i , adj = c(0.5,0.5), col='black' , cex=textsize/3) }
		if (i == 20)             { text(0.09,  bottom+0.01,  paste('>=',(maxSNPcnt/20) * i ,sep='') , adj = c(0.5,0.5), col='black' , cex=textsize/3) }
	}

}

#--------------------------------------------------------------------------------------------------------
plot.inner.labels <- function(inner.labels,textsize,minDeltaPos)
{
	inner.labels <- inner.labels[inner.labels[,2] >=1,]
	idx <- sort(inner.labels[,1],index.return=T)$ix
	inner.labels <- inner.labels[idx,]
	print(inner.labels)
	lastplottedpos <- -1e10
	for (i in 1:nrow(inner.labels))
	{
		if ((inner.labels[i,1] - lastplottedpos) >= minDeltaPos) { circos.text(inner.factor*inner.labels[i,1],inner.start-0.3,   as.character( inner.labels[i,2] )  ,facing='clockwise',niceFacing=TRUE,cex=0.8*textsize,adj=c(1.0,0.5),col='grey30') ; lastplottedpos <- inner.labels[i,1] }
	}
}
#--------------------------------------------------------------------------------------------------------

plot.segment <- function(idx, STRAIGHT , PLOTPOSITIONS,textsize,segment.color)
{
	if (length(idx) == 0)
	{
		return(NULL)
	}
	col.intensity <- match$SNPcnt[idx] / maxSNPcnt
	col.intensity[col.intensity < 0.0] <- 0.0
	col.intensity[col.intensity > 1.0] <- 1.0
	circos.rect(outer.factor*match$srcStart[idx],rep(outer.start,length(idx)),outer.factor*match$srcEnd[idx],rep(outer.end,length(idx)),col=hsv(outer.col,1,1.0-col.intensity),border=hsv(outer.col,1,0.5))
	circos.rect(inner.factor*match$dstStart[idx],rep(inner.start,length(idx)),inner.factor*match$dstEnd[idx],rep(inner.end,length(idx)),col=hsv(inner.col,1,1.0-col.intensity),border=hsv(inner.col,1,0.5))

	ignore <- which(match$origin[idx] == 2)
	if (length(ignore) > 0)
	{
		print(match[idx[ignore],])
		idx <- idx[-ignore]
	}
	circos.segments(outer.factor*match$srcStart[idx], rep(outer.start,length(idx)), inner.factor*match$dstStart[idx], rep(inner.end,length(idx)), straight=STRAIGHT, col=segment.color)
	circos.segments(outer.factor*match$srcEnd[idx]  , rep(outer.start,length(idx)), inner.factor*match$dstEnd[idx]  , rep(inner.end,length(idx)), straight=STRAIGHT, col=segment.color)
	inner.labels <- NULL
	if (PLOTPOSITIONS == "REV")
	{
		circos.text(outer.factor*match$srcStart[idx],rep(outer.end + 0.8,length(idx)),  as.character( outer.len - match$srcStart[idx])   ,facing='clockwise',niceFacing=TRUE,cex=textsize,adj=c(1.0,0.5),col='grey30')
		inner.labels <- cbind(match$dstStart[idx],match$dstStart[idx])
	}
	if (PLOTPOSITIONS == "FWD")
	{
		circos.text(outer.factor*match$srcStart[idx],rep(outer.end + 0.8,length(idx)),   as.character( match$srcStart[idx])  ,facing='clockwise',niceFacing=TRUE,cex=textsize,adj=c(1.0,0.5),col='grey30')
		inner.labels <- cbind(match$dstStart[idx],match$dstStart[idx])
	}
	return(inner.labels)
}

#--------------------------------------------------------------------------------------------------------
# MAIN
#--------------------------------------------------------------------------------------------------------


library('circlize')
options(scipen=13)


maxSNPcnt <- 20
notComparableIndelLenThreshold <- 10000
minLabelDeltaChr1 <- 7000
minLabelDeltaChr2 <- 6000

pdf(file=paste('compareStrains.CP028827.lbl',minLabelDeltaChr1,'.',minLabelDeltaChr2,'.pdf',sep=''),width=30,height=15)
layout(matrix(1:2,ncol=2))

outer.col <- 1.00
inner.col <- 0.10

# -- Chr 1

outer.len <- 3015094
inner.len <- 2975504

outer.factor <- 1.0
if (outer.len >= inner.len)
{
	inner.factor <- outer.len / inner.len
} else {
	inner.factor <- inner.len / outer.len
}


match <- LoadFile('Align_A1552_onto_N16961.chr1.maxdel50000.txt','Align_A1552_onto_N16961.chr1.maxdel512.txt',TRUE,NULL)
match <- cbind(match,delta=match$insCnt - match$delCnt)

match$srcStart <- outer.len - match$srcStart
match$srcEnd <- outer.len - match$srcEnd
tmp <- match$srcEnd
match$srcEnd <- match$srcStart
match$srcStart <- tmp


annot <- read.table('coordinates_islands.txt',sep="\t",header=T)
idx <- which(annot$Chr == 1)
annot$startPos[idx] <- outer.len - annot$startPos[idx]
annot$endPos[idx] <- outer.len - annot$endPos[idx]
tmp <- annot$endPos[idx]
annot$endPos[idx] <- annot$startPos[idx]
annot$startPos[idx] <- tmp

par(mar=c(2,2,2,2))
margins <- c(-0.9,0.9)

circos.clear()
circos.par("start.degree" = 90, canvas.xlim = margins, canvas.ylim = margins,
        gap.degree = 0, cell.padding = c(0.0, 0.0, 0.0, 0.0), track.margin = c(0.001, 0.001),track.margin = c(0.001, 0.001), track.height = 0.1)
#circos.initialize(factors = 'a', xlim = c(1,max(outer.len,inner.len)))
circos.initialize(factors = 'a', xlim = c(1,outer.len))
circos.par(track.margin = c(0, convert_height(25, "mm")))
circos.track(ylim = c(0, 1), track.height = convert_height(20, "mm"),bg.border='white')

inner.start <- -1.5
inner.end <- -1.0

outer.start <- 0.0
outer.end   <- 0.5

textsize <- 2
positions.textsize <- 0.7


circos.rect(outer.factor*1,outer.start,outer.factor*outer.len,outer.end,col='white',border=hsv(outer.col,1,0.5))
circos.rect(inner.factor*1,inner.start,inner.factor*inner.len,inner.end,col='white',border=hsv(inner.col,1,0.5))


idx <- which(match$orientation == 'REV' & match$SNPcnt != -1 & abs(match$delta) < notComparableIndelLenThreshold)
inner.labels <- plot.segment( idx , TRUE ,'REV' ,positions.textsize,'black')
idx <- which(match$orientation == 'FWD' & match$SNPcnt != -1 & abs(match$delta) < notComparableIndelLenThreshold)
inner.labels <- rbind(inner.labels,plot.segment( idx , TRUE ,'REV' ,positions.textsize,'blue'))
plot.inner.labels(inner.labels,positions.textsize,minLabelDeltaChr1) 

plotAnnot(which(annot$Chr == 1 & annot$NAME != "LAT-1"),0.0)
plotAnnot(which(annot$Chr == 1 & annot$NAME == "LAT-1"),0.5)


text(0.0,  0.15, "I", adj = c(0.5,0.5) , cex=2.0*textsize)
text(0.0,  0.0, "CP028894 (A1552)", adj = c(0.5,0.5), col=hsv(outer.col,1,1.0) , cex=0.75*textsize)
text(0.0,  0.05, "CP028827 (N16961)", adj = c(0.5,0.5), col=hsv(inner.col,1,1.0) , cex=0.75*textsize)


PlotLegend()


# -- Chr 2

outer.len <- 1070374
inner.len <- 1072331

outer.factor <- 1.0
if (outer.len >= inner.len)
{
	inner.factor <- outer.len / inner.len
} else {
	inner.factor <- inner.len / outer.len
}

match <- LoadFile('Align_A1552_onto_N16961.chr2.maxdel50000.txt','Align_A1552_onto_N16961.chr2.maxdel512.txt',FALSE,NULL)
match <- cbind(match,delta=match$insCnt - match$delCnt)
match$dstStart[which(match$srcStart == 340001 & match$chr == 'chr2')] <- 341440  # patch because of gap in integron just in 10Kb junction

par(mar=c(5,5,5,5))
margins <- c(-1.1,1.1)
circos.clear()
circos.par("start.degree" = 90, canvas.xlim = margins, canvas.ylim = margins,
        gap.degree = 0, cell.padding = c(0.0, 0.0, 0.0, 0.0), track.margin = c(0.001, 0.001),track.margin = c(0.001, 0.001), track.height = 0.1)
circos.initialize(factors = 'a', xlim = c(1,outer.len))
circos.par(track.margin = c(0, convert_height(25, "mm")))
circos.track(ylim = c(0, 1), track.height = convert_height(20, "mm"),bg.border='white')

inner.start <- -1.5
inner.end <- -1.0

outer.start <- 0.0
outer.end   <- 0.5

textsize <- 2
positions.textsize <- 0.7

circos.rect(outer.factor*1,outer.start,outer.factor*outer.len,outer.end,col='white',border=hsv(outer.col,1,0.5))
circos.rect(inner.factor*1,inner.start,inner.factor*inner.len,inner.end,col='white',border=hsv(inner.col,1,0.5))


idx <- which(match$orientation == 'FWD' & match$SNPcnt != -1 & abs(match$delta) < notComparableIndelLenThreshold)
inner.labels <- plot.segment( idx , TRUE ,'FWD' ,positions.textsize,'black')
print(inner.labels)
plot.inner.labels(inner.labels,positions.textsize,minLabelDeltaChr2)

plotAnnot(which(annot$Chr == 2),0.0)

text(0.0,  0.15, "II", adj = c(0.5,0.5) , cex=2.0*textsize)
text(0.0,  0.0, "CP028895 (A1552)", adj = c(0.5,0.5), col=hsv(outer.col,1,1.0) , cex=0.75*textsize)
text(0.0,  0.05, "CP028828 (N16961)", adj = c(0.5,0.5), col=hsv(inner.col,1,1.0) , cex=0.75*textsize)


PlotLegend()

# ---

dev.off()

