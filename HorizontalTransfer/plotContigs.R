#--------------------------------------------------------------------------------------------------------
#
# PURPOSE: produce a plot detailing the transferred regions for a given clone
#
# VERSION: 1
#
# COPYRIGHT: Nicolas Guex & Christian Iseli, SIB - Swiss Institute of Bioinformatics
#
#
# DEPENDENCES: library('circlize')
# INPUT:       list of contigs (.lst) covering donor and acceptor, and transfered (.trf and .lim) elements
# OUTPUT:      pdf file with plot, and its associated text file containing the retained transferred elements
#
#
# REMARK: this script is not generalized
#         (e.g. heavily specific for the figure, with several hardcoded values and heavy use of global variables)
#         It is provided solely with the intent to be able to reproduce the figure using the associated script cmd.sh
#
#--------------------------------------------------------------------------------------------------------

PlotAcceptorCoveredRegions <- function(dataS,fillcolor,addBoundaryLabel)
{

	if (nrow(dataS) > 0)
	{
		delta <- dataS[2:nrow(dataS),2] - dataS[1:(nrow(dataS)-1),3]
		prevrightFlank <- -100000
		for (i in 1:nrow(dataS))
		{
			leftFlank  <- dataS[i,2]
			rightFlank <- dataS[i,3]
			nextleftFlank <- dataS[i+1,2]
			leftlabeloffset <-  leftFlank - prevrightFlank
			rightlabeloffset <-  rightFlank - nextleftFlank
			if (is.na(rightlabeloffset)) { rightlabeloffset <- 100000}
			if (rightlabeloffset < 2000) { rightlabeloffset <- 2000 } else { rightlabeloffset <- 0 }
			if (leftlabeloffset < 2000) { leftlabeloffset <- 2000 } else { leftlabeloffset <- 0 }

				circos.rect(leftFlank ,-halftrackheight ,rightFlank ,  halftrackheight + 0.0 , col=fillcolor,border='NA')
					if ((rightFlank - leftFlank) < 8000)
					{
						if (addBoundaryLabel)
						{
							circos.text(leftFlank+leftlabeloffset  ,lblPos,   paste(leftFlank,rightFlank,sep="- ")  ,facing='clockwise',niceFacing=TRUE,cex=htpos.fontsize,adj=c(1.0,0.5),col='black')
						}
					} else {
						if (addBoundaryLabel && leftFlank >= 1000) # avoid plotting label near origin
						{
							circos.text(leftFlank+leftlabeloffset  ,lblPos,   leftFlank  ,facing='clockwise',niceFacing=TRUE,cex=htpos.fontsize,adj=c(1.0,0.5),col='black')
							circos.segments(leftFlank  ,0.5,   leftFlank+leftlabeloffset , 1.0 , straight=TRUE, col='black')
						}
						if (addBoundaryLabel && rightFlank <= template.len-1000)
						{
							circos.text(rightFlank-rightlabeloffset  ,lblPos,   rightFlank  ,facing='clockwise',niceFacing=TRUE,cex=htpos.fontsize,adj=c(1.0,0.5),col='black')
							circos.segments(rightFlank  ,0.5,   rightFlank-rightlabeloffset , 1.0 , straight=TRUE, col='black')
						}
					}
			prevrightFlank <- rightFlank
		}
	}
}
#-----------------------------------------------------------------------------------

PlotDonorCoveredRegions <- function(dataS,fillcolor,addBoundaryLabel,chr)
{
	color <- 'darkgreen'
	minOffset <- 3500
	summary <- NULL
	if (nrow(dataS) > 0)
	{
		dataS <- dataS[dataS$len >0,]
		delta <- dataS[2:nrow(dataS),2] - dataS[1:(nrow(dataS)-1),3]
		prevrightFlank <- -100000
		for (i in 1:nrow(dataS))
		{
			leftFlank  <- dataS[i,2] + dataS[i,'fistSnpOffset']
			rightFlank <- leftFlank + dataS[i,'len']
			nextleftFlank <- dataS[i+1,2] + dataS[i+1,'fistSnpOffset']
			leftlabeloffset <-  leftFlank - prevrightFlank
			rightlabeloffset <-  rightFlank - nextleftFlank
			if (is.na(rightlabeloffset)) { rightlabeloffset <- 100000}
			if (rightlabeloffset < minOffset) { rightlabeloffset <- minOffset } else { rightlabeloffset <- 0 }
			if (leftlabeloffset < minOffset) { leftlabeloffset <- minOffset } else { leftlabeloffset <- 0 }

			circos.rect(leftFlank ,-halftrackheight ,rightFlank ,  halftrackheight + 0.0 , col=fillcolor,border=fillcolor)
			boundaries <- cbind(sample=sample,chr=chr,leftFlank=leftFlank,rightFlank=rightFlank,len=rightFlank-leftFlank+1)
			summary <- rbind(summary,boundaries)
			
			if (dataS[i,'len'] < 8000)
			{
				if (addBoundaryLabel)
				{
					circos.text(leftFlank+leftlabeloffset  ,lblPos,   paste(leftFlank,rightFlank,sep="- ")  ,facing='clockwise',niceFacing=TRUE,cex=htpos.fontsize,adj=c(1.0,0.5),col=color)
				}
			} else {
				if (addBoundaryLabel && leftFlank >= 1000) # avoid plotting label near origin
				{
					circos.text    (leftFlank+leftlabeloffset  ,lblPos,   leftFlank  ,facing='clockwise',niceFacing=TRUE,cex=htpos.fontsize,adj=c(1.0,0.5),col=color)
					circos.segments(leftFlank  ,0.5,   leftFlank+leftlabeloffset , 1.0 , straight=TRUE, col=color)
				}
				if (addBoundaryLabel && rightFlank <= template.len-1000)
				{
					circos.text    (rightFlank-rightlabeloffset  ,lblPos,   rightFlank  ,facing='clockwise',niceFacing=TRUE,cex=htpos.fontsize,adj=c(1.0,0.5),col=color)
					circos.segments(rightFlank  ,0.5,   rightFlank-rightlabeloffset , 1.0 , straight=TRUE, col=color)
				}
			}
			prevrightFlank <- rightFlank
		}
	}
	return(summary)
}
#-----------------------------------------------------------------------------------

PlotNonOriginTransferBoundaries <- function(data)
{
	prevrightFlank <- -100000
	minOffset <- 3500
	
	if (nrow(data) > 0)
	{

		idx <- sort(data[,"acceptorFirstSNP"],index.return=T)$ix
		data <- data[idx,]
		data <- data[data$len >0,]

		for (i in 1:nrow(data))
		{
			color <- 'darkgreen'
			
			leftFlank  <- data[i,"acceptorFirstSNP"]
			rightFlank <- data[i,"acceptorLastSNP"]
			nextleftFlank <- data[i+1,"acceptorFirstSNP"]
			leftlabeloffset <-  leftFlank - prevrightFlank
			rightlabeloffset <-  rightFlank - nextleftFlank
			if (is.na(rightlabeloffset)) { rightlabeloffset <- 100000}

			if (rightlabeloffset < minOffset) { rightlabeloffset <- minOffset } else { rightlabeloffset <- 0 }
			if (leftlabeloffset < minOffset) { leftlabeloffset <- minOffset } else { leftlabeloffset <- 0 }

			circos.rect(leftFlank ,-halftrackheight ,rightFlank ,  0.0 , col=color,border='NA')

			if ((rightFlank - leftFlank) < 8000)
			{
				circos.text(leftFlank  ,lblPos,   paste(leftFlank,rightFlank,sep="- ")  ,facing='clockwise',niceFacing=TRUE,cex=htpos.fontsize,adj=c(1.0,0.5),col=color)
			} else
			{
				circos.text(leftFlank+leftlabeloffset  ,lblPos,   leftFlank  ,facing='clockwise',niceFacing=TRUE,cex=htpos.fontsize,adj=c(1.0,0.5),col=color)
				circos.segments(leftFlank  ,0.5,   leftFlank+leftlabeloffset , 1.0 , straight=TRUE, col=color)

				circos.text(rightFlank-rightlabeloffset  ,lblPos, rightFlank  ,facing='clockwise',niceFacing=TRUE,cex=htpos.fontsize,adj=c(1.0,0.5),col=color)
				circos.segments(rightFlank  ,0.5,   rightFlank-rightlabeloffset , 1.0 , straight=TRUE, col=color)
			}
			prevrightFlank <- rightFlank
		}
	}
}
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#
#                                   MAIN
#
#  NOTE:  expects input   in a directory T?/
#                 outputs in a directory plots_T?/
#
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

options(scipen=13)
options(width=150)
library('circlize')

args=commandArgs(TRUE)
sample <- args[1]
whichTemplate <- as.character(args[2])
SNPratio <- as.numeric(args[3])
TEMPLATE <- "0"

if (identical(whichTemplate, "1"))
{
	TEMPLATE <- "1"
	template.Chr1.len <- 2955400
	template.Chr2.len <- 1096954
	RESISTANCEA2.present <- F
	RESISTANCE1.present <- F
	RESISTANCE2.present <- T
	RESISTANCE2.start <- 117375
	RESISTANCE2.end <- 118852
	RESISTANCE2.text <- 'VCA0107::FRT-Kan-FRT'
}

if (identical(whichTemplate, "2"))
{
	TEMPLATE <- "2"
	template.Chr1.len <- 2955400
	template.Chr2.len <- 1096951
	RESISTANCEA2.present <- F
	RESISTANCE1.present <- F
	RESISTANCE2.present <- T
	RESISTANCE2.start <- 721147
	RESISTANCE2.end <- 722619
	RESISTANCE2.text <- 'VCA0747::FRT-Kan-FRT'
}

if (identical(whichTemplate, "A"))
{
	TEMPLATE <- "A"
	template.Chr1.len <- 2955400
	template.Chr2.len <- 1097967
	RESISTANCEA2.present <- F
	RESISTANCE1.present <- F
	RESISTANCE2.present <- T
	RESISTANCE2.start <- c(117375,722163)
	RESISTANCE2.end <- c(118387,723635)
	RESISTANCE2.text <- c('VCA0107::FRT-Cm-FRT','VCA0747::FRT-Kan-FRT')
}

if (identical(whichTemplate, "B"))
{
	TEMPLATE <- "B"
	template.Chr1.len <- 2955810
	template.Chr2.len <- 1097967
	RESISTANCEA2.present <- T
	RESISTANCEA2.start <- 119187
	RESISTANCEA2.end <- 120199
	RESISTANCEA2.text <- 'VCA0107::FRT-Cm-FRT'
	RESISTANCE1.present <- T
	RESISTANCE1.start <- 1062346
	RESISTANCE1.end <- 1063345
	RESISTANCE1.text <- 'VC1879-Amp'
	RESISTANCE2.present <- T
	RESISTANCE2.start <- c(117375,722163)
	RESISTANCE2.end <- c(118387,723635)
	RESISTANCE2.text <- c('VCA0107::FRT-Cm-FRT','VCA0747::FRT-Kan-FRT')
}

if (identical(whichTemplate, "C1") || identical(whichTemplate, "C2"))
{
	TEMPLATE <- "C"
	template.Chr1.len <- 2956374
	template.Chr2.len <- 1095478
	RESISTANCEA2.present <- F
	RESISTANCE1.present <- T
	RESISTANCE1.start <- c(592268,1062910)
	RESISTANCE1.end <- c(593746,1063909)
	RESISTANCE1.text <- c('VC2338::FRT-Kan-FRT','VC1879-Amp')
	RESISTANCE2.present <- F
}

if (identical(whichTemplate, "D"))
{
	TEMPLATE <- "D"
	template.Chr1.len <- 2956374
	template.Chr2.len <- 1096494
	RESISTANCEA2.present <- F
	RESISTANCE1.present <- T
	RESISTANCE1.start <- c(592268,1062910)
	RESISTANCE1.end <- c(593746,1063909)
	RESISTANCE1.text <- c('VC2338::FRT-Kan-FRT','VC1879-Amp')
	RESISTANCE2.present <- T
	RESISTANCE2.start <- 117375
	RESISTANCE2.end <- 118387
	RESISTANCE2.text <- 'VCA0107::FRT-Cm-FRT'
}

if (identical(whichTemplate, "E"))
{
	TEMPLATE <- "E"
	template.Chr1.len <- 2955400
	template.Chr2.len <- 1096954
	RESISTANCEA2.present <- F
	RESISTANCE1.present <- F
	RESISTANCE2.present <- T
	RESISTANCE2.start <- 117375
	RESISTANCE2.end <- 118850
	RESISTANCE2.text <- 'VCA0107::FRT-Kan-FRT'
}

if (identical(TEMPLATE, "0"))
{
	print("Unknown template")
	quit(save="no",status=1,runLast=TRUE)
}



template.color <- 'white'
template.bordercolor <- 'red'
halftrackheight <- 0.33
axisFontSize <- 1.333

testsize <- 3.0
lblPos <- 1.9

# set acceptor chr lengths

Chr1.len  <- 3015094
Chr2.len  <- 1070371		#  was 1070374  before correction for sequencing errors.
if (identical(whichTemplate, "B"))
{
	Chr2.len  <- 1071387   # unlike all other acceptor, this one contains an insertion K7
}


resistance.textsize <- 2.0
resistance.color <- 'skyblue'

plotStart <- 1.2  # do not change
templateLen.fontsize <- 1.1
htpos.fontsize <- 1.0


# --- read contigs

contigs <- read.table(paste('T',whichTemplate,'/',sample,'.lst',sep=''),header=F,sep="\t")
colnames(contigs) <- c('V1','V2','chr','contigStart','contigLen','contigEnd')
transfer <- read.delim(paste('T',whichTemplate,'/',sample,'.trf',sep=''),header=F,sep=",",colClasses=c('character','character','numeric','numeric','numeric','numeric','numeric','numeric','numeric','character','character','numeric','character','character'))
ignored <- which(transfer$V7 == -1 | is.na(transfer$V7))
ignoredStartPos <- NULL
if (length(ignored > 0)) {
	ignoredStartPos <- transfer$V3[ignored]
	transfer <- transfer[-ignored,]
}
minimaltransfer <- read.delim(paste('T',whichTemplate,'/',sample,'.lim',sep=''),header=F,sep="\t",colClasses=c('numeric','character','numeric','numeric','numeric','numeric','numeric','numeric','numeric'))
colnames(minimaltransfer) <- c('sample','chr','donorAlignStart','donorAlignLen','acceptorFirstSNP','acceptorLastSNP','len','acceptorAlignStart','fistSnpOffset')

tbl <- table(minimaltransfer$donorAlignStart == transfer$V3 & minimaltransfer$donorAlignLen == transfer$V4)
if (length(tbl) != 1) { print(paste("MERGE CHR1 ERROR",sample)) ; print(transfer) ; print(minimaltransfer) ; quit() }

transfer <- cbind(transfer,minimaltransfer[,5:9])
transfer1 <- transfer[which((transfer$V4 / transfer$V7 < SNPratio | transfer$V8 >= 100 )  & transfer$V2 == "ch1"),]
transfer2 <- transfer[which((transfer$V4 / transfer$V7 < SNPratio | transfer$V8 >= 100)  & transfer$V2 == "ch2"),]


dataT1 <- contigs[grep('Don',contigs$V1),]
dataT1 <- dataT1[dataT1$V2=='Don',]


if (length(ignoredStartPos) > 0) { print(dataT1); print(paste("Ignoring dataT1 startpos",ignoredStartPos,"Of Sample",sample)) ;  dataT1 <- dataT1[-which(dataT1$contigStart %in% ignoredStartPos),] ; print(dataT1); }
tbl <- table(minimaltransfer$donorAlignStart == dataT1$contigStart & minimaltransfer$donorAlignLen == dataT1$contigLen)
if (length(tbl) != 1) { print(paste("MERGE dataT1 ERROR",sample)) ; print(contigs) ; print(dataT1) ; print(minimaltransfer) ; quit() }
dataT1 <- cbind(dataT1,minimaltransfer[,7:9])
dataT1 <- dataT1[which(dataT1$chr=='ch1'),c(3,4,6,7,8,9)]
dataT1 <- dataT1[which(dataT1$contigStart %in% transfer1$V3),]


dataA1 <- contigs[grep('Acc',contigs$V1),]
dataA1 <- dataA1[dataA1$V2=='Acc',]
dataA1 <- dataA1[which(dataA1$chr=='ch1'),c(3,4,6)]

dataT2 <- contigs[grep('Don',contigs$V1),]
dataT2 <- dataT2[dataT2$V2=='Don',]

if (length(ignoredStartPos) > 0) { print(dataT2); print(paste("Ignoring dataT2 startpos",ignoredStartPos,"Of Sample",sample)) ;  dataT2 <- dataT2[-which(dataT2$contigStart %in% ignoredStartPos),] ; print(dataT2); }

tbl <- table(minimaltransfer$donorAlignStart == dataT2$contigStart & minimaltransfer$donorAlignLen == dataT2$contigLen)
if (length(tbl) != 1) { print(paste("MERGE dataT2 ERROR",sample)) ; quit() }

dataT2 <- cbind(dataT2,minimaltransfer[,7:9])
dataT2 <- dataT2[which(dataT2$chr=='ch2'),c(3,4,6,7,8,9)]
dataT2 <- dataT2[which(dataT2$contigStart %in% transfer2$V3),]


dataA2 <- contigs[grep('Acc',contigs$V1),]
dataA2 <- dataA2[dataA2$V2=='Acc',]
dataA2 <- dataA2[which(dataA2$chr=='ch2'),c(3,4,6)]

# --- save retained transferred events

retained <- rbind(transfer1,transfer2)
retained[,10] <- paste(retained[,10],'>',retained[,11])
retained[,13] <- paste(retained[,13],'>',retained[,14])
retained <- retained[,-c(1,11,14)]
colnames(retained) <- c('chr','DonorStart','DonorLen','DonorSNPCnt','DonorIndelLen','AcceptorSNPCnt','AcceptorIndelLen','AcceptorStartPos','AcceptorStartVCF','AcceptorEndPos','AcceptorEndVCF')
retained <- cbind(sample=sample,retained)
retained <- cbind(retained,TransferLen=retained$AcceptorEndPos - retained$AcceptorStartPos + 1)

# --- plot

retainedTransfers <- NULL
pdf(file=paste('plots_T',whichTemplate,'/',sample,'_horizontal.pdf',sep=''),width=30,height=30)
layout(matrix(1:4,ncol=2,byrow=T))


	# -- Chr 1 of A1552


	template.len <- Chr1.len

	margins <- c(-1.0,1.0)
	par(mar=c(2,2,2,2))
	circos.clear()
	circos.par("start.degree" = 90, canvas.xlim = margins, canvas.ylim = margins,
			gap.degree = 0, cell.padding = c(0.0, 0.0, 0.0, 0.0), track.margin = c(0.001, 0.001),track.margin = c(0.001, 0.001), track.height = 0.1)
	circos.initialize(factors = 'a', xlim = c(template.len,1))
	circos.par(track.margin = c(0, convert_height(10, "mm")))
	circos.track(ylim = c(0, 1), track.height = convert_height(20, "mm"),bg.border='white')

	text(0.0,  0.3, sample, adj = c(0.5,0.5) , cex=testsize*1.5)
	text(0.0,  0.15, "Chromosome I", adj = c(0.5,0.5) , cex=testsize)
	text(0.0,  0.05, paste("A1552"), adj = c(0.5,0.5),  , cex=testsize)
	text(0.0,  -0.03, "based on CP028894 (A1552)", adj = c(0.5,0.5), , cex=0.5*testsize)
	circos.rect(1  ,-halftrackheight , template.len, -halftrackheight-0.03 ,col=template.bordercolor,border='NA')
	circos.rect(1  ,halftrackheight , template.len, halftrackheight+0.03 ,col=template.bordercolor,border='NA')

	idx <- NULL
	if (nrow(transfer1) > 0)
	{
		for (i in 1: nrow(dataA1))
		{
			for (j in 1: nrow(transfer1))
			{
				if (dataA1[i,2] >= transfer1$acceptorAlignStart[j] && dataA1[i,3] <= transfer1$acceptorAlignStart[j]+transfer1$len[j])
				{
					idx <- c(idx,i)
				}
			}
		}
	}

	if (length(idx) > 0)
	{
		print(paste("REMOVE",length(idx),"FROM ACCEPTOR chr1",sample))
		dataA1 <- dataA1[-idx,]
	}

	idx <- NULL
	if (nrow(transfer1) > 0)
	{
		for (i in 1: nrow(transfer1))
		{
			for (j in 1: nrow(dataA1))
			{
				if (transfer1$acceptorAlignStart[i] >= dataA1[j,2] && transfer1$acceptorAlignStart[i]+transfer1$fistSnpOffset[i]+transfer1$len[i] <= dataA1[j,3])
				{
					idx <- c(idx,i)
				}
			}
		}
	}

	if (length(idx) > 0)
	{
		print(paste("REMOVE",length(idx),"FROM DONOR chr1",sample))
		transfer1 <- transfer1[-idx,]
	}

	PlotAcceptorCoveredRegions(dataA1,'red',FALSE)


	# plot origin spanning and mark it as done by putting a len of zero.

	idx1 <- which(transfer1[,"acceptorLastSNP"] >= template.len - 1000)
	idx2 <- which(transfer1[,"acceptorFirstSNP"] <= 1000)
	if (length(idx1) > 0 && length(idx2 > 0))
	{
		leftFlank <- transfer1[idx1,"acceptorFirstSNP"]
		circos.rect(leftFlank ,-halftrackheight ,template.len ,  0.0 , col='darkgreen',border='NA')
		circos.rect(1 ,-halftrackheight ,  transfer1[idx2,"acceptorLastSNP"]  , 0.0, col='darkgreen',border='NA')
		circos.text(leftFlank+5000  ,lblPos,   paste(leftFlank,transfer1[idx2,"acceptorLastSNP"],sep="- ")  ,facing='clockwise',niceFacing=TRUE,cex=htpos.fontsize,adj=c(1.0,0.5),col='darkgreen')
		transfer1[idx1,"len"] <- 0
		transfer1[idx2,"len"] <- 0
	}

	PlotNonOriginTransferBoundaries(transfer1)

	circos.track(ylim = c(0, 1), track.height = convert_height(1, "mm"),bg.border='white')
	first_label_height <- 1 ;last_label_height <-1 ; circos.axis(major.at = seq(template.len,1,by=-50000),labels.cex=axisFontSize,labels.facing='clockwise' , direction='inside')
	lblseq <- seq(5000,template.len,by=50000);circos.text(lblseq  ,rep(-1.0,length(lblseq)),  lblseq ,facing='clockwise',niceFacing=TRUE,cex=axisFontSize,adj=c(1.0,0.5),col='black')
	for (l in lblseq) { circos.segments(l  ,0.0,   l , 1.0 , straight=TRUE, col='black') }
	lblseq <- seq(15000,template.len,by=50000); for (l in lblseq) { circos.segments(l  ,0.0,   l , 0.5, straight=TRUE, col='black') }
	lblseq <- seq(25000,template.len,by=50000); for (l in lblseq) { circos.segments(l  ,0.0,   l , 0.5, straight=TRUE, col='black') }
	lblseq <- seq(35000,template.len,by=50000); for (l in lblseq) { circos.segments(l  ,0.0,   l , 0.5, straight=TRUE, col='black') }
	lblseq <- seq(45000,template.len,by=50000); for (l in lblseq) { circos.segments(l  ,0.0,   l , 0.5, straight=TRUE, col='black') }


	# -- Chr 2 of A1552


	template.len <- Chr2.len
	margins <- c(-1.33,1.33 )
	par(mar=c(5,5,5,5))
	circos.clear()
	circos.par("start.degree" = 90, canvas.xlim = margins, canvas.ylim = margins,
			gap.degree = 0, cell.padding = c(0.0, 0.0, 0.0, 0.0), track.margin = c(0.001, 0.001),track.margin = c(0.001, 0.001), track.height = 0.1)
	circos.initialize(factors = 'a', xlim = c(plotStart,template.len))
	circos.par(track.margin = c(0, convert_height(10, "mm")))
	circos.track(ylim = c(0, 1), track.height = convert_height(20, "mm"),bg.border='white')

	text(0.0,  0.3, sample, adj = c(0.5,0.5) , cex=testsize*1.5)
	text(0.0,  0.15, "Chromosome II", adj = c(0.5,0.5) , cex=testsize)
	text(0.0,  0.05, paste("A1552"), adj = c(0.5,0.5),  , cex=testsize)
	text(0.0,  -0.04, "based on CP028895 (A1552)", adj = c(0.5,0.5), , cex=0.5*testsize)
	circos.rect(1  ,-halftrackheight , template.len, -halftrackheight-0.03 ,col=template.bordercolor,border='NA')
	circos.rect(1  ,halftrackheight , template.len, halftrackheight+0.03 ,col=template.bordercolor,border='NA')

	idx <- NULL
	if (nrow(transfer2) > 0)
	{
		for (i in 1: nrow(dataA2))
		{
			for (j in 1: nrow(transfer2))
			{
				if (dataA2[i,2] >= transfer2$acceptorAlignStart[j] && dataA2[i,3] <= transfer2$acceptorAlignStart[j]+transfer2$len[j])
				{
					idx <- c(idx,i)
				}
			}
		}
	}

	if (length(idx) > 0)
	{
		print(paste("REMOVE",length(idx),"FROM ACCEPTOR chr2",sample))
		dataA2 <- dataA2[-idx,]
	}

	idx <- NULL
	if (nrow(transfer2) > 0)
	{
		for (i in 1: nrow(transfer2))
		{
			for (j in 1: nrow(dataA2))
			{
				if (transfer2$acceptorAlignStart[i] >= dataA2[j,2] && transfer2$acceptorAlignStart[i]+transfer2$fistSnpOffset[i]+transfer2$len[i] <= dataA2[j,3])
				{
					idx <- c(idx,i)
				}
			}
		}
	}

	if (length(idx) > 0)
	{
		print(paste("REMOVE",length(idx),"FROM DONOR chr2",sample))
		transfer2 <- transfer2[-idx,]
	}

	PlotAcceptorCoveredRegions(dataA2,'red',FALSE)

	# plot origin spanning and mark it as done by putting a len of zero.

	idx1 <- which(transfer2[,"acceptorLastSNP"] >= template.len - 1000)
	idx2 <- which(transfer2[,"acceptorFirstSNP"] <= 1000)
	if (length(idx1) > 0 && length(idx2 > 0))
	{
		leftFlank <- transfer2[idx1,"acceptorFirstSNP"]
		circos.rect(leftFlank ,-halftrackheight ,template.len ,  0.0 , col='darkgreen',border='NA')
		circos.rect(1 ,-halftrackheight ,  transfer2[idx2,"acceptorLastSNP"]  , 0.0, col='darkgreen',border='NA')
		circos.text(leftFlank+5000  ,lblPos,   paste(leftFlank,transfer2[idx2,"acceptorLastSNP"],sep="- ")  ,facing='clockwise',niceFacing=TRUE,cex=htpos.fontsize,adj=c(1.0,0.5),col='darkgreen')
		transfer2[idx1,"len"] <- 0
		transfer2[idx2,"len"] <- 0
	}

	PlotNonOriginTransferBoundaries(transfer2)

	if (RESISTANCEA2.present) {
	  for (i in 1:length(RESISTANCEA2.start)) {
	    circos.rect(RESISTANCEA2.start[i], -halftrackheight , RESISTANCEA2.end[i], halftrackheight ,col=resistance.color,border=resistance.color)
	    circos.text(RESISTANCEA2.start[i]  ,halftrackheight+1.8,   RESISTANCEA2.text[i]  ,facing='inside',niceFacing=TRUE,cex=resistance.textsize,adj=c(0.5,0.5),col=resistance.color)
	  }
	}
	circos.track(ylim = c(0, 1), track.height = convert_height(1, "mm"),bg.border='white')
	circos.axis(major.at = seq(0,template.len,by=50000),labels.cex=axisFontSize,labels.facing='clockwise' , direction='inside')


	# -- Chr 1 of Template

	template.color <- 'white'
	template.bordercolor <- 'darkgreen'

	template.len <- template.Chr1.len

	margins <- c(-1.0,1.0)
	par(mar=c(2,2,2,2))
	circos.clear()
	circos.par("start.degree" = 90, canvas.xlim = margins, canvas.ylim = margins,
			gap.degree = 0, cell.padding = c(0.0, 0.0, 0.0, 0.0), track.margin = c(0.001, 0.001),track.margin = c(0.001, 0.001), track.height = 0.1)
	circos.initialize(factors = 'a', xlim = c(1,template.len))
	circos.par(track.margin = c(0, convert_height(10, "mm")))
	circos.track(ylim = c(0, 1), track.height = convert_height(20, "mm"),bg.border='white')

	text(0.0,  0.3, sample, adj = c(0.5,0.5) , cex=testsize*1.5)
	text(0.0,  0.15, "Chromosome I", adj = c(0.5,0.5) , cex=testsize)
	text(0.0,  0.05, paste("Template",TEMPLATE), adj = c(0.5,0.5),  , cex=testsize)
	text(0.0,  -0.03, "based on CP028892 (Sa5Y)", adj = c(0.5,0.5), , cex=0.5*testsize)
	circos.rect(1  ,-halftrackheight , template.len, -halftrackheight-0.03 ,col=template.bordercolor,border='NA')
	circos.rect(1  ,halftrackheight , template.len, halftrackheight+0.03 ,col=template.bordercolor,border='NA')

	dataT1 <- dataT1[which(dataT1$contigStart %in% transfer1$V3),]

	# plot origin spanning and mark it as done by putting a len of zero.

	idx1 <- which(dataT1[,"contigEnd"] >= template.len - 1000)
	idx2 <- which(dataT1[,"contigStart"] <= 1000)
	if (length(idx1) > 0 && length(idx2 > 0))
	{
		leftFlank <- dataT1[idx1,"contigStart"] + dataT1[idx1,"fistSnpOffset"]
		rightFlank <- dataT1[idx2,"contigStart"] + dataT1[idx2,"fistSnpOffset"] + dataT1[idx2,"len"]
		circos.rect(leftFlank ,-halftrackheight ,template.len ,  halftrackheight , col='darkgreen',border='NA')
		circos.rect(1 ,-halftrackheight ,  rightFlank  , halftrackheight, col='darkgreen',border='NA')
		circos.text(leftFlank+5000  ,lblPos,   paste(leftFlank,rightFlank,sep="- ")  ,facing='clockwise',niceFacing=TRUE,cex=htpos.fontsize,adj=c(1.0,0.5),col='darkgreen')
		dataT1[idx1,"len"] <- 0
		dataT1[idx2,"len"] <- 0
		boundaries <- cbind(sample=sample,chr='chr1',leftFlank=leftFlank,rightFlank=rightFlank,len=(template.len-leftFlank)+(rightFlank-1)+1)
		retainedTransfers <- rbind(retainedTransfers,boundaries)
	}

	retainedTransfers <- rbind(retainedTransfers,PlotDonorCoveredRegions(dataT1,'darkgreen',TRUE,'chr1'))

	# plot resistance K7

	if (RESISTANCE1.present) {
	  for (i in 1:length(RESISTANCE1.start)) {
	    circos.rect(RESISTANCE1.start[i], -halftrackheight , RESISTANCE1.end[i], halftrackheight ,col=resistance.color,border=resistance.color)
	    circos.text(RESISTANCE1.start[i]  ,halftrackheight+1.8,   RESISTANCE1.text[i]  ,facing='inside',niceFacing=TRUE,cex=resistance.textsize,adj=c(0.5,0.5),col=resistance.color)
	  }
	}
	circos.track(ylim = c(0, 1), track.height = convert_height(5, "mm"),bg.border='white')
	circos.axis(major.at = seq(0,template.len,by=50000),labels.cex=axisFontSize,labels.facing='clockwise' , direction='inside')

	# -- Chr 2  of Template

	template.len <- template.Chr2.len

	margins <- c(-1.33,1.33 )
	par(mar=c(5,5,5,5))
	circos.clear()
	circos.par("start.degree" = 90, canvas.xlim = margins, canvas.ylim = margins,
			gap.degree = 0, cell.padding = c(0.0, 0.0, 0.0, 0.0), track.margin = c(0.001, 0.001),track.margin = c(0.001, 0.001), track.height = 0.1)
	circos.initialize(factors = 'a', xlim = c(plotStart,template.len))
	circos.par(track.margin = c(0, convert_height(10, "mm")))
	circos.track(ylim = c(0, 1), track.height = convert_height(20, "mm"),bg.border='white')

	text(0.0,  0.3, sample, adj = c(0.5,0.5) , cex=testsize*1.5)
	text(0.0,  0.15, "Chromosome II", adj = c(0.5,0.5) , cex=testsize)
	text(0.0,  0.05, paste("Template",TEMPLATE), adj = c(0.5,0.5),  , cex=testsize)
	text(0.0,  -0.04, "based on CP028893 (Sa5Y)", adj = c(0.5,0.5), , cex=0.5*testsize)
	circos.rect(1  ,-halftrackheight , template.len, -halftrackheight-0.03 ,col=template.bordercolor,border='NA')
	circos.rect(1  ,halftrackheight , template.len, halftrackheight+0.03 ,col=template.bordercolor,border='NA')

	dataT2 <- dataT2[which(dataT2$contigStart %in% transfer2$V3),]

	# plot origin spanning and mark it as done by putting a len of zero.

	idx1 <- which(dataT2[,"contigEnd"] >= template.len - 1000)
	idx2 <- which(dataT2[,"contigStart"] <= 1000)
	if (length(idx1) > 0 && length(idx2 > 0))
	{
		leftFlank <- dataT2[idx1,"contigStart"] + dataT2[idx1,"fistSnpOffset"]
		rightFlank <- dataT2[idx2,"contigStart"] + dataT2[idx2,"fistSnpOffset"] + dataT2[idx2,"len"]
		circos.rect(leftFlank ,-halftrackheight ,template.len ,  halftrackheight , col='darkgreen',border='NA')
		circos.rect(1 ,-halftrackheight ,  rightFlank  , halftrackheight, col='darkgreen',border='NA')
		circos.text(leftFlank+5000  ,lblPos,   paste(leftFlank,rightFlank,sep="- ")  ,facing='clockwise',niceFacing=TRUE,cex=htpos.fontsize,adj=c(1.0,0.5),col='darkgreen')
		dataT2[idx1,"len"] <- 0
		dataT2[idx2,"len"] <- 0
		boundaries <- cbind(sample=sample,chr='chr2',leftFlank=leftFlank,rightFlank=rightFlank,len=(template.len-leftFlank)+(rightFlank-1)+1)
		retainedTransfers <- rbind(retainedTransfers,boundaries)
	}

	retainedTransfers <- rbind(retainedTransfers,PlotDonorCoveredRegions(dataT2,'darkgreen',TRUE,'chr2'))

	# plot resistance K7

	if (RESISTANCE2.present) {
	  for (i in 1:length(RESISTANCE2.start)) {
	    circos.rect(RESISTANCE2.start[i], -halftrackheight , RESISTANCE2.end[i], halftrackheight ,col=resistance.color,border=resistance.color)
	    circos.text(RESISTANCE2.start[i]  ,halftrackheight+1.8,   RESISTANCE2.text[i]  ,facing='inside',niceFacing=TRUE,cex=resistance.textsize,adj=c(0.5,0.5),col=resistance.color)
	  }
	}
	circos.track(ylim = c(0, 1), track.height = convert_height(1, "mm"),bg.border='white')
	circos.axis(major.at = seq(0,template.len,by=50000),labels.cex=axisFontSize,labels.facing='clockwise' , direction='inside')



dev.off()

write.table(retainedTransfers,file=paste('plots_T',whichTemplate,'/',sample,'_transferred.txt',sep=''),sep="\t",quote=F,row.names=F)


