#--------------------------------------------------------------------------------------------------------
#
# PURPOSE: produce a summary plot with the transferred regions of all clones of a given experiment
#
# VERSION: 1
#
# COPYRIGHT: Nicolas Guex & Christian Iseli, SIB - Swiss Institute of Bioinformatics
#
#
# DEPENDENCES: library('circlize')
# INPUT:       list of contigs (.lst) covering donor and acceptor, and transfered (.trf) elements
# OUTPUT:      pdf file with plot, and its associated text file containing the retained transferred elements
#
#
# REMARK: this script is not generalized
#         (e.g. heavily specific for the figure, with several hardcoded values and heavy use of global variables)
#         It is provided solely with the intent to be able to reproduce the figure using the associated script cmd.sh
#
#--------------------------------------------------------------------------------------------------------

PlotCoveredRegions <- function(dataS,fillcolor,addBoundaryLabel,ringOffset)
{
	totSize <- 0
	circos.rect(0 ,ringOffset-halftrackheight ,template.len ,  ringOffset+halftrackheight + 0.0 , col='grey',border='grey')
	if (nrow(dataS) > 0)
	{
		for (i in 1:nrow(dataS))
		{
				if (dataS$rightFlank[i] > dataS$leftFlank[i])
				{
					circos.rect(dataS$leftFlank[i] ,ringOffset-5*halftrackheight ,dataS$rightFlank[i],  ringOffset+5*halftrackheight , col=fillcolor,border=fillcolor)
				} else {
					circos.rect(dataS$leftFlank[i] ,ringOffset-5*halftrackheight ,template.len,  ringOffset+5*halftrackheight , col=fillcolor,border=fillcolor)
					circos.rect(1                  ,ringOffset-5*halftrackheight ,dataS$rightFlank[i],  ringOffset+5*halftrackheight , col=fillcolor,border=fillcolor)
				}
				totSize <- totSize + dataS$len[i]
		}
	}
	return(cbind(Count=nrow(dataS),TransferSize=totSize))
}
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#
#                                   MAIN
#
#  NOTE:  expects input   in a directory plots_T?/
#                 outputs plots and text summary in current directory
#
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------


options(scipen=13)
library('circlize')

args=commandArgs(TRUE)
whichTemplate <- as.character(args[1])
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


filelist <- list.files(path = paste('./T',whichTemplate,sep=''), pattern = '.trf')
results <- NULL

fn <- paste('./T',whichTemplate,'_summary.pdf',sep='')
pdf(file=fn,width=30,height=15)
layout(matrix(1:2,ncol=2,byrow=T))

# plot page 1

	# -- Chr 1 of Template -----------------------------------------------------

	template.color <- 'white'
	template.bordercolor <- 'red'
	halftrackheight <- 0.01
	axisFontSize <- 1.333

	testsize <- 1.5
	maxY <- 25


	resistance.textsize <- 1.5
	resistance.color <- 'skyblue'

	plotStart <- 1.2  # do not change
	templateLen.fontsize <- 1.1
	htpos.fontsize <- 0.9

	template.color <- 'white'
	template.bordercolor <- 'darkgreen'

	template.len <- template.Chr1.len
	marginFactor <- 2.5
	margins <- c(-marginFactor,marginFactor)
	par(mar=c(2,2,2,2))
	circos.clear()
	circos.par("start.degree" = 90, canvas.xlim = margins, canvas.ylim = margins,
			gap.degree = 0, cell.padding = c(0.0, 0.0, 0.0, 0.0), track.margin = c(0.001, 0.001),track.margin = c(0.001, 0.001), track.height = 0.1)
	circos.initialize(factors = 'a', xlim = c(1,template.len))
	circos.par(track.margin = c(0, convert_height(10, "mm")))
	circos.track(ylim = c(0, 1), track.height = convert_height(20, "mm"),bg.border='white')


	text(0.0,  0.15, "Chr I", adj = c(0.5,0.5) , cex=testsize)
	text(0.0,  0.03, paste("Template",TEMPLATE), adj = c(0.5,0.5),  , cex=testsize)
	text(0.0,  -0.07, "based on CP028892 (Sa5Y)", adj = c(0.5,0.5), , cex=0.4*testsize)

	circos.track(ylim = c(0, 1), track.height = convert_height(5, "mm"),bg.border='white')
	circos.axis(h=maxY,major.at = seq(0,template.len,by=100000),labels.cex=axisFontSize,labels.facing='clockwise' , direction='outside')


	# plot resistance K7

	if (RESISTANCE1.present) {
	  for (i in 1:length(RESISTANCE1.start)) {
		circos.rect(RESISTANCE1.start[i], 0 , RESISTANCE1.end[i], maxY+2 ,col=resistance.color,border=resistance.color)
		circos.text(RESISTANCE1.start[i]  ,maxY+6,   RESISTANCE1.text[i]  ,facing='inside',niceFacing=TRUE,cex=resistance.textsize,adj=c(0.5,0.5),col=resistance.color)
	  }
	}

	# ---- loop to plot each clone coverage
	step <- maxY / length(filelist)
	for (filenum in 1:length(filelist))
	{
		sample <- unlist(strsplit(filelist[filenum],'.trf'))
		f <- (maxY - filenum * step)

		# --- read transferred regions

		retainedTransfers <- read.table(paste('plots_T',whichTemplate,'/',sample,'_transferred.txt',sep=''),header=T,sep="\t")
		retainedTransfers <- retainedTransfers[retainedTransfers$chr == 'chr1',]
		temp <- PlotCoveredRegions(retainedTransfers,'darkgreen',TRUE, f )
		temp <- cbind(Sample=sample,Chr='Chr1',temp)
		results <- rbind(results,temp)
	} # for


	# -- Chr 2 of Template -----------------------------------------------------

	maxY <- 16
	template.len <- template.Chr2.len

	margins <- c(-marginFactor,marginFactor)
	par(mar=c(2,2,2,2))
	circos.clear()
	circos.par("start.degree" = 90, canvas.xlim = margins, canvas.ylim = margins,
			gap.degree = 0, cell.padding = c(0.0, 0.0, 0.0, 0.0), track.margin = c(0.001, 0.001),track.margin = c(0.001, 0.001), track.height = 0.1)
	circos.initialize(factors = 'a', xlim = c(plotStart,template.len))
	circos.par(track.margin = c(0, convert_height(10, "mm")))
	circos.track(ylim = c(0, 1), track.height = convert_height(20, "mm"),bg.border='white')

	text(0.0,  0.15, "Chr II", adj = c(0.5,0.5) , cex=testsize)
	text(0.0,  0.03, paste("Template",TEMPLATE), adj = c(0.5,0.5),  , cex=testsize)
	text(0.0,  -0.07, "based on CP028893 (Sa5Y)", adj = c(0.5,0.5), , cex=0.4*testsize)

	circos.track(ylim = c(0, 1), track.height = convert_height(5, "mm"),bg.border='white')
	circos.axis(h = maxY, major.at = seq(0,template.len,by=50000),labels.cex=axisFontSize,labels.facing='clockwise' , direction='outside')

	# plot resistance K7

	if (RESISTANCE2.present) {
	  for (i in 1:length(RESISTANCE2.start)) {
	    circos.rect(RESISTANCE2.start[i], 0, RESISTANCE2.end[i], maxY+2 ,col=resistance.color,border=resistance.color)
	    circos.text(RESISTANCE2.start[i]  ,maxY+6 ,   RESISTANCE2.text[i]  ,facing='inside',niceFacing=TRUE,cex=resistance.textsize,adj=c(0.5,0.5),col=resistance.color)
	  }
	}


	# ---- loop to plot each clone coverage
	step <- maxY / length(filelist)
	for (filenum in 1:length(filelist))
	{

		sample <- unlist(strsplit(filelist[filenum],'.trf'))
		f <- (maxY - filenum * step)

		# --- read transferred regions

		retainedTransfers <- read.table(paste('plots_T',whichTemplate,'/',sample,'_transferred.txt',sep=''),header=T,sep="\t")
		retainedTransfers <- retainedTransfers[which(retainedTransfers$chr == 'chr2'),]
		temp <- PlotCoveredRegions(retainedTransfers,'darkgreen',TRUE,f)
		temp <- cbind(Sample=sample,Chr='Chr2',temp)
		results <- rbind(results,temp)

	} # for

	# output summary of transfers

	fn <- paste('./T',whichTemplate,'_summary.txt',sep='')
	write.table(results,file=fn,sep="\t",row.names=F,quote=F)

# plot page 2

	layout(matrix(1:2,ncol=2))
	par(mar=c(8,8,8,8))
	x <- read.table(file=fn,header=T)
	x1 <- x[x$Chr == "Chr1",]
	x2 <- x[x$Chr == "Chr2",]
	mx <- max(c(x1$Count,x2$Count))
	data <- rbind( table(factor(x1$Count,levels=0:mx)) , table(factor(x2$Count,levels=0:mx)) )
	barplot(data,beside=T,las=1,legend=c('Chr I','Chr II'),xlab='Transfers',ylab='Frequency')
	boxplot(x1$TransferSize+x2$TransferSize,ylab='Total Transferred Length')

dev.off()


