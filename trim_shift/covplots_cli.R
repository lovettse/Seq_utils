require(ggplot2)
require(scales)
setwd(".")

args <- commandArgs(trailingOnly = TRUE)

input_list <- read.table(args[1], header=FALSE, stringsAsFactors=FALSE)
binSize <- as.numeric(args[2])
cols <- as.numeric(args[3])
chrom <- args[4]

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <-  plotlist
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

slidingwindowplot <- function(windowsize, inputseq)
{
  starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
  n <- length(starts)
  chunkbps <- numeric(n)
  chunkstats<- numeric(n)
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)] 
    chunkmean <- mean(chunk)
    chunkstdv<-sd(chunk)
    chunkbps[i] <- chunkmean
    chunkstats[i]<-chunkstdv
  }
  return (list(starts,chunkbps,chunkstats))
}

pl <- list()
out <- matrix()


for (i in 1:dim(input_list)[1]) {
  inputFile = input_list[i,]
  print(inputFile)
  column.types <- c("character", "numeric", "numeric")
  all.data <- read.csv(inputFile, header=FALSE, sep="\t",colClasses=column.types)
  temp <- all.data
  myvector_all<-as.vector(as.matrix(all.data[3]))
  windowAll<-slidingwindowplot(binSize,myvector_all)
  df<-data.frame(windowAll[[1]],windowAll[[2]],windowAll[[3]])
  colname<-c("x","mean","sd")
  colnames(df)<-colname
  maxy <- max(all.data[3])
  
  eb <- aes(ymax=mean+sd,ymin=mean-sd)
  
  pl[[length(pl)+1]] <- ggplot(data = df, aes(x = x, y = mean)) + geom_line(colour="#0066CC",size=0.5) + geom_ribbon(eb, alpha=.5, colour=NA,fill="#0066CC")+theme_bw()+xlab(unlist(strsplit(inputFile,"_"))[1])+ylab("")+scale_x_continuous(expand = c(0,0), labels=comma)+scale_y_continuous(limits = c(0, maxy))+theme(axis.text.x=element_text(angle=90, hjust=0.5))
}

pdf(file=paste(chrom,"covplots.pdf",sep="_"))
multiplot(plotlist = pl,file=paste(chrom,"covplots.pdf",sep="_"),cols=cols)
dev.off()
