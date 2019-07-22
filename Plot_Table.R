library(openxlsx)
library(ggplot2)
library(dplyr)
library(nortest)
library(clusterSim)
library(zoo)
library(qwraps2)

ExcelFile <- file.choose()

df1 <- read.xlsx(ExcelFile, sheet = 1, skipEmptyRows = FALSE)

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
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
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

find <- function(df1, p1, p2){
  
  # p1 <- parametro do tipo de sinal
  # p2 <- parametro da característica
  
  # p1 <- "GYR1.X"
  # p2 <- "PeakAmp"
  
 Nevents <- unique(df1$Event)
 index1 <- grep(p1, df1$Signal)
 index2 <- grep(p2, df1$Feature)

 indexf2 <- Reduce(intersect, list(index1,index2))

 idxt.f2 <- lapply(Nevents, function(x){Reduce(intersect, list(indexf2, grep(x, df1$Event)))})
 
 value <- lapply(idxt.f2, function(x){df1$Value[x]})

  return(value)
}

# Multiple boxplots ----------------------------------------------------------------------------------

Mult.boxplot <- function(df1, sensor, feature){ 
  
  Nsensor <- length(sensor)
  s <- sensor
  lit <- vector("list", Nsensor)
  
  for(i in 1:Nsensor){
    
    k <- find(df1, s[i], feature)
    d <- data.frame(PeakAmp = unlist(k), Events = rep(c(1:length(k)),times = sapply(k,length)))
    lit[[i]] <- ggplot(d, aes(x = as.factor(Events), y = PeakAmp)) +  geom_boxplot(fill="skyblue") + xlab("Event") + ylab("Peak Amplitude") +
      ggtitle(paste("Plot of", sep = " ", s[i], collapse = NULL))
    
  }
  
  return(multiplot(plotlist = lit, cols=4))
   
}

s1 <- c("GYR1.X", "GYR1.Y", "GYR1.Z", "ACC1.X", "ACC1.Y", "ACC1.Z", "MAG1.X", "MAG1.Y", "MAG1.Z", "EMG1")
s2 <- c("GYR2.X", "GYR2.Y", "GYR2.Z", "ACC2.X", "ACC2.Y", "ACC2.Z", "MAG2.X", "MAG2.Y", "MAG2.Z", "EMG2")
s3 <- c("GYR3.X", "GYR3.Y", "GYR3.Z", "ACC3.X", "ACC3.Y", "ACC3.Z", "MAG3.X", "MAG3.Y", "MAG3.Z")
s4 <- c("GYR4.X", "GYR4.Y", "GYR4.Z", "ACC4.X", "ACC4.Y", "ACC4.Z", "MAG4.X", "MAG4.Y", "MAG4.Z")

Mult.boxplot(df1, s1, "PeakAmp")
Mult.boxplot(df1, s2, "PeakAmp")
Mult.boxplot(df1, s3, "PeakAmp")
Mult.boxplot(df1, s4, "PeakAmp")

Mult.boxplot(df1, s1, "ValleyAmp")
Mult.boxplot(df1, s2, "ValleyAmp")
Mult.boxplot(df1, s3, "ValleyAmp")
Mult.boxplot(df1, s4, "ValleyAmp")

Mult.boxplot(df1, s1, "PeakTime")
Mult.boxplot(df1, s2, "PeakTime")
Mult.boxplot(df1, s3, "PeakTime")
Mult.boxplot(df1, s4, "PeakTime")

Mult.boxplot(df1, s1, "ValleyTime")
Mult.boxplot(df1, s2, "ValleyTime")
Mult.boxplot(df1, s3, "ValleyTime")
Mult.boxplot(df1, s4, "ValleyTime")

Mult.boxplot(df1, s1, "PPI")
Mult.boxplot(df1, s2, "PPI")
Mult.boxplot(df1, s3, "PPI")
Mult.boxplot(df1, s4, "PPI")

Mult.boxplot(df1, s1, "VVI")
Mult.boxplot(df1, s2, "VVI")
Mult.boxplot(df1, s3, "VVI")
Mult.boxplot(df1, s4, "VVI")

Mult.boxplot(df1, s1, "PVI")
Mult.boxplot(df1, s2, "PVI")
Mult.boxplot(df1, s3, "PVI")
Mult.boxplot(df1, s4, "PVI")


# Multiple histograms ----------------------------------------------------------------------------------

Mult.histograms <- function(df1, sensor, feature){ 
  
  Nsensor <- length(sensor)
  s <- sensor
  lit <- vector("list", Nsensor)
  
  par(mfrow=c(3, 3))
 
  for(i in 1:Nsensor){
    
    k <- find(df1, s[i], feature)
    d <- data.frame(PeakAmp = unlist(k), Events = rep(c(1:length(k)),times = sapply(k,length)))
    lit[[i]] <- hist(d$PeakAmp, main=paste("Plot of", sep = " ", s[i], collapse = NULL))
    
  }
  
  return(lit)
  
}

Mult.histograms(df1, s1, "PeakAmp")

# Multiple Density -------------------------------------------------------------------------------------

Mult.density <- function(df1, sensor, feature){ 
  
  Nsensor <- length(sensor)
  s <- sensor
  lit <- vector("list", Nsensor)
  
  par(mfrow=c(3, 3))
  
  for(i in 1:Nsensor){
    
    k <- find(df1, s[i], feature)
    d <- data.frame(PeakAmp = unlist(k), Events = rep(c(1:length(k)),times = sapply(k,length)))
    dens <- density(d$PeakAmp)
    plot(dens, type="n", main=paste("Plot of", sep = " ", s[i], collapse = NULL))
    lit[[i]] <- polygon(dens, col="skyblue", border="gray")
    
  }
  
  return(lit)
  
}

Mult.density(df1, s1, "PeakAmp")

# Mean and standard deviation of interval time-----------------------------------------------

Mean.sd.Time <- function(df1, sensor){ 
  
  Nsensor <- length(sensor)
  s <- sensor
  
  feature <- c("PPI", "VVI", "PVI")
  LPPI <- vector("list", Nsensor)
  LVVI <- vector("list", Nsensor)
  LPVI <- vector("list", Nsensor)
  
  for(i in 1:Nsensor){ # Select the intervals times
    
    LPPI[[i]] <- unlist(find(df1, s[i], feature[1]))
    LVVI[[i]] <- unlist(find(df1, s[i], feature[2]))
    LPVI[[i]] <- unlist(find(df1, s[i], feature[3]))
    
  }
  
  l1 <- vector("list", 1)
  l2 <- vector("list", 1)
  l3 <- vector("list", 1)
  
  for(i in 1:Nsensor){ # Calculating the mean and standard deviation of sensors
    
      l1[[1]][[i]] <- paste(round(mean(LPPI[[i]]), 3), round(sd(LPPI[[i]]), 3), sep=" ± ") 
      l2[[1]][[i]] <- paste(round(mean(LVVI[[i]]), 3), round(sd(LVVI[[i]]), 3), sep=" ± ") 
      l3[[1]][[i]] <- paste(round(mean(LPVI[[i]]), 3), round(sd(LPVI[[i]]), 3), sep=" ± ")
  }
  
  # Create the excel file data:
  PPI <- data.frame(Sensor = s)
  VVI <- data.frame(Sensor = s)
  PVI <- data.frame(Sensor = s)
  
  PPI$PPI.Mean.SD <- unlist(l1[[1]])
  VVI$VVI.Mean.SD <- unlist(l2[[1]])
  PVI$PVI.Mean.SD <- unlist(l3[[1]])
  
  Int <-as.matrix(merge(PPI, merge(VVI, PVI)))
  
  return(Int)
}

Table.MSD1 <- Mean.sd.Time(df1, s1)
Table.MSD2 <- Mean.sd.Time(df1, s2)
Table.MSD3 <- Mean.sd.Time(df1, s3)
Table.MSD4 <- Mean.sd.Time(df1, s4)






