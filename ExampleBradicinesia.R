
# Libraries ---------------------------------------------------------------
source("Thresholding.r")
source("TremsenToolbox.r")
library(openxlsx)
library(stringi)
library(grid)
library(gridExtra)
library(dplyr)


# Load TREMSEN file -------------------------------------------------------

# Choose a TREMSEN (.txt) file 
# strFileName <- file.choose() # If data are available, select the file

strFileName <- file.choose() # If data are available, select the file

df <- LoadTREMSENFile(strFileName) # Load TREMSEN file

## Plotting the signal
dfplot1<-data.frame(df$X.Time., df$X.A1.X.)
dygraph(dfplot1)%>%dyRangeSelector()

# Entrada dos parâmetros do indivíduo:
Directory=strFileName
Person="AFP"
Gender="M"
Age="x"
Condition="Parkinson"
Member="Upper"
Assignment=2  #Extensão, flexão, pronação, supinação e pinça

Ident <- data.frame(Directory, Person, Gender, Age, Condition, Member, Assignment)

# Detrend and remove DC offset --------------------------------------------
df.nonlineardetrended <- nonLineardetrendTremsenData(df) #Remove both linear and nonlinear trends

# Event detection  ---------------------------------------------------

# Using detrend dataframe for event detection
dftmp<-df.nonlineardetrended
eventDetect = detectEvent(dftmp)

channel <- as.list(df.nonlineardetrended) # Converte o dataframe df.nonlineardetrend em uma lista



# Create labels and add to PULSE.LABEL ----------------------------------------------------

Labels <- createLabels(eventDetect)
df.nonlineardetrended$X.PULSE.LABEL <- Labels
channel <- as.list(df.nonlineardetrended)


# Filter signal -----------------------------------------------------------
# Define butterworth filter to remove high frequency noise
Fs <- 1/0.02
Fn <- Fs/2
f <- 2
per <-f/Fn
bf <- butter(3, per, type="low") # Filter out high frequency noise

# Filtrando o canais dos dados(menos a primeira e as últimas quatro colunas do dataframe)
ff <- rapply(channel[-c(1,40,41,42,43)], function(x){filtfilt(bf, x)}, how = "list", classes = "numeric")
ff <- as.data.frame(ff)
ff <- cbind(Time=channel[[1]], ff)

#write.xlsx(ff, file = "DataFilt.2.xlsx", col.names = TRUE)


# Partition signal ----------------------------------------------- 
sep <- findSep(eventDetect)

PartitionSignal <- function(df, ff, sep){
  
  # Separando os dataframes ff e df em canais filtrados e originais:
  Filt.channel1 <- ff[grep(pattern=".*1",colnames(ff))] # Seleciona todos sensores do canal 1 - Simbolo .*1 seleciona todas colunas do df que possuem o número 1
  Filt.channel2 <- ff[grep(pattern=".*2",colnames(ff))]
  Filt.channel3 <- ff[grep(pattern=".*3",colnames(ff))]
  Filt.channel4 <- ff[grep(pattern=".*4",colnames(ff))]

  
  Orig.channel1 <- df[grep(pattern=".*1",colnames(df))]
  Orig.channel2 <- df[grep(pattern=".*2",colnames(df))]
  Orig.channel3 <- df[grep(pattern=".*3",colnames(df))]
  Orig.channel4 <- df[grep(pattern=".*4",colnames(df))]

  
  PartSig <- function(Filt.sensor, Orig.sensor){
    
    partChannel <- rapply(Filt.sensor, function(x){partitionChannel(x, sep)}, how = "list", classes = "numeric")
    partOrignal <- rapply(Orig.sensor, function(x){partitionChannel(x, sep)}, how = "list", classes = "numeric")
    
    part <- list(FiltSignal=partChannel, OrigSignal=partOrignal)
    
    return(part)
  }
  
# Canais separados em eventos:
# Partitions: lista para salvar os canais analisados
# Na primeira lista da partition terá duas outras listas: uma será a separação em eventos do sinal filtrado e a outra do sinal original
# Na partChannel terá eventos de um vetor do sinal filtrado, por exemplo o G1.X
  
 
  partitions <- vector("list", 4)
  partitions[[1]] <- PartSig(Filt.channel1, Orig.channel1)
  partitions[[2]] <- PartSig(Filt.channel2, Orig.channel2)
  partitions[[3]] <- PartSig(Filt.channel3, Orig.channel3)
  partitions[[4]] <- PartSig(Filt.channel4, Orig.channel4)
  
  return(partitions)
  
}

PartitionTime <- function(ff, sep){
  
  Vector.Time <- ff$Time
  partTime <- partitionChannel(Vector.Time, sep)
  
  
  # Encontrando tempo incial e final de cada evento no sinal anlisado:
  t.begin <- rep(NA, length(partTime))
  t.end <- rep(NA, length(partTime))
  
  for(i in 1:length(partTime)){
    for(j in 1:length(partTime[[i]])){
      
      if(j == 1){
        t.begin[i] <- partTime[[i]][[j]]
      }
      if(j == length(partTime[[i]])){
        t.end[i] <- partTime[[i]][[j]]
      }
      
    }
  }
  
  event.tempo <- vector("list", 2)
  event.tempo[[1]] <- t.begin
  event.tempo[[2]] <- t.end
  
  return(event.tempo)
}

Events <- PartitionSignal(df, ff, sep)
Event.Time <- PartitionTime(ff, sep)

# Finding peaks and valleys -----------------------------------------------

# Função para listar os picos e vales de um canal
Find.Peaks.Valley <- function(Events1, Events2){
  
  
  # X <-  Events[[1]][[1]]
  # XOriginal <- Events[[1]][[2]]
  
  X <- Events1
  XOriginal <- Events2
  
  deltaT = 0.02 #Time between data points
  tsig<-c(0:(length(X)-1))*deltaT
  
  
  xMaxIdx<- vector("list", length(X))
  xMaxValue <- vector("list", length(X))
  xMinIdx<- vector("list", length(X))
  xMinValue <- vector("list", length(X))
  xMaxtime <- vector("list", length(X))
  xMintime <- vector("list", length(X))
  
  
  PeakValley <- rapply(X, findPeakValley, how = "list", classes = "numeric")
  # Função findPeakValley retorna uma lista com xMaxIdx, xMaxValue, xMinIdx e xMinValue respectivamente
  
# PeakValley: lista para salvar picos e vales do canais analisados
# Dentro dessa lista, existirá listas para cada um dos sensores ou sinal EMG estudado. 
# Nesse exemplo, serão os sensores: G1.X, G1.Y, G1.Z, A1.X, A1.Y, A1.Z, M1.X, M1.Y, M1.Z
# Para cada um desses sensores, existirá listas com os eventos detectados com a função PartitionSignal
# Dentro de cada evento existirá os picos e vales dectectados com a função findPeakValley
  
  for(i in 1:length(X)){
    for(j in 1:length(X[[1]])){
      
      xMaxIdx[[i]][[j]] <- PeakValley[[i]][[j]][[1]]
      xMaxValue[[i]][[j]] <- PeakValley[[i]][[j]][[2]]
      xMinIdx[[i]][[j]] <-PeakValley[[i]][[j]][[3]]
      xMinValue[[i]][[j]] <- PeakValley[[i]][[j]][[4]]
      xMaxtime[[i]][[j]] <- ff$Time[PeakValley[[i]][[j]][[1]]]
      xMintime[[i]][[j]] <- ff$Time[PeakValley[[i]][[j]][[3]]]
      
    }
  }
  
  
  
  Loc <- vector("list", 6)
  Loc[[1]] <- xMaxIdx # Indíce do picos
  Loc[[2]] <- xMinIdx # Indíce vales
  Loc[[3]] <- xMaxValue # Amplitude picos
  Loc[[4]] <- xMinValue # Amplitude Vales
  Loc[[5]] <- xMaxtime # Tempo de ocorrência de picos
  Loc[[6]] <- xMintime # Tempo de ocorrência de vales
    
  return(Loc)
}   

# Picos e vales dos canais 1 à 4 e dos sinais EMG1 e EMG2
# Primeiro parâmetro é o sinal filtrado e o segundo o sinal bruto
PeakValley_c1 <- Find.Peaks.Valley(Events[[1]][[1]], Events[[1]][[2]])
PeakValley_c2 <- Find.Peaks.Valley(Events[[2]][[1]], Events[[2]][[2]])
PeakValley_c3 <- Find.Peaks.Valley(Events[[3]][[1]], Events[[3]][[2]])
PeakValley_c4 <- Find.Peaks.Valley(Events[[4]][[1]], Events[[4]][[2]])


# Excel files ----------------------------------------------------------

Inter <- function(xMaxIdx,xMinIdx){
  
  Intervals <- vector("list", length(xMaxIdx[[1]]))
  
  for(i in 1:length(xMaxIdx)){
    Intervals[[i]] <- vector("list", length(xMaxIdx[[1]]))
  }
  
  for(i in 1:length(xMaxIdx)){
    for(j in 1:length(xMaxIdx[[i]])){
      
      Intervals[[i]][[j]] <- findInterval(xMaxIdx[[i]][[j]], xMinIdx[[i]][[j]])
      
    }
  }
  
  return(Intervals) 
  
}

# Intervalos de tempos entre picos, vales e entre picos/vales de cada um dos canais e sinais EMG:
I1 <- Inter(PeakValley_c1[[1]], PeakValley_c1[[2]])
I2 <- Inter(PeakValley_c2[[1]], PeakValley_c2[[2]])
I3 <- Inter(PeakValley_c3[[1]], PeakValley_c3[[2]])
I4 <- Inter(PeakValley_c4[[1]], PeakValley_c4[[2]])


Interval.Peak.Valley <- function(Intervals){
  
  # Intervals <- I1
  
  xInterPeakInterval <- vector("list", length(Intervals))
  xInterValleyInterval <- vector("list", length(Intervals))
  xCrossInterval <- vector("list", length(Intervals))
  
  for(i in 1:length(Intervals)){
    for(j in 1:length(Intervals[[1]])){
      if(is.null(Intervals[[i]])){
        break
      }
      else{
        xInterPeakInterval[[i]][[j]] <- Intervals[[i]][[j]][[1]]
        xInterValleyInterval[[i]][[j]] <- Intervals[[i]][[j]][[2]]
        xCrossInterval[[i]][[j]] <- Intervals[[i]][[j]][[3]]
      }
    }
  }
  
  Inter.PV <- vector("list", 3)
  Inter.PV[[1]] <- xInterPeakInterval
  Inter.PV[[2]] <- xInterValleyInterval
  Inter.PV[[3]] <- xCrossInterval
  
   return(Inter.PV)

}

# Listas divididas em tempos entre picos[[1]], vales[[2]] e picos/vales[[3]] nos canais de 1 à 4 e nos sinais EMG:
InterPV.1 <- Interval.Peak.Valley(I1) # CANAL 1
InterPV.2 <- Interval.Peak.Valley(I2) # CANAL 2
InterPV.3 <- Interval.Peak.Valley(I3) # CANAL 3
InterPV.4 <- Interval.Peak.Valley(I4) # CANAL 4


AmpTime <- function(PeakValley, InterPV, Ident, channel, Event.Time){

   AT <- function(PeakValley, InterPV, Ident, sensor, channel, Event.Time){ 

     Name.signal <- case_when(
       sensor == 1 ~ paste("GYR", sep = "", channel, ".X", collapse = NULL),
       sensor == 2 ~ paste("GYR", sep = "", channel, ".Y", collapse = NULL),
       sensor == 3 ~ paste("GYR", sep = "", channel, ".Z", collapse = NULL),
       sensor == 4 ~ paste("ACC", sep = "", channel, ".X", collapse = NULL),
       sensor == 5 ~ paste("ACC", sep = "", channel, ".Y", collapse = NULL),
       sensor == 6 ~ paste("ACC", sep = "", channel, ".Z", collapse = NULL),
       sensor == 7 ~ paste("MAG", sep = "", channel, ".X", collapse = NULL),
       sensor == 8 ~ paste("MAG", sep = "", channel, ".Y", collapse = NULL),
       sensor == 9 ~ paste("MAG", sep = "", channel, ".Z", collapse = NULL),
       sensor == 10 ~ paste("EMG", sep = "", channel, collapse = NULL),
       TRUE ~ as.character(sensor)
     )
     
  # Valores tempo inicio e fim
  # Criando dataframes para armazenar os tempos de picos, vales e picos/vales:
     t <- length(unlist(Event.Time[[1]]))
     df.begin <- data.frame(Count=seq(from=1, to=t, by=1))
     df.end <- data.frame(Count=seq(from=1, to=t, by=1))
     
     beginEvent <- seq(from=1, to=length(Event.Time[[1]]), by=1)
     endEvent <- seq(from=1, to=length(Event.Time[[2]]), by=1)
     
     df.begin$Signal <- Name.signal
     df.end$Signal <- Name.signal
     
     df.begin$Event <- beginEvent
     df.end$Event <- endEvent
     
     df.begin$Value <-  unlist(Event.Time[[1]])
     df.end$Value <- unlist(Event.Time[[2]])
     
     df.begin$Feature<-'beginTime'  #Feature
     df.end$Feature<-'endTime'
  
     evt.Time <-rbind(df.begin, df.end)
     evt.Time$Count <- NULL

  datSave <- function(v1, v2, c1, c2, sensor){
    
    n1 <- length(unlist(v1[[sensor]]))
    n2 <- length(unlist(v2[[sensor]]))
    
    df1 <- data.frame(Count=seq(from=1, to=n1, by=1))
    df2 <- data.frame(Count=seq(from=1, to=n2, by=1))
    
    v1Event <- vector("list", length(v1[[sensor]]))
    v2Event <- vector("list", length(v2[[sensor]]))
    
    for(i in 1:length(v1[[sensor]])){
      v1Event[[i]] <- rep(i, length(v1[[sensor]][[i]]))
      v2Event[[i]] <- rep(i, length(v2[[sensor]][[i]]))
    }
    
    df1$Signal <- Name.signal
    df2$Signal <- Name.signal
    
    df1$Event <- as.numeric(as.character(unlist(v1Event)))
    df2$Event <- as.numeric(as.character(unlist(v2Event)))
    
    df1$Value <-  as.numeric(as.character(unlist(v1[[sensor]])))
    df2$Value <-  as.numeric(as.character(unlist(v2[[sensor]])))
    
    df1$Feature<-c1  #Feature
    df2$Feature<-c2
    
    dft <-rbind(df1, df2)
    dft$Count <- NULL
    
    return(dft)
    
  }
  
  AmptPeakValley <- datSave(PeakValley[[3]], PeakValley[[4]], 'PeakAmp', 'ValleyAmp', sensor) # Amplitudes picos e vales
  
  TPV <- datSave(PeakValley[[5]], PeakValley[[6]], 'PeakTime', 'ValleyTime', sensor) # Tempo dos picos e vales
  
  
  # Valores Intervalos de Tempo
    # Criando dataframes para armazenar os tempos de picos, vales e picos/vales:
    npt <- length(unlist(InterPV[[1]][[sensor]]))
    nvt <- length(unlist(InterPV[[2]][[sensor]]))
    nct <- length(unlist(InterPV[[3]][[sensor]]))
    
    dfPeak <- data.frame(Count=seq(from=1, to=npt, by=1))
    dfValley <- data.frame(Count=seq(from=1, to=nvt, by=1)) 
    dfCross <- data.frame(Count=seq(from=1, to=nct, by=1)) 
    
    # Adicionando colunas eventos e distância aos dataframes:
    
    PeakEvent <- vector("list", length(InterPV[[1]]))
    ValleyEvent <- vector("list", length(InterPV[[2]]))
    CrossEvent <- vector("list", length(InterPV[[3]]))
    
    for(i in 1:length(InterPV[[1]][[1]])){
      # [[1]][[sensor]][[1]] -> Definição intervalo de tempo - sensor - evento
      PeakEvent[[i]] <- rep(i, length(InterPV[[1]][[sensor]][[i]]))
      ValleyEvent[[i]] <- rep(i, length(InterPV[[2]][[sensor]][[i]]))
      CrossEvent[[i]] <- rep(i, length(InterPV[[3]][[sensor]][[i]]))
    }
    
    # Definindo quais trechos ocorreu os eventos:
    dfPeak$Signal <- Name.signal
    dfValley$Signal <- Name.signal
    dfCross$Signal <- Name.signal
    
    dfPeak$Event <- unlist(PeakEvent)
    dfValley$Event <- unlist(ValleyEvent)
    dfCross$Event <- unlist(CrossEvent)
    
    dfPeak$Value <-  unlist(InterPV[[1]][[sensor]])
    dfValley$Value <- unlist(InterPV[[2]][[sensor]])
    dfCross$Value <- unlist(InterPV[[3]][[sensor]])
    
    #dfPeak$dist <- rep(NA,length(xInterPeakInterval))
    dfPeak$Feature<-'PPI'
    dfValley$Feature<-'VVI'
    dfCross$Feature<-'PVI'
    
    IntTime <-rbind(dfPeak, dfValley, dfCross)
    IntTime$Count <- NULL
    d1 <- as.matrix(evt.Time)
    d2 <- as.matrix(AmptPeakValley)
    d3 <- as.matrix(TPV)
    d4 <- as.matrix(IntTime)
    d5 <- as.data.frame(rbind(d1, d2, d3, d4))
    d5 <- cbind(Ident, d5)
    
    return(d5)
   }
   
   a <- length(PeakValley[[1]])
   
   at <- vector("list", a)
   
   for(i in 1:a){
     at[[i]] <- AT(PeakValley, InterPV, Ident, i, channel, Event.Time) 
   }
   
   at<-do.call(rbind, at)
  
   return(at)
  
}

ch1 <- as.matrix(AmpTime(PeakValley_c1, InterPV.1, Ident, 1, Event.Time))
ch2 <- as.matrix(AmpTime(PeakValley_c2, InterPV.2, Ident, 2, Event.Time))
ch3 <- as.matrix(AmpTime(PeakValley_c3, InterPV.3, Ident, 3, Event.Time))
ch4 <- as.matrix(AmpTime(PeakValley_c4, InterPV.4, Ident, 4, Event.Time))
f5 <- as.data.frame(rbind(ch1, ch2, ch3, ch4))
f5 <- type.convert(f5)  # Corrige problemas de escrever numeros como caracteres

write.xlsx(f5, file="Amp.Time.AFP2.xlsx")



