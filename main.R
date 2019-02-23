#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Master's Thesis - Remote Sensing                      %
# Environmental Engineering - ISA/UL - Lisbon, Portugal %
# (c) 2014 by Jonas Schmedtmann & Manuel Campagnolo     %
#                                                       %
# MAIN SCRIPT                                           %
#                                                       %
# Implements the program logic using functions from     %
# functions.R. All steps typically found in remote      %
# sensing can be found here.                            %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Processing steps:
#   01. DATA AQUISITION AND PROCESSING
#   02. PARCEL DATA SELECTION
#   03. EXPLORATORY DATA ANALYSIS
#   04. VARIABLE SELECTION
#   05. 06. 07. SELECTION/TRAINING/VALIDATION OF CLASSIFIERS


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 1. DATA AQUISITION AND PROCESSING                                                   ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Loading all fucntions and initializing program
source('functions.R')
init()

# Defining study area, clipping, and projecting
coordsArea <- cbind(AREA.X,AREA.Y)
cdg <- carregaRecortaDadosEspaciais(coordsArea, PROJ4.UTM)
plot(cdg$area);plot(cdg$parc2005, add=T)

# Loading and correcting images
rm(todasImagens)
todasImagens <- constroiListaImagensLandsat(landsatPath=CAMINHO.LANDSAT,
                                            areaEstudo=cdg$area,
                                            prefixo="CORR_14.08",
                                            ano=2005,
                                            corrige=TRUE)

# Building a list holding all data
rm(listaDados)
listaDados <- constroiListaDados(ano=2005)

# Getting a list of all parcels (data.frame)
listaTodasParcelasIniciais <- constroiTodasParcelas()


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2. PARCEL DATA SELECTION                                                            ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Excluding some crop codes that don't make sense to include (e.g. too difficult to classify)
codExclusao  <-  c(87,88,666)
listaTodasParcelas <- listaTodasParcelasIniciais[!(listaTodasParcelasIniciais$cultura %in% codExclusao),]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Determining which crops, together, occupy the majority of study area (tested for 90%, 95%, 98% of area). This is limit the amount of crops to anazlize

# Getting total parcel area of all crops
dados <- data.table(listaTodasParcelas)
areasCulturas <- dados[,list(area=sum(area),numParc=length(area)),by=cultura]
areasCulturas <- areasCulturas[order(areasCulturas$area, decreasing = TRUE),]
areasCulturas <- cbind(areasCulturas,cumsum(areasCulturas$area))
areaTotal <- sum(areasCulturas$area)

# Visualizing crops for the 3 area thresholds
plot(1:nrow(areasCulturas),areasCulturas$V2)
abline(h=areaTotal*0.9,col='blue')
abline(h=areaTotal*0.95,col='orange')
abline(h=areaTotal*0.98,col='red')

# Selecting 95% of area as a good candidate
limite<-0.95
cultInfl <- areasCulturas[areasCulturas$V2 < areaTotal*limite,]
cultInfl <- cultInfl[!cultInfl$cultura == 27]
fraccaoInfl <- sum(cultInfl$area)/areaTotal
cultInfluentes <- cultInfl$cultura
length(cultInfluentes)

# Number or remaining crops
nrow(areasCulturas)-13

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Selecting the most influential crops (95% of total area) and reclassifying parcels to a standard 1-12 crop code to make analysis more sytraightforward
listaTodasParcelas <- listaTodasParcelas[listaTodasParcelas$cultura %in% cultInfluentes,]

novasClasses <- as.data.frame(cbind(cultInfluentes, c(1,2,3,4,5,6,7,8,5,9,10,11,12)))
colnames(novasClasses) <- c('cultura','novaClasse')
nClasses <- length(table(novasClasses$novaClasse))

for(i in 1:length(listaTodasParcelas$cultura))
  listaTodasParcelas$cultura[i] <- novasClasses$novaClasse[which(novasClasses$cultura==listaTodasParcelas$cultura[i])]


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 3. EXPLORATORY DATA ANALYSIS                                                        ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##### ANOVA to determine if crop and parcel affect reflectance values


#%%%%%%%%%%%%%%%%%%%%%%%%%%
# APPROACH 1)
dadosANOVA <- constroiDadosANOVA(data=4, banda=4, dimAmostra=11852/2)
table(dadosANOVA$cultura)
plot(dadosANOVA$reflectancias, dadosANOVA$cultura)

nrow(dadosANOVA)
length(table(dadosANOVA$parcela))

load('finalListaDados2005_03.12.Robj')
load('dadosANOVA_18.04.2015.Robj')



dadosANOVA.d1.b1 <- constroiDadosANOVA(data=1, banda=1, dimAmostra=11852)
dadosANOVA.d1.b2 <- constroiDadosANOVA(data=1, banda=2, dimAmostra=11852)
dadosANOVA.d1.b3 <- constroiDadosANOVA(data=1, banda=3, dimAmostra=11852)
dadosANOVA.d1.b4 <- constroiDadosANOVA(data=1, banda=4, dimAmostra=11852)
dadosANOVA.d1.b5 <- constroiDadosANOVA(data=1, banda=5, dimAmostra=11852)
dadosANOVA.d1.b6 <- constroiDadosANOVA(data=1, banda=6, dimAmostra=11852)

dadosANOVA.d2.b1 <- constroiDadosANOVA(data=2, banda=1, dimAmostra=11852)
dadosANOVA.d2.b2 <- constroiDadosANOVA(data=2, banda=2, dimAmostra=11852)
dadosANOVA.d2.b3 <- constroiDadosANOVA(data=2, banda=3, dimAmostra=11852)
dadosANOVA.d2.b4 <- constroiDadosANOVA(data=2, banda=4, dimAmostra=11852)
dadosANOVA.d2.b5 <- constroiDadosANOVA(data=2, banda=5, dimAmostra=11852)
dadosANOVA.d2.b6 <- constroiDadosANOVA(data=2, banda=6, dimAmostra=11852)

dadosANOVA.d3.b1 <- constroiDadosANOVA(data=3, banda=1, dimAmostra=11852)
dadosANOVA.d3.b2 <- constroiDadosANOVA(data=3, banda=2, dimAmostra=11852)
dadosANOVA.d3.b3 <- constroiDadosANOVA(data=3, banda=3, dimAmostra=11852)
dadosANOVA.d3.b4 <- constroiDadosANOVA(data=3, banda=4, dimAmostra=11852)
dadosANOVA.d3.b5 <- constroiDadosANOVA(data=3, banda=5, dimAmostra=11852)
dadosANOVA.d3.b6 <- constroiDadosANOVA(data=3, banda=6, dimAmostra=11852)

dadosANOVA.d4.b1 <- constroiDadosANOVA(data=4, banda=1, dimAmostra=11852)
dadosANOVA.d4.b2 <- constroiDadosANOVA(data=4, banda=2, dimAmostra=11852)
dadosANOVA.d4.b3 <- constroiDadosANOVA(data=4, banda=3, dimAmostra=11852)
dadosANOVA.d4.b4 <- constroiDadosANOVA(data=4, banda=4, dimAmostra=11852)
dadosANOVA.d4.b5 <- constroiDadosANOVA(data=4, banda=5, dimAmostra=11852)
dadosANOVA.d4.b6 <- constroiDadosANOVA(data=4, banda=6, dimAmostra=11852)

dadosANOVA.d5.b1 <- constroiDadosANOVA(data=5, banda=1, dimAmostra=11852)
dadosANOVA.d5.b2 <- constroiDadosANOVA(data=5, banda=2, dimAmostra=11852)
dadosANOVA.d5.b3 <- constroiDadosANOVA(data=5, banda=3, dimAmostra=11852)
dadosANOVA.d5.b4 <- constroiDadosANOVA(data=5, banda=4, dimAmostra=11852)
dadosANOVA.d5.b5 <- constroiDadosANOVA(data=5, banda=5, dimAmostra=11852)
dadosANOVA.d5.b6 <- constroiDadosANOVA(data=5, banda=6, dimAmostra=11852)

dadosANOVA.d6.b1 <- constroiDadosANOVA(data=6, banda=1, dimAmostra=11852)
dadosANOVA.d6.b2 <- constroiDadosANOVA(data=6, banda=2, dimAmostra=11852)
dadosANOVA.d6.b3 <- constroiDadosANOVA(data=6, banda=3, dimAmostra=11852)
dadosANOVA.d6.b4 <- constroiDadosANOVA(data=6, banda=4, dimAmostra=11852)
dadosANOVA.d6.b5 <- constroiDadosANOVA(data=6, banda=5, dimAmostra=11852)
dadosANOVA.d6.b6 <- constroiDadosANOVA(data=6, banda=6, dimAmostra=11852)

#%%%%%%%%%%%%%%%%%%%%%%%%%%
# APPROACH 2) Hierarchical ANOVA

# df$x é o factor dominante: cultura
# df$y é o factor subordinado: parcela
# df$z é a resposta (reflectância)

decomposition.hierarq.anova<-function(df, data, banda)
{
  names(df) <- c("z","x","y")
  
  # clean NAs
  df<-df[!is.na(df$z) & !is.na(df$x) & !is.na(df$y),]
  sqa<-0
  sqre<-0
  crops.nparcels<-c()
  great.mean<-mean(df$z)
  for (levC in levels(df$x)) #crops
  {
    pixels.diff<-c()
    crop.mean<-mean(df$z[df$x==levC])
    crop.length<-sum(df$x==levC);#print(crop.length)
    crop.nparcels<-length(unique(as.character(df$y[df$x==levC])))
    crops.nparcels<-c(crops.nparcels,crop.nparcels)
    sqa<-sqa+crop.length*(crop.mean-great.mean)^2 #
    for (levP in unique(as.character(df$y[df$x==levC]))) # parcels
    {
      pixels.diff<-df$z[df$x==levC & df$y==levP]-mean(df$z[df$x==levC & df$y==levP])
      sqre<-sqre+sum(pixels.diff^2)
    }
  }
  
  #print(crops.nparcels) #b_i's
  N<-length(df$z)
  gla<-(nlevels(df$x)-1) # a-1
  glb<-(sum(crops.nparcels-1)) #\sum (b_i -1)
  gle<-(N-sum(crops.nparcels)) #n-\sum b_i
  sqt<-var(df$z)*(N-1)
  sqb<-sqt-sqa-sqre
  qma<-sqa/gla
  qmb<-sqb/glb
  qmre<-sqre/gle
  print(paste(gla,round(sqa,5),round(qma,5),round(qma/qmre,3)))
  print(paste(glb,round(sqb,5),round(qmb,5),round(qmb/qmre,3)))
  print(paste(gle,round(sqre,5),round(qmre,5)))
  return(list(culturas=c(gla,sqa,qma),parcelas=c(glb,sqb,qmb),residuos=c(gle,sqre,qmre), data=data, banda=banda))
}

decomp.d1.b1 <- decomposition.hierarq.anova(dadosANOVA.d1.b1, 1, 1)
decomp.d1.b2 <- decomposition.hierarq.anova(dadosANOVA.d1.b2, 1, 2)
decomp.d1.b3 <- decomposition.hierarq.anova(dadosANOVA.d1.b3, 1, 3)
decomp.d1.b4 <- decomposition.hierarq.anova(dadosANOVA.d1.b4, 1, 4)
decomp.d1.b5 <- decomposition.hierarq.anova(dadosANOVA.d1.b5, 1, 5)
decomp.d1.b6 <- decomposition.hierarq.anova(dadosANOVA.d1.b6, 1, 6)

decomp.d2.b1 <- decomposition.hierarq.anova(dadosANOVA.d2.b1, 2, 1)
decomp.d2.b2 <- decomposition.hierarq.anova(dadosANOVA.d2.b2, 2, 2)
decomp.d2.b3 <- decomposition.hierarq.anova(dadosANOVA.d2.b3, 2, 3)
decomp.d2.b4 <- decomposition.hierarq.anova(dadosANOVA.d2.b4, 2, 4)
decomp.d2.b5 <- decomposition.hierarq.anova(dadosANOVA.d2.b5, 2, 5)
decomp.d2.b6 <- decomposition.hierarq.anova(dadosANOVA.d2.b6, 2, 6)

decomp.d3.b1 <- decomposition.hierarq.anova(dadosANOVA.d3.b1, 3, 1)
decomp.d3.b2 <- decomposition.hierarq.anova(dadosANOVA.d3.b2, 3, 2)
decomp.d3.b3 <- decomposition.hierarq.anova(dadosANOVA.d3.b3, 3, 3)
decomp.d3.b4 <- decomposition.hierarq.anova(dadosANOVA.d3.b4, 3, 4)
decomp.d3.b5 <- decomposition.hierarq.anova(dadosANOVA.d3.b5, 3, 5)
decomp.d3.b6 <- decomposition.hierarq.anova(dadosANOVA.d3.b6, 3, 6)

decomp.d4.b1 <- decomposition.hierarq.anova(dadosANOVA.d4.b1, 4, 1)
decomp.d4.b2 <- decomposition.hierarq.anova(dadosANOVA.d4.b2, 4, 2)
decomp.d4.b3 <- decomposition.hierarq.anova(dadosANOVA.d4.b3, 4, 3)
decomp.d4.b4 <- decomposition.hierarq.anova(dadosANOVA.d4.b4, 4, 4)
decomp.d4.b5 <- decomposition.hierarq.anova(dadosANOVA.d4.b5, 4, 5)
decomp.d4.b6 <- decomposition.hierarq.anova(dadosANOVA.d4.b6, 4, 6)

decomp.d5.b1 <- decomposition.hierarq.anova(dadosANOVA.d5.b1, 5, 1)
decomp.d5.b2 <- decomposition.hierarq.anova(dadosANOVA.d5.b2, 5, 2)
decomp.d5.b3 <- decomposition.hierarq.anova(dadosANOVA.d5.b3, 5, 3)
decomp.d5.b4 <- decomposition.hierarq.anova(dadosANOVA.d5.b4, 5, 4)
decomp.d5.b5 <- decomposition.hierarq.anova(dadosANOVA.d5.b5, 5, 5)
decomp.d5.b6 <- decomposition.hierarq.anova(dadosANOVA.d5.b6, 5, 6)

decomp.d6.b1 <- decomposition.hierarq.anova(dadosANOVA.d6.b1, 6, 1)
decomp.d6.b2 <- decomposition.hierarq.anova(dadosANOVA.d6.b2, 6, 2)
decomp.d6.b3 <- decomposition.hierarq.anova(dadosANOVA.d6.b3, 6, 3)
decomp.d6.b4 <- decomposition.hierarq.anova(dadosANOVA.d6.b4, 6, 4)
decomp.d6.b5 <- decomposition.hierarq.anova(dadosANOVA.d6.b5, 6, 5)
decomp.d6.b6 <- decomposition.hierarq.anova(dadosANOVA.d6.b6, 6, 6)


listaANOVA <- list(decomp.d1.b1, decomp.d1.b2, decomp.d1.b3, decomp.d1.b4, decomp.d1.b5, decomp.d1.b6,
                   decomp.d2.b1, decomp.d2.b2, decomp.d2.b3, decomp.d2.b4, decomp.d2.b5, decomp.d2.b6,
                   decomp.d3.b1, decomp.d3.b2, decomp.d3.b3, decomp.d3.b4, decomp.d3.b5, decomp.d3.b6,
                   decomp.d4.b1, decomp.d4.b2, decomp.d4.b3, decomp.d4.b4, decomp.d4.b5, decomp.d4.b6,
                   decomp.d5.b1, decomp.d5.b2, decomp.d5.b3, decomp.d5.b4, decomp.d5.b5, decomp.d5.b6,
                   decomp.d6.b1, decomp.d6.b2, decomp.d6.b3, decomp.d6.b4, decomp.d6.b5, decomp.d6.b6)






df<-data.frame(z=runif(10),x=as.factor(c(rep("a",5), rep("b",5))),y=as.factor(c(rep("A",2), rep("B",2),rep("C",2),rep("D",4))) )
j1 <- decomposition.hierarq.anova(df,1,1)
j2 <- decomposition.hierarq.anova(df,2,2)
j3 <- decomposition.hierarq.anova(df,2,3)

j <- list(j1,j2,j3)
names(j) <- c("a", "b", "c")



#------------------------------------------------------------------
# OLD VERSION WITH 1 FACTOR ONLY

#df$x is the factor (parcel) and df$y is the reflectance

decomposition.one.way.anova <- function(df, data, banda)
{
  means.parcels <- c()
  n.pixels <- c()
  great.mean <- mean(df$reflectancias,na.rm=TRUE)
  
  for (lev in levels(df$parcela))
  {
    means.parcels <- c(means.parcels,mean(df$reflectancias[df$parcela==lev],na.rm=TRUE))
    n.pixels <- c(n.pixels,length(df$reflectancias[!is.na(df$parcela) & df$parcela==lev]))
  }
  
  sqf <- sum(n.pixels*(means.parcels-great.mean)^2)
  sqt <- var(df$reflectancias)*(length(df$reflectancias[!is.na(df$reflectancias)])-1)
  qmf <- sqf/(nlevels(df$parcela)-1)
  qmre <- (sqt-sqf)/(length(df$reflectancias[!is.na(df$reflectancias)])-nlevels(df$parcela))
  #print(paste("QMF=",qmf))
  #print(paste("QMRE=",qmre))
  #print(paste("F=",qmf/qmre))
  return(c(data, banda, qmf, qmre))
}

decomp.d1.b1 <- decomposition.one.way.anova(dadosANOVA.d1.b1, 1, 1)
decomp.d1.b2 <- decomposition.one.way.anova(dadosANOVA.d1.b2, 1, 2)
decomp.d1.b3 <- decomposition.one.way.anova(dadosANOVA.d1.b3, 1, 3)
decomp.d1.b4 <- decomposition.one.way.anova(dadosANOVA.d1.b4, 1, 4)
decomp.d1.b5 <- decomposition.one.way.anova(dadosANOVA.d1.b5, 1, 5)
decomp.d1.b6 <- decomposition.one.way.anova(dadosANOVA.d1.b6, 1, 6)

decomp.d2.b1 <- decomposition.one.way.anova(dadosANOVA.d2.b1, 2, 1)
decomp.d2.b2 <- decomposition.one.way.anova(dadosANOVA.d2.b2, 2, 2)
decomp.d2.b3 <- decomposition.one.way.anova(dadosANOVA.d2.b3, 2, 3)
decomp.d2.b4 <- decomposition.one.way.anova(dadosANOVA.d2.b4, 2, 4)
decomp.d2.b5 <- decomposition.one.way.anova(dadosANOVA.d2.b5, 2, 5)
decomp.d2.b6 <- decomposition.one.way.anova(dadosANOVA.d2.b6, 2, 6)

decomp.d3.b1 <- decomposition.one.way.anova(dadosANOVA.d3.b1, 3, 1)
decomp.d3.b2 <- decomposition.one.way.anova(dadosANOVA.d3.b2, 3, 2)
decomp.d3.b3 <- decomposition.one.way.anova(dadosANOVA.d3.b3, 3, 3)
decomp.d3.b4 <- decomposition.one.way.anova(dadosANOVA.d3.b4, 3, 4)
decomp.d3.b5 <- decomposition.one.way.anova(dadosANOVA.d3.b5, 3, 5)
decomp.d3.b6 <- decomposition.one.way.anova(dadosANOVA.d3.b6, 3, 6)

decomp.d4.b1 <- decomposition.one.way.anova(dadosANOVA.d4.b1, 4, 1)
decomp.d4.b2 <- decomposition.one.way.anova(dadosANOVA.d4.b2, 4, 2)
decomp.d4.b3 <- decomposition.one.way.anova(dadosANOVA.d4.b3, 4, 3)
decomp.d4.b4 <- decomposition.one.way.anova(dadosANOVA.d4.b4, 4, 4)
decomp.d4.b5 <- decomposition.one.way.anova(dadosANOVA.d4.b5, 4, 5)
decomp.d4.b6 <- decomposition.one.way.anova(dadosANOVA.d4.b6, 4, 6)

decomp.d5.b1 <- decomposition.one.way.anova(dadosANOVA.d5.b1, 5, 1)
decomp.d5.b2 <- decomposition.one.way.anova(dadosANOVA.d5.b2, 5, 2)
decomp.d5.b3 <- decomposition.one.way.anova(dadosANOVA.d5.b3, 5, 3)
decomp.d5.b4 <- decomposition.one.way.anova(dadosANOVA.d5.b4, 5, 4)
decomp.d5.b5 <- decomposition.one.way.anova(dadosANOVA.d5.b5, 5, 5)
decomp.d5.b6 <- decomposition.one.way.anova(dadosANOVA.d5.b6, 5, 6)

decomp.d6.b1 <- decomposition.one.way.anova(dadosANOVA.d6.b1, 6, 1)
decomp.d6.b2 <- decomposition.one.way.anova(dadosANOVA.d6.b2, 6, 2)
decomp.d6.b3 <- decomposition.one.way.anova(dadosANOVA.d6.b3, 6, 3)
decomp.d6.b4 <- decomposition.one.way.anova(dadosANOVA.d6.b4, 6, 4)
decomp.d6.b5 <- decomposition.one.way.anova(dadosANOVA.d6.b5, 6, 5)
decomp.d6.b6 <- decomposition.one.way.anova(dadosANOVA.d6.b6, 6, 6)

resultado <- rbind(decomp.d1.b1, decomp.d1.b2, decomp.d1.b3, decomp.d1.b4, decomp.d1.b5, decomp.d1.b6,
      decomp.d2.b1, decomp.d2.b2, decomp.d2.b3, decomp.d2.b4, decomp.d2.b5, decomp.d2.b6,
      decomp.d3.b1, decomp.d3.b2, decomp.d3.b3, decomp.d3.b4, decomp.d3.b5, decomp.d3.b6,
      decomp.d4.b1, decomp.d4.b2, decomp.d4.b3, decomp.d4.b4, decomp.d4.b5, decomp.d4.b6,
      decomp.d5.b1, decomp.d5.b2, decomp.d5.b3, decomp.d5.b4, decomp.d5.b5, decomp.d5.b6,
      decomp.d6.b1, decomp.d6.b2, decomp.d6.b3, decomp.d6.b4, decomp.d6.b5, decomp.d6.b6)
resultado <- cbind(resultado, resultado[,3]/resultado[,4], sqrt(resultado[,4]))
colnames(resultado) <- c("Data", "Banda", "QMF", "QMRE", "F", "RMSE")







#For latex output
round(resultado[1:6,5:6][,"F"], 1)

res.d1 <- paste0(round(resultado[1:6,5:6][,'F'], 1), ", ", round(resultado[1:6,5:6][,'RMSE'], 3))
res.d2 <- paste0(round(resultado[7:12,5:6][,'F'], 1), ", ", round(resultado[7:12,5:6][,'RMSE'], 3))
res.d3 <- paste0(round(resultado[13:18,5:6][,'F'], 1), ", ", round(resultado[13:18,5:6][,'RMSE'], 3))
res.d4 <- paste0(round(resultado[19:24,5:6][,'F'], 1), ", ", round(resultado[19:24,5:6][,'RMSE'], 3))
res.d5 <- paste0(round(resultado[25:30,5:6][,'F'], 1), ", ", round(resultado[25:30,5:6][,'RMSE'], 3))
res.d6 <- paste0(round(resultado[31:36,5:6][,'F'], 1), ", ", round(resultado[31:36,5:6][,'RMSE'], 3))

res.final <- cbind(res.d1, res.d2, res.d3, res.d4, res.d5, res.d6)

res.final.xtab <- xtable(res.final)
align(res.final.xtab) <- rep("l", 6)
print.xtable(res.final.xtab, booktabs=T, include.rownames=F)


save(resultado, file="resultadoDecompANOVA.object")

RMSE <- sqrt()

rm(dadosANOVA.d6.b1)
rm(dadosANOVA.d6.b2)
rm(dadosANOVA.d6.b3)
rm(dadosANOVA.d6.b4)
rm(dadosANOVA.d6.b5)
rm(dadosANOVA.d6.b6)



#%%%%%%%%%%%%%%%%%%%%%%%%%%
# APPROACH 3) 

MEGAdadosANOVA <- constroiDadosANOVA(data=5, banda=4, dimAmostra=11582)
somaDeQuadrados(MEGAdadosANOVA)
somaDeQuadrados(dadosANOVA)



#%%%%%%%%%%%%%%%%%%%%%%%%%%
# APPROACH 4)

# Test if the crop affects combinations of image/date and NDVIs
FValues <- c()
pValues <- c()
efeito <- c()
for(i in 7:length(listaTodasParcelas))
{
  parcelas.aov <- aov(listaTodasParcelas[,i] ~ as.factor(listaTodasParcelas$cultura))
  FV <- summary(parcelas.aov)[[1]][["F value"]][[1]]
  pV <- summary(parcelas.aov)[[1]][["Pr(>F)"]][[1]]
  
  FValues <- c(FValues,FV)
  pValues <- c(pValues,pV)
  if(pV <= 0.05) ef <- 1 else ef <- 0
  efeito <- c(efeito,ef)
}
nn <- colnames(listaTodasParcelas[,7:length(listaTodasParcelas)])
resultadoANOVA <- data.frame(cbind(nn,FValues,pValues,efeito))
resultadoANOVA$FValues <- as.numeric(as.character(resultadoANOVA$FValues))
resultadoANOVA <- resultadoANOVA[order(-resultadoANOVA$FValues),]; resultadoANOVA




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4. VARIABLE SELECTION                                                               ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4.1 WHOLE DATASET

# Prparing data for the classifier in the next step. Only spectral signatures for the first 6 dates
dadosClassificadores <- listaTodasParcelas[,c(5,7:42)]
dadosClassificadores$cultura <- as.factor(dadosClassificadores$cultura)
dadosClassificadores[,-1][dadosClassificadores[,-1] > 1] <- 1


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 4.2 PARTIAL DATASET (DATASET 2) (This is where variable selection happens)

# Using trim.matrix from subselect package
dadosTrim <- dadosClassificadores
criterioExclusaoVars <- 0.02
classTrim <- trim.matrix(cor(dadosTrim[,-1]), criterioExclusaoVars);classTrim

classTrim$names.discarded == classTrimOLD$names.discarded

varsRetirar <- classTrim$numbers.discarded+1  #+1 because in the original data frame, col 1 is for the crop
varsSobram <- 1:length(dadosTrim)
varsSobram <- varsSobram[! varsSobram %in% varsRetirar]

# Preparing DATASET 2 (after variable selection)
dadosClassificadoresSub <- dadosTrim[,varsSobram]
dadosClassificadoresSub$cultura <- as.factor(dadosClassificadoresSub$cultura)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5. 6. 7. ESCOLHA/TREINO/VALIDACAO DOS CLASSIFICADORES                               ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5.1 KNN

#%%%%%%%%%%%%%%%%%%%%%%%%
# DATASET 1

# Tuning
resTune.comp.knn <- matrix(nrow=21, ncol=0)
for(i in 1:10)
{
  print(i)
  KNN.tune.comp <- tune.knn(dadosClassificadores[,-1], dadosClassificadores[,1], k=1:20)
  k <- KNN.tune.comp[1][[1]][1,1]
  erros <- KNN.tune.comp$performances$error
  resTune.comp.knn <- cbind(resTune.comp.knn, c(k, erros))
}
# RESULT: use k=7, because it produces the lowest error

# Cross validation
KNN.comp.cruz <- validacaoCruzada(n = 10, tipo = 'KNN', dados = dadosClassificadores, lambda = 0.8, k = 7)
KNN.comp.cruz$result$correcTot


#%%%%%%%%%%%%%%%%%%%%%%%%
#DATASET 2

# Tuning
resTune.sub.knn <- matrix(nrow=21, ncol=0)
for(i in 1:10)
{
  print(i)
  KNN.tune.sub <- tune.knn(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], k=1:20)
  k <- KNN.tune.sub[1][[1]][1,1]
  erros <- KNN.tune.sub$performances$error
  resTune.sub.knn <- cbind(resTune.sub.knn, c(k, erros))
}
# RESULT: use k=7, because it produces the lowest error

# Cross validation
KNN.sub.cruz <- validacaoCruzada(n = 10, tipo = 'KNN', dados = dadosClassificadoresSub, lambda = 0.8, k = 7)

KNN.sub.cruz$result$correcTot


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 5.2 SVM

#%%%%%%%%%%%%%%%%%%%%%%%%
# DATASET 1

# Tuning
rGamma <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4)
rCost <- c(0.5, 1, 1.5, 2, 2.5, 3)
resTune.comp.v1 <- matrix(ncol = length(rCost), nrow = length(rGamma))

for(i in 1:length(rGamma))
  for(j in 1:length(rCost))
  {
    print(paste0("i = ", i, " and j = ", j))
    res <- SVM.tune.comp <- tune.svm(dadosClassificadores[,-1], dadosClassificadores[,1], gamma = rGamma[i], cost = rCost[j])
    
    resTune.comp.v1[i,j] <- res$best.performance
  }
colnames(resTune.comp.v1) <- rCost
rownames(resTune.comp.v1) <- rGamma
min(resTune.comp.v1)


# Cross validation
SVM.comp.cruz <- validacaoCruzada(n = 10, tipo = 'SVM', dados = dadosClassificadores, lambda = 0.8, gamma = 0.4, cost = 1.5)

SVM.comp.cruz$result$correcTot


#%%%%%%%%%%%%%%%%%%%%%%%%
#DATASET 2

# Tuning
rGamma <- c(0.4, 0.6, 0.8, 1, 1.2, 1.4)
rCost <- c(0.5, 1, 1.5, 2, 2.5, 3)
resTune.sub.v1 <- matrix(ncol = length(rCost), nrow = length(rGamma))

for(i in 1:length(rGamma))
  for(j in 1:length(rCost))
  {
    print(paste0("i = ", i, " and j = ", j))
    res <- SVM.tune.sub <- tune.svm(dadosClassificadoresSub[,-1], dadosClassificadoresSub[,1], gamma = rGamma[i], cost = rCost[j])

    resTune.sub.v1[i,j] <- res$best.performance
  }
colnames(resTune.sub.v1) <- rCost
rownames(resTune.sub.v1) <- rGamma
min(resTune.sub.v1)


# Cross validation
SVM.sub.cruz <- validacaoCruzada(n = 10, tipo = 'SVM', dados = dadosClassificadoresSub, lambda = 0.8, gamma = 0.4, cost = 2.5)



#%%%%%%%%%%%%%%%%%%%%%%%%
# METHOD CALIBRATION
# q_j estimation. Follow the estimaQ function in functions.R
SVM.sub.qi <- estimaQ(classificador = SVM.sub.cruz$classificador, lambdas = c(0.6, 0.7, 0.8, 0.9, 0.95))



# CALIBRATION ASSESSMENT
# Follow the assessmentMetodo function in functions.R
SVM.sub.assess <- assessmentMetodo(classificador = SVM.sub.cruz$classificador,  qis = SVM.sub.qi)




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Final data structures: output for paper

#%%%%%%%%%%%%%%%%%%%%%%%%
# Data.frame containing all parcels and data from calibration mehtod 
comp1 <- SVM.sub.analise$resultClass
comp1$parcela <- as.character(comp1$parcela)
comp1 <- comp1[order(comp1$parcela, decreasing = TRUE),]

comp2 <- listaTodasParcelas
comp2$id <- as.character(comp2$id)
comp2 <- comp2[order(comp2$id, decreasing = TRUE),]

igual <- comp1$verdade == comp1$class1
igual[igual == TRUE]  <- 'Igual'
igual[igual == FALSE]  <- 'Diferente'

qj <- c()
for(i in 1:nrow(comp1)) qj[i] <- SVM.sub.analise$l0.8$qi[comp1$class1[i]]

resultado <- comp1$prob1 >= qj
resultado[resultado == TRUE]  <- 'Aceite'
resultado[resultado == FALSE]  <- 'Rejeitado'

dadosFinais <- cbind(comp1, igual, qj, resultado, area=comp2$area)

head(dadosFinais)
