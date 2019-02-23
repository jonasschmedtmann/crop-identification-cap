#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Master's Thesis - Remote Sensing                      %
# Environmental Engineering - ISA/UL - Lisbon, Portugal %
# (c) 2014 by Jonas Schmedtmann & Manuel Campagnolo     %
#                                                       %
# FUNCTIONS SCRIPT                                      %
#                                                       %
# This script contains the functions used by the main   %
# sript. Each function is documented (in portuguese).   %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#1. Funcao init:
#   input:  - NA
#   output: - NA
init  <-  function()
{
  #Carregamento de packages
  library(maptools)
  library(raster)
  library(sp)
  library(rgeos)
  library(rgdal)
  library(proj4)
  library(xlsx)
  library(MASS)
  library(data.table)
  library(e1071)
  library(class)
  #library(tree)
  #library(party)
  #library(rpart)
  library(subselect)
  library(plotrix)
  library(kknn)
  library(tikzDevice)
  library(xtable)
  library(caret)
  
  #Definicao de codigos de culturas
  codigos2005 <<- read.delim("culturas2005.txt",header=T)
  codigosNovos2005 <<- read.delim("culturasNovas2005.txt",header=T)
  
  #Definicao de constantes
  source('constants.R')
  cat("\014")
}


#2. Funcao carregaRecortaDadosEspaciais:
#   input:  - coordsAreaEstudo: matriz de coordenadas X,Y que vai ser a area de estudo
#           - crs: string de projeccao da classe CRS
#   output: - lista com os diferentes cdg's transformados e recortados pela area de estudo (apenas as parelas)
carregaRecortaDadosEspaciais <- function(coordsAreaEstudo, crs)
{
  #Ler os dados dos ficheiros shp
  p2005 <- readShapePoly(paste0(CAMINHO.SHP,"/2005/2005_PARC"), proj4string=PROJ4.IGEOE)
  c <- readShapePoly(paste0(CAMINHO.SHP,"/Enquadramentos PT IGEO/CONCELHOS"), proj4string=PROJ4.ETRS)
  
  #Transformar cdg's para sistema de coordenadas dado no argumento crs
  p2005Trans <- spTransform(p2005, crs)
  cTrans <- spTransform(c, crs)
  
  #Construir area de estudo a partir da matriz de entrada
  areaEstudo <- SpatialPolygons(list(Polygons(list(Polygon(coordsAreaEstudo)), "1")),proj4string=crs)
  
  #Recortar as parcelas pela area de estudo  
  p2005inter <- gWithin(p2005Trans, areaEstudo, byid=TRUE)
  p2005clip <- subset(p2005Trans,as.vector(p2005inter))
  p2005clip@data$ID <- 1:length(p2005clip@data$ID)  
  
  out <- list(parc2005=p2005clip, area=areaEstudo, concelhos=cTrans)
  return(out)
}


#4. Funcao corrigeLeConjuntoImagensLandsat:
#   input:  - conjuntoPath: caminho para o conjunto de imagens em formato raster a corrigir
#           - areaEstudo:   cdg com poligono pelo qual as imagens serao cortadas
#           - prefixo:      string que e juntada ao nome original da imagem, na imagem corrigida
#           - corrige:      FALSE se nao for para corrigir as imagens, mas sim ler imagens ja corrigidas
#   output: - lista de imagens corrigidas com metodo DOS4 em formato raster, essas imagens tambem sao escritas em
#             ficheiros
corrigeLeConjuntoImagensLandsat <- function(conjuntoPath, areaEstudo=FALSE, prefixo, corrige=TRUE)
{
  out <- list()
  splitConjuntoPath <- strsplit(conjuntoPath, "/")
  codigoImagens <- splitConjuntoPath[[1]][length(splitConjuntoPath[[1]])]
  
  #O QUE FAZER SE FOR PARA CORRIGIR OS DADOS (DEFAULT)
  if (corrige == TRUE)
  {
    #Metadados
    metadataFile <- paste0(conjuntoPath,"/",codigoImagens,"_MTL.txt")
    metadata <- read.table(metadataFile, nrows = 186)
    
    resultadosDOS4 <- data.frame(Imagem=NA,k=NA,tau=NA,T_v=NA,T_z=NA,E_down=NA,L_p=NA,ro_min=NA,ro_max=NA,ro_mean=NA)
    
    #Corrige cada imagem do conjunto dado
    for(j in 1:length(BANDAS.LANDSAT))
    {
      print(paste0('A processar imagem ', j, ' de ', length(BANDAS.LANDSAT)))
      #Abrir GeoTIFF
      nomeImagem <- paste0(codigoImagens,"_B",BANDAS.LANDSAT[j],".TIF")
      ficheiroImagem <- paste0(conjuntoPath,"/",nomeImagem)
      imagem <- raster(ficheiroImagem)
      
      #if (areaEstudo == FALSE)
      if(3 > 23)
        imagemArea <- imagem
      else
      {
        #Fazer clip do GEOTIFF pela area de estudo
        imagemArea <- crop(imagem, areaEstudo)
        
        #Mascarar valores que estao fora da area de estudo como NA
        imagemArea <- mask(x=imagemArea, mask=areaEstudo, updatevalue=NA)
      }
      
      #Converter valores de 0 (FAIXAS PRETAS) para NA
      imagemArea[imagemArea == 0] <- NA
      
      #1) Radiancia no sensor
      #A banda 7 e diferente porque existem duas bandas 6
      if (BANDAS.LANDSAT[j] == 7)
      {
        Lmax <- as.numeric(as.character(metadata$V3[83]))
        Lmin <- as.numeric(as.character(metadata$V3[84]))
        Qmax <- as.numeric(as.character(metadata$V3[103]))
        Qmin <- as.numeric(as.character(metadata$V3[104]))
      } else {
        Lmax <- as.numeric(as.character(metadata$V3[69+2*BANDAS.LANDSAT[j]-2]))
        Lmin <- as.numeric(as.character(metadata$V3[69+2*BANDAS.LANDSAT[j]-1]))
        Qmax <- as.numeric(as.character(metadata$V3[89+2*BANDAS.LANDSAT[j]-2]))
        Qmin <- as.numeric(as.character(metadata$V3[89+2*BANDAS.LANDSAT[j]-1]))
      }
      
      LSat <- Lmin+(Lmax-Lmin)*((imagemArea-Qmin)/(Qmax-Qmin))
      
      #2) Correccao atmosferica - DOS4
      #Estimacao iterativa de T_v, T_z, E_down, L_p
      #LSatMin <- cellStats(LSat, min)
      LSatMin <- LSat@data@min
      doyStr <- substr(as.character(metadata$V3[5]),14,16)
      doy <- as.numeric(doyStr)
      d <- (1-0.017^2)/(1+0.017*cos(365.25/360*doy*(pi/180)))
      E0 <- ESOLAR.ETM7[j]/d^2
      tetaElev <- as.numeric(as.character(metadata$V3[62]))
      tetaz <- 2*pi/360*(90-tetaElev)
      tetav <- 0
      
      #tetav <- 2*pi/360*(90-tetaElev)
      #tetaz <- 0
      
      Tv <- 1
      Tz <- 1
      EDown <- 0
      tau <- 0
      k <- 0
      repeat
      {
        Lp <- LSatMin-0.01*(E0*cos(tetaz)*Tz+EDown)*Tv/pi
        EDown <- pi*Lp
        oldTau <- tau
        tau <- -cos(tetaz)*log(1-((4*pi*(LSatMin-0.01*(E0*cos(tetaz)*Tz+EDown)*Tv/pi))/(E0*cos(tetaz))))
        Tv <- exp(-tau/cos(tetav))
        Tz <- exp(-tau/cos(tetaz))
        k <- k+1
        if(abs(oldTau-tau)<0.000001){
          break
        }
      }
      if (Tv > 1) Tv <- 1
      if (Tz > 1) Tz <- 1
      if (EDown < 0) EDown <- 0
      if (Lp < 0) Lp <- 0
      
      #3) Reflectancia (objecto da classe Raster contendo imagem com reflectancia)
      reflect <- (pi*(LSat-Lp))/(Tv*(E0*cos(tetaz)*Tz+EDown))   
      reflect[reflect<0] <- 0
      reflect[reflect>1] <- 1
      
      #Resultados para data.frame resultadosDOS4
      #roMin <- cellStats(reflect, min)
      #roMax <- cellStats(reflect, max)
      #roMean <- cellStats(reflect, mean)
      
      reflectMat <- as.matrix(reflect)
      roMin <- min(reflectMat, na.rm = T)
      roMax <- max(reflectMat, na.rm = T)
      roMean <- mean(reflectMat, na.rm = T)
      
      resultadosDOS4 <- rbind(resultadosDOS4,c(nomeImagem,k,tau,Tv,Tz,EDown,Lp,roMin,roMax,roMean))
      #resultadosDOS4 <- resultadosDOS4[-1,]
      
      #Fazer output do raster
      writeRaster(reflect,paste0(conjuntoPath,"/",prefixo,"_",nomeImagem), format="GTiff", overwrite=T)
      outItemName <- paste0("b",BANDAS.LANDSAT[j])
      out[[outItemName]] <- reflect
    }
    
    #Devolver resultados
    write.xlsx(resultadosDOS4, paste0(conjuntoPath,"/",prefixo,"_output.xlsx"), "Resultados DOS4", showNA=FALSE)
    
  }
  else
  {
    #O QUE FAZER SE FOR APENAS PARA LER OS DADOS *JA CORRIGIDOS* E NAO PARA CORRIGIR
    for(j in 1:length(BANDAS.LANDSAT))
    {
      ficheiroImagem <- paste0(conjuntoPath,"/",prefixo,"_",codigoImagens,"_B",BANDAS.LANDSAT[j],".TIF")
      reflect <- raster(ficheiroImagem)
      outItemName <- paste0("b",BANDAS.LANDSAT[j])
      out[[outItemName]] <- reflect
    }
  }
  
  return(out)
}


leImagemCDR <- function(conjuntoPath, areaEstudo=FALSE)
{
  out <- list()
  codigoImagens <- substr(list.files(conjuntoPath, full.names = FALSE)[2], 1, 21)
  
  for(j in 1:length(BANDAS.LANDSAT))
  {
    #Abrir imagem
    pathImagem <- paste0(conjuntoPath,"/",codigoImagens,"_sr_band",BANDAS.LANDSAT[j],".tif")
    imagem <- raster(pathImagem)
    
    #Fazer clip do raster pela area de estudo
    imagemArea <- crop(imagem, areaEstudo)
    
    #Mascarar valores que estao fora da area de estudo como NA
    imagemArea <- mask(x=imagemArea, mask=areaEstudo, updatevalue=NA)
  
    #Converter valores de 0 (FAIXAS PRETAS) para NA
    #imagemArea[imagemArea == 0] <- NA
    
    
    outItemName <- paste0("b",BANDAS.LANDSAT[j])
    out[[outItemName]] <- imagemArea/10000
  }
  
  return(out)
}


#5. Funcao constroiListaImagensLandsat:
#   input:  - landsatPath: caminho para a directoria que contem todas as datas de Landsat: nesta directoria, as 
#                         pastas para cada imagem TEM que estar ordenadas cronoligicamente
#           - areaEstudo:  cdg com poligono pelo qual as imagens serao cortadas
#           - prefixo:     string que e juntada ao nome original da imagem, na imagem corrigida
#           - ano:        ano para o qual se pretende obter a lista de imagens
#           - corrige: FALSE se nao for para corrigir as imagens, mas sim ler imagens ja corrigidas
#   output: - lista com todas as datas e imagens Landsat em forma de matriz, incluindo detalhes da imagem
constroiListaImagensLandsat <- function(landsatPath, areaEstudo, prefixo, ano=ANO, corrige=FALSE)
{
  conjuntos <- list.files(landsatPath, full.names = FALSE)
  doys <- 0
  out <- list()
  
  #Leitura dos dados a partir dos nomes das directorias
  for(i in 1:length(conjuntos))
    doys[i] <- as.numeric(substr(conjuntos[i], 14, 16))
  
  strAno <- paste0("a", ano)
  out[[strAno]] <- list()
  
  #DATAS E IMAGENS
  for(j in 1:length(doys))
  {
    #if(ano != todosAnos[j]) next
    doy <- doys[j]
    strDoy <- paste0("d", doy)
    #out[[strAno]][[strDoy]][["DR"]] <- corrigeLeConjuntoImagensLandsat(paste0(landsatPath,"/",conjuntos[j]), areaEstudo, prefixo, corrige)
    
    out[[strAno]][[strDoy]][["DR"]] <- leImagemCDR(paste0(landsatPath,"/",conjuntos[j]), areaEstudo)
    
    
    
    
    
    #Obter dados de coordenadas min e max
    img <- out[[strAno]][[strDoy]][["DR"]][["b1"]]
    xmin <- xmin(img)
    ymax <- ymax(img)
    cellsize <- res(img)[1]
    TAMANHO.CELULA <<- cellsize
    out[[strAno]][[strDoy]][["info"]] <- list("xmin"=xmin, "ymax"=ymax,"ymax2"=ymax, "cellsize"=cellsize)
    
    #obter coordenadas de todos os pontos
    coordRas <- coordinates(img)
    lin <- nrow(img)
    col <- ncol(img)
    xCoord <- matrix(coordRas[,1],nrow=lin,ncol=col,byrow=T)
    yCoord <- matrix(coordRas[,2],nrow=lin,ncol=col,byrow=T)
    out[[strAno]][[strDoy]][["X"]] <- xCoord
    out[[strAno]][[strDoy]][["Y"]] <- yCoord
    
    #Finalmente, converter as imagens a matrizes para mais tarde facilitar o processamento
    for (k in 1:length(out[[strAno]][[strDoy]][["DR"]]))
      out[[strAno]][[strDoy]][["DR"]][[k]] <- as.matrix(out[[strAno]][[strDoy]][["DR"]][[k]])
  }
  
  return(out)
}


#6. Funcao constroiListaDados:
#   input:  - ano: ano em analise
#   output: - lista com toda a informacao disponivel sobre parcelas e dados de reflectancias, tendo esta forma:
#               lista$nNIFAP$aANO$pPARCELA$DR$dDATA$bBANDA
#                                         $atributos
#                                         $coordsPoly
constroiListaDados <- function(ano=ANO)
{
  todosNIFAP <- c(cdg$parc2005@data$NIFAP)#, cdg$parc2004@data$NIFAP, cdg$parc2005@data$NIFAP)
  unicosNIFAP <- sort(unique(todosNIFAP))
  listaGlobal <- list()
    
  #Construcao da lista
  for(i in 1:length(unicosNIFAP))
  {
    print(paste0("i = ",i))
    print(proc.time())
    #NIFAP
    NIFAP <- unicosNIFAP[i]
    strNIFAP <- paste0("n",NIFAP)
    listaGlobal[[strNIFAP]] <- list()
    
    #for(j in 1:1)
    for(j in 1:length(ano))
    {
      #ANO
      ano <- ano[j]
      strAno <- paste0("a", ano)
      listaGlobal[[strNIFAP]][[strAno]] <- list()
      
      #PARCELA
      esteCDG <- cdg[[j]]
      select <- esteCDG@data[esteCDG$NIFAP==NIFAP,]
      if (nrow(select) == 0) next
      IDs <- as.vector(as.numeric(select$ID))
      parcelas <- as.vector(select$COD_CRUZAM)
      culturas <- as.vector(select$COD_CULTUR)
      sequeiro <- as.vector(select$COD_SEQ_RE)
      
      for(k in 1:length(parcelas))
      {
        #print(paste0("k = ",k))
        #print(proc.time())
        parcela <- parcelas[k]
        strParcela <- paste0("p", parcela)
        listaGlobal[[strNIFAP]][[strAno]][[strParcela]] <- list()
        esteID <- IDs[k]
        if (sequeiro[k] == "S"){
          esteSeq <- TRUE
        } else {
          esteSeq <- FALSE}
        
        #ATRIBUTOS
        estaArea <- esteCDG@polygons[[esteID]]@Polygons[[1]]@area
        listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["atributos"]] <- list(ID=esteID,
                                                                             CULTURA=culturas[k],
                                                                             AREA=estaArea,
                                                                             SEQEIRO=esteSeq)
        
        #COORDENADAS DAS PARCELAS
        estasCoordenadas <- coordinates(esteCDG@polygons[[esteID]]@Polygons[[1]])
        listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["coordsPoly"]] <- estasCoordenadas
        
        #DR EM DIFERENTES DATAS
        listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["DR"]] <- list()
        
        #Definicao do rectangulo para extrair com base nas coordenadas da parcela
        minX <- min(estasCoordenadas[,1])
        maxX <- max(estasCoordenadas[,1])
        minY <- min(estasCoordenadas[,2])
        maxY <- max(estasCoordenadas[,2])
        
        for(l in 1:length(todasImagens[[j]]))
        {
          #DATAS
          datas <- names(todasImagens[[j]])
          estaData <- datas[[l]]
          numBandas <- length(todasImagens[[j]][[l]][["DR"]])
          #listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["DR"]][[estaData]] <- array(dim=c(1,1,numBandas))
          
          #BANDAS
          for(m in 1:numBandas)
          {
            #print(m);print(proc.time())
            estaImg <- todasImagens[[j]][[l]][["DR"]][[m]]
            
            #Extrair valores de interesse da matriz
            minXglobal <- todasImagens[[j]][[l]][["info"]][["xmin"]]
            maxYglobal <- todasImagens[[j]][[l]][["info"]][["ymax"]]
            cellsize <- todasImagens[[j]][[l]][["info"]][["cellsize"]]
            
            selectMinX <<- floor((minX-minXglobal)/cellsize)
            selectMaxX <<- ceiling((maxX-minXglobal)/cellsize)
            selectMinY <<- floor((maxYglobal-maxY)/cellsize)
            selectMaxY <<- ceiling((maxYglobal-minY)/cellsize)
            
            if (selectMaxX >= 0 && selectMaxY >= 0 && selectMinX >= 0 && selectMinY >= 0)
              imgContent <- estaImg[selectMinY:selectMaxY,selectMinX:selectMaxX]
            
            #Criar um novo array na primeira iteracao
            if (m==1)
              listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["DR"]][[estaData]] <- array(dim=c(nrow(imgContent),ncol(imgContent),numBandas+2))
            
            listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["DR"]][[estaData]][,,m] <- imgContent            
          }
          
          if (selectMaxX >= 0 && selectMaxY >= 0 && selectMinX >= 0 && selectMinY >= 0)
          {
            listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["DR"]][[estaData]][,,m+1] <- todasImagens[[j]][[l]][["X"]][selectMinY:selectMaxY,selectMinX:selectMaxX]
            listaGlobal[[strNIFAP]][[strAno]][[strParcela]][["DR"]][[estaData]][,,m+2] <- todasImagens[[j]][[l]][["Y"]][selectMinY:selectMaxY,selectMinX:selectMaxX]
          }
        }
      }
    }
  }
  return(listaGlobal)
}


#XX. Funcao obtemConteudoParcela
#   input:  - parcela:        parcela sob forma de lista
#           - apenasInterior: TRUE/FALSE a indicar se apenas devem ser considerados os pixeis no interior da parc
#           - excluiVizinhos: TRUE/FALSE que pode seleccionar apenas pixeis que nao se toquem
#   output: - lista que para cada data tem uma matriz que tem uma coluna para cada banda, contendo todas as
#             reflectancias
obtemConteudoParcela <- function(parcela, apenasInterior=TRUE, excluiVizinhos=FALSE)
{
  out <- list()
  poly.x <- parcela$coordsPoly[,1]
  poly.y <- parcela$coordsPoly[,2]
  
  datas <- names(parcela$DR)
  for(i in 1:length(datas))
  {
    #out[[datas[i]]] <- list()
    #bandas <- names(parcela$DR[[datas[i]]])
    numBandas <- length(parcela$DR[[datas[i]]][1,1,])
    for(j in 1:(numBandas-2))  #-2 para nao incluir as "bandas" X e Y
    {
      #O que fazer se queremos excluir os vizinhos
      if(excluiVizinhos==FALSE)
      {
        #Obter coordenadas dos CENTROS dos pixeis
        img.x <- as.vector(parcela$DR[[i]][,,numBandas-1])
        img.y <- as.vector(parcela$DR[[i]][,,numBandas])
        
        #Obter a banda em questao
        bandaVec <- as.vector(parcela$DR[[i]][,,j])
      }
      else
      {
        #Obter coordenadas dos CENTROS dos pixeis
        mat.x <- parcela$DR[[i]][,,numBandas-1]
        img.x <- as.vector(mat.x[seq(1,nrow(mat.x),2),seq(1,ncol(mat.x),2)])
        
        mat.y <- parcela$DR[[i]][,,numBandas]
        img.y <- as.vector(mat.y[seq(1,nrow(mat.y),2),seq(1,ncol(mat.y),2)])
        
        #Obter a banda em questao
        mat <- parcela$DR[[i]][,,j]
        matSemViz <- mat[seq(1,nrow(mat),2),seq(1,ncol(mat),2)]
        bandaVec <- as.vector(matSemViz)
      }
      
      #O que fazer quando nao interessa se o pixel esta todo dentro da parcela ou nao
      if(apenasInterior==FALSE)
      {
        dentro <- point.in.polygon(img.x, img.y, poly.x, poly.y)
        conteudo <- bandaVec[dentro>0]
        
        #Criar uma nova matriz na primeira iteracao
        if(j==1)
          out[[datas[i]]] <- matrix(nrow=length(conteudo),ncol=numBandas-2)
        
        out[[datas[i]]][,j] <- conteudo
      }
      #O que fazer quando queremos apenas os pixeis com 100% dentro da parcela
      else
      {
        img.x.1 <- img.x-TAMANHO.CELULA/2
        img.y.1 <- img.y+TAMANHO.CELULA/2
        dentro.1 <- point.in.polygon(img.x.1, img.y.1, poly.x, poly.y)
        img.x.2 <- img.x+TAMANHO.CELULA/2
        img.y.2 <- img.y+TAMANHO.CELULA/2
        dentro.2 <- point.in.polygon(img.x.2, img.y.2, poly.x, poly.y)
        img.x.3 <- img.x+TAMANHO.CELULA/2
        img.y.3 <- img.y-TAMANHO.CELULA/2
        dentro.3 <- point.in.polygon(img.x.3, img.y.3, poly.x, poly.y)
        img.x.4 <- img.x-TAMANHO.CELULA/2
        img.y.4 <- img.y-TAMANHO.CELULA/2
        dentro.4 <- point.in.polygon(img.x.4, img.y.4, poly.x, poly.y)
        dentroFinal <- dentro.1*dentro.2*dentro.3*dentro.4
        conteudo <- bandaVec[dentroFinal>0]
        
        #Criar uma nova matriz na primeira iteracao
        if(j==1)
          out[[datas[i]]] <- matrix(nrow=length(conteudo),ncol=numBandas-2)
        
        out[[datas[i]]][,j] <- conteudo
      }
    }
  }
  return(out)
}


#XX. Funcao parcNDVI
#   input:  - parcela:        parcela sob forma de lista
#           - apenasInterior: TRUE/FALSE a indicar se apenas devem ser considerados os pixeis no interior da parc
#   output: - lista que para cada data tem o NDVI de cada pixel na parcela
parcNDVI <- function(parcela, apenasInterior=TRUE)
{
  out <- list()
  conteudo <- obtemConteudoParcela(parcela,apenasInterior)
  datas <- names(conteudo)
  
  for(i in 1:length(datas))
  {
    NDVI <- (conteudo[[i]][,4]-conteudo[[i]][,3])/(conteudo[[i]][,4]+conteudo[[i]][,3])
    out[[datas[i]]] <- NDVI
  }
  return(out)
}


#XX. Funcao parcNIFAP
#   input:  - nomeParc: string na forma 'pXXXXXXXXXXXXX.XXX.XX' com nome da parcela cujo NIFAP se pretende saber 
#   output: - NIFAP ao qual pertence a parcela introduzida
parcNIFAP<-function(nomeParc)
{
  for(i in 1:length(listaDados))
    for(j in 1:length(listaDados[[i]]))
      for(k in 1:length(listaDados[[i]][[j]]))
      {
        if (names(listaDados[[i]][[j]])[[k]] == nomeParc)
          return(names(listaDados)[[i]])
      }
}


#XX. Funcao constroiDadosANOVA
#   input:  - ano:        ano em analise
#           - data:       data em analise
#           - banda:      banda espectral em analise
#           - dimAmostra: numero de parcelas que se pretendem analisar
#   output: - data.frame com dados organizados para se efectuar uma ANOVA hierarquizada
constroiDadosANOVA <- function(data, banda, dimAmostra)
{
  dadosANOVA <- matrix(nrow=0,ncol=3)
  set.seed(0)    #Depois e para tirar
  amostraParcs <- sort(sample(1:length(finalListaDados),dimAmostra))
  
  for(i in amostraParcs)
  {
    parc <- finalListaDados[[i]]
    reflectancias <- obtemConteudoParcela(parc,apenasInterior=TRUE,excluiVizinhos=TRUE)[[data]][,banda]
          
    #Passar ao ciclo seguinte se houver um unico NA na parcela
    if(any(is.na(reflectancias))) next
      
    cultura <- rep(parc$atributos$CULTURA,length(reflectancias))
    parcela <- rep(names(finalListaDados)[[i]],length(reflectancias))
    novaParc <- cbind(reflectancias,cultura,parcela)
    dadosANOVA <- rbind(dadosANOVA,novaParc)
    }
  
  #Converter resultado para data.frame com formatos correctos
  dadosANOVA <- as.data.frame(dadosANOVA)
  dadosANOVA$reflectancias <- as.numeric(as.vector(dadosANOVA$reflectancias))
  dadosANOVA$cultura <- as.factor(dadosANOVA$cultura)
  dadosANOVA$parcela <- as.factor(dadosANOVA$parcela)
  
  return(dadosANOVA)
}


somaDeQuadrados <- function(dados)
{
  # Calculo SQT
  SQT <- var(dados$reflectancias)*(length(dados$reflectancias)-1)
  
  # Calculo SQA
  SQA <- 0
  grandeMedia <- mean(dados$reflectancias)
  for(i in 1:12)
  {
    x <- dados[dados$cultura == i,]
    ni <- nrow(x)
    if(nrow(x) == 0) next
    media <- mean(x$reflectancias)
    
    SQA <- SQA + ni*(media - grandeMedia)^2
  }
  
  # Calculo SQRE
  rle <- rle(as.vector(dados$parcela))$lengths
  rleCont <- rep(1:length(rle), rle)
  dados <- cbind(dados, rle=rleCont)
  medias <- tapply(dados$reflectancias, rleCont, mean)
  SQRE <- sum((dados$reflectancias - medias[dados$rle])^2)
  
  # Calculo SQAB
  SQAB <- SQT-SQA-SQRE
  
  return(list(SQT=SQT, SQA=SQA, SQAB=SQAB, SQRE=SQRE))
}




#XX. Funcao constroiTodasParcelas
#   input:  - excluiVazias: TRUE/FALSE a indicar se as parcelas com NA devem ser excluidas
#   output: - data.frame com dados relevantes sobre todas as parcelas existentes na lista de dados
constroiTodasParcelas <- function(excluiVazias=TRUE)
{
  count <- 1
  for(i in 1:length(listaDados))
  #for(i in 1)
    for(j in 1:length(listaDados[[i]]))
      for(k in 1:length(listaDados[[i]][[j]]))
      {
        print(paste0("count = ",count))
        count <- count+1
        
        #Obter dados base para parcela
        parc <- listaDados[[i]][[j]][[k]]
        id <- names(listaDados[[i]][[j]])[k]
        nifap <- names(listaDados)[i]
        ano <- names(listaDados[[i]])[j]
        area <- parc$atributos$AREA/10000
        cultura <- parc$atributos$CULTURA
        seq <- parc$atributos$SEQEIRO
        
        #Obter reflectancias para diferentes datas - media e desvio padrao
        reflectancias <- obtemConteudoParcela(parc,apenasInterior=TRUE)
        matRef <- matRefNomes <- matrix(nrow=length(reflectancias),ncol=ncol(reflectancias[[1]]))
        for(m in 1:length(reflectancias))
          for(n in 1:ncol(reflectancias[[1]]))
          {  
            matRef[m,n] <- mean(reflectancias[[m]][,n], na.rm=T)
            matRefNomes[m,n] <- paste0("B",n,"_d",m)
          }
        
        vecRef <- as.vector(matRef)
        vecRefNomes <- as.vector(matRefNomes)
        
        #Obter NDVI's para diferentes datas
        NDVIs <- parcNDVI(parc)
        vecNDVI <- c()
        for(l in 1:length(NDVIs))
          vecNDVI[l] <- mean(NDVIs[[l]], na.rm=T)
        
        #Saltar esta parcela se houver algum NA e se isso for desejado
        if(excluiVazias == TRUE && (any(is.na(vecRef)) || any(is.na(vecNDVI)) || 0 %in% vecRef || -1 %in% vecNDVI)) next
        
        #Construir linha para data.frame
        if(i==1 && j==1 && k==1)
        {
          nDatas <- length(names(NDVIs))
          nomesCol <- c("id","NIFAP","ano","area","cultura","seq",vecRefNomes,paste0("NDVI_d",1:nDatas))
          out <- matrix(nrow=0,ncol=length(nomesCol))
          colnames(out) <- nomesCol
        }
        linha <- c(id,nifap,ano,as.numeric(area),cultura,seq,vecRef,vecNDVI)
        
        out <- rbind(out,linha)
        
      }
  
  #Construir data.frame com numeros
  out <- as.data.frame(out)
  rownames(out) <- out$id
  out$area <- as.numeric(as.character(out$area))
  out$cultura <- as.numeric(as.character(out$cultura))
  
  for(i in 7:length(out))
    out[,i] <- as.numeric(as.character(out[,i]))
  
  return(out)
}


#XX. Funcao treinaValida
#   input:  - ver 'validacaoCruzada'
#   output: - ver 'validacaoCruzada'
treinaValida <- function(tipo, treino, valid, k=NA, gamma=NA, cost=NA)
{
  #Treino do classificador
  if(tipo == 'KNN')
  {
    classificador <- kknn(cultura ~ ., treino, valid[,-1], k=k, kernel='rectangular')
    probs <- classificador$prob
  }
  else if(tipo == 'SVM')
  {
    classificador <- svm(cultura ~ ., data=treino, probability=TRUE, kernel='radial', gamma=gamma, cost=cost)
    probs <- attr(predict(classificador, valid[,-1], probability=T), "prob")
  }
  
  #Obtencao das probabilidades principais e secundarias 
  probs1 <- as.vector(apply(probs, 1, max))
  probs2 <- as.vector(apply(probs, 1, function(x) max(x[x!=max(x)])))
  
  #Obtencao das classificacoes primarias e secundarias associadas as probabilidades
  class1 <- class2 <- c() 
  for (i in 1:nrow(probs))
  {
    class1[i] <- as.numeric(names(which(probs[i,] == probs1[i])))
    class2[i] <- as.numeric(names(which(probs[i,] == probs2[i])))
  }
  
  #Data.frame com os resultados obtidos
  result <- data.frame(parcela=rownames(valid),
                       verdade=valid[,1],
                       class1=class1,
                       prob1=probs1)
  
  #Matriz de erro e percentagem de classificacoes correctas
  matErro <- as.data.frame.matrix(table(class = factor(class1, levels=levels(valid[,1])), verdade = valid[,1]))  
  correcTot <- cc(matErro)
  
  #User's e producer's accuracy
  userA <- diag(as.matrix(matErro))/rowSums(matErro)
  prodA <- diag(as.matrix(matErro))/colSums(matErro)
  names(userA) <- names(prodA) <- names(matErro)
  
  #Construcao do output
  out<-list(#classificador=classificador,
            result=result,
            matErro=matErro,
            #matErroProb=matErroProb,
            correcTot=correcTot,
            userA=userA,
            prodA=prodA,
            probs=round(probs,3))
}


#XX. Funcao validacaoCruzada
#   input:  - n:           cs
#           - tipo:        'KNN' ou 'SVM'
#           - dados:       totalidade dos dados de parcelas (daqui serao seleccionados treino e validacao)
#           - lambda:      nivel de 'confianca' desejado
#           - tune:        TRUE/FALSE conforme se deseja fazer o tuning
#           - k:           KNN: numero de classes
#           - gamma:       SVM: parametro gamma
#           - cost:        SVM: parametro de custo
#   output: lista com resultados da validacao cruzada efectuada, contendo 2 elementos:
#           - lista resultTodos (resultados em cada uma das n validacoes 'folds')
#           - lista result (media dos resultados mais importantes de resultTodos):
#               - correcTot:    fraccao total de parcelas bem classificadas
#               - userA:        user's accuracy
#               - prodA:        producer's accuracy
#               - k:            numero de classes para classificador KNN
#               - gamma:        para classificador SVM
#               - cost:         para classificador SMV
#               - qi:           limites de probabilidade P(wi|x) que garantem lambda para classe i
validacaoCruzada <- function(n, tipo, dados, lambda, k=NA, gamma=NA, cost=NA)
{
  #Estas sao as particoes que FICAM DE FORA em cada uma das validacoes
  resultTodos <- result <- list()
  set.seed(23)
  particoes <- createFolds(dados$cultura, k = n)
  
  #Fazer n validacoes
  for(i in 1:n)
  {
    print(paste0('A processar validacao ', i, ' de ', n, '...'))
    
    #Seleccao dos dados a partir dos folds criados
    dadosTrein <- dados[-particoes[[i]],]
    dadosValid <- dados[particoes[[i]],]
    
    #Treina e valida o classificador
    classificador <- treinaValida(tipo, dadosTrein, dadosValid, k=k, gamma=gamma, cost=cost)
    
    #Resultados desta validacao
    resultTodos[[paste0('valid',i)]] <- classificador
  }
  
  #Juncao das classificacoes de cada fold
  classificador <- NA
  for(i in 1:n)
    classificador <- rbind(classificador, resultTodos[[i]]$result)
  classificador <- classificador[-1,]
  
  #Construcao dos resultados de output
  nomes <- names(resultTodos$valid1$userA)
  tmp <- unlist(resultTodos, F)
  result[["OA"]] <- as.vector(sapply("correcTot", function(x) rowMeans(na.rm = T, do.call(cbind, tmp[grep(x, names(tmp))]))))
  result[["UA"]] <- as.vector(sapply("userA", function(x) rowMeans(na.rm = T, do.call(cbind, tmp[grep(x, names(tmp))]))))
  names(result[["UA"]]) <- nomes
  result[["PA"]] <- as.vector(sapply("prodA", function(x) rowMeans(na.rm = T, do.call(cbind, tmp[grep(x, names(tmp))]))))
  names(result[["PA"]]) <- nomes
  
  print("Fim")
  return(list(result=result, resultTodos=resultTodos, classificador=classificador))
}


#XX. Funcao estimaQ
#   input:  - ver 'validacaoCruzada'
#   output: - ver 'validacaoCruzada'
estimaQ <- function(classificador, lambdas)
{
  c <- length(table(classificador$verdade))
  
  # Calculate qi limits for different lambas (confidence level)
  qis <- matrix(ncol = length(lambdas), nrow = c)
  for(j in 1:length(lambdas))
  {    
    print(paste0('qi = ', lambdas[j]))
    qi <- c()
    
    # For each crop
    for(i in 1:c)
    {
      RqInicial <- classificador[classificador$class1 == i,]
      RqInicial <- RqInicial[order(RqInicial$prob1, decreasing = F),]
      
      for(k in 1:nrow(RqInicial))
      {
        Rq <- RqInicial[k:nrow(RqInicial),]
        UA <- nrow(Rq[Rq$verdade == Rq$class1,])/nrow(Rq)
        
        if(UA >= lambdas[j])
        {
          qi[i] <- RqInicial[k, 'prob1']
          break
        }
        if(k == nrow(RqInicial)) {qi[i] <- 1; break}
      }
    }
    qis[,j] <- qi
  }
  colnames(qis) <- lambdas
  rownames(qis) <- COD_CULT
  return(qis)
}




#XX. Funcao assessmentMetodo
#   input:  - class:  o classificador resultado duma dada validacao cruzada efectuada
#   output: - ...
assessmentMetodo <- function(classificador, qis)
{
  out <- list()
  lambdas <- colnames(qis)
  for(i in 1:length(lambdas))
  {
    #Seleccao apenas das classificacoes aceites em funcao de lambda
    qi <- qis[,i]
    classAceites <- classificador[classificador$prob1 >= qi[classificador$class1],]
    
    #Matriz de erro e accuracy ststistics
    matErro <- as.data.frame.matrix(table(class = factor(classAceites$class1, levels=levels(classAceites$verdade)), verdade = classAceites$verdade))
    OA <- cc(matErro)
    userA <- diag(as.matrix(matErro))/rowSums(matErro)
    prodA <- diag(as.matrix(matErro))/colSums(matErro)
    
    #Percentagem de parcelas com classificacao aceites
    percAcTot <- nrow(classAceites)/nrow(dadosClassificadores)
    
    #Isto e a versao com a classificacao
    percAcCult <- table(factor(classAceites$class1, levels=levels(classAceites$verdade)))/table(classificador$class1)
    
    out[[paste0('l', lambdas[i])]] <- list(qi=qi, matErro=matErro, OA=OA, UA=userA, PA=prodA, OACP=percAcTot, ACP=percAcCult)
  }
  return(out)
}





#Pequenas funcoes de ajuda
squareTable <- function(x,y) {
  commonLevels <- unique(c(levels(x), levels(y)))
  x <- factor(x, levels = commonLevels)
  y <- factor(y, levels = commonLevels)
  return(table(x,y))
}
cc <- function(m) return(sum(diag(as.matrix(m)))/sum(m))
condensaMatriz <- function(m, ind) {
  x <- y <- seq(0,0,nrow(m))
  m <- as.matrix(m)
  for(i in 1:length(ind))
    x <- x + as.vector(m[ind[i],])
  m[ind[1],] <- x
  
  for(i in 1:length(ind))
    y <- y + as.vector(m[,ind[i]])
  m[,ind[1]] <- y
  
  m<-m[-ind[-1],-ind[-1]]
  print(cc(m))
  return(m)
}
convertCoord <- function(coord){
  deg <- floor(coord)
  x <- (coord - deg) * 60
  min <- floor(x)
  y <- (x - min) * 60
  sec <- floor(y)
  return(paste0(deg, "d ", min, "' ", sec, "''"))
}
jsort.vector <- function(x, sort.by.vector) {
  if(length(sort.by.vector) == length(x)) {
    y <- x[sort.by.vector[1]]
    for(i in 2:length(sort.by.vector))
      y <- cbind(y, x[sort.by.vector[i]])
    return(as.vector(y))
  } else {
    print("Erro")
  } 
}
jsort.data.frame.col <- function(x, sort.by.vector) {
  if(length(sort.by.vector) == ncol(x)) {
    y <- x[,sort.by.vector[1]]
    for(i in 2:length(sort.by.vector))
      y <- cbind(y, x[,sort.by.vector[i]])
    colnames(y) <- jsort.vector(colnames(x), sort.by.vector)
    return(y)
  } else {
    print("Erro")
  } 
}
jsort.data.frame.row <- function(x, sort.by.vector) {
  if(length(sort.by.vector) == nrow(x)) {
    y <- x[sort.by.vector[1],]
    for(i in 2:length(sort.by.vector))
      y <- rbind(y, x[sort.by.vector[i],])
    #rownames(y) <- jsort.vector(rownames(x), sort.by.vector)
    return(y)
  } else {
    print("Erro")
  } 
}


NOVO.trim.matrix <- function (mat, tolval) 
{
  p <- dim(mat)[2]
  matindices <- 1:p
  mat.eig <- eigen(mat, symmetric = TRUE)
  discard <- rep(FALSE, p)
  newmat <- mat
  newmatindices <- matindices
  k <- c()
  while (mat.eig$values[p]/mat.eig$values[1] < tolval) {
    k <- c(k, mat.eig$values[1]/mat.eig$values[p])
    int <- as.numeric(newmatindices[order(abs(mat.eig$vectors[,p]), decreasing = TRUE)[1]])
    
    discard[int] <- TRUE
    newmat <- mat[!discard, !discard]
    newmatindices <- matindices[!discard]
    p <- p - 1
    mat.eig <- eigen(newmat, symmetric = TRUE)
  }
  size <- dim(newmat)[2]
  output <- list(newmat, as.numeric(matindices[discard]), colnames(mat)[discard], 
                 size, k)
  names(output) <- c("trimmedmat", "numbers.discarded", "names.discarded", 
                     "size", "k")
  output
}



