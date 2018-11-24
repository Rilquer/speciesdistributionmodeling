#Workflow do Script
##1. Calcular matriz de dist?ncia geogr?fica entre os pontos, ou importar essa matriz
##de outro programa.
##2. Calcular para cada localidade quantas localidades est?o a menos de x km de dist?ncia,
##salvando em arquivo "quant".
##3. Ordenar as localidades por aquelas que est?o mais pr?ximas de outras (ordenar quant
##de modo decrescente)
##5. Excluir a primeira localidade do arquivo "quant" ordenado do arquivo de entrada
##original.
##6. Repetir de 1 a 6, at? que cada localidade esteja pr?xima de apenas uma outra.
##6. Criar um arquivo com os pares de localidades e, guidando-se por ele, excluir
##uma localidade de cada par, de modo aleat?rio.

##Instalando pacotes necess?rios
install.packages('sp','rgal','raster','rgeos')
##Definir caminho do arquivo


##Criar fun??o para filtrar os pontos. Rodar primeiro o c?digo daqui at? a linha 115.
geodel <- function(x,m) {

#Carregando pacotes necess?rios
library(sp)
library(rgdal)
library(raster)
library(rgeos)

#Criando fun??o que calcula a dist?ncia
distcalc <- function(x,m) { 
  ##Define quais as colunas que tem as coordenadas
  coordinates(x) <- c('long', 'lat')
  ##Define o sistema de coordenadas
  proj4string(x) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
  ##Transforma o data.frame em um objeto SpatialPoints.
  locs <- spTransform(x, CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
  ##Verifica quais locaidades distam entre si menos do que a dist?ncia a definida
  ##(argumento dist), retornando TRUE para tais localidades.
  dist<-gWithinDistance(locs, dist = m, byid = TRUE)
  return(dist)
}

#Puxando tabela.

##Puxar arquivo txt, formato long-lat, coordenadas decimas, primeira linha como
##cabe?alho
coord<-read.table(x,header=T)
##Transforma arquivo txt em um dataframe
coord <- data.frame(coord) #Criando arquivo data frame de coordenadas originais
points <- data.frame(coord) #Criando arquivo data frame para trabalho.
nloci<-nrow(coord)
cat("Calculando matriz de dist?ncia entre ",nloci," pontos...")

#Parte 1. Calcular a matriz de dist?ncia, fazer arquivo cont e remover localidades
#com mais cont at? que n?o haja localidades pareadas com mais de uma.
nremoved=0
repeat {
  distmat<-distcalc(points,m) ##Calcula matriz de dist?ncia
  diag(distmat) <- NA ##Define a diagona como NA.
  nloc<-ncol(distmat)
  quant<-c()
  for (i in 1:nloc)
  {
    q<-c(i,(nrow(subset(distmat,distmat[,i]==TRUE))))
    quant<-rbind(quant,q)
  }
  colnames(quant)<-c("loc","quant") #Nomeando colunas
  quant<-data.frame(quant) ##Transformando em data.frame
  quant<-quant[order(-quant$quant),] ##Ordenando de modo descrescente
  ##Testando se todos os valores na coluna quant s?o menores ou iguais a 1. Se sim,
  ##parar o loop. Se n?o, pegar a localidade com maior quantidade de pares pr?ximos
  ##e remover do data frame original, utilizado no c?lculo da dist?ncia.
  if (all(quant$quant<=1)) {
    break  
  } else { 
    remove<-quant[1,1]
    points<-points[-remove,]
    nremoved<-nremoved+1
    cat(".")
    }
}

##Caso o else tenha sido ativado, o loop n?o termina e a dist?ncia ? calculada novamente.
cat(".\n")
cat("\n")
cat("Removendo pontos...")

#Parte 2. Criar arquivos de pares, para remover uma localidade de cada par.
pares<-c()
##Retira a parte triangular inferior da matriz e salva em matpairs.
distmat[lower.tri(distmat, diag=TRUE)]<-NA
nloc<-nrow(points)

##Loop for para criar arquivo pares.
for (i in 2:nloc) {
  for (j in 1:(i-1)) {
    if (distmat[j,i]==TRUE) {
      p<-c(i,j)
      pares<-rbind(pares,p)
    }
  }
}
colnames(pares)<-c("loc1","loc2") #Nomeando colunas
prow<-nrow(pares) #Pegando n?mero de linhas do objeto pares
vecpares<-c()
for (i in 1:prow) {
  colremove<-sample(c(1:2),1)
  vecpares<-rbind(vecpares,pares[i,colremove])
  nremoved=nremoved+1
  cat(".")
}
points<-points[-vecpares,]
write.table(points,paste(substr(x,1,nchar(x)-4),'_',m,'m_filtered.txt',sep=''),row.names=F,col.names = F,quote=F,sep='\t')
cat("\n")
cat(paste("Arquivo filtrado com sucesso. Foram removidos ",nremoved," pontos. As ",nloci-nremoved," coordenadas restantes \n foram salvas no arquivo '",substr(x,1,nchar(x)-4),"_",m,"_filtered.txt'.",sep=''))
}
##Rodar primeiro at? aqui. Depois, rodar a linha abaixo, colocando o nome do
##seu arquivo e a dist?ncia m?nima em metros.

geodel('cyanocorax_cyanopogon.txt',10000)
