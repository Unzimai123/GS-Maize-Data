#carregando pacotes necessarios
library(tidyverse)
library(kableExtra)
library(AGHmatrix)

#carregando banco de dados
info <- read.csv(file = 'info.csv', header = T, sep = ',')
mkr <- read.csv(file = 'markers.csv', header = T, sep = ',')
set.seed(3)
mkr <- mkr[,sample(ncol(mkr), 200)]
data <- cbind(info,mkr)
data <- data[sample(nrow(data),5000),]

#exluindo NAs
duble_ALT <- which(grepl(",", data$ALT))
data <- data[-duble_ALT,]
mkr <- data[,9:208]

#codificando os valores
for (i in 1:ncol(mkr)){
  mkr[,i] <-  replace(mkr[,i], mkr[,i] == "0/0", 2)
  mkr[,i] <-  replace(mkr[,i], mkr[,i] == "1/0", 1) 
  mkr[,i] <-  replace(mkr[,i], mkr[,i] == "0/1", 1)
  mkr[,i] <-  replace(mkr[,i], mkr[,i] == "1/1", 0)
}

#porcentagem de NA
count_na <- sum(is.na(mkr))
count_snp <- (ncol(mkr)*nrow(mkr))
porcent_na <- count_na/count_snp

#MAF
ac <- numeric() 
an <- numeric()

for (i in 1:nrow(mkr)){
  ac[i] <- 0
  an[i] <- 0
  for (j in 1:ncol(mkr)){
    if (!is.na(mkr[i,j])){
      if (mkr[i,j]== 2){
      ac[i] <- ac[i] + 0
      an[i] <- an[i] + 2
    }else if (mkr[i,j] == 1){
      ac[i] <- ac[i] + 1
      an[i] <- an[i] + 2
    }else if (mkr[i,j] == 0){
      ac[i] <- ac[i] + 2
      an[i] <- an[i] + 2
    }
    }
  }
}

MAF <- vector()

for (i in 1:length(ac)){
  if (ac[i] < 0.5*an[i]){
    MAF[i] <- ac[i]/an[i]
  }
  else{
    MAF[i] <- (an[i]-ac[i])/an[i]
  }
}

data <- cbind(data[,1:7],AC = ac, AN = an, MAF = MAF, data[,9:208])

#Calculando Gmatrix (VanRaden)
SNP <- t(mkr)
colnames(SNP) <- c(data$ID)
rownames(SNP) <- c(colnames(mkr))
SNP <- as.numeric(as.matrix(mkr))
SNP <- matrix(SNP, nrow = ncol(mkr), ncol = nrow(mkr), byrow = TRUE)

coef_par <- Gmatrix(SNPmatrix = SNP, method = "VanRaden", missingValue = NA)

heatmap(coef_par)
