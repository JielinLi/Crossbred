#input
DIR <- ""
library(purrr)
library(plyr)
library(expm)
library(MASS)
library(blupADC)
source(paste0(DIR,"functions.r"))

#-------------------------------------------------------------------------------------------------------------------
# population parameters
#-------------------------------------------------------------------------------------------------------------------
nRep     <- 2 
nBreeds  <- Npop <- 3
nIndPop  <- c(4500, 4500, 4500)  ## The number of elements must equal nBreeds(base population)
nSires0  <- c(2250, 2250, 2250)
nDams0   <- c(2250, 2250, 2250)
NSires   <- c(30,45,60)  ## The number of selected sires; If it is truncation selection(next generation)
NDams    <- c(300,450,600)  ## The number of selected dams;  If it is truncation selection
NAni     <- c(2400,4500,7200) ## The total number of animals in the population (offspring)
nGen     <- 16            ## The number of generations
#nGen  <- c(5,10,20); nGen[ibreed]
#-------------------------------------------------------------------------------------------------------------------
# genome paramters
#-------------------------------------------------------------------------------------------------------------------
nChr    <- 18
chrl    <- 120


#-------------------------------------------------------------------------------------------------------------------
# Read stored haplotypes
#-------------------------------------------------------------------------------------------------------------------
hap.all <- list()

for (ibreed in 1:nBreeds) {
  setwd(paste0(DIR, "breed", ibreed))
  chr.all <- list()
  for (ichr in 1:nChr) {
    chr.all[[ichr]] <- read.table(paste0("chr", ichr, ".txt")); print(dim(chr.all[[ichr]]))
  }
  hap.all[[ibreed]] <- chr.all
}

#check dimension
dim(hap.all[[ibreed]][[ichr]]) #haplotypes are stored as ibreed - ichr


#-------------------------------------------------------------------------------------------------------------------
# trait parameters
#-------------------------------------------------------------------------------------------------------------------
## Three breeds, 2 traits for each.
## we know additive matrix (genomic matrix), residual matrix and dominance matrix
nTraits  <- 6

a.dis    <- "normal"       ## The distribution of the allelic effects; either "gamma" or even
nqtl     <- 200            ## The number of QTL per chromosome
nloci    <- 4092         ## The number of loci per chromosome(marker+QTL)

#variance covariance matrix of additive effects
gList <- list()
gList[[1]] <- matrix(c(1.5, 0.5, 0.5, 1), byrow = TRUE, ncol = nTraits/nBreeds)  
gList[[2]] <- matrix(c(2,0.7,0.7,1), byrow=TRUE, ncol = nTraits/nBreeds)  
gList[[3]] <- matrix(c(2.5,0.9,0.9,2.5), byrow=TRUE, ncol = nTraits/nBreeds)
genomicMatrix <- matrix(0,nrow=nTraits,ncol=nTraits)
for (ibreed in 1:nBreeds) {
  genomicMatrix[(2*(ibreed-1)+1):(2*ibreed),(2*(ibreed-1)+1):(2*ibreed)] <- gList[[ibreed]]
}

#variance covariance matrix of residuals
resList <- list()
resList[[1]] <- matrix(c(0.8, 0, 0, 5.2), byrow = TRUE, ncol = nTraits/nBreeds)  
resList[[2]] <- matrix(c(2.5, 0, 0, 4.5), byrow=TRUE, ncol = nTraits/nBreeds)  
resList[[3]] <- matrix(c(3.2, 0, 0, 13.2), byrow=TRUE, ncol = nTraits/nBreeds)
resMatrix <- matrix(0,nrow=nTraits,ncol=nTraits)
for (ibreed in 1:nBreeds) {
  resMatrix[(2*(ibreed-1)+1):(2*ibreed),(2*(ibreed-1)+1):(2*ibreed)] <- resList[[ibreed]]
}

#variance covariance matrix of dominance effects
domList <- list()
domList[[1]] <- matrix(c(1.3, 0, 0, 1.3), byrow = TRUE, ncol = nTraits/nBreeds)  
domList[[2]] <- matrix(c(1.3, 0, 0, 1.3), byrow=TRUE, ncol = nTraits/nBreeds)  
domList[[3]] <- matrix(c(1.3, 0, 0, 1.3), byrow=TRUE, ncol = nTraits/nBreeds)
domMatrix <- matrix(0,nrow=nTraits,ncol=nTraits)
for (ibreed in 1:nBreeds) {
  domMatrix[(2*(ibreed-1)+1):(2*ibreed),(2*(ibreed-1)+1):(2*ibreed)] <- domList[[ibreed]]
}

dMean   <- rep(0.192,nTraits)     ##mean of dominance degree for each QTL
dSd     <- rep(0.312,nTraits)     ##sa of dominance degree for each QTL


#variance covariance matrix of dominance effects
for (ibreed in 1:nBreeds) {
  domMatrix[(2*(ibreed-1)+1):(2*ibreed),(2*(ibreed-1)+1):(2*ibreed)] <- domList[[ibreed]]
}

ecoValueBreed1 <- matrix(c(0.7,0.3,0,0,0,0),nrow=1) 
ecoValueBreed2 <- matrix(c(0,0,0.5,0.5,0,0),nrow=1) 
ecoValueBreed3 <- matrix(c(0,0,0,0,0.4,0.6),nrow=1) 
ecoValue <- list(ecoValueBreed1,ecoValueBreed2,ecoValueBreed3)
newecoValueBreed1 <- matrix(c(0.7,0.3),nrow=1) 
newecoValueBreed2 <- matrix(c(0.5,0.5),nrow=1) 
newecoValueBreed3 <- matrix(c(0.4,0.6),nrow=1) 
newecoValue <- list(newecoValueBreed1,newecoValueBreed2,newecoValueBreed3)
#genome
#-------------------------------------------------------------------------------------------------------------------
# sample QTL and markers
#-------------------------------------------------------------------------------------------------------------------
LS 	     <- NAni/NDams    ## Litter size 
mating   <- NDams/NSires  ## The number of dams per sire
d        <- seq(0,(chrl-(chrl/nloci)),length.out=nloci)  ## evenly distributed 
BreedDat <- list()

ibreed <- 1 ##redundant, only for testing 

for (ibreed in 1:nBreeds){   ## loop over breeds
  NAnim <- nIndPop[ibreed]
  Nmrk <- NULL      ## an empty vector to store the number of markers in each chromosome, in the current the number was given 
  Wmrk <- NULL      ## an empty vector to store which ones are the markers
  Wqtl <- NULL      ## an empty vector to store which ones are the qtl
  Wneu <- NULL
  chromosome <- list()
  
  for (ichr in 1:nChr) {               ## loop over chromosomes
    dat  <- hap.all[[ibreed]][[ichr]]  ## put all the chromosomes to a list
    wpqtl <- 1:(ncol(dat)/2);wpqtl
    wqtl    <- wpqtl[floor(seq(1,length(wpqtl),length.out=nqtl))]
    wmrk    <- setdiff(wpqtl,wqtl) ;length(wmrk)          ## which ones are marker (nmrk)
    nmrk    <- nloci-nqtl
    sel.L   <- as.numeric(c(wqtl,wmrk))                   ## which ones are the selected loci (qtl+mrk)
    nkept   <- length(sel.L)                              ## the number of qtl+mrk
    dkept   <- d[sel.L[order(sel.L)]]                     ## the locations of the selected loci on the chromosome
    
    crowmrk <- wcol(wmrk)[order(wcol(wmrk))]              ## which columns in dat are for markers
    colnames(dat)[crowmrk] <- NameMrk(nmrk)               ## name markers
    crowqtl <- wcol(wqtl)[order(wcol(wqtl))]              ## which columns in dat are for qtl
    colnames(dat)[crowqtl] <- Nameqtl(length(wqtl))       ## name qtl
    
    
    Animals <- dat[,wcol(sel.L)[order(wcol(sel.L))]]  ## the haplotypes for the base population
    
    
    ## This part of the code can be re-written
    DATA     <- cbind(rep(dkept,each=2),t(Animals))  
    cqtl     <- which(DATA[,1] %in% d[wqtl]==TRUE)
    wqtl     <- cqtl[is.even(cqtl)]/2                     ## which columns are qtl in the new data (the data with only selected loci)
    
    cmrk     <- which(DATA[,1] %in% d[wmrk]==TRUE)         
    wmrk     <- cmrk [is.even(cmrk)]/2                    ## which columns are markers in the new data
    
    ## which columns are qtl in the new data (the data with only selected loci)
    
    dqtl     <- dkept[wqtl]                               ## the locations of qtl in the new data
    dmrk     <- dkept[wmrk]                               ## the locations of marker in the new data
    
    
    locqtl <- data.frame(chr=rep(ichr,times=length(wqtl)),number=1:length(wqtl), dqtl)   ## The location data of qtl 
    locmrk <- data.frame(chr=rep(ichr,times=nmrk),number=1:nmrk, dmrk)                   ## The location data of mrk
    
    
    qtlFreq  <- getFreq(Animals[,cqtl])                   ## the frequency of allele 1 at qtl
    mrkFreq  <- getFreq(Animals[,cmrk])                   ## the frequency of allele 1 at marker
    
    
    ## writing all the information to a list 
    temp <- list(Animals=Animals,nkept=nkept,dkept=dkept,locqtl=locqtl,locmrk=locmrk,
                 qtlFreq=qtlFreq,mrkFreq=mrkFreq)
    
    chromosome[[ichr]] <- temp                            ## write the results to the list chromosome
    
    Nmrk <- c(Nmrk,nmrk)   ## can be changed (neglected)
    Wmrk <- c(Wmrk,wmrk)
    Wqtl <- c(Wqtl,wqtl)
    
    
    cat(paste("chr ",ichr," is done \n"))
  }                                                    ## loop over chromosomes
  
  Animals <- chromosome[[1]]$Animals
  dkept <- list()
  wmrk <- wqtl <- wneu <- nkept <- qtlFreq <- mrkFreq <- LOCqtl <- LOCmrk <- NULL
  for (i in 1:nChr){
    if (i != 1) {
      Animals    <- cbind(Animals,chromosome[[i]]$Animals)
    }
    wqtl       <- c(wqtl,chromosome[[i]]$wqtl)
    wmrk       <- c(wmrk,chromosome[[i]]$wmrk)
    wneu       <- c(wneu,chromosome[[i]]$wneu)
    nkept      <- c(nkept,chromosome[[i]]$nkept)
    dkept[[i]] <- chromosome[[i]]$dkept
    qtlFreq    <- c(qtlFreq,chromosome[[i]]$qtlFreq)
    mrkFreq    <- c(mrkFreq,chromosome[[i]]$mrkFreq)
    
    LOCqtl <- rbind(LOCqtl,chromosome[[i]]$locqtl)
    LOCmrk <- rbind(LOCmrk,chromosome[[i]]$locmrk)
    
    
  }
  eqtl <- 0
  eeqtl <- 0
  for (i in 1:nChr){
    eeqtl1 <- sum(eeqtl,nkept[i])
    eeqtl  <- eeqtl1
    eqtl   <- c(eqtl,eeqtl)
  }
  
  add.term_qtl <- c(rep(eqtl[-(nChr+1)],each=nqtl))
  Wqtl <- Wqtl + add.term_qtl
  
  add.term_mrk <- rep(eqtl[-(nChr+1)],c(Nmrk))
  Wmrk <- Wmrk + add.term_mrk
  
  
  cqtl <- c(matrix(wcol(Wqtl),byrow=T,nrow=2))              ## This funcation can be re-written
  cmrk <- c(matrix(wcol(Wmrk),byrow=T,nrow=2))
  
  
  colnames(Animals)[cmrk] <- NameMrk(nmrk*nChr)             ## give the colnames for markers
  colnames(Animals)[cqtl] <- Nameqtl(nqtl*nChr)             ## give the colnames for qtl
  
  
  
  each <- eqtl*2
  cut.it <- list()
  for (i in 2:(nChr+1)){
    res <- (each[i-1]+1):each[i]
    cut.it[[i-1]] <- res
  }                                                                 
  
  sort <- rank(runif(nIndPop[ibreed],0,1))
  BASE <- data.frame(breed=ibreed,Gen=0,ID=1:nIndPop[ibreed],sex=sample(c(1,2),size=nIndPop[ibreed],
                                                                        prob=c(NSires[ibreed]/(NSires[ibreed]+NDams[ibreed]),
                                                                               NDams[ibreed]/(NSires[ibreed]+NDams[ibreed])),replace=TRUE),#sex= c(1,2), 
                     Animals)#[sort <=NAnim,])
  
  BreedDat[[ibreed]] <- list(chromosome=chromosome, 
                             BASE=BASE, 
                             cut.it=cut.it, 
                             #qtlEffects=aMatrix,
                             cqtl=cqtl,
                             cmrk=cmrk,
                             Wqtl=Wqtl,
                             Wmrk=Wmrk) 
  
  
  cat(paste("breed ",ibreed," is done \n"))
}  ## loop over breeds


#-------------------------------------------------------------------------------------------------------------------
# sample additive effects  
#-------------------------------------------------------------------------------------------------------------------
baseHapQTLBreeds <- NULL  ##QTL haplotypes
for (ibreed in 1:nBreeds){
  baseHapQTLBreeds <- rbind(baseHapQTLBreeds,BreedDat[[ibreed]]$BASE[,cqtl+4]);lala(baseHapQTLBreeds)
}

hap <- baseHapQTLBreeds  ##transfer QTL haplotypes to hap
QTL <- abs(hap[,is.eve(1:ncol(hap))]+hap[,is.odd(1:ncol(hap))]-4)  ##QTL gentoypes

set.seed(1:5)  
if (nTraits==1) {  ##If one trait, d was not simulated due to time. Can add later 
  ## additive
  add <- rowSums(t(t(QTL)*a1))
  TempVar <- var(add) ;TempVar
  ## scaling
  scale <- sqrt(genomicMatrix/TempVar)
  sa1   <- a1*scale
  saMatrix  <- matrix(sa1,nrow=length(sa1))
  baseAdd  <- rowSums(t(t(QTL)*saMatrix[,itrait]))
  var(baseAdd)
  
} else { ## If multiple traits
  #-------------------------------------------------------------------------------------------------------------------
  # sample additive effects  (before scaling)
  #-------------------------------------------------------------------------------------------------------------------
  aMatrix <- mvrnorm(n = nqtl*nChr, rep(0,nTraits), 
                     Sigma=genomicMatrix, 
                     empirical = FALSE);var(aMatrix[,1]);var(aMatrix[,2]);var(aMatrix)
  
  baseAdd <- matrix(nrow=nrow(QTL),ncol=nTraits)
  for (itrait in 1:nTraits) {
    add <- rowSums(t(t(QTL)*aMatrix[,itrait]))
    baseAdd[,itrait] <- add
    TempVar <- var(add) ;TempVar
  }
  
  #-------------------------------------------------------------------------------------------------------------------
  # sample additive effects  (after scaling)   - saMatrix, baseAdd
  #-------------------------------------------------------------------------------------------------------------------  
  #saMatrix =read.table("saMatrix")
  V_beforeScaling <- solve(var(baseAdd))
  S <- genomicMatrix%*%V_beforeScaling; S
  S <- sqrtm(S);S
  saMatrix <- t(S%*%t(aMatrix))
  
  #overwritten by new baseAdd after scaling
  baseAdd <- matrix(nrow=nrow(QTL),ncol=nTraits)
  for (itrait in 1:nTraits) {
    add <- rowSums(t(t(QTL)*saMatrix[,itrait]))
    baseAdd[,itrait] <- add
    TempVar <- var(add) ;TempVar
  }
  ## test the the variance and covariance for additive effects after scaling
  round(var(baseAdd),2)
  
  #-------------------------------------------------------------------------------------------------------------------
  # sample dominance effects (before scaling)   
  #-------------------------------------------------------------------------------------------------------------------  
  dMatrix <- aMatrix
  #  overDominance <- rep(NA,nTraits)   ## remember to report this result 
  for (itrait in 1:nTraits) {
    domDegree <- rnorm(nrow(aMatrix),mean=dMean[itrait],sd=dSd[itrait])
    dMatrix[,itrait] <- domDegree*abs(aMatrix[,itrait])
    #    overDominance[itrait] <- length(which(abs(dMatrix[,itrait])>abs(aMatrix[,itrait])))
  }
  
  QTLDom <- (abs(abs(QTL-1)-1))
  
  baseDom <- matrix(nrow=nrow(QTL),ncol=nTraits)
  for (itrait in 1:nTraits) {
    tempD <- rowSums(t(t(QTLDom)*dMatrix[,itrait]))
    baseDom[,itrait] <- tempD
    
  }
  TempVar <- var(baseDom) ;TempVar
  #-------------------------------------------------------------------------------------------------------------------
  # sample dominance effects (after scaling)   - sdMatrix, baseDom
  #-------------------------------------------------------------------------------------------------------------------  
  #sdMatrix=read.table("sdMatrix")
  V_beforeScaling <- solve(var(baseDom))
  SDom <- domMatrix%*%V_beforeScaling; SDom
  SDom <- sqrtm(SDom);SDom
  sdMatrix <- t(SDom%*%t(dMatrix))
  
  #overwritten by new baseDom after scaling
  baseDom <- matrix(nrow=nrow(QTL),ncol=nTraits)
  for (itrait in 1:nTraits) {
    tempD <- rowSums(t(t(QTLDom)*sdMatrix[,itrait]))
    baseDom[,itrait] <- tempD
    
  }  
  
  round(var(baseDom),2)
  #-------------------------------------------------------------------------------------------------------------------
  # sample residuals for each trait
  #-------------------------------------------------------------------------------------------------------------------  
  baseRes <- baseDom 
  for (itrait in 1:nTraits) {
    baseRes[,itrait]  <- rnorm(n=nrow(QTL),mean=0,sd=sqrt(diag(resMatrix)[itrait]))
  }
  
} ## loop for number of traits 

basePhen <- baseAdd+baseDom+baseRes
nIndPop0 <- c(0,nIndPop)
for (ibreed in 1:nBreeds){
  TakeRow <- (nIndPop0[ibreed]+1):(nIndPop0[ibreed]+nIndPop0[ibreed+1])
  BreedDat[[ibreed]]$BASE <- data.frame(BreedDat[[ibreed]]$BASE[,1:4],baseAdd[TakeRow,],baseDom[TakeRow,],
                                        baseRes[TakeRow,],basePhen[TakeRow,],BreedDat[[ibreed]]$BASE[,-c(1:4)])  
}

NameIt <- function(x,y="add"){
  return(paste0(y, rep(1:x)))
}

for (ibreed in 1:nBreeds){
  pointer <- 5
  colnames(BreedDat[[ibreed]]$BASE)[pointer:(pointer-1+nTraits)] <- NameIt(nTraits,"add")
  pointer <- pointer+nTraits
  colnames(BreedDat[[ibreed]]$BASE)[pointer:(pointer-1+nTraits)] <- NameIt(nTraits,"dom")
  pointer <- pointer+nTraits
  colnames(BreedDat[[ibreed]]$BASE)[pointer:(pointer-1+nTraits)] <- NameIt(nTraits,"res")
  pointer <- pointer+nTraits
  colnames(BreedDat[[ibreed]]$BASE)[pointer:(pointer-1+nTraits)] <- NameIt(nTraits,"phen")
}

BreedDat[[ibreed]]$BASE[1:6,1:34]

# cor(BreedDat[[1]]$BASE$add1,BreedDat[[1]]$BASE$add2)
# cor(BreedDat[[2]]$BASE$add1,BreedDat[[2]]$BASE$add2)
# cor(BreedDat[[3]]$BASE$add1,BreedDat[[3]]$BASE$add2)
# 
# 
# 
# cor(BreedDat[[1]]$BASE$add1,BreedDat[[1]]$BASE$add2)
# cor(BreedDat[[2]]$BASE$add1,BreedDat[[2]]$BASE$add2)
# cor(BreedDat[[3]]$BASE$add1,BreedDat[[3]]$BASE$add2)
write.table(saMatrix,"saMatrix",quote=F)
write.table(sdMatrix,"sdMatrix",quote=F)
write.table(aMatrix,"aMatrix",quote=F)
write.table(cqtl,"cqtl",quote=F)
write.table(cmrk,"cmrk",quote=F)
write.table(Wmrk,"Wmrk",quote=F)
write.table(Wqtl,"Wqtl",quote=F)
write.table(cut.it,"cut.it",quote=F)
write.table(domDegree,"domDegree",quote=F)
#write.table(BreedDat,"BreedDat",quote=F)
#selection


SelMeth    <- "pesudoGS"
MatingMeth <- "RANDOM"


#source(paste0(DIR,"input.r"))        ## The input 
#source(paste0(DIR,"genome.r"))       ## The genome

##---------------------------------------------------------------------------------------
## Base population
##---------------------------------------------------------------------------------------
Dat.list <- list()
rr <- 1
for(rr in 1){
  #ibreed <- 1
  #for (ibreed in 3:nBreeds) {
  ibreed=2
  newnewpedigree=NULL
  newsnp.code=NULL
  newQTL=NULL
  newpheno=NULL
  id=NULL
  b <- ibreed
  NAnim <- nIndPop[ibreed]
  BASE <- BreedDat[[ibreed]]$BASE   ## BASE is to store haplotypes for each breed
  BASE.tag <- temp <- BASE[,-c(1:(4+4*nTraits))];lala(temp) ## To track IBD
  BASE.tag[,c(is.odd(1:ncol(temp)))] <-  as.numeric(rownames(BASE))
  BASE.tag[,c(is.eve(1:ncol(temp)))] <- -as.numeric(rownames(BASE))
  BASE.tag <- data.frame(BASE[,c(1:(4+4*nTraits))], BASE.tag); lala(BASE.tag)
  pedigree <- data.frame(BASE[,c(1:(4+4*nTraits))], SireID=0, DamID=0);lala(pedigree) ##only the pedigree and phenotypes
  Dat.list [[ibreed]] <- list(Dat=BASE, Dat.tag=BASE.tag, pedigree=pedigree)  ## put haplotypes, IBD data and pedigree in to a list for all breeds
  #str(Dat.list [[ibreed]])
  
  
  DAT <- list() ## To store all the data including haplotypes, pedigree, inbreeding and accuracy
  GG  <- list() ## To store the pedigree from generation 0 to the current generation for each breed (need to change)
  #g <- 2
  for (g in (0:nGen)){
    if(g==0) {
      Dat     <- Dat.list[[ibreed]]$Dat ;head(Dat)                          ## haplotype data
      Dat.tag <- Dat.list[[ibreed]]$Dat.tag ;head(Dat.tag)                  ## inbreeding data
      pedigree <- data.frame(BASE[,c(1:4)],SireID=0,DamID=0)                ##,tbv=stbv,res=res0,phen=phen0)    
    }
    ##----------------------------------------------------------------------------------------------
    ## calculate true breeding values and phenotypes
    cqtl    <- BreedDat[[ibreed]]$cqtl                         ## columns for qtl
    cmrk    <- BreedDat[[ibreed]]$cmrk
    if (g==0) {
      hap    <- (Dat[,-c(1:(4+4*nTraits))])[,cqtl] 
    } else {
      
      hap    <- (Dat[,-c(1:4)])[,cqtl] 
    }
     ## qtl data
    QTL    <- abs(hap[,is.eve(1:ncol(hap))]+hap[,is.odd(1:ncol(hap))]-4) ;dim(QTL) ## genotypes: 0 1 2 
    QTLDom <- (abs(abs(QTL-1)-1)) ## genotypes 0 1 0 
    Add <- matrix(nrow=nrow(QTL),ncol=nTraits)
    for (itrait in 1:nTraits) {
        add <- rowSums(t(t(QTL)*saMatrix[,itrait]))
        Add[,itrait] <- add
      }
      sum(genomicMatrix-var(Add));#should be around 0!
      
      Add <- matrix(nrow=nrow(QTL),ncol=nTraits)
      for (itrait in 1:nTraits) {
        add <- rowSums(t(t(QTL)*saMatrix[,itrait]))
        Add[,itrait] <- add
      }
      sum(genomicMatrix-var(Add));#should be around 0!
      
      Dom <- matrix(nrow=nrow(QTL),ncol=nTraits)
      for (itrait in 1:nTraits) {
        dom <- rowSums(t(t(QTLDom)*sdMatrix[,itrait]))
        Dom[,itrait] <- dom
      }
      sum(domMatrix-var(Dom));#should be around 0!
      
      Res <- Dom 
      for (itrait in 1:nTraits) {
        Res[,itrait]  <- rnorm(n=nrow(QTL),mean=0,sd=sqrt(diag(resMatrix)[itrait]))
      }
      sum(resMatrix-var(Res))
      
      Phen <- Add+Dom+Res
    
    
    if (g==0) {
      GGx <- data.frame(BASE[,c(1:4)],SireID=0,DamID=0,Add,Dom,Res,Phen)
      colnames(GGx)[-c(1:6)] <- c(NameIt(nTraits,"add"),NameIt(nTraits,"dom"),
                                  NameIt(nTraits,"res"),NameIt(nTraits,"phen"))
    } else {
      GGx     <- data.frame(pedigree,Add,Dom,Res,Phen)
      colnames(GGx)[-c(1:6)] <- c(NameIt(nTraits,"add"),NameIt(nTraits,"dom"),
                                  NameIt(nTraits,"res"),NameIt(nTraits,"phen"))
      lala(GGx)
    }
    ## phenotype of each animal
    if(g==0) {
      Ped.all <-  data.frame(pedigree)  ## pedigree data for each breed
    }
    
    ##----------------------------------------------------------------------------------------------
    ## calculate allele frequencies and homozygosities at QTL and markers
    ##----------------------------------------------------------------------------------------------
    if (g==0) {
      pqtl <- getFreq((Dat[,-c(1:(4+4*nTraits))])[,cqtl])
      pmrk <- getFreq((Dat[,-c(1:(4+4*nTraits))])[,cmrk])	
      
      homQtl <- t(getHetero((Dat[,-c(1:(4+4*nTraits))])[,cqtl])) 
      homMrk <- t(getHetero((Dat[,-c(1:(4+4*nTraits))])[,cmrk])) 
    } else {
      pqtl <- getFreq((Dat[,-c(1:4)])[,cqtl]) 
      pmrk <- getFreq((Dat[,-c(1:4)])[,cmrk])	
      
      homQtl <- t(getHetero((Dat[,-c(1:4)])[,cqtl])) 
      homMrk <- t(getHetero((Dat[,-c(1:4)])[,cmrk])) 
      
    }
    
    
    ##----------------------------------------------------------------------------------------------
    ## calculate true inbreeding
    ##----------------------------------------------------------------------------------------------
    if (g==0) {
      F.IBD   <- findIBD.dir(Dat.tag[,-c(1:(4+4*nTraits))])   ## the true inbreeding at each locus
    } else {
      F.IBD   <- findIBD.dir(Dat.tag[,-c(1:4)])
    }
    F.qtl   <- F.IBD[BreedDat[[ibreed]]$Wqtl]   ## the true inbreeding at qtl
    F.mrk   <- F.IBD[BreedDat[[ibreed]]$Wmrk]   ## the true inbreeding at mrk
    setwd(paste0(DIR,"inbreeding"))
    write.table(F.qtl,paste0("Rep",rr,"-pedigree-breed",ibreed,"qtltrueinbreeding",g,".txt"))
    write.table(F.mrk,paste0("Rep",rr,"-pedigree-breed",ibreed,"mrktrueinbreeding",g,".txt"))
    ##----------------------------------------------------------------------------------------------
    ## calculate pedigree inbreeding
    ##----------------------------------------------------------------------------------------------
    GG[[g+1]] <- GGx	
    ppp <- NULL	  
    for (gg in (1:(g+1))){
      ppp <- rbind(ppp,GG[[gg]][,c(3,6,5)]) ; #9 ,2
    }  
    Fped.ind  <- calcInbreeding(ppp)[GGx$ID]
    Fped      <- mean(Fped.ind)
    setwd(paste0(DIR,"inbreeding"))
    write.table(Fped.ind,paste0("Rep",rr,"-pedigree-breed",ibreed,"pedigreeinbreeding",g,".txt"))
    
    ##----------------------------------------------------------------------------------------------
    ## snp data is extracted if you need to do genomic selection
    ##----------------------------------------------------------------------------------------------
    if(g==0) {
      snp.data <- Dat[,-c(1:(4+4*nTraits))] 
    } else {
      snp.data <- Dat[,-c(1:4)] 
    }
    
    snp.uncode <- snp.data[,cmrk]  #(haplotype)
    snp.code   <- codeGenotype(snp.uncode) ; lala(snp.code) ## marker genotypes 
    
    ##----------------------------------------------------------------------------------------------
    ## For simplicity, we select on aggregate additive value
    ##----------------------------------------------------------------------------------------------     
    setwd(DIR)
    newpedigree=cbind(pedigree$ID,pedigree$SireID,pedigree$DamID)
    newnewpedigree=rbind(newnewpedigree,newpedigree)
    turn=matrix(c(1:nrow(newnewpedigree)))
    newnewnewpedigree=cbind(newnewpedigree,turn)
    pheno=cbind(GGx$ID,GGx[,ibreed*2+23],GGx[,ibreed*2+24])
    newpheno=rbind(newpheno,pheno)
    mean=matrix(1,nrow(newpheno))
    fixpheno=cbind(newpheno[,1],mean,newpheno[,2],newpheno[,3])
    write.table(fixpheno,"phenotype.txt",col.names = F,row.names = F,quote = F)
    num <-seq(from=1,to=NAni[ibreed],by=5)
    if(g == 10){
      newsnp.code=snp.code[num,]
      newsnp.code=apply(newsnp.code, 2, as.integer)
      rownames(newsnp.code)=GGx$ID[num]
      newQTL=QTL[num,]
      rownames(newQTL)=c(GGx$ID[num])
    }
    if(g > 10){
      newnewsnp.code=snp.code[num,]
      newnewsnp.code=apply(newnewsnp.code, 2, as.integer)
      rownames(newnewsnp.code)=c(GGx$ID[num])
      newnewQTL=QTL[num,]
      rownames(newnewQTL)=c(GGx$ID[num])
      newsnp.code=rbind(newsnp.code,newnewsnp.code)
      newQTL=rbind(newQTL,newnewQTL)
    }
    if(g < 4){
      AggregateAddTrue <- ecoValue[[ibreed]]%*%t(Add); AggregateAddTrue
      GGx$GEBV <- t(AggregateAddTrue)
    }
    if(g>3&g<10){
      
      write.table(newnewnewpedigree,"pedigree.txt",col.names = F,row.names = F,quote = F)
      run_DMU(
        phe_col_names=c("Id","mean","Trait1","Trait2"), 
        target_trait_name=list(c("Trait1"),(c("Trait2"))),
        fixed_effect_name=list(c("mean"),c("mean")),    
        random_effect_name=list(c("Id"),c("Id")),       
        covariate_effect_name=NULL,                             
        genetic_effect_name="Id",	                 
        phe_path=DIR,                          
        phe_name="phenotype.txt",                   
        integer_n=2,                                 
        analysis_model="PBLUP_A",                    
        dmu_module="dmuai",                          
        relationship_path=DIR,                 
        relationship_name="pedigree.txt",          
        output_result_path=DIR                  
      )
      ebv=read.table("corrected_phe_dmu.txt")
      ebv=ebv[,c(1,2,6)]
      id=GGx$ID
      newebv=ebv[ebv[,1]%in%id,]
      newebv=newebv[,-1]
      newebv=apply(newebv,2,as.numeric)
      AggregateAddTrue <- ecoValue[[ibreed]]%*%t(Add); AggregateAddTrue
      GGx$GEBV <- t((newecoValue[[ibreed]])%*%t(newebv))
    }
    setwd(DIR)
    if(g > 9){
      #newnewnewpedigree=newnewpedigree[(nrow(newnewpedigree)*0.7):nrow(newnewpedigree),]
      write.table(newnewnewpedigree,"pedigree.txt",col.names = F,row.names = F,quote = F)
      #newsnp.code=as.matrix(newsnp.code)
      kinship_result=cal_kinship(
        input_data_numeric=newsnp.code,          #provided hapmap data object
        kinship_type=c("G_A"),    #type of dominance effect
        inbred_type=c("Homozygous"),      #type of inbreeding coefficients
        output_matrix_type="col_three",
        output_matrix_path=DIR,
        output_matrix_name="G_A_col_three.txt",
        return_result=FALSE)
      run_DMU(
        phe_col_names=c("Id","mean","Trait1","Trait2"),
        target_trait_name=list(c("Trait1"),c("Trait2")),                                
        fixed_effect_name=NULL,     #¹Ì¶¨Ð§Ó¦Ãû³Æ
        random_effect_name=list(c("Id"),c("Id")),       
        covariate_effect_name=NULL,                             
        genetic_effect_name="Id",	                  
        phe_path=DIR,                          
        phe_name="phenotype.txt",                    
        integer_n=2,
        analysis_model="SSBLUP_A",                   
        dmu_module="dmu5",                          
        dmu_algorithm_code=2,
        provided_prior_file_path="PBLUP_A_Trait1_Trait2",
        provided_prior_file_name="PAROUT",
        relationship_path=DIR,               
        relationship_name=c("pedigree.txt","G_A_col_three.txt_G_A_col_three.txt"),   
        output_result_path=DIR                  
      )
      
      
      id=matrix(GGx$ID)
      gebv=read.table("corrected_phe_dmu.txt")
      gebv=gebv[,c(1,2,6)]
      id=GGx$ID
      newgebv=gebv[gebv[,1]%in%id,]
      newgebv=newgebv[,-1]
      newgebv=apply(newgebv,2,as.numeric)
      AggregateAddTrue <- ecoValue[[ibreed]]%*%t(Add); AggregateAddTrue
      GGx$GEBV <- t((newecoValue[[ibreed]])%*%t(newgebv))
    }
    Acc <- cor(t(AggregateAddTrue),GGx$GEBV);Acc
    cat(Acc)
    setwd(DIR)
    write.table(GGx,paste0("Rep",rr,"-pedigree-breed",ibreed,"generation",g,".txt"))
    if(ibreed==1){
      if(g==16){
        write.table(newnewnewpedigree,"breed1ped")
        write.table(newsnp.code,"breed110-16geno")
        write.table(newQTL,"breed110-16QTL")
      }
    }
    if(ibreed==2){
      if(g==16){
        write.table(newnewnewpedigree,"breed2ped")
        write.table(newsnp.code,"breed210-16geno")
        write.table(newQTL,"breed210-16QTL")
      }
    }
    if(ibreed==3){
      if(g==16){
        write.table(newnewnewpedigree,"breed3ped")
        write.table(newsnp.code,"breed310-16geno")
        write.table(newQTL,"breed310-16QTL")
      }
    }
    #Acc <- cor(t(AggregateAddTrue),GGx$GEBV);Acc    ## calculate the accuracy; Now the accuracy is 1
    
    #select on GEBV
    Males <- subset(GGx,sex==1)
    Sires <- Males[order(Males$GEBV, decreasing = T), 3][1:NSires[ibreed]];Sires
    Sires <- sample(Sires,length(Sires))
    
    Dams <- subset(GGx,sex==2)
    Dams <- Dams[order(Dams$GEBV, decreasing = T), 3][1:NDams[ibreed]];Dams
    Dams <- sample(Dams,length(Dams))      
    
    DamID <- rep(Dams, each = LS[ibreed])
    length(unique(DamID))
    DS    <- floor(NDams[ibreed] / NSires[ibreed])
    SireID <- rep(Sires, each = DS * LS[ibreed])
    
    
    extra <- (NAni[ibreed]-length(SireID))/LS[ibreed]
    if (extra > 0) {
      rest <- sample(SireID,size=extra)
      SireID <- c(SireID,rep(rest,each=LS[ibreed]))
    }
    
    pedigree  <- data.frame(breed=ibreed,Gen=g+1,ID=c((NAnim+g*NAni[ibreed]+1):(NAnim+NAni[ibreed]+g*NAni[ibreed])),
                            sex=sample(c(1,2),size=NAni[ibreed],prob=c(0.5,0.5),replace=TRUE),SireID=c(SireID),
                            DamID)
    
    
    parents <- as.numeric(t(cbind(pedigree[, 5], pedigree[, 6])))
    ##----------------------------------------------------------------------------------------------
    
    ## adjust the id of the parents  (check later)      
    if (g==0){
      parents.adj <- parents
    } else{
      subB1 <- NAnim+(g-1)*sum(NAni[b])
      parents.adj <- parents-subB1
    }
    Dat.next      <- NULL
    Dat.tag.next  <- NULL
    ##----------------------------------------------------------------------------------------------
    ## sampling the offspring    
    for (wchr in 1:nChr) {
      if (g==0) {
        parhap       <- displayHap(as.matrix((Dat[,-c(1:(4+4*nTraits))])[,cut.it[[wchr]]][parents.adj,]))  ## there is a problem with cut.it 
        parhap.tag   <- displayHap(as.matrix(Dat.tag[,-c(1:(4+4*nTraits))] [,cut.it[[wchr]]][parents.adj,]))
      } else {
        parhap       <- displayHap(as.matrix((Dat[,-c(1:(4))])[,cut.it[[wchr]]][parents.adj,]))  ## there is a problem with cut.it 
        parhap.tag   <- displayHap(as.matrix(Dat.tag[,-c(1:(4))] [,cut.it[[wchr]]][parents.adj,]))
      }
      
      nrec    <- rpois((2*NAni[b]),1)
      nkept   <- nmrk+nqtl 
      recpos  <- chopos(nrec, nkept, (2*NAni[b]))
      outt    <- recom2((2*NAni[b]),recpos,parhap,parhap.tag)		
      temp.untag    <- outt[[2]]
      temp.tag      <- outt[[1]]
      Dat.next     <- cbind(Dat.next,temp.untag)
      Dat.tag.next <- cbind(Dat.tag.next,temp.tag)
      cat(paste("chr ",wchr," is done for breed ",b, "\n"))
    }
    Dat.next <- data.frame(Dat.next)
    colnames(Dat.next)[c(cqtl,cmrk)] <- c(Nameqtl(length(cqtl)/2),NameMrk(length(cmrk)/2))
    colnames(Dat.tag.next) <- colnames(Dat.next)
    
    setwd(DIR)
    write.table(Dat.next,paste0("Rep",rr,"-Haplotypes-breed",ibreed,"generation",g+1,".txt"))
    ##----------------------------------------------------------------------------------------------
    ## update the haplotypes
    Dat     <- data.frame(pedigree[,1:4],Dat.next)                      ## haplotype data
    Dat.tag <- data.frame(pedigree[,1:4],Dat.tag.next)             ## inbreeding data
    ##----------------------------------------------------------------------------------------------
    
    cat(paste("Breed",b,"- Generation",g,"is done \n"))
    
    
    cat(paste("print results of additive effects \n"))
    if (ibreed==1) {
      cat(paste("meanAdd1 - Generation",g,"is ",mean(GGx$add1) ," \n"))
      cat(paste("varAdd1 - Generation",g,"is ",var(GGx$add1) ," \n"))
    }
    
    if (ibreed==2) {
      cat(paste("meanAdd3 - Generation",g,"is ",mean(GGx$add3) ," \n"))
      cat(paste("varAdd3 - Generation",g,"is ",var(GGx$add3) ," \n"))
    }
    
    if (ibreed==3) {
      cat(paste("meanAdd5 - Generation",g,"is ",mean(GGx$add5) ," \n"))
      cat(paste("varAdd15 - Generation",g,"is ",var(GGx$add5) ," \n"))
    }
    
    
    cat(paste("print results of dominance effects \n"))
    
    if (ibreed==1) {
      cat(paste("meanDom1 - Generation",g,"is ",mean(GGx$dom1) ," \n"))
      cat(paste("varDom1 - Generation",g,"is ",var(GGx$dom1) ," \n"))
    }
    
    if (ibreed==2) {
      cat(paste("meanDom3 - Generation",g,"is ",mean(GGx$dom3) ," \n"))
      cat(paste("varDom3 - Generation",g,"is ",var(GGx$dom3) ," \n"))
    }
    
    if (ibreed==3) {
      cat(paste("meanDom5 - Generation",g,"is ",mean(GGx$dom5) ," \n"))
      cat(paste("varDom5 - Generation",g,"is ",var(GGx$dom5) ," \n"))
    }
    
    cat(paste("-------generation",g,"is done ------------------------- \n"))  
    ## Some output
  } #g
  cat(paste("-------breed",ibreed,"is done ------------------------- \n"))
  #} #breed  
  cat(paste("-------rep",rr,"is done ------------------------- \n"))
}#rep
write.table(saMatrix,"saMatrix",quote=F)
write.table(sdMatrix,"sdMatrix",quote=F)
write.table(aMatrix,"aMatrix",quote=F)
write.table(cqtl,"cqtl",quote=F)
write.table(cmrk,"cmrk",quote=F)
write.table(Wmrk,"Wmrk",quote=F)
write.table(Wqtl,"Wqtl",quote=F)
write.table(cut.it,"cut.it",quote=F)
write.table(domDegree,"domDegree",quote=F)
write.table(BreedDat,"BreedDat",quote=F)
