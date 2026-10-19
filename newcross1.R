DIR <- ""
library(purrr)
library(plyr)
library(expm)
library(MASS)
library(hiblup)
library(data.table)
source(paste0(DIR,"functions.r"))
setwd(DIR)

nChr    <- 18
chrl    <- 120

nqtl     <- 200            ## The number of QTL per chromosome
nloci    <- 4092         ## The number of loci per chromosome(marker+QTL)
nmrk    <- nloci-nqtl

a=read.table('Rep1-pedigree-breed3generation10.txt')
bb=fread('Rep1-Haplotypes-breed3generation10.txt',data.table=F)
c=read.table('Rep1-pedigree-breed2generation10.txt')
dd=fread('Rep1-Haplotypes-breed2generation10.txt',data.table=F)

nRep     <- 2 
nBreeds  <- Npop <- 3
nIndPop  <- c(4500, 4500, 4500)  ## The number of elements must equal nBreeds(base population)
nSires0  <- c(2250, 2250, 2250)
nDams0   <- c(2250, 2250, 2250)
NSires   <- c(60,45,30)  ## The number of selected sires; If it is truncation selection(next generation)
NDams    <- c(600,450,300)  ## The number of selected dams;  If it is truncation selection
NAni     <- c(7200,4500,2400) ## The total number of animals in the population (offspring)
#aMatrix  <- read.table("/home/txiang/jlli/newnewa/newaMatrix")
saMatrix     <-    read.table("saMatrix")
cut.it   <- read.table("cut.it")
Wqtl   <- read.table("Wqtl")
cqtl   <- read.table("cqtl")
Wqtl=as.matrix(Wqtl)
Wqtl=as.vector(Wqtl)
cqtl=as.matrix(cqtl)
cqtl=as.vector(cqtl)
saMatrix=cbind(saMatrix[,5],saMatrix[,6],saMatrix[,3],saMatrix[,4],saMatrix[,1],saMatrix[,2])

LSc=c(3,3)
ibreed=4
g=1
NAnic=c(1500,3000)
NDamsc=c(500,1000)
NSiresc=c(50,100)

#variance covariance matrix of dominance effects
nTraits=6
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
##计算纯种育种值


ebv1=cbind(a$ID,a$sex,a$GEBV)
ebv2=cbind(c$ID,c$sex,c$GEBV)
#ebv3=cbind(e$ID,e$sex,e$GEBV)
## parents of purebreed and crossbreed from breed1
ibreed=1
Males1 <- subset(ebv1,ebv1[,2]==1)
Sires1 <- Males1[order(Males1[,3], decreasing = T),1 ][1:NSires[ibreed]] ##breed1纯种
c1sires <- Males1[order(Males1[,3], decreasing = T),1 ][1:(NSires[ibreed]+NSiresc[ibreed])]##杂种中breed1父本,假设20个
Dams1 <- subset(ebv1,ebv1[,2]==2)
females1 <- Dams1[order(Dams1[,3], decreasing = T), 1][1:NDams[ibreed]]
c2dams <- Dams1[order(Dams1[,3], decreasing = T), 1][1:(NDams[ibreed]+NDamsc[ibreed+1])]##杂种中breed1母本,假设40个
## parents of purebreed and crossbreed from breed2
ibreed=2
Males2 <- subset(ebv2,ebv2[,2]==1)
Sires2 <- Males2[order(Males2[,3], decreasing = T),1 ][1:NSires[ibreed]] ##breed2纯种
c2sires2 <- Males2[order(Males2[,3], decreasing = T),1 ][1:(NSires[ibreed]+NSiresc[ibreed])]##杂种中breed2父本,假设20个
Dams2 <- subset(ebv2,ebv2[,2]==2)
females2 <- Dams2[order(Dams2[,3], decreasing = T), 1][1:NDams[ibreed]]
c1dams2 <- Dams2[order(Dams2[,3], decreasing = T), 1][1:(NDams[ibreed]+NDamsc[ibreed-1])]##杂种中breed2母本,假设40个
## parents of purebreed and crossbreed from breed3


##----------------------------------------------------------------------------------------------
##crossbreed pedigree
LSc=c(3,3)
ibreed=4
g=1
NAnic=c(1500,3000)
NDamsc=c(500,1000)
NSiresc=c(50,100)
c1dams2=sample(c1dams2,NDamsc[1])
c2dams=sample(c2dams,NDamsc[2])
cDamID1 <- rep(c1dams2, each = 3)
cDamID2 <- rep(c2dams, each = 3)
length(unique(cDamID1))
DS    <- floor(NDamsc[ibreed-3] / NSiresc[ibreed-3])
c1sires=sample(c1sires,NSiresc[1])
c2sires2=sample(c2sires2,NSiresc[2])
cSireID1 <- rep(c1sires, each = DS*3)
cSireID2 <- rep(c2sires2, each = DS*3)
sex1 <- matrix(1,NAnic[ibreed-3],1)
sex1[sample(1:NAnic[ibreed-3],0.5*NAnic[ibreed-3]),] <- 2
sex2 <- matrix(1,NAnic[ibreed-2],1)
sex2[sample(1:NAnic[ibreed-2],0.5*NAnic[ibreed-2]),] <- 2
cross1pedigree  <- data.frame(breed=ibreed,Gen=g,ID=c(((g-1)*NAnic[ibreed-3]+1):((g-1)*NAnic[ibreed-3]+NAnic[ibreed-3])),
                              sex=sex1,Sirebreed=1,dambreed=2,SireID=c(cSireID1),
                              DamID=c(cDamID1))
cross2pedigree  <- data.frame(breed=ibreed,Gen=g,ID=c((g*NAnic[ibreed-3]+1):(g*NAnic[ibreed-3]+NAnic[ibreed-2])),
                              sex=sex2,Sirebreed=2,dambreed=1,SireID=c(cSireID2),
                              DamID=c(cDamID2))
cparents1 <- cbind(cross1pedigree[, 7],cross1pedigree[, 8])
cparents2 <- cbind(cross2pedigree[, 7],cross2pedigree[, 8])
g=10
NAnim <- 4500
subB1 <- NAnim+(g-1)*sum(NAni[1])
subB2 <- NAnim+(g-1)*sum(NAni[2])
newcparents1 <- cbind(as.matrix(cparents1[,1]-subB1),as.matrix(cparents1[,2]-subB2))
newcparents2 <- cbind(as.matrix(cparents2[,1]-subB2),as.matrix(cparents2[,2]-subB1))
c1parents.adj <- as.numeric(t(cbind(newcparents1[, 1], newcparents1[, 2])))
c2parents.adj <- as.numeric(t(cbind(newcparents2[, 1], newcparents2[, 2])))
bBASE.tag1 <- bb[1:3600,]
bBASE.tag2 <- bb[3601:7200,]
bBASE.tag1[,c(is.odd(1:ncol(bb)))] <-  as.numeric(rownames(bBASE.tag1))
bBASE.tag2[,c(is.odd(1:ncol(bb)))] <-  as.numeric(rownames(bBASE.tag2))
bBASE.tag1[,c(is.eve(1:ncol(bb)))] <- -as.numeric(rownames(bBASE.tag1))
bBASE.tag2[,c(is.eve(1:ncol(bb)))] <- -as.numeric(rownames(bBASE.tag2))
bBASE.tag=rbind(bBASE.tag1,bBASE.tag2)
dBASE.tag <- dd
dBASE.tag[,c(is.odd(1:ncol(dd)))] <-  as.numeric(rownames(dd))
dBASE.tag[,c(is.eve(1:ncol(dd)))] <- -as.numeric(rownames(dd))

#c1sire.adj <- c1sires-subB1
#c1dams2.adj <- c1dams2-subB1
#c2sire2.adj <- c2sires2-subB1
#c2dams.adj <- c2dams-subB1
#c1sireshap <- b[c(c1sire.adj),]
#c2damshap <- b[c(c2dams.adj),]
#c2sireshap <- d[c(c2sire2.adj),]
#c1damshap <- d[c(c1dams2.adj),]
#data3.tag <- matrix(0, nrow = NAnic[b], ncol = nkept*2)
#newc1sireshap <- displayHap(c1sireshap)
#newc2damshap <- displayHap(c2damshap)
#newc2sireshap <- displayHap(c2sireshap)
#newc1damshap <- displayHap(c1damshap)
##----------------------------------------------------------------------------------------------
##sire from breed1 and dam from breed2 
b=1
data1=matrix(0,length(c1parents.adj),ncol(bb))
data1.tag=matrix(0,length(c1parents.adj),ncol(bBASE.tag))
bb=as.matrix(bb)
dd=as.matrix(dd)
bBASE.tag=as.matrix(bBASE.tag)
dBASE.tag=as.matrix(dBASE.tag)
data1[seq(1,nrow(data1),2),]=bb[c1parents.adj[seq(1,length(c1parents.adj),2)],]
data1[seq(0,nrow(data1),2),]=dd[c1parents.adj[seq(0,length(c1parents.adj),2)],]
data1.tag[seq(1,nrow(data1.tag),2),]=bBASE.tag[c1parents.adj[seq(1,length(c1parents.adj),2)],]
data1.tag[seq(0,nrow(data1.tag),2),]=dBASE.tag[c1parents.adj[seq(0,length(c1parents.adj),2)],]

cDat.next      <- NULL
cDat.tag.next  <- NULL
cDat.next2      <- NULL
cDat.tag.next2  <- NULL
for (wchr in 1:nChr) {
  cparhap       <- displayHap(as.matrix(data1[,cut.it[[wchr]]]))  ## there is a problem with cut.it 
  cparhap.tag   <- displayHap(as.matrix(data1.tag[,cut.it[[wchr]]]))
  
  
  nrec    <- rpois((2*NAnic[b]),1)
  nkept   <- nmrk+nqtl 
  recpos  <- chopos(nrec, nkept, (2*NAnic[b]))
  outt    <- recom2c((2*NAnic[b]),recpos,cparhap,cparhap.tag)		
  ctemp.untag    <- outt[[2]]
  ctemp.tag      <- outt[[1]]
  ctemp.untag2    <- outt[[4]]
  ctemp.tag2      <- outt[[3]]
  cDat.next     <- cbind(cDat.next,ctemp.untag)
  cDat.tag.next <- cbind(cDat.tag.next,ctemp.tag)
  cDat.next2     <- cbind(cDat.next2,ctemp.untag2)
  cDat.tag.next2 <- cbind(cDat.tag.next2,ctemp.tag2)
  cat(paste("chr ",wchr," is done for breed ",b, "\n"))
}
cDat.next <- data.frame(cDat.next)
cDat.next2 <- data.frame(cDat.next2)
write.table(cDat.next2,"Rep1-Haplotypes-c1generation1.txt",col.names = F,row.names = F,quote = F)
#colnames(cDat.next)[c(Wqtl,Wmrk)] <- c(Nameqtl(length(Wqtl)/2),NameMrk(length(Wmrk)/2))
#colnames(cDat.tag.next) <- colnames(cDat.next)
##----------------------------------------------------------------------------------------------
##sire from breed2 and dam from breed1
data2=matrix(0,length(c2parents.adj),ncol(bb))
data2.tag=matrix(0,length(c2parents.adj),ncol(bBASE.tag))

data2[seq(1,nrow(data2),2),]=dd[c2parents.adj[seq(1,length(c2parents.adj),2)],]
data2[seq(0,nrow(data2),2),]=bb[c2parents.adj[seq(0,length(c2parents.adj),2)],]
data2.tag[seq(1,nrow(data2.tag),2),]=dBASE.tag[c2parents.adj[seq(1,length(c2parents.adj),2)],]
data2.tag[seq(0,nrow(data2.tag),2),]=bBASE.tag[c2parents.adj[seq(0,length(c2parents.adj),2)],]

c2Dat.next      <- NULL
c2Dat.tag.next  <- NULL
c2Dat.next2      <- NULL
c2Dat.tag.next2  <- NULL
b=2
for (wchr in 1:nChr) {
  c2parhap       <- displayHap(as.matrix(data2[,cut.it[[wchr]]]))  ## there is a problem with cut.it 
  c2parhap.tag   <- displayHap(as.matrix(data2.tag[,cut.it[[wchr]]]))
  
  
  nrec    <- rpois((2*NAnic[b]),1) 
  nkept   <- nmrk+nqtl 
  recpos  <- chopos(nrec, nkept, (2*NAnic[b]))
  outt    <- recom2c((2*NAnic[b]),recpos,c2parhap,c2parhap.tag)		
  c2temp.untag    <- outt[[2]]
  c2temp.tag      <- outt[[1]]
  c2temp.untag2    <- outt[[4]]
  c2temp.tag2      <- outt[[3]]
  c2Dat.next     <- cbind(c2Dat.next,c2temp.untag)
  c2Dat.tag.next <- cbind(c2Dat.tag.next,c2temp.tag)
  c2Dat.next2     <- cbind(c2Dat.next2,c2temp.untag2)
  c2Dat.tag.next2 <- cbind(c2Dat.tag.next2,c2temp.tag2)
  cat(paste("chr ",wchr," is done for breed ",b, "\n"))
}
c2Dat.next <- data.frame(c2Dat.next)
c2Dat.next2 <- data.frame(c2Dat.next2)
write.table(c2Dat.next2,"Rep1-Haplotypes-c2generation1.txt",col.names = F,row.names = F,quote = F)


qtlhap1=cDat.next[,Wqtl]  ##haplotype of qtl
qtlhap2=c2Dat.next[,Wqtl]
qtlhap1.code=codeGenotype(cDat.next2[,cqtl])
qtlhap2.code=codeGenotype(c2Dat.next2[,cqtl])

newqtlhap1=abs(qtlhap1-2) ##transfer 1 to 1 ,2 to 0
newqtlhap2=abs(qtlhap2-2)
##----------------------------------------------------------------------------------------------
##additive effect of crossbreed
##additive effect of crossbreed
hap1.1=newqtlhap1[is.odd(1:nrow(newqtlhap1)),]
hap1.2=newqtlhap1[is.eve(1:nrow(newqtlhap1)),]
hap2.1=newqtlhap2[is.odd(1:nrow(newqtlhap2)),]
hap2.2=newqtlhap2[is.eve(1:nrow(newqtlhap2)),]
hap1.1=as.matrix(hap1.1)
hap1.2=as.matrix(hap1.2)
hap2.1=as.matrix(hap2.1)
hap2.2=as.matrix(hap2.2)
##trait1
t1addhap1.1=hap1.1%*%(as.matrix(saMatrix[,1]))
t1addhap1.2=hap1.2%*%(as.matrix(saMatrix[,3]))
c1trait1add=t1addhap1.1+t1addhap1.2
t1addhap2.1=hap2.1%*%(as.matrix(saMatrix[,3]))
t1addhap2.2=hap2.2%*%(as.matrix(saMatrix[,1]))
c2trait1add=t1addhap2.1+t1addhap2.2
##trait2
t2addhap1.1=hap1.1%*%(as.matrix(saMatrix[,2]))
t2addhap1.2=hap1.2%*%(as.matrix(saMatrix[,4]))
c1trait2add=t2addhap1.1+t2addhap1.2
t2addhap2.1=hap2.1%*%(as.matrix(saMatrix[,4]))
t2addhap2.2=hap2.2%*%(as.matrix(saMatrix[,2]))
c2trait2add=t2addhap2.1+t2addhap2.2

c1traitadd=matrix(0,nrow(c1trait1add),2)
c1traitadd[,1]=c1trait1add
c1traitadd[,2]=c1trait2add
c2traitadd=matrix(0,nrow(c2trait1add),2)
c2traitadd[,1]=c2trait1add
c2traitadd[,2]=c2trait2add


##----------------------------------------------------------------------------------------------
##dominance of C(AB) before scaling
nTraitsc1=2
QTLc1 <- abs(qtlhap1[is.eve(1:nrow(qtlhap1)),]+qtlhap1[is.odd(1:nrow(qtlhap1)),]-4)  ##QTL
QTLc2 <- abs(qtlhap2[is.eve(1:nrow(qtlhap2)),]+qtlhap2[is.odd(1:nrow(qtlhap2)),]-4)  ##QTL
QTLDomc1 <- (abs(abs(QTLc1-1)-1))
QTLDomc2 <- (abs(abs(QTLc2-1)-1))

sacMatrix=cbind(as.matrix((saMatrix[,1]+saMatrix[,3])/2),as.matrix((saMatrix[,2]+saMatrix[,4])/2))
sdcMatrix <- sacMatrix
for (itrait in 1:nTraitsc1) {
  domcDegree <- rnorm(nrow(sacMatrix),mean=0.193,sd=0.312)
  sdcMatrix[,itrait] <- domcDegree*abs(sacMatrix[,itrait])
}

write.table(sdcMatrix,"sdcMatrix",col.names=F,row.names=F,quote=F)
write.table(domcDegree,"domcDegree",col.names=F,row.names=F,quote=F)
#sdcMatrix=read.table("sdcMatrix")
#sdc3Matrix=read.table("sdc3Matrix")
##----------------------------------------------------------------------------------------------
##dominance of C(AB)
Domc1 <- matrix(nrow=nrow(QTLc1),ncol=nTraitsc1)
for (itrait in 1:nTraitsc1) {
  domc1 <- rowSums(t(t(QTLDomc1)*sdcMatrix[,itrait]))
  Domc1[,itrait] <- domc1
}


Domc2 <- matrix(nrow=nrow(QTLc2),ncol=nTraitsc1)
for (itrait in 1:nTraitsc1) {
  domc2 <- rowSums(t(t(QTLDomc2)*sdcMatrix[,itrait]))
  Domc2[,itrait] <- domc2
}
##res of C(AB) 

resMatrixc=matrix(c(1.4,0,0,4.7),2,2)
Resc1 <- Domc1 
for (itrait in 1:nTraitsc1) {
  Resc1[,itrait]  <- rnorm(n=nrow(QTLc1),mean=0,sd=sqrt(diag(resMatrixc)[itrait]))
}
sum(resMatrixc-var(Resc1))

Resc2 <- Domc2 
for (itrait in 1:nTraitsc1) {
  Resc2[,itrait]  <- rnorm(n=nrow(QTLc2),mean=0,sd=sqrt(diag(resMatrixc)[itrait]))
}
sum(resMatrixc-var(Resc2))

Phenc1 <- c1traitadd+Domc1+Resc1
Phenc2 <- c2traitadd+Domc2+Resc2
newphec1=cbind(c1traitadd,Domc1,Resc1,Phenc1)
newphec2=cbind(c2traitadd,Domc2,Resc2,Phenc2)
colnames(newphec1)[1:8]=c("add1","add2","dom1","dom2","res1","res2","phe1","phe2")
colnames(newphec2)[1:8]=c("add1","add2","dom1","dom2","res1","res2","phe1","phe2")
ped1=cbind(cross1pedigree,newphec1)
ped2=cbind(cross2pedigree,newphec2)
write.table(ped1,"Rep1-pedigree-c1generation1.txt",col.names = T,row.names = F,quote = F)
write.table(ped2,"Rep1-pedigree-c2generation1.txt",col.names = T,row.names = F,quote = F)


