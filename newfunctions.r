
cat("#######################################################\n")
cat("##                   libraries                         \n")
cat("#######################################################\n")

x <- c("pedigree", "doBy", "rrBLUP", "BLR","functional")
lapply (x,require,character.only=T)

cat("#######################################################\n")
cat("##                   functions                         \n")
cat("#######################################################\n")

##--------------------------------------------------------------
## print the first 10 rows
lala  <- function(data, x=10, y=10) {
  if(is.vector(data)==TRUE){
    if (length(data) < 10) {
      return(data)
    } else {
      return(data[1:10])
    } 
  }
  if(class(data)=="matrix"|class(data)=="data.frame") {
    if (nrow(data) <= 10 & ncol(data) <= 10) {
      return(data[1:nrow(data),1:ncol(data)])
    }
    if (nrow(data) <= 10 & ncol(data) >= 10) {
      return(data[1:nrow(data),1:10])
    }
    if (nrow(data) >= 10 & ncol(data) <= 10) {
      return(data[1:10,1:ncol(data)])
    }
    if (nrow(data) >=10 & ncol(data)>=10) {
      return(data[1:10,1:10])
    }
  }
}
##--------------------------------------------------------------
## take the odd or even columns                
is.odd <- function(x) x %% 2 != 0 
is.eve <- function(x) x %% 2 == 0

is.even <- function(x) x %% 2 == 0
##--------------------------------------------------------------
## get the frequency of allele "1"
getFreq <- function(dat) {
  geno <- abs (dat[,is.odd(1:ncol(dat))] + 
                 dat[,is.eve(1:ncol(dat))]-4)                ## haplotypes are coverted to genotypes
  freq <- apply(geno, 2, function(x) sum(x)/(length(x)*2))   ## calculate the frequency of each column of genotypes
  return(freq)
}

getFreq3 <- function(x){
  x <- as.matrix(x)
  if (ncol(x)> 2) {
    qtl1 <- x[,seq(from=1,to=ncol(x)-1,by=2)]
    qtl2 <- x[,seq(from=2,to=ncol(x),by=2)]
    qtl3 <- abs(qtl1+qtl2-4)
    getFreq2 <- function(try){
      res <-(length(try[try==2])*2+length(try[try==1]))/(2*length(try))
      return(res)}
    result <- apply(qtl3,2,FUN=getFreq2)
    return(result)
  }
  if (ncol(x)==2){
    qtl1 <- x[,1]
    qtl2 <- x[,2]
    qtl3 <- abs(qtl1+qtl2-4)
    
    getFreq2 <- function(try){
      res <-(length(try[try==2])*2+length(try[try==1]))/(2*length(try))
      return(res)}
    result <- getFreq2(qtl3)
    return(result)
  }
}

##-------------------------------------------------------------- 
## get minor allele frequencies
getMinor <- function(dat){
  Freq  <- matrix(getFreq(dat),nrow=1)
  Minor <- apply(Freq,2,function(x)min(x,(1-x)))
  return(Minor)  
}

##--------------------------------------------------------------
## mainly used to capture the columns for marker and qtl 
wcol    <- function(x)  c((2*x-1),2*x)            
##--------------------------------------------------------------
## get the true breeding values for each qtl
get.tbv <- function(effdata,a){
  nqtl    <- length(a)
  qtleff  <- matrix(0, nrow = nrow(effdata), ncol = nqtl)
  for (j in 1:(nqtl)) {
    gc <- chooserow(effdata, c(2 * j - 1), 2 * j)[[2]]
    gr <- chooserow(effdata, c(2 * j - 1), 2 * j)[[1]]
    qtleff[gr, j] <- c(rep(a[j], times = gc[1]), 
                       rep(0, times = gc[2]), rep(0, times = gc[3]), 
                       rep(-a[j], times = gc[4]))
  }
  return(qtleff)
}
##--------------------------------------------------------------
## name qtl and Mrk  
NameMrk <- function (x){
  paste("M", "_", rep(1:x,each = 2), 
        "_", rep(1:2, times = x), sep = "")
}
			   
Nameqtl <- function (x){
  paste("qtl", "_", rep(1:x,each = 2), 
        "_", rep(1:2, times = x), sep = "")
}

NameNeu <- function (x){
  paste("neu", "_", rep(1:x,each = 2), 
        "_", rep(1:2, times = x), sep = "")
} 

NameMrkloci <- function(x){
  paste("M", "_", rep(1:x),sep="")
}

Nameqtlloci <- function(x){
  paste("qtl", "_", rep(1:x),sep="")
}

NameNeuloci <- function(x){
  paste("neu", "_", rep(1:x),sep="")
}
##--------------------------------------------------------------
## find IBD
findIBD.dir <- function (x,cal=c("col","row")) {
  cal <- match.arg(cal,choices=c("col","row"),several.ok =FALSE)
  if (cal %in% c("col","row")){
    yy <- x[,is.odd(1:ncol(x))]==x[,is.eve(1:ncol(x))]
    yyy <- yy[,1]
    if (cal=="col"){
      res <- apply(yy,2,FUN=function(yyy) length(which(yyy=="TRUE"))/length(yyy))
    } else {
      res <- apply(yy,1,FUN=function(yyy) length(which(yyy=="TRUE"))/length(yyy))
    }
    return(res)
  } 
}

##--------------------------------------------------------------
codeGenotype <- function(x) {
  col.odd  <- is.odd(1:ncol(x))		
  col.eve  <- is.eve(1:ncol(x))
  snp.code <- 4-(x[,col.odd]+x[,col.eve])
  return(snp.code)
}
##--------------------------------------------------------------
codeGenotype2 <- function(x) {
  col.odd  <- is.odd(1:ncol(x))		
  col.eve  <- is.eve(1:ncol(x))
  snp.code <- (x[,col.odd]+x[,col.eve])-2
  return(snp.code)
}
##--------------------------------------------------------------
## display the haplotypes
displayHap <- function(parhap){ 
    fill    <- matrix(NA,nrow=2*nrow(parhap),ncol=ncol(parhap)/2)
    fill[is.odd(1:nrow(fill)),]  <-  parhap[,is.odd(1:ncol(parhap))]
    fill[is.even(1:nrow(fill)),] <-  parhap[,is.even(1:ncol(parhap))]
    return(fill)
}

##--------------------------------------------------------------
## choose the postions of recombinations		  
chopos   <- function(x,nloci,NAnim){
  uu <- list()
  for (i in 1:NAnim){
    if (x[i] ==0){uu[[i]] <- 0
    }
    else{uu[[i]]<-sample(1:nloci,size=x[i],replace=T)
    }
  }
  return(uu)
}

##--------------------------------------------------------------	
## do recombinations
recom  <- function (N,recpos,parhap1){
  res <- list()
  resu <- c()
  for (i in 1:N){
    nkept <- ncol(parhap1)
    pos <- recpos[[i]]
    pos <- pos[order(pos)]
    
    if (is.odd(nrec[i])==T){
      if (nrec[i]!=1){
        res1 <- NULL
        pos1 <- 0:pos[1]
        for (j in seq(from=2,to=(nrec[i]-1),by=2)){
          pos2 <- c(pos[j]:pos[j+1])
          res1 <- c(res1,pos2)
        }
        resuc<- c(pos1,res1)
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu<- resuc 
        }
        if (samp==2){
          resu<- setdiff(1:nkept,resuc) 
        }
      }
      
      if (nrec[i]==1){
        pos1 <- pos[1]
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu <- c(1:pos1)
        }
        if (samp==2){
          resu <- c(pos1:nkept)
        }
      }
    }
    
    if (is.odd(nrec[i])==F) {
      
      if (nrec[i]!=2 && nrec[i]!=0){
        res1 <- NULL
        pos1 <- 0:pos[1]
        posn <- pos[nrec[i]]:nkept
        for (j in seq(from=2,to=(nrec[i]-2),by=2)){
          pos2 <- c(pos[j]:pos[j+1])
          res1 <- c(res1,pos2)
        }
        resuc <- c(pos1,res1,posn)
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu<- resuc 
        }
        if (samp==2){
          resu<- setdiff(1:nkept,resuc) 
        }
      }
      
      if (nrec[i]==0){
        resu <- 0
      }
      
      if (nrec[i]==2){
        resuc<- c(pos[1]:pos[2])
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu<- resuc 
        }
        if (samp==2){
          resu<- setdiff(1:nkept,resuc) 
        }
      }
    }
    res[[i]] <- resu
  }
  
  Geno1 <- Geno <- parhap1
  for (i in 1:N) {
    Geno[2 * i - 1, res[[i]]] <- Geno1[2 * i, res[[i]]]
    Geno[2 * i, res[[i]]] <- Geno1[2 * i - 1, res[[i]]]
  }
  
  rownames(Geno) <- rep(rownames(parhap1)/2,each=2)
  seg <- c()
  for (j in 1:N){
    seg[j] <- sample((2 * j - 1):(2 * j), 
                     size = 1, prob = c(0.5, 0.5), replace = T)
  }
  Genonext <- Geno[seg,]
  
  #N <- N/2   ##NAni=200
  data3 <- matrix(0, nrow = NAni, ncol = nkept*2)  #dim(data3) <- 2*nloci,NAnim
  
  data3[, is.odd(1:ncol(data3))]  <- Genonext[is.odd(1:nrow(Genonext)),]
  data3[,is.even(1:ncol(data3))]  <- Genonext[is.even(1:nrow(Genonext)),]
  
  return(data3)
}






##-----------------------------------------------------------------------------------------------------------------------
## Long time ago.... needs to work on them again
##-----------------------------------------------------------------------------------------------------------------------

##take cols based on locus ID
##find the cols by index
getq <- function(X,Y){ 
  output <- Y[,X]
  return(output)}


chooserow <- function(data, x, y) {
  row11 <- intersect(which(data[, x] == 1), which(data[, y] == 1))
  row12 <- intersect(which(data[, x] == 1), which(data[, y] == 2))
  row21 <- intersect(which(data[, x] == 2), which(data[, y] == 1))
  row22 <- intersect(which(data[, x] == 2), which(data[, y] == 2))
  rowlala <- c(row11, row12, row21, row22)
  rown <- c(length(row11), length(row12), length(row21), length(row22))
  out <- list(rowlala, rown)
  out
}




Gblup.Huiming <- function(formula,M,lambda,ped){
  data <- ped 
  data$ID <- factor(data$ID)
  depVar <- match(all.vars(formula)[1], colnames(data))
  indVars <- match(all.vars(formula)[-1], colnames(data))
  if (length(indVars) > 0) 
    data <- data[apply(as.matrix(data[, indVars]), 1, function(x) all(!is.na(x))), 
                 ]
  M <- M[match(data$ID, rownames(M)), ]
  
  X <- Matrix(model.matrix(formula, data))
  Y <- Matrix(model.frame(formula, data)[, 1])
  Z <- Matrix(model.matrix(data[, depVar] ~ data$ID))
  ## Ginv <- calcG(M, data = data[, match(all.vars(formula)[1], 
  ##     colnames(data))], solve = TRUE)
  
  
  M <- M[, apply(M, 2, function(x) length(unique(x)) > 1)]
  dd <- apply(M, 2, paste, collapse = "")
  M <- M[, !duplicated(dd)]
  dd <- apply(M, 1, paste, collapse = "")
  dd <- as.numeric(factor(dd, levels = unique(dd)))
  M <- M[unique(dd), ]
  p <- apply(M, 2, mean)/2
  D <- Diagonal(ncol(M), 1/(ncol(M) * 2 * p * (1 - p)))
  M <- Matrix(M - 1)
  G <- as(M %*% D %*% t(M), "dgCMatrix")
  #if (solve) 
  diag(G) <- diag(G)+0.001
  G <- solve(G)
  dd <- diag(ncol(G))[, dd]
  G <- t(dd) %*% G %*% dd
  Ginv <- G
  xtx <- crossprod(X)
  xtz <- crossprod(X, Z)
  ztx <- crossprod(Z, X)
  ztzginv <- crossprod(Z) + lambda * Ginv
  nr <- nrow(xtx) + nrow(ztx)
  LHS <- Matrix(0, ncol = nr, nrow = nr)
  LHS[1:nrow(xtx), 1:ncol(xtx)] <- xtx
  LHS[nrow(xtx) + seq(1:nrow(ztx)), 1:ncol(ztx)] <- ztx
  LHS[1:nrow(xtz), ncol(xtx) + seq(1, ncol(xtz))] <- xtz
  LHS[nrow(xtz) + seq(1:nrow(ztzginv)), ncol(ztx) + seq(1, 
                                                        ncol(ztzginv))] <- ztzginv
  RHS <- Matrix(0, ncol = 1, nrow = nr)
  xty <- crossprod(X, Y)
  RHS[1:nrow(xty), 1] <- xty
  zty <- crossprod(Z, Y)
  RHS[nrow(xty) + seq(1, nrow(zty))] <- zty
  sol <- solve(LHS, RHS)
  return(sol[-1])
}



getIBD.col <- function(data.in){
  data.in <- as.matrix(data.in)
  R <- nrow(data.in)
  C <- ncol(data.in)/2
  IBD <- c()
  for (i in 1:C){
    eIBD <- data.in[,(2*i-1):(2*i)]	   
    findIBD <- sum((table(eIBD)/(2*R))**2)
    IBD     <- c(IBD,findIBD)
  }
  return(round(IBD,digits=4))
}

getIBD.row <- function(data.in){
  data.in <- displayHap(data.in)
  R <- nrow(data.in)/2
  C <- ncol(data.in)
  IBD <- c()
  for (i in 1:R){
    eIBD <- data.in[(2*i-1):(2*i),]	   
    findIBD <- length(which(eIBD[1,]==eIBD[2,]))
    IBD     <- c(IBD,findIBD)
  }
  return(round(IBD,digits=4))
}


getr2 <- function(r2.data1){
  x11      <- length(which(c(r2.data1[,c(1,2)]+
                               r2.data1[,c(3,4)])==2))/(2*nrow(r2.data1))
  p1 <- getFreq3(r2.data1)[1]
  q1 <- getFreq3(r2.data1)[2]	  
  r2 <- ((x11-p1*q1))**2/(p1*(1-p1)*q1*(1-q1))
  return(r2)
}


transAG <- function(snp.matrix){
  
  bases <- c("A","T","C","G")
  nsnps <- ncol(snp.matrix)/2
  snp.base.matrix <- matrix(nrow=nrow(snp.matrix),
                            ncol=ncol(snp.matrix))
  col.index <- matrix(1:(4*nsnps),nrow=nsnps,ncol=2,byrow=TRUE)
  base.index <- matrix(nrow=nsnps,ncol=2)
  
  for (i in 1:nsnps){
    base.index[i,] <- bases[sample(1:4,2,replace=FALSE)]
    replace.index  <-which(snp.matrix[,col.index[i,]]==1)
    replace.index2 <-which(snp.matrix[,col.index[i,]]==2)
    snp.base.matrix[,col.index[i,]] <- 
      replace(snp.matrix[,col.index[i,]],
              replace.index,base.index[i,][1])
    snp.base.matrix[,col.index[i,]] <- 
      replace(snp.base.matrix[,col.index[i,]],
              replace.index2,base.index[i,][2])
  }
  return(snp.base.matrix)	
}	





trans12 <- function(snp.code.tag){
  c.row <- nrow(snp.code.tag) 
  c.col <-ncol(snp.code.tag) 
  test <- matrix(nrow=c.row,ncol=c.col)
  test[,is.odd(1:c.col)] <- 
    snp.code.tag[,is.odd(1:c.col)]/snp.code.tag[,is.even(1:c.col)]
  test[,is.even(1:c.col)] <- 
    snp.code.tag[,is.even(1:c.col)]/snp.code.tag[,is.odd(1:c.col)]
  test1 <- test[,is.odd(1:c.col)]
  test1[test[,is.odd(1:c.col)]!=1] <- 1
  test2 <- test[,is.even(1:c.col)]
  test2[test[,is.even(1:c.col)]!=1] <- 2
  
  test[,is.odd(1:c.col)] <- test1
  test[,is.even(1:c.col)] <- test2
  return(test)
}




#This is to organize the results from ROH
#DK PC
###############################################################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This is does the summary; it's not easy to understand...
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun= function(xx, col, na.rm) {
                   c( N    = length2(xx[,col], na.rm=na.rm),
                      mean = mean   (xx[,col], na.rm=na.rm),
                      sd   = sd     (xx[,col], na.rm=na.rm)
                   )
                 },
                 measurevar,
                 na.rm
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean"=measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}


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

limdig <- function(x) round(x,digits=4)


reform <- function (x){
  reformed <- matrix(0,nrow=2*nrow(x),ncol=ncol(x)/2)
  for (i in 1:(ncol(x)/2)){
    reformed[,i] <- matrix(cbind(x[,(2*i-1)],x[,(2*i)]),ncol=1,byrow=F)
  }
  return(reformed)
}  

findIBD <- function(tryit){
  tryit.r <- reform(tryit)
  R <- nrow(tryit)
  findIBD <- function(x)sum((table(x)/(2*R))**2)
  ress <- apply(tryit.r,2,findIBD)
  return(ress)
}
############################################################################
#### new functions
recom2  <- function (N,recpos,parhap.untag, parhap.tag){
  res    <- list()
  resu   <- c()
  for (i in 1:N){
    nkept <- ncol(parhap.untag)
    pos <- recpos[[i]]
    pos <- pos[order(pos)]
    
    if (is.odd(nrec[i])==T){
      if (nrec[i]!=1){
        res1 <- NULL
        pos1 <- 0:pos[1]
        for (j in seq(from=2,to=(nrec[i]-1),by=2)){
          pos2 <- c(pos[j]:pos[j+1])
          res1 <- c(res1,pos2)
        }
        resuc<- c(pos1,res1)
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu<- resuc 
        }
        if (samp==2){
          resu<- setdiff(1:nkept,resuc) 
        }
      }
      
      if (nrec[i]==1){
        pos1 <- pos[1]
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu <- c(1:pos1)
        }
        if (samp==2){
          resu <- c(pos1:nkept)
        }
      }
    }
    
    if (is.odd(nrec[i])==FALSE) {
      
      if (nrec[i]!=2 && nrec[i]!=0){
        res1 <- NULL
        pos1 <- 0:pos[1]
        posn <- pos[nrec[i]]:nkept
        for (j in seq(from=2,to=(nrec[i]-2),by=2)){
          pos2 <- c(pos[j]:pos[j+1])
          res1 <- c(res1,pos2)
        }
        resuc <- c(pos1,res1,posn)
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu<- resuc 
        }
        if (samp==2){
          resu<- setdiff(1:nkept,resuc) 
        }
      }
      
      if (nrec[i]==0){
        resu <- 0
      }
      
      if (nrec[i]==2){
        resuc<- c(pos[1]:pos[2])
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu<- resuc 
        }
        if (samp==2){
          resu<- setdiff(1:nkept,resuc) 
        }
      }
    }
    res[[i]] <- resu
  }
  
  Geno1 <- Geno.tag <- parhap.tag
  for (i in 1:N) {
    Geno.tag[2 * i - 1, res[[i]]] <- Geno1[2 * i, res[[i]]]
    Geno.tag[2 * i, res[[i]]] <- Geno1[2 * i - 1, res[[i]]]
  }
  rownames(Geno.tag) <- rep(rownames(parhap.tag)/2,each=2)
  
  Geno2 <- Geno.untag <- parhap.untag
  for (i in 1:N) {
    Geno.untag[2 * i - 1, res[[i]]] <- Geno2[2 * i, res[[i]]]
    Geno.untag[2 * i, res[[i]]] <- Geno2[2 * i - 1, res[[i]]]
  }
  rownames(Geno.untag) <- rep(rownames(parhap.untag)/2,each=2)
  
  seg <- c()
  for (j in 1:N){
    seg[j] <- sample((2 * j - 1):(2 * j), 
                     size = 1, prob = c(0.5, 0.5), replace = T)
  }
  Genonext.tag   <- Geno.tag[seg,]
  Genonext.untag <- Geno.untag[seg,]
  
  
  data3.tag <- matrix(0, nrow = NAni[b], ncol = nkept*2)  #dim(data3) <- 2*nloci,NAnim
  
  data3.tag[,is.odd(1:ncol(data3.tag))]  <- Genonext.tag[is.odd(1:nrow(Genonext.tag)),]
  data3.tag[,is.eve(1:ncol(data3.tag))]  <- Genonext.tag[is.eve(1:nrow(Genonext.tag)),]
  
  
  data3.untag <- matrix(0, nrow = NAni[b], ncol = nkept*2)  #dim(data3) <- 2*nloci,NAnim
  
  data3.untag[,is.odd(1:ncol(data3.untag))]<-Genonext.untag[is.odd(1:nrow(Genonext.untag)),]
  data3.untag[,is.eve(1:ncol(data3.untag))]<-Genonext.untag[is.eve(1:nrow(Genonext.untag)),]
  return(list(tag=data3.tag, untag=data3.untag))
}




findIBD.exp <- function(tryit,cal=c("col","row")){
  cal <- match.arg(cal,choices=c("col","row"),several.ok =FALSE)
  if (cal %in% c("col","row")){
    tryit.r <- reform(tryit)
    R <- nrow(tryit)
    C <- ncol(tryit)
    if (cal=="col"){
      findIBD <- function(x)sum((table(x)/(length(x)))**2)
      ress <- apply(tryit.r,2,findIBD)
    } else {
      findIBD <- function(x)sum((table(x)/(length(x)))**2)
      ress <- apply(tryit,1,findIBD)
    }
    return(ress)
  }
}
getFreq <- getFreq3 
########################################################################
## This is the method from VanRaden 2008
#########################################################################
calM_G2 <- function (snp.code,BASE,cmrk) {
  
  M <- snp.code-1
  pm <- getFreq(BASE[,cmrk+3])
  Pm <- 2*(pm-0.5)
  x <- pm*(1-pm)
  P <- matrix(rep(Pm,times=nrow(M)),nrow=nrow(M),byrow=T)
  Z <- M-P
  G <- (Z%*%t(Z))/(2*sum(x))
  FVR2 <-  rd(mean(diag(G)-1))
  G2 <- diag(G)
  return(list(diag(G),FVR2))
}



########################################################################
## This is the method from Yang et al 2010 Nature
########################################################################
calM_G3 <- function(snp.code,BASE,cmrk){
  x <- snp.code
  p <- getFreq(BASE[,cmrk+3])
  p <- matrix(rep(p,times=nrow(x)),nrow=nrow(x),byrow=T)
  G <- rd(((x**2)-(1+2*p)*x+2*p**2)/(2*p*(1-p)))
  FVR3 <- rd(mean(G))  
  return(list(rd(rowMeans(G)),FVR3,rd(G)))
}

########################################################################
## This is the method from VanRaden 2008 (regression)
########################################################################
calM_G <- function (snp.code,pp,NAni) {
  res <- (snp.code-1)%*%(t(snp.code-1))
  pa2 <- (nrow(pp)-NAni+1):nrow(pp)
  pa1  <- setdiff(1:nrow(pp),pa2)
  A <- makeA(ped=pp,which=rep(TRUE,nrow(pp)))
  A <- read.table("A.txt")
  nInd <- length(pa1)+length(pa2)
  Ar <- matrix(0,nrow = nInd,ncol = nInd)
  Ar[as.matrix(A[,1:2])] <- A[,3]
  VA <- Ar[pa2,pa2]
  lala(VA)
  mod = lm(c(res) ~ 1+c(VA))
  g0 <- mod$coefficient[1]
  g1 <- mod$coefficient[2]
  g0
  g1
  M1 <- matrix(rep(1,length(pa2)),ncol=1)
  M11 <- M1%*%t(M1)
  Gmatrix <- (res-g0*(M11))/g1
  FVR1 <- mean(diag(Gmatrix))-1
  return(FVR1)
}
#####################################################################
calM_G1.1 <- function (snp.code,pp,NAni) {
  M <- snp.code-1
  res <- M%*%(t(M))
  pa2 <- (nrow(pp)-NAni+1):nrow(pp)
  pa1  <- setdiff(1:nrow(pp),pa2)
  A <- makeA(ped=pp,which=rep(TRUE,nrow(pp)))
  A <- read.table("A.txt")
  nInd <- length(pa1)+length(pa2)
  Ar <- matrix(0,nrow = nInd,ncol = nInd)
  Ar[as.matrix(A[,1:2])] <- A[,3]
  A <- Ar[pa2,pa2]
  x1 <- sum(A)
  x2 <- sum(A**2)
  x3 <- sum(res)
  x4 <- sum(res*A)
  y <- matrix(c(x3,x4),nrow=2)
  x <- matrix(c(NAni**2,x1,x1,x2),nrow=2)
  pro <- solve(x)%*%y
  g0 <- pro[1]
  g1 <- pro[2]
  M1 <- matrix(rep(1,length(pa2)),ncol=1)
  M11 <- M1%*%t(M1)
  Gmatrix <- (res-g0*(M11))/g1
  FVR1 <- mean(diag(Gmatrix))-1
  return(FVR1)
}


calkin <- function(M,ret=c("rl","rlAB","sm","sim_smin")) {
  ret <- match.arg(ret,choices=c("rl","rlAB","sm","sim_smin"),
                   several.ok=FALSE)		
  if (ret %in% c("sm", "sm-smin")) {
    M <- snp.code
    marker <- M - (max(M, na.rm = TRUE) - 1)
    m <- ncol(marker)
    d <- 1 - (tcrossprod(marker) + m)/(2 * m)
    #if (ret == "sm-smin") {
    s <- 1 - d
    smin <- min(s, na.rm = TRUE)
    f <- (s - smin)/(1 - smin)
    mean(f)
  } else {f <- 1 - d
  kmat <- 2 * f
  U <- kmat
  }	
  return(mean(diag(U)-1))		
}

makeG <- function(M,BASE,cmrk){
  # M <- snp.code
  p <- getFreq(BASE[,cmrk+3])
  NN <- nrow(M)
  PP <- ncol(M)
  Gmatrix <- matrix(NA,nrow=NN,ncol=NN)
  pq2 <- 2*p*(1-p)
  Z1 <- sweep(M,2,2*p)
  Z1 <- (Z1/sqrt(pq2))
  Z1 <- as.matrix(Z1)
  Gmatrix <- (Z1%*%t(Z1))/PP
  col2 <- rep(1:NN ,rep(NN :1))
  RR <- cbind(1:NN,NN)
  col1 <- unlist(apply(RR,1,function(x) return(x[1]:x[2])))
  pedi <- cbind(col1=col1,col2=col2)
  p <- matrix(rep(p,times=nrow(M)),nrow=nrow(M),byrow=T)
  GGG <- rd(((M**2)-(1+2*p)*M+2*p**2)/(2*p*(1-p)))
  diag(Gmatrix) <- 1+rowMeans(GGG)
  G3 <- Gmatrix[lower.tri(Gmatrix, diag = TRUE)]
  N <-0
  if (g>=1) {
    if (g==1) {
      N <- 0
    } else {
      N <- pedigree[1, 2]-NAni-1
    }
  }
  pedi <- pedi+N 
  finalG <- cbind(pedi,G3)
  write.table(finalG,file=paste0("Gmatrix.txt"),row.names=FALSE,col.names=FALSE)
  return(TRUE)
}



makeGdis <- function(M,BASE,cmrk,wm,pedigree){
  # M <- snp.code
#  M <- (codeGenotype(snp.uncode.all))[,wm]
  p <- getFreq(BASE[,cmrk+4])
  p <- p[wm]
  NN <- nrow(M)
  PP <- ncol(M)
  Gmatrix <- matrix(NA,nrow=NN,ncol=NN)
  pq2 <- 2*p*(1-p)
  Z1 <- sweep(M,2,2*p)
  Z1 <- (Z1/sqrt(pq2))
  Z1 <- as.matrix(Z1)
  Gmatrix <- (Z1%*%t(Z1))/PP
  col2 <- rep(1:NN ,rep(NN :1))
  RR <- cbind(1:NN,NN)
  col1 <- unlist(apply(RR,1,function(x) return(x[1]:x[2])))
  pedi <- cbind(col1=col1,col2=col2)
  p <- matrix(rep(p,times=nrow(M)),nrow=nrow(M),byrow=T)
  GGG <- rd(((M**2)-(1+2*p)*M+2*p**2)/(2*p*(1-p)))
  diag(Gmatrix) <- 1+rowMeans(GGG)
  G3 <- Gmatrix[lower.tri(Gmatrix, diag = TRUE)]
  N <- nrow(pedigree)
  FULL <-data.frame(pedi,G3)
  #start <- which(pedi[,2]==nrow(snp.uncode.uncode)-nrow(pedigree)+1)[1]
  start <- which(pedi[,2]==nrow(snp.uncode)-nrow(pedigree)+1)[1]
  end <- nrow(pedi)
  FinalG <- FULL[start:end,]
  N <-0
  write.table(FinalG,file=paste0("Gmatrix.txt"),row.names=FALSE,col.names=FALSE)
  return(TRUE)
}




makeGdis <- function(M,BASE,cmrk,wm,pedigree){
  # M <- snp.code
  #  M <- (codeGenotype(snp.uncode.all))[,wm]
  p <- getFreq(BASE[,cmrk+4])
  p <- p[wm]
  NN <- nrow(M)
  PP <- ncol(M)
  Gmatrix <- matrix(NA,nrow=NN,ncol=NN)
  pq2 <- 2*p*(1-p)
  Z1 <- sweep(M,2,2*p)
  Z1 <- (Z1/sqrt(pq2))
  Z1 <- as.matrix(Z1)
  Gmatrix <- (Z1%*%t(Z1))/PP
  col2 <- rep(1:NN ,rep(NN :1))
  RR <- cbind(1:NN,NN)
  col1 <- unlist(apply(RR,1,function(x) return(x[1]:x[2])))
  pedi <- cbind(col1=col1,col2=col2)
  p <- matrix(rep(p,times=nrow(M)),nrow=nrow(M),byrow=T)
  GGG <- rd(((M**2)-(1+2*p)*M+2*p**2)/(2*p*(1-p)))
  diag(Gmatrix) <- 1+rowMeans(GGG)
  G3 <- Gmatrix[lower.tri(Gmatrix, diag = TRUE)]
  N <-0
  if (g>=1) {
    if (g==1) {
      N <- NAnim
    } else {
      N <- pedigree[1, 3]-1#-NAni-1
    }
  }
  pedi <- pedi+N 
  finalG <- cbind(pedi,G3)
  write.table(finalG,file=paste0("Gmatrix.txt"),row.names=FALSE,col.names=FALSE)
  return(TRUE)
}




rd <- function(x){
  round(x,digits=4)
}

getHetero <- function(x) {
  x <- as.matrix(x)
  if (ncol(x) > 2) {
    qtl1 <- x[, seq(from = 1,
                    to = ncol(x) - 1,
                    by = 2)]
    qtl2 <- x[, seq(from = 2,
                    to = ncol(x),
                    by = 2)]
    qtl3 <- abs(qtl1 + qtl2 - 4)
    getH <- function(try) {
      res <-(length(try[try==1]))/(length(try))
      return(res)
    }
    result <- apply(qtl3, 2, FUN = getH)
    return(result)
  }
  if (ncol(x) == 2) {
    qtl1 <- x[, 1]
    qtl2 <- x[, 2]
    qtl3 <- abs(qtl1 + qtl2 - 4)
    
    getH <- function(try) {
      res <-(length(try[try==1]))/(length(try))
      return(res)
    }
    result <- getH(qtl3)
    return(result)
  }
}

makeCorrelatedColumn <- function(tbv,rho=0.6) {
  n     <- length(tbv)          # length of vector
  #rho   <- 0.6                   # desired correlation = cos(angle)
  theta <- acos(rho)             # corresponding angle
  x1    <- tbv        # fixed given data
  x2    <- rnorm(n, 2, 0.5)      # new random data
  X     <- cbind(x1, x2)         # matrix
  Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)
  
  Id   <- diag(n)                               # identity matrix
  Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
  P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
  x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
  Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
  Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
  
  x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
  cor(x1, x)  
  return(x)
  
  
}


getEff <- function(QTL,a1){

  AA <- 2*a1
  Aa <- a1
  aa <- rep(0,times=length(a1))
  nqtl    <- ncol(QTL)
  nind    <- nrow(QTL)
  qtlBV  <- matrix(0, nrow = nrow(QTL), ncol = ncol(QTL))
  for (j in 1:nqtl) {
    for (i in 1:nind){
      if(QTL[i,j]==2) {qtlBV[i,j] <- AA[j]}
      if(QTL[i,j]==1) {qtlBV[i,j] <- Aa[j]}
      if(QTL[i,j]==0) {qtlBV[i,j] <- aa[j]}
    }
  }
  
  return(qtlBV)
}
getEff2 <- function(QTL,alpha1,alpha2,qtlFreq){
  qq <- 1-qtlFreq
  
  AA <- 2*alpha1
  Aa <- alpha1+alpha2
  aa <- 2*alpha2
  nqtl    <- ncol(QTL)
  nind    <- nrow(QTL)
  qtlBV  <- matrix(0, nrow = nrow(QTL), ncol = ncol(QTL))
  for (j in 1:nqtl) {
    for (i in 1:nind){
      if(QTL[i,j]==2) {qtlBV[i,j] <- AA[j]}
      if(QTL[i,j]==1) {qtlBV[i,j] <- Aa[j]}
      if(QTL[i,j]==0) {qtlBV[i,j] <- aa[j]}
    }
  }
  
  return(qtlBV)
}

recom2c  <- function (N,recpos,parhap.untag, parhap.tag){
  res    <- list()
  resu   <- c()
  for (i in 1:N){
    nkept <- ncol(parhap.untag)
    pos <- recpos[[i]]
    pos <- pos[order(pos)]
    
    if (is.odd(nrec[i])==T){
      if (nrec[i]!=1){
        res1 <- NULL
        pos1 <- 0:pos[1]
        for (j in seq(from=2,to=(nrec[i]-1),by=2)){
          pos2 <- c(pos[j]:pos[j+1])
          res1 <- c(res1,pos2)
        }
        resuc<- c(pos1,res1)
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu<- resuc 
        }
        if (samp==2){
          resu<- setdiff(1:nkept,resuc) 
        }
      }
      
      if (nrec[i]==1){
        pos1 <- pos[1]
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu <- c(1:pos1)
        }
        if (samp==2){
          resu <- c(pos1:nkept)
        }
      }
    }
    
    if (is.odd(nrec[i])==FALSE) {
      
      if (nrec[i]!=2 && nrec[i]!=0){
        res1 <- NULL
        pos1 <- 0:pos[1]
        posn <- pos[nrec[i]]:nkept
        for (j in seq(from=2,to=(nrec[i]-2),by=2)){
          pos2 <- c(pos[j]:pos[j+1])
          res1 <- c(res1,pos2)
        }
        resuc <- c(pos1,res1,posn)
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu<- resuc 
        }
        if (samp==2){
          resu<- setdiff(1:nkept,resuc) 
        }
      }
      
      if (nrec[i]==0){
        resu <- 0
      }
      
      if (nrec[i]==2){
        resuc<- c(pos[1]:pos[2])
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu<- resuc 
        }
        if (samp==2){
          resu<- setdiff(1:nkept,resuc) 
        }
      }
    }
    res[[i]] <- resu
  }
  
  Geno1 <- Geno.tag <- parhap.tag
  for (i in 1:N) {
    Geno.tag[2 * i - 1, res[[i]]] <- Geno1[2 * i, res[[i]]]
    Geno.tag[2 * i, res[[i]]] <- Geno1[2 * i - 1, res[[i]]]
  }
  rownames(Geno.tag) <- rep(rownames(parhap.tag)/2,each=2)
  
  Geno2 <- Geno.untag <- parhap.untag
  for (i in 1:N) {
    Geno.untag[2 * i - 1, res[[i]]] <- Geno2[2 * i, res[[i]]]
    Geno.untag[2 * i, res[[i]]] <- Geno2[2 * i - 1, res[[i]]]
  }
  rownames(Geno.untag) <- rep(rownames(parhap.untag)/2,each=2)
  
  seg <- c()
  for (j in 1:N){
    seg[j] <- sample((2 * j - 1):(2 * j), 
                     size = 1, prob = c(0.5, 0.5), replace = T)
  }
  Genonext.tag   <- Geno.tag[seg,]
  Genonext.untag <- Geno.untag[seg,]

  data3.tag <- matrix(0, nrow = NAnic[b], ncol = nkept*2)  #dim(data3) <- 2*nloci,NAnim
  
  data3.tag[,is.odd(1:ncol(data3.tag))]  <- Genonext.tag[is.odd(1:nrow(Genonext.tag)),]
  data3.tag[,is.eve(1:ncol(data3.tag))]  <- Genonext.tag[is.eve(1:nrow(Genonext.tag)),]
  
  
  data3.untag <- matrix(0, nrow = NAnic[b], ncol = nkept*2)  #dim(data3) <- 2*nloci,NAnim
  
  data3.untag[,is.odd(1:ncol(data3.untag))]<-Genonext.untag[is.odd(1:nrow(Genonext.untag)),]
  data3.untag[,is.eve(1:ncol(data3.untag))]<-Genonext.untag[is.eve(1:nrow(Genonext.untag)),]
  
  return(list(tag=Genonext.tag, untag=Genonext.untag,tag2=data3.tag, untag2=data3.untag))
}

recom2cc  <- function (N,recpos,parhap.untag, parhap.tag){
  res    <- list()
  resu   <- c()
  for (i in 1:N){
    nkept <- ncol(parhap.untag)
    pos <- recpos[[i]]
    pos <- pos[order(pos)]
    
    if (is.odd(nrec[i])==T){
      if (nrec[i]!=1){
        res1 <- NULL
        pos1 <- 0:pos[1]
        for (j in seq(from=2,to=(nrec[i]-1),by=2)){
          pos2 <- c(pos[j]:pos[j+1])
          res1 <- c(res1,pos2)
        }
        resuc<- c(pos1,res1)
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu<- resuc 
        }
        if (samp==2){
          resu<- setdiff(1:nkept,resuc) 
        }
      }
      
      if (nrec[i]==1){
        pos1 <- pos[1]
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu <- c(1:pos1)
        }
        if (samp==2){
          resu <- c(pos1:nkept)
        }
      }
    }
    
    if (is.odd(nrec[i])==FALSE) {
      
      if (nrec[i]!=2 && nrec[i]!=0){
        res1 <- NULL
        pos1 <- 0:pos[1]
        posn <- pos[nrec[i]]:nkept
        for (j in seq(from=2,to=(nrec[i]-2),by=2)){
          pos2 <- c(pos[j]:pos[j+1])
          res1 <- c(res1,pos2)
        }
        resuc <- c(pos1,res1,posn)
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu<- resuc 
        }
        if (samp==2){
          resu<- setdiff(1:nkept,resuc) 
        }
      }
      
      if (nrec[i]==0){
        resu <- 0
      }
      
      if (nrec[i]==2){
        resuc<- c(pos[1]:pos[2])
        samp <- sample(c(1,2),1)
        if (samp==1){
          resu<- resuc 
        }
        if (samp==2){
          resu<- setdiff(1:nkept,resuc) 
        }
      }
    }
    res[[i]] <- resu
  }
  
  Geno1 <- Geno.tag <- parhap.tag
  for (i in 1:N) {
    Geno.tag[2 * i - 1, res[[i]]] <- Geno1[2 * i, res[[i]]]
    Geno.tag[2 * i, res[[i]]] <- Geno1[2 * i - 1, res[[i]]]
  }
  rownames(Geno.tag) <- rep(rownames(parhap.tag)/2,each=2)
  
  Geno2 <- Geno.untag <- parhap.untag
  for (i in 1:N) {
    Geno.untag[2 * i - 1, res[[i]]] <- Geno2[2 * i, res[[i]]]
    Geno.untag[2 * i, res[[i]]] <- Geno2[2 * i - 1, res[[i]]]
  }
  rownames(Geno.untag) <- rep(rownames(parhap.untag)/2,each=2)
  
  seg <- c()
  for (j in 1:N){
    seg[j] <- sample((2 * j - 1):(2 * j), 
                     size = 1, prob = c(0.5, 0.5), replace = T)
  }
  Genonext.tag   <- Geno.tag[seg,]
  Genonext.untag <- Geno.untag[seg,]
  
  data3.tag <- matrix(0, nrow = NAnicc[b], ncol = nkept*2)  #dim(data3) <- 2*nloci,NAnim
  
  data3.tag[,is.odd(1:ncol(data3.tag))]  <- Genonext.tag[is.odd(1:nrow(Genonext.tag)),]
  data3.tag[,is.eve(1:ncol(data3.tag))]  <- Genonext.tag[is.eve(1:nrow(Genonext.tag)),]
  
  
  data3.untag <- matrix(0, nrow = NAnicc[b], ncol = nkept*2)  #dim(data3) <- 2*nloci,NAnim
  
  data3.untag[,is.odd(1:ncol(data3.untag))]<-Genonext.untag[is.odd(1:nrow(Genonext.untag)),]
  data3.untag[,is.eve(1:ncol(data3.untag))]<-Genonext.untag[is.eve(1:nrow(Genonext.untag)),]
  return(list(tag=Genonext.tag, untag=Genonext.untag,tag2=data3.tag, untag2=data3.untag))
}

cat(paste("All the functions have been read :) \n",sep=""))

