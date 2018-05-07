## functions required to run the classifier

library(preprocessCore)
library(randomForest)
library(gsubfn)

ReadFromCSV=function(DirPath, mapdir, idxsort=NULL){
  DirNames=dir(DirPath, "*.csv")
  dir.create(file.path(DirPath, "RdataParts"))
  Stains=sapply(c("CDX2", "HTR2B", "ZEB1", "FRMD6", "ker"), 
                function(x) as.character(strapply(DirNames, x, ignore.case=T)))
  Stains=sapply(1:nrow(Stains), function(x) Stains[x, which(Stains[x, ]!="NULL")])
  Stains=toupper(Stains)
  
  temp=list()
  for (i in 1:length(DirNames)){
    # print(sprintf("%g", i))
    tempA=read.table(sprintf("%s/%s", DirPath, DirNames[i]), sep="\t", header=T, stringsAsFactors =F)
    tempA=tempA[ ,-7]
    tt=as.character(strapply(tempA$ImageNo, "[R|_][0-9]+[C|_][0-9]+_[_|A-Z]"))
    tt=strapply(tt, "[0-9]+")
    t2=sapply(1:length(tt), function(x) sprintf("R%sC%s", tt[[x]][1], tt[[x]][2]))
    tempA$ImageNo=t2
    tempA$Slide=ifelse(regexpr("[0-9]{6}", DirNames[i])==-1, i, unlist(strapply(DirNames[i], "[0-9]{6}")))
    colnames(tempA)=c("ImageNo", paste(Stains[i], colnames(tempA)[-1], sep="_"))
    temp[[i]]=tempA
    names(temp)[[i]]=unlist(tempA[1,7])
  }
  MFiles=dir(mapdir, "*.csv")
  xslide=grep("*lide_", MFiles)
  Slide_No=read.csv(sprintf("%s/%s", mapdir, MFiles[xslide]), header=T, stringsAsFactors =F)
  if (length(MFiles)==nrow(Slide_No)+1){
    ## FOr Cairo and AMC mappings
    for (i in 1:nrow(Slide_No)){
      # print(sprintf("%g", i))
      idx=match(as.character(Slide_No[i, 4:8]), names(temp))
      if (length(na.omit(idx))!=5){
        next 
      }
      tempB=merge(temp[[idx[1]]], temp[[idx[2]]], by.x="ImageNo", by.y="ImageNo", sort=T)
      for (j in 3:5){
        tempB=merge(tempB, temp[[idx[j]]], by.x="ImageNo", by.y="ImageNo")
      }
      tempB$Row=strapply(tempB$ImageNo, "[R][0-9]+")
      tempB$Col=strapply(tempB$ImageNo, "C[0-9]+")
      tempB$Row=substr(tempB$Row, 2, 3)   
      tempB$Col=substr(tempB$Col, 2, 3)
      idx2=grep(Slide_No[ i,1], MFiles)
      mapping=read.csv(sprintf('%s/%s', mapdir, MFiles[idx2]), header=F)

      Maps=sapply(1:nrow(tempB), function(x) as.character(mapping[as.numeric(tempB$Row[x]),as.numeric(tempB$Col[x])]))
      tempB$Pat=Maps
      save(tempB, file=sprintf("%s/RdataParts/data_%g.RData", DirPath, i))
    }
  } else {
    ## Leiden Mappings
    for (i in 1:nrow(Slide_No)){
      # merge the lists with the same mapping template
      idx=match(as.character(Slide_No[i, 4:8]), names(temp))
      temp2=list()
      for (j in 1:5){
        tempB=temp[[idx[j]]]
        tempB$Row=strapply(tempB$ImageNo, "R[0-9]+")
        tempB$Row=substr(tempB$Row, 2, 3)
        tempB$Col=strapply(tempB$ImageNo, "C[0-9]+")
        tempB$Col=substr(tempB$Col, 2, 3)
        idx2=grep(sprintf("%s_%s",Slide_No[ i,1], names(temp)[[idx[j]]]), MFiles)
        mapping=read.csv(sprintf('%s%s', mapdir, MFiles[idx2]), header=F)
        Maps=as.character(sapply(1:nrow(tempB), function(x) mapping[as.numeric(tempB$Row[x]),as.numeric(tempB$Col[x])]))
        tempB$Pat=Maps
        temp2[[j]]=tempB      
      }
      names(temp2)=c("CDX2","HTR2B", "FRMD6", "Panker", "ZEB1")
      Out=MergeLeiden(temp2, Slide_No[ i,1])
      tempB=Out$data
      Pattemp=strsplit(tempB$Pat, " ")
      tempB$Pat=sapply(1:nrow(tempB), function(x) Pattemp[[x]][1])
      tempB$Suffix=sapply(1:nrow(tempB), function(x) Pattemp[[x]][length(Pattemp[[x]])])
      ExcSamp=Out$excl
      save(tempB, ExcSamp, file=sprintf("%sRdataParts/data_%g.RData", DirPath, i))
    }
  }
}


MergeLeiden=function(Stains, Slide){
  ## functions
  MergeFun=function(S1, S2, x1, ...){
    merge(S1[S1$Pat %in% names(x1), ], 
          S2[S2$Pat %in% names(x1), ],
          by.x="Pat", by.y="Pat",...)}
  WrapFn=function(Stain, name){
    xa=sapply(1:length(Stain), function(x) which(Stain[[x]]$Pat==name))
    if (class(xa)=="list"){
      idx=sapply(1:length(Stain), function(x) cbind(as.numeric(Stain[[x]]$Row[xa[[x]]]), 
                                                    as.numeric(Stain[[x]]$Col[xa[[x]]])))
      t1=idx[[1]]  
      for (i in 2:length(idx)){
        t1=rbind(t1, idx[[i]])
      }
      t2=t1[duplicated(t1), ]
      t4=unique(t1)
      
      xb=list()
      for (i in 1:nrow(t4)){
        xb[[i ]]=sapply(1:5, function(x) xa[[x]][which(apply(idx[[x]], 1, identical, t4[i, ]))])
      }  
      t3=sapply(xb, function(x) length(unlist(x)))
      sv=which(t3==5)
      if (length(sv)>0){
        xa=xb[[sv]]
        idx=matrix(0, nrow=2, ncol=2)
      } else {
        idx=matrix(c(1,2,3,4), nrow=2, ncol=2)
      }
    }
    else if (class(xa)=="matrix"){
      idx=sapply(1:length(Stain), function(x) cbind(as.numeric(Stain[[x]]$Row[xa[ ,x]]), 
                                                    as.numeric(Stain[[x]]$Col[xa[ ,x]])))    
    }
    if (sum(idx-idx[ ,1])==0){
      if (class(xa)=="matrix"){
        A1=cbind(Stain[[1]][xa[ ,1], ], Stain[[2]][xa[,2], ], Stain[[3]][xa[,3], ],
                 Stain[[4]][xa[,4], ],Stain[[5]][xa[,5], ])}
      else {
        A1=cbind(Stain[[1]][xa[1], ], Stain[[2]][xa[2], ], Stain[[3]][xa[3], ],
                 Stain[[4]][xa[4], ],Stain[[5]][xa[5], ])  
      }
      A1=A1[ ,-grep("Pat", colnames(A1))]
      A1=cbind(paste(name, 1:nrow(A1), sep="_"), A1)
      A1
    }}
  
  ## estimate number of times a name appears
  PatNames=unlist(sapply(1:5, function(x) Stains[[x]]$Pat))
  tt=table(PatNames)
  # merge samples which are present in all stains
  x1=which(tt==5)  
  A1=MergeFun(Stains[[1]], Stains[[2]], x1)
  B1=MergeFun(Stains[[3]], Stains[[4]], x1)
  C1=MergeFun(A1, B1, x1)
  D1=MergeFun(C1, Stains[[5]], x1)
  # merge sample with more than 2 core
  x2=names(tt[which(tt>5&tt<15)])
  x2=x2[grep("EZ", x2)]
  if (length(x2)>0){
    E1=lapply(1:length(x2), function(x) WrapFn(Stains, x2[x]))
    D2=E1[[1]]
    
    if (length(E1)>1){
      for (i in 2:length(E1)){
        D2=rbind(D2, E1[[i]])
      }}
    
    if (!is.null(D2)){
      colnames(D2)=colnames(D1)
      D1=rbind(D1, D2)
    }}
  # merge together
  ## Change the Info in the matrix
  D1$Sample=paste(D1$Pat, Slide, sep="_")
  rownames(D1)=D1$Sample
  D1=D1[ ,-c(grep("Row.", colnames(D1)), grep("Col.", colnames(D1)), 
             grep("ImageNo.", colnames(D1)))]
  ExclSamp=names(tt[round(tt/5)*5!=tt])
  Out=list(data=D1, excl=ExclSamp)
  Out
}



LoadFiles=function(DirPath, rename=F, exclsamps=F){
  pathID=file.path(DirPath, "RdataParts")
  DirNames=dir(pathID, "*RData")
  load(sprintf("%s/%s", pathID, DirNames[1]))  
  CombinedMat=tempB
  NormFeature=function(CombinedMat, feature, normfeature){    
    t1=grep(feature, colnames(CombinedMat))
    t2=grep(paste("KER", feature, sep="_"), colnames(CombinedMat))
    tt=sapply(t1, function(x) CombinedMat[ ,x]/CombinedMat[ ,t2])
    # set an upper threshold of 5?
    tt[tt>2.5]=2.5
    StNames=unlist(strsplit(colnames(CombinedMat)[t1], "_"))[seq(1, 10, 2)]
    StNames=paste(StNames, paste(feature, "norm", sep="."), sep="_")
    colnames(tt)=StNames
    tt}
  if (exclsamps==F & length(DirNames)>1) {
    for (i in 2:length(DirNames)){
      load(sprintf("%s/%s", pathID, DirNames[i]))
      CombinedMat=rbind(CombinedMat, tempB)  
    }
    if (rename==T){
      temp=paste("Sample",CombinedMat$ImageNo,"_", CombinedMat$Slide, sep="_")
      CombinedMat$ImageNo=temp
    }   
  }
  else if (exclsamps==T & length(DirNames)>1) {
    exclSamps=ExcSamp
    for (i in 2:length(DirNames)){
      load(sprintf("%s/%s", pathID, DirNames[i]))
      CombinedMat=rbind(CombinedMat, tempB)
      exclSamps=c(exclSamps, ExcSamp)
    }
  }
  t1=grep("StAreaFrac", colnames(CombinedMat))
  t2=grep("StainInt", colnames(CombinedMat))    
  NewFeat=sapply(1:5, function(x) CombinedMat[ ,t1[x]]*CombinedMat[ , t2[x]])
  StNames=unlist(strsplit(colnames(CombinedMat)[t1], "_"))[seq(1, 10, 2)]    
  StNames=paste(StNames, "Brown.total", sep="_")
  colnames(NewFeat)=StNames
  CombinedMat=cbind(CombinedMat, NewFeat)
  # Add features
  features=c("StAreaFrac", "Brown.total")
  NormFeat=lapply(1:2, function(x) NormFeature(CombinedMat, features[x]))    
  CombinedMat=cbind(CombinedMat, NormFeat[[1]], NormFeat[[2]])
  
  if (exclsamps==T){
    ret=list(data=CombinedMat, exc= exclSamps)
    ret
  }else {
    return(CombinedMat)}
  
}

PreProcess=function(data, classes, ID='Train', rm.range=c(0.2, 0.6), na.rm=T){
  # Determine unique number of patients
  rmID=c("x", "\\.", "Lever", "Milt", "Colon", "\\?")
  x2=unlist(sapply(rmID, function(x) grep(x, data$Pat) ))
  if (length(x2)>0){
    data=data[-x2, ]
  } 
  print(sprintf("Number of Patients with clinical data= %g", length(unique(classes$ID))))
  print(sprintf('Number of MSI patients= %g', length(which(classes$MSI=="MSI"))))
  msi1=which(data$Pat %in% toupper(classes[classes$MSI=="MSI","ID"]))
  d2=data[msi1, ]
  msiNA=which(data$Pat %in% toupper(classes[is.na(classes$MSI),"ID"]))
  print(sprintf('Number of patients with unknown MSI = %g', length(which(is.na(classes$MSI)))))
  print(sprintf('Number of patients to be classified CCS1/3= %g', length(which(classes$MSI=="MSS"))))
  mssPat=which(data$Pat %in% toupper(classes[classes$MSI=="MSS","ID"]))
  data=data[mssPat, ]
  print(sprintf('Number of patients with TMA data = %g', length(unique(data$Pat))))
  print(sprintf('Number of cores belonging to these patients = %g', length(mssPat)))
  # check samples with image fraction outside the given range
  x1=which(data[ ,grep("ImFrac", colnames(data))]>rm.range[2] | data[ ,grep("ImFrac", colnames(data))]<rm.range[1], arr.ind=T)
  print(sprintf('Number available cores after removing poor quality samples = %g', nrow(data[-(x1[ ,1]), ]) ))
  # tempA=length(unique(data$Pat))-length(unique(data$Pat[-x1]))
  print(sprintf('Number of remaining patients to classify = %g', length(unique(data[-(x1[ ,1]) , "Pat"]))))
  
  if (length(x1)>0){
    data=data[-x1[ ,1], ]
    #  data=rbind(data, d2)
  }
  #save PatID in sep matrix
  PatID=data[ , c("Row", "Col", "Pat", "CDX2_Slide", "ImageNo")]  
  # remove this info from the data frame
  data=data[ , -c(grep("Row", colnames(data)), grep("Col", colnames(data)), 
                  grep("Pat", colnames(data)), grep("Slide", colnames(data)))]
  x1=c(grep("KER_Im", colnames(data)), grep("KER_StAreaFrac.norm", colnames(data)), grep("KER_Brown.total.norm", colnames(data)))
  data=data[ ,-x1]
 
  switch (ID,
          'Train'={
            matchID=sapply(tolower(PatID$Pat), function(x) classes[match(x, classes[ ,1]),3])
            PatID$predClass=unlist(matchID)
            print(sprintf('Rm cores with conflicting MSI and CCS2 stat= %g', length(which(PatID$predClass==2))))
            x1=which(PatID$predClass==2)
            if (length(x1)>0){
              PatID=PatID[-x1, ]
              data=data[-x1, ]
            }
            matchID=sapply(tolower(PatID$Pat), function(x) classes[match(x, classes[ ,1]),2])
            PatID$MSI=(matchID)
            print(sprintf("Final remaining patients = %g", length(unique(PatID$Pat))))
            print(sprintf("Final remaining cores = %g", length(PatID$Pat)))
          },
          'Test'={
            matchID=sapply(tolower(PatID$Pat), function(x) classes[match(x, classes[ ,1]),2])
            PatID$MSI=(matchID)
            print(sprintf("Final remaining patients = %g", length(unique(PatID$Pat))))
            print(sprintf("Final remaining cores = %g", length(PatID$Pat)))
          })  
  
  ret=list(data=data, PatID=PatID)
}

MedScale=function(data){
xmin=apply(data, 2, function(x) quantile(x, 0.025, na.rm=T))
xmax=apply(data, 2, function(x) quantile(x, 0.975, na.rm=T))
temp=sapply(1:ncol(data), function(x) (data[ ,x]-xmin[x])/xmax[x], simplify=T)
temp[temp>1.5]=1.5
rownames(temp)=rownames(data)
colnames(temp)=colnames(data)
temp
}


Compress.Scores=function(Probs, Cairo.h, meth=1, prob=0.6){
  MProbs=sapply(1:nrow(Probs), function(x) ifelse(Probs[x, 2]<prob, 1, 2))
  names(MProbs)=rownames(Probs)
  MProbs=unlist(MProbs)
  # browser()
  uni.samp <- unique(c(Cairo.h))
  report.co <- matrix(0, length(uni.samp), 3)
  for(i in 1:length(uni.samp)) {
    # print(sprintf("%g", i))
    inds <- which(Cairo.h==uni.samp[i])
    this.pred <- as.character(MProbs[inds])
    report.co[i, 1] <- sum(this.pred==1, na.rm = T)
    report.co[i, 2] <- sum(this.pred==2, na.rm = T)
    #if (sum(this.pred==2)>0){
    if (report.co[i, 2]>=report.co[i, 1]){
      report.co[i, 3] <- mean(Probs[inds[which(this.pred==2)], 2])
    } else {
      report.co[i, 3]<-mean(Probs[inds[which(this.pred==1)], 1])
    }
  }
  rownames(report.co) <- uni.samp
  colnames(report.co) <- c("CCS1", "CCS3", "Probability")
  
  pred1 <- rep(NA, nrow(report.co))
  if (meth==1){
    pred1[report.co[, "CCS3"] >= 1] <- 3
    pred1[report.co[, "CCS3"] < 1] <- 1}
  if (meth==2){
    pred1[report.co[, "CCS3"] >= report.co[, "CCS1"]] <- 3
    pred1[report.co[, "CCS3"] < report.co[, "CCS1"]] <- 1
  }
  if (meth==3){
    pred1 <- rep(1, nrow(report.co))
    pred1[as.numeric(report.co[, "Probability"])>=prob]=3
  }
  report=cbind(unique(c(Cairo.h)), report.co, pred1)
  colnames(report)=c("Pat", "ccs1.count", "ccs3.count", "certainty", "CCS")
  rownames(report)=unique(c(Cairo.h))
  #idx=match(MSIpos, rownames(report))
  #report[idx, 5]=2
  report
}
