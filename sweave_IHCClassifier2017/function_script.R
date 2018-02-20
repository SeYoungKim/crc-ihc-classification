#library("gsubfn")
library(e1071)
library(nnet)
library(preprocessCore)
library(randomForest)
library(ROCR)
library(rpart)
library(survival)
# remove ggplot2 if not necessary
library(ggplot2)
library(reshape2)
library(knitr)
library(RColorBrewer)
library(scatterplot3d)
library(MASS)
library(clinfun)
library(irr)
library(heatmap.plus)
library(psych)

subCols <- c("skyblue", "palegreen2", "darkblue","dodgerblue2","#E31A1C", "green4",
             "#6A3D9A","#FF7F00", "black","gold1",
             "skyblue2","#FB9A99", "palegreen2",
             "#CAB2D6", "#FDBF6F", "gray70", "khaki2",
             "maroon","orchid1","deeppink1","blue1","steelblue4",
             "darkturquoise","green1","yellow4","yellow3",
             "darkorange4","brown")
palette(subCols)

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
    tempA=read.table(sprintf("%s%s", DirPath, DirNames[i]), sep="\t", header=T, stringsAsFactors =F)
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
  Slide_No=read.csv(sprintf("%s%s", mapdir, MFiles[xslide]), header=T, stringsAsFactors =F)
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
      mapping=read.csv(sprintf('%s%s', mapdir, MFiles[idx2]), header=F)
      Maps=sapply(1:nrow(tempB), function(x) as.character(mapping[as.numeric(tempB$Row[x]),as.numeric(tempB$Col[x])]))
      tempB$Pat=Maps
      save(tempB, file=sprintf("%sRdataParts/data_%g.RData", DirPath, i))
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

PreProcess=function(data, classes, ID=NULL, rm.range=c(0.2, 0.6), na.rm=T, KER.rm=T){
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
  if (KER.rm==T){
    data=data[ ,-grep("KER", colnames(data))]}
  else {
    x1=c(grep("KER_Im", colnames(data)), grep("KER_StAreaFrac.norm", colnames(data)), grep("KER_Brown.total.norm", colnames(data)))
    data=data[ ,-x1]
  }
  switch (ID,
          'Lei'={
            rmID=c("!", "\\*", "&", "%", "\\$")
            x2=unlist(sapply(rmID, function(x) grep(x, data$Suffix)))
            data=data[-x2, ]
            PatID=PatID[-x2, ]
            PatID$Suffix=data$Suffix
            data=data[ , -c(grep("Suffix", colnames(data)), grep("Sample", colnames(data)))]
            matchID=sapply(tolower(PatID$Pat), function(x) classes[match(x, classes[ ,1]),2])
            t2=sapply(matchID, function(x) length(x))
            PatID$MSI=matchID
            print(sprintf("Final remaining patients = %g", length(unique(PatID$Pat))))
            print(sprintf("Final remaining cores = %g", length(PatID$Pat)))
          },
          'AMC'={
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
          'Cai'={
            matchID=sapply(tolower(PatID$Pat), function(x) classes[match(x, classes[ ,1]),2])
            PatID$MSI=(matchID)
            print(sprintf("Final remaining patients = %g", length(unique(PatID$Pat))))
            print(sprintf("Final remaining cores = %g", length(PatID$Pat)))
          },{
            print(sprintf("Final remaining patients = %g", length(unique(PatID$Pat))))
            print(sprintf("Final remaining cores = %g", length(PatID$Pat)))
          } )  
  
  ret=list(data=data, PatID=PatID)
}

EditCairo2Clin=function(PatInfo){
  PatInfo$EGFR.GCN=factor(PatInfo$EGFR.GCN) 
  levels(PatInfo$EGFR.GCN)=c("0", "1", "2", NA)
  PatInfo$BRAF=factor(PatInfo$BRAF)
  levels(PatInfo$BRAF)=c("0", "1", NA)
  PatInfo$Sex=factor(PatInfo$Sex*-1, labels=c("F", "M"))
  PatInfo$Arm=factor(PatInfo$Arm, labels=c("cntl", "wCetux"))
  PatInfo$PTEN=factor(PatInfo$PTEN, levels=c("0", "1", NA))
  PatInfo$KRAS_01=factor(PatInfo$KRAS_01, levels=c("0", "1", NA))
  PatInfo$HER2.GCN=factor(PatInfo$HER2.GCN, levels=c("0", "1", NA))
  PatInfo$tumorPC[grep("NA", PatInfo$tumorPC)]=NA
  levels(PatInfo$Response)=c("CR", NA, "PD", "PR", "SD", NA)
  PatInfo$Response=factor(PatInfo$Response, levels=c("CR", "PR", "SD", "PD"))
  PatInfo$Age=as.numeric(as.character(PatInfo$Age))
  PatInfo$Age[PatInfo$Age<0]=NA
  PatInfo$AgeCat=cut(PatInfo$Age, c(0, 50, 60, 70, 100), c("<50", "50-60", "60-70", "70+"))
  PatInfo
}

Cairo1Edit=function(data){
  colnames(data)[c(6:8, 12,14,16)]=c("Response.First", "Response.Second",
                                     "Response.Third", "PFS1.days", "PFS2.days",
                                     "PFS3.days")
  levels(data$Geslacht) <- list("F"="v", "M"="m", "NA"=".")
  levels(data$Geslacht)=c("F", "M", NA)
  levels(data$treatment.Arm)=list("ArmA"="Sequential chemotherapy (Arm A)", 
                                      "ArmB"="Combination chemotherapy (Arm B)")
  data$Age2=cut(data$Age, breaks=c(0, 50, 60, 70, 120), 
                labels=c("0-50","50-60", "60-70", "70+"))
  data$OS_MONTH=data$OS_DAY/30.25
  levels(data$Response.First)=c(NA, "CR", NA, NA, NA, "SD", NA, "PD", "PR")
  data$Response.First=factor(data$Response.First, levels=c("CR", "PR", "SD", "PD"))
  levels(data$Response.Second)=c(NA, "CR", NA, NA, NA, "SD", NA, "PD", "PR")
  data$Response.Second=factor(data$Response.Second, levels=c("CR", "PR", "SD", "PD"))
  levels(data$Response.Third)=c(NA,  NA,  NA, "SD", NA, "PD", "PR")
  data$Response.Third=factor(data$Response.Third, levels=c("CR", "PR", "SD", "PD"))
  #OverallBestResponse
  data$OverallBestResponse=factor(paste(data$Response.First, data$Response.Second, sep=""))
  levels(data$OverallBestResponse)=c("CR", "CR", "CR", "CR", "CR", NA, "PD", "PR",
                                     "SD", "PD", "PD", "PR", "SD", "CR", "PR", "PR", 
                                     "PR", "PR", "CR", "SD", "SD", "PR", "SD")
  data$OverallBestResponse=factor(paste(data$OverallBestResponse, data$Response.Third, sep=""))
  levels(data$OverallBestResponse)=c("CR", "CR", "CR",NA, "PD", "PD", "PR",
                                     "SD", "PR", "PR", "PR", "PR", "SD", "SD", "PR","SD")
  data$OverallBestResponse=factor(data$OverallBestResponse, levels=c("CR", "PR", "SD", "PD"))
  
  data 
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

PathComparison=function(inputData, classInfo, destxt,...){
  Lvls=lapply(2:5, function(x) sort(unique(na.omit(classInfo[ ,x]))))
  nLvl=sapply(Lvls, function(x) length(x))
  boxplot(NA, xlim=c(0.5, sum(nLvl)+0.5), ..., yaxt='n', horizontal=T)
  idx=toupper(unlist((strsplit(colnames(inputData), "_")))) [seq(1, 10, 2)]
  #browser()
  idx2=toupper(colnames(classInfo))
  idx3=na.omit(match(idx2, idx))
  idx2=which(!is.na(match(idx2, idx)))
  MVal=min((...)[2]*0.9, quantile(max(inputData), 0.98, na.rm=T), na.rm=T)
  colVals=list(colorRampPalette(c("darkblue", "skyblue"))(nLvl[1]),
            colorRampPalette(c("darkgreen", "palegreen"))(nLvl[2]),
            colorRampPalette(c("darkgoldenrod1", "khaki1"))(nLvl[3]),
            colorRampPalette(c("darkred", "pink"))(nLvl[4]))
  atIdx=list(c(1:nLvl[1]), c(sum(nLvl[1],1):sum(nLvl[1:2])), 
              c(sum(nLvl[1:2], 1):sum(nLvl[1:3])), c(sum(nLvl[1:3],1):sum(nLvl)))
#     "skyblue","blue", "darkblue", "palegreen", "forestgreen", "darkgreen",
#             "khaki1", "yellow", "darkgoldenrod1", "pink", "red", "darkred")
#browser()
  for (i in 1:4){
    boxplot(inputData[, idx3[i]]~classInfo[ ,idx2[i]], horizontal=T,
            add=T,col=colVals[[i]], at=atIdx[[i]], xaxt='n', yaxt='n')
  #  browser()
    for (j in 1:(length(Lvls[[i]])-1)){
      temp=Lvls[[i]]
      ind1=which(classInfo[ ,idx2[i]]==temp[j])
      ind2=which(classInfo[ ,idx2[i]]==temp[j+1])
      t1=wilcox.test(inputData[c(ind1, ind2), idx3[i]]~classInfo[c(ind1, ind2) ,idx2[i]])
        if (t1$p.val<0.05){
        text(MVal, atIdx[[i]][j]+0.5, sprintf("%g", round(t1$p.val*1000)/1000), cex=0.7)
      }
    }
  }
  axis(2, at=unlist(atIdx), c(paste(rep(colnames(classInfo)[2:5], nLvl), unlist(Lvls))), las=2)
  title(main=sprintf("%s", destxt))
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

ContTable=function(tab, title, chisqtest=F,ylabL="GE classifier", scaleV=T){
  library("RColorBrewer")
  if (chisqtest==T){
    a1=chisq.test(tab)
    tit2=sprintf("Chisq = %g", round(a1$p.value*100)/100)
  } else {
    tit2=" "
  }
  nr=nrow(tab)
  nc=ncol(tab)
  if (chisqtest==T){
    l1=(a1$observed-a1$expected)
    l1[which(l1<0, arr.ind=T)]=0
    image(l1, col=brewer.pal(9, "BuGn"), xaxt="n", yaxt="n",
          xlab="TMA classifier", ylab=ylabL,
          main=sprintf("%s %s", title, tit2))
  }else{
    if (scaleV==T){
  image(scale(tab), col=brewer.pal(9, "BuGn"), xaxt="n", yaxt="n",
        xlab="TMA classifier", ylab=ylabL,
        main=sprintf("%s %s", title, tit2))
  }else{
    image(tab, col=brewer.pal(9, "BuGn"), xaxt="n", yaxt="n",
          xlab="TMA classifier", ylab=ylabL,
          main=sprintf("%s %s", title, tit2))   
    }}
  xval=seq(0, 1, 1/(nr-1))
  yval=seq(0, 1, 1/(nc-1))
  axis(1, at=xval, rownames(tab))
  axis(2, at=yval, colnames(tab))
  for(i in 1:nr){
    for (j in 1:nc){
      text(xval[i], yval[j], tab[i,j], cex=1.5)
    }
  }
}

CollapseMethod=function(m1, m2, m3, data){
  idx=which(is.na(m1[ ,1]))
  if (length(idx)>0){
    C2_RPa=Compress.Scores(m1[-idx, ], data[-idx],  1)
    C2_RPb=Compress.Scores(m1[-idx, ], data[-idx],  3, 0.6)}
  else{
    C2_RPa=Compress.Scores(m1, data, 1)
    C2_RPb=Compress.Scores(m1, data, 3)
  }
  m1out=as.numeric(C2_RPa[ ,5])+as.numeric(C2_RPb[ ,5])
  idx=which(is.na(m2[ ,1]))
  if (length(idx)>0){
    C2_NPa=Compress.Scores(m2[-idx, ], data[-idx],  1)
    C2_NPb=Compress.Scores(m2[-idx, ], data[-idx],  3, 0.6)
  }else{
    C2_NPa=Compress.Scores(m2, data,1)
    C2_NPb=Compress.Scores(m2, data,3, 0.6)}
  
  m2out=as.numeric(C2_NPa[ ,5])+as.numeric(C2_NPb[ ,5])
  idx=which(is.na(m3[ ,1]))
  if (length(idx)>0){
    C2_QPa=Compress.Scores(m3[-idx, ], data[-idx],1)
    C2_QPb=Compress.Scores(m3[-idx, ], data[-idx], 3, 0.6)
  }else{
    C2_QPa=Compress.Scores(m3, data, 1)
    C2_QPb=Compress.Scores(m3, data, 3, 0.6)
  }
  m3out=as.numeric(C2_QPa[ ,5])+as.numeric(C2_QPb[ ,5])
  PatList=unique(c(C2_RPa[ ,1], C2_NPa[ ,1], C2_QPa[ ,1]))
  PatList=data.frame(PatList, range=NA, z=NA, qn=NA)
  x1=match(C2_RPa[ ,1], PatList$PatList)
  PatList$range[x1]=m1out
  x1=match(C2_NPa[ ,1], PatList$PatList)
  PatList$z[x1]=m2out
  x1=match(C2_QPa[ ,1], PatList$PatList)
  PatList$qn[x1]=m3out
  PatList
}

LeidenCleanUp=function(LeiClin){
  LeiClin$ID=gsub(" ", "", LeiClin$ID)
  # change the DFS to 60 months
  LeiClin$DFS_month60=LeiClin$DFS_months
  LeiClin$DFS_event60=LeiClin$DFS_event
  x1=which(LeiClin$DFS_month60>60)
  LeiClin$DFS_event60[x1[LeiClin$DFS_event60[x1]==1]]=0
  LeiClin$DFS_month60[x1]=60
  # change the OS to 60 months 
  LeiClin$OS_month60=LeiClin$DFS_months
  LeiClin$OS_event60=LeiClin$DFS_event
  x1=which(LeiClin$OS_month60>60)
  LeiClin$OS_event60[x1[LeiClin$OS_event60[x1]==1]]=0
  LeiClin$OS_month60[x1]=60
  # Clean up TNM staging system
  LeiClin$TNM=factor(LeiClin$TNM)
  levels(LeiClin$TNM)=c("1", "2", "3", "4", NA)
  levels(LeiClin$AgeCat)=list("0-50"="0-50", "50-60"="50-60",
                              "60-70"="60-70", "70+"=">70")
  LeiClin
}

ClassifyPrediction=function(tdat, N=100, ptype, meth='none', nclass=2){
  Accuracy=rep(NA, N)
  Labels=matrix(NA, nrow(tdat), N)#, dimnames=list(rownames(tdat), 1:N))
  ImpMat=matrix(NA, ncol(tdat)-1, N) #, dimnames=list(rownames(tdat), 1:N))
  Err.Rate=rep(NA,N)
  if (nclass==3){
    LabelsC2=matrix(NA, nrow(tdat), N)
    LabelsC3=matrix(NA, nrow(tdat), N)
  }
  rownames(Labels)=rownames(tdat)
  nfeat=ncol(tdat)-1
  rownames(ImpMat)=colnames(tdat)[1:nfeat]
  for (i in 1:N){
    Idx=sample(1:nrow(tdat), round(nrow(tdat)/3*2))
    TrainSet=tdat[Idx, ]
    TestSet=tdat[-Idx, ]
    switch(meth,
           'none'={},
           'QN'={
             nlength=nrow(TestSet)
             TS=lapply(1:nfeat, function(x)
               normalize.quantiles(cbind(TrainSet[ ,x], TestSet[ ,x])))
             TrainSet=data.frame(sapply(1:nfeat, function(x) TS[[x]][ ,1]))
             TestSet=data.frame(sapply(1:nfeat, function(x) TS[[x]][ ,2]))
             TrainSet$CCS=tdat$CCS[Idx]
             TestSet=TestSet[ 1:nlength, ]
             TestSet$CCS=tdat$CCS[-Idx]
             colnames(TrainSet)=colnames(tdat)
             colnames(TestSet)=colnames(tdat)
           },
           'z'={
             TrainSet=data.frame(scale(TrainSet[ ,-(nfeat+1)]))
             TestSet=data.frame(scale(TestSet[ ,-(nfeat+1)]))
             TrainSet$CCS=tdat$CCS[Idx]
             TestSet$CCS=tdat$CCS[-Idx]
             colnames(TrainSet)=colnames(tdat)
             colnames(TestSet)=colnames(tdat)},
           'range'={
             TrainSet=as.data.frame(MedScale(TrainSet[ ,-(nfeat+1)]))
             TestSet=as.data.frame(MedScale(TestSet[ ,-(nfeat+1)]))
             TrainSet$CCS=tdat$CCS[Idx]
             TestSet$CCS=tdat$CCS[-Idx]
             colnames(TrainSet)=colnames(tdat)
             colnames(TestSet)=colnames(tdat)}
    )
    # Train the classifier
    switch(ptype,
           'RF'={
             rf1 = randomForest(CCS ~ ., data=TrainSet, importance=TRUE, ntree=1000)
             ImpMat[ ,i]=importance(rf1, type=1, scale=F)
             Err.Rate[i]=mean(rf1$err.rate[ ,1])
             rfpred =predict(rf1, TestSet[ , -(nfeat+1)], type="response") 
             rfpredB =predict(rf1, TestSet[ , -(nfeat+1)], type="prob") 
           },
           'SVM'={
             rf1 = svm(CCS ~ ., data=TrainSet, type="C-classification", kernel="linear", probability=T)    
             rfpred =predict(rf1, TestSet[ , -ncol(TestSet)], probability=T) 
             rfpredB =attr(rfpred, "prob")
             rfpredB=rfpredB[ ,c("1", "3")]
           },
           'NN'={
             rf1 = nnet(CCS ~ ., data=TrainSet, decay=0.01, size=2, maxit=200) 
             rfpred =predict(rf1, TestSet[ , -ncol(TestSet)], 'class') 
             rfpredB =predict(rf1, TestSet[ , -ncol(TestSet)]) 
           },
           'rpart'={
             rf1 = rpart(CCS ~ ., data=TrainSet, method="class") 
             rfpred =predict(rf1, TestSet[ , -ncol(TestSet)], type="class") 
             rfpredB =predict(rf1, TestSet[ , -ncol(TestSet)]) 
             temp=rf1$variable.importance
             x1=match(names(temp), rownames(ImpMat))
             ImpMat[x1,i]=temp
           },
           'LR'={
             rf1=glm(factor(CCS)~., data=TrainSet,  family = binomial(link="logit"))
             rfpredB=predict(rf1, TestSet[ , -ncol(TestSet)], type="response") 
             rfpred=rfpredB>=0.5
             rfpredB=cbind(1-rfpredB, rfpredB)
             x1=match(names(coef(rf1)), rownames(ImpMat))
             x2=which(is.na(x1))
             ImpMat[x1[-x2], i]=coefficients(rf1)[-x2]
           }
    )
    # predict outcome
    Acc=table(rfpred, factor(TestSet[ ,ncol(TestSet)]))
    Accuracy[i]=sum(diag(Acc))/sum(Acc)
    Labels[-Idx, i]=rfpredB[ ,1]
    if (nclass>2){
      LabelsC2[-Idx, i]=rfpredB[,2]
      LabelsC3[-Idx, i]=rfpredB[,3]
    }
  }
  if (nclass>2){
    ret=list(Accuracy=Accuracy, Labels=Labels, ImpMat=ImpMat, LabelsC2=LabelsC2, LabelsC3=LabelsC3)
  }else{
    ret=list(Accuracy=Accuracy, Labels=Labels, ImpMat=ImpMat, Err.Rate=Err.Rate)}
  ret  
}



VisualiseTable=function(tab, title){
  cnames=colnames(tab)[ncol(tab):1]
  tab=sapply(1:ncol(tab), function(x) as.numeric(factor(tab[ ,x])))
  image(tab[ ,ncol(tab):1], col=palette(), xaxt='n', yaxt='n', main=title)
  axis(2, at=seq(0, 1, length=ncol(tab)), labels=cnames, las=1)
}

FindAUC=function(pred, true){
  x1=prediction(pred, true)
  perf <- performance(x1, "auc") 
  unlist(perf@y.values)
}

  

CoreConsistencyPlot=function(mA, mB, mC, mD, mE,trueS){
  par(mfrow = c(5,1),
      oma = c(5,4,0,0) + 0.1,
      mar = c(0,0,1,1) + 0.1)
  mA=rowMeans(mA$Labels, na.rm=T)
  mA=cbind(mA, 1-mA)
  mB=rowMeans(mB$Labels, na.rm=T)
  mB=cbind(mB, 1-mB)
  mC=rowMeans(mC$Labels, na.rm=T)
  mC=cbind(1-mC, mC)
  mE=rowMeans(mE$Labels, na.rm=T)
  mE=cbind( mE, 1-mE)
  IDX=order(trueS[ 1, ], mA[ ,1])
  barplot(t1[,IDX ], col=c(1, 3))
  barplot(t(mA[IDX, ]), col=c(1, 3))
  barplot(t(mB[IDX, ]), col=c(1, 3))
  barplot(t(mC[IDX, ]), col=c(1, 3))
  
  if (is.null(mD)==F){
    mD=rowMeans(mD$Labels, na.rm=T)
    mD=cbind(mD, 1-mD)
    barplot(t(mD[IDX, ]), col=c(1, 3))
  }
  barplot(t(mE[IDX, ]), col=c(1, 3))
  
}

PatConsistencyPlot=function(mA, mB, mC, mD, mE, PatScore, method, val){
  par(mfrow=c(1,1))
  AUCscores=rep(NA, 4)
  mA=rowMeans(mA$Labels, na.rm=T)
  AUCscores[1]=FindAUC(1-mA, PatScore[ ,2])
  mA=cbind(mA, 1-mA)
  mA=Compress.Scores(mA, as.character(PatScore[ ,1]), method, val)
  mB=rowMeans(mB$Labels, na.rm=T)
  AUCscores[2]=FindAUC(1-mB, PatScore[ ,2])
  mB=cbind(mB, 1-mB)
  mB=Compress.Scores(mB,  as.character(PatScore[ ,1]), method, val)
  mC=rowMeans(mC$Labels, na.rm=T)
  AUCscores[3]=FindAUC(mC, PatScore[ ,2])
  mC=cbind(1-mC, mC)
  mC=Compress.Scores(mC,  as.character(PatScore[ ,1]), method, val)
  if (is.null(mD)==F){
    mD=rowMeans(mD$Labels, na.rm=T)
    AUCscores[4]=FindAUC(1-mD, PatScore[ ,2])
    mD=cbind(mD, 1-mD)
    mD=Compress.Scores(mD, as.character(PatScore[ ,1]),method, val)
  }
  mE=rowMeans(mE$Labels, na.rm=T)
  AUCscores[5]=FindAUC(1-mE, PatScore[ ,2])
  mE=cbind( mE, 1-mE)
  mE=Compress.Scores(mE, as.character(PatScore[ ,1]),method, val)
  
  idx1=match(rownames(mA), PatScore[ ,1])
  TrueClass=PatScore[idx1, 2]
  names(TrueClass)=unique(PatScore[idx1, 1])
  x1=order(TrueClass)
  if (is.null(mD)==F){
    tempA=cbind( as.numeric(mA[x1 ,5]), as.numeric(mB[x1 ,5]), as.numeric(mC[ x1,5]), 
                 as.numeric(mD[x1 ,5]), as.numeric(mE[x1, 5]),as.numeric(TrueClass[x1]))
  } else{
    tempA=cbind(mA[x1 ,5], mB[x1 ,5], mC[ x1,5], TrueClass[x1])   
  }
  image(tempA, xaxt='n', yaxt='n', col=c("lightblue","darkblue" ))
  # report the accuracy
  if (is.null(mD)==F){
    NCorrect=c(sum(diag(table(TrueClass, mA[ ,5]))),
               sum(diag(table(TrueClass, mB[ ,5]))),
               sum(diag(table(TrueClass, mC[ ,5]))),
               sum(diag(table(TrueClass, mD[ ,5]))),
               sum(diag(table(TrueClass, mE[ ,5]))),
               length(TrueClass))
    axis(2, at=seq(0, 1, length=6), paste(c("RF", "SVM", "NN", "CART", "LR", "True"),
                                         NCorrect), las=1)
  }
  AUCscores
}

getConds<-function(tree, VarNameList){
  #store all conditions into a list

  #start by the terminal nodes and find previous conditions
  id.leafs<-which(tree$status==-1)
  conds<-matrix(NA, nrow=11, ncol=length(VarNameList))
  colnames(conds)=VarNameList
  j<-0
  condst=list()
  for(i in id.leafs){
    j<-j+1
    if (j>11){
      break
    }
    prevConds<-prevCond(tree,i)
    condst[[j]]<-paste(prevConds$var, prevConds$val)
    conds[j, match(prevConds$var, VarNameList)]=prevConds$val
    while(prevConds$id>1){
      prevConds<-prevCond(tree,prevConds$id)
      condst[[j]]<-paste(condst[[j]]," & ",prevConds$var, prevConds$val)
      if (is.na(conds[j, match(prevConds$var, VarNameList)])){
      conds[j, match(prevConds$var, VarNameList)]=prevConds$val}
      else{
        print('help!')
        x1=conds[j, match(prevConds$var, VarNameList)]
        x2=prevConds$val
        t1=paste(substr(x1, 1,1), substr(x2, 1, 1), sep="")
        switch(t1, 
               "<<"={ t2=min(substr(x1, 3,6), substr(x2, 3,6))
               if (t2==1.5){t2=1}else{t2=paste("<", t2)} },
               ">>"={ t2=max(substr(x1, 3,6), substr(x2, 3,6)) 
               if (t2==2.5){t2=3}else{t2=paste(">", t2)}},
               "<>"={
                if (abs(as.numeric(substr(x2, 3,6))-as.numeric(substr(x1, 3,6)))==1){
                  t2=2
                } else{
                 t2=paste(substr(x2, 3,6),"-", substr(x1, 3,6))} },
               "><"={
                 if (abs(as.numeric(substr(x2, 3,6))-as.numeric(substr(x1, 3,6)))==1){
                   t2=2
                 } else {t2=paste(substr(x1, 3,6),"-", substr(x2, 3,6))}},
               "3>"={t2=3},
               ">3"={t2=3},
               "1<"={t2=1},
               "<1"={t2=1})
        conds[j, match(prevConds$var, VarNameList)]=t2
      }
      if(prevConds$id==1){
        condst[[j]]<-paste(condst[[j]]," => ",tree$prediction[i])
        conds[j,1]=tree$prediction[i]
        break()
      }
    }
    
  }
  
  return(conds)
}

prevCond<-function(tree,i){

  if(i %in% tree$right_daughter){
    id<-which(tree$right_daughter==i)
    var=tree$split_var[id]
    if (tree$split_point[id]==2.5){
      val=3
    } else{
    val<-paste(">",tree$split_point[id])}
  }
  if(i %in% tree$left_daughter){
    id<-which(tree$left_daughter==i)
    var=tree$split_var[id]
    if (tree$split_point[id]==1.5){
      val=1
    } else{
      val<-paste("<",tree$split_point[id])}
  }
  
  return(list(var=var, val=val,id=id))
}

#remove spaces in a word
collapse<-function(x){
  x<-sub(" ","_",x)
  
  return(x)
}

PermuteTrees=function(tdat, N=100){
  Acc=rep(NA, N)
  CutVals=list()
  for (i in 1:N){
    t1=which(tdat$CCS==3)
    t2=which(tdat$CCS==1)
  l1=c(sample(t1, round(length(t1)*2/3)), sample(t2, round(length(t2)*2/3)))
  AMCclass3=rpart(factor(CCS)~., data=tdat[l1, ], method="class", 
                  control=rpart.control(minsplit=10, cp=0.01))
  l2=table(predict(AMCclass3, tdat[-l1, ])[ ,1]<0.5, tdat[-l1, "CCS"])
  Acc[i]=sum(diag(l2))/sum(l2)
  
  ## look at the variables used in the tree
  VarUsed=AMCclass3$frame$var
  VarUsed2=unique(VarUsed)
  VarUsed2=setdiff(VarUsed2, "<leaf>")
  AMCclass3$splits[grep(VarUsed2[1], rownames(AMCclass3$splits)), 3:4]
  
  ## create a table of the variables used and their associated values:
  temp=lapply(1:length(VarUsed2), function(x) (AMCclass3$splits[c(grep(VarUsed2[x], rownames(AMCclass3$splits)),1), 1:4]))
  temp2=do.call(rbind, temp) #}, error=function(e) {browser()})
  temp2=temp2[!duplicated(temp2), ]
  temp2=temp2[temp2[ ,1]!=0, ]
  
  CutVals[[i]]=temp2
  }
  ret=list(Acc=Acc, CutOff=CutVals)
  ret
}