#preprocess cohort data

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

ContTable=function(tab, title, chisqtest=F,ylabL="GE classifier", scaleV=T,
                   xlabL="TMA classifier"){
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
          xlab=xlabL, ylab=ylabL,
          main=sprintf("%s %s", title, tit2))
  }else{
    if (scaleV==T){
      image(scale(tab), col=brewer.pal(9, "BuGn"), xaxt="n", yaxt="n",
            xlab=xlabL, ylab=ylabL,
            main=sprintf("%s %s", title, tit2))
    }else{
      image(tab, col=brewer.pal(9, "BuGn"), xaxt="n", yaxt="n",
            xlab=xlabL, ylab=ylabL,
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


AllPatList=function(x, Marker="Panker"){
  #x is the directory with all the files
  dirIn=dir(path=x, pattern="*.csv")
  mapidx=grep("lide_no", dirIn)
  mapin=read.csv(paste(x,dirIn[mapidx], sep=""))
  rmaps=setdiff(dirIn, dirIn[mapidx])
  
  ra=regmatches(rmaps, regexpr("_[0-9]+.csv", rmaps))
  if (nchar(ra[1])==7){
    ridx=gsub(".csv","", rmaps)
  }else if (nchar(ra[1])==11){
    ridx=strsplit(rmaps, ra)
    ridx=unlist(ridx)
    
    rmaps2=sapply(mapin[ ,match(Marker, colnames(mapin))], function(x) grep(x, rmaps))
    rmaps=rmaps[rmaps2]
    ridx=ridx[rmaps2]
  }
  
  
  temp=data.frame()
  for (i in 1:length(rmaps)){
    a1=read.csv(paste(x, rmaps[i], sep=""), header=F, stringsAsFactors = F)
    nR=nrow(a1)
    nC=ncol(a1)
    a1=stack(a1)
    a1$Row=rep(c(1:nR), nC)
    a1$Col=rep(c(1:nC), each=nR)
    a1$SlideName=ridx[i]
    a1$SlideNo=match(ridx[i], mapin$TMA_slide_Number)
    a1$Aperio=mapin[match(ridx[i], mapin$TMA_slide_Number), match(Marker, colnames(mapin))]
    temp=rbind(temp, a1)
  }
  temp
}
