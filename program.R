x<-'TTCATA'
y<-'TGCTCGTA'
global_align<-function(x,y,match,mismatch,gap){
  m<-matrix(0,nrow = (nchar(x)+1), ncol = (nchar(y)+1))
  n<-matrix(0,nrow = (nchar(x)+1), ncol = (nchar(y)+1))
  l<-matrix(0,nrow = (nchar(x)+1), ncol = (nchar(y)+1))
  for(j in 2:ncol(m)){
    m[1,j]<-(m[1,j-1]+gap)
    n[1,j]<-'←'
    l[1,j]<-'left'
  }
  for(i in 2:nrow(m)){
    m[i,1]<-(m[i-1,1]+gap)
    n[i,1]<-'↑'
    l[i,1]<-'up'
  }
  for(i in 2:nrow(m)){
    for(j in 2:ncol(m)){
      xl<-unlist(strsplit(x,split = ""))
      yl<-unlist(strsplit(y,split = ""))
      if(xl[i-1]==yl[j-1]){s1<-m[i-1,j-1]+match} else {s1<-m[i-1,j-1]+mismatch}
      s2<-(m[i-1,j]+gap)
      s3<-(m[i,j-1]+gap)
      m[i,j]<-max(s1,s2,s3)
      if(s1==m[i,j]){n[i,j]<-'↖';l[i,j]<-'diag'} 
      if(s2==m[i,j]){n[i,j]<-'↑';l[i,j]<-'up'}
      if(s3==m[i,j]){n[i,j]<-'←';l[i,j]<-'left'}
      if(m[i,j]==s1&&m[i,j]==s2){n[i,j]<-paste0('↖','↑');l[i,j]<-'diag_up'}
      if(m[i,j]==s1&&m[i,j]==s3){n[i,j]<-paste0('↖','←');l[i,j]<-'diag_left'}
      if(m[i,j]==s2&&m[i,j]==s3){n[i,j]<-paste0('↑','←');l[i,j]<-'up_left'}
      if(m[i,j]==s1&&m[i,j]==s2&&m[i,j]==s3){n[i,j]<-paste0('↖','↑','←');l[i,j]<-'diag_up_left'}
    }
  }
  cat("Sequence x:",x,"\n")
  cat("Sequence y:",y,"\n")
  cat("Scoring system:",match,"for match,",mismatch,"for mismatch,",gap,"for gap","\n","\n")
  cat("Dynamic programming matrix 1:","\n")
  print(m)
  cat("Dynamic programming matrix 2:","\n")
  print(n)
  cat("Dynamic programming matrix 3:","\n")
  print(l)
  writeLines(paste0("\n","Alignment:"))
  lx<-unlist(strsplit(x,split = ""))
  ly<-unlist(strsplit(y,split = ""))
  L2=length(lx)
  L1=length(ly)
  D<-c('left','up','diag')
  if(l[L2+1,L1+1]==D[1]){ 
    RX2<-c(ly[L1]);RX1<-c("-");i=L2+1;j=L1 
  } else {
    if(l[L2+1,L1+1]==D[2]){
      RX2<-c("-");RX1<-c(lx[L2]);i=L2;j=L1+1
    } else { 
      RX2<-c(ly[L1]);RX1<-c(lx[L2]);i=L2;j=L1 
    } 
  } 
  while((i>1)&&(j>1)){ 
    # browser() 
    if(l[i,j]==D[1]){ 
      RX2<-c(ly[j-1],RX2);RX1<-c("-",RX1);j=j-1 
    } 
    else if(l[i,j]==D[2]){ 
      RX2<-c("-",RX2);RX1<-c(lx[i-1],RX1);i=i-1 
    } 
    else {RX2<-c(ly[j-1],RX2);RX1<-c(lx[i-1],RX1);j=j-1;i=i-1} 
  }
  RX3<-c()
  for(r in 1:length(RX1)){
    if(RX1[r]==RX2[r]){RX3<-c(RX3,"|")}
    else {RX3<-c(RX3," ")}
  }
  #hamming_distance<-(max(nchar(x),nchar(y))-length(grep("[|]",RX3)))
  if(length(lx)>length(ly)){lmer<-ly;s<-lx} else {lmer<-lx;s<-ly}
  hamming_distance<-numeric(length(s)-length(lmer)+1)
  for(i in 1:(length(s)-length(lmer)+1)){
    ss<-s[i:(i+length(lmer)-1)]
    for(j in 1:length(lmer)){
      if(lmer[j]!=ss[j]){hamming_distance[i]<-(hamming_distance[i]+1)}
    }
  }
  hamming_distance<-min(hamming_distance)
  cat(" x: ",RX1,"\n","   ",RX3,"\n","y: ",RX2,"\n")
  writeLines(paste0("\n","#1 score: ",m[L2+1,L1+1],"\n","#2 hamming-distance: ",hamming_distance))
}
#writeLines(paste0(paste0(unlist(strsplit(ss,split = "\n"))[1],unlist(strsplit(ss2,split = "\n"))[1]),"\n",paste0(unlist(strsplit(ss,split = "\n"))[2],unlist(strsplit(ss2,split = "\n"))[2]),"\n",paste0(unlist(strsplit(ss,split = "\n"))[3],unlist(strsplit(ss2,split = "\n"))[3])))
z<-global_align(x,y,5,-2,-6)



###################################################################################################
#read data 
data<-read.table("blast1.txt",header=T,fill=F) 
blastdata<-as.matrix(data) 
X1= blastdata[,1] 
X2= blastdata[,2] 
L1=length(X1) 
for(i in 1:length(X1)){ 
  if(X1[i]=="-") 
    L1=L1-1 
} 
L2=length(X2) 
for(i in 1:length(X2)){ 
  if(X2[i]=="-") 
    L2=L2-1 
} 
SCORE<-matrix(0,nrow=L2+1,ncol=L1+1)#score matrix 
DIRE<-matrix(0,nrow=L2+1,ncol=L1+1)#record direction for later traceback 
S<-c(10,-2,-5)#match,dismatch,space 
D<-c('left','up','diag') 

#fill score matrix and record direction in DIRE 
for(i in 1:L2+1){ 
  SCORE[i,1]=S[3]*(i-1) 
  DIRE[i,1]=D[2] 
} 
for(j in 1:L1+1){ 
  SCORE[1,j]=S[3]*(j-1) 
  DIRE[1,j]=D[1] 
} 
for(i in 2:(L2+1)){ 
  for(j in 2:(L1+1)){ 
    left=SCORE[i,j-1]+S[3] 
    up=SCORE[i-1,j]+S[3] 
    if(X1[j-1]==X2[i-1]) 
      diag=SCORE[i-1,j-1]+S[1] 
    else 
      diag=SCORE[i-1,j-1]+S[2] 
    max=max(left,up,diag) 
    SCORE[i,j]=max 
    if(left==max) DIRE[i,j]=D[1] 
    else if(up==max) DIRE[i,j]=D[2] 
    else if(diag==max) DIRE[i,j]=D[3] 
  } 
} 

#print out alignment 
if(DIRE[L2+1,L1+1]==D[1]){ 
  RX1<-c(X1[L1]);RX2<-c("-");i=L2+1;j=L1 
} else {if(DIRE[L2+1,L1+1]==D[2]){ 
  RX1<-c("-");RX2<-c(X2[L2]);i=L2;j=L1+1}else{ 
    RX1<-c(X1[L1]);RX2<-c(X2[L2]);i=L2;j=L1 
  } 
} 

#traceback 
while((i>1)&&(j>1)){ 
  # browser() 
  if(DIRE[i,j]==D[1]){ 
    RX1<-c(X1[j-1],RX1);RX2<-c("-",RX2);j=j-1 
  } 
  else if(DIRE[i,j]==D[2]){ 
    RX1<-c("-",RX1);RX2<-c(X2[i-1],RX2);i=i-1 
  } 
  else {RX1<-c(X1[j-1],RX1);RX2<-c(X2[i-1],RX2);j=j-1;i=i-1} 
} 
#printout 
cat(" 2 seq: ","\n"," X1: ",X1,"\n"," X2: ",X2,"\n","alignment","\n"," RX1: ",RX1,"\n"," RX2: ",RX2,"\n"," score:", SCORE[L2+1,L1+1])
###################################################################################################


