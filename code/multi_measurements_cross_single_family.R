#!/usr/bin/env Rscript
name_snp1 = commandArgs(trailingOnly=TRUE)[1]
name_snp2 = commandArgs(trailingOnly=TRUE)[2]
err=as.numeric(commandArgs(trailingOnly=TRUE)[3])
rep=as.numeric(commandArgs(trailingOnly=TRUE)[4])
#miss=as.numeric(commandArgs(trailingOnly=TRUE)[5])
#m=as.numeric(commandArgs(trailingOnly=TRUE)[6])
miss=0
m=100
library(data.table)
MM <- function(N1,N2,data){
  a=na.omit(data)
#####Jacquard 1-9
  j1=length(intersect(which(a[,1]==0),which(a[,2]==0)))
  j2=length(intersect(which(a[,1]==1),which(a[,2]==0)))
  j3=length(intersect(which(a[,1]==2),which(a[,2]==0)))
  j4=length(intersect(which(a[,1]==0),which(a[,2]==1)))
  j5=length(intersect(which(a[,1]==1),which(a[,2]==1)))
  j6=length(intersect(which(a[,1]==2),which(a[,2]==1)))
  j7=length(intersect(which(a[,1]==0),which(a[,2]==2)))
  j8=length(intersect(which(a[,1]==1),which(a[,2]==2)))
  j9=length(intersect(which(a[,1]==2),which(a[,2]==2)))
  n1=j5
  n2=j3+j7
  n3=j2+j5+j8
  n4=j4+j5+j6
  n5=min(c(n3,n4))
  kin=(n1-2*n2)/(n3+n4)
  kin0=0.5-(0.25*(n3+n4-2*n1)+n2)/n5
  ibs0=(j3+j7)/nrow(a)
  ibs1=(j2+j4+j6+j8)/nrow(a)
  ibs2=(j1+j5+j9)/nrow(a)
  ibs02=ibs0+ibs2
  ibs12=ibs1+ibs2
  ibs01=ibs0+ibs1
  results=c(kin,kin0,ibs0,ibs1,ibs2,ibs02,ibs12,ibs01,j1/nrow(a),j2/nrow(a),j3/nrow(a),j4/nrow(a),j5/nrow(a),j6/nrow(a),j7/nrow(a),j8/nrow(a),j9/nrow(a))
  return(results)
}

table=as.matrix(read.table("sim_table_multi_measurements_balance.txt",head=T))

data_snp1=data.frame(fread(name_snp1,head=T)[,-(1:6)])
data_snp2=data.frame(fread(name_snp2,head=T)[,-(1:6)])
  for(k in 1:nrow(table)){
	temp1=gsub("-",'.',table[k,1])
        temp1=paste("X0",temp1,sep='_')
        N1=which(colnames(data_snp1)==temp1)
	temp2=gsub("-",'.',table[k,2])
        temp2=paste("X0",temp2,sep='_')
        N2=which(colnames(data_snp2)==temp2)
        data_snp=cbind(data_snp1[,N1],data_snp2[,N2])
        table[k,5:ncol(table)]=MM(N1,N2,data_snp)
  }
table[,1]=paste("fam",rep,table[,1],sep='_')
table[,2]=paste("fam",rep,table[,2],sep='_')
colnames(table)=c("individual_1","individual_2","pedigree-relationship","pedigree-degree","kin","kin0","ibs0","ibs1","ibs2","ibs02","ibs12","ibs01","j1","j2","j3","j4","j5","j6","j7","j8","j9")
write.table(table,paste("miss",miss,"_sim_table_multi_measurements_err",err,"_",rep,"_",m,".txt",sep=''),quote=F,row.names = F,col.names = T)



