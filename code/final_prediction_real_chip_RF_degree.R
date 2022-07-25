#!/usr/bin/env Rscripti
list=c(0,0.01,0.03,0.05,0.07,0.1)
#list=c(0.3,0.4,0.5)
results=matrix(NA,9,length(list))
for(err in 1:length(list)){
name1=paste("merge_EUR_miss0_sim_table_multi_measurements_err",list[err],"_1_100.txt",sep='')
#name2="real_chip_table_multi_measurements.txt"
name2="real_chip_table_multi_measurements_GSA.txt"
library(data.table)
library(ranger)
data=as.data.frame(fread(name1,head=T))
base=c(4,c(1,10,16,12,15)+4)
#base=c(4,5)
x=data.frame(data[,c(base)])
x$pedigree.degree=factor(x$pedigree.degree)
  data1=as.data.frame(fread(paste(name2,sep=''),head=T))
  x1=data.frame(data1[,c(base)])
  x1$pedigree.degree=factor(x1$pedigree.degree)
  rf <- ranger(pedigree.degree ~ ., data = x)
            myY_P = predict(rf,x1)$predictions
K=data1[,5:8]
for(i in 1:nrow(K)){
    if(K[i,1]<1/(2^1.5) & K[i,1]>=1/(2^2.5))
      K[i,3]=0.25
    else if (K[i,1]<1/(2^2.5) & K[i,1]>=1/(2^3.5))
      K[i,3]=0.125
    else if (K[i,1]<1/(2^3.5) & K[i,1]>=1/(2^4.5))
      K[i,3]=0.0625
    else if (K[i,1]<1/(2^4.5) & K[i,1]>=-1)
      K[i,3]=0
  }
for(i in 1:nrow(K)){
    if(K[i,2]<1/(2^1.5) & K[i,2]>=1/(2^2.5))
      K[i,4]=0.25
    else if (K[i,2]<1/(2^2.5) & K[i,2]>=1/(2^3.5))
      K[i,4]=0.125
    else if (K[i,2]<1/(2^3.5) & K[i,2]>=1/(2^4.5))
      K[i,4]=0.0625
    else if (K[i,2]<1/(2^4.5) & K[i,2]>=-1)
      K[i,4]=0
  }
r=cbind(as.character(x1[,1]),as.character(myY_P),K[,3:4])
results[1,err]=length(which(r[1:28,2]==r[1:28,1]))/28
results[2,err]=length(which(r[1:28,3]==r[1:28,1]))/28
results[3,err]=length(which(r[1:28,4]==r[1:28,1]))/28
results[4,err]=length(which(r[29:56,2]==r[29:56,1]))/28
results[5,err]=length(which(r[29:56,3]==r[29:56,1]))/28
results[6,err]=length(which(r[29:56,4]==r[29:56,1]))/28
results[7,err]=length(which(r[57:84,2]==r[57:84,1]))/28
results[8,err]=length(which(r[57:84,3]==r[57:84,1]))/28
results[9,err]=length(which(r[57:84,4]==r[57:84,1]))/28
}
#write.table(results,paste("real_chip_table_results_4M_1.txt",sep=''),quote=F,row.names = F,col.names = T)
write.table(results,paste("real_chip_table_results_GSA.txt",sep=''),quote=F,row.names = F,col.names = T)
 
