#!/usr/bin/env Rscript

j=as.numeric(commandArgs(trailingOnly=TRUE)[1])
name1=commandArgs(trailingOnly=TRUE)[2]
name2=commandArgs(trailingOnly=TRUE)[3]
miss=as.numeric(commandArgs(trailingOnly=TRUE)[4])
fam=as.numeric(commandArgs(trailingOnly=TRUE)[5])
library(data.table)
library(ranger)

data=as.data.frame(fread(name1,head=T))
rep=1
cycles=10
base=c(4,c(1,10,16,12,15)+4)#column number of selected features
x=data.frame(data[,c(base)])
x$pedigree.degree=factor(x$pedigree.degree)
Pre=matrix(NA,nrow=nrow(x), ncol=2)
Pre[,1]=1:nrow(x)
results=matrix(NA,nrow=cycles, ncol=rep)
  data1=as.data.frame(fread(paste(name2,sep=''),head=T))
  x1=data.frame(data1[,c(base)])
  x1$pedigree.degree=factor(x1$pedigree.degree)
  for(i in 1:rep){
    s=data.frame(matrix(sample(1:(floor(nrow(x)/cycles)*cycles)), ncol=cycles))
    for(r in 1:cycles){
      test=s[,r]
    rf <- ranger(pedigree.degree ~ ., data = x[-test,])
            myY_P = predict(rf,x1[test,])$predictions
    Pre[test,2]=as.character(myY_P)
results[r,i]=cor(x1[test,1],as.character(myY_P))
    }
  }
 # print(mean(results[,1]))
results=as.data.frame(results)
colnames(results)=c("index","accuracy")
  write.table(results,paste("EUR_miss",miss,"_final_degree_RF_err",j,"_1_100.txt",sep=''),quote=F,row.names = F,col.names = T)
#write.table(Pre,paste("EUR_miss",miss,"_final_degree_err",j,"_prediction.txt",sep=''),quote=F,row.names = F,col.names = T)


