#!/usr/bin/env Rscript

j=as.numeric(commandArgs(trailingOnly=TRUE)[1])
name1=commandArgs(trailingOnly=TRUE)[2]
name2=commandArgs(trailingOnly=TRUE)[3]
miss=0
library(data.table)
library(ranger)

data=as.data.frame(fread(name1,head=T))
rep=1
cycles=10
base=4
#rest save all the top features from of different error rates
rest=c(1,16,12,15,10,3,11,13,2,4,6,7)+4
results=matrix(NA,nrow=cycles, ncol=length(rest))
  data1=as.data.frame(fread(paste(name2,sep=''),head=T))
for(i in 1:length(rest)){
      base=c(base,rest[i])
x=data.frame(data[,c(base)])
x$pedigree.degree=factor(x$pedigree.degree)
x1=data.frame(data1[,c(base)])
  x1$pedigree.degree=factor(x1$pedigree.degree)
    s=data.frame(matrix(sample(1:(floor(nrow(x)/cycles)*cycles)), ncol=cycles))
    for(r in 1:cycles){
      test=s[,r]
    rf <- ranger(pedigree.degree ~ ., data = x[-test,])
            myY_P = predict(rf,x1[test,])$predictions
results[r,i]=length(which(x1[test,1]==as.character(myY_P)))/length(test)
    }
  }
 # print(mean(results[,1]))
results1=as.data.frame(matrix(NA,nrow=2, ncol=length(rest)))
for(i in 1:length(rest)){
	results1[1,i]=paste("add_feature_",rest[i]-4,sep='')
	results1[2,i]=mean(results[,i])
}
  write.table(results1,paste("EUR_err",j,"_best_final_degree_RF.txt",sep=''),quote=F,row.names = F,col.names = F)



