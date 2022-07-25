#!/usr/bin/env Rscript

j=as.numeric(commandArgs(trailingOnly=TRUE)[1])
name1=commandArgs(trailingOnly=TRUE)[2]
name2=commandArgs(trailingOnly=TRUE)[3]
library(data.table)
library(ranger)

data=as.data.frame(fread(name1,head=T))
index=c(which(data[,3]=="PO"),which(data[,3]=="FS"))

data=data[index,]
rep=1
cycles=10
base=c(3)
rest=c(1,2,5,3,4,7,8,6,9,13,15,17,10,14,16)+4

results=matrix(NA,nrow=cycles, ncol=length(rest))
  data1=as.data.frame(fread(paste(name2,sep=''),head=T))
data1=data1[index,]
for(i in 1:length(rest)){
	base=c(base,rest[i])
	x=data.frame(data[,c(base)])
x$pedigree.relationship=factor(x$pedigree.relationship)
  x1=data.frame(data1[,c(base)])
  x1$pedigree.relationship=factor(x1$pedigree.relationship)
    s=data.frame(matrix(sample(1:(floor(nrow(x)/cycles)*cycles)), ncol=cycles))
    for(r in 1:cycles){
      test=s[,r]
    rf <- ranger(pedigree.relationship ~ ., data = x[-test,])
            myY_P = predict(rf,x1[test,])$predictions
results[r,i]=length(which(x1[test,1]==as.character(myY_P)))/length(test)

    }
    }
results1=as.data.frame(matrix(NA,nrow=2, ncol=length(rest)))
for(i in 1:length(rest)){
	results1[1,i]=paste("add_feature_",rest[i]-4,sep='')
	results1[2,i]=mean(results[,i])
}
write.table(results1,paste("EUR_best_final_RF_table_relationship_degree1_err",j,".100fam.txt",sep=''),quote=F,row.names = F,col.names = F)
