#!/usr/bin/env Rscript
err=as.numeric(commandArgs(trailingOnly=TRUE)[1])#genotyping error rate
fam=as.numeric(commandArgs(trailingOnly=TRUE)[2])#number of simlated families
miss=as.numeric(commandArgs(trailingOnly=TRUE)[3])#missing SNPs rate
library(data.table)
library(ranger)
library(aricode)
data=as.data.frame(fread(paste("/home/meng/Ped-sim/merge_EUR_miss",miss,"_sim_table_multi_measurements_err0_1_",fam,".txt",sep=''),head=T))
table=matrix(NA,20,ncol(data)-4)
rep=1
cycles=10
n_method=0
r0=0
base=c(4)
rest=5:ncol(data)-4
for(k in 1:17){
  if(length(rest)>0){
    for(p in rest){
      # x=data.frame(data[,-c(1,2,3,(p+6))])
      x=data.frame(data[,c(base,(p+4))])
  x$pedigree.degree=factor(x$pedigree.degree)
      for(j in err){
        results=matrix(NA,nrow=cycles, ncol=rep)
        data1=as.data.frame(fread(paste("/home/meng/Ped-sim/merge_EUR_miss",miss,"_sim_table_multi_measurements_err",j,"_1_",fam,".txt",sep=''),head=T))
       x1=data.frame(data1[,c(base,(p+4))])
  x1$pedigree.degree=factor(x1$pedigree.degree)
	for(i in 1:rep){
          s=data.frame(matrix(sample(1:(floor(nrow(x)/cycles)*cycles)), ncol=cycles))
          for(r in 1:cycles){
            test=s[,r]
rf <- ranger(pedigree.degree ~ ., data = x[-test,])
            myY_P = predict(rf,x1[test,])$predictions
#results[r,i]=ARI(x1[test,1],as.character(myY_P))
results[r,i]=length(which(x1[test,1]==as.character(myY_P)))/length(myY_P)
	  }
        }
        table[k,p]=mean(colMeans(results))
      }
    }
    #if((max(table[k,])-r0)>0.01){
      index=which(table[k,]==max(na.omit(table[k,])))[1]
      base=c(base,(index+4))
      rest=rest[-which(rest==index)]
    #}
    table[18,k]=max(na.omit(table[k,]))
    table[19,k]=index
    table[20,k]=min(na.omit(table[k,]))
  }
}
write.table(table,paste("/home/meng/Ped-sim/EUR_miss",miss,"_model_selection_RF_table_degree_noBase_err",err,".",fam,"fam.txt",sep=''),quote=F,row.names = F,col.names = F)

