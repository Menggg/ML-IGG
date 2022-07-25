#!/usr/bin/env Rscript
err=as.numeric(commandArgs(trailingOnly=TRUE)[1])
fam=as.numeric(commandArgs(trailingOnly=TRUE)[2])
miss=as.numeric(commandArgs(trailingOnly=TRUE)[3])
library(data.table)
library(ranger)
library(aricode)
data=as.data.frame(fread(paste("merge_EUR_miss",miss,"_sim_table_multi_measurements_err0_1_",fam,".txt",sep=''),head=T))
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
        data1=as.data.frame(fread(paste("merge_EUR_miss",miss,"_sim_table_multi_measurements_err",j,"_1_",fam,".txt",sep=''),head=T))
       x1=data.frame(data1[,c(base,(p+4))])
  x1$pedigree.degree=factor(x1$pedigree.degree)
	for(i in 1:rep){
          s=data.frame(matrix(sample(1:(floor(nrow(x)/cycles)*cycles)), ncol=cycles))
          for(r in 1:cycles){
            test=s[,r]
system(paste("mkdir -p degree/relationship_miss_",miss,"_err_",j,"_rep_",r,sep=''))
system(paste("cp SVM.py degree/relationship_miss_",miss,"_err_",j,"_rep_",r,sep=''))
setwd(paste("degree/relationship_miss_",miss,"_err_",j,"_rep_",r,sep=''))
  write.csv(x1[test,],"data1.csv",quote=F,row.names = F)
write.csv(x[-test,],"data.csv",quote=F,row.names = F)
 write.csv(x[-test,1],"label.csv",quote=F,row.names = F)
system("python3 SVM.py data.csv data1.csv label.csv .")
myY_P=read.csv("SVM_pred.csv")[,1]
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
write.table(table,paste("EUR_miss",miss,"_model_selection_svm_table_degree_noBase_err",err,".",fam,"fam.txt",sep=''),quote=F,row.names = F,col.names = F)

