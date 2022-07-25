#!/usr/bin/env Rscript
err=as.numeric(commandArgs(trailingOnly=TRUE)[1])
fam=as.numeric(commandArgs(trailingOnly=TRUE)[2])
miss=as.numeric(commandArgs(trailingOnly=TRUE)[3])
library(data.table)
data=as.data.frame(fread(paste("merge_EUR_miss",miss,"_sim_table_multi_measurements_err0_1_",fam,".txt",sep=''),head=T))
index=c(which(data[,3]=="GP"),which(data[,3]=="HS"),which(data[,3]=="Uncle"))
data=data[index,]
data1=as.data.frame(fread(paste("merge_EUR_miss",miss,"_sim_table_multi_measurements_err",err,"_1_",fam,".txt",sep=''),head=T))
        data1=data1[index,]
table=matrix(NA,20,ncol(data)-7)
rep=1
cycles=10
n_method=0
r0=0
base=c(3,5:7)
rest=c(8:ncol(data))-7
for(k in 1:14){
  if(length(rest)>0){
    for(p in rest){
      # x=data.frame(data[,-c(1,2,3,(p+6))])
      x=data.frame(data[,c(base,(p+7))])
      x$pedigree.relationship=factor(x$pedigree.relationship)
      for(j in err){
        results=matrix(NA,nrow=cycles, ncol=rep)
        # x1=data.frame(data1[,-c(1,2,3,(p+6))])
        x1=data.frame(data1[,c(base,(p+7))])
        x1$pedigree.relationship=factor(x1$pedigree.relationship)
        for(i in 1:rep){
          s=data.frame(matrix(sample(1:(floor(nrow(x)/cycles)*cycles)), ncol=cycles))
          for(r in 1:cycles){
            test=s[,r]
	  #temp=aggregate(. ~ pedigree.relationship, x[-test,], mean)
	  system(paste("mkdir -p relationship/relationship_miss_",miss,"_err_",j,"_rep_",r,sep=''))
system(paste("cp SVM.py relationship/relationship_miss_",miss,"_err_",j,"_rep_",r,sep=''))
setwd(paste("relationship/relationship_miss_",miss,"_err_",j,"_rep_",r,sep=''))
  write.csv(x1[test,],"data1.csv",quote=F,row.names = F)
     write.csv(x[-test,],"data.csv",quote=F,row.names = F)
 write.csv(x[-test,1],"label.csv",quote=F,row.names = F)
system("python3 SVM.py data.csv data1.csv label.csv .")
myY_P=read.csv("SVM_pred.csv")[,1]
temp1=read.csv("data1.csv")
results[r,i]=length(which(temp1[,1]==as.character(myY_P)))/nrow(temp1)
	  }
        }
        table[k,p]=mean(colMeans(results))
      }
    }
    #if((max(table[k,])-r0)>0.01){
      index=which(table[k,]==max(na.omit(table[k,])))[1]
      base=c(base,(index+7))
      rest=rest[-which(rest==index)]
    #}
    table[18,k]=max(na.omit(table[k,]))
    table[19,k]=index
    table[20,k]=min(na.omit(table[k,]))
  }
}
write.table(table,paste("EUR_miss",miss,"_model_selection_svm_table_relationship_degree2_err",err,".",fam,"fam.txt",sep=''),quote=F,row.names = F,col.names = F)

