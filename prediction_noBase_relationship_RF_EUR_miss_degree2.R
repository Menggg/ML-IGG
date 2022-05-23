#!/usr/bin/env Rscript
err=as.numeric(commandArgs(trailingOnly=TRUE)[1])#genotyping error rate
fam=as.numeric(commandArgs(trailingOnly=TRUE)[2])#number of simlated families
miss=as.numeric(commandArgs(trailingOnly=TRUE)[3])#missing SNPs rate
library(data.table)
library(ranger)
library(aricode)
data=as.data.frame(fread(paste("/home/meng/Ped-sim/merge_EUR_miss",miss,"_sim_table_multi_measurements_err0_1_",fam,".txt",sep=''),head=T))
data_pre=as.data.frame(fread(paste("EUR_miss",miss,"_final_degree_err",err,"_prediction.txt",sep=''),head=F))
index=which(data_pre[,2]=="0.125")
data=data[index,]
data1=as.data.frame(fread(paste("/home/meng/Ped-sim/merge_EUR_miss",miss,"_sim_table_multi_measurements_err",err,"_1_",fam,".txt",sep=''),head=T))
        data1=data1[index,]
table=matrix(NA,20,ncol(data)-4)
rep=1
cycles=10
n_method=0
r0=0
base=c(3)
rest=c(5:ncol(data))-4
for(k in 1:17){
  if(length(rest)>0){
    for(p in rest){
      # x=data.frame(data[,-c(1,2,3,(p+6))])
      x=data.frame(data[,c(base,(p+4))])
      x$pedigree.relationship=factor(x$pedigree.relationship)
      for(j in err){
        results=matrix(NA,nrow=cycles, ncol=rep)
        # x1=data.frame(data1[,-c(1,2,3,(p+6))])
        x1=data.frame(data1[,c(base,(p+4))])
        x1$pedigree.relationship=factor(x1$pedigree.relationship)
        for(i in 1:rep){
          s=data.frame(matrix(sample(1:(floor(nrow(x)/cycles)*cycles)), ncol=cycles))
          for(r in 1:cycles){
            test=s[,r]
	  index=c(which(x[-test,]$pedigree.relationship=="GP"),which(x[-test,]$pedigree.relationship=="HS"),which(x[-test,]$pedigree.relationship=="Uncle"))
	  temp=x[-test,][index,]
	  temp$pedigree.relationship=factor(temp$pedigree.relationship)
	  rf <- ranger(pedigree.relationship ~ ., data = temp)
            myY_P = predict(rf,x1[test,])$predictions
#results[r,i]=ARI(temp1,as.character(temp2))
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
write.table(table,paste("/home/meng/Ped-sim/EUR_miss",miss,"_model_selection_RF_table_relationship_degree2_err",err,".",fam,"fam.txt",sep=''),quote=F,row.names = F,col.names = F)

