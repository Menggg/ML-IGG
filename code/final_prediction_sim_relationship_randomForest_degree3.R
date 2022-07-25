#!/usr/bin/env Rscript

j=as.numeric(commandArgs(trailingOnly=TRUE)[1])
name1=commandArgs(trailingOnly=TRUE)[2]
name2=commandArgs(trailingOnly=TRUE)[3]
miss=as.numeric(commandArgs(trailingOnly=TRUE)[4])
library(data.table)


data=as.data.frame(fread(name1,head=T))
index=c(which(data[,3]=="GGP"),which(data[,3]=="Great-uncle"),which(data[,3]=="First-cousin"),which(data[,3]=="Half-uncle"))
data=data[index,]
rep=1
cycles=10
#base=c(3,5,7,11,14,18:20)
base=c(3,c(10,13,15,16,9,11,12,1,2,3,14)+4)
#base=c(3,5:21)
x=data.frame(data[,c(base)])
x$pedigree.relationship=factor(x$pedigree.relationship)

results=matrix(NA,nrow=cycles, ncol=rep)
  data1=as.data.frame(fread(paste(name2,sep=''),head=T))
data1=data1[index,]
  x1=data.frame(data1[,c(base)])
  x1$pedigree.relationship=factor(x1$pedigree.relationship)
    s=data.frame(matrix(sample(1:(floor(nrow(x)/cycles)*cycles)), ncol=cycles))
    for(r in 1:cycles){
      test=s[,r]
    system(paste("mkdir -p relationship/relationship_miss_",miss,"_err_",j,"_rep_",r,sep=''))
system(paste("cp randomforest_ranger.R relationship/relationship_miss_",miss,"_err_",j,"_rep_",r,sep=''))
setwd(paste("relationship/relationship_miss_",miss,"_err_",j,"_rep_",r,sep=''))
  write.csv(x1[test,],"data1.csv",quote=F,row.names = F)
     write.csv(x[-test,],"data.csv",quote=F,row.names = F)
system("/usr/bin/Rscript randomforest_ranger.R")
    }

#write.table(results,paste("EUR_miss",miss,"_final_RF_table_relationship_degree3_err",j,".100fam.txt",sep=''),quote=F,row.names = F,col.names = F)
