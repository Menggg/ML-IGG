#!/usr/bin/env Rscript

j=as.numeric(commandArgs(trailingOnly=TRUE)[1])
name1=commandArgs(trailingOnly=TRUE)[2]
name2=commandArgs(trailingOnly=TRUE)[3]
miss=as.numeric(commandArgs(trailingOnly=TRUE)[4])
library(data.table)


data=as.data.frame(fread(name1,head=T))
index=c(which(data[,3]=="GP"),which(data[,3]=="HS"),which(data[,3]=="Uncle"))
data=data[index,]
rep=1
cycles=10
#base=c(3,5,7,11,14,18:20)
base=c(3,c(10,13,15,12,16,9,14,4,17,11,1:3,7)+4)
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
    #rf <- randomForest(pedigree.relationship ~ ., data = x[-test,], importance = TRUE, proximity = TRUE)
            #myY_P = predict(rf,x1[test,])
    #results[r,1]=length(which(x1[test,1]==as.character(myY_P)))/nrow(x1)

    }

#write.table(results,paste("EUR_miss",miss,"_final_RF_table_relationship_degree2_err",j,".100fam.txt",sep=''),quote=F,row.names = F,col.names = F)
