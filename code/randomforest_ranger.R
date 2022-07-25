library(ranger)
x=read.csv("data.csv",head=T)
x$pedigree.relationship=factor(x$pedigree.relationship)
x1=read.csv("data1.csv",head=T)
x1$pedigree.relationship=factor(x1$pedigree.relationship)
rf <- ranger(pedigree.relationship ~ ., data = x)
            myY_P = predict(rf,x1)$predictions

myY_P=as.data.frame(myY_P)
colnames(myY_P)[1]="predicted_relationship_type"
write.csv(myY_P,"RF_pred.csv",quote=F,row.names = F)

