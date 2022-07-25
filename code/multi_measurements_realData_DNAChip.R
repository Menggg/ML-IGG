#convert vcf file to SNP major text file .traw
system("plink1.9 --vcf 50ng.vcf --export A-transpose --snps-only --autosome --out 50ng")
system("plink1.9 --vcf 500pg.vcf --export A-transpose --snps-only --autosome --out 500pg")
system("plink1.9 --vcf 100pg.vcf --export A-transpose --snps-only --autosome --out 100pg")

GSA=read.table("GSA-24v3-0_A1.csv.autosomes.chr-pos.tsv",head=F)
GSA[,3]=paste(GSA[,1],GSA[,2],sep=':')
colnames(GSA)[3]="X.C.M.x"

MM <- function(N1,N2,data){
  a=na.omit(data[,c(N1,N2)])
  n1=length(intersect(which(a[,1]==1),which(a[,2]==1)))
  n2=length(intersect(which(a[,1]*a[,2]==0),which(a[,1]+a[,2]==2)))
  n3=length(which(a[,1]==1))
  n4=length(which(a[,2]==1))
  n5=min(c(n3,n4))
  kin=(n1-2*n2)/(n3+n4)
  kin0=0.5-(0.25*(n3+n4-2*n1)+n2)/n5
  j1=length(intersect(which(a[,1]==0),which(a[,2]==0)))
  j2=length(intersect(which(a[,1]==1),which(a[,2]==0)))
  j3=length(intersect(which(a[,1]==2),which(a[,2]==0)))
  j4=length(intersect(which(a[,1]==0),which(a[,2]==1)))
  j5=length(intersect(which(a[,1]==1),which(a[,2]==1)))
  j6=length(intersect(which(a[,1]==2),which(a[,2]==1)))
  j7=length(intersect(which(a[,1]==0),which(a[,2]==2)))
  j8=length(intersect(which(a[,1]==1),which(a[,2]==2)))
  j9=length(intersect(which(a[,1]==2),which(a[,2]==2)))
  
  ibs0=(j3+j7)/nrow(a)
  ibs1=(j2+j4+j6+j8)/nrow(a)
  ibs2=(j1+j5+j9)/nrow(a)
  ibs02=ibs0+ibs2
  ibs12=ibs1+ibs2
  ibs01=ibs0+ibs1
  
  results=c(kin,kin0,ibs0,ibs1,ibs2,ibs02,ibs12,ibs01,j1/nrow(a),j2/nrow(a),j3/nrow(a),j4/nrow(a),j5/nrow(a),j6/nrow(a),j7/nrow(a),j8/nrow(a),j9/nrow(a))
  return(results)
}

r=matrix(NA,1000,21)
k=1
ID=c("500pg","100pg")
relationship=c("UN","GP","UN","UN","GP","PO","GP","GP","UN","UN","GP","PO","GP","PO","UN","first-cousin","PO","FS","UN","UN","UN","PO","PO","UN","UN","uncle-nephew","first-cousin","PO")
degree=c(0,0.125,0,0,0.125,0.25,0.125,0.125,0,0,0.125,0.25,0.125,0.25,0,0.0625,0.25,0.25,0,0,0,0.25,0.25,0,0,0.125,0.0625,0.25)
data1=read.table("50ng.traw",head=T)
data1[,3]=paste(data1[,2],data1[,4],sep=':')
data1=merge(GSA,data1,by.x="X.C.M.x",by.y="X.C.M.x")
#data2=merge(GSA,data2,by.x="SNP",by.y="SNP")
data1=data1[,-(1:8)]
for(j in 1:length(ID)){
	q=1
	data2=read.table(paste(ID[j],".traw",sep=''),head=T)
	data2[,3]=paste(data2[,2],data2[,4],sep=':')
data2=merge(GSA,data2,by.x="X.C.M.x",by.y="X.C.M.x")
	data2=data2[,-(1:8)]
for(i in 1:(ncol(data1)-1)){
                for(p in (i+1):ncol(data1)){
		data=cbind(data1[,i],data2[,p])
r[k,5:21]=MM(1,2,data)
r[k,3:4]=c(ID[j],degree[q])
r[k,1:2]=c(paste(colnames(data1)[i],sep=''),paste(colnames(data2)[p],sep=''))
k=k+1
q=q+1
}
}
}
r=na.omit(r)
colnames(r)=c("individual_1","individual_2","pedigree-relationship","pedigree-degree","kin","kin0","ibs0","ibs1","ibs2","ibs02","ibs12","ibs01","j1","j2","j3","j4","j5","j6","j7","j8","j9")
write.table(r,paste("real_chip_table_multi_measurements_GSA.txt",sep=''),quote=F,row.names = F,col.names = T)
