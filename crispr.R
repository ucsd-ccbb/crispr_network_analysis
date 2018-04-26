C1<-read.table("../T0_1/T0_1_GeneCounts.tsv",sep="\t",header=FALSE)
C2<-read.table("../T0_2/T0_2_GeneCounts.tsv",sep="\t",header=FALSE)
C3<-read.table("../T0_3/T0_3_GeneCounts.tsv",sep="\t",header=FALSE)
T1<-read.table("../T21_1/T21_1_GeneCounts.tsv",sep="\t",header=FALSE)
T2<-read.table("../T21_2/T21_2_GeneCounts.tsv",sep="\t",header=FALSE)
T3<-read.table("../T21_3/T21_3_GeneCounts.tsv",sep="\t",header=FALSE)
X<-data.frame(C1[,2],C2[,2],C3[,2],T1[,2],T2[,2],T3[,2])

data<-data.matrix(X)
genes<-(as.character(C1[,1]))
genes[which(genes=="5-Mar")]<-"Mar5"
genes[which(genes=="5-Sep")]<-"Sep5" #takes care of excel problem

rownames(data)<-genes #unique IDs
ngenes<-length(genes)

total<-apply(data,2,sum)
cdata<-t(t(data)/total)

mint<-min(total)
ndata<-cdata*mint #roughly normalized data, only for detection call
detected<-apply(ndata[,1:3],1,mean)>32 | apply(ndata[,4:6],1,mean)>32

gdata<-data[detected,] #only detected genes
gdata[gdata==0]<-1 #smallest measurable amount other than 0
total<-apply(gdata,2,sum)
cdata<-t(t(gdata)/total)
xdata<-(cdata[,4:6]/cdata[,1:3])

medx<-apply(xdata,2,median)
xdata<-t(t(xdata)/medx) #adjusted to median 1

lxdata<-log2(xdata)
mlxdata<-apply(lxdata,1,mean)

#y<-sqrt(xdata) #I will do a ONE TAIL TEST for suppressors (KO causes growth)
#y<-1/sqrt(xdata) #I will do a ONE TAIL TEST for enhancers (KO causes decay)
biolog<-function(x) ifelse(x>1,log(x)+1,x)



#*********************CRUCIAL***********************
y<-biolog(1/xdata) #for enhancers
#y<-biolog(xdata) #for suppressors
#looking for large y
v<-c(1,1,1)/sqrt(3)
s<-function(x) { #define s-statistic
   x111<-sum(v*x)*v #along 111
   xp<-x-x111 #perpendicular to 111
   u<-sum(x111*x111)-sum(xp*xp)
   return(u)
}
ys<-apply(y,1,s)

#permutation null:
ng<-nrow(y)
nmc<-1000 #this many samples
sampleys<-matrix(0,nrow=ng,ncol=nmc)
sampley<-y
for (i in 1:nmc) {
sampley[,2]<-sample(sampley[,2])
sampley[,3]<-sample(sampley[,3])
sampleys[,i]<-apply(sampley,1,s)
}

e0<-ecdf(-sampleys)
e<-ecdf(-ys)
pperm<-e0(-ys) #permutation p-value
pperm[pperm==0]<-1/(ng*nmc)
fdr<-pperm/e(-ys) #empirical FDR
o<-order(pperm)

res<-data.frame(lxdata,mlxdata,pperm,fdr)
colnames(res)<-c("alphaT_1","alphaT_2","alphaT_3","mean alphaT","permutation p","permutation FDR")
write.table(res[o,],file="2D_enhancers_p_fdr.txt",sep="\t",row.names=TRUE,col.names=NA,quote=FALSE)
