library(dplyr)
library(VGAM)
args <- commandArgs(trailingOnly=TRUE)
filename <- args[1]
outfilename <- args[2]
fsize <- args[3]
betaNLL <- function(params,data){
  a<-params[1]
  b<-params[2]
  #negative log likelihood for beta
  return(-sum(dbeta(data,shape1=a, shape2=b, log=TRUE)))
}

calculate.bb.pvalue<-function(x2){
  a1<-x2$Coverage*x2$Max.Non.ref.Allele.Frequency
  b1<-x2$Coverage
  m<-mean(a1/b1)
  v<-var(a1/b1)
  a0<-m*(m * (1-m) / v-1 )
  b0<-(1-m)*(m * (1-m) / v-1 )
  params0=c(a0,b0)
  fit <- nlm(betaNLL,params0,a0/b0)
  a<-fit$estimate[1]
  b<-fit$estimate[2]
  print(a)
  print(b)
  pval<-1-pbetabinom.ab(a1,b1,a,b)
  return(pval)
}

x <- read.table(filename,sep="\t",header=TRUE)
x2 <- x %>% filter(Consensus.group.size==fsize & Coverage >= 100 & !Name=="")
pval<-calculate.bb.pvalue(x2)
fdr<-p.adjust(pval,method="fdr")
results<-cbind(x2,pval,fdr)
write.table(results[fdr<0.05,],outfilename,sep='\t',row.names=FALSE,quote=FALSE)

