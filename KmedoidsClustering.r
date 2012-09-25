#####################Parameters##################################
input.monitors.file="c:/Users/tkhachik/Desktop/mon.csv"
load.file="c:/Users/tkhachik/Desktop/load.csv"
output.file="c:/Users/tkhachik/Desktop/clusters.csv"
######################Initialization#############################
require(cluster)
threshold=0.975
monitors<- read.csv(file=input.monitors.file,head=TRUE,sep=",")
loading<- read.csv(file=load.file,head=TRUE,sep=",")
monitors.names<- names(monitors)
loading.names<- names(loading)
tstmploc<- which(names(monitors)=="Timestamp")
if (length(tstmploc)>0) {monitors<- monitors[-tstmploc];monitors.names<- monitors.names[-tstmploc];}
nobsm<- min(sapply(monitors,length))
nobsl<- min(sapply(loading,length))
nobs<- min(nobsm,nobsl)
##############Removing sigma=0 columns###########################
msigma<- mat.or.vec(length(monitors),0)
lsigma<- mat.or.vec(length(loading),0)
for(i in 1:(length(monitors)+length(loading))) 
{if (length(monitors)>=i) {msigma[i]=sd(monitors[[i]][seq(1,nobs)],na.rm=TRUE)}
 if (length(monitors)<i) {lsigma[i-length(monitors)]=sd(loading[[i-length(monitors)]][seq(1,nobs)],na.rm=TRUE)}
}
if (length(which(msigma==0))>0) {monitors0<- monitors[which(msigma==0)]; monitors<- monitors[-which(msigma==0)];monitors.names<- monitors.names[-which(msigma==0)]}
if (length(which(lsigma==0))>0) {loading0<- loading[which(msigma==0)];loading<- loading[-which(lsigma==0)];loading.names<- loading.names[-which(lsigma==0)]}
##############Combining monitors and load########################
monitors.std<- mat.or.vec(length(monitors)+length(loading),nobs)
param.names<- append(monitors.names,loading.names)
for(i in 1:(length(monitors)+length(loading))) 
{if (length(monitors)>=i) {monitors.std[i,]<- as.vector(monitors[[i]][seq(1,nobs)])}
 if (length(monitors)<i) {monitors.std[i,]<- as.vector(loading[[i-length(monitors)]][seq(1,nobs)])}
}
##############Calculating coreeletion matrix and PCA#############
monitors.std<- t(monitors.std)
moncor<- cor(monitors.std)
mpcc<- princomp(moncor,cor=TRUE)
var.mpcc<- mpcc$sdev^2
cum.mpcc<- cumsum(var.mpcc/sum(var.mpcc))
ncluster.cor=max(which(cum.mpcc<threshold))
#clara.cor<-  clara(moncor, ncluster.cor)
##plot(clara.cor,color=TRUE, shade=TRUE, lines=1,labels=5)
#################################################################
n<- length(param.names)
CorrMat<- mat.or.vec(n,n+1)
CorrMat[,1] <- param.names
for(i in 1:n){
  CorrMat[,1+i] <- moncor[,i]
}
M <- moncor
kk <- ncluster.cor
clarax <-  clara(M, kk)
cl_clusinfo<-clarax$clusinfo
plot(clarax,color=TRUE, shade=TRUE, lines=1,labels=5)
caa <- clarax$clust 
length(clarax)
cl_clustering<-clarax$clustering
cl_medoids <- clarax$medoids
cl_diss <-clarax$diss
dissMat<-as.matrix(cl_diss)
clustersSize<- cl_clusinfo[,1]
ClusterVector<- vector(mode = "list", length =kk)
for(i in 1:kk ){
  ClusterVector[i]<-list(which(cl_clustering==i))
}
########################################################
ClusterDistance <- function(XMat,medoid,clusterInd){
  nn<-length(clusterInd)
  distvector<-as.vector(nn)
  for(i  in 1:nn){
    aa  <- as.matrix(rbind(XMat[medoid,],XMat[clusterInd[i],]))
    distvector[i] <- dist(aa)
  }
  return (distvector)
}

DistanceMatrix <- vector(mode = "list", length =kk)

for(j in 1:kk ){
  DistanceMatrix[[j]] <- list(ClusterDistance(M,clarax$i.med[j],ClusterVector[[j]]))
}

clusterInd_Distance <- sapply(1:kk,function(i) data.frame(ClusterVector[[i]],DistanceMatrix[[i]]))

##############################

OrderedSingleCluster <- function(CorrMat,clusterInd_Distance,N_Clust){
  #for(j in 1:kk )
  j<-N_Clust
  n<- clustersSize[j]
  #  print(n)
  if(n != 1){
    clmat <- mat.or.vec(nr=3,nc=n)
    (ii<- order(clusterInd_Distance[[2*j]],clusterInd_Distance[[2*j-1]]))
    aa <- rbind(clusterInd_Distance[[2*j-1]],clusterInd_Distance[[2*j]])[,ii]
    bb<- as.matrix(CorrMat[aa[1,],1])
    comment(clmat[1,])<- "cluster"
    comment(clmat[2,])<- "index"
    comment(clmat[3,])<- "distance"
    
    clmat[1,]<- t(bb)
    clmat[2,]<- aa[1,]
    clmat[3,]<- aa[2,]
  }
  else
  {
    clmat <- mat.or.vec(nr=3,nc=n)
    (ii<- order(clusterInd_Distance[[2*j]],clusterInd_Distance[[2*j-1]]))
    aa <- rbind(clusterInd_Distance[[2*j-1]],clusterInd_Distance[[2*j]])
    bb<- CorrMat[aa[1,],1]
    
    clmat[1]<- bb
    clmat[2]<- aa[1,]
    clmat[3]<- aa[2,]
    
  }
  return (clmat)
}

############################################################

OrderedClusters <- function(CorrMat,clusterInd_Distance){
  for(j in 1:kk ){
    
    
    n_cl<- clustersSize[j]
    
    if(n_cl != 1){
      clmat <- mat.or.vec(nr=3,nc=n_cl)
      (ii<- order(clusterInd_Distance[[2*j]],clusterInd_Distance[[2*j-1]]))
      aa <- rbind(clusterInd_Distance[[2*j-1]],clusterInd_Distance[[2*j]])[,ii]
      bb<- as.matrix(CorrMat[aa[1,],1])
      
      clmat[1,]<- t(bb)
      clmat[2,]<- aa[1,]
      clmat[3,]<- aa[2,]
    }
    else
    {
      clmat <- mat.or.vec(nr=3,nc=n_cl)
      (ii<- order(clusterInd_Distance[[2*j]],clusterInd_Distance[[2*j-1]]))
      aa <- rbind(clusterInd_Distance[[2*j-1]],clusterInd_Distance[[2*j]])
      bb<- as.matrix(CorrMat[aa[1,],1])
      
      clmat[1]<- (bb)
      clmat[2]<- aa[1,]
      clmat[3]<- aa[2,]
    }
    
    
    
    write.table((clmat), file = output.file, append = TRUE, quote = FALSE, sep = ", ",
                eol = "\n", na = "NA", dec = ".", row.names = c(paste("cluster ",j)  ,"index","dist"),
                col.names = FALSE, qmethod = c( "escape"),
                fileEncoding = "")
    
    
  }
  
  
  
  return (clmat)
  
}
constmon <- t(names(monitors0))
write.table(constmon, file = output.file,append = TRUE, quote = FALSE, sep = ", ",
            eol = "\n", na = "NA", dec = ".", row.names = c("constant metrics " ),
            col.names = FALSE, qmethod = c( "escape"),
            fileEncoding = "")
vvv <- OrderedClusters(CorrMat,clusterInd_Distance)

