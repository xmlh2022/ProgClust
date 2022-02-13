# R functions used
speccluster<-function(data,method='leiden'){
  
  if(method=='leiden')
    clust<-leiden(tem)
  if(method=='spec'){
    
    te<-data/(apply(data,1,max)+1e-10)
    te_n <- t(t(te)/(sqrt(apply(te^2,2,sum))+1e-10))
    tem<-as.matrix(t(te_n))%*%as.matrix(te_n)
    
    #mymethod
    
    for(i in 1:dim(tem)[1]){
      tem[i,i]<-0
    }

    #clust<-speccalt(tem,maxk=min(dim(tem)[1],7))
    means<-c(sum(tem)/(dim(tem)[1]^2-dim(tem)[1]))
    clusts<-c(rep(1,dim(tem)[1]))
    clustb<-rep(means,dim(tem)[1])
    for ( k in 2:(min(10,dim(tem)[1]-1))){
      clust<-speccalt(tem,k=k)
      clusts<-c(clusts,clust)
      add_means<-c()
      for (i in 1:k){
        smean<-sum(tem[which(clust==i),which(clust==i)])
        num<-(length(which(clust==i)))^2-length(which(clust==i))
        mean_<-smean/(num+1e-10)
        add_means<-c(add_means,(mean_*length(which(clust==i))-sum(clustb[which(clust==i)]))/sum(clustb[which(clust==i)]))
        if(mean_!=0)
          clustb[which(clust==i)]<-mean_
      }

      means<-c(means,max(add_means)+sum(add_means[which(add_means<0)]*(table(clust)[which(add_means<0)]/table(clust)[which(add_means==max(add_means))][1])))
      if(means[k]<=k_crit){
        ek<-k-1
        clust<-clusts[(1+(ek-1)*dim(tem)[1]):(ek*dim(tem)[1])]
        return (clust)
        break
      }
    }
    clust<-clusts[(1+(k-1)*dim(tem)[1]):(k*dim(tem)[1])]
    return (clust)
  }
  
  if(method=='apclust'){
    #data<-data+1
    #te<-data
    te<-data/(apply(data,1,max)+1e-10)
    te_n <- t(t(te)/(sqrt(apply(te^2,2,sum))+1e-10))
    tem<-as.matrix(t(te_n))%*%as.matrix(te_n)
    
    clust<-apcluster(tem)
    clust<-clust@idx
  }
  return (list(clust,tem))
}

kmeanscl<-function(m,verb){
  pca_Fano<-calcul.pca(t(m),50)
  gapStatistic<-clusGap(pca_Fano$pca,FUN=kmeans,K.max=min(dim(pca_Fano$pca)[1]-1,K.max),spaceH0="original",d.power = 2,iter.max=50,verbose = verb)
  k<-with(gapStatistic,maxSE(Tab[,"gap"],Tab[,"SE.sim"],method="Tibs2001SEmax",SE.factor=1))
  # k=4
  seeds<-sample(1:1000,20,replace=FALSE)
  errors<-c()
  for (seed in seeds){
    set.seed(seed)
    MaxFanoClusters_pca<-kmeans(pca_Fano$pca,k,iter.max=50,algorithm="MacQueen")
    errors<-c(errors, MaxFanoClusters_pca$tot.withinss)
  }
  seed<-seeds[which.min(errors)]
  set.seed(seed)
  MaxFanoClusters_pca<-kmeans(pca_Fano$pca,k,iter.max=50,algorithm="MacQueen")
  
  P_F<-MaxFanoClusters_pca$cluster #Fano clustering result
  return (P_F)
}

filt <- function(tem,m_filt,rec,crit){
  rr_i<-c()
  for(i in 1:length(table(rec))){
    rma<-tem[which(rec==names(table(rec))[i]),which(rec==names(table(rec))[i])]
    rr_i<-c(rr_i,(sum(rma)-table(rec)[i])/(table(rec)[i]^2-table(rec)[i]+1e-6))
  }
  if(length(table(rec))==1){
    unselected<-c()
    return (unselected)
  }
  rr_o<-c()
  for(i in 1:length(table(rec))){
    rma<-tem[which(rec==names(table(rec))[i]),which(rec!=names(table(rec))[i])]
    rr_o<-c(rr_o,mean(rma))
  }
  
  
  unselected<-c()
  for(i in which(length(table(rec))==1)){
    if(rr_o[i]<=o_crit)
      unselected<-c(unselected,which(rec==names(table(rec)[i])))
  }
  for(i in which(rr_i/(rr_o+1e-10)>=crit & length(table(rec))!=1)){
    unselected<-c(unselected,which(rec==names(table(rec)[i])))
  }
  return (unselected)
}

filt2 <- function(data){
  if (dim(data)[2]<=2)
    return (1:dim(data)[2])
  cellsum<-apply(data,2,sum)
  scellsum<-sort(cellsum,decreasing = F)
  flag<-filt_c
  select<-which(cellsum<=flag)
  select
}

Fanocluster<-function(data,label,flag=FALSE,deep,Gene_be,verb){
  if(is.null(dim(data))[1])
    return(label)
  if(dim(data)[2]==2)
    return(c(label,label))
  if(verb)
    print('num of cells')
    print(dim(data)[2])
  Fano=apply(data,1,function(x) var(x)/(mean(x)+1e-10))
  FanoDecreasing<-sort(Fano,decreasing=TRUE)
  maxFano=names(FanoDecreasing[1:1000])
  Gene_be = unique(c(Gene_be,maxFano))
  m<-data[Gene_be,]
  if(verb)
    print('dim of fano-based matrix')
    print(dim(m))
  
  if(dim(data)[2]<=50){
    P_F_n<-speccluster(m,'spec')
  }
  else
    P_F_n<-kmeanscl(m,verb)
  if(verb)
    print('fano_based result')
    print(table(P_F_n))
  if(length(table(P_F_n))==1 & !flag ){
    return(rep(label,length(P_F_n)))
  }
  flag=FALSE
  label_next<-label
  P_F <- rep(0,dim(data)[2])
  P_F_p<-P_F_n
  for(i in 1:length(table(P_F_n))){
    if(length(which(P_F_n==i))>2){
      if(length(which(P_F_n==i))>50)
        selected<-Ginicluster(data[,which(P_F_n==i)],verb=verb)
      else
        selected<-rep(1,length(which(P_F_n==i)))
      if(sum(selected)!=length(selected))
        flag = TRUE
      else{
        if(length(table(P_F_n))==1)
          return (rep(label,length(P_F_n)))
      }
      if(deep>1){
        labeli<-Fanocluster(data[,which(P_F_n==i)][,which(selected==1)],label_next,flag,deep-1,Gene_be,verb=verb)
        P_F[which(P_F_n==i)][which(selected==1)]=labeli
      }
      else{
        labeli<-label_next
        P_F[which(P_F_n==i)][which(selected==1)]=labeli
      }
      P_F_p[which(P_F_n==i)][which(selected==0)]<-0
    }
    else{
      labeli<-label_next
      P_F[which(P_F_n==i)]=labeli
    }
    label_next<-max(labeli)+1
  }
  # if(length(table(P_F_p[which(P_F_p!=0)]))>1){
  #   ExprM.norCounts<-m
  #   for (i in 1:dim(m)[1]){
  #     ExprM.norCounts[i,]<-.normalize_by_gene(m[i,])
  #   }
    # dir.create(paste("F:/LH/my_method/result/myresult/real_mouse/tree",label,'_',deep,sep=""), showWarnings = FALSE)
    # save(m,file=paste("F:/LH/my_method/result/myresult/real_mouse/tree",label,'_',deep,'/m.RData',sep=""))
    # save(P_F_p,file=paste("F:/LH/my_method/result/myresult/real_mouse/tree",label,'_',deep,'/P_F_p.RData',sep=""))
    # DE(m,P_F_p,paste("F:/LH/my_method/result/myresult/real_mouse/tree",label,'_',deep,"/",sep=""))
    # difgenep(ExprM.norCounts,P_F_p,paste("F:/LH/my_method/result/myresult/real_mouse/tree",label,'_',deep,"/",sep=""))
  # }
  return (P_F)
}

Ginicluster<-function(data,verb){
  ExprM.RawCounts.filter<-data
  Genelist.top_pvalue<-ginfitting(ExprM.RawCounts.filter)
  if(verb)
    print('num of gini genes')
    print(length(Genelist.top_pvalue))
  if(length(Genelist.top_pvalue)==0)
    return (rep(1,dim(m)[2]))
  else{
    m <- ExprM.RawCounts.filter[Genelist.top_pvalue,]
    m_filt <- filt2(m)
    lst<-speccluster(m,method='apclust')
    P_G<-lst[[1]]
    tem<-lst[[2]]
    if (verb)
      print('num of gini genes')
      print(table(P_G))
    unselected<-filt(tem,m_filt,P_G,crit)
    selected<-rep(1,dim(m)[2])
    # selected[m_filt][unselected]<-0
    selected[unselected]<-0
    selected[m_filt]<-1
    if(sum(selected)==0){
      selected[which(P_G==names(table(P_G)[which(table(P_G)==max(table(P_G)))])[1])]<-1
    }
    if(verb){
      print('gini_based result')
      print(table(selected))
    }
    return(selected)
  }
}

igclust<-function(data,deep=2,verb=FALSE){
  # te<-data/(apply(data,1,max)+1e-10)
  # te_n <- t(t(te)/(sqrt(apply(te^2,2,sum))+1e-10))
  # tem2<-as.matrix(t(te_n))%*%as.matrix(te_n)
  # for(i in 1:dim(tem2)[1]){
  #   tem2[i,i]<-0
  # }
  finalCluster<-rep(0,dim(data)[2])
  i<-1
  while(sum(which(finalCluster==0))!=0){
    print(paste('clustering tree ',i,'...'))
    if(sum(which(finalCluster==0))>1){
      finalCluster[which(finalCluster==0)]<-Fanocluster(data[,which(finalCluster==0)],i,flag=TRUE,deep=deep,c(),verb=verb)
      
      i<-max(finalCluster)+1
    }
    else
      finalCluster[which(finalCluster==0)]<-i
  }
  # ExprM.norCounts<-data
  # for (i in 1:dim(data)[1]){
  #   ExprM.norCounts[i,]<-.normalize_by_gene(data[i,])
  # }
  # dir.create("F:/LH/my_method/result/myresult/real_mouse/final", showWarnings = FALSE)
  # save(data,file="F:/LH/my_method/result/myresult/real_mouse/final/m.RData")
  # save(finalCluster,file="F:/LH/my_method/result/myresult/real_mouse/final/P_F_p.RData")
  # DE(data,finalCluster,"F:/LH/my_method/result/myresult/real_mouse/final/")
  # difgenep(ExprM.norCounts,finalCluster,"F:/LH/my_method/result/myresult/real_mouse/final/")
  return(finalCluster)
} 


H<-function(x,w){
  t<-c()
  for(i in 1:length(table(x))){
    t<-c(t,sum(w[which(x==names(table(x)[i]))]))
  }
  hx<-0
  for(i in 1:length(t)){
    p<-t[i]/sum(t)
    hx<-hx-p*log(p)
  }
  hx
}

MI <- function(x,y,w){
  t<-matrix(0,nrow=length(table(x)),ncol=length(table(y)))
  for(i in 1:length(table(x)))
    for(j in 1:length(table(y)))
      t[i,j]<-sum(w[which(x==names(table(x)[i]) & y==names(table(y)[j]))])
    hxy<-0
    for(i in 1:length(table(x)))
      for(j in 1:length(table(y))){
        if(t[i,j]!=0){
          pij<-t[i,j]/sum(t)
          pi<-sum(t[i,])/sum(t)
          pj<-sum(t[,j])/sum(t)
          hxy<-hxy+pij*log(pij/(pi*pj))
        }
      }
    hxy      
}

wnmi<-function(finalCluster,label){
  ww<-1/as.numeric(table(label))
  w<-rep(1,length(label))
  for( i in 1:length(table(label))){
    w[which(label==as.numeric(names(table(label))[i]))]<-ww[i]
  }
  return(MI(finalCluster,label,w)*2/(H(finalCluster,w)+H(label,w)))
}

wnmi_rare<-function(finalCluster,label,rare_label){
  nas<-c(1:length(table(finalCluster)))
  for(kk in names(table(finalCluster[which(label==rare_label)]))){
    nas[as.numeric(kk)]<-0
  }
  finalCluster_r<-finalCluster
  for(ii in nas){
    if(ii!=0){
      finalCluster_r[finalCluster_r==ii]<-0
    }
  }
  label_r<-label
  label_r[label_r!=rare_label]<-0
  
  nmi_r<-wnmi(finalCluster_r,label_r)
  return(nmi_r)
}


wpurity_rare<-function(finalCluster,label,rare_label){
  nas<-c(1:length(table(finalCluster)))
  for(kk in names(table(finalCluster[which(label==rare_label)]))){
    nas[as.numeric(kk)]<-0
  }
  finalCluster_r<-finalCluster
  for(ii in nas){
    if(ii!=0){
      finalCluster_r[finalCluster_r==ii]<-0
    }
  }
  label_r<-label
  label_r[label_r!=rare_label]<-0
  pu_r<-wpurity(finalCluster_r,label_r)
  return(pu_r)
}

wpurity<-function(finalCluster,label){
  ta<-table(finalCluster,label)
  ww<-1/as.numeric(table(label))
  w<-rep(1,length(label))
  pi<-c()
  for(i in 1:dim(ta)[1]){
    wp<-c()
    for(j in 1:dim(ta)[2])
      wp<-c(wp,ww[j]*ta[i,j])
    pi<-c(pi,max(wp))
  }
  return (sum(pi)/length(table(label)))
}



.normalize_by_gene<-function(v) {
  m <- log(1+v)
  m<-m-mean(m)
  m<-m/sd(m)
  m
}


difgenep<-function(ExprM.RawCounts.filter,finalCluster,de_dir){
  Genes<-c()
  for(i in 1:100){
    file = paste(de_dir,i,"_lrTest_Sig.csv",sep='')
    if (file_test("-f", file)) {
      x<-read.csv(file, header=T, na.strings=c("NA"))
      if (length(x$Gene)>=max(2,round(15/length(table(finalCluster[which(finalCluster!=0)])))))
        Genes<-c(Genes,x$Gene[1:max(2,round(15/length(table(finalCluster[which(finalCluster!=0)]))))])
      else
        Genes<-c(Genes,x$Gene)
    }
  }
  Genes<-unique(Genes)
  # for(i in 1:dim(table(finalCluster[which(finalCluster!=0)]))){
  #   if(i == 1){
  #     testt<-apply(ExprM.RawCounts.filter[,which(finalCluster==i)],1,mean)
  #   }
  #   else{
  #     testt<-cbind(testt,apply(ExprM.RawCounts.filter[,which(finalCluster==i)],1,mean))
  #   }
  # }
  
  for(i in 1:dim(table(finalCluster))){
    if(i == 1){
      testt<-ExprM.RawCounts.filter[,which(finalCluster==i)]
    }
    else{
      testt<-cbind(testt,ExprM.RawCounts.filter[,which(finalCluster==i)])
    }
  }
  
  print(length(Genes))
  colnames(testt)<-NULL
  testt<-testt[Genes,]

  # for(i in 1:(dim(testt)[2]-2)){
  #   dis<-c()
  #   for(j in (i+1):dim(testt)[2]){
  #     dis<-c(dis,sum((te[,i]*te[,j])/((sum(te[,i]^2))^0.5*(sum(te[,j]^2))^0.5)))
  #   }
  #   to<-which(dis==min(dis))
  #   to2<-testt[,i+1]
  #   testt[,i+1]<-testt[,i+to]
  #   testt[,i+to]<-to2
  # }
  # ce<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
  # ce<-c(1,8,2,3,13,11,4,5,10,12,7,14,6,9)
  # GENE<-c('Cdh5','Pecam1','Ednrb','Gm13889','Ins1','Ins2','G6pc2','Igf1r','Ptprn','Mt1','Mt2','Clu','Krt8','Krt18','S100a11','Ncl','Rps3','Fth1','Pla2g7','B2m','Coro1a','Ttr','Gcg','Pyy','Sst')
  # testt<-testt[GENE,ce]
  # rownames(testt)<-GENE

  test2<-melt(as.matrix(testt),varnames = c("GENE","Clust"))
  
  ggplot(data = test2, aes(x=Clust, y=GENE, fill=value)) + 
    geom_tile()+scale_fill_gradient2(low = "blue", high = "red", mid = "white",midpoint = 0)+scale_x_discrete(limits = c(1:length(table(finalCluster[which(finalCluster!=0)]))))+ggtitle("ISClust")+theme(plot.title = element_text(hjust = 0.5)) 
  ggsave(paste(de_dir,"hm",".png",sep =""),plot = last_plot())
}


calcul.gini = function(x, unbiased = TRUE, na.rm = FALSE){
  if (!is.numeric(x)){
    warning("'x' is not numeric; returning NA")
    return(NA)
  }
  if (!na.rm && any(na.ind = is.na(x)))
    stop("'x' contain NAs")
  if (na.rm)
    x = x[!na.ind]
  n = length(x)
  mu = mean(x)
  N = if (unbiased) n * (n - 1) else n * n
  ox = x[order(x)]
  dsum = drop(crossprod(2 * 1:n - n - 1,  ox))
  dsum / (mu * N)
}

jaccard <- function(m) {
  ## common values:
  A = tcrossprod(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)
  
  ## only non-zero values of common
  Aim = A[im]
  
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )
  
  return( J )
}


Mean.in.log2space=function(x) {   #input is a vetor in log2 space, the return is their mean in log2 space.
  return(log2(mean(2^(x)-1)+1))
}

calcul.pca = function(x,n){
  use_genes <- which(colSums(x) > 0)
  m <- x[,use_genes]
  bc_tot <- rowSums(m)
  median_tot <- median(bc_tot)
  m <- sweep(m, 1, median_tot/bc_tot, '*')
  m[is.na(m)]<-0
  m <- log(1+m)
  m <- sweep(m, 2, colMeans(m), '-')
  ppk<-propack.svd(as.matrix(m),neig=n)
  pca<-t(ppk$d*t(ppk$u))
  list(ppk=ppk,pca=pca, m=m,use_genes=use_genes)
}


jaccard_dist_large_matrix <- function(m) { #use when matrix larger than 65,536 (2^16) cells
  A = tcrossprod(m)
  
  r<-35000 #65,536 was max possible before, now 70000 is- can increase r to increase this
  # approach is to split up matrix into 4, then recombine
  im = which(A[1:r,1:r] > 0,arr.ind=TRUE)
  im2 = which(A[(r+1):dim(A)[1],(r+1):dim(A)[1]] > 0,arr.ind=TRUE)
  im2[,1] = im2[,1]+r
  im2[,2] = im2[,2]+r
  im3 = which(A[1:r,(r+1):dim(A)[1]] > 0,arr.ind=TRUE)
  im3[,2]=im3[,2]+r
  im4 = which(A[(r+1):dim(A)[1],1:r] > 0,arr.ind=TRUE)
  im4[,1]=im4[,1]+r
  im_all<-rbind(im,im2,im3,im4)

  b = rowSums(m) 
  
  ## only non-zero values of common
  Aim = A[im_all]
  #print("line 4")
  
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im_all[,1],
    j = im_all[,2],
    x = Aim / (b[im_all[,1]] + b[im_all[,2]] - Aim),
    dims = dim(A)
  )
  
  J1<-J[1:r,1:r]
  oneMinusJ1<-1-J1
  index.cell.zero.1<-index.cell.zero[which(index.cell.zero<(r+1))]
  oneMinusJ1[index.cell.zero.1,index.cell.zero.1]=0
  
  J2<-J[(r+1):dim(A)[1],(r+1):dim(A)[1]]
  oneMinusJ2<-1-J2
  index.cell.zero.2<-index.cell.zero[which(index.cell.zero>r)]
  oneMinusJ2[index.cell.zero.2-r,index.cell.zero.2-r]=0
  
  J3<-J[1:r,(r+1):dim(A)[1]]
  oneMinusJ3<-1-J3
  oneMinusJ3[index.cell.zero.1,index.cell.zero.2-r]=0
  
  J4<-J[(r+1):dim(A)[1],1:r]
  oneMinusJ4<-1-J4
  oneMinusJ4[index.cell.zero.2-r,index.cell.zero.1]=0
  
  #combine:
  
  cell.cell.jaccard.distance<-matrix(rep(0,dim(A)[1]**2),ncol=dim(A)[1],nrow=dim(A)[1])
  cell.cell.jaccard.distance[1:r,1:r]<-as.vector(oneMinusJ1)
  cell.cell.jaccard.distance[(r+1):dim(A)[1],(r+1):dim(A)[1]]<-as.vector(oneMinusJ2)
  cell.cell.jaccard.distance[1:r,(r+1):dim(A)[1]]<-as.vector(oneMinusJ3)
  cell.cell.jaccard.distance[(r+1):dim(A)[1],1:r]<-as.vector(oneMinusJ4)
  
  return(cell.cell.jaccard.distance)
}

#for MAST DE
#additional functions
Mean.in.log2space=function(x,pseudo.count) {
  return(log2(mean(2^(x)-pseudo.count)+pseudo.count))
}

stat.log2=function(data.m, group.v, pseudo.count){
  #data.m=data.used.log2
  log2.mean.r <- aggregate(t(data.m), list(as.character(group.v)), function(x) Mean.in.log2space(x,pseudo.count))
  log2.mean.r <- t(log2.mean.r)
  colnames(log2.mean.r) <- paste("mean.group",log2.mean.r[1,], sep="")
  log2.mean.r = log2.mean.r[-1,]
  log2.mean.r = as.data.frame(log2.mean.r)
  log2.mean.r = unfactor(log2.mean.r)  #from varhandle
  log2.mean.r[,1] = as.numeric(log2.mean.r[,1])
  log2.mean.r[,2] = as.numeric(log2.mean.r[,2])
  log2_foldchange = log2.mean.r$mean.group1-log2.mean.r$mean.group0
  results = data.frame(cbind(log2.mean.r$mean.group0,log2.mean.r$mean.group1,log2_foldchange))
  colnames(results) = c("log2.mean.group0","log2.mean.group1","log2_fc")
  rownames(results) = rownames(log2.mean.r)
  return(results)
}


####### function m.auc  ######
#install.packages("ROCR")
library("ROCR")
v.auc = function(data.v,group.v) {
  prediction.use=prediction(data.v, group.v, 0:1)
  perf.use=performance(prediction.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  return(auc.use)
}
m.auc=function(data.m,group.v) {
  AUC=apply(data.m, 1, function(x) v.auc(x,group.v))
  AUC[is.na(AUC)]=0.5
  return(AUC)
  
}  
####### function m.auc END ######

ginfitting<-function(ExprM.RawCounts.filter){
  # .2 calculate (may also add more statistics in the Future 10.08.2015)
  if(gini.bi==0){
    gini = apply(as.data.frame(ExprM.RawCounts.filter), 1, function(x){calcul.gini(as.numeric(x)) } )    #theoretically, gini have very low chance to have a 1 value
    GiniIndex = as.data.frame(cbind(1:dim(ExprM.RawCounts.filter)[1], gini))
  } else {
    GiniIndex1 <- as.data.frame(apply(ExprM.RawCounts.filter, 1, function(x){calcul.gini(as.numeric(x)) } ) )
    GiniIndex2 <- as.data.frame(apply(ExprM.RawCounts.filter+0.00001, 1, function(x){calcul.gini(as.numeric(1/x)) } ) ) #bi directional
    GiniIndex  <- cbind(GiniIndex1, GiniIndex2)
    colnames(GiniIndex)=c("gini1","gini2")
    GiniIndex$gini2_sign = 0 - GiniIndex$gini2;
    GiniIndex$gini = apply(GiniIndex, 1, max)
    GiniIndex <- na.omit(GiniIndex)
    GiniIndex$gini_sign = GiniIndex$gini
    for(genei in 1:dim(GiniIndex)[1])
    {
      GiniIndex[genei, 5] = ifelse(  GiniIndex[genei, 1] > GiniIndex[genei,2], "up-regulation", "down-regulation") 
    }
    #dim(GiniIndex) 
    write.table(GiniIndex, file=paste("results/", exprimentID,"_bi-directional.GiniIndexTable.csv", sep=""), sep=",", row.names = TRUE,  col.names = TRUE,  quote = FALSE) 
  }
  
  Maxs          = apply(ExprM.RawCounts.filter,1,max)
  Means         = apply(ExprM.RawCounts.filter,1,mean)
  log2.Maxs     = log2(Maxs+0.1)
  ExprM.Stat1   = as.data.frame(cbind(Maxs,GiniIndex$gini,log2.Maxs))
  colnames(ExprM.Stat1) = c("Maxs","Gini","log2.Maxs")
  ExprM.Stat1 = ExprM.Stat1[ExprM.Stat1$log2.Maxs>log2.expr.cutoffl & ExprM.Stat1$log2.Maxs<=log2.expr.cutoffh ,]  # is this necessary?
  log2.Maxs = ExprM.Stat1$log2.Maxs
  Gini      = ExprM.Stat1$Gini
  Maxs      = ExprM.Stat1$Maxs
  
  
  # .3 fitting in max-gini space 
  Gini.loess.fit        = loess(Gini~log2.Maxs, span=span, degree=1)
  Normlized.Gini.Score  = Gini.loess.fit$residuals   #residuals = Gini - Gini.fitted
  Gini.fitted           = Gini.loess.fit$fitted    
  ExprM.Stat1           = as.data.frame(cbind(ExprM.Stat1[,c("Maxs","Gini", "log2.Maxs")], Normlized.Gini.Score, Gini.fitted))
  colnames(ExprM.Stat1) = c("Maxs","Gini","log2.Maxs", "Norm.Gini", "Gini.fitted")
  ### remove 25% of first round outlier genes, do second round loess
  Gini.loess.fit.residual = residuals(Gini.loess.fit)                               
  thresh.outlier = quantile(Gini.loess.fit.residual[Gini.loess.fit.residual>0], outlier_remove) 
  id.genes.loess.fit = which(Gini.loess.fit.residual < thresh.outlier)               
  id.outliers.loess.fit = which(Gini.loess.fit.residual >= thresh.outlier)          
  log2.Maxs.genes = log2.Maxs[id.genes.loess.fit]                                   
  log2.Maxs.outliers = log2.Maxs[id.outliers.loess.fit]                            
  Gini.loess.fit.2 = loess(Gini[id.genes.loess.fit]~log2.Maxs[id.genes.loess.fit], span=span, degree = 1)
  Gini.loess.fit.2.predict = predict(Gini.loess.fit.2)  
  
  #plot second round fit
  Gini.loess.fit.2.x.y = cbind(log2.Maxs.genes,Gini.loess.fit.2.predict)
  Gini.loess.fit.2.x.y.uniq = as.data.frame(unique(Gini.loess.fit.2.x.y))
  Gini.loess.fit.2.x.y.uniq = Gini.loess.fit.2.x.y.uniq[order(Gini.loess.fit.2.x.y.uniq[,1]),]
  log2.Maxs.genes.sorted = log2.Maxs.genes[order(log2.Maxs.genes)]                   
  Gini.loess.fit.2.predict.sorted = Gini.loess.fit.2.predict[order(log2.Maxs.genes)] 
  #using Gini.loess.fit.2 as model, predict gini value for those outlier which are not used for build model.
  #for each max in outliers set, find the id of max value which is most close in fitted data set
  loc.outliers = apply(matrix(log2.Maxs.outliers),1,function(x){
    if(x<max(log2.Maxs.genes.sorted)){
      return(which(log2.Maxs.genes.sorted>=x)[1])
    }else{
      return(which.max(log2.Maxs.genes.sorted))
    }})                
  #check the results
  outlier_max_in_fit <- cbind(log2.Maxs.outliers, loc.outliers, log2.Maxs.genes.sorted[loc.outliers])
  
  #based on Gini.loess.fit.2, predict outliers which was not used for fitting
  Gini.outliers.predict = apply(cbind(seq(length(log2.Maxs.outliers)),log2.Maxs.outliers),1,function(x){
    id = x[1]
    value = x[2]
    if(value == log2.Maxs.genes.sorted[loc.outliers[id]]){
      return(as.numeric(Gini.loess.fit.2.x.y.uniq[which(Gini.loess.fit.2.x.y.uniq$log2.Maxs.genes>=value)[1],2]))
    }else{
      if(loc.outliers[id]>1){
        return(Gini.loess.fit.2.predict.sorted[loc.outliers[id]-1]+(Gini.loess.fit.2.predict.sorted[loc.outliers[id]]-Gini.loess.fit.2.predict.sorted[loc.outliers[id]-1])*(value-log2.Maxs.genes.sorted[loc.outliers[id]-1])/(log2.Maxs.genes.sorted[loc.outliers[id]]-log2.Maxs.genes.sorted[loc.outliers[id]-1]))
      }else{
        return(Gini.loess.fit.2.predict.sorted[2]-(Gini.loess.fit.2.predict.sorted[2]-Gini.loess.fit.2.predict.sorted[1])*(log2.Maxs.genes.sorted[2]-value)/(log2.Maxs.genes.sorted[2]-log2.Maxs.genes.sorted[1]))
      }
    }
  })
  
  #plot outliers predict results
  outliers.precit.x.y.uniq = as.data.frame(unique(cbind(log2.Maxs.outliers, Gini.outliers.predict)))
  #plot(outliers.precit.x.y.uniq)
  #plot whole fit2 
  colnames(outliers.precit.x.y.uniq) = colnames(Gini.loess.fit.2.x.y.uniq)
  Gini.loess.fit.2.full.x.y.uniq = rbind(Gini.loess.fit.2.x.y.uniq, outliers.precit.x.y.uniq)
  #plot(Gini.loess.fit.2.full.x.y.uniq)
  
  #calcualte Normlized.Gini.Score2
  Normlized.Gini.Score2                        = rep(0,length(Gini.loess.fit.residual))               
  Normlized.Gini.Score2[id.genes.loess.fit]    = residuals(Gini.loess.fit.2)                         
  Normlized.Gini.Score2[id.outliers.loess.fit] = Gini[id.outliers.loess.fit] - Gini.outliers.predict 
  
  Gini.fitted2           = Gini - Normlized.Gini.Score2         
  ExprM.Stat1            = as.data.frame(cbind(ExprM.Stat1[,c("Maxs","Gini", "log2.Maxs", "Gini.fitted", "Norm.Gini" )], Gini.fitted2, Normlized.Gini.Score2))
  colnames(ExprM.Stat1)  = c("Maxs","Gini","log2.Maxs", "Gini.fitted","Norm.Gini",  "Gini.fitted2", "Norm.Gini2")
  Gini.pvalue            = pnorm(-abs(scale(ExprM.Stat1$Norm.Gini2, center=TRUE,scale=TRUE)))
  ExprM.Stat2            = cbind(ExprM.Stat1, Gini.pvalue)  #first time use ExprM.Stat2
  
  #for each measurement, first ranked by themself.
  # .4 identify High Gini Genes with Norm.Gini
  #ExprM.Stat2         = ExprM.Stat2[order(ExprM.Stat2$Norm.Gini2,decreasing=T),]
  #Genelist.HighNormGini = rownames(ExprM.Stat2[intersect(which(ExprM.Stat2$Norm.Gini2 > sort(ExprM.Stat2$Norm.Gini2,decreasing = T)[50]),which(!is.nan(ExprM.Stat2$Norm.Gini2))),])  # cutoff approach, use a cutoff to choose top genes. 
  #length(Genelist.HighNormGini)
  
  ExprM.Stat2         = ExprM.Stat2[order(ExprM.Stat2$Gini.pvalue),]
  Genelist.top_pvalue = rownames(ExprM.Stat2[intersect(which(ExprM.Stat2$Gini.pvalue <= max(Gini.pvalue_cutoff,sort(ExprM.Stat2$Gini.pvalue[which(ExprM.Stat2$Norm.Gini2 > 0)],decreasing = F)[30])),which(ExprM.Stat2$Norm.Gini2 > 0)),])
  #Genelist.top_pvalue = rownames(ExprM.Stat2[intersect(which(ExprM.Stat2$Gini.pvalue <= max(Gini.pvalue_cutoff)),which(ExprM.Stat2$Norm.Gini2 > 0)),])
  length(Genelist.top_pvalue)
  return(Genelist.top_pvalue)
}
