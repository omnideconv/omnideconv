#solve using OLS, constrained such that cell type numbers>0
#NEEDED
solveOLS<-function(S,B){
  D<-t(S)%*%S
  d<-t(S)%*%B
  A<-base::cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  sc <- norm(D,"2")
  solution<-quadprog::solve.QP(D/sc,d/sc,A,bzero)$solution
  names(solution)<-colnames(S)
  #print(round(solution/sum(solution),5))
  return(solution/sum(solution))
}

#return cell number, not proportion
#do not print output
#NEEDED
solveOLSInternal<-function(S,B){
  D<-t(S)%*%S
  d<-t(S)%*%B
  A<-base::cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  sc <- norm(D,"2")
  solution<-quadprog::solve.QP(D/sc,d/sc,A,bzero)$solution
  names(solution)<-colnames(S)
  return(solution)
}

#solve using WLS with weights dampened by a certain dampening constant
#NEEDED
solveDampenedWLS<-function(S,B){
  #first solve OLS, use this solution to find a starting point for the weights
  solution<-solveOLSInternal(S,B)
  #now use dampened WLS, iterate weights until convergence
  iterations<-0
  changes<-c()
  #find dampening constant for weights using cross-validation
  j<-findDampeningConstant(S,B,solution)
  change<-1
  while(change>.01 & iterations<1000){
    newsolution<-solveDampenedWLSj(S,B,solution,j)
    #decrease step size for convergence
    solutionAverage<-base::rowMeans(base::cbind(newsolution,matrix(solution,nrow = length(solution),ncol = 4)))
    change<-norm(as.matrix(solutionAverage-solution))
    solution<-solutionAverage
    iterations<-iterations+1
    changes<-c(changes,change)
  }
  #print(round(solution/sum(solution),5))
  return(solution/sum(solution))
}

#solve WLS given a dampening constant
#NEEDED
solveDampenedWLSj<-function(S,B,goldStandard,j){
  multiplier<-1*2^(j-1)
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsDampened<-wsScaled
  wsDampened[base::which(wsScaled>multiplier)]<-multiplier
  W<-base::diag(wsDampened)
  D<-t(S)%*%W%*%S
  d<- t(S)%*%W%*%B
  A<-base::cbind(diag(dim(S)[2]))
  bzero<-c(rep(0,dim(S)[2]))
  sc <- norm(D,"2")
  solution<-quadprog::solve.QP(D/sc,d/sc,A,bzero)$solution
  names(solution)<-colnames(S)
  return(solution)
}

#find a dampening constant for the weights using cross-validation
#NEEDED
findDampeningConstant<-function(S,B,goldStandard){
  solutionsSd<-NULL
  #goldStandard is used to define the weights
  sol<-goldStandard
  ws<-as.vector((1/(S%*%sol))^2)
  wsScaled<-ws/min(ws)
  wsScaledMinusInf<-wsScaled
  #ignore infinite weights
  if(max(wsScaled)=="Inf"){
    wsScaledMinusInf<-wsScaled[-base::which(wsScaled=="Inf")]
  }
  #try multiple values of the dampening constant (multiplier)
  #for each, calculate the variance of the dampened weighted solution for a subset of genes
  for (j in 1:base::ceiling(log2(max(wsScaledMinusInf)))){
    multiplier<-1*2^(j-1)
    wsDampened<-wsScaled
    wsDampened[which(wsScaled>multiplier)]<-multiplier
    solutions<-NULL
    seeds<-c(1:100)
    for (i in 1:100){
      base::set.seed(seeds[i]) #make nondeterministic
      subset<-sample(length(ws),size=length(ws)*0.5) #randomly select half of gene set
      #solve dampened weighted least squares for subset
      fit = lm (B[subset] ~ -1+S[subset,],weights=wsDampened[subset])
      sol<-fit$coef*sum(goldStandard)/sum(fit$coef)
      solutions<-base::cbind(solutions,sol)
    }
    solutionsSd<-base::cbind(solutionsSd,base::apply(solutions,1,sd))
  }
  #choose dampening constant that results in least cross-validation variance
  j<-which.min(base::colMeans(solutionsSd^2))
  return(j)
}

#NEEDED
solveSVR<-function(S,B){
  #scaling
  ub=max(c(as.vector(S),B)) #upper bound
  lb=min(c(as.vector(S),B)) #lower bound
  Bs=((B-lb)/ub)*2-1
  Ss=((S-lb)/ub)*2-1

  #perform SVR
  model<-e1071::svm(Ss,Bs, nu=0.5,scale = TRUE, type = "nu-regression",kernel ="linear",cost = 1)
  coef <- t(model$coefs) %*% model$SV
  coef[base::which(coef<0)]<-0
  coef<-as.vector(coef)
  names(coef)<-colnames(S)
  #print(round(coef/sum(coef),5))
  return(coef/sum(coef))
}

#perform DE analysis using Seurat
#when path = NULL, the generated files in the processes will not be saved and output.
DEAnalysis<-function(scdata,id,path){
  list.names <- unique(id)
  list_de_group<-as.list(rep(0, length(list.names)))

  exprObj<-Seurat::CreateSeuratObject(raw.data=as.data.frame(scdata), project = "DE")
  exprObj2<-Seurat::SetIdent(exprObj,ident.use=as.vector(id))
  #print("Calculating differentially expressed genes:")
  for (i in unique(id)){
    de_group <- Seurat::FindMarkers(object=exprObj2, ident.1 = i, ident.2 = NULL,
                            only.pos = TRUE, test.use = "bimod")

    index<-which(list.names==i)
    list_de_group[[index]]<-de_group

    if(!is.null(path)){
      save(de_group,file=paste(path,"/de_",i,".RData",sep=""))
    }
  }
  return(list_de_group)
}



#build signature matrix using genes identified by DEAnalysis()
#when path = NULL, the generated files in the processes will not be saved and output.
buildSignatureMatrixUsingSeurat<-function(scdata,id,path,diff.cutoff=0.5,pval.cutoff=0.01){

  #perform differential expression analysis
  list_de_groups <- DEAnalysis(scdata,id,path)

  numberofGenes<-c()
  for (i in unique(id)){
    name_list<-unique(id)
    index<-which(name_list==i)
    de_group <- list_de_groups[[index]]

    DEGenes<-rownames(de_group)[intersect(which(de_group$p_val_adj<pval.cutoff),which(de_group$avg_logFC>diff.cutoff))]
    nonMir = base::grep("MIR|Mir", DEGenes, invert = T)
    base::assign(paste("cluster_lrTest.table.",i,sep=""),de_group[which(rownames(de_group)%in%DEGenes[nonMir]),])
    numberofGenes<-c(numberofGenes,length(DEGenes[nonMir]))
  }

  #need to reduce number of genes
  #for each subset, order significant genes by decreasing fold change, choose between 50 and 200 genes
  #choose matrix with lowest condition number
  conditionNumbers<-c()
  for(G in 50:200){
    Genes<-c()
    j=1
    for (i in unique(id)){
      if(numberofGenes[j]>0){
        temp<-base::paste("cluster_lrTest.table.",i,sep="")
        temp<-get(temp)
        temp<-temp[order(temp$p_val_adj,decreasing=TRUE),]
        Genes<-c(Genes,(rownames(temp)[1:min(G,numberofGenes[j])]))
      }
      j=j+1
    }
    Genes<-unique(Genes)
    #make signature matrix
    ExprSubset<-scdata[Genes,]
    Sig<-NULL
    for (i in unique(id)){
      Sig<-base::cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
    }
    colnames(Sig)<-unique(id)
    conditionNumbers<-c(conditionNumbers,kappa(Sig))
  }
  G<-which.min(conditionNumbers)+min(49,numberofGenes-1) #G is optimal gene number
  #
  Genes<-c()
  j=1
  for (i in unique(id)){
    if(numberofGenes[j]>0){
      temp<-base::paste("cluster_lrTest.table.",i,sep="")
      temp<-get(temp)
      temp<-temp[order(temp$p_val_adj,decreasing=TRUE),]
      Genes<-c(Genes,(rownames(temp)[1:min(G,numberofGenes[j])]))
    }
    j=j+1
  }
  Genes<-unique(Genes)
  ExprSubset<-scdata[Genes,]
  Sig<-NULL
  for (i in unique(id)){
    Sig<-base::cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig)<-unique(id)

  if(!is.null(path)){
    save(Sig,file=paste(path,"/Sig.RData",sep=""))
  }

  return(Sig)
}

##alternative differential expression method using MAST

#functions for DE

Mean.in.log2space=function(x,pseudo.count) {
  return(log2(mean(2^(x)-pseudo.count)+pseudo.count))
}

stat.log2=function(data.m, group.v, pseudo.count){
  #data.m=data.used.log2
  log2.mean.r <- stats::aggregate(t(data.m), list(as.character(group.v)), function(x) Mean.in.log2space(x,pseudo.count))
  log2.mean.r <- t(log2.mean.r)
  colnames(log2.mean.r) <- paste("mean.group",log2.mean.r[1,], sep="")
  log2.mean.r = log2.mean.r[-1,]
  log2.mean.r = as.data.frame(log2.mean.r)
  log2.mean.r = varhandle::unfactor(log2.mean.r)  #from varhandle
  log2.mean.r[,1] = as.numeric(log2.mean.r[,1])
  log2.mean.r[,2] = as.numeric(log2.mean.r[,2])
  log2_foldchange = log2.mean.r$mean.group1-log2.mean.r$mean.group0
  results = data.frame(base::cbind(log2.mean.r$mean.group0,log2.mean.r$mean.group1,log2_foldchange))
  colnames(results) = c("log2.mean.group0","log2.mean.group1","log2_fc")
  rownames(results) = rownames(log2.mean.r)
  return(results)
}

v.auc = function(data.v,group.v) {
  prediction.use=ROCR::prediction(data.v, group.v, 0:1)
  perf.use=ROCR::performance(prediction.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  return(auc.use)
}
m.auc=function(data.m,group.v) {
  AUC=apply(data.m, 1, function(x) v.auc(x,group.v))
  AUC[is.na(AUC)]=0.5
  return(AUC)

}

#perform DE analysis using MAST
#when path = NULL, the generated files in the processes will not be saved and output.
#NEEDED
DEAnalysisMAST<-function(scdata,id,path){
  list.names <- unique(id)
  list_lrTest.table<-as.list(rep(0, length(list.names)))

  pseudo.count = 0.1
  data.used.log2   <- log2(scdata+pseudo.count)
  colnames(data.used.log2)<-make.unique(colnames(data.used.log2))
  diff.cutoff=0.5
  for (i in unique(id)){
    cells.symbol.list2     = colnames(data.used.log2)[which(id==i)]
    cells.coord.list2      = match(cells.symbol.list2, colnames(data.used.log2))
    cells.symbol.list1     = colnames(data.used.log2)[which(id != i)]
    cells.coord.list1      = match(cells.symbol.list1, colnames(data.used.log2))
    data.used.log2.ordered  = cbind(data.used.log2[,cells.coord.list1], data.used.log2[,cells.coord.list2])
    group.v <- c(rep(0,length(cells.coord.list1)), rep(1, length(cells.coord.list2)))
    #ouput
    log2.stat.result <- stat.log2(data.used.log2.ordered, group.v, pseudo.count)
    Auc <- m.auc(data.used.log2.ordered, group.v)
    bigtable <- data.frame(cbind(log2.stat.result, Auc))

    DE <- bigtable[bigtable$log2_fc >diff.cutoff,]
    dim(DE)
    if(dim(DE)[1]>1){
      data.1                 = data.used.log2[,cells.coord.list1, drop=F]
      data.2                 = data.used.log2[,cells.coord.list2, drop=F]
      genes.list = rownames(DE)
      log2fold_change        = cbind(genes.list, DE$log2_fc)
      colnames(log2fold_change) = c("gene.name", "log2fold_change")
      counts  = as.data.frame(cbind( data.1[genes.list,], data.2[genes.list,] ))
      groups  = c(rep("Cluster_Other", length(cells.coord.list1) ), rep(i, length(cells.coord.list2) ) )
      groups  = as.character(groups)
      data_for_MIST <- as.data.frame(cbind(rep(rownames(counts), dim(counts)[2]), reshape::melt(counts),rep(groups, each = dim(counts)[1]), rep(1, dim(counts)[1] * dim(counts)[2]) ))
      colnames(data_for_MIST) = c("Gene", "Subject.ID", "Et", "Population", "Number.of.Cells")
      vbeta = data_for_MIST
      vbeta.fa <-MAST::FromFlatDF(vbeta, idvars=c("Subject.ID"),
                            primerid='Gene', measurement='Et', ncells='Number.of.Cells',
                            geneid="Gene",  cellvars=c('Number.of.Cells', 'Population'),
                            phenovars=c('Population'), id='vbeta all')
      vbeta.1 <- subset(vbeta.fa,Number.of.Cells==1)
      # .3 MAST
      head(SummarizedExperiment::colData(vbeta.1))
      zlm.output <- MAST::zlm(~ Population, vbeta.1, method='bayesglm', ebayes=TRUE)
      show(zlm.output)
      coefAndCI <- summary(zlm.output, logFC=TRUE)
      zlm.lr <- MAST::lrTest(zlm.output, 'Population')
      zlm.lr_pvalue <- reshape::melt(zlm.lr[,,'Pr(>Chisq)'])
      zlm.lr_pvalue <- zlm.lr_pvalue[which(zlm.lr_pvalue$test.type == 'hurdle'),]



      lrTest.table <-  merge(zlm.lr_pvalue, DE, by.x = "primerid", by.y = "row.names")
      colnames(lrTest.table) <- c("Gene", "test.type", "p_value", paste("log2.mean.", "Cluster_Other", sep=""), paste("log2.mean.",i,sep=""), "log2fold_change", "Auc")
      cluster_lrTest.table <- lrTest.table[rev(order(lrTest.table$Auc)),]

      #. 4 save results

      index<-which(list.names==i)
      list_lrTest.table[[index]]<-cluster_lrTest.table

      if(!is.null(path)){
        write.csv(cluster_lrTest.table, file=paste(path,"/",i,"_lrTest.csv", sep=""))
        save(cluster_lrTest.table, file=paste(path,"/",i,"_MIST.RData", sep=""))}
    }
  }

  return(list_lrTest.table)
}

#build signature matrix using genes identified by DEAnalysisMAST()
#when path = NULL, the generated files in the processes will not be saved and output.

buildSignatureMatrixMAST<-function(scdata,id, path, diff.cutoff=0.5,pval.cutoff=0.01){

  id <- gsub(" ","_",id)

  #compute differentially expressed genes for each cell type
  list_cluster_table <- DEAnalysisMAST(scdata,id,path)

  #for each cell type, choose genes in which FDR adjusted p-value is less than 0.01 and the estimated fold-change
  #is greater than 0.5
  numberofGenes<-c()
  for (i in unique(id)){
    name_list<-unique(id)
    index<-which(name_list==i)
    cluster_lrTest.table = list_cluster_table[[index]]
    pvalue_adjusted<-p.adjust(cluster_lrTest.table[,3], method = "fdr", n = length(cluster_lrTest.table[,3]))
    cluster_lrTest.table<-cbind(cluster_lrTest.table,pvalue_adjusted)
    DEGenes<-cluster_lrTest.table$Gene[intersect(which(pvalue_adjusted<pval.cutoff),which(cluster_lrTest.table$log2fold_change>diff.cutoff))]
    nonMir = grep("MIR|Mir", DEGenes, invert = T)  # because Mir gene is usually not accurate
    assign(paste("cluster_lrTest.table.",i,sep=""),cluster_lrTest.table[which(cluster_lrTest.table$Gene%in%DEGenes[nonMir]),])
    numberofGenes<-c(numberofGenes,length(DEGenes[nonMir]))
  }


  #need to reduce number of genes
  #for each subset, order significant genes by decreasing fold change, choose between 50 and 200 genes
  #for each, iterate and choose matrix with lowest condition number
  conditionNumbers<-c()
  for(G in 50:200){
    Genes<-c()
    j = 1
    for (i in unique(id)){
      if(numberofGenes[j]>0){
        temp<-paste("cluster_lrTest.table.",i,sep="")
        temp<-get(temp)
        temp<-temp[order(temp$log2fold_change,decreasing=TRUE),]
        Genes<-c(Genes,varhandle::unfactor(temp$Gene[1:min(G,numberofGenes[j])]))
      }
      j=j+1
    }
    Genes<-unique(Genes)
    #make signature matrix
    ExprSubset<-scdata[Genes,]
    Sig<-NULL
    for (i in unique(id)){
      Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
    }
    colnames(Sig)<-unique(id)
    conditionNumbers<-c(conditionNumbers,kappa(Sig))
  }
  G<-which.min(conditionNumbers)+min(49,numberofGenes-1)
  Genes<-c()
  j=1
  for (i in unique(id)){
    if(numberofGenes[j]>0){
      temp<-paste("cluster_lrTest.table.",i,sep="")
      temp<-get(temp)
      temp<-temp[order(temp$log2fold_change,decreasing=TRUE),]
      Genes<-c(Genes,varhandle::unfactor(temp$Gene[1:min(G,numberofGenes[j])]))
    }
    j=j+1
  }
  Genes<-unique(Genes)
  ExprSubset<-scdata[Genes,]
  Sig<-NULL
  for (i in unique(id)){
    Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }

  colnames(Sig)<-unique(id)
  rownames(Sig)<-Genes

  if(!is.null(path)){
    save(Sig,file=paste(path,"/Sig.RData",sep=""))}


  return(Sig)
}
