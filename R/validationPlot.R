#' Plots to decide number of clusters
#'
#' The function generates multiple plots that help users to decide number of clusters or MetICA components.
#'
#' @param M1 The entire list object generated from the function MetICA(X, pcs = 15...)
#' @param cluster_index Boolean object. TRUE if geometric indices for clusterinig quality needs to be calculted. The calculation provides an additional criterion for cluster number selection but can be time-consuming.
#'
#' @return Return an overall graph with multiple plots. If cluster_index = TRUE, also return a list object containing clustering quality score.
#'
#' @export

validationPlot<-function(M1,cluster_index=T){

  memory.size(max=2000000000)
  options(warn=-1)

  # This function helps to find the best number of clusters by visualizing

  new_component_size=sapply(M1$Stage3$Kurt_history,length) # How many new components each time
  cluster_list=list()
  for (i in 1:M1$Stage3$max.cluster){cluster_list[[i]]=rep(i,new_component_size[i])}

  # Evaluate newly added clusters when increasing cluster number

  par(mfrow=c(1,2))

  plot(unlist(cluster_list),unlist(M1$Stage3$tn_history),col="blue",pch=19,xlab="Number of clusters",ylab="Number of estimates",font.lab=2,font.axis=2,bty="n")
  box(lwd=2)
  title("Total number of estimates in each cluster")

  plot(unlist(cluster_list),unlist(M1$Stage3$Kurt_history),col="blue",pch=19,xlab="Number of clusters",ylab="Kurtosis",font.lab=2,font.axis=2,bty="n")
  box(lwd=2)
  title("Kurtosis of cluster centers (components)")

  plot(unlist(cluster_list),unlist(M1$Stage3$bn_history),col="blue",pch=19,xlab="Number of clusters",ylab="Proportion",font.lab=2,font.axis=2,bty="n")
  box(lwd=2)
  polygon(c(-1,100,100,-1),c(0.45,0.45,1,1),col=rgb(0.8, 0.9, 0.8,0.5), border = "grey",lty="dashed")
  title("Proportion of bootstrap estimates in each cluster")

  plot(unlist(cluster_list),unlist(M1$Stage3$boot_eval),col="blue",pch=19,xlab="Number of clusters",ylab="Bootstrap distance",font.lab=2,font.axis=2,bty="n")
  box(lwd=2)
  polygon(c(-1,100,100,-1),c(0,0,0.8,0.8),col=rgb(0.8, 0.9, 0.8,0.5), border = "grey",lty="dashed")
  title("Cluster stability against bootstrapping")

  plot(unlist(cluster_list),unlist(M1$Stage3$dispersion_history),col="blue",pch=19,xlab="Number of clusters",ylab="Dispersion index",font.lab=2,font.axis=2,bty="n")
  box(lwd=2)
  polygon(c(-1,100,100,-1),c(0,0,0.8,0.8),col=rgb(0.8, 0.9, 0.8,0.5), border = "grey",lty="dashed")
  title("Dispersion index")

  if (cluster_index){
    message("Generating geometric index... please be patient!")
    results=R_index(M1)
    plot(results$x,results$y,col="blue",pch=19,xlab="Number of clusters",ylab="R index",font.lab=2,font.axis=2,bty="n")
    box(lwd=2)
    title("Geometric index of clusters")
    suggested=which.min(results$y)
    message("Suggested cluster number according to the geometric index is:",suggested)
    return(results)
  }
  else{return("Please set cluster_index=T if you cannont make a decision!")}

  par(mfrow=c(1,1))
}

R_index<-function(M1){

  S_sum=M1$Stage1$S
  trends=M1$Stage2$trends
  max.cluster=M1$Stage3$max.cluster

  # Recalculate distance object:

  if (ncol(S_sum)>4000){ # Faster correlation matrix for big data
    R=bigcor(S_sum,size= 2000,method="spearman",verbose = F) # Big correlation matrix
    R=R[1:nrow(R), 1:ncol(R)]}
  else{ # Smaller data size
    R=cor(S_sum,method="spearman")}

  if (trends){D=(1-R)/2} # Disimilariy Matrix for time dependent data [0,1]
  else {D=1-abs(R)} # Disimilariy Matrix for no-time dependent data [0,1]

  # Calculate score
  score_list=c()
  for (nbcluster in 2:max.cluster){
    cluster=cutree(M1$Stage2$clusterObj,nbcluster) # CUtting tree
    cluster_list=1:nbcluster
    cluster_scores=c()
    for (i in cluster_list){
      index_in=which(cluster==i) # Index of elenments in cluster i
      Cm_in=length(index_in) # Nb of elements in the cluster
      Var_in=1/(Cm_in^2)*sum(D[index_in,index_in]) # Variance inside the cluster
      Var_out=c()
      for (j in cluster_list[-i]){
        index_out=which(cluster==j)
        Cm_out=length(index_out) # Nb of elements in cluster j
        Var_out=c(Var_out,1/(Cm_in*Cm_out)*sum(D[index_in,index_out]))} # Variance between cluster i and j
      Var_out=min(Var_out) # Take the minimal variance
      cluster_scores=c(cluster_scores,Var_in/Var_out) # Score of the evaluated cluster
    }
    score_list=c(score_list,mean(cluster_scores))}
  return(list(x=2:max.cluster,y=score_list))
}

