#' MetICA simulations on metabolomics data
#'
#' The main function for MetICA simulation on a sample × variables (n × p) metabolomics data matrix.
#' @param X A numeric matrix obtained from metabolomics expriments. Its dimension should be n (samples) × p (metabolic features), either centered or not. Notice: scaled data is not recommended for MetICA!
#' @param pcs Number of principal components used to whiten the data before ICA, also number of components estimated in each IPCA run. It should be at least 3. Its value can be modified after the function is launched and percentage of variance explained is calculated. We recommend that PCA whitening should keep at least 80 percent of total variance.
#' @param max_iter Number of IPCA iterations. It should be at least 50 to provide reliable results. More than 500 runs can lead to long computational time. To avoid computer memory issues, the total number of estimates (pcs × max_iter) must be under 20 000.
#' @param boot.prop Proportion of samples replaced in bootstrap iterations (when X is resampled). It should not exceed 0.4.
#' @param max.cluster The number of clusters in HCA of estimated components is evaluated from 2 to max.cluster. Its value can be modified in the function if one cluster contains fewer than 30 estimates.
#' @param trends Boolean variable. TRUE if your observations are time-dependent (e.g. blood samples taken over a period of time from a patient).
#' @param verbose Boolean variable. If TRUE the completion of each stage of the algorithm will be reminded by a message.
#'
#' @details MetICA is a three-stage algorithm:
#' \itemize{
#' \item{Stage_1 We performs multiple IPCA iterations on data matrix X. IPCA from mixOmics package combines the advantages of both PCA and ICA, and is therefore dedicated to large biological datasets. For all IPCA iterations, the initial un-mixing matrices are randomly generated from gaussian distributions. For half of our simulations, the initial data matrix X is re-sampled through bootstrapping. We use parallel-based FastICA for the first half of IPCA runs and deflation-based for the rest.}
#' \item{Stage_2 Estimated components from all IPCA runs are then submitted to hierarchical cluster analysis (HCA). The spearman distance metric is used since the objective of MetICA is to group estimated components that have similar shapes or trends.}
#' \item{Stage_3 One general difficulty of ICA on metabolomics data analysis is the selection of number of component. In MetICA, the problem becomes the choice of cluster number since components generated are the geometric centers of clusters. By increasing the number of clusters in HCA, geometric indexes of each cluster are calculated to reflect the convergence of IPCA algorithm. The stability of clusters against bootstrapping is evaluated as well as the kurtosis of cluster centers.}
#'}
#'
#' @references A. Hyvarinen and E. Oja, Independent Component Analysis: Algorithms and Applications, Neural Networks (2000) vol. 13 no. 4-5
#' @references Fangzhou Yao, Jeff Coquery and Kim-Anh Le Cao, Independent Principal Component Analysis for biologically meaningful dimension reduction of large biological data sets, BMC Bioinformatics (2012) Vol. 13 no. 24
#' @references Youzhong Liu, Kirill Smirnov, Marianna Lucio, Regis D. Gougeon, Herve Alexandre and Philippe Schmitt-Kopplin, MetICA: independent component analysis for high-resolution mass-spectrometry based non-targeted metabolomics, BMC Bioinformatics (2016) Vol. 17 no. 114
#'
#' @return A model (a list object) that contains results from each stage of the simulation
#' \itemize{
#'  \item{Stage1$S Estimated scores from IPCA runs. Data matrix contains n rows (samples) and pcs × max_iter columns (estimated components)}
#'  \item{Stage1$A Estimated loadings from IPCA runs. Data matrix contains p rows (metabolic features) and pcs × max_iter columns (estimated components)}
#'  \item{Stage1$boot_id pcs × max_iter vector indicating bootstrap/without bootstrap iterations.}
#'  \item{Stage1$type_id pcs × max_iter vector indicating deflation/parallel iterations.}
#'  \item{Stage2$clusterObj A object of class hclust from HCA analysis of estimated components.}
#'  \item{Stage2$trends Boolean object defined by users.}
#'  \item{Stage3$max.cluster Maximal number of clusters evaluated.}
#'  \item{Stage3$S_history List object of length max.cluster. Each element of the list (from 2nd element on) is a data matrix with n rows (samples) and nb columns (number of clusters) that represents center of each generated cluster.}
#'  \item{Stage3$A_history List object of length max.cluster. Each element of the list is a data matrix with p columns (metabolic features) and nb columns (number of clusters) that represents center of each generated cluster.}
#'  \item{Stage3$kurt_history List object of length max.cluster. Each element of the list is a numeric vector that describes the kurtosis of each cluster center obtained.}
#'  \item{Stage3$tn_history Number of estimates in each cluster. }
#'  \item{Stage3$bn_history Proportion of bootstrap estimates in each cluster.}
#'  \item{Stage3$dispersion_history Dispersion index of each cluster. Lower index (closer to 0) indicates a compact cluster thus a good algorithm convergence.}
#'  \item{Stage3$boot_eval Bootstrap stability index of each cluster. Lower index (closer to 0) means that there's no big difference between bootstrap and no bootstrap estimates, and the cluster is stable towards bootstrapping.}
#' }
#'
#' @author Youzhong Liu, \email{Youzhong.Liu@uantwerpen.be}
#'
#' @examples
#' data(bacteria_peptides)
#' # Perform 100 IPCA simulations on centered metabolomics data:
#' M1=MetICA(bacteria_peptides$X,pcs = 20,max_iter = 100,boot.prop = 0.3,max.cluster = 40,trends = T)
#' # Generate validation plots along with geometric index calculation to help decide number of clusters
#' validationPlot(M1)
#' # According to the validation, we choose 10 components
#' M2=MetICA_extract_model(M1,10,tops=7)
#'
#' @export


MetICA<-function(X, pcs = 15, max_iter = 400, boot.prop = 0.3, max.cluster = 20, trends = T, verbose=T){

  t=as.numeric(Sys.time())
  set.seed((t - floor(t))*1e8->seed)
  memory.size(max=2000000000)
  options(warn=-1)

  ##############################
  ### Check function inputs:
  ##############################

  if (missing(X) || (!is.matrix(X))) {
    stop("Please enter a profile matrix (samples x features)!")}

  if (pcs<3){
    message("Too few pcs for whitening: reset pcs to 3")
    pcs=3}

  checked=F
  while (!checked){
    prep_model<-ipca(X,pcs) # Precheck the model
    tot_var = sum(prep_model$explained_variance)*100
    message(pcs," principal components explained about ",round(tot_var,2),"% of variance")
    switch(menu(c("Yes", "No"), title="Do you want to whiten the data with this number of principal components ?"),
    Yes={
      checked=T},
    No={
      pcs=readline("Please enter a new number: ")
      pcs=as.integer(pcs)}
    )}

  if (max_iter < 50){
    message("'Too few iterations: reset max_iter to 50")
    max_iter=50}

  if (pcs*max_iter>20000){
    message("Out-of-memory problem: reset max_iter to ", floor(20000/pcs))
    max_iter=floor(20000/pcs)}

  if (pcs*max_iter>10000){
    message("The simulation can be very long, please have a coffee and take a rest!")
  }

  if (boot.prop>0.4){
    message("No more than 40% samples should be replaced: reset boot.prop to 0.4")
    boot.prop = 0.4}

  nbt = round(nrow(X)*boot.prop) # Number of samples replaced in bootstrap simuations
  if (nbt<3){
    message("At least 3 samples should be replaced: reset boot.prop to ",round(3/nrow(X),2))
    nbt = 3}

  if (max.cluster<3){
    message("Too few cluster: reset max.cluster to 3")
    max.cluster=3}


  ##############################
  ### 1 Random iterations:
  ##############################

  Sys.sleep(1)
  if (verbose){print("Start simulating random components...")}

  W_sum=c()  # Matrix storing all simulated mixing matrices
  W0_sum=c() # Matrix storing the initial demixing matrices
  A_sum=c() # Matrix storing the loading matrices
  S_sum=c() # Matrix storing the scoring matrices
  boot_id=c() # Whether a component simulated is from bootstrapped
  type_id=c()

  for (i in 1:max_iter){

    if (verbose){print(paste0("Iteration: ",i))}

    if (i>floor(max_iter/2)){
      type='deflation'
      type_id=c(type_id,"d")}
    else{
      type='parallel'
      type_id=c(type_id,"p")} # Half as deflation, half as parallel

    w.init=matrix(rnorm(pcs*pcs),pcs,pcs) # Initialization
    if (i%%2==0){     # One bootstrapped + one without
      Xb=X[bootstrap_generator(nrow(X),nbt),]
      rownames(Xb)=rownames(X)
      boot_id=c(boot_id,rep("B",pcs))
      wines.ica <- try(sipca(Xb, ncomp=pcs, mode=type, max.iter=500,w.init = w.init),silent=T)}
    else{
      Xb=X
      boot_id=c(boot_id,rep("O",pcs))
      wines.ica <- try(sipca(Xb, ncomp=pcs, mode=type, max.iter=500, w.init = w.init),silent=T)} # Label components from bootstrapped simulations

    if (class(wines.ica)!="try-error"){
      S_sum=cbind(S_sum,wines.ica$x) # Matrix storing the scoring matrices
      W_sum=rbind(W_sum,wines.ica$unmixing) # Storage of demixing matrix
      W0_sum=rbind(W0_sum,w.init) # Storage of initial demixing matrix
      A_sum=cbind(A_sum,wines.ica$loadings$X)} # Storage of loading matrix
    else{message("Simulation failed")}
}

  colnames(S_sum)=paste0("EST",1:ncol(S_sum))
  colnames(A_sum)=paste0("EST",1:ncol(A_sum))
  Stage1=list(S=S_sum,A=A_sum,boot_id=boot_id,type_id=type_id)

  if (verbose){print(paste0("First stage completed: ",pcs*max_iter," sources generated!"))}

  ##############################
  ### 2 Clustering:
  ##############################

  Sys.sleep(1)
  if (verbose){print("Start clustering random components...")}

  if (ncol(S_sum)>4000){ # Faster correlation matrix for big data
      R=bigcor(S_sum,size= 2000,method="spearman",verbose = F) # Big correlation matrix
      R=R[1:nrow(R), 1:ncol(R)]}
    else{ # Smaller data size
      R=cor(S_sum,method="spearman")}

  if (trends){dist_R_matrix=(1-R)/2} # Disimilariy Matrix for time dependent data [0,1]
  else {dist_R_matrix=1-abs(R)} # Disimilariy Matrix for no-time dependent data [0,1]

  dist_R=as.dist(dist_R_matrix)
  colnames(dist_R_matrix)=rownames(dist_R_matrix)=paste0("EST",1:ncol(S_sum))
  clusterObj <- hclust(dist_R, method="average") # Hierarchical clustering

 Stage2=list(clusterObj=clusterObj,trends=trends)
 if (verbose){print("Second stage completed: cluster object created!")}

 ##############################
 ### 3 Cluster aggregation:
 ##############################

  Sys.sleep(1)
  if (verbose){print("Start aggregating clusters...")}

  S_center_summary=list()
  A_center_summary=list()
  kurtosis_summary=list()
  total_number_summary=list()
  boot_prop_summary=list()
  divergence_summary=list()
  boot_eval_summary=list()
  defl_prop_summary=list()

  gcc_index = find_center(dist_R_matrix,rep(1,ncol(S_sum))) # Global center
  gcc_distance=dist_R_matrix[gcc_index,] # Distance of all points to global center
  nb_cluster = 2

  while (nb_cluster<=max.cluster){

    cluster <- cutree(clusterObj,nb_cluster)

    # Find cluster centers:
    cc_index=find_center(dist_R_matrix,cluster) # Which simulated estimates are cluster centers

    S_center=S_sum[,cc_index] # Cluster center scores
    A_center=A_sum[,cc_index] # Cluster center loadings

    # Find newly added clusters:
    if (nb_cluster>2){
      new_comp=which(!(colnames(S_center) %in% colnames(S_center_summary[[nb_cluster-1]]))) # The new component(s) found
      old_comp=which(colnames(S_center_summary[[nb_cluster-1]]) %in% colnames(S_center))       # The older component
      S_center=cbind(S_center[,new_comp],S_center_summary[[nb_cluster-1]][,old_comp])
      A_center=cbind(A_center[,new_comp],A_center_summary[[nb_cluster-1]][,old_comp])}
    else{
      new_comp=c(1,2)
      old_comp=c()}

    # Evaluate new clusters:

    total_number=c() # Number of total estimates in the newly added cluster
    boot_prop=c() # Proportion of bootstrap estimates in the newly added cluster
    divergence=c() # Geometry index of each cluster
    dis_boot_no_boot=c() # Compare bootstrapped and no bootstrapped estimates
    defl_prop=c()

    #for (nc in new_comp){
    for (nc in new_comp){
      total_number=c(total_number,sum(cluster==nc))
      boot_prop=c(boot_prop,sum(cluster[boot_id=="B"]==nc)/sum(cluster==nc))
      defl_prop=c(defl_prop,sum(cluster[type_id=="d"]==nc)/sum(cluster==nc))

      cc_distance=dist_R_matrix[cc_index[nc],] # Distances from members of cluster to its center
      divergence<-c(divergence,mean(cc_distance[cluster==nc])/mean(gcc_distance[cluster==nc]))

      boot_members=no_boot_members=rep(1,length(cluster))
      boot_index=intersect(which(cluster==nc),which(boot_id=="B")) # Index of bootstrapped members in the cluster
      no_boot_index=intersect(which(cluster==nc),which(boot_id=="O")) # Index of non-bootstrapped members in the cluster

      if ((length(boot_index)>10) & (length(no_boot_index)>10)){ # At least 10 to evaluate!
        boot_members[boot_index]=2
        no_boot_members[no_boot_index]=2
        boot_center=find_center(dist_R_matrix,boot_members)[2]
        no_boot_center=find_center(dist_R_matrix,no_boot_members)[2]
        mean_distance=mean(dist_R_matrix[cluster==nc,cluster==nc]) # Average distance between points in the cluster
        dis_boot_no_boot=c(dis_boot_no_boot,dist_R_matrix[boot_center,no_boot_center]/mean_distance) # Distance between boot center-no boot center
      }
      else {dis_boot_no_boot=c(dis_boot_no_boot,NA)}
    }

    if (nb_cluster>2){
    total_number=c(total_number,total_number_summary[[nb_cluster-1]][old_comp])
    boot_prop=c(boot_prop,boot_prop_summary[[nb_cluster-1]][old_comp])
    divergence=c(divergence,divergence_summary[[nb_cluster-1]][old_comp])
    dis_boot_no_boot=c(dis_boot_no_boot,boot_eval_summary[[nb_cluster-1]][old_comp])
    defl_prop=c(defl_prop,defl_prop_summary[[nb_cluster-1]][old_comp])}

    S_center_summary[[nb_cluster]]=S_center
    A_center_summary[[nb_cluster]]=A_center
    kurtosis_summary[[nb_cluster]]=kurtosis(S_center)
    total_number_summary[[nb_cluster]]=total_number
    boot_prop_summary[[nb_cluster]]=boot_prop
    divergence_summary[[nb_cluster]]=divergence
    boot_eval_summary[[nb_cluster]]=dis_boot_no_boot
    defl_prop_summary[[nb_cluster]]=defl_prop

    if (!all(total_number>30)){ # At least 30 estimates in each cluster
      max.cluster=nb_cluster
      message("Too few estimates in new cluster(s): max.cluster reset to: ",nb_cluster)}

    nb_cluster=nb_cluster+1}

  Stage3=list(max.cluster=max.cluster,S_history=S_center_summary,A_history=A_center_summary,Kurt_history=kurtosis_summary,defl_history=defl_prop_summary,
              tn_history=total_number_summary,bn_history=boot_prop_summary,dispersion_history=divergence_summary,boot_eval=boot_eval_summary)

  if (verbose){print("Third stage completed: cluster aggregated and validated!")}

return(list(Stage1=Stage1,Stage2=Stage2,Stage3=Stage3))}


bootstrap_generator <- function(n,c){

  # This function create bootstrap indices
  # n: total number of samples, in data matrix the indices of samples should be 1:n
  # c: number of samples taken for bootstrapping

  boot=1:n
  boot_chosen=sample(boot,c) # Indices chosen to be replaced
  boot_used=sample(setdiff(boot,boot_chosen),c) # Indices chosen to replace boot_chosen
  boot[boot_chosen]=boot_used

  return(boot)
}


find_center <-function(D, cluster){

  # This function calculates index of cluster centers based on distance matrix
  # cluster: vector e.g 1,2,1,3,2,2...

  nb_cluster=length(unique(cluster))
  ci_index=c()
  for (p in 1:nb_cluster){
    cl=which(cluster==p) # Indices of estimates that belong to cluster p
    D_cl=D[cl,cl] # Distance matrix for estimates belonging to this cluster
    ci=which.min(apply(D_cl,1,sum))  # Which estimates has minimal distance to other points in the clusteR
    ci_index=c(ci_index,cl[ci])}

  return(ci_index)
}

