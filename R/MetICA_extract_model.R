#' Extracting MetICA model
#'
#' This function extracts a MetICA model with user-decided component numbers.
#'
#' @param M1 The entire list object generated from the function MetICA(X, pcs = 15...)
#' @param ics Number of clusters (MetICA components) decided by users according to validationPlot
#' @param tops Integer or NULL object. If integer, the tops "best" components based on each evaluation criterion will be plotted.
#'
#' @return MetICA model (a list object) whose number of components is decided by users.
#' \itemize{
#'  \item{"S"}{Scores of MetICA components. Data matrix contains n rows (samples) and ics columns (components)}
#'  \item{"A"}{Loadings of MetICA components. Data matrix contains p rows (metabolic features) and ics columns (components)}
#'  \item{"eval$kurtosis"}{Vector of length ics containing Kurtosis measure of each MetICA component.}
#'  \item{"eval$cluster_size"}{Number of IPCA estimates in each corresponding cluster.}
#'  \item{"eval$divergence_index"}{Divergence indexes of clusters.}
#'  \item{"eval$boot_prop"}{Proportion of bootstrap estimates in each cluster.}
#'  \item{"eval$boot_stability"}{Stability of each cluster against bootstrapping.}
#' }
#'
#' @export
#'
MetICA_extract_model<-function(M1,ics,tops=ics){

  t=as.numeric(Sys.time())
  set.seed((t - floor(t))*1e8->seed)
  memory.size(max=2000000000)
  options(warn=-1)

  # Check function inpus:
  S_sum=M1$Stage1$S
  trends=M1$Stage2$trends
  max.cluster=M1$Stage3$max.cluster
  boot_id=M1$Stage1$boot_id

  if (is.null(S_sum)){
    stop("Please enter a valid MetICA object !")}

  if (missing(ics)){
    stop("Please enter the selected number of clusters according to validationPlot !")}

  if (ics>max.cluster || ics<2){
    stop("Number of clusters should between 2 and ", M1$Stage3$max.cluster)}

  if (is.numeric(tops)){
    if (tops>ics){tops=ics}}

  # Extract the selected model:

  M2=list()

  M2$S=M1$Stage3$S_history[[ics]]
  AAA=M1$Stage3$A_history[[ics]]

  M2$A=AAA
  norm=apply(AAA,1,function(x) sqrt(sum(x^2)))
  for (i in 1:ncol(AAA)){AAA[,i]=AAA[,i]/norm}
  M2$A1=AAA

  colnames(M2$S)=colnames(M2$A)=colnames(M2$A1)=paste0("IC",1:ics)

  kurtosis=M1$Stage3$Kurt_history[[ics]]
  total_number=M1$Stage3$tn_history[[ics]]
  boot_prop=M1$Stage3$bn_history[[ics]]
  divergence=M1$Stage3$dispersion_history[[ics]]
  dis_boot_no_boot=M1$Stage3$boot_eval[[ics]]

  names(kurtosis)=names(total_number)=names(boot_prop)=names(divergence)=names(dis_boot_no_boot)=paste0("IC",1:ics)
  M2$eval=list(kurtosis=kurtosis,cluster_size=total_number,divergence_index=divergence,boot_prop=boot_prop,boot_stability=dis_boot_no_boot)

  # Visualize and rank each cluster:

  if (!is.null(tops)){
    par(mfrow=c(1,1))

    kurtosis2=kurtosis[order(kurtosis,decreasing = T)[1:tops]]
    barplot(kurtosis2,font = 2,ylab="Component kurtosis", cex.lab = 1.3)

    divergence2=divergence[order(divergence)[1:tops]]
    barplot(divergence2,font = 2,ylab="Cluster dispersion index", cex.lab = 1.3)

    boot_prop2=boot_prop[order(abs(boot_prop-0.5))[1:tops]]
    barplot(boot_prop2,font = 2,ylab="Proportion of bootstrap estimates", cex.lab = 1.3)

    boot_stability2=dis_boot_no_boot[order(dis_boot_no_boot)[1:tops]]
    barplot(boot_stability2,font = 2,ylab="Distance bootstrap-non-bootstrap", cex.lab = 1.3)}

  return(M2)
}
