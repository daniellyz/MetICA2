#' Bacteria exo-metabolome example data
#'
#' Exo-metabolomic profiles of two wine bacteria strains BETA and VP41. The experiment included 6 different growth media and 3 biological replicates. Samples were taken at different day points throughout bacterial fermentation and measured on maXis ESI–QTOF-MS coupled with reversed-phase LC.
#'
#' @docType data
#'
#' @usage data(bacteria_peptides)
#'
#' @format A two-element list object:
#' \itemize{
#'   \item features: metabolic features detected and extracted. Each feature is represented by a unique ID, a mass (m/z) and a retention time in minute (RT)
#'   \item X: data matrix with 159 rows (samples) and 8332 columns (metabolic features). This is the example data matrix for MetICA simulations.
#' }
#'
#' @references Youzhong Liu, Sara Forcisi, Marianna Lucio, Mourad Harir, Florian Bahut, Magali Deleris-Bou, Sibylle Krieger-Weber, Regis D. Gougeon, Herve Alexandre and Philippe Schmitt-Kopplin, Digging into the low molecular weight peptidome with the OligoNet web server, Scientific reports (2017) Vol. 7 no. 11692
#'
"bacteria_peptides"

