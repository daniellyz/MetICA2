#' Yeast exo-metabolome example data
#'
#' Exo-metabolomic profiles of 15 yeast commercial strains for winemaking. The experiment included 3 biological replicates. Samples were taken at the end of alcoholic fermentationand measured on direct-infusion FT-ICR-MS.
#'
#' @docType data
#'
#' @usage data(yeast_metabolome)
#'
#' @format A two-element list object:
#' \itemize{
#'   \item features: metabolic features detected and extracted. Each feature is represented by a unique ID and a mass (m/z).
#'   \item X: data matrix with 45 rows (samples) and 2700 columns (metabolic features). This is the example data matrix for MetICA simulations.
#' }
#'
#' @references Youzhong Liu, Sara Forcisi, Mourad Harir, Magali Deleris-Bou, Sibylle Krieger-Weber, Marianna Lucio, Cedric Longin, Claudine Degueurce, Regis D. Gougeon, Philippe Schmitt-Kopplin and Herve Alexandre, New molecular evidence of wine yeast-bacteria interaction unraveled by non-targeted exometabolomic profiling, Metabolomics (2016) Vol. 12 no. 69

"yeast_metabolome"
