#' Convert contrasts
#'
#' @param HR A hazard ratio.
#' @param scale A scale parameter.
#'
#' @return RMST difference.
#' @export
#'
#' @examples
#' convert_contrasts(HR=42, scale=1)
convert_contrasts <- function(HR, scale){
  print("hello world")
  RMSTD <- HR*scale+10
  return(RMSTD)
}
