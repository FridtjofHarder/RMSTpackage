#' Converts contrast to RMSTD or RMSTR
#'
#' @param HR A scalar specifying the hazard ratio
#'
#' @return A scalaer specifying either RMSTD or RMSTR
#' @export
#'
#' @examples
#' RMSTD <- convert_contrast(2)
convert_contrast <- function(x){
  print("this will later allow converting a contrast to RMSTR or RMSTD. For now, it only prints a test number and return a test number, which is double the input")
  print(2*x)
  return(2*x)
}
