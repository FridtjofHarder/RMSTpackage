#' Compare sample size calculations
#'
#' @param input1 first input
#' @param input2 second input
#'
#' @return sum of input1 and input2
#' @export
#'
#' @examples
#' sum_of_inputs <- compare_ss_calculations(1, 42)
compare_ss_calculations <- function(input1, input2){
  print("this function will later allow comparing sample sizes calculated by different methods.
        for now it only ouputs input+input2")
  return(input1+input2)

}
