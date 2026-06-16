#' Creates report for function convert_contrast_ph
#'
#' Requires rmarkdown.
#'
#' @param x A list created by the function \code{convert_contrast_ph()}.
#' @param output_file Output file name.
#' @param output_dir Output directory.
#'
#' @returns A character string with the path to the generated PDF report, returned invisibly.
#' @export
#'
#' @examples
#' # convert contrasts using convert_contrast_ph() and create a report on its results
#' if (requireNamespace("rmarkdown", quietly = TRUE)) {
#'   results <- convert_contrast_ph(scale_trmt = 10, scale_ctrl = 6, tau = 4,
#'                                  percentile = 80, plot_curves = TRUE)
#'   report_convert_contrasts_ph(results)
#' }
report_convert_contrasts_ph <- function(x, output_file = "report_convert_contrasts_ph.pdf", output_dir = getwd()) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' must be installed.")
  }

  template <- system.file("reports", "report_convert_contrasts_ph.Rmd", package = "RMSTpackage")

  if (template == "") {
    stop("Report template not found in package.")
  }

  rmarkdown::render(
    input = template,
    output_file = output_file,
    output_dir = output_dir,
    params = list(result = x),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )
}
