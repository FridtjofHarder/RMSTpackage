#' Launches function convert_contrast_ph in Shiny
#'
#' Run using: app_convert_contrast_ph()
#'
#' @returns A Shiny app object. Call this function to run the app.
#' @export
#' @examples
#' if (interactive()) {
#'   app_convert_contrast_ph()
#' }
app_convert_contrast_ph <- function() {

  ui <- shiny::fluidPage(
    shiny::titlePanel("convert_contrast_ph"),

    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::h4("Input / output parameters"),

        shiny::numericInput("scale_trmt", "scale trmt", value = NA_real_),
        shiny::numericInput("scale_ctrl", "scale ctrl", value = NA_real_),
        shiny::numericInput("shape",      "shape",      value = 1),
        shiny::numericInput("parameterization", "parameterization", value = 1),

        shiny::numericInput("RMSTD", "RMSTD", value = NA_real_),
        shiny::numericInput("RMSTR", "RMSTR", value = NA_real_),

        shiny::numericInput("tau",  "tau",  value = NA_real_),
        shiny::numericInput("HR",   "hazard ratio (HR)", value = NA_real_),

        shiny::numericInput("median_diff",     "median difference",     value = NA_real_),
        shiny::numericInput("percentile_diff", "percentile difference", value = NA_real_),
        shiny::numericInput("percentile",      "percentile",            value = 50),
        shiny::numericInput("survival_diff",   "survival difference at tau", value = NA_real_),

        shiny::checkboxInput("plot_curves", "Show plot", value = TRUE),

        shiny::actionButton("go", "Update")
      ),

      shiny::mainPanel(
        shiny::h4("Results"),
        shiny::verbatimTextOutput("results"),

        shiny::h4("Survival curves"),
        shiny::plotOutput("surv_plot", height = "350px")
      )
    )
  )

  server <- function(input, output, session) {

    clean_input <- function(x) {
      if (is.null(x)) return(NULL)
      if (length(x) == 0L) return(NULL)
      if (is.na(x)) return(NULL)
      x
    }

    res <- shiny::eventReactive(input$go, {
      convert_contrast_ph(
        scale_trmt       = clean_input(input$scale_trmt),
        scale_ctrl       = clean_input(input$scale_ctrl),
        shape            = clean_input(input$shape),
        parameterization = clean_input(input$parameterization),
        RMSTD            = clean_input(input$RMSTD),
        RMSTR            = clean_input(input$RMSTR),
        tau              = clean_input(input$tau),
        HR               = clean_input(input$HR),
        median_diff      = clean_input(input$median_diff),
        percentile_diff  = clean_input(input$percentile_diff),
        percentile       = clean_input(input$percentile),
        survival_diff    = clean_input(input$survival_diff),
        plot_curves      = input$plot_curves
      )
    })

    output$results <- shiny::renderPrint({
      shiny::req(res())
      tmp <- res()
      if (!is.null(tmp$plot)) tmp$plot <- NULL
      tmp
    })

    output$surv_plot <- shiny::renderPlot({
      shiny::req(res())
      if (!is.null(res()$plot)) {
        grDevices::replayPlot(res()$plot)
      }
    }, res = 96)

    shiny::observeEvent(res(), {
      out <- res()

      if (!is.null(out$`scale trmt`)) {
        shiny::updateNumericInput(session, "scale_trmt",
                                  value = out$`scale trmt`
        )
      }
      if (!is.null(out$`scale ctrl`)) {
        shiny::updateNumericInput(session, "scale_ctrl",
                                  value = out$`scale ctrl`
        )
      }
      if (!is.null(out$shape)) {
        shiny::updateNumericInput(session, "shape",
                                  value = out$shape
        )
      }
      if (!is.null(out$parameterization)) {
        shiny::updateNumericInput(session, "parameterization",
                                  value = out$parameterization
        )
      }
      if (!is.null(out$RMSTD)) {
        shiny::updateNumericInput(session, "RMSTD",
                                  value = out$RMSTD
        )
      }
      if (!is.null(out$RMSTR)) {
        shiny::updateNumericInput(session, "RMSTR",
                                  value = out$RMSTR
        )
      }
      if (!is.null(out$tau)) {
        shiny::updateNumericInput(session, "tau",
                                  value = out$tau
        )
      }
      if (!is.null(out$`hazard ratio`)) {
        shiny::updateNumericInput(session, "HR",
                                  value = out$`hazard ratio`
        )
      }
      if (!is.null(out$`median difference`)) {
        shiny::updateNumericInput(session, "median_diff",
                                  value = out$`median difference`
        )
      }
      if (!is.null(out$`percentile difference`)) {
        shiny::updateNumericInput(session, "percentile_diff",
                                  value = out$`percentile difference`
        )
      }
      if (!is.null(out$percentile)) {
        shiny::updateNumericInput(session, "percentile",
                                  value = out$percentile
        )
      }
      if (!is.null(out$`survival difference at tau`)) {
        shiny::updateNumericInput(session, "survival_diff",
                                  value = out$`survival difference at tau`
        )
      }
    })
  }

  shiny::shinyApp(ui, server)
}
