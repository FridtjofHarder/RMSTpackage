#' Launches function convert_contrast_ph in Shiny
#'
#' Run using: app_convert_contrast_ph()
#'
#' @return A Shiny app object. Call this function to run the app.
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
        shiny::numericInput("shape", "shape", value = 1),
        shiny::numericInput("parameterisation", "parameterisation", value = 1),
        shiny::numericInput("RMSTD", "RMSTD", value = NA_real_),
        shiny::numericInput("RMSTR", "RMSTR", value = NA_real_),
        shiny::numericInput("tau", "tau", value = NA_real_),
        shiny::numericInput("HR", "hazard ratio (HR)", value = NA_real_),
        shiny::numericInput("median_diff", "median difference", value = NA_real_),
        shiny::numericInput("percentile_diff", "percentile difference", value = NA_real_),
        shiny::numericInput("percentile", "percentile", value = 50),
        shiny::numericInput("survival_diff", "survival difference at tau", value = NA_real_),
        shiny::checkboxInput("plot_curves", "Show plot", value = TRUE),
        shiny::actionButton("go", "Update")
      ),
      shiny::mainPanel(
        shiny::h4("Survival curves"),
        shiny::plotOutput("surv_plot", height = "350px")
      )
    )
  )

  server <- function(input, output, session) {
    # Tracks the last value the app wrote to each contrast field. Used by
    # clean_contrast() to distinguish app-generated values from user-typed ones.
    app_written <- shiny::reactiveValues(
      RMSTD           = NA,
      RMSTR           = NA,
      HR              = NA,
      median_diff     = NA,
      percentile_diff = NA,
      survival_diff   = NA
    )

    clean_input <- function(x) {
      if (is.null(x) || length(x) == 0L || is.na(x)) return(NULL)
      x
    }

    # Like clean_input, but also returns NULL when the value matches what the
    # app itself wrote last time (i.e. it is an output being echoed, not a
    # fresh user entry). This prevents computed contrasts from being fed back
    # as inputs on subsequent runs.
    clean_contrast <- function(x, field) {
      val <- clean_input(x)
      if (is.null(val)) return(NULL)
      last <- app_written[[field]]
      if (!is.na(last) && isTRUE(all.equal(val, last))) return(NULL)
      val
    }

    res <- shiny::eventReactive(input$go, {
      convert_contrast_ph(
        scale_trmt       = clean_input(input$scale_trmt),
        scale_ctrl       = clean_input(input$scale_ctrl),
        shape            = clean_input(input$shape),
        parameterisation = clean_input(input$parameterisation),
        RMSTD            = clean_contrast(input$RMSTD,           "RMSTD"),
        RMSTR            = clean_contrast(input$RMSTR,           "RMSTR"),
        tau              = clean_input(input$tau),
        HR               = clean_contrast(input$HR,              "HR"),
        median_diff      = clean_contrast(input$median_diff,     "median_diff"),
        percentile_diff  = clean_contrast(input$percentile_diff, "percentile_diff"),
        percentile       = clean_input(input$percentile),
        survival_diff    = clean_contrast(input$survival_diff,   "survival_diff"),
        plot_curves      = input$plot_curves
      )
    })

    output$surv_plot <- shiny::renderPlot(
      {
        shiny::req(res())
        if (!is.null(res()$plot)) {
          grDevices::replayPlot(res()$plot)
        }
      },
      res = 96
    )

    shiny::observeEvent(res(), {
      out <- res()

      # Write all computed values back to the input fields so the user can see
      # them. For contrast fields, also record the written value in app_written
      # so clean_contrast() can treat them as non-inputs on the next run.
      shiny::updateNumericInput(session, "scale_trmt",      value = out$`scale trmt`)
      shiny::updateNumericInput(session, "scale_ctrl",      value = out$`scale ctrl`)
      shiny::updateNumericInput(session, "shape",           value = out$shape)
      shiny::updateNumericInput(session, "parameterisation",value = out$parameterisation)
      shiny::updateNumericInput(session, "tau",             value = out$tau)
      shiny::updateNumericInput(session, "percentile",      value = out$percentile)

      shiny::updateNumericInput(session, "RMSTD",           value = out$RMSTD)
      app_written$RMSTD           <- out$RMSTD

      shiny::updateNumericInput(session, "RMSTR",           value = out$RMSTR)
      app_written$RMSTR           <- out$RMSTR

      shiny::updateNumericInput(session, "HR",              value = out$`hazard ratio`)
      app_written$HR              <- out$`hazard ratio`

      shiny::updateNumericInput(session, "median_diff",     value = out$`median difference`)
      app_written$median_diff     <- out$`median difference`

      shiny::updateNumericInput(session, "percentile_diff", value = out$`percentile difference`)
      app_written$percentile_diff <- out$`percentile difference`

      shiny::updateNumericInput(session, "survival_diff",   value = out$`survival difference at tau`)
      app_written$survival_diff   <- out$`survival difference at tau`
    })
  }

  shiny::shinyApp(ui, server)
}
