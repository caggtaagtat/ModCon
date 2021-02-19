# Run shiny app
startModConApp <- function() {
  appDir <- system.file("shinyApp","app.R", package = "ModCon")
  if (appDir == "") {
    stop("Could not find the shiny directory.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}

