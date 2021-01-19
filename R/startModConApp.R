# Run shiny app
startVarConApp <- function() {
  appDir <- system.file("extdata","app.R", package = "ModCon")
  if (appDir == "") {
    stop("Could not find the shiny directory.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}