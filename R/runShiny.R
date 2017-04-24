#' @export
runShiny <- function() {
  appDir <- system.file("shiny-examples", "myapp", package = "ANOVAreplication")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `ANOVAreplication`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

