#' @name ClinicalTrialInterimAnalysis

#' @title Clinical Trial Interim Report

#' @description This function will run a Shiny application

#' @export



  ClinicalTrialInterimAnalysis <- function() {

  appDir <- system.file("shiny-examples", "myapp", package = "CinicalTrialInterimAnalysis")

  if (appDir == "") {

    stop("Could not find example directory. Try re-installing `CinicalTrialInterimAnalysis`.", call. = FALSE)

  }



  shiny::runApp(paste(appDir,'/app.R',sep=''), launch.browser =T,host = getOption( "127.0.0.1"))

}
