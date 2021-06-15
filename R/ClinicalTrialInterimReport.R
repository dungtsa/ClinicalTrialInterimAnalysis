
#' @name ClinicalTrialInterimReport

#' @title Clinical Trial Interim Report

#' @description This function will run a Shiny application

#' @export



ClinicalTrialInterimReport <- function() {

  appDir <- system.file("shiny-examples", "myapp", package = "StatisticalClinicalSignificance")

  if (appDir == "") {

    stop("Could not find example directory. Try re-installing `StatisticalClinicalSignificance`.", call. = FALSE)

  }



  shiny::runApp(paste(appDir,'/app.R',sep=''), launch.browser =T,host = getOption( "127.0.0.1"))

}
