#' @title Plot AIC and LRT-statistics values from `cpobj` object
#'
#' @description Create a plot of AIC or Likelihood ratio test statistic values
#'   for the estimation procedure. If there are two cutpoints, a Contour-plot
#'   and an Index-plot can be generated.
#' @name cp_value_plot
#' @param cpobj list, contains a vector of AIC values (AIC_values) and
#'   Likelihood ratio test statistic values (LRT_values) of the estimating
#'   procedure
#' @param plotvalues character, either `AIC` or `LRT`. Either the AIC or
#'   LRT values are displayed. Default is `AIC`.
#' @param dp.plot numeric, digits for the AIC values and LRT values.
#'   Default is `2`.
#' @param show_limit logical, if `TRUE` the minimum AIC value is shown in the
#'   plot if `plotvalues = "AIC"`, and the maximum LRT value is shown if
#'   `plotvalues = "LRT"`
#' @param plottype2cp character, either `contour` or `index`. Default is
#'   `contour`. This option is available only when searching for two cutpoints.
#'   Index plots and contour plots can be selected. Index plots display all AIC
#'   or LRT values from the estimation process as a scatter plot. Contour plots
#'   are shown in the RStudio viewer and illustrate the two potential cutpoints
#'   along with the corresponding AIC or LRT values. Index plots that do not
#'   show extreme values suggest that there may not be any actual cutpoints in
#'   the data. Contour plots provide an opportunity to explore whether there
#'   might be other potential cutpoints with similar AIC or LRT values. The
#'   smaller the `bandwidth` (minimum group size per group), the more precise
#'   and meaningful the contour plots can be interpreted.
#' @returns Plots the AIC- or LRT-values, derived from the estimation procedure.
#' @examples
#' \dontrun{
#' # Example 1
#' # Plot AIC-values and potential cutpoints of the estimation process
#'
#' # Create AIC values:
#' AIC_values <- c(1950:1910, 1910:1920, 1920:1880, 1880:1920)
#' AIC_values <- round(AIC_values + rnorm(length(AIC_values),
#'                    mean = 0, sd = 5), digits = 2)
#'
#' # Create a cutpoint variable:
#' cpvariable_values <- matrix(NA, nrow = length(AIC_values), ncol = 2)
#' cpvariable_values[ ,1] <- c(1:length(AIC_values))
#'
#' # Create a cutpoint object (cpobj):
#' cpobj <- list(AIC_values        = AIC_values,
#'               nb_of_cp          = 1,
#'               cpvariable_values = cpvariable_values,
#'               cpvarname         = "Cutpoint variable"
#'               )
#'
#' cp_value_plot(cpobj, plotvalues = "AIC", dp.plot = 2, show_limit = TRUE)
#'
#' # Example 2
#' # Splines plot based on data1
#' # The data set data1 is included in this package
#' cpobj <- cp_est(
#'   cpvarname    = "biomarker",
#'   covariates   = c("covariate_1", "covariate_2"),
#'   data         = data1,
#'   nb_of_cp     = 2,
#'   plot_splines = TRUE,
#' )
#' # Example 3
#' # Contour plot based on data1
#' # The data set data1 is included in this package
#' cpobj <- cp_est(
#'    cpvarname    = "biomarker",
#'    covariates   = c("covariate_1", "covariate_2"),
#'    data         = data1,
#'    nb_of_cp     = 2,
#'    plot_splines = FALSE,
#' )
#' cp_value_plot(cpobj, plotvalues = "AIC", plottype2cp = "contour")
#' }
#' @importFrom graphics plot title legend
#' @importFrom utils globalVariables
#' @importFrom plotly plot_ly
#' @export
#'
#' @seealso [cp_est()] for main function of the package, [cp_splines_plot()]
#'   for penalized spline plots
NULL
cp_value_plot <- function(cpobj,
                          plotvalues  = "AIC",
                          dp.plot     = 2,
                          show_limit  = TRUE,
                          plottype2cp = "contour"
                          )  {

   if (!is.list(cpobj)) {
      stop("Cutpoint object cpobj must be a list")
   }

   if (!is.character(plotvalues)) {
      stop("plotvalues must be a character")
   }
   plotvalues <- toupper(plotvalues)
   if (plotvalues != "AIC" && plotvalues != "LRT") {
      stop("plotvalues must be either 'AIC' or 'LRT'")
   }

   if (cpobj$nb_of_cp == 2){

      if (!is.character(plottype2cp)) {stop("plottype2cp must be a character")}
      plottype2cp <- tolower(plottype2cp)
      if (plottype2cp != "contour" && plottype2cp != "index") {
         stop("plottype2cp must be either 'contour' or 'index'")
      }
   }

   if (cpobj$nb_of_cp != 1 && cpobj$nb_of_cp != 2) {
      stop("nb_of_cp must be either 1 or 2")
   }
   if (cpobj$nb_of_cp == 1) {
      if (!is.numeric(cpobj$cpvariable_values[ ,1])) {
         stop("cpvariable_values must be numeric")}
      if (length(cpobj$cpvariable_values[ ,1]) == 0) {
         stop("cpvariable_values must not be empty")}
      if (length(cpobj$cpvariable_values[ ,1]) < 2) {
         stop("cpvariable_values must have at least 2 values")}
   }
   if (cpobj$nb_of_cp == 2) {
      if (!is.numeric(cpobj$cpvariable_values[ ,1])) {
         stop("cpvariable_values must be numeric")}
      if (length(cpobj$cpvariable_values[ ,1]) == 0) {
         stop("cpvariable_values must not be empty")}
      if (length(cpobj$cpvariable_values[ ,1]) < 2) {
         stop("cpvariable_values must have at least 2 values")}
      if (!is.numeric(cpobj$cpvariable_values[ ,2])) {
         stop("cpvariable_values must be numeric")}
      if (length(cpobj$cpvariable_values[ ,2]) == 0) {
         stop("cpvariable_values must not be empty")}
      if (length(cpobj$cpvariable_values[ ,2]) < 2) {
         stop("cpvariable_values must have at least 2 values")}
   }

   if (!is.numeric(dp.plot))
      stop("dp.plot must be numeric")
   if (dp.plot %% 1 != 0)
      stop("dp.plot must be an integer")
   if (dp.plot < 0)
      stop("dp.plot should be 0 or greater than 0")
   if (dp.plot > 19)
      stop("dp.plot must be smaller than 20")


   if (!is.logical(show_limit))
      stop("show_limit must be logical (TRUE or FALSE)")


   if (!is.character(cpobj$cpvarname)) {
      stop("cpvarname must be a character")
   }
   cpvarname <- cpobj$cpvarname


   #' Checks and assignments if plotvalues == AIC
   if (plotvalues == "AIC") {

         if (!is.vector(cpobj$AIC_values)) {
            stop("AIC values must be a vector")}
         if (length(cpobj$AIC_values) == 0) {
            stop("AIC values must not be empty")}
         if (length(cpobj$AIC_values) < 2) {
            stop("AIC values must have at least 2 plot_values")}
         if (any(cpobj$AIC_values < 0)) {
            stop("AIC values must be greater than 0")}
         if (!is.numeric(cpobj$AIC_values)) {
            stop("AIC values must be numeric")}

         plot_values       <- round(cpobj$AIC_values, digits = dp.plot)
         value_name        <- "AIC"
         extrem_value      <- min(plot_values)
         extrem_value_name <- "min"
         plot_color        <- "#6a93b0"
         legend_position   <- "topright"

      } else {

   #' Checks and assignments if plotvalues == LRT
         if (!is.vector(cpobj$LRT_values)) {
            stop("LRT values must be a vector")}
         if (length(cpobj$LRT_values) == 0) {
            stop("LRT values must not be empty")}
         if (length(cpobj$LRT_values) < 2) {
            stop("LRT values must have at least 2 values")}
         if (any(cpobj$LRT_values < 0)) {
            stop("LRT values must be greater than 0")}
         if (!is.numeric(cpobj$LRT_values)) {
            stop("LRT values must be numeric")}

         plot_values       <- round(cpobj$LRT_values, digits = dp.plot)
         value_name        <- "LRT"
         extrem_value      <- max(plot_values)
         extrem_value_name <- "max"
         plot_color        <- "blue"
         legend_position   <- "bottomright"
      }


   #' Create plot if number of cutpoints is 1
   if(cpobj$nb_of_cp == 1) {

      plot(cpobj$cpvariable_values[ ,1], plot_values,
           col  = plot_color, pch = 19, cex = 0.5,
           xlab = paste0(cpvarname),
           ylab = paste0(value_name, " - values"),
           main = paste0("Plot of ", value_name,
                         " - values of the estimation process"))

      if (show_limit == TRUE) {

         legend(legend_position, legend =
                   paste0(extrem_value_name, " ", value_name, ": ",
                          extrem_value
                   ),
                cex = 0.75
         )
      }
   }

   #' Create plot if number of cutpoints is 2 and plottype2cp == "index"
   if(cpobj$nb_of_cp == 2 && plottype2cp == "index") {

      plot(plot_values, col = plot_color, pch = 19, cex = 0.5,
           xlab = "Index",
           ylab = paste0(value_name, " - values"),
           main = paste0("Plot of ", value_name,
                         " - values of the estimation process"))

      if (show_limit == TRUE) {

         legend(legend_position, legend =
                   paste0(extrem_value_name, " ", value_name, ": ",
                          extrem_value
                   ),
                cex = 0.75
         )
      }
   }


   #' Create plot if number of cutpoints is 2 and plottype2cp == "index"
   if(cpobj$nb_of_cp == 2 && plottype2cp == "contour") {

      cp_plot <- plotly::plot_ly(
         x = cpobj$cpvariable_values[ ,1],
         y = cpobj$cpvariable_values[ ,2],
         z = plot_values,
         type = "contour"
      )

      # Define X-axis title
      xaxis_title <- paste0(
         "Potential values of '",cpvarname,"' for cutpoint 1\nhigh ", value_name,
         " values: light colors; low ", value_name,
         " values: dark colours; ",
         extrem_value_name, " ", value_name, ": ",
         extrem_value
      )

      # Define Y-axis title
      yaxis_title <- paste0("Potential values of '", cpvarname,"' for cutpoint 2")

      cp_plot <- cp_plot %>%
         plotly::layout(
         title  = paste0(
            "Contour plot - coloured ", value_name,
            " value areas for potential cutpoints"
            ),
         yaxis  = list(title = yaxis_title),
         xaxis  = list(title = xaxis_title),
         legend = list(title=value_name)
         )

      print(cp_plot)
   }

   return(invisible())
}
