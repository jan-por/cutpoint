#' @title Plot AIC values
#'
#' @description Create plot of AIC values of estimating procedure
#' @name aicplot
#' @param cpobj list, contains a vector of AIC values (AIC_values) of the
#'   estimating procedure
#' @param dp.aic numeric, digits for the AIC values
#' @returns returns a plot which shows the AIC values of the estimating
#'   procedure
#' @examples
#' # Create a vector with AIC values
#' AIC_values <- c(1950:1910, 1910:1920, 1920:1880, 1880:1920)
#' AIC_values <- round(AIC_values + rnorm(length(AIC_values),
#'                     mean = 0, sd = 5), digits = 2)
#'
#' cpobj <- list(AIC_values = AIC_values)
#' aicplot(cpobj, dp.aic = 2)
#' @importFrom graphics plot title legend
#' @importFrom utils globalVariables
#' @export
#'
#' @seealso \code{\link{est.cutpoint}}

aicplot <-
   function(cpobj, dp.aic = 2) {

      if (!is.list(cpobj)) {
         stop("Cutpoint object cpobj must be a list")
      }

      AIC_values <- cpobj$AIC_values

      if (!is.vector(AIC_values)) {
         stop("AIC_values must be a vector")
      }

      if (!is.numeric(dp.aic))
         stop("dp.aic must be numeric")
      if (dp.aic %% 1 != 0)
         stop("dp.aic must be an integer")
      if (dp.aic < 0)
         stop("dp.aic should be 0 or greater than 0")
      if (dp.aic > 19)
         stop("dp.aic must be smaller than 20")

      AIC_values   <- round(AIC_values, digits = dp.aic)
      smallest_AIC <- min(AIC_values)

      plot(AIC_values)
      title(main = "AIC-values of the estimating prozess")

      legend("bottomleft",
             legend = paste("min AIC: ", smallest_AIC),
             cex = 0.75)

      return(invisible())
   }
