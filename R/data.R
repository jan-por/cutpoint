#' Dataset for testing the cutpoint estimating function: cp_est
#'
#' A dataset containing data for testing the estimating of one or two cutpoints
#'
#' @format ## "data1"
#' A data frame with 100 rows and 5 variables:
#' \describe{
#'   \item{biomarker}{numeric from 1 to 257}
#'   \item{covariate_1}{numeric, from 4.25 to 12.33, with effect of the cutpoint
#'    of the biomarker}
#'   \item{covariate_2}{numeric, from 465 to 1205, with no or small effect of
#'   the cutpoint of the biomarker}
#'   \item{time}{numeric, from 3 to 328}
#'   \item{event}{numeric, 0 or 1}
#' }
#' @source Self-generated example data
#' @name data1
#' @usage data(data1)
#' @docType data
#' @author Jan Porthun
#' @examples
#' data(data1)
#'
"data1"
NULL


#' Dataset for testing the ushape argument of cp_est function
#'
#' A dataset containing data for testing the ushape argument of cp_est function.
#'
#' @format ## "data2_ushape"
#' A data frame with 200 rows and 4 variables:
#' \describe{
#'   \item{biomarker}{numeric from 1e-04 to 4.7}
#'   \item{covariate_1}{numeric, from 8.07e-05 to 1.90}
#'   \item{time}{numeric, from 0.002 to 5.09}
#'   \item{event}{numeric, 0 or 1}
#' }
#' @source Self-generated example data
#' @name data2_ushape
#' @usage data(data2_ushape)
#' @docType data
#' @author Jan Porthun
#' @examples
#' data(data2_ushape)
#'
"data2_ushape"
NULL

