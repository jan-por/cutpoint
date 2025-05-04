#' @title Summarise cutpoint estimation
#' @description Writes the summary of the cutpoint estimation to the console.
#' @name cp_estsum
#' @param cpobj list, contains variables for `cp_estsum` function
#' @param verbose logical value: if `TRUE` the summary of the cutpoint
#'   estimation is writing to the console. Default is `TRUE`.
#' @returns Summary of the cutpoint estimation.
#' @examples
#' # Example
#' # Writes the summary to the console
#' # The data set data1 is included in this package
#' cpobj <- cp_est(
#'   cpvarname    = "biomarker",
#'   covariates   = c("covariate_1", "covariate_2"),
#'   data         = data1,
#'   nb_of_cp     = 2,
#'   plot_splines = FALSE,
#'   verbose      = FALSE
#' )
#' cp_estsum <- function(cpobj, verbose = TRUE)
#' @export
#'
#' @seealso [cp_est()] for main function of the package.
NULL
cp_estsum <- function( cpobj,
                       verbose = TRUE){


   #' Check if cpobj is a list
   if (!is.list(cpobj)) {
      stop("Cutpoint object (cpobj) for function cp_estsum must be a list")
   }

   if (!is.logical(verbose))
      stop("verbose for function cp_estsum must be logical (TRUE or FALSE)")

   #' Check cpobj - variables

   if (!is.numeric(cpobj$bandwith))
      stop("bandwith must be numeric for function cp_estsum")
   if (cpobj$bandwith < 0.05 | cpobj$bandwith > 0.3)
      stop("bandwith must be between 0.05 and 0.3 for function cp_estsum")

   if (!is.logical(cpobj$change_bw))
      stop("change_bw must be logical (TRUE or FALSE) for function cp_estsum")

   if (!is.numeric(cpobj$cp)) stop("cp must be numeric for function cp_estsum")

   if (!is.character(cpobj$cpvarname)) {
      stop("cpvarname must be a character for function cp_estsum")
   }


   if (!is.vector(cpobj$cpvariable_original)) {
      stop("cpvariable_original must be a vector for function cp_estsum")}
   if (length(cpobj$cpvariable_original) == 0) {
      stop("cpvariable_original must not be empty for function cp_estsum")}
   if (length(cpobj$cpvariable_original) < 2) {
      stop("cpvariable_original must have at least 2 values for fun cp_estsum")}
   if (!is.numeric(cpobj$cpvariable_original)) {
      stop("cpvariable_original must be numeric for function cp_estsum")}

   if (!is.logical(cpobj$cpvar_strata))
      stop("cpvar_strata for fun cp_estsum must be logical (TRUE or FALSE)")

   if (!is.null(cpobj$covariates) & !is.character(cpobj$covariates)) {
      stop("covariates must be a character vector or NULL for fun cp_estsum")
   }

   if (!is.character(cpobj$est_type))
      stop("est_type must be a character ('AIC' or 'LRT') for fun cp_estsum")
   cpobj$est_type <- toupper(cpobj$est_type)
   if ((cpobj$est_type != "AIC") && (cpobj$est_type != "LRT"))
      stop("est_type must be 'AIC' or 'LRT' for function cp_estsum")

   if (!is.numeric(cpobj$nb_of_cp))
      stop("nb_of_cp must be numeric for function cp_estsum")
   if (cpobj$nb_of_cp != 1 &
       cpobj$nb_of_cp != 2)
      stop("nb_of_cp must be 1 or 2 for function cp_estsum")

   if (!is.numeric(cpobj$nrm_start))
      stop("nrm_start must be numeric for function cp_estsum")

   if (!is.logical(cpobj$sample_yes))
      stop("sample_yes for function cp_estsum must be logical (TRUE or FALSE)")

   if (!is.logical(cpobj$symtails))
      stop("symtails for function cp_estsum must be logical (TRUE or FALSE)")

   if (!is.logical(cpobj$ushape))
      stop("ushape for function cp_estsum must be logical (TRUE or FALSE)")


   #' Extracting necessary variables from cpobj
   bandwith            <- cpobj$bandwith
   change_bw           <- cpobj$change_bw
   cp                  <- cpobj$cp
   cpvarname           <- cpobj$cpvarname
   cpvariable_original <- cpobj$cpvariable_original
   cpvar_strata        <- cpobj$cpvar_strata
   covariates          <- cpobj$covariates
   est_type            <- cpobj$est_type
   nb_of_cp            <- cpobj$nb_of_cp
   nrm_start           <- cpobj$nrm_start
   sample_yes          <- cpobj$sample_yes
   symtails            <- cpobj$symtails
   ushape              <- cpobj$ushape


   #' Get the counts and percentage of groups in relation to the original data
   lcpvo <- length(cpvariable_original)

   if(nb_of_cp == 1){
      x            <- which(cpvariable_original == cp[1])
      nbbygroup1   <- x[length(x)]
      percbygroup1 <- nbbygroup1 / lcpvo
      percbygroup1 <- round(percbygroup1*100,1)
      nbbygroup2   <- lcpvo - nbbygroup1
      percbygroup2 <- nbbygroup2 / lcpvo
      percbygroup2 <- round(percbygroup2*100,1)
      rm(x)
   } else {
      x            <- which(cpvariable_original == cp[1])
      nbbygroup1   <- x[length(x)]
      percbygroup1 <- nbbygroup1 / lcpvo
      percbygroup1 <- round(percbygroup1*100,1)
      x            <- which(cpvariable_original == cp[2])
      nbbygroup2   <- x[length(x)] - nbbygroup1
      percbygroup2 <- nbbygroup2 / lcpvo
      percbygroup2 <- round(percbygroup2*100,1)
      nbbygroup3   <- lcpvo - nbbygroup1 - nbbygroup2
      percbygroup3 <- nbbygroup3 / lcpvo
      percbygroup3 <- round(percbygroup3*100, 1)
      rm(x)
   }

   min_gr_size   <- round((nrm_start * bandwith), 0)
   bandwith_perc <- bandwith * 100

   #' Writing summary to the console if verbose == TRUE ------------------------

   if (verbose == TRUE){

   if (sample_yes == TRUE){
      cat("--------------------------------------------------------------------\n")
      cat("! Because the number of observations in the original dataset is >1000\n")
      cat("  a random sample of 1000 observations is used for estimating the cutpoints\n")
   }
   cat("--------------------------------------------------------------------\n")
   cat("SETTINGS:\n")
   cat(" Cutpoint-variable                    = ", cpvarname, "\n")
   cat(" Number of cutpoints   (nb_of_cp)     = ", nb_of_cp, "\n")
   if (change_bw == TRUE) {
      cat(" Min. group size in %  (bandwith) = ", bandwith,
          " (was set to 0.1 because ushape is TRUE)\n")
   } else {
      cat(" Min. group size in %  (bandwith)     = ", bandwith, "\n")}
   cat(" Estimation type       (est_type)     = ", est_type, "\n")
   cat(" CP-variable as strata (cpvar_strata) = ", cpvar_strata, "\n")
   cat(" Symmetric tails       (symtails)     = ", symtails,
       "  (is set to FALSE if nb_of_cp = 1)\n")
   cat(" Cutpoints for u-shape (ushape)       = ", ushape,
       "  (is set to FALSE if nb_of_cp = 1)\n")
   cat("--------------------------------------------------------------------\n")
   if (is.null(covariates)) {
      cat("No covariates were selected\n")
   } else {
      cat("Covariates or factors are:\n")
      cat(" ",covariates, "\n") }
   cat("--------------------------------------------------------------------\n")
   cat("Minimum group size is ", min_gr_size,
       " (", bandwith_perc, "% of sample size in original dataset, N = ",
       nrm_start,  ")\n", sep = "")
   cat("--------------------------------------------------------------------\n")
   cat("Number of cutpoints searching for:", nb_of_cp, "\n")

   if (nb_of_cp == 1) {

      cat("Cutpoint:", cpvarname, "\u2264", cp[1], "\n")
      cat("-----------------------------------------------------------------\n")
      cat("Group size in relation to valid data of",cpvarname ,
          "in original data set\n")
      cat(" Total:   N = ", lcpvo, " (100%)\n", sep = "")
      cat(" Group A: n = ", nbbygroup1, " (", percbygroup1,
          "%)\n", sep = "")
      cat(" Group B: n = ", nbbygroup2, " (", percbygroup2,
          "%)\n", sep = "")
   }

   if (nb_of_cp == 2) {

      cat(" 1.Cutpoint:", cpvarname, "\u2264", cp[1], "\n")
      cat(" 2.Cutpoint:", cpvarname, "\u2264", cp[2], "\n")
      cat("-----------------------------------------------------------------\n")
      cat("Group size in relation to valid data of", cpvarname ,
          "in original data set:\n")

      if(ushape == FALSE) {

         cat(" Total:   N = ", lcpvo, " (100%)\n", sep = "")
         cat(" Group A: n = ", nbbygroup1, " (",
             percbygroup1, "%)\n", sep = "")
         cat(" Group B: n = ", nbbygroup2, " (", percbygroup2, "%)\n", sep = "")
         cat(" Group C: n = ", nbbygroup3, " (", percbygroup3, "%)\n", sep = "")

      } else {

         cat(" Total:                N = ", lcpvo, " (100%)\n", sep = "")
         cat(" Group A (lower part): n = ", nbbygroup1, " (",
             percbygroup1, "%)\n", sep = "")
         cat(" Group B:              n = ", nbbygroup2, " (",
             percbygroup2, "%)\n", sep = "")
         cat(" Group A (upper part): n = ", nbbygroup3, " (",
             percbygroup3, "%)\n", sep = "")
      }
   }

   # End: Writing summary to console--------------------------------------------

} # End: if (verbose == TRUE)

   return(invisible())

} # End: cp_estsum
