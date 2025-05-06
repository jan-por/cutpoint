#' @title Estimate cutpoints in a multivariable setting for survival data
#'
#' @description One or two cutpoints of a metric variable are estimated using
#'   either the AIC (Akaike Information Criterion) or the LRT
#'   (Likelihood-Ratio Test statistic) within a multivariable Cox proportional
#'   hazards model. These cutpoints are used to create two or three groups with
#'   different survival probabilities.
#' @description The cutpoints are estimated by dichotomising the variable of
#'   interest, which is then incorporated into the Cox regression model. The
#'   cutpoint of this variable is the value at which the AIC reaches its lowest
#'   value or the LRT statistic achieves its maximum for the corresponding
#'   Cox-regression model.
#' @description This process occurs within a multivariable framework, as other
#'   covariates and/or factors are considered during the search for the
#'   cutpoints. Cutpoints can also be estimated when the variable of interest
#'   shows a U-shaped or inverse U-shaped relationship to the hazard ratio of
#'   time-to-event data. The argument `symtail` facilitates the estimation of two
#'   cutpoints, ensuring that the two outer tails represent groups of equal size.
#' @name cp_est
#' @param cpvarname character, the name of the variable for which the cutpoints
#'   are estimated.
#' @param time character, this is the follow-up time.
#' @param event character, the status indicator, normally 0=no event, 1=event
#' @param covariates character vector with the names of the covariates and/ or
#'   factors. If no covariates are used, set `covariates = NULL`.
#' @param data a data.frame, contains the following variables:
#' * variable which is dichotomized
#' * follow-up time
#' * event (status indicator)
#' * covariates and/or cofactors
#' @param nb_of_cp numeric, number of cutpoints to be estimated (1 or 2). The
#'   default is: `nb_of_cp = 1`. The other option is `nb_of_cp = 2`.
#' @param bandwith numeric, minimum group size per group in percent of the total
#'   sample size, `bandwith` must be between 0.05 and 0.30, default is 0.1
#'   If `ushape = TRUE`, `bandwidth` must be at least 0.1.
#' @param est_type character, the method used to estimate the cutpoints. The
#'   default is 'AIC' (Akaike information criterion). The other options is 'LRT'
#'   (likelihood ratio test statistic)
#' @param cpvar_strata logical value: if `FALSE`, The dichotomised variable
#'   serves as covariate in the Cox-regression model for cutpoint determination.
#'   If `TRUE`, the dichotomised variable is included as a strata in the
#'   Cox-regression model to determine the cutpoint rather than as a covariate.
#'   Default is `FALSE`.
#' @param ushape logical value: if `TRUE`, the cutpoints are estimated under the
#'   assumtion that the spline plot shows a U-shaped form or a inverted U-shaped
#'   curve. Default is `FALSE`.
#' @param symtails logical value: if `TRUE`, the cutpoints are estimated with
#'   symmetric tails. If `nb_of_cp = 1`, symtails is set to `FALSE`. Default is
#'   `FALSE`.
#' @param dp numeric, number of decimal places the cutpoints are rounded to.
#'   Default is `dp = 2`.
#' @param plot_splines logical value: if `TRUE`, a penalized spline plot is
#'   created. Default is `TRUE`.
#' @param all_splines logical value: if `TRUE`, The plot shows splines with
#'   different degrees of freedom. This may help determine whether
#'   misspecification or overfitting occurs. Default is `TRUE`.
#' @param print_res logical value: if `TRUE` the function prints the summary
#'   of the cutpoint estimation to the console. Default is `TRUE`.
#' @param verbose logical value: if `TRUE` the function prints the approximate
#'   remaining process-time and other information to the console. If `FALSE`, no
#'   information will be printed to the console, including the summary of the
#'   cutpoint estimation. Default is `TRUE`.
#' @references Govindarajulu, U., & Tarpey, T. (2020). Optimal partitioning for
#'   the proportional hazards model. Journal of Applied Statistics, 49(4),
#'   968–987. https://doi.org/10.1080/02664763.2020.1846690
#' @returns Returns the `cpobj` object with cutpoints and the characteristics
#'   of the formed groups.
#' @examples
#' \donttest{
#' # Example 1:
#' # Estimate two cutpoints of the variable biomarker.
#' # The dataset data1 is included in this package and contains
#' # the variables time, event, biomarker, covariate_1, and covariate_2.
#' cpobj <- cp_est(
#'   cpvarname  = "biomarker",
#'   covariates = c("covariate_1", "covariate_2"),
#'   data       = data1,
#'   nb_of_cp   = 2,
#'   plot_splines = FALSE
#'   )
#'
#' # Example 2:
#' # Searching for cutpoints, if the variable shows a U-shaped or
#' # inverted U-shaped relationship to the hazard ratio.
#' # The dataset data2_ushape is included in this package and contains
#' # the variables time, event, biomarker, and cutpoint_1.
#' cpobj <- cp_est(
#'   cpvarname  = "biomarker",
#'   covariates = c("covariate_1"),
#'   data       = data2_ushape,
#'   nb_of_cp   = 2,
#'   bandwith   = 0.2,
#'   ushape     = TRUE,
#'   plot_splines = FALSE
#'   )
#'   }
#' @importFrom survival coxph
#' @importFrom survival Surv
#' @importFrom stats AIC complete.cases median quantile rnorm as.formula
#'                   update.formula
#' @importFrom utils globalVariables
#' @export
#'
#' @seealso [cp_splines_plot()] for penalized spline plots, [cp_value_plot()]
#'   for Value plots and Index plots
NULL
cp_est <- function(cpvarname,
            time         = "time",
            event        = "event",
            covariates   = NULL,
            data         = data,
            nb_of_cp     = 1,
            bandwith     = 0.1,
            est_type     = "AIC",
            cpvar_strata = FALSE,
            ushape       = FALSE,
            symtails     = FALSE,
            dp           = 2,
            plot_splines = TRUE,
            all_splines  = TRUE,
            print_res    = TRUE,
            verbose      = TRUE
           ) {

#' Verify that the input is correct
#' Check if the cutpoint variable is in the data
   if (!all(cpvarname %in% colnames(data))) {
      stop("cpvarname must be included in data, check the name of cpvarname")
   }
   if (!is.character(cpvarname)) {
      stop("cpvarname must be a character")
   }
   if (length(unique(data[ ,cpvarname])) == 1) {
      stop("cpvarname must have more than one unique value")
   }
   if (length(unique(data[ ,cpvarname])) < 3) {
       stop("cpvarname must have more than two unique values")
   }

#' Check if the time variable is in the data
   if (!is.character(time)) {
       stop("time must be a character")
   }
   if (!all(is.numeric(data[ ,time]) | is.na(data[ ,time]))) {
       stop("time must be numeric or NA")
   }
   if (length(unique(data[ ,time])) == 1) {
       stop("time must have more than one unique value")
   }
   if (length(unique(data[ ,time])) < 3) {
       stop("time must have more than two unique values")
   }

   if (any(data[ ,time] <= 0)){
       message("\nPlease note: For at least one observation the follow-up time is ",
       "\u2264"," 0\n this can lead to an error message of the Cox regression\n")
   }

   if (!is.character(event)) {
       stop("event must be a character")
   }
   if (length(unique(data[ ,event])) != 2) {
        stop("event must have two unique values")
   }

   if (!is.null(covariates) & !is.character(covariates)) {
      stop("covariates must be a character vector or NULL")
   }
   if (!all(covariates %in% colnames(data))) {
      stop("all covariates must be included in data, check the names of the covariates")
   }

   if (!is.data.frame(data)) {
      stop("data must be a data.frame")
   }

   if (!is.numeric(nb_of_cp))
      stop("nb_of_cp must be numeric")
   if (nb_of_cp != 1 &
       nb_of_cp != 2)
      stop("nb_of_cp must be 1 or 2")

   if (!is.numeric(bandwith))
      stop("bandwith must be numeric")
   if (bandwith < 0.05 |
       bandwith > 0.30)
      stop("bandwith must be between 0.05 and 0.30")

   if (!is.character(est_type))
         stop("est_type must be a character ('AIC' or 'LRT')")
      est_type <- toupper(est_type)
   if ((est_type != "AIC") && (est_type != "LRT"))
         stop("est_type must be 'AIC' or 'LRT'")

   if (!is.logical(cpvar_strata))
         stop("cpvar_strata must be logical (TRUE or FALSE)")

   if (!is.logical(ushape))
      stop("ushape must be logical (TRUE or FALSE)")

   if (!is.logical(symtails))
      stop("symtails must be logical (TRUE or FALSE)")

   if (!is.numeric(dp))
      stop("dp must be numeric")
   if (dp %% 1 != 0)
      stop("dp must be an integer")
   if (dp < 0)
      stop("dp must be 0 or greater than 0")
   if (dp > 19)
      stop("dp must be smaller than 20")

   if (!is.logical(plot_splines))
      stop("plot_splines must be logical (TRUE or FALSE)")

   if (!is.logical(all_splines))
      stop("all_splines must be logical (TRUE or FALSE)")

   if (!is.logical(print_res))
         stop("print_res must be logical (TRUE or FALSE)")

   if (!is.logical(verbose))
      stop("verbose must be logical (TRUE or FALSE)")


   # "cpvarname" is used for labelling and "cpvar" for calculations
   cpvar <- cpvarname

   #' data frame includes only the variables that should be part of the model
   cpdata <- data[, c(cpvar, time, event, covariates)]

   #' Keep the original cutpoint variable
   cpvariable_original <- sort(data[ ,cpvar])


   #' nrm = Number of rows in cpdata before possible changes
   nrm_start  <- nrow(cpdata)

   colnames(cpdata)[1:3] <- c("cpvar", "time", "event")

   #' Observations are deleted for missing values of the cutpoint variable
   cpdata <- cpdata[complete.cases(cpdata$cpvar), ]

   #' If there are more than 1000 observations, a random sample with 1000
   #'     observations is created
   sample_yes <- FALSE

   if (nrow(cpdata) > 1000) {
      cpdata <- cpdata[sample(nrow(cpdata), size=1000), ]
      sample_yes <- TRUE
   }

   #' Data frame is sorted by cpvar in ascending order
   cpdata <- cpdata[order(cpdata$cpvar), ]
   cpvar <- cpdata$cpvar


   #' nrm = Number of rows in cpdata after removing observations with
   #'     missing values in cpvar
   nrm  <- nrow(cpdata)


      #' If only one cutpoint is estimated, symtails and ushape is set to FALSE
   if (nb_of_cp == 1) {
      symtails <- FALSE
      ushape   <- FALSE
   }

   #' If ushape is TRUE and bandwith < 0.1, then bandwith is set to 0.1 and
   #' change_bw is set to TRUE
   change_bw <- FALSE
   if (ushape == TRUE && bandwith < 0.1) {
      bandwith <- 0.1
      change_bw <- TRUE    # Info for output
   }

   #' If there are no covariates, constant_var is defined as a vector of 1
   if (is.null(covariates)) { constant_var   <- rep(1, nrm) }


   # Formula (FML) for Cox-regression from vector of covariates:

   if (cpvar_strata == FALSE) {
      if (!is.null(covariates)) {
         FML <- as.formula(paste0('~ cpvariable_dicho +',
                                  paste(covariates, collapse = "+")))
      } else {
         FML <- as.formula(paste0('~ cpvariable_dicho + constant_var')) }

   # if cpvar_strata == TRUE:
   } else {

      if (!is.null(covariates)) {
         FML <- as.formula(paste0('~ strata(cpvariable_dicho) +',
                                  paste(covariates, collapse = "+")))
      } else {
         FML <- as.formula(paste0('~ strata(cpvariable_dicho) + constant_var'))
      }
   }


   #' Numbers of observations which should at least remain in row
   m.perm <- factors_combine(bandwith, nb_of_cp, nrm, symtails)


   #' If ushape is TRUE then create new m.perm with 2 categories, as u-shape only
   #'     has 2 categories
   if (ushape == TRUE) {
      m.perm[m.perm == 3] <- 1
   }


   loop_nr <- 0

   #' Number of variants of the cpvar, the number of times the loop must
   #'     be run through
   nbr.m.perm <- nrow(m.perm)

   #' Vector for AIC values
   AIC_values <- rep(NA, nbr.m.perm)

   #' Vector for Likelihood ratio test (LRT) values
   LRT_values <- rep(NA, nbr.m.perm)

   #' Matrix for cpvariable values
   cpvariable_values <- matrix(NA, nrow = nbr.m.perm, ncol = 2)


   #' After 5%, remaining time is communicated, for this timefactor is defined
   timefactor <- 0.05
   nbr.m.perm.time <- round(nbr.m.perm * timefactor, 0)

   #' Start time measurement
   ptm <- proc.time()


   for (i in 1:(nbr.m.perm)) {

      loop_nr <- loop_nr + 1

      cpvariable_dicho <- m.perm[i, ]

      result.cox <- survival::coxph(update.formula(Surv(time, event)~., FML),
                                  data = cpdata, eps = 1e-09, toler.inf = 1e-03)

      #' Extraction of AIC value
      AIC_values[i] <- AIC(result.cox)

      #' Extraction of the likelihood ratio chi2 -test
      LRT_values[i] <- (summary(result.cox))$logtest["test"]

      # Assigning the corresponding cpvariable values
      cpvariable_values[i, 1] <- cpvar[(rle(cpvariable_dicho)$lengths[1])]

      if (nb_of_cp == 2) { cpvariable_values[i, 2] <-
                           cpvar[(sum(rle(cpvariable_dicho)$lengths[1:2]))]
      }

      rm(result.cox)

      # Show the user the approx. remaining time for estimation if verbose==TRUE

      if (nbr.m.perm.time == loop_nr) {
         tm <- proc.time() - ptm

         if (verbose == TRUE) {
            message("\nApprox. remaining time for estimation: ",
                round((tm[3] * (1 / timefactor)
            ), 0), " seconds \n")
         }

         ptm <- proc.time()
      }


   } #' End: for (i in 1:nbr.m.perm)

   rm(FML)
   rm(i)
   rm(timefactor)
   rm(tm)
   rm(nbr.m.perm.time)


   #' if ushape==TRUE then create new m.perm with 3 categories, because
   #'    has only 2 categories
   if (ushape == TRUE) {
      m.perm <- factors_combine(bandwith, nb_of_cp, nrm, symtails)
   }

   #' Extract the cutpoints from the cpvar variable
   #'

   if(est_type == "AIC") {

      # Cutpoint is at position minAIC_row_nb of, all those less than or equal
      # ... to Cutpoint (minAIC_row_nb) belong to the first group
      minAIC_row_nb <- which.min(AIC_values)
      cp1_position  <- sum(m.perm[minAIC_row_nb, ] == 1)
      cp2_position  <- sum(m.perm[minAIC_row_nb, ] <= 2)

      } else {

      # Cutpoint is at position maxLRT_row_nb of , all those less than or equal
      # ...to Cutpoint (minLRT_row_nb) belong to the first group
      maxLRT_row_nb <- which.max(LRT_values)
      cp1_position  <- sum(m.perm[maxLRT_row_nb, ] == 1)
      cp2_position  <- sum(m.perm[maxLRT_row_nb, ] <= 2)
   }

   # Generate vector with cutpoints

   if(nb_of_cp == 1){
      cp <- c(NA)
      cp[1] <- cpvar[cp1_position]
      names(cp[1]) <- "CP1"

   } else {
      cp <- c(NA, NA)
      cp[1] <- cpvar[cp1_position]
      names(cp[1]) <- "CP1"
      cp[2] <- cpvar[cp2_position]
      names(cp[2]) <- "CP2"


      # check if cp[1] and cp[2] are not NA
      if (!is.na(cp[1]) && (!is.na(cp[2]))) {

         # check if cp[1] is smaller than cp[2], if yes, switch the values
         if(cp[2] < cp[1]){ cpx   <- cp[1]
                            cp[1] <- cp[2]
                            cp[2] <- cpx
                            rm(cpx)
         }
      }
   }


   returnlist <- list(
      bandwith            = bandwith,
      change_bw           = change_bw,
      cp                  = cp,
      cpdata              = cpdata,
      cpvarname           = cpvarname,
      cpvariable_original = cpvariable_original,
      cpvariable_values   = cpvariable_values,
      cpvar_strata        = cpvar_strata,
      covariates          = covariates,
      est_type            = est_type,
      nb_of_cp            = nb_of_cp,
      nrm_start           = nrm_start,
      dp                  = dp,
      sample_yes          = sample_yes,
      symtails            = symtails,
      ushape              = ushape,
      AIC_values          = AIC_values,
      LRT_values          = LRT_values
   )


   # Print summary of estimation process if verbose == TRUE:--------------------
   if (verbose == TRUE & print_res == TRUE) {
      cp_estsum (returnlist)
   }

   # Create splines plot
   if (plot_splines == TRUE) {
      cp_splines_plot(returnlist, show_splines = all_splines )
   }

   return(returnlist)

} # End: cp_est <- function
