#' @title Estimates cutpoints
#'
#' @description One or two cutpoints of a variable are estimated by the AIC
#'   criterion in a Cox proportional hazards model. The cutpoints are estimated
#'   by dichotomizing the variable. The cutpoints with the lowest AIC value are
#'   chosen.
#' @name est_cutpoint
#'
#' @param cpvarname character, name of the variable for which the cutpoints are
#'   estimated
#' @param time character, this is the follow-up time
#' @param event character, the status indicator, normally 0=no event, 1=event
#' @param covariates character vector with the names of the covariates or/ and
#'   factors. If there are no covariates set: covariates = NULL
#' @param data a data.frame, contains the following variables: variable which is
#'   dichotomized, follow-up time, event (status indicator) and the covariates
#' @param nb_of_cp numeric, number of cutpoints to be estimated
#' @param bandwith numeric, minimum group size per group in percent of the total
#'   sample size, bandwith must be between 0.05 and 0.30, default is 0.10
#' @param ushape logical value: if TRUE, the cutpoints are estimated under the
#'   assumtion that the spline plot shows a u-form
#' @param symtails logical value: if TRUE, the cutpoints are estimated with
#'   symmetric tails. If nb_of_cp=1, symtails is set to FALSE
#' @param dp numeric, number of decimal places the cutpoints are rounded to
#' @param plot_splines logical value: if TRUE, a penalized spline plot is
#'   created
#' @param all_splines logical value: if TRUE, the plot shows splines with
#'   different degrees of freedom. This may help identify if overfitting occurs.
#' @references Govindarajulu, U., & Tarpey, T. (2020). Optimal partitioning for
#'   the proportional hazards model. Journal of Applied Statistics, 49(4),
#'   968–987. https://doi.org/10.1080/02664763.2020.1846690
#' @returns returns the estimated cutpoints and the characteristics of the
#'   groups in relation to the original data set
#' @importFrom survival coxph
#' @importFrom survival Surv
#' @importFrom stats AIC complete.cases median quantile rnorm
#' @importFrom utils globalVariables
#' @export
#'
est_cutpoint <-
function(cpvarname,
         time         = "time",
         event        = "event",
         covariates   = NULL,
         data         = data,
         nb_of_cp     = 1,
         bandwith     = 0.1,
         ushape       = FALSE,
         symtails     = FALSE,
         dp           = 2,
         plot_splines = TRUE,
         all_splines  = TRUE) {

   #' Check if the input is correct
   #' -------------------------------------------------------------------------
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
      cat("\nPlease note: For at least one observation the follow-up time is","\u2264","0\n")
      cat("this can lead to an error message of the Cox regression\n")
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

   # "cpvarname" is used for labelling and "biomarker" for calculations
   biomarker <- cpvarname

   #' Data frame contains only the variables that must be included in the model
   cpdata <- data[, c(biomarker, time, event, covariates)]

   #' nrm = Number of rows in cpdata before possible changes
   nrm_start  <- nrow(cpdata)

   colnames(cpdata)[1:3] <- c("biomarker", "time", "event")

   #' Observations are deleted for missing values of the cutpoint variable
   cpdata <- cpdata[complete.cases(cpdata$biomarker), ]

   #' If there are more than 1000 observations, a random sample with 1000
   #'     observations is created
   sample_yes <- FALSE

   if (nrow(cpdata) > 1000) {
      cpdata <- cpdata[sample(nrow(cpdata), size=1000), ]
      sample_yes <- TRUE
   }

   #' Data frame is sorted by biomarker in ascending order
   cpdata <- cpdata[order(cpdata$biomarker), ]
   biomarker <- cpdata$biomarker

   cov_ <- cbind(cpdata)

   #' Remove elements from cov_ (time, event, biomarker)
   cov_ <- cov_[!names(cov_) %in% c("time", "event", "biomarker")]

   cov_ <- as.matrix(cov_)

   #' nrm = Number of rows in cpdata after removing observations with
   #'     missing values in biomarker
   nrm  <- nrow(cpdata)

   #if (is.null(covariates) && nrm < 200 && ushape == TRUE) {
   #   cat(" \n")
   #   cat("! With few observations and without the use of covariates, a warning \n")
   #   cat("may occur: - Loglik converges before ... may be infinite. - \n")
   #   cat("The test that is triggered to generate this warning is very sensitive.")
   #   cat(" \n")
   #}

   #' If there are no covariates, cov is defined as a constant
   if (is.null(covariates)) {cov_ <- rep(1, nrm)}

   #' If only one cutpoint is estimated, symtails and ushape is set to FALSE
   if (nb_of_cp == 1) {
      symtails <- ushape <- FALSE
   }

   # Prevent infinite of Cox regression in case of no covariate and ushape==TRUE
   if (is.null(covariates) && ushape == TRUE) {

      # Calculate SD on basis of the median of the biomarker:

      median_biomarker <- median(biomarker, na.rm = TRUE)

      # Prevent division of 0
      ifelse (median_biomarker == 0, median_biomarker <- 0.1, median_biomarker)
      if (median_biomarker > 0) {sd_value <- median_biomarker / 1e+8} else {
         sd_value <- abs( 1/(median_biomarker * (1e+8)))}

      # Set cov_ as infinite preventer
      cov_ <- rnorm(nrm, sd = sd_value)
   }

   #' Numbers of observations which should at least remain in line
   m.perm <- combine_factors(bandwith, nb_of_cp, nrm, symtails)

   #' If ushape is TRUE then create new m.perm with 2 categories, as u-shape only
   #'     has 2 categories
   if (ushape == TRUE) {
      m.perm[m.perm == 3] <- 1
   }

   #' If ushape is TRUE and bandwith < 0.1, then bandwith is set to 0.1
   change_bw <- FALSE
   if (ushape == TRUE && bandwith < 0.1) {
      bandwith <- 0.1
      change_bw <- TRUE
   }

   loop_nr <- 0

   #' Number of variants of the biomarker, the number of times the loop must
   #'     be run through
   nbr.m.perm <- nrow(m.perm)

   #' Vector for AIC values
   AIC_values <- rep(NA, nbr.m.perm)

   #' After 5%, remaining time is communicated for this timefactor is defined
   timefactor <- 0.05
   nbr.m.perm.time <- round(nbr.m.perm * timefactor, 0)

   #' Start time measurement
   ptm <- proc.time()


   for (i in 1:(nbr.m.perm)) {
      loop_nr <- loop_nr + 1

      biomarker_dicho <- as.factor(m.perm[i, ])

      result.cox <-
         survival::coxph(Surv(time, event) ~ cov_ + biomarker_dicho, data =
                            cpdata)

      AIC_values[i] <- AIC(result.cox)

      rm(result.cox)

      #' Show user approx. remaining time
      if (nbr.m.perm.time == loop_nr) {
         tm <- proc.time() - ptm
         cat("\nApprox. remaining time for estimation in seconds:", round((tm[3] * (1 / timefactor)
         ), 0), "\n")
         ptm <- proc.time()
      }

   } #' End: for (i in 1:nbr.m.perm)

   rm(i)
   rm(timefactor)
   rm(tm)
   rm(nbr.m.perm.time)


   #' if ushape==TRUE then create new m.perm with 3 categories, as ushape
   #'    has only 2 categories
   if (ushape == TRUE) {
      m.perm <- combine_factors(bandwith, nb_of_cp, nrm, symtails)
   }


   #' Cutpoint is at position minAIC_row_nb of , all those less than or equal to
   #'    Cutpoint (minAIC_row_nb) belong to the first group
   minAIC_row_nb <- which.min(AIC_values)
   cp1_position <- sum(m.perm[minAIC_row_nb, ] == 1)
   cp2_position <- sum(m.perm[minAIC_row_nb, ] <= 2)

   #' Generate vector with cutpoints
   cp <- c(NA, NA)


   cp[1] <- biomarker[cp1_position]
   names(cp[1]) <- "CP1"

   cp[2] <- biomarker[cp2_position]
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


   #' Output:------------------------------------------------------------------

   #' Get the counts of each group
   group_counts <- table(m.perm[minAIC_row_nb, ])

   #' Calculate the total count
   total_count <- sum(group_counts)

   #' Get the percentage
   percbygroup <- group_counts / total_count

   #' Calculate numbers per group for original dataset
   nbbygroup <- round(percbygroup * nrm_start, 0)

   cat("--------------------------------------------------------------------\n")
   if (sample_yes == TRUE){
   cat("! Because the number of obervations in original dataset is >1000\n")
   cat("  a random sample of 1000 obervations is used for estimating the cutpoints\n")
   }
   cat("--------------------------------------------------------------------\n")
   cat("SETTINGS:\n")
   cat(" Cutpoint-variable                = ", cpvarname, "\n")
   cat(" Number of cutpoints   (nb_of_cp) = ", nb_of_cp, "\n")
   if (change_bw == TRUE) {
      cat(" Min. group size in %  (bandwith) = ", bandwith,
          " (was set to 0.1 because ushape is TRUE)\n")
      } else {
   cat(" Min. group size in %  (bandwith) = ", bandwith, "\n")}
   cat(" Symmetric tails       (symtails) = ", symtails,
       "  (is set FALSE if nb_of_cp = 1)\n")
   cat(" Cutpoints for u-shape (ushape)   = ", ushape,
       "  (is set FALSE if nb_of_cp = 1)\n")
   cat("--------------------------------------------------------------------\n")
   if (is.null(covariates)) {
      cat("No covariates were selected\n")
   } else {
      cat("Covariates or factors are:\n")
      cat(" ",covariates, "\n") }
   cat("--------------------------------------------------------------------\n")
   cat("Minimum group size is ", round((nrm_start * bandwith), 0),
     " (", bandwith * 100, "% of sample size in original dataset, N = ", nrm_start,  ")\n", sep = "")
   cat("--------------------------------------------------------------------\n")
   cat("Number of Cutpoints searching for:", nb_of_cp, "\n")

   if (nb_of_cp == 1) {

      cat("Cutpoint:", cpvarname, "\u2264", cp[1], "\n")
      cat("-----------------------------------------------------------------\n")
      cat("Group size in relation to the original data set\n")
      cat(" Total:   N = ", nrm_start, "\n", sep = "")
      cat(" Group 1: n = ", nbbygroup[1], " (", round(percbygroup[1]*100,1), "%)\n", sep = "")
      cat(" Group 2: n = ", nbbygroup[2], " (", round(percbygroup[2]*100,1), "%)\n", sep = "")

      cp <- cp[-2]
   }

   if (nb_of_cp == 2) {

      cat(" 1.Cutpoint:", cpvarname, "\u2264", cp[1], "\n")
      cat(" 2.Cutpoint:", cpvarname, "\u2264", cp[2], "\n")
      cat("-----------------------------------------------------------------\n")
      cat("Group size in relation to the original data set\n")
      cat(" Total:   N = ", nrm_start, "\n", sep = "")
      cat(" Group 1: n = ", nbbygroup[1], " (", round(percbygroup[1]*100,1), "%)\n", sep = "")
      cat(" Group 2: n = ", nbbygroup[2], " (", round(percbygroup[2]*100,1), "%)\n", sep = "")
      cat(" Group 3: n = ", nbbygroup[3], " (", round(percbygroup[3]*100,1), "%)\n", sep = "")

   }

   # End: Output --------------------------------------------------------------


   returnlist <- list(
      cp = cp,
      cpdata = cpdata,
      cpvarname = cpvarname,
      nb_of_cp = nb_of_cp,
      dp = dp,
      AIC_values = AIC_values
   )

   #' Plot Splines
   if (plot_splines == TRUE) {
      splines_plot(returnlist, show_splines = all_splines)
   }

   return(returnlist)

} # End: est_cutpoint <- function
