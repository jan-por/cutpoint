#' @title Estimates cutpoints
#'
#' @description One or two cutpoints of a variable are estimated by the AIC
#'   criterion in a Cox proportional hazards model. The cutpoints are estimated
#'   by dichotomizing the variable and testing the log-likelihood ratio test.
#'   The cutpoints with the lowest AIC value are chosen.
#' @name est.cutpoint
#'
#' @param cpvarname character, name of the variable for which the cutpoints are
#'     estimated
#' @param time character, name of the time variable
#' @param event character, name of the event variable
#' @param cofactors character, vector with the names of the covariates or/ and
#'   factors
#' @param data data.frame, data set with the variables
#' @param nbofcp numeric, number of cutpoints searching for
#' @param bandwith numeric, minimum group size in percent of the total sample
#'   size, bandwith must be between 0 and 0.3
#' @param ushape logical, if TRUE the cutpoints are estimated with three
#'   categories
#' @param symtails logical, if TRUE the cutpoints are estimated with symmetric
#'   tails
#' @param dp numeric, number of decimal places the cutpoints are rounded to
#' @param plot_splines logical, if TRUE creates pspline plot
#' @param all_splines logical, if TRUE all splines are shown
#'
#' @returns returns the estimated cutpoints and the characteristics of the groups
#' @export
#'
#' @import utils
utils::globalVariables(c("stats", "complete.cases", "AIC"))

est.cutpoint <-
function(cpvarname,
         time = "time",
         event = "event",
         cofactors,
         data = data,
         nbofcp = 1,
         bandwith = 0.1,
         ushape     = FALSE,
         symtails   = FALSE,
         dp         = 2,
         plot_splines = TRUE,
         all_splines = TRUE) {

   if (!is.character(cpvarname)) {
      stop("time must be a character")
   }

   #if (length(unique(data$cpvarname)) == 1) {
   #   stop("cpvarname must have more than one unique value")
   #}
   #if (length(unique(data$cpvarname)) < 3) {
   #   stop("cpvarname must have more than two unique values")
   #}

   if (!is.character(time)) {
      stop("time must be a character")
   }
   #if (length(unique(data$time)) == 1) {
   #   stop("time must have more than one unique value")
   #}
   #if (length(unique(data$time)) < 3) {
   #   stop("time must have more than two unique values")
   #}

   if (!is.character(event)) {
      stop("event must be a character")
   }
   # if (length(unique(data$event)) != 2) {
   #    stop("event must have two unique values")
   # }
   # if (length(unique(data$event)) < 2) {
   #    stop("event must have more than one unique value")
   # }

   if (!is.character(cofactors)) {
      stop("cofactors must be a character")
   }

   if (!is.data.frame(data)) {
      stop("data must be a data.frame")
   }

   if (!is.numeric(nbofcp))
      stop("nbofcp must be numeric")
   if (nbofcp != 1 &
       nbofcp != 2)
      stop("nbofcp must be 1 or 2")

   if (!is.numeric(bandwith))
      stop("bandwith must be numeric")
   if (bandwith < 0 |
       bandwith > 0.3)
      stop("bandwith must be between 0 and 0.3")

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
   cpdata <- data[, c(biomarker, time, event, cofactors)]

   #' nrm = Number of rows in cpdata before possible changes
   nrm_start  <- nrow(cpdata)

   colnames(cpdata)[1:3] <- c("biomarker", "time", "event")

   #' Observations are deleted for missing values of the cutpoint variable
   cpdata <- cpdata[complete.cases(cpdata$biomarker), ]

   #' If there are more than 1000 observations, a random sample with 1000 observations is created
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

   #' If only one cutpoint is estimated, symtails and ushape is set to FALSE
   if (nbofcp == 1) {
      symtails <- ushape <- FALSE
   }

   #' Numbers of observations which should at least remain in line
   m.perm <- combine.factors(bandwith, nbofcp, nrm, symtails)

   #' If ushape==TRUE then create new m.perm with 3 categories, as ushape only
   #'     has 2 categories
   if (ushape == TRUE) {
      m.perm[m.perm == 3] <- 1
   }

   loop_nr <- 0

   #' Number of variants of the biomarker, the number of times the loop must
   #'     be run through
   nbr.m.perm <- nrow(m.perm)

   #' vector for AIC values
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
         cat("Approx. remaining time in seconds:", round((tm[3] * (1 / timefactor)
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
      m.perm <- combine.factors(bandwith, nbofcp, nrm, symtails)
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
   cat(" Number of cutpoints   (nbofcp) = ", nbofcp, "\n")
   cat(" Min. group size in %  (bandwith) = ", bandwith, "\n")
   cat(" Symmetric tails       (symtails) = ", symtails,
       "   (is set FALSE if nbofcp = 1)\n")
   cat(" Cutpoints for u-shape (ushape)   = ", ushape,
       "   (is set FALSE if nbofcp = 1)\n")
   cat("--------------------------------------------------------------------\n")
   cat("Cofactors are:\n")
   cat(" ",cofactors, "\n")
   cat("--------------------------------------------------------------------\n")
   cat("Minimum group size for original dataset is ", round((nrm_start * bandwith), 0),
     " (", bandwith * 100, "% of sample size, N = ", nrm_start,  ")\n", sep = "")
   cat("--------------------------------------------------------------------\n")
   cat("Number of Cutpoints searching for:", nbofcp, "\n")

   if (nbofcp == 1) {

      cat("Cutpoint for", cpvarname, "=", cp[1], "\n")
      cat("-----------------------------------------------------------------\n")
      cat("Group sizes for original dataset\n")
      cat(" Total:   N = ", nrm_start, "\n", sep = "")
      cat(" Group 1: n = ", nbbygroup[1], " (", round(percbygroup[1]*100,1), "%)\n", sep = "")
      cat(" Group 2: n = ", nbbygroup[2], " (", round(percbygroup[2]*100,1), "%)\n", sep = "")

      cp <- cp[-2]
   }

   if (nbofcp == 2) {

      cat(" 1.Cutpoint for", cpvarname, "=", cp[1], "\n")
      cat(" 2.Cutpoint for", cpvarname, "=", cp[2], "\n")
      cat("-----------------------------------------------------------------\n")
      cat("Group sizes of the dataset\n")
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
      nbofcp = nbofcp,
      dp = dp,
      AIC_values = AIC_values
   )

   #' Plot Splines
   if (plot_splines == TRUE) {
      pspline.plot(returnlist, show.splines = all_splines)
   }

   return(returnlist)

} # End: est.cutpoint <- function
