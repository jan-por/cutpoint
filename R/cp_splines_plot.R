#' @title Plot penalized smoothing splines and shows cutpoints from `cpobj` object
#'
#' @description Create penalized smoothing splines plot with different degrees
#'    of freedom and shows the cutpoints of the dichotomised variable.
#' @name cp_splines_plot
#' @param cpobj list, contains variables for pspline plot:
#' * `nb_of_cp` (number of cutpoints)
#' * `cp` (contain one or two cutpoint/s)
#' * `dp` (digits for plot)
#' * `cpvarname` (name of the variable for that the cutpoints are estimated)
#' * `cpdata` a data frame, contains the following variables: a variable that is
#'   dichotomized, `time` (follow-up time), `even`t (status indicator),
#'   `covariates` (a vector with the names of the covariates and/or factors))
#' @param show_splines logical, if `TRUE`, The plot shows splines with
#'   different degrees of freedom. This may help determine whether
#'   misspecification or overfitting occurs.
#' @param adj_splines logical, if `TRUE`, the splines are adjusted for the
#'   covariates. Default is `TRUE`.
#' @returns Plots penalized smoothing splines and shows the cutpoints.
#' @examples
#' cpvar <- rnorm(100, mean = 100, sd = 10)
#' time <- seq(1, 100, 1)
#' event <- rbinom(100, 1, 0.5)
#' datf <- data.frame(time, event, cpvar)
#' plot_splines_list <- list(cpdata = datf, nb_of_cp = 1, cp = 95, dp = 2,
#'     cpvarname = "Biomarker")
#' cp_splines_plot(plot_splines_list)
#' @importFrom stats quantile
#' @importFrom survival coxph Surv pspline
#' @importFrom graphics abline legend lines
#' @importFrom utils globalVariables
#' @export
#'
#' @seealso [cp_est()] for main function of the package, [cp_value_plot()]
#'   for Value plots and Index plots
NULL
cp_splines_plot <- function(cpobj, show_splines = TRUE, adj_splines = TRUE) {

      #' Check if cpobj is a list
      if (!is.list(cpobj)) {
         stop("Cutpoint object (cpobj) must be a list")
      }

      #' Extract necessary variables from cpobj
      nb_of_cp   <- cpobj$nb_of_cp
      cp         <- cpobj$cp
      dp         <- cpobj$dp
      cpdata     <- cpobj$cpdata
      cpvarname  <- cpobj$cpvarname
      covariates <- cpobj$covariates

      #' Check variables

      if (!is.numeric(nb_of_cp))
         stop("nb_of_cp must be numeric")
      if (nb_of_cp != 1 &
          nb_of_cp != 2)
         stop("nb_of_cp must be 1 or 2")

      if (!is.numeric(cp)) stop("cp must be numeric")

      if (!is.numeric(dp)) stop("dp must be numeric")
      if (dp %% 1 != 0)
         stop("dp must be an integer")
      if (dp < 0)
         stop("dp must be 0 or greater than 0")
      if (dp > 19)
         stop("dp must be smaller than 20")

      if (!is.logical(show_splines))
         stop("show_splines must be logical (TRUE or FALSE)")

      if (!is.logical(adj_splines))
         stop("show_splines must be logical (TRUE or FALSE)")

      if (!("time" %in% names(cpdata))) {
         stop("time must be a column in cpdata")
      }
      if (length(unique(cpdata$time)) == 1) {
         stop("time must have more than one unique value")
      }
      if (length(unique(cpdata$time)) < 3) {
         stop("time must have more than two unique values")
      }

      if (!("event" %in% names(cpdata))) {
         stop("event must be a column in cpdata")
      }
      if (length(unique(cpdata$event)) != 2) {
         stop("event must have two unique values")
      }
      if (length(unique(cpdata$event)) < 2) {
         stop("event must have more than one unique value")
      }

      if (!("cpvar" %in% names(cpdata))) {
         stop("cpvar must be a column in cpdata")
      }
      if (length(unique(cpdata$cpvar)) == 1) {
         stop("cpvar must have more than one unique value")
      }
      if (length(unique(cpdata$cpvar)) < 3) {
         stop("cpvar must have more than two unique values")
      }

      if (!is.character(cpvarname)) {
         stop("cpvarname must be a character")
      }


      if (!is.null(covariates) & !is.character(covariates)) {
         stop("covariates must be a character vector or NULL")
      }
      if (!all(covariates %in% colnames(cpdata))) {
         stop("all covariates must be included in cpdata, check the names of the
              covariates")
      }

      cpvar <- cpdata$cpvar

      #' Get quantiles of cpvar
      q <- quantile(cpvar, na.rm = TRUE)

      #' Get optimal degree of freedom
      # tfit <- survival::coxph(
      #    formula = Surv(time, event) ~ survival::pspline(
      #       x = cpvar,
      #       df = 0,
      #       caic = TRUE,
      #       plot = FALSE
      #    ) + age + sex,
      #    data = cpdata
      # )


      #' Define formula (FML) for Cox-reg. with cpvar and vector of covariates:

      #' Without covariates splines cannot be adjusted
      if(is.null(covariates)) { adj_splines <- FALSE }

      if(adj_splines == TRUE) {
         FML <- as.formula(paste0(' ~ survival::pspline(x = cpvar, df = 0,
                                  caic = TRUE, plot = FALSE) +',
                                  paste(covariates, collapse = "+")))
      } else {
         #' If adj_splines == FALSE
         #' Define formula for Cox-regression with cpvar only
         FML <- as.formula(paste0(' ~ survival::pspline(x = cpvar, df = 0,
                                  caic = TRUE, plot = FALSE)'))
      }

      try(tfit <- survival::coxph(update.formula(Surv(time, event)~., FML),
                                  data = cpdata))

      degfr_optimal <- round(tfit$df, 1)

      #' Define degree of freedom used for termplot
      degfr <- c(5, 4, 3, 2)

      tempcolors <- c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6","black")

      #' Define main text for plot
      if(adj_splines == TRUE) {
         main_text_part1 <- "Covariate-adjusted splines plot with "
      }
      else {
         main_text_part1 <- "Splines plot with "
      }

      if (show_splines ==  TRUE) {
         main_text <- paste0(main_text_part1,
                             "different degrees of freedom (df)"
                            )
      }
      else {
         main_text <- paste0(main_text_part1,
                             "optimal degrees of freedom (df = ",
                             degfr_optimal, ")"
                            )
      }

      #' Visualization: pspline Plot - termplot -----------------------
      termplot(
         tfit,
         terms = 1,
         se       = TRUE,
         col.term = "black",
         col.se   = "black",
         lwd.term = 2,
         font.lab = 2,
         xlabs    = paste("Cutpoint variable: ", cpvarname),
         ylabs    = "log relative hazard",
         main     = main_text,
         sub      = paste0("Mean: ",
            round(mean(cpvar), 2),
            " (dashed, red line),  Q1: ",
            round(q[2], 2),
            " (grey),  Median: ",
            round(q[3], 2),
            " (grey),  Q3: ",
            round(q[4], 2) ,
            " (grey)"
         )
      ) # End: termplot


      #' Show splines with different degrees of freedom and add legend
      if (show_splines ==  TRUE) {

         for (i in 1:length(degfr)) {
            #' Formula (FML) for Cox-reg. with cpvar and vector of covariates:

            if(adj_splines == TRUE) {

               #' Formula (FML) for Cox-reg. with cpvar and vector of covariates
               FML <- as.formula(paste0(' ~ survival::pspline(x = cpvar, df = ',
                                        degfr[i],', caic = TRUE) +',
                                        paste(covariates, collapse = "+")
                                        )
                                 )
            } else {
               #' If adj_splines == FALSE
               #' Define formula for Cox-regression with cpvar only
               FML <- as.formula(paste0(' ~ survival::pspline(x = cpvar, df = ',
                                        degfr[i],', caic = TRUE)'
                                        )
                                 )
            }

            try(tfit <- survival::coxph(update.formula(Surv(time, event)~., FML),
                                        data = cpdata
                                        )
                )


            temp <- termplot(tfit, se = FALSE, terms = 1, plot = FALSE)
            lines(temp$cpvar$x,
                  temp$cpvar$y,
                  col = tempcolors[i],
                  lwd = 2)
            rm(temp)
         }

         legend(
            "bottomright",
            cex = 0.8,
            paste0("df=", c(degfr[1:length(degfr)],
               paste((degfr_optimal[[length(degfr_optimal)]]), "(optimal)"
            ))),
            lty = 1,
            col = tempcolors[1:(length(degfr)+1)],
            lwd = 2
         )

      } # End: if (show_splines ==  TRUE)

      #' Add lines for mean and quantiles of cpvar
      abline(
         v = c(q[2], q[3], q[4], mean(cpvar)),
         col = c("grey", "grey", "grey", "red"),
         lty = c("dotted", "dotted", "dotted", "dashed")
      )

      if (nb_of_cp == 1) {
         abline(v = cp[1], col = "red", lwd = 2)
      }

      if (nb_of_cp == 2) {
         abline(v = cp[1], col = "red", lwd = 2)
         abline(v = cp[2], col = "red", lwd = 2)
      }

      #' Show cutpoints as legend in plot

      #' Define title for legend
      if (nb_of_cp == 1) { cptext <- "Cutpoint:  " } else {
         cptext <- "Cutpoints:  " }

      legend(
            "bottomleft",
            title = cptext,
            cex = 0.8,
            paste("\u2264", round((cp[1]), dp)),
            lty = 1,
            col = tempcolors[1:(length(degfr)+1)],
            lwd = 2
         )

      legend(
         "bottomleft",
         title = cptext,
         cex = 0.8,
         if (nb_of_cp == 1) {paste("\u2264", round((cp[1]), dp))} else {
            paste("\u2264", round((cp[1]),dp), "and", "\u2264",
                  round((cp[2]), dp)) },
         lty = 1,
         col = tempcolors[1:(length(degfr)+1)],
         lwd = 2
      )


      rm(FML)

      return(invisible())

   } # End
