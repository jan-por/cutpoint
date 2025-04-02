#' @title Plot psplines
#'
#' @description Create pspline plot with different degrees of freedom and
#'   shows the cutpoints of the biomarker
#' @name pspplot
#' @param cpobj list contains variables for pspline plot
#' @param show.splines logical, if TRUE shows the splines with different DF
#' @returns returns a plot with the pspline and shows the cutpoints
#' @examples
#' biomarker <- rnorm(100, mean = 100, sd = 10)
#' time <- seq(1, 100, 1)
#' event <- rbinom(100, 1, 0.5)
#' datf <- data.frame(time, event, biomarker)
#' plot_splines_list <- list(cpdata = datf, nbofcp = 1, cp = 95, dp = 2,
#'     cpvarname = "Biomarker")
#' pspplot(plot_splines_list)
#' @importFrom stats quantile
#' @importFrom survival coxph Surv
#' @importFrom graphics abline legend lines mtext
#' @importFrom utils globalVariables
#' @export
#'
#' @keywords cutpoint pspline plot visualization termplot
#'
pspplot <-
   function(cpobj, show.splines = TRUE) {

      if (!is.list(cpobj)) {
         stop("Cutpoint object (cpobj) must be a list")
      }

      nbofcp    <- cpobj$nbofcp
      cp        <- cpobj$cp
      dp        <- cpobj$dp
      cpdata    <- cpobj$cpdata
      cpvarname <- cpobj$cpvarname

      if (!is.numeric(nbofcp))
         stop("nbofcp must be numeric")
      if (nbofcp != 1 &
          nbofcp != 2)
         stop("nbofcp must be 1 or 2")

      if (!is.numeric(cp)) stop("cp must be numeric")

      if (!is.numeric(dp)) stop("dp must be numeric")
      if (dp %% 1 != 0)
         stop("dp must be an integer")
      if (dp < 0)
         stop("dp must be 0 or greater than 0")
      if (dp > 19)
         stop("dp must be smaller than 20")

      if (!is.logical(show.splines))
         stop("show.splines must be logical (TRUE or FALSE)")

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

      if (!("biomarker" %in% names(cpdata))) {
         stop("biomarker must be a column in cpdata")
      }
      if (length(unique(cpdata$biomarker)) == 1) {
         stop("biomarker must have more than one unique value")
      }
      if (length(unique(cpdata$biomarker)) < 3) {
         stop("biomarker must have more than two unique values")
      }

      if (!is.character(cpvarname)) {
         stop("cpvarname must be a character")
      }

      biomarker <- cpdata$biomarker

      #' Get quantiles of biomarker
      q <- quantile(biomarker, na.rm = TRUE)

      #' Get optimal degree of freedom
      tfit <- survival::coxph(
         formula = Surv(time, event) ~ pspline(
            x = biomarker,
            df = 0,
            caic = TRUE,
            plot = FALSE
         ) ,
         data = cpdata
      )
      degfr_optimal <- round(tfit$df, 1)

      #' Define degree of freedom used for termplot
      degfr <- c(5, 4, 3, 2)

      tempcolors <- c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6","black")


      if (show.splines ==  TRUE) {
         main_text <- "Splines with different degrees of freedom (df)"
      }
      else {
         main_text <-
            paste0("Splines with optimal degrees of freedom (df = ",
                   degfr_optimal,
                   ")")
      }


      termplot(
         tfit,
         se       = TRUE,
         col.term = "black",
         col.se   = "black",
         lwd.term = 2,
         font.lab = 2,
         xlabs     = paste("Cutpoint variable: ", cpvarname),
         ylabs     = "log relative hazard",
         main     = main_text,
         sub      = paste0("Mean: ",
            round(mean(biomarker), 2),
            " (dashed, red line),  Q1: ",
            round(q[2], 2),
            " (grey),  Median: ",
            round(q[3], 2),
            " (grey),  Q3: ",
            round(q[4], 2) ,
            " (grey)"
         )
      ) # End: termplot



      if (show.splines ==  TRUE) {
         for (i in 1:length(degfr)) {
            try(tfit <- survival::coxph(
               formula = Surv(time, event) ~ pspline(
                  x = biomarker,
                  df = degfr[i],
                  caic = TRUE
               ),
               data = cpdata
            ))


            temp <- termplot(tfit, se = FALSE, plot = FALSE)
            lines(temp$biomarker$x,
                  temp$biomarker$y,
                  col = tempcolors[i],
                  lwd = 2)
            rm(temp)
         }

         legend(
            "topright",
            cex = 0.8,
            paste0("df=", c(degfr[1:length(degfr)],
               paste((degfr_optimal[[length(degfr_optimal)]]), "(optimal)"
            ))),
            lty = 1,
            col = tempcolors[1:(length(degfr)+1)],
            lwd = 2
         )

      } # End: if (show.splines ==  TRUE)

      abline(
         v = c(q[2], q[3], q[4], mean(biomarker)),
         col = c("grey", "grey", "grey", "red"),
         lty = c("dotted", "dotted", "dotted", "dashed")
      )

      if (nbofcp == 1) {
         abline(v = cp[1], col = "red", lwd = 2)
      }

      if (nbofcp == 2) {
         abline(v = cp[1], col = "red", lwd = 2)
         abline(v = cp[2], col = "red", lwd = 2)
      }

      legend(
         "bottomleft",
         title = "Cutpoint:  ",
         cex = 0.8,
         if (nbofcp == 1) {paste("\u2264", cp[1])} else {
            paste("\u2264", cp[1], "and", "\u2264", cp[2]) } ,
         lty = 1,
         col = tempcolors[1:(length(degfr)+1)],
         lwd = 2
      )

      return(invisible())

   } # End: Visualization: pspline Plot - termplot -----------------------
