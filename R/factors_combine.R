#' @title Combine Factors
#'
#' @description Intern function, used for creation of a matrix with all factor
#'   combinations of the cutpoint-variable
#' @param bandwith numeric, determines the minimum size per group of the
#'   dichitomised variable
#' @param nb_of_cp numeric, number of cutpoints searching for
#' @param nrm numeric, number of rows in cpdata after removing observations
#'   with missing values in biomarker
#' @param symtails logical, if TRUE the tails of the dichotomised variable
#'   are symmetrical
#' @returns Returns all factor combinations of the dichotomized variable.

factors_combine <-
   function(bandwith = 0.1, nb_of_cp = 1, nrm, symtails = FALSE) {

   if (!is.numeric(bandwith))
      stop("bandwith must be numeric")
   if (bandwith < 0 |
       bandwith > 0.3)
      stop("bandwith must be between 0 and 0.3")

   if (!is.numeric(nb_of_cp))
      stop("nb_of_cp must be numeric")
   if (nb_of_cp != 1 &
       nb_of_cp != 2)
      stop("nb_of_cp must be 1 or 2")

   if (!is.numeric(nrm))
      stop("nrm must be numeric")
   if (nrm %% 1 != 0)
      stop("nrm must be an integer")
   if (nrm < 1)
      stop("nrm must be greater than 0")
   if (nrm > 1000)
      stop("nrm must be smaller than 1000")

   if (!isTRUE(symtails) &
       !isFALSE(symtails))
      stop("symtails must be logical TRUE or FALSE")


   #' Creation of a matrix with all factor combinations of the dichotomised
   #'   variable
   m.perm.fun <- as.matrix(RcppAlgos::comboGeneral((nb_of_cp + 1),
                                                   nrm,
                                                   repetition = TRUE))

   #' Delete rows with less than ‘bandwith’ percent of the number of columns
   #' # convert percentages into number
   bandwith_nb <- round(ncol(m.perm.fun) * bandwith, 0)
   if (bandwith_nb < 1)
      bandwith_nb <- 1
   # as this variable is used for deletion, bandwith_nb must be reduced by 1
   bandwith_nb <- bandwith_nb - 1

   #' If nb_of_cp == 1
   #' How often does the 1 occur per line
   y <- rowSums(m.perm.fun == 1)
   m.perm.fun <- cbind(m.perm.fun, y)
   rm(y)

   #' Delete rows with less than ‘bandwith’ percent of the number of columns
   m.perm.fun <- m.perm.fun[!(m.perm.fun[, (ncol(m.perm.fun))] %in% 0:bandwith_nb), ]
   m.perm.fun <- m.perm.fun[, -seq(ncol(m.perm.fun), ncol(m.perm.fun))]

   #' If nb_of_cp == 2
   y <- rowSums(m.perm.fun == 2)
   m.perm.fun <- cbind(m.perm.fun, y)
   rm(y)

   #' Delete all rows with less than bandwith per cent of the nb. of columns
   m.perm.fun <- m.perm.fun[!(m.perm.fun[, (ncol(m.perm.fun))] %in% 0:bandwith_nb), ]
   m.perm.fun <- m.perm.fun[, -seq(ncol(m.perm.fun), ncol(m.perm.fun))]


   #' How often does the 2 occur per line
   if (nb_of_cp == 2) {
      y <- rowSums(m.perm.fun == 3)
      m.perm.fun <- cbind(m.perm.fun, y)
      rm(y)

      m.perm.fun <- m.perm.fun[!(m.perm.fun[, (ncol(m.perm.fun))] %in% 0:bandwith_nb), ]
      m.perm.fun <- m.perm.fun[, -seq(ncol(m.perm.fun), ncol(m.perm.fun))]
   }


   #' If m.perm.fun should be symmetrical
   if (symtails == TRUE) {
      y1 <- rowSums(m.perm.fun == 1)
      y3 <- rowSums(m.perm.fun == 3)
      y <- y1 - y3
      m.perm.fun <- cbind(m.perm.fun, y)
      rm(y1)
      rm(y3)
      rm(y)
      # The length of the tails can differ by a maximum of 1.
      m.perm.fun <- m.perm.fun[(m.perm.fun[, (ncol(m.perm.fun))] %in% 0:1), ]
      m.perm.fun <- m.perm.fun[, -seq(ncol(m.perm.fun), ncol(m.perm.fun))]
   }


   return(m.perm.fun)


   } #' End: factors_combine <- function( --------------------------------------------------------------
