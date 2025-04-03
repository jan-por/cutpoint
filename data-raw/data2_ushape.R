## code to prepare `cp_ushape_data` dataset goes here


#Generate cp_ushape_data as example
n <- 200
x <- rnorm(n)
y <- 0.5 * x^2
z <- 0.3 * x^2
covariate_1 <- z
log_relative_hazard <- y + z
#binary_variable <- rbinom(n, 1, 0.5)
survival_time <- rweibull(n, shape = 1, scale = exp(-log_relative_hazard))
censoring_time <- rexp(n, rate = 0.05)
censoring_indicator <- ifelse(survival_time < censoring_time, 1, 0)
time <- pmin(survival_time, censoring_time)
event <- ifelse(survival_time < censoring_time, 1, 0)
cp_ushape_data <- data.frame(biomarker = x, covariate_1, time, event)

# Check how much events
unique_event <- 200-(sum(event))
unique_event
