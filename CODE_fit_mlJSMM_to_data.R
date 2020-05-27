##################################################################################
##              Multiple landscapes Joint species movement modelling            ##
##                                   Fit Model                                  ##
##################################################################################

## code from the manuscript "Forest and connectivity loss drive changes in movement behavior of bird species"
## authors: "Danielle Leal Ramos, Marco Aurélio Pizo, Milton Cezar Ribeiro, Rafael Souza Cruz, Juan Manuel Morales, Otso Ovaskainen"
## date: "May, 2020"

## Required packages

if (!require("MCMCpack")) install.packages("MCMCpack") # for function riwish
if (!require("mvtnorm")) install.packages("mvtnorm") # for functions dmvnorm and rmvnorm
library(Matrix)
library(doParallel) # for parallel computing

## Required functions

trunca <- function(x) min(max(x, 10^-5), 10^5) # required function to run the JSMM function below
source("MLjsmm.r") # this is the JSMM adapted for multiple-landscapes data analysis

## Load dataset

load("multi_landscape_real_bird_dataset.Rdata")

## Fitting the model to spatial aspects of the movement ##

source("loglik_steps.r") # log-likelihood function for spatial movement parameters

# cl <- makeCluster(5)
# registerDoParallel(cl)

# tpm <- foreach(repl = 1:5, .export = c("loglik", "trunca", "jsmm.mcmc"), .packages = c('Matrix','mvtnorm','MCMCpack')) %dopar%
# {

output <- jsmm.mcmc(loglikelihood = loglik, data = data, n.iter = 100, n.adapt.iter = 20, 
                    n.thin = 2, rotate = TRUE)
# }
# stopCluster(cl)

save(output, data, file = "posteriors_steps_100it.Rdata")

## Fitting the model to temporal aspects of the movement ##

source("loglik_time.r") # log-likelihood function for temporal movement parameters

# cl <- makeCluster(5)
# registerDoParallel(cl)

# tpm <- foreach(repl = 1:5, .export = c("loglik", "trunca", "jsmm.mcmc"), .packages = c('Matrix','mvtnorm','MCMCpack')) %dopar%
# {

output <- jsmm.mcmc(loglikelihood = loglik, data = data, n.iter = 100, n.adapt.iter = 20, 
                    n.thin = 2, rotate = TRUE)
# }
# stopCluster(cl)

save(output, data, file = "posteriors_time_100it.Rdata")
