# Assignment 4
# Pia Elbe

cat("\014")  # Clear console
rm(list=ls()) # Delete all variables
graphics.off() # Close all open plots

# getting JAGS ...I didn't end up using JAGS. I used Stan.

# install.packages("rjags")
# install.packages("runjags")

# library("rjags")
# library("runjags")

library("rstan")
library("shinystan")

############################################################
# Task A. 1.
# getting the data set up for figure 6.4.

setwd("~/Desktop/bayesian_assignments/assignment_4")
source("DBDA2E-utilities.R")
source("BernBeta.R")

# specify the prior (page 138)

t = 0.5                 # prior mode (could also use m for mean)
n = 500                 # effective prior sample size
a = t*(n-2) + 1         # beta shape parameter a
b = (1-t)*(n-2) + 1     # beta shape parameter b
a
b
Prior = c(a,b)          # prior as a vector with the two shape parameters.

# specify the data (page 139)

N = 20   #total number of flips
z = 17   #ntotal number of heads
Data = c(rep(0,N-z),rep(1,z))  # convert N and z into vectors of 0's and 1's.

# recreating the left side of figure 6.4 page 135.

openGraph(width=5, height=7)
posterior = BernBeta( priorBetaAB=Prior, Data=Data, plotType="Bars", showCentTend = "Mode", showHDI = TRUE, showpD = FALSE)
saveGraph(file="GraphforFig6_4", type="jpg")

#Using the Data in figure 6.4 left side, making a list.
N <- 20
z <- 17
y <- c(rep(0,N-z),rep(1,z))
stan_1_data <- list(N=N, y=y)
stan_1_data

hist(y)

# compile model in stan (see the stan file)
model <- stan_model("figure1_left_model.stan")

#pass data to stan and run model
parallel::detectCores() # make sure you have the cores on your machine to run the chains in parallel.
options(mc.cores=4) # run 4 parallel cores for 4 chains, speeds the process
fit <- sampling(model, list(N=N, y=y), iter=2000, chains=4)

#diagnose
print(fit)
# by default, stan uses half of the iterations for warmup. this adds up to 400 posterior samples.
# the sampler is NUTS (no u turn sampling), apparently a type of hamiltonian.
#Rhat close to one means it did a good job of finding all of the distribution.

#We can look at a density plot:
stan_dens(fit) + xlim(0,1)

#graph

parameters <- extract(fit) #in the global environment, check parameters has 3 entries: a, b, and lp__ (low probability)
hist(parameters$a) #I get an error that 'x' must be numeric. 

#with shinystan
launch_shinystan(fit) #this will load a new window in the default browser.
#in shinystan you can look up the diagnostics and plots of estimates.

###################################################################
# Task A. 2.
#a. y = [1,0,1,1,0,1,1,1,0,1,1,1,1,1]

#let's define the data again.
N <- 14
y <- c(1,0,1,1,0,1,1,1,0,1,1,1,1,1)
stan_1_data <- list(N=N, y=y)
stan_1_data

hist(y)

# I'll just use the same model again.
model <- stan_model("figure1_left_model.stan")

#pass data to stan and run model
options(mc.cores=4) # run 4 parallel cores for 4 chains, speeds the process
fit2 <- sampling(model, list(N=N, y=y), iter=2000, chains=4)

#diagnose
print(fit2)

#We can look at a density plot:
stan_dens(fit2) + xlim(0,1)

#graph
parameters2 <- extract(fit2)

launch_shinystan(fit2) #this will load a new window in the default browser.

#uncertainty, credibility intervals
library(bayestestR)
library(dplyr)
library(ggplot2)

describe_posterior(
  fit2,
  centrality = "mean",
  dispersion = FALSE,
  ci = 0.95,
  ci_method = "hdi",
  test = c("p_direction", "rope"),
  rope_range = c(0,1),
  rope_ci = 0.95,
  keep_iterations = FALSE
)

# Compute HDI and ETI (equal tailed)
ci_hdi <- ci(fit2, method = "HDI")
ci_hdi
ci_eti <- ci(fit2, method = "ETI")
ci_eti

##########################################
#2.b
#Let's make another model for z = [1,0,0,0,0,0,0,1,1,0].

N <- 10
y <- c(1,0,0,0,0,0,0,1,1,0)
stan_2_data <- list(N=N, z=z)
stan_2_data

hist(y)

# Now the model again.
model <- stan_model("figure1_left_model.stan")

#pass data to stan and run model
options(mc.cores=4) # run 4 parallel cores for 4 chains, speeds the process
fit3 <- sampling(model, list(N=N, y=y), iter=2000, chains=4)

#diagnose
print(fit3)

#We can look at a density plot:
stan_dens(fit3) + xlim(0,1)

#to compare:
stan_dens(fit2) + xlim(0,1)

#graph
parameters3 <- extract(fit3)

launch_shinystan(fit3) #this will load a new window in the default browser.

ci_hdi <- ci(fit3, method = "HDI")
ci_hdi
ci_eti <- ci(fit3, method = "ETI")
ci_eti











#Most of the below code is wrong for the task. My first attempt. Does not matter anymore.
#Please disregard.

#################################################
# now with JAGS. taken from Jags-ExampleScript.R
#################################################
library("runjags")
fileNameRoot="Jags-Script" # For output file names.

# load the data:

y = Data$y
Ntotal = length(y)
dataList = list(
  y = y ,
  Ntotal = Ntotal
)

# Load the data:
myData = read.csv("z15N50.csv") # Read data file; must be in curr. work. dir.
y = myData$y        # The y values are in the column named y.
Ntotal = length(y)  # Compute the total number of flips.
dataList = list(    # Put the information into a list.
  y = y ,
  Ntotal = Ntotal 
)

# Define the model:
modelString = "
model {
  for ( i in 1:Ntotal ) {
    y[i] ~ dbern( theta )
  }
  theta ~ dbeta( 1 , 1 )
}
" 
writeLines( modelString , con="TEMPmodel.txt" )

# Initialize the chains based on MLE of data.
# Option: Use single initial value for all chains:
#  thetaInit = sum(y)/length(y)
#  initsList = list( theta=thetaInit )
# Option: Use function that generates random values for each chain:
initsList = function() {
  resampledY = sample( y , replace=TRUE )
  thetaInit = sum(resampledY)/length(resampledY)
  thetaInit = 0.001+0.998*thetaInit # keep away from 0,1
  return( list( theta=thetaInit ) )
}

# Run the chains:
jagsModel = jags.model( file="TEMPmodel.txt" , data=dataList , inits=initsList , 
                        n.chains=3 , n.adapt=500 )
update( jagsModel , n.iter=500 )
codaSamples = coda.samples( jagsModel , variable.names=c("theta") ,
                            n.iter=3334 )
save( codaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") )

# Examine the chains:

# Convergence diagnostics:
diagMCMC( codaObject=codaSamples , parName="theta" )
saveGraph( file=paste0(fileNameRoot,"ThetaDiag") , type="eps" )
# Posterior descriptives:
openGraph(height=3,width=4)
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( codaSamples[,"theta"] , main="theta" , xlab=bquote(theta) )
saveGraph( file=paste0(fileNameRoot,"ThetaPost") , type="eps" )
# Re-plot with different annotations:
plotPost( codaSamples[,"theta"] , main="theta" , xlab=bquote(theta) , 
          cenTend="median" , compVal=0.5 , ROPE=c(0.45,0.55) , credMass=0.90 )
saveGraph( file=paste0(fileNameRoot,"ThetaPost2") , type="eps" )

#################################################
# now with STAN
#################################################
library("rstan")
setwd("~/Desktop/bayesian_assignments/assignment_4")

# The data for 2.a
N <- 14
Y = c(1,0,1,1,0,1,1,1,0,1,1,1,1,1)

# The data for 2.b
N <- 100
Y <- rnorm(N, 0.5, sd=1)
hist(Y)

# compile model
model <- stan_model("first_model.stan")

#pass data to stan and run model
parallel::detectCores() # make sure you have the cores on your machine to run the chains in parallel.
options(mc.cores=4) # run 4 parallel cores for 4 chains, speeds the process
fit <- sampling(model, list(N=N, Y=Y), iter=200, chains=4)

#diagnose
print(fit)
# by default, stan uses half of the iterations for warmup. this adds up to 400 posterior samples.
# the sampler is NUTS (no u turn sampling).

#graph

parameters <- extract(fit) #in the global environment, check parameters has 3 entries: a, b, and lp__ (low probability)
hist(parameters$a)

library("shinystan")
launch_shinystan(fit) #this will load a new window in the default browser

##################
# The data for 2.b¨
##################
N <- 10
z <- c(1,0,0,0,0,0,0,1,1,0)
hist(z)
#the probability should be quite low give that there are more 0's in this sample compared to sample y.

#Let's make another model for parameter z with stan.
model2 <- stan_model("second_model.stan")

#pass data to stan and run model
parallel::detectCores() # make sure you have the cores on your machine to run the chains in parallel.
options(mc.cores=4) # run 4 parallel cores for 4 chains, speeds the process
fit2 <- sampling(model2, list(N=N, z=z), iter=500, chains=4)
# I needed to add more iterations, because it didn't work the first time and got stuck.

#diagnose
print(fit2)
#Rhat looks okay since close to 1.

#graph

parameters2 <- extract(fit2) #in the global environment, check parameters has 3 entries: a, b, and lp__ (low probability)
hist(parameters2$a)

#The distributions look quite different. Now how to test if they are different or not?
#we have to take the difference between the probability distributions of y and z. if the answer is close to zero, they are the same.
#2.b.2

#This is WRONG. because there is no reason to change the priors. we just want to compare the probability distributions.
#Let's make a 3rd model for the difference between y and z, called delta.
N <- 10  #not sure how to deal with the size...
delta <- setdiff(Y, z) #to find the difference between two vectors in r.

#This returns with delta as an empty numeric. that can't be right...
#This was incorrect! Now to find the difference between two probability distributions:

delta <- setdiff(parameters$a, parameters2$a)

#2.b.3
hist(delta)

#Okay that looks better. Now we can see in the histogram that the distribution is not centered around zero.

#Central tendency
mode(delta)
mean(delta)
#The mode is 1, which means there are more 1's in the difference than zeros. This means that more than 50% of the values are different between the two distributions.
#The mean is 0.33 which means that there is a non-zero chance that the distributions are different.

#credible interval

#Let's put it into ShinyStan. I will treat delta as a prior, because I don't know how to treat it as a posterior.

#Let's make another model for parameter delta with stan.
model_delta <- stan_model("delta_model.stan")

#pass data to stan and run model
fitdelta <- sampling(model_delta, list(N=N, delta=delta), iter=500, chains=1)
N=396
# I needed to add more iterations, because it didn't work the first time and got stuck.

#diagnose
print(fitdelta)
#Rhat looks okay since close to 1.

#graph

parametersDelta <- extract(fitdelta) #in the global environment, check parameters has 3 entries: a, b, and lp__ (low probability)
hist(parametersDelta$a)


