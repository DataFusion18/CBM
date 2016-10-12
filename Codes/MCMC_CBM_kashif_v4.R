# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (2016)

# This code carries out Bayesian calibration for 4 variables (allocation fractions: "k","af","as","sf") on 
# daily time scale (e.g. 30 days) to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot)
# The algorithm only considers the instantaneous time frame (Not the previous/cumulative pool amount)

rm(list=ls())
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM/Results")

# install.packages("mvtnorm")
library(mvtnorm)
chainLength = 1000 # Setting the length of the Markov Chain to be generated

# Generate synthetic data for Cstorage,Cleaf,Cstem,Croot with Mean and SD
time = 1:30 # Time in days
# Function to gerenate normally distributed random numbers with mean and SD
rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }
GPP = rnorm2(length(time),15,3) # Generate GPP data with mean=15, sd=3
Rd = rnorm2(length(time),4,0.8) # Generate Rd data with mean=4, sd=0.8
# Generate Cstorage.true to create a measurement sets of Cleaf, Cstem, Croot 
Cstorage.true = rnorm2(length(time),7.5,2) # Generate Cstorage data with mean=7.5, sd=2
Y = 0.3

# Generate random parameter sets for synthetic Cleaf,Cstem,Croot data generation
k.true = rnorm2(length(time),0.55,0.15) # Generate true 'k' values with mean=0.55, sd=0.15 
af.true = rnorm2(length(time),1/7,0.05) # Generate true 'af' values with mean=1/7, sd=0.05 
as.true = rnorm2(length(time),3/7,0.15) # Generate true 'as' values with mean=3/7, sd=0.15
sf.true = rnorm2(length(time),1/30,1/100) # Generate true 'sf' values with mean=1/30, sd=1/100

Cleaf.true <- matrix(, nrow=length(time), ncol=1)
Cstem.true <- matrix(, nrow=length(time), ncol=1)
Croot.true <- matrix(, nrow=length(time), ncol=1)
for (i in 1:length(GPP)){
  Cleaf.true[i] = (af.true[i] * (1-Y) * k.true[i]*Cstorage.true[i]) / (1 + sf.true[i])
  Cstem.true[i] = as.true[i] * (1-Y) * k.true[i]*Cstorage.true[i]
  Croot.true[i] = (1 - af.true[i] - as.true[i]) * (1-Y) * k.true[i]*Cstorage.true[i]
}
# Calculate SD of data sets
sd.Cstorage = sd(Cstorage.true)
sd.Cleaf = sd(Cleaf.true)
sd.Cstem = sd(Cstem.true)
sd.Croot = sd(Croot.true)

# Setting lower and upper bounds of the prior parameter pdf, and starting point of the chain
no.param = length(time)
no.var = 4 # variables are k,af,as,sf
param.k <- matrix(c(0.1,0.55,1) , nrow=no.param, ncol=3, byrow=T) 
param.af <- matrix(c(0,1/7,0.35) , nrow=no.param, ncol=3, byrow=T) 
param.as <- matrix(c(0.2,3/7,0.75) , nrow=no.param, ncol=3, byrow=T) 
param.sf <- matrix(c(0,1/30,1/10) , nrow=no.param, ncol=3, byrow=T) 
param = data.frame(param.k,param.af,param.as,param.sf)
names(param) <- c("k_min", "k", "k_max", "af_min", "af", "af_max","as_min", "as", "as_max", "sf_min", "sf", "sf_max")
pMinima <- param[ ,c("k_min", "af_min", "as_min", "sf_min")]
pMaxima <- param[ ,c("k_max", "af_max", "as_max", "sf_max")]

# Initialising the log likelihood 
logli <- matrix(, nrow=length(time), ncol=1) 

# Defining the variance-covariance matrix for proposal generation
vcovProposal = diag( (0.1*(pMaxima-pMinima)) ^2 ) 
# vcovProposal = diag( (0.05*(pMaxima$k_max-pMinima$k_min)) ^2 )

pValues <- param[ ,c("k","af","as","sf")] # Start point of the chain
pChain <- matrix(0, nrow=chainLength, ncol = no.param*no.var+1) # Initialising the chain

# Prior probability density
# logPrior0 <- sum(log(dunif(pValues, min=pMinima, max=pMaxima))) 
prior.dist = vector("list", no.var)
for (i in 1:no.var)
{prior.dist[i] = list(log(dunif(pValues[ , i], pMinima[ , i], pMaxima[ , i])))}
logPrior0 <- sum(unlist(prior.dist))

# Defining the model to iteratively calculate Cstorage, Cleaf, Cstem, Croot
model <- function (GPP,Rd,Cleaf.true,Cstem.true,Croot.true,Y,k,af,as,sf) {
  Cstorage <- matrix(, nrow=length(GPP), ncol=1)
  Cleaf <- matrix(, nrow=length(GPP), ncol=1)
  Cstem <- matrix(, nrow=length(GPP), ncol=1)
  Croot <- matrix(, nrow=length(GPP), ncol=1)
  output = data.frame(matrix(ncol = 4, nrow = length(GPP)))
  names(output)[1:4] <- c("Cstorage", "Cleaf", "Cstem", "Croot") 
  for (i in 1:length(GPP)){
    output$Cstorage[i] = GPP[i] - Rd[i]*(Cleaf.true[i] + Cstem.true[i] + Croot.true[i]) / (1 + k[i]) 
    output$Cleaf[i] = (af[i] * (1-Y) * k[i]*output$Cstorage[i]) / (1 + sf[i])
    output$Cstem[i] = as[i] * (1-Y) * k[i]*output$Cstorage[i]
    output$Croot[i] = (1 - af[i] - as[i]) * (1-Y) * k[i]*output$Cstorage[i]
  }
  return(output)}

# Calculating model outputs for the start point of the chain and then the likelihood
output = model(GPP,Rd,Cleaf.true,Cstem.true,Croot.true,Y,pValues$k,pValues$af,pValues$as,pValues$sf)
for (i in 1:length(time)) 
{logli[i] = - 0.5*((output$Cleaf[i] - Cleaf.true[i])/sd.Cleaf)^2
- 0.5*((output$Cstem[i] - Cstem.true[i])/sd.Cstem)^2 - 0.5*((output$Croot[i] - Croot.true[i])/sd.Croot)^2 
- log(sd.Cleaf) - log(sd.Cstem) - log(sd.Croot)}
# {logli[i] = -0.5*((output$Cstorage[i] - Cstorage.true[i])/sd.Cstorage)^2 - log(sd.Cstorage) - 0.5*((output$Cleaf[i] - Cleaf.true[i])/sd.Cleaf)^2
# - 0.5*((output$Cstem[i] - Cstem.true[i])/sd.Cstem)^2 - 0.5*((output$Croot[i] - Croot.true[i])/sd.Croot)^2 
# - log(sd.Cleaf) - log(sd.Cstem) - log(sd.Croot)} 
logL0 <- sum(logli) # Log likelihood
pChain[1,] <- c(pValues$k,pValues$af,pValues$as,pValues$sf,logL0) # Assign the first parameter set with log likelihood

# Calculating the next candidate parameter vector, as a multivariate normal jump away from the current point
for (c in (2 : chainLength))
{candidatepValues = matrix(ncol = no.var, nrow = no.param)
for (i in 1:no.var)
{candidatepValues[ , i] = rmvnorm(n=1, mean=pValues[ , i],
                                  sigma=diag(vcovProposal[i],no.param)) }
candidatepValues = data.frame(candidatepValues)
names(candidatepValues) <- c("k", "af", "as", "sf")

# Reflected back to generate another candidate value
reflectionFromMin = pmin( unlist(matrix(0,nrow=no.param,ncol=no.var)), unlist(candidatepValues-pMinima) )
reflectionFromMax = pmax( unlist(list(rep(0, no.param))), unlist(candidatepValues-pMaxima) )
candidatepValues = candidatepValues - 2 * reflectionFromMin - 2 * reflectionFromMax 

# Calculating the prior probability density for the candidate parameter vector
if (all(candidatepValues>pMinima) && all(candidatepValues<pMaxima)){
  uni.dist = vector("list", no.var)
  for (i in 1:no.var)
  {uni.dist[i] = list(dunif(candidatepValues[ , i], pMinima[ , i], pMaxima[ , i]))}
  # uni.dist[1+(i-1)*no.param : i*no.param] = list(dunif(candidatepValues[ , i], pMinima[ , i], pMaxima[ , i]))}
  Prior1 <- prod(unlist(uni.dist))
} else {Prior1 <- 0}

# Calculating the outputs for the candidate parameter vector, log likelihood
if (Prior1 > 0){
  out.cand = model(GPP,Rd,Cleaf.true,Cstem.true,Croot.true,Y,
                   candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf)
  for (i in 1:length(time))
  {logli[i] = -0.5*((out.cand$Cleaf[i] - Cleaf.true[i])/sd.Cleaf)^2
  - 0.5*((out.cand$Cstem[i] - Cstem.true[i])/sd.Cstem)^2 - 0.5*((out.cand$Croot[i] - Croot.true[i])/sd.Croot)^2 
  - log(sd.Cleaf) - log(sd.Cstem) - log(sd.Croot)} 
  # {logli[i] = -0.5*((out.cand$Cstorage[i] - Cstorage.true[i])/sd.Cstorage)^2 - log(sd.Cstorage) - 0.5*((out.cand$Cleaf[i] - Cleaf.true[i])/sd.Cleaf)^2
  # - 0.5*((out.cand$Cstem[i] - Cstem.true[i])/sd.Cstem)^2 - 0.5*((out.cand$Croot[i] - Croot.true[i])/sd.Croot)^2 
  # - log(sd.Cleaf) - log(sd.Cstem) - log(sd.Croot)} 
  logL1 <- sum(logli)
  # Calculating the logarithm of the Metropolis ratio
  logalpha <- (log(Prior1)+logL1) - (logPrior0+logL0) 
  # Accepting or rejecting the candidate vector
  if ( log(runif(1, min = 0, max =1)) < logalpha ){ 
    pValues <- candidatepValues
    logPrior0 <- log(Prior1)
    logL0 <- logL1
  }
}

pChain[c,] <- c(pValues$k,pValues$af,pValues$as,pValues$sf,logL0)
# pChain[c,1:no.param] <- pValues
# pChain[c,no.param*no.var+1] <- logL0
}

# Store the final parameter set values
param.set = colMeans(pChain[ , 1:(no.param*no.var)])
param.final = data.frame(matrix(ncol = no.var, nrow = no.param))
names(param.final) <- c("k", "af", "as", "sf")
param.final$k = param.set[1:no.param]
param.final$af = param.set[(1+no.param):(2*no.param)]
param.final$as = param.set[(1+2*no.param):(3*no.param)]
param.final$sf = param.set[(1+3*no.param):(4*no.param)]
# Calculate final output set from the predicted parameter set
output.final = model(GPP,Rd,Cleaf.true,Cstem.true,Croot.true,Y,
                     param.final$k,param.final$af,param.final$as,param.final$sf)

# Get the cumulative sums over the length of time
output.csum <- cumsum(output.final)
Cstorage.csum = cumsum(Cstorage.true)
Cleaf.csum = cumsum(Cleaf.true)
Cstem.csum = cumsum(Cstem.true)
Croot.csum = cumsum(Croot.true)
# Acceptance rate of the chain
nAccepted = length(unique(pChain[,1]))
acceptance = (paste(nAccepted, "out of ", chainLength, "candidates accepted ( = ",
                    round(100*nAccepted/chainLength), "%)"))
print(acceptance)

# mp <- apply(pChain, 2, mean)
# print(mp)
# pCovMatrix <- cov(pChain)
# print(pCovMatrix)

# Find the correlation coefficient between original measurements and predictions
corrMatrix.1 = cor(Cstorage.csum,output.csum$Cstorage)
t1 = (paste("Correlation Coefficient, r2 of original Cstorage measurements and predictions: ", corrMatrix.1*corrMatrix.1))
corrMatrix.2 = cor(Cleaf.csum,output.csum$Cleaf)
t2 = (paste("Correlation Coefficient, r2 of original Cleaf measurements and predictions: ", corrMatrix.2*corrMatrix.2))
corrMatrix.3 = cor(Cstem.csum,output.csum$Cstem)
t3 = (paste("Correlation Coefficient, r2 of original Cstem measurements and predictions: ", corrMatrix.3*corrMatrix.3))
corrMatrix.4 = cor(Croot.csum,output.csum$Croot)
t4 = (paste("Correlation Coefficient, r2 of original Croot measurements and predictions: ", corrMatrix.4*corrMatrix.4))
print(t1)
print(t2)
print(t3)
print(t4)

# Plot few accepted parameter values over time to find whether the chain converged
par(mfrow=c(1,1))                    
# png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/MCMC_CBM_v4_outcomes/Parameter_1.png")
plot(pChain[,no.param/2], main="Parameter 1 at time=10", xlab="Chain length", ylab="Parameter 1")
# dev.off()

# png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/MCMC_CBM_v4_outcomes/Parameter_2.png")
plot(pChain[,(no.param+no.param/2)], main="Parameter 2 at time=10", xlab="Chain length", ylab="Parameter 1")
# dev.off()

# png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/MCMC_CBM_v4_outcomes/Parameter_end.png")
plot(pChain[,no.param*no.var-no.param/2], main="Last Parameter at time=10", xlab="Chain length", ylab="Last Parameter")
# dev.off()

# plot(pChain[,1],pChain[,2])  
# title('Parameter 1 vs Parameter 2')
# plot(Cstorage.csum,output.csum$Cstorage) 
# lines(Cstorage.csum,Cstorage.csum)
# title('Original Cstorage measurements vs predictions')

# Plot original measurements vs predictions for Cstorage, Cleaf, Cstem, Croot
# png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/MCMC_CBM_v4_outcomes/Cstorage.png")
plot(time,Cstorage.csum, main="Original Cstorage measurements vs predictions", 
     xlab="Time", ylab="Cstorage", pch=20, col="red") 
lines(time,output.csum$Cstorage,lwd = 2,col="green")
legend('topleft', c("Measurements", "Predictions") , 
       lty=1, col=c('red','green'), bty='n', cex=0.75)
# dev.off()

# png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/MCMC_CBM_v4_outcomes/Cleaf.png")
plot(time,Cleaf.csum, main="Original Cleaf measurements vs predictions", 
     xlab="Time", ylab="Cleaf", pch=20, col="red") 
lines(time,output.csum$Cleaf,lwd = 2,col="green")
legend('topleft', c("Measurements", "Predictions") , 
       lty=1, col=c('red','green'), bty='n', cex=0.75)
# dev.off()

# png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/MCMC_CBM_v4_outcomes/Cstem.png")
plot(time,Cstem.csum, main="Original Cstem measurements vs predictions", 
     xlab="Time", ylab="Cstem", pch=20, col="red") 
lines(time,output.csum$Cstem,lwd = 2,col="green")
legend('topleft', c("Measurements", "Predictions") , 
       lty=1, col=c('red','green'), bty='n', cex=0.75)
# dev.off()

# png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/MCMC_CBM_v4_outcomes/Croot.png")
plot(time,Croot.csum, main="Original Croot measurements vs predictions", 
     xlab="Time", ylab="Croot", pch=20, col="red") 
lines(time,output.csum$Croot,lwd = 2,col="green")
legend('topleft', c("Measurements", "Predictions") , 
       lty=1, col=c('red','green'), bty='n', cex=0.75)
# dev.off()

# plot(Cleaf.true, output.final$Cleaf, main="Scatterplot Example", 
#      xlab="True Cleaf", ylab="Estimated Cleaf", pch=19)
# plot(Cleaf.true, output.final$Cleaf, main="Scatterplot Example", 
#      xlab="True Cleaf", ylab="Estimated Cleaf", pch=19)
# abline(lm(output.final$Cleaf~Cleaf.true), col="red") # regression line (y~x)
# plot(Cstem.true, output.final$Cstem, main="Scatterplot Example", 
#      xlab="True Cstem", ylab="Estimated Cstem", pch=19)
# abline(lm(output.final$Cstem~Cstem.true), col="red") # regression line (y~x)
# 
# lines(lowess(Cleaf.true, output.final$Cleaf), col="blue") # lowess line (x,y)
# 
# 
# 
# 
# 
# 
# t=paste("Parameter Number / True Parameter Values / Estimated Parameter Values:");
# print(t)
# for (i in 1:length(time)){
#   t = paste (i, ": ", k.true[i], "    ", format(round(param.final[i], 2),, nsmall = 2) );
#   print(t)
# }  

