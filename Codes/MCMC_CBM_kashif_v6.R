# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (October 2016)
# k.mahmud@westernsydney.edu.au

# This code carries out Bayesian calibration for 4 variables (allocation fractions: "k","af","as","sf") on 
# daily time scale (e.g. 120 days) to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot)

##############################
# Version = v6: Removing some data points randomly to test the code with gap data
##############################
rm(list=ls())
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM/Results")

# install.packages("mvtnorm")
library(mvtnorm) # Creates candidate parameter vector as a multivariate normal jump away from the current candidate
chainLength = 1000 # Setting the length of the Markov Chain to be generated

# Generate synthetic data for Cstorage,Cleaf,Cstem,Croot with Mean and SD
time = 1:30 # Time in days

# Function to gerenate normally distributed random numbers with mean and SD
rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }
GPP.data = rnorm2(length(time),1,0.4) # Generate GPP data with mean=15, sd=3  # Units g C d-1
GPP.data = GPP.data[order(GPP.data)]
Rd.data = rnorm2(length(time),0.04,0.008) # Generate Rd data with mean=0.04, sd=0.008  # Units g C g-1 C d-1
Y = 0.3  # Allocation fraction to growth respiration (Fromm literature: Y ~ 0.3)

# Generate random parameter sets for synthetic Cleaf,Cstem,Croot data generation
k.data = rnorm2(length(time),0.5,0.05) # Generate data 'k' values with mean=0.5, sd=0.05 
af.data = rnorm2(length(time),1/3,0.1) # Generate data 'af' values with mean=1/3, sd=0.1 
as.data = rnorm2(length(time),1/3,0.1) # Generate data 'as' values with mean=1/3, sd=0.1
sf.data = rnorm2(length(time),1/50,1/500) # Generate data 'sf' values with mean=1/50, sd=1/500

Cstorage.data <- Cleaf.data <- Croot.data <- Cstem.data <- c()
Cstorage.data[1] <- 0.5 # Units g C d-1
Cleaf.data[1] <- 2 # Units g C d-1
Cstem.data[1] <- 1.5 # Units g C d-1
Croot.data[1] <- 1 # Units g C d-1

# Defining the model to iteratively calculate Cstorage, Cleaf, Cstem, Croot
model <- function (GPP,Rd,Cstorage,Cleaf,Cstem,Croot,Y,k,af,as,sf) {
  for (i in 2:length(GPP)){
    Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rd[i-1]*(Cleaf[i-1] + Croot[i-1] + Cstem[i-1]) - k[i-1]*Cstorage[i-1]
    Cleaf[i] <- Cleaf[i-1] + k[i-1]*Cstorage[i-1]*af[i-1]*(1-Y) - sf[i-1]*Cleaf[i-1]
    Cstem[i] <- Cstem[i-1] + k[i-1]*Cstorage[i-1]*as[i-1]*(1-Y) 
    Croot[i] <- Croot[i-1] + k[i-1]*Cstorage[i-1]*(1-af[i-1]-as[i-1])*(1-Y) 
  }
  output = data.frame(Cstorage,Cleaf,Cstem,Croot)
  return(output)}

# Generating the synthetic data sets
output.data = model(GPP.data,Rd.data,Cstorage.data,Cleaf.data,Cstem.data,Croot.data,Y,k.data,af.data,as.data,sf.data)
names(output.data) = c("Cstorage.data","Cleaf.data","Cstem.data","Croot.data")

# Removing some data points randomly
library(purrr)
output.gap.data = map_df(output.data, function(x) {x[sample(c(TRUE, NA), prob = c(0.2, 0.8), size = length(x), replace = TRUE)]})
output.gap.data[1,] = output.data[1,]
output.data = output.gap.data

# Plotting the synthetic data sets
matplot(output.data[ , 2:ncol(output.data)],type = c("b"),pch=1,col = 1:3,xlab="Days",ylab="gC",main="Different Carbon pool measurements") #plot
legend("topleft", legend = c("Cleaf.data","Cstem.data","Croot.data"), col=1:3, pch=0.75) # optional legend

# Initialize SD of data sets
sd.Cleaf = 0.1
sd.Cstem = 0.1
sd.Croot = 0.1

# Setting lower and upper bounds of the prior parameter pdf, and starting point of the chain
no.param = length(time)
no.var = 4 # variables are k,af,as,sf
param.k <- matrix(c(0.1,0.5,1) , nrow=no.param, ncol=3, byrow=T) 
param.af <- matrix(c(0.1,1/3,0.6) , nrow=no.param, ncol=3, byrow=T) 
param.as <- matrix(c(0.1,1/3,0.6) , nrow=no.param, ncol=3, byrow=T) 
param.sf <- matrix(c(0,1/50,1/25) , nrow=no.param, ncol=3, byrow=T) 
param = data.frame(param.k,param.af,param.as,param.sf)
names(param) <- c("k_min", "k", "k_max", "af_min", "af", "af_max","as_min", "as", "as_max", "sf_min", "sf", "sf_max")
pMinima <- param[ ,c("k_min", "af_min", "as_min", "sf_min")]
pMaxima <- param[ ,c("k_max", "af_max", "as_max", "sf_max")]
pValues <- param[ ,c("k","af","as","sf")] # Starting point of the chain
pChain <- matrix(0, nrow=chainLength, ncol = no.param*no.var+1) # Initialising the chain


# Defining the variance-covariance matrix for proposal generation
vcovProposal = diag( (0.5*(pMaxima-pMinima)) ^2 ) 


# Find the Prior probability density
prior.dist = vector("list", no.var)
for (i in 1:no.var)
{prior.dist[i] = list(log(dunif(pValues[ , i], pMinima[ , i], pMaxima[ , i])))}
logPrior0 <- sum(unlist(prior.dist))


# Calculating model outputs for the starting point of the chain
Cstorage <- Cleaf <- Croot <- Cstem <- c()
Cstorage[1] <- 0.5
Cleaf[1] <- 2
Cstem[1] <- 1.5
Croot[1] <- 1
output = model(GPP.data,Rd.data,Cstorage,Cleaf,Cstem,Croot,Y,pValues$k,pValues$af,pValues$as,pValues$sf)


# Calculating the log likelihood of starting point of the chain
logli <- matrix(0, nrow=length(time), ncol = 1) # Initialising the logli
for (i in 1:length(time)) {
  if (!is.na(output.data$Cleaf.data[i]))
  {logli[i] = - 0.5*((output$Cleaf[i] - output.data$Cleaf.data[i])/sd.Cleaf)^2 - log(sd.Cleaf)}
  if (!is.na(output.data$Cstem.data[i]))
  {logli[i] = logli[i] - 0.5*((output$Cstem[i] - output.data$Cstem.data[i])/sd.Cstem)^2 - log(sd.Cstem)}
  if (!is.na(output.data$Croot.data[i]))
  {logli[i] = logli[i] - 0.5*((output$Croot[i] - output.data$Croot.data[i])/sd.Croot)^2 - log(sd.Croot)}
}
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
  Cstorage <- Cleaf <- Croot <- Cstem <- c()
  Cstorage[1] <- 0.5
  Cleaf[1] <- 2
  Cstem[1] <- 1.5
  Croot[1] <- 1
  out.cand = model(GPP.data,Rd.data,Cstorage,Cleaf,Cstem,Croot,Y,
                   candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf)
  
  logli <- matrix(0, nrow=length(time), ncol = 1) # Initialising the logli
  for (i in 1:length(time)) {
    if (!is.na(output.data$Cleaf.data[i]))
    {logli[i] = - 0.5*((out.cand$Cleaf[i] - output.data$Cleaf.data[i])/sd.Cleaf)^2 - log(sd.Cleaf)}
    if (!is.na(output.data$Cstem.data[i]))
    {logli[i] = logli[i] - 0.5*((out.cand$Cstem[i] - output.data$Cstem.data[i])/sd.Cstem)^2 - log(sd.Cstem)}
    if (!is.na(output.data$Croot.data[i]))
    {logli[i] = logli[i] - 0.5*((out.cand$Croot[i] - output.data$Croot.data[i])/sd.Croot)^2 - log(sd.Croot)}
  }
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
output.final = model(GPP.data,Rd.data,Cstorage,Cleaf,Cstem,Croot,Y,
                     param.final$k,param.final$af,param.final$as,param.final$sf)


# Plotting the Measured(data) vs Modelled Plant Carbon pools for comparison
output.final$time = time
# output.data$GPP.Cum.data = cumsum(GPP.data)
# output.data$Resp.Cum.data = cumsum(Rd.data * (output.data$Cleaf.data + output.data$Cstem.data + output.data$Croot.data))
output.data$time = time
names(output.final) = c("Cstorage.modelled","Cleaf.modelled","Cstem.modelled","Croot.modelled","time")
melted.output = melt(output.final, id.vars="time")
melted.data = melt(output.data[ , c(2:ncol(output.data))], id.vars="time")

# png(file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/Measured_vs_modelled_Cpools.png")
ggplot(melted.data, aes(x = time, y = value, group = variable, colour=factor(variable))) + 
  geom_point(pch=15) +
  geom_line(data = melted.output, aes(x = time, y = value, group = variable, colour=factor(variable))) + 
  xlab("Days") +
  ylab("Plant Carbon pool (gC)") +
  ggtitle("Measured (points) vs Modelled (lines) Plant Carbon pools")
# dev.off()

# Plotting the parameter sets over time
param.final$time = time
param.final$ar = 1 - param.final$af - param.final$as
melted.param = melt(param.final, id.vars="time")
# png(file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/Allocation_fractions_over_time.png")
ggplot() + 
  geom_line(data = melted.param, aes(x = time, y = value, group = variable, colour=factor(variable))) + 
  xlab("Days") +
  ylab("Parameters") +
  ggtitle("Modelled allocation fractions")
# dev.off()

# Acceptance rate of the chain
nAccepted = length(unique(pChain[,1]))
acceptance = (paste(nAccepted, "out of ", chainLength, "candidates accepted ( = ",
                    round(100*nAccepted/chainLength), "%)"))
print(acceptance)


# Find the correlation coefficients (r2) between original measurements and predictions 
corrMatrix.1 = cor(output.final$Cleaf[!is.na(output.data$Cleaf.data)],output.data$Cleaf.data[!is.na(output.data$Cleaf.data)])
t1 = (paste("Correlation Coefficient, r2 of original Cleaf measurements and predictions: ", corrMatrix.1*corrMatrix.1))
corrMatrix.2 = cor(output.final$Cstem[!is.na(output.data$Cstem.data)],output.data$Cstem.data[!is.na(output.data$Cstem.data)])
t2 = (paste("Correlation Coefficient, r2 of original Cstem measurements and predictions: ", corrMatrix.2*corrMatrix.2))
corrMatrix.3 = cor(output.final$Croot[!is.na(output.data$Croot.data)],output.data$Croot.data[!is.na(output.data$Croot.data)])
t3 = (paste("Correlation Coefficient, r2 of original Croot measurements and predictions: ", corrMatrix.3*corrMatrix.3))
print(t1)
print(t2)
print(t3)


# Find the standard error (SE) between original measurements and predictions 
se <- function(data,meas) { se = SE = 0
  for (i in 1:length(meas)) {
    N = sum(!is.na(data))
    if (!is.na(data[i]))
    {se[i] = (data[i] - meas[i])^2 / N }
  }
  SE = sum(se,na.rm=TRUE) 
  return(SE)}
se1 = se(output.data$Cleaf.data,output.final$Cleaf)
t1 = (paste("Standard error (SE) of Cleaf measurements: ", se1))
se2 = se(output.data$Cstem.data,output.final$Cstem)
t2 = (paste("Standard error (SE) of Cstem measurements: ", se2))
se3 = se(output.data$Croot.data,output.final$Croot)
t3 = (paste("Standard error (SE) of Croot measurements: ", se3))
print(t1)
print(t2)
print(t3)




