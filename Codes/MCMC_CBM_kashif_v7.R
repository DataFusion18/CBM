# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (October 2016)
# k.mahmud@westernsydney.edu.au

# This code carries out Bayesian calibration for 4 variables (allocation fractions: "k","af","as","sf") on 
# daily time scale (e.g. 120 days) to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot)

##############################
# Version = v7: First trial with soil manipulation pot experiment data
##############################
rm(list=ls())
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM/Data_files")

# install.packages("mvtnorm")
library(mvtnorm) # Creates candidate parameter vector as a multivariate normal jump away from the current candidate
chainLength = 1000 # Setting the length of the Markov Chain to be generated
vol = 5

# Import daily GPP, daily Rd
GPP.data = read.csv("GPP.csv") # Units gC d-1
GPP.data = subset(GPP.data,volume==vol) # Consider only free seedling to start with
names(GPP.data)[3] = "GPP"
Rd.data = read.csv("Rd.csv") # Units g C g-1 plant d-1
Rd.data = subset(Rd.data,volume==vol)
Y = 0.3 # Carbon allocation fraction to growth respiration (Fromm literature: Y ~ 0.3)
time = 1:nrow(GPP.data) # Time in days

# Import weekly Cleaf, weekly Cstem, initial/harvest Croot data with Mean and SD
Cleaf.data = read.csv("Cleaf_weekly_data.csv") # Units gC d-1
Cleaf.data = subset(Cleaf.data,volume==vol)
Cstem.data = read.csv("Cstem_weekly_data.csv") # Units gC d-1
Cstem.data = subset(Cstem.data,volume==vol)
Croot.data = read.csv("Croot_twice_data.csv") # Units gC d-1
Croot.data = subset(Croot.data,volume==vol)

# Merge all GPP, Rd, Cleaf, Cstem, Croot data
data = merge(GPP.data,Rd.data, all = TRUE)
data = merge(data,Cleaf.data, all = TRUE)
data = merge(data,Cstem.data, all = TRUE)
data = merge(data,Croot.data, all = TRUE)
names(data)[4:ncol(data)] = c("Rd","Cleaf","Cleaf_SD","Cstem","Cstem_SD","Croot","Croot_SD")

# From Duan's experiment for TNC partitioning to tree organs
# Leaf TNC/Leaf DW =  0.1401421; Stem TNC/Stem DW =  0.0453869; Root TNC/Root DW =  0.02154037
Cstorage = c()
Cstorage[1] <- data$Cleaf[1] * 0.1401421 + data$Cstem[1] * 0.0453869 + data$Croot[1] * 0.02154037 

# Plotting the data sets
matplot(data[ , c("Cleaf","Cstem","Croot")],type = c("b"),pch=1,col = 1:3,xlab="Days",ylab="gC",main="Different Carbon pool measurements") #plot
legend("topleft", legend = c("Cleaf.data","Cstem.data","Croot.data"), col=1:3, pch=0.75) # optional legend


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


# Setting lower and upper bounds of the prior parameter pdf, and starting point of the chain
no.param = length(time)
no.var = 4 # variables are k,af,as,sf
param.k <- matrix(c(0.4,0.7,1) , nrow=no.param, ncol=3, byrow=T) 
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
vcovProposal = diag( (0.1*(pMaxima-pMinima)) ^2 ) 


# Find the Prior probability density
prior.dist = vector("list", no.var)
for (i in 1:no.var)
{prior.dist[i] = list(log(dunif(pValues[ , i], pMinima[ , i], pMaxima[ , i])))}
logPrior0 <- sum(unlist(prior.dist))


# Calculating model outputs for the starting point of the chain
Cleaf <- Croot <- Cstem <- c()
Cleaf[1] <- data$Cleaf[1]
Cstem[1] <- data$Cstem[1]
Croot[1] <- data$Croot[1]
output = model(data$GPP,data$Rd,Cstorage,Cleaf,Cstem,Croot,Y,pValues$k,pValues$af,pValues$as,pValues$sf)


# Calculating the log likelihood of starting point of the chain
logli <- matrix(0, nrow=length(time), ncol = 1) # Initialising the logli
for (i in 1:length(time)) {
  if (!is.na(data$Cleaf[i]))
  {logli[i] = - 0.5*((output$Cleaf[i] - data$Cleaf[i])/data$Cleaf_SD[i])^2 - log(data$Cleaf_SD[i])}
  if (!is.na(data$Cstem[i]))
  {logli[i] = logli[i] - 0.5*((output$Cstem[i] - data$Cstem[i])/data$Cstem_SD[i])^2 - log(data$Cstem_SD[i])}
  if (!is.na(data$Croot[i]))
  {logli[i] = logli[i] - 0.5*((output$Croot[i] - data$Croot[i])/data$Croot_SD[i])^2 - log(data$Croot_SD[i])}
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
  Cleaf <- Croot <- Cstem <- c()
  Cleaf[1] <- data$Cleaf[1]
  Cstem[1] <- data$Cstem[1]
  Croot[1] <- data$Croot[1]
  out.cand = model(data$GPP,data$Rd,Cstorage,Cleaf,Cstem,Croot,Y,
                   candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf)
  
  logli <- matrix(0, nrow=length(time), ncol = 1) # Initialising the logli
  for (i in 1:length(time)) {
    if (!is.na(data$Cleaf[i]))
    {logli[i] = - 0.5*((out.cand$Cleaf[i] - data$Cleaf[i])/data$Cleaf_SD[i])^2 - log(data$Cleaf_SD[i])}
    if (!is.na(data$Cstem[i]))
    {logli[i] = logli[i] - 0.5*((out.cand$Cstem[i] - data$Cstem[i])/data$Cstem_SD[i])^2 - log(data$Cstem_SD[i])}
    if (!is.na(data$Croot[i]))
    {logli[i] = logli[i] - 0.5*((out.cand$Croot[i] - data$Croot[i])/data$Croot_SD[i])^2 - log(data$Croot_SD[i])}
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
output.final = model(data$GPP,data$Rd,Cstorage,Cleaf,Cstem,Croot,Y,
                     param.final$k,param.final$af,param.final$as,param.final$sf)


# Plotting the Measured(data) vs Modelled Plant Carbon pools for comparison
output.final$time = time
# output.data$GPP.Cum.data = cumsum(GPP.data)
# output.data$Resp.Cum.data = cumsum(Rd.data * (output.data$Cleaf.data + output.data$Cstem.data + output.data$Croot.data))
data$time = time
names(output.final) = c("Cstorage.modelled","Cleaf.modelled","Cstem.modelled","Croot.modelled","time")
melted.output = melt(output.final, id.vars="time")
melted.data = melt(data[ , c("Cleaf","Cstem","Croot","time")], id.vars="time")

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
corrMatrix.1 = cor(output.final$Cleaf[!is.na(data$Cleaf)],data$Cleaf[!is.na(data$Cleaf)])
t1 = (paste("Correlation Coefficient, r2 of original Cleaf measurements and predictions: ", corrMatrix.1*corrMatrix.1))
corrMatrix.2 = cor(output.final$Cstem[!is.na(data$Cstem)],data$Cstem[!is.na(data$Cstem)])
t2 = (paste("Correlation Coefficient, r2 of original Cstem measurements and predictions: ", corrMatrix.2*corrMatrix.2))
corrMatrix.3 = cor(output.final$Croot[!is.na(data$Croot)],data$Croot[!is.na(data$Croot)])
t3 = (paste("Correlation Coefficient, r2 of original Croot measurements and predictions: ", corrMatrix.3*corrMatrix.3))
print(t1)
print(t2)
print(t3)


