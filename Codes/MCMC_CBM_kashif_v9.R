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
library(reshape2)
library(ggplot2)
chainLength = 1000 # Setting the length of the Markov Chain to be generated
vol = 1000

# Import daily GPP, daily Rd
GPP.data = read.csv("GPP.csv") # Units gC d-1
GPP.data = subset(GPP.data,volume==vol) # Consider only free seedling to start with
names(GPP.data)[3] = "GPP"
GPP.data$GPP = GPP.data$GPP
Rd.data = read.csv("Rd.csv") # Units g C g-1 plant d-1
Rd.data = subset(Rd.data,volume==vol)
tnc.data = read.csv("tnc_fortnightly_data.csv") # Units g plant-1
Sleaf.data = tnc.data = subset(tnc.data,volume==vol)

# Import weekly Cleaf, weekly Cstem, initial/harvest Croot data with Mean and SD
Mleaf.data = read.csv("Cleaf_weekly_data.csv") # Units gC d-1
Mleaf.data = subset(Mleaf.data,volume==vol)
Mstem.data = read.csv("Cstem_weekly_data.csv") # Units gC d-1
Mstem.data = subset(Mstem.data,volume==vol)
Mroot.data = read.csv("Croot_twice_data.csv") # Units gC d-1
Mroot.data = subset(Mroot.data,volume==vol)
############## Interpolate weekly root mass (considering root growth similar to the stem growth)
# root_to_stem = Mroot.data$rootmass[nrow(Mroot.data)] / Mstem.data$stemmass[nrow(Mstem.data)]
# mult = seq(1, root_to_stem, (root_to_stem-1)/(nrow(Mstem.data)-1))
# # mult = seq(Mroot.data$rootmass[1], Mroot.data$rootmass[nrow(Mroot.data)], (Mroot.data$rootmass[nrow(Mroot.data)] -
#                                                                                       # Mroot.data$rootmass[1]) / (nrow(Mstem.data)-1))
# mult_SD = seq(Mroot.data$rootmass_SD[1], Mroot.data$rootmass_SD[nrow(Mroot.data)], (Mroot.data$rootmass_SD[nrow(Mroot.data)] -
#                                                                                       Mroot.data$rootmass_SD[1]) / (nrow(Mstem.data)-1))
# Mroot.data = data.frame(Mstem.data$volume,Mstem.data$Date,Mstem.data$stemmass * mult,mult_SD)
# # Mroot.data = data.frame(Mstem.data$volume,Mstem.data$Date,mult,mult_SD)
# names(Mroot.data)[1:ncol(Mroot.data)] = c("volume","Date","rootmass","rootmass_SD")
##############

# Merge all GPP, Rd, Cleaf, Cstem, Croot data
data = merge(GPP.data,Rd.data, all = TRUE)
data = merge(data,Sleaf.data, all = TRUE)
data = merge(data,Mleaf.data, all = TRUE)
data = merge(data,Mstem.data, all = TRUE)
data = merge(data,Mroot.data, all = TRUE)
names(data)[4:ncol(data)] = c("Rd","Sleaf","Sleaf_SD","Mleaf","Mleaf_SD","Mstem","Mstem_SD","Mroot","Mroot_SD")
data[ , c(7:ncol(data))] = data[ , c(7:ncol(data))] * 0.65 # Unit conversion: gDM to gC

# Plotting the data sets
matplot(data[ , c("Mleaf","Mstem","Mroot")],type = c("b"),pch=1,col = 1:3,xlab="Days",ylab="gC",main="Different Carbon pool measurements") #plot
legend("topleft", legend = c("Mleaf.data","Mstem.data","Mroot.data"), col=1:3, pch=0.75) # optional legend


# # Estimate the Y value from total C coming in (GPP) and total C going out (sum of all Cpools, respiration)
# C.in = gpp.sum = sum(data$GPP)
# resp.sum = sum(data$Rd) * (mean(data$Mleaf,na.rm = TRUE) + mean(data$Mstem,na.rm = TRUE) + mean(data$Mroot,na.rm = TRUE))
# Cpool.sum = data$Mleaf[nrow(data)] + data$Mstem[nrow(data)] + data$Mroot[nrow(data)]
# Clit = sum(param.final$sf) * (mean(data$Mleaf,na.rm = TRUE))
# Y = 1 - (resp.sum + (Cpool.sum + Clit)) / C.in

Y = 0.357 # Fixed value of Y

# Defining the model to iteratively calculate Cstorage, Cleaf, Cstem, Croot, Sleaf, Sstem, Sroot
model <- function (GPP,Rd,Mleaf,Mstem,Mroot,Y,k,af,as,sf) {
  Cstorage = Sleaf = Sstem = Sroot = c()
  
  # From Duan's experiment for TNC partitioning to tree organs
  # Leaf TNC/Leaf DW =  0.1401421; Stem TNC/Stem DW =  0.0453869; Root TNC/Root DW =  0.02154037
  Sleaf[1] = Mleaf[1] / 0.65 * 0.1401421
  Sstem[1] = Mstem[1] / 0.65 * 0.0453869
  Sroot[1] = Mroot[1] / 0.65 * 0.02154037
  Cstorage[1] <- Sleaf[1] + Sstem[1] + Sroot[1] 
  
  Cleaf <- Croot <- Cstem <- c()
  Cleaf[1] <- Mleaf[1] - Sleaf[1]
  Cstem[1] <- Mstem[1] - Sstem[1]
  Croot[1] <- Mroot[1] - Sroot[1]
  
  for (i in 2:length(GPP)){
    Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1]) - k[i-1]*Cstorage[i-1]
    Sleaf[i] <- Cstorage[i] * 0.75 # 75% of storage goes to leaf (Duan's experiment)
    Sstem[i] <- Cstorage[i] * 0.16 # 16% of storage goes to stem (Duan's experiment)
    Sroot[i] <- Cstorage[i] * 0.09 # 9% of storage goes to root (Duan's experiment)
    
    Cleaf[i] <- Cleaf[i-1] + k[i-1]*Cstorage[i-1]*af[i-1]*(1-Y[i-1]) - sf[i-1]*Cleaf[i-1]
    Cstem[i] <- Cstem[i-1] + k[i-1]*Cstorage[i-1]*as[i-1]*(1-Y[i-1])
    Croot[i] <- Croot[i-1] + k[i-1]*Cstorage[i-1]*(1-af[i-1]-as[i-1])*(1-Y[i-1])
    
    Mleaf[i] <- Cleaf[i] + Sleaf[i]
    Mstem[i] <- Cstem[i] + Sstem[i]
    Mroot[i] <- Croot[i] + Sroot[i]
    
    # Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rd[i-1]*(Cleaf[i-1] + Croot[i-1] + Cstem[i-1]) - k[i-1]*Cstorage[i-1]
    # Cleaf[i] <- Cleaf[i-1] + k[i-1]*Cstorage[i-1]*af[i-1]*(1-Y) - sf[i-1]*Cleaf[i-1]
    # Cstem[i] <- Cstem[i-1] + k[i-1]*Cstorage[i-1]*as[i-1]*(1-Y)
    # Croot[i] <- Croot[i-1] + k[i-1]*Cstorage[i-1]*(1-af[i-1]-as[i-1])*(1-Y)
  }
  output = data.frame(Cstorage,Mleaf,Mstem,Mroot,Sleaf)
  return(output)}


# Setting lower and upper bounds of the prior parameter pdf, and starting point of the chain
no.param = length(GPP.data$Date)
no.var = 4 # variables are k,af,as,sf,Y
param.k <- matrix(c(0.1,0.4,0.8) , nrow=no.param, ncol=3, byrow=T) 
param.af <- matrix(c(0.1,0.45,0.7) , nrow=no.param, ncol=3, byrow=T) 
param.as <- matrix(c(0.1,0.25,0.5) , nrow=no.param, ncol=3, byrow=T) 
param.sf <- matrix(c(0,1/100,1/50) , nrow=no.param, ncol=3, byrow=T) 
# param.Y <- matrix(c(0.2,0.35,0.5) , nrow=no.param, ncol=3, byrow=T) 
param = data.frame(param.k,param.af,param.as,param.sf)
names(param) <- c("k_min", "k", "k_max", "af_min", "af", "af_max","as_min", "as", "as_max", 
                  "sf_min", "sf", "sf_max")
pMinima <- param[ ,c("k_min", "af_min", "as_min", "sf_min")]
pMaxima <- param[ ,c("k_max", "af_max", "as_max", "sf_max")]
pValues <- param[ ,c("k","af","as","sf")] # Starting point of the chain
# param = data.frame(param.k,param.af,param.as,param.sf,param.Y)
# names(param) <- c("k_min", "k", "k_max", "af_min", "af", "af_max","as_min", "as", "as_max", 
#                   "sf_min", "sf", "sf_max","Y_min", "Y", "Y_max")
# pMinima <- param[ ,c("k_min", "af_min", "as_min", "sf_min", "Y_min")]
# pMaxima <- param[ ,c("k_max", "af_max", "as_max", "sf_max", "Y_max")]
# pValues <- param[ ,c("k","af","as","sf","Y")] # Starting point of the chain
pChain <- matrix(0, nrow=chainLength, ncol = no.param*no.var+1) # Initialising the chain


# Defining the variance-covariance matrix for proposal generation
vcovProposal = diag( (0.1*(pMaxima-pMinima)) ^2 ) # The higher the coefficient, the lower the acceptance rate with better matching


# Find the Prior probability density
prior.dist = vector("list", no.var)
for (i in 1:no.var)
{prior.dist[i] = list(log(dunif(pValues[ , i], pMinima[ , i], pMaxima[ , i])))}
logPrior0 <- sum(unlist(prior.dist))


# Calculating model outputs for the starting point of the chain
Mleaf = Mstem = Mroot = c()
Mleaf[1] <- data$Mleaf[1]
Mstem[1] <- data$Mstem[1]
Mroot[1] <- data$Mroot[1]
output = model(data$GPP,data$Rd,Mleaf,Mstem,Mroot,Y,pValues$k,pValues$af,pValues$as,pValues$sf)
# output = model(data$GPP,data$Rd,Mleaf,Mstem,Mroot,pValues$Y,pValues$k,pValues$af,pValues$as,pValues$sf)


# Calculating the log likelihood of starting point of the chain
logli <- matrix(0, nrow=length(GPP.data$Date), ncol = 1) # Initialising the logli
for (i in 1:length(GPP.data$Date)) {
  if (!is.na(data$Mleaf[i]))
  {logli[i] = - 0.5*((output$Mleaf[i] - data$Mleaf[i])/data$Mleaf_SD[i])^2 - log(data$Mleaf_SD[i])}
  if (!is.na(data$Mstem[i]))
  {logli[i] = logli[i] - 0.5*((output$Mstem[i] - data$Mstem[i])/data$Mstem_SD[i])^2 - log(data$Mstem_SD[i])}
  if (!is.na(data$Mroot[i]))
  {logli[i] = logli[i] - 0.5*((output$Mroot[i] - data$Mroot[i])/data$Mroot_SD[i])^2 - log(data$Mroot_SD[i])}
  if (!is.na(data$Sleaf[i]))
  {logli[i] = logli[i] - 0.5*((output$Sleaf[i] - data$Sleaf[i])/data$Sleaf_SD[i])^2 - log(data$Sleaf_SD[i])}
}
logL0 <- sum(logli) # Log likelihood
pChain[1,] <- c(pValues$k,pValues$af,pValues$as,pValues$sf,logL0) # Assign the first parameter set with log likelihood
# pChain[1,] <- c(pValues$k,pValues$af,pValues$as,pValues$sf,pValues$Y,logL0) # Assign the first parameter set with log likelihood


# Calculating the next candidate parameter vector, as a multivariate normal jump away from the current point
for (c in (2 : chainLength))
{candidatepValues = matrix(ncol = no.var, nrow = no.param)
for (i in 1:no.var)
{candidatepValues[ , i] = rmvnorm(n=1, mean=pValues[ , i],
                                  sigma=diag(vcovProposal[i],no.param)) }
candidatepValues = data.frame(candidatepValues)
names(candidatepValues) <- c("k", "af", "as", "sf")
# names(candidatepValues) <- c("k", "af", "as", "sf","Y")


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
  Mleaf = Mstem = Mroot = c()
  Mleaf[1] <- data$Mleaf[1]
  Mstem[1] <- data$Mstem[1]
  Mroot[1] <- data$Mroot[1]
  out.cand = model(data$GPP,data$Rd,Mleaf,Mstem,Mroot,Y,
                   candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf)
  # out.cand = model(data$GPP,data$Rd,Mleaf,Mstem,Mroot,candidatepValues$Y,
  #                  candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf)
  
  logli <- matrix(0, nrow=length(GPP.data$Date), ncol = 1) # Initialising the logli
  for (i in 1:length(GPP.data$Date)) {
    if (!is.na(data$Mleaf[i]))
    {logli[i] = - 0.5*((out.cand$Mleaf[i] - data$Mleaf[i])/data$Mleaf_SD[i])^2 - log(data$Mleaf_SD[i])}
    if (!is.na(data$Mstem[i]))
    {logli[i] = logli[i] - 0.5*((out.cand$Mstem[i] - data$Mstem[i])/data$Mstem_SD[i])^2 - log(data$Mstem_SD[i])}
    if (!is.na(data$Mroot[i]))
    {logli[i] = logli[i] - 0.5*((out.cand$Mroot[i] - data$Mroot[i])/data$Mroot_SD[i])^2 - log(data$Mroot_SD[i])}
    if (!is.na(data$Sleaf[i]))
    {logli[i] = logli[i] - 0.5*((out.cand$Sleaf[i] - data$Sleaf[i])/data$Sleaf_SD[i])^2 - log(data$Sleaf_SD[i])}
  }
  logL1 <- sum(logli)
  # Calculating the logarithm of the Metropolis ratio
  logalpha <- (log(Prior1)+logL1) - (logPrior0+logL0) 
  # Accepting or rejecting the candidate vector
  # browser()
  if ( log(runif(1, min = 0, max =1)) < logalpha ){ 
    pValues <- candidatepValues
    logPrior0 <- log(Prior1)
    logL0 <- logL1
  }
}
pChain[c,] <- c(pValues$k,pValues$af,pValues$as,pValues$sf,logL0)
# pChain[c,] <- c(pValues$k,pValues$af,pValues$as,pValues$sf,pValues$Y,logL0)
}

# Store the final parameter set values
param.set = colMeans(pChain[ , 1:(no.param*no.var)])
param.final = data.frame(matrix(ncol = no.var, nrow = no.param))
names(param.final) <- c("k", "af", "as", "sf")
param.final$k = param.set[1:no.param]
param.final$af = param.set[(1+no.param):(2*no.param)]
param.final$as = param.set[(1+2*no.param):(3*no.param)]
param.final$sf = param.set[(1+3*no.param):(4*no.param)]
# param.final$Y = param.set[(1+4*no.param):(5*no.param)]

# Calculate final output set from the predicted parameter set
Mleaf = Mstem = Mroot = c()
Mleaf[1] <- data$Mleaf[1]
Mstem[1] <- data$Mstem[1]
Mroot[1] <- data$Mroot[1]
output.final = model(data$GPP,data$Rd,Mleaf,Mstem,Mroot,Y,
                     param.final$k,param.final$af,param.final$as,param.final$sf)
# output.final = model(data$GPP,data$Rd,Mleaf,Mstem,Mroot,param.final$Y,
#                      param.final$k,param.final$af,param.final$as,param.final$sf)

# Plotting the Measured(data) vs Modelled Plant Carbon pools for comparison
output.final$Date = data$Date
# output.data$GPP.Cum.data = cumsum(GPP.data)
# output.data$Resp.Cum.data = cumsum(Rd.data * (output.data$Cleaf.data + output.data$Cstem.data + output.data$Croot.data))
# data$time = time
names(output.final) = c("Cstorage.modelled","Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled","Date")
melted.output = melt(output.final, id.vars="Date")
melted.data = melt(data[ , c("Mleaf","Mstem","Mroot","Sleaf","Date")], id.vars="Date")

# png(file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/Measured_vs_modelled_Cpools.png")
ggplot(melted.data, aes(x = Date, y = value, group = variable, colour=factor(variable))) + 
  geom_point(pch=15) +
  geom_line(data = melted.output, aes(x = Date, y = value, group = variable, colour=factor(variable))) + 
  xlab("Days") +
  ylab("Plant Carbon pool (gC)") +
  ggtitle("Measured (points) vs Modelled (lines) Plant Carbon pools")
# dev.off()

# Plotting the parameter sets over time
param.final$Date = data$Date
param.final$ar = 1 - param.final$af - param.final$as
melted.param = melt(param.final, id.vars="Date")
# png(file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/Allocation_fractions_over_time.png")
ggplot() + 
  geom_line(data = melted.param, aes(x = Date, y = value, group = variable, colour=factor(variable))) + 
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
corrMatrix.1 = cor(output.final$Mleaf[!is.na(data$Mleaf)],data$Mleaf[!is.na(data$Mleaf)])
t1 = (paste("Correlation Coefficient, r2 of original Mleaf measurements and predictions: ", corrMatrix.1*corrMatrix.1))
corrMatrix.2 = cor(output.final$Mstem[!is.na(data$Mstem)],data$Mstem[!is.na(data$Mstem)])
t2 = (paste("Correlation Coefficient, r2 of original Mstem measurements and predictions: ", corrMatrix.2*corrMatrix.2))
corrMatrix.3 = cor(output.final$Mroot[!is.na(data$Mroot)],data$Mroot[!is.na(data$Mroot)])
t3 = (paste("Correlation Coefficient, r2 of original Mroot measurements and predictions: ", corrMatrix.3*corrMatrix.3))
print(t1)
print(t2)
print(t3)


# Find the standard error (SE) between original measurements and predictions 
se <- function(data,res) { se = SE = 0
for (i in 1:length(res)) {
  
  if (!is.na(data[i]))
  {se[i] = (data[i] - res[i])^2 }
}
N = sum(!is.na(data))
SE = sum(se,na.rm=TRUE) / N / mean(data,na.rm=TRUE)
return(SE)}

se1 = se(data$Mleaf,output.final$Mleaf)
t1 = (paste("Normalized Standard error (SE) of Mleaf measurements: ", se1))
se2 = se(data$Mstem,output.final$Mstem)
t2 = (paste("Normalized Standard error (SE) of Mstem measurements: ", se2))
se3 = se(data$Mroot,output.final$Mroot)
t3 = (paste("Normalized Standard error (SE) of Mroot measurements: ", se3))
print(t1)
print(t2)
print(t3)


# Validating total C coming in (GPP) vs total C going out (sum of all Cpools, respiration and !!!storage pool!!!)
C.in = gpp.sum = sum(data$GPP)

resp.sum = sum(data$Rd) * (mean(data$Mleaf,na.rm = TRUE) + mean(data$Mstem,na.rm = TRUE) + mean(data$Mroot,na.rm = TRUE))
# Cstorage.sum = sum(output.final$Cstorage.modelled)
Cpool.sum = data$Mleaf[nrow(data)] + data$Mstem[nrow(data)] + data$Mroot[nrow(data)]
Clit = sum(param.final$sf) * (mean(data$Mleaf,na.rm = TRUE))
C.out = resp.sum + (Cpool.sum + Clit) / (1-Y)
Y = 1- (resp.sum + (Cpool.sum + Clit)) / C.in
# C.out = resp.sum + Cpool.sum + Clit + Cstorage.sum
cat('C.in = ',C.in)
cat('C.out = ',C.out)

# plot(data$Mleaf,col="red")
# lines(output.final$Mleaf,col="green")
# lines(output.final$Sleaf,col="blue")
lines(param.final$k,col="red")
lines(param.final$af,col="green")


