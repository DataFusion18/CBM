# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (October 2016)
# k.mahmud@westernsydney.edu.au

# This code carries out Bayesian calibration for 4 variables (allocation fractions: "k","af","as","sf") on 
# daily time scale (e.g. 120 days) to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot)

##############################
# Version = v11: MCMC with soil manipulation pot experiment data for all treatments (including the free seedling), 
# This version considers either daily/weekly/monthly/just one parameter set for 5 variables ("k","Y","af","as","sf")
# So we can set the parameters for various time frames
# Also calculate the MCMC SDs for all parameters at different time steps
##############################
rm(list=ls())
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM/Data_files")

# install.packages("mvtnorm")
library(mvtnorm) # Creates candidate parameter vector as a multivariate normal jump away from the current candidate
library(reshape2)
library(ggplot2)
chainLength = 1500 # Setting the length of the Markov Chain to be generated
vol = 1000

# Import daily GPP, daily Rd
GPP.data = read.csv("GPP.csv") # Units gC d-1
GPP.data = subset(GPP.data,volume==vol) # Consider only free seedling to start with
names(GPP.data)[3] = "GPP"
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


# Merge all GPP, Rd, Cleaf, Cstem, Croot data
data = merge(GPP.data,Rd.data, all = TRUE)
data = merge(data,Sleaf.data, all = TRUE)
data = merge(data,Mleaf.data, all = TRUE)
data = merge(data,Mstem.data, all = TRUE)
data = merge(data,Mroot.data, all = TRUE)
names(data)[4:ncol(data)] = c("Rd","Sleaf","Sleaf_SD","Mleaf","Mleaf_SD","Mstem","Mstem_SD","Mroot","Mroot_SD")
data[ , c(7:ncol(data))] = data[ , c(7:ncol(data))] * 0.65 # Unit conversion: gDM to gC

# Reducing the measurements uncertainty (SDs) to fit the data perfectly
# data[,c("Sleaf_SD","Mleaf_SD","Mstem_SD","Mroot_SD")] = data[,c("Sleaf_SD","Mleaf_SD","Mstem_SD","Mroot_SD")] / 10000

# Plotting the data sets
par(mfrow=c(1,1))
matplot(data[ , c("Mleaf","Mstem","Mroot")],type = c("b"),pch=1,col = 1:3,xlab="Days",ylab="gC",main="Different Carbon pool measurements") #plot
legend("topleft", legend = c("Mleaf.data","Mstem.data","Mroot.data"), col=1:3, pch=0.75) # optional legend


# Defining the model to iteratively calculate Cstorage, Cleaf, Cstem, Croot, Sleaf, Sstem, Sroot
model <- function (GPP,Rd,j,Mleaf,Mstem,Mroot,Y,k,af,as,sf) {
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
  for (i in 2:length(GPP)) {
    Cstorage[i] <- Cstorage[i-1] + GPP[i-1] - Rd[i-1]*(Mleaf[i-1] + Mroot[i-1] + Mstem[i-1]) - k[(i-1)-(j[i-1])]*Cstorage[i-1]
    Sleaf[i] <- Cstorage[i] * 0.75 # 75% of storage goes to leaf (Duan's experiment)
    Sstem[i] <- Cstorage[i] * 0.16 # 16% of storage goes to stem (Duan's experiment)
    Sroot[i] <- Cstorage[i] * 0.09 # 9% of storage goes to root (Duan's experiment)
    
    Cleaf[i] <- Cleaf[i-1] + k[(i-1)-(j[i-1])]*Cstorage[i-1]*af[(i-1)-(j[i-1])]*(1-Y[(i-1)-(j[i-1])]) - sf[(i-1)-(j[i-1])]*Cleaf[i-1]
    Cstem[i] <- Cstem[i-1] + k[(i-1)-(j[i-1])]*Cstorage[i-1]*as[(i-1)-(j[i-1])]*(1-Y[(i-1)-(j[i-1])])
    Croot[i] <- Croot[i-1] + k[(i-1)-(j[i-1])]*Cstorage[i-1]*(1-af[(i-1)-(j[i-1])]-as[(i-1)-(j[i-1])])*(1-Y[(i-1)-(j[i-1])])
    
    Mleaf[i] <- Cleaf[i] + Sleaf[i]
    Mstem[i] <- Cstem[i] + Sstem[i]
    Mroot[i] <- Croot[i] + Sroot[i]
  }
  output = data.frame(Cstorage,Mleaf,Mstem,Mroot,Sleaf)
  return(output)
}


# Setting lower and upper bounds of the prior parameter pdf, and starting point of the chain
no.var = 5 # variables are k,Y,af,as,sf
no.param.par.var = 4 # parameter count per variable based on time frame 
param.vary = ceiling(nrow(data)/no.param.par.var) # How many days the parameter set remain unchanged (weekly = 7; monthly = 30; just one parameter = nrow(data))
# param.vary = 30 # How many days the parameter set remain unchanged (weekly = 7; monthly = 30; just one parameter = nrow(data))
no.param = ceiling(nrow(data)/param.vary) # number of parameter set for the whole duration of experiment (121 days)
j = c()
j[1] = 0
i = seq(1,nrow(data),1)
j[i] = i - ceiling(i/param.vary)*1  # j is for parameter settings for various time frames

param.k <- matrix(c(0,0.45,1) , nrow=no.param, ncol=3, byrow=T) 
param.Y <- matrix(c(0.2,0.3,0.4) , nrow=no.param, ncol=3, byrow=T) 
param.af <- matrix(c(0,0.45,0.7) , nrow=no.param, ncol=3, byrow=T) 
param.as <- matrix(c(0,0.17,0.5) , nrow=no.param, ncol=3, byrow=T) 
param.sf <- matrix(c(0,1/50,1/25) , nrow=no.param, ncol=3, byrow=T) 
param = data.frame(param.k,param.Y,param.af,param.as,param.sf)
names(param) <- c("k_min","k","k_max","Y_min","Y","Y_max","af_min","af","af_max","as_min","as","as_max","sf_min","sf","sf_max")
pMinima <- param[ ,c("k_min","Y_min","af_min","as_min","sf_min")]
pMaxima <- param[ ,c("k_max","Y_max","af_max","as_max","sf_max")]
pValues <- param[ ,c("k","Y","af","as","sf")] # Starting point of the chain
pChain <- matrix(0, nrow=chainLength, ncol = no.param*no.var+1) # Initialising the chain


# Defining the variance-covariance matrix for proposal generation
vcov = (0.01*(pMaxima-pMinima))^2
vcovProposal =  vcov[1,] # The higher the coefficient, the lower the acceptance rate with better matching


# Find the Prior probability density
prior.dist = vector("list", no.var)
for (i in 1:no.var) {
  prior.dist[i] = list(log(dunif(pValues[ , i], pMinima[ , i], pMaxima[ , i])))
}
logPrior0 <- sum(unlist(prior.dist))


# Calculating model outputs for the starting point of the chain
Mleaf = Mstem = Mroot = c()
Mleaf[1] <- data$Mleaf[1]
Mstem[1] <- data$Mstem[1]
Mroot[1] <- data$Mroot[1]
output = model(data$GPP,data$Rd,j,Mleaf,Mstem,Mroot,pValues$Y,pValues$k,pValues$af,pValues$as,pValues$sf)


# Calculating the log likelihood of starting point of the chain
logli <- matrix(0, nrow=length(GPP.data$Date), ncol = 1) # Initialising the logli
for (i in 1:length(GPP.data$Date)) {
  if (!is.na(data$Mleaf[i])) {
    logli[i] = - 0.5*((output$Mleaf[i] - data$Mleaf[i])/data$Mleaf_SD[i])^2 - log(data$Mleaf_SD[i])
  }
  if (!is.na(data$Mstem[i])) {
    logli[i] = logli[i] - 0.5*((output$Mstem[i] - data$Mstem[i])/data$Mstem_SD[i])^2 - log(data$Mstem_SD[i])
  }
  if (!is.na(data$Mroot[i])) {
    logli[i] = logli[i] - 0.5*((output$Mroot[i] - data$Mroot[i])/data$Mroot_SD[i])^2 - log(data$Mroot_SD[i])
  }
  if (!is.na(data$Sleaf[i])) {
    logli[i] = logli[i] - 0.5*((output$Sleaf[i] - data$Sleaf[i])/data$Sleaf_SD[i])^2 - log(data$Sleaf_SD[i])
  }
}
logL0 <- sum(logli) # Log likelihood
pChain[1,] <- c(pValues$k,pValues$Y,pValues$af,pValues$as,pValues$sf,logL0) # Assign the first parameter set with log likelihood


# Calculating the next candidate parameter vector, as a multivariate normal jump away from the current point
for (c in (2 : chainLength)) {
  candidatepValues = matrix(ncol = no.var, nrow = no.param)
  for (i in 1:no.var) {
    candidatepValues[ , i] = rmvnorm(n=1, mean=pValues[ , i],
                                     sigma=diag(vcovProposal[i],no.param)) 
  }
  candidatepValues = data.frame(candidatepValues)
  names(candidatepValues) <- c("k","Y","af","as","sf")
  
  
  # Reflected back to generate another candidate value
  reflectionFromMin = pmin( unlist(matrix(0,nrow=no.param,ncol=no.var)), unlist(candidatepValues-pMinima) )
  reflectionFromMax = pmax( unlist(list(rep(0, no.param))), unlist(candidatepValues-pMaxima) )
  candidatepValues = candidatepValues - 2 * reflectionFromMin - 2 * reflectionFromMax 
  
  
  # Calculating the prior probability density for the candidate parameter vector
  if (all(candidatepValues>pMinima) && all(candidatepValues<pMaxima)){
    uni.dist = vector("list", no.var)
    for (i in 1:no.var) {
      uni.dist[i] = list(log(dunif(candidatepValues[ , i], pMinima[ , i], pMaxima[ , i])))
    }
    logPrior1 <- sum(unlist(uni.dist))
    Prior1 = 1
  } else {
    Prior1 <- 0
  }
  
  
  # Calculating the outputs for the candidate parameter vector, log likelihood
  if (Prior1 > 0) {
    Mleaf = Mstem = Mroot = c()
    Mleaf[1] <- data$Mleaf[1]
    Mstem[1] <- data$Mstem[1]
    Mroot[1] <- data$Mroot[1]
    out.cand = model(data$GPP,data$Rd,j,Mleaf,Mstem,Mroot,candidatepValues$Y,
                     candidatepValues$k,candidatepValues$af,candidatepValues$as,candidatepValues$sf)
    
    logli <- matrix(0, nrow=length(GPP.data$Date), ncol = 1) # Initialising the logli
    for (i in 1:length(GPP.data$Date)) {
      if (!is.na(data$Mleaf[i])) {
        logli[i] = - 0.5*((out.cand$Mleaf[i] - data$Mleaf[i])/data$Mleaf_SD[i])^2 - log(data$Mleaf_SD[i])}
      if (!is.na(data$Mstem[i])) {
        logli[i] = logli[i] - 0.5*((out.cand$Mstem[i] - data$Mstem[i])/data$Mstem_SD[i])^2 - log(data$Mstem_SD[i])
      }
      if (!is.na(data$Mroot[i])) {
        logli[i] = logli[i] - 0.5*((out.cand$Mroot[i] - data$Mroot[i])/data$Mroot_SD[i])^2 - log(data$Mroot_SD[i])
      }
      if (!is.na(data$Sleaf[i])) {
        logli[i] = logli[i] - 0.5*((out.cand$Sleaf[i] - data$Sleaf[i])/data$Sleaf_SD[i])^2 - log(data$Sleaf_SD[i])
      }
    }
    logL1 <- sum(logli)
    # Calculating the logarithm of the Metropolis ratio
    logalpha <- (logPrior1+logL1) - (logPrior0+logL0) 
    # Accepting or rejecting the candidate vector
    if ( log(runif(1, min = 0, max =1)) < logalpha ) { 
      pValues <- candidatepValues
      logPrior0 <- logPrior1
      logL0 <- logL1
    }
  }
  pChain[c,] <- c(pValues$k,pValues$Y,pValues$af,pValues$as,pValues$sf,logL0)
}
# Discard the first 500 iterations for Burn-IN in MCMC
pChain <- pChain[501:nrow(pChain),]

# Store the final parameter set values
param.set = colMeans(pChain[ , 1:(no.param*no.var)])
param.SD = apply(pChain[ , 1:(no.param*no.var)], 2, sd)
param.final = data.frame(matrix(ncol = (no.var+1)*2, nrow = no.param))
names(param.final) <- c("k","Y","af","as","ar","sf","k_SD","Y_SD","af_SD","as_SD","ar_SD","sf_SD")
param.final$k = param.set[1:no.param]
param.final$Y = param.set[(1+no.param):(2*no.param)]
param.final$af = param.set[(1+2*no.param):(3*no.param)]
param.final$as = param.set[(1+3*no.param):(4*no.param)]
param.final$sf = param.set[(1+4*no.param):(5*no.param)]
param.final$ar = 1 - param.final$af - param.final$as
param.final$k_SD = param.SD[1:no.param]
param.final$Y_SD = param.SD[(1+no.param):(2*no.param)]
param.final$af_SD = param.SD[(1+2*no.param):(3*no.param)]
param.final$as_SD = param.SD[(1+3*no.param):(4*no.param)]
param.final$sf_SD = param.SD[(1+4*no.param):(5*no.param)]
param.final$ar_SD = with(param.final, (af_SD*af_SD + as_SD*as_SD)^0.5)
# Calculate final output set from the predicted parameter set
Mleaf = Mstem = Mroot = c()
Mleaf[1] <- data$Mleaf[1]
Mstem[1] <- data$Mstem[1]
Mroot[1] <- data$Mroot[1]
output.final = model(data$GPP,data$Rd,j,Mleaf,Mstem,Mroot,param.final$Y,
                     param.final$k,param.final$af,param.final$as,param.final$sf)


# Plotting the Measured (data) vs Modelled Plant Carbon pools for comparison
output.final$Date = data$Date
names(output.final) = c("Cstorage.modelled","Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled","Date")
melted.output = melt(output.final[,c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled","Date")], id.vars="Date")
melted.data = melt(data[ , c("Mleaf","Mstem","Mroot","Sleaf","Date")], id.vars="Date")
melted.data$Date = as.Date(melted.data$Date)
melted.output$Date = as.Date(melted.output$Date)

# png(file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/Measured_vs_modelled_Cpools_monthly.png")
par(mfrow=c(1,1))
ggplot(melted.data, aes(x = Date, y = value, group = variable, colour=factor(variable))) +
  geom_point(shape = 1, size = 1, stroke = 1.25) +
  geom_line(data = melted.output, aes(x = Date, y = value, group = variable, colour=factor(variable))) + 
  ylab("Plant Carbon pool (gC)") +
  theme(axis.text.x=element_text(angle=50, size=10, vjust=0.5)) + 
  ggtitle("Measured (circles) vs Modelled (lines) Plant Carbon pools \n for monthly parameter set") +
  theme(legend.title = element_text(colour="chocolate", size=10, face="bold")) +
  scale_color_discrete(name="C pools") +
  annotate("text", x = melted.output$Date[20], y = max(output$Mstem,na.rm = TRUE), size = 3, 
           label = paste("Mean k = ", round(mean(param.final[,1]), 3), "\nMean Y = ", round(mean(param.final[,2]), 3), 
                         "\nMean af = ", round(mean(param.final[,3]), 3), "\nMean as = ", round(mean(param.final[,4]), 3), 
                         "\nMean ar = ", round(mean(param.final[,7]), 3), "\nMean sf = ",round(mean(param.final[,5]), 3), "\nChain length = ", chainLength))
# dev.off()


# Plotting the parameter sets over time
param.final$Date = data$Date[seq(1,nrow(data),param.vary)]
# param.final$ar = 1 - param.final$af - param.final$as
# melted.param = melt(param.final[,c("k","Y","af","as","ar","sf","Date")], id.vars="Date")
melted.param1 = melt(param.final[,c("k","Y","af","as","ar","sf","Date")], id.vars="Date")
melted.param2 = melt(param.final[,c("k_SD","Y_SD","af_SD","as_SD","ar_SD","sf_SD","Date")], id.vars="Date")
melted.param = data.frame(melted.param1$Date, melted.param1$variable, melted.param1$value, melted.param2$value)
names(melted.param) = c("Date","variable","Parameter","Parameter_SD")
melted.param$Date = as.Date(melted.param$Date)

# png(file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/Allocation_fractions_over_time_monthly.png")
pd <- position_dodge(3) # move the overlapped errorbars horizontally
ggplot(data = melted.param, aes(x = Date, y = Parameter, group = variable, colour=factor(variable))) +
  geom_line(position=pd) +
  geom_errorbar(data = melted.param, aes(ymin=Parameter-Parameter_SD, ymax=Parameter+Parameter_SD), width=5, position=pd) +
  geom_point(position=pd, size=1.5, shape=21, stroke=1.25, fill="white") + # 21 is filled circle
  xlab("Days") +
  ylab("Parameters") +
  ggtitle("Modelled allocation fractions") +
  scale_colour_hue(name="Parameter",    # Legend label, use darker colors
                   l=40) +                    # Use darker colors, lightness=40
  scale_y_continuous(breaks=0:10*0.1)  # Set tick every 0.1
  theme_bw() +
  theme(legend.justification=c(1,1),
        legend.position=c(1.1,1.1)) # Position legend in bottom right
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
se <- function(data,res) { 
  se = SE = 0
  for (i in 1:length(res)) {
    if (!is.na(data[i])) {
      se[i] = (data[i] - res[i])^2 
    }
  }
  N = sum(!is.na(data))
  SE = sum(se,na.rm=TRUE) / N / (mean(data,na.rm=TRUE)*mean(data,na.rm=TRUE))
  return(SE)
}

se1 = se(data$Mleaf,output.final$Mleaf)
t1 = (paste("Normalized Standard error (SE) of Mleaf measurements: ", se1))
se2 = se(data$Mstem,output.final$Mstem)
t2 = (paste("Normalized Standard error (SE) of Mstem measurements: ", se2))
se3 = se(data$Mroot,output.final$Mroot)
t3 = (paste("Normalized Standard error (SE) of Mroot measurements: ", se3))
print(t1)
print(t2)
print(t3)


# Plotting all parameter time series seperately with a moving average for the overall trend
ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
n = 10
# png(file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/Parameters_over_time_monthly.png")
par(mfrow=c(2,3))
plot(param.final$k,type='l',col="red",main="Utilization coefficient, k",xlab="Days")
# lines(ma(param.final$k,n),type='l',col="black")

plot(param.final$Y,type='l',col="chocolate",main="Allocation fraction to Biomass, Y",xlab="Days")
# lines(ma(param.final$Y,n),type='l',col="black")

plot(param.final$af,type='l',col="green",main="Allocation fraction to foliage, af",xlab="Days")
# lines(ma(param.final$af,n),type='l',col="black")

plot(param.final$as,type='l',col="blue",main="Allocation fraction to stem, as",xlab="Days")
# lines(ma(param.final$as,n),type='l',col="black")

plot(param.final$ar,type='l',col="magenta",main="Allocation fraction to root, ar",xlab="Days")
# lines(ma(param.final$ar,n),type='l',col="black")

plot(param.final$sf,type='l',col="purple",main="Foliage tunrover rate, sf",xlab="Days")
# lines(ma(param.final$sf,n),type='l',col="black")
# dev.off()


# Plotting all parameter whole iterations for Day 1 only to check the convergance
# png(file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/Parameter_iterations_monthly.png")
par(mfrow=c(2,3))
plot(pChain[,1],col="red",main="Utilization coefficient at Day 1",xlab="Iterations",ylab="k")
plot(pChain[,1+no.param],col="green",main="Alloc frac to Biomass at Day 1",xlab="Iterations",ylab="Y")
plot(pChain[,1+2*no.param],col="magenta",main="Alloc frac to foliage at Day 1",xlab="Iterations",ylab="af")
plot(pChain[,1+3*no.param],col="blue",main="Alloc frac to stem at Day 1",xlab="Iterations",ylab="as")
plot(pChain[,1+4*no.param],col="green",main="Foliage turnover at Day 1",xlab="Iterations",ylab="sf")
plot(pChain[,1+5*no.param],col="magenta",main="Log-likelihood at Day 1",xlab="Iterations",ylab="Log-likelihood")
# dev.off()


# Display the final mean parameter values
colMeans(param.final[ , c(1:6)])

# Calcualte AIC and BIC to find the most accurate model for best balance between model fit and complexity
k = 2 # k = 2 for the usual AIC
npar = no.param*no.var # npar = total number of parameters in the fitted model
aic = -2*sum(pChain[,ncol(pChain)]) + k*npar

n = sum(!is.na(data$Sleaf)) + sum(!is.na(data$Mleaf)) + sum(!is.na(data$Mstem)) + sum(!is.na(data$Mroot))
k = log(n) # n being the number of observations for the so-called BIC
bic = -2*sum(pChain[,ncol(pChain)]) + k*npar


