# Carbon balance model 
# Developed by Kashif Mahmud and Belinda Medlyn (November 2016)
# k.mahmud@westernsydney.edu.au

# This code carries out Bayesian calibration for 4 variables (allocation fractions: "k","af","as","sf") on 
# daily time scale (e.g. 120 days) to estimate Carbon pools (Cstorage,Cleaf,Cstem,Croot)

##############################
# Version = v15: MCMC with soil manipulation pot experiment data for all treatments (including the free seedling), 
# This version considers either daily/weekly/monthly/just one parameter set for 5 variables ("k","Y","af","as","sf")
# So we can set the parameters for various time frames
# Also calculates the MCMC SDs for all parameters with different time frames, and also the LogLi, AIC, BIC, time (measures 
# for best model selection) to select the best parameter set
# Finally save the figures in Github/results folder
##############################

rm(list=ls())
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM/Data_files")


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


# install.packages("mvtnorm")
library(mvtnorm) # Creates candidate parameter vector as a multivariate normal jump away from the current candidate
library(reshape2)
library(ggplot2)
chainLength = 1500 # Setting the length of the Markov Chain to be generated
no.var = 5 # variables to be modelled are: k,Y,af,as,sf
vol = c(15,35) # test run
no.param.par.var = c(1,3) # test run
# GPP.data.raw = read.csv("GPP.csv") # Units gC d-1
# vol = unique(GPP.data.raw$volume)[order(unique(GPP.data.raw$volume))]
# no.param.par.var = c(1,2,3,4,5,6,9) # temporal parameter count per variable
param.mean = data.frame(matrix(ncol = no.var+1, nrow = length(no.param.par.var)*length(vol)))
names(param.mean) = c("k","Y","af","as","ar","sf")
aic.bic = data.frame(matrix(ncol = 4, nrow = length(no.param.par.var)*length(vol)))
time = data.frame(no.param=rep(no.param.par.var,length(vol)),
                  start.time=numeric(length(no.param.par.var)*length(vol)),
                  end.time=numeric(length(no.param.par.var)*length(vol)),
                  time.taken=numeric(length(no.param.par.var)*length(vol)))
q = 0 # Indicates the iteration number


# Import daily GPP, daily Rd
GPP.data.raw = read.csv("GPP.csv") # Units gC d-1
Rd.data.raw = read.csv("Rd.csv") # Units g C g-1 plant d-1
tnc.data.raw = read.csv("tnc_fortnightly_data.csv") # Units g plant-1

# Import weekly Cleaf, weekly Cstem, initial/harvest Croot data with Mean and SD
Mleaf.data.raw = read.csv("Cleaf_weekly_data.csv") # Units gC d-1
Mstem.data.raw = read.csv("Cstem_weekly_data.csv") # Units gC d-1
Mroot.data.raw = read.csv("Croot_twice_data.csv") # Units gC d-1


for (z in 1:length(no.param.par.var)) {
  for (v in 1:length(vol)) {
    GPP.data = subset(GPP.data.raw,volume==vol[v]) # Consider only free seedling to start with
    names(GPP.data)[3] = "GPP"
    Rd.data = subset(Rd.data.raw,volume==vol[v])
    Sleaf.data = tnc.data = subset(tnc.data.raw,volume==vol[v])
    Mleaf.data = subset(Mleaf.data.raw,volume==vol[v])
    Mstem.data = subset(Mstem.data.raw,volume==vol[v])
    Mroot.data = subset(Mroot.data.raw,volume==vol[v])
    
    
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
    
    
    # Initialize few output data files
    q = q + 1
    time$start.time[q] <- Sys.time()
    param.vary = ceiling(nrow(data)/no.param.par.var[z]) # How many days the parameter set remain unchanged (weekly = 7; monthly = 30; just one parameter = nrow(data))
    no.param = ceiling(nrow(data)/param.vary) # number of parameter set for the whole duration of experiment (121 days)
    j = c()
    j[1] = 0
    i = seq(1,nrow(data),1)
    j[i] = i - ceiling(i/param.vary)*1  # j is for parameter settings for various time frames
    
    
    # Setting lower and upper bounds of the prior parameter pdf, and starting point of the chain
    param.k <- matrix(c(0,0.45,1) , nrow=no.param, ncol=3, byrow=T) 
    param.Y <- matrix(c(0.1,0.3,0.5) , nrow=no.param, ncol=3, byrow=T) 
    param.af <- matrix(c(0,0.45,0.7) , nrow=no.param, ncol=3, byrow=T) 
    param.as <- matrix(c(0,0.2,0.5) , nrow=no.param, ncol=3, byrow=T) 
    param.sf <- matrix(c(0,0.05,0.1) , nrow=no.param, ncol=3, byrow=T) 
    
    # param.k <- matrix(c(0,0.45,1) , nrow=no.param, ncol=3, byrow=T) 
    # param.Y <- matrix(c(0.2,0.3,0.4) , nrow=no.param, ncol=3, byrow=T) 
    # param.af <- matrix(c(0,0.45,0.7) , nrow=no.param, ncol=3, byrow=T) 
    # param.as <- matrix(c(0,0.17,0.5) , nrow=no.param, ncol=3, byrow=T) 
    # param.sf <- matrix(c(0,0.02,0.04) , nrow=no.param, ncol=3, byrow=T) 
    
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
    
    
    # Plotting the Measured (data) vs Modelled Plant Carbon pools for plotting and comparison
    output.final$Date = data$Date
    names(output.final) = c("Cstorage.modelled","Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled","Date")
    melted.output = melt(output.final[,c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled","Date")], id.vars="Date")
    melted.Cstorage = output.final[,c("Cstorage.modelled","Date")]
    melted.data = melt(data[ , c("Mleaf","Mstem","Mroot","Sleaf","Date")], id.vars="Date")
    melted.data$Date = as.Date(melted.data$Date)
    melted.error = melt(data[ , c("Mleaf_SD","Mstem_SD","Mroot_SD","Sleaf_SD","Date")], id.vars="Date")
    melted.error$Date = as.Date(melted.error$Date)
    melted.error$volume = as.factor(vol[v])
    melted.error$parameter = melted.data$value
    melted.output$Date = as.Date(melted.output$Date)
    melted.data$volume = as.factor(vol[v])
    melted.output$volume = as.factor(vol[v])
    melted.output$no.param = as.factor(no.param.par.var[z])
    melted.Cstorage$Date = as.Date(melted.Cstorage$Date)
    melted.Cstorage$volume = as.factor(vol[v])
    melted.Cstorage$no.param = as.factor(no.param.par.var[z])
    
    
    # Plotting the parameter sets over time
    param.final$Date = data$Date[seq(1,nrow(data),param.vary)]
    melted.param1 = melt(param.final[,c("k","Y","af","as","ar","sf","Date")], id.vars="Date")
    melted.param2 = melt(param.final[,c("k_SD","Y_SD","af_SD","as_SD","ar_SD","sf_SD","Date")], id.vars="Date")
    melted.param = data.frame(melted.param1$Date, melted.param1$variable, melted.param1$value, melted.param2$value)
    names(melted.param) = c("Date","variable","Parameter","Parameter_SD")
    melted.param$Date = as.Date(melted.param$Date)
    melted.param$volume = vol[v]
    melted.param$no.param = as.factor(no.param.par.var[z])
    
    # Set working directory for saving figures
    setwd("/Users/kashifmahmud/WSU/ARC_project/CBM/Results")
    
    # Plotting C pools over time for individual volume and No. of parameter
    pd <- position_dodge(3) # move the overlapped errorbars horizontally
    # p1 = ggplot(melted.data, aes(x = Date, y = value, group = variable, colour=factor(variable))) +
    #   geom_point(shape = 1, size = 1, stroke = 1.25) +
    #   geom_line(data = melted.output, aes(x = Date, y = value, group = variable, colour=factor(variable))) + 
    #   ylab("Plant Carbon pool (gC)") +
    #   theme(axis.text.x=element_text(angle=50, size=10, vjust=0.5)) + 
    #   ggtitle("Measured (circles) vs Modelled (lines) Plant Carbon pools") +
    #   theme(legend.title = element_text(colour="chocolate", size=10, face="bold")) +
    #   scale_color_discrete(name="C pools") +
    #   annotate("text", x = melted.output$Date[20], y = max(output$Mstem,na.rm = TRUE), size = 3, 
    #            label = paste("Mean k = ", round(mean(param.final[,1]), 3), "\nMean Y = ", round(mean(param.final[,2]), 3), 
    #                          "\nMean af = ", round(mean(param.final[,3]), 3), "\nMean as = ", round(mean(param.final[,4]), 3), 
    #                          "\nMean ar = ", round(mean(param.final[,7]), 3), "\nMean sf = ",round(mean(param.final[,5]), 3), "\nChain length = ", chainLength))
    # p1
    p1 = ggplot(melted.error, aes(x=Date, y=parameter, colour=variable, group=variable)) + 
      geom_errorbar(aes(ymin=parameter-value, ymax=parameter+value,colour=variable, group=variable), width=5) +
      geom_line(data = melted.output, aes(x = Date, y = value, group = variable, colour=variable)) + 
      geom_point(shape = 1, size = 1, stroke = 1.25) +
      ylab("Plant Carbon pool (gC)") +
      ggtitle("Measured (circles) vs Modelled (lines) Plant Carbon pools") +
      labs(colour="C pools") +
      theme(legend.title = element_text(colour="chocolate", size=10, face="bold")) +
      annotate("text", x = melted.output$Date[20], y = max(output$Mstem,na.rm = TRUE), size = 3, 
                          label = paste("Mean k = ", round(mean(param.final[,1]), 3), "\nMean Y = ", round(mean(param.final[,2]), 3),
                                        "\nMean af = ", round(mean(param.final[,3]), 3), "\nMean as = ", round(mean(param.final[,4]), 3),
                                        "\nMean ar = ", round(mean(param.final[,7]), 3), "\nMean sf = ",round(mean(param.final[,5]), 3), "\nChain length = ", chainLength-500))
    p1
    ggsave(p1,filename=paste("Measured_vs_Modelled_Carbon_pools_vol_",vol[v],"_par_",no.param.par.var[z],".png",sep=""))
    
    
    # Plotting Allocation fractions over time for individual volume and No. of parameter
    pd <- position_dodge(3) # move the overlapped errorbars horizontally
    p2 = ggplot(data = melted.param, aes(x = Date, y = Parameter, group = variable, colour=factor(variable))) +
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
    p2
    ggsave(p2,filename=paste("Allocation_fractions_over_time_vol_",vol[v],"_par_",no.param.par.var[z],".png",sep=""))
    
    # Storing the summary of data, outputs, Cstorage, parameters
    if (q == 1) {
      summary.data = melted.data
      summary.error = melted.error
      summary.output = melted.output
      summary.Cstorage = melted.Cstorage
      summary.param = melted.param
    }
    if (q > 1) {
      summary.output = rbind(summary.output,melted.output)
      summary.Cstorage = rbind(summary.Cstorage,melted.Cstorage)
      summary.param = rbind(summary.param,melted.param)
      if (z == 1) {
        summary.data = rbind(summary.data,melted.data)
        summary.error = rbind(summary.error,melted.error)
      }
    }
    
    
    # Display the Acceptance rate of the chain
    nAccepted = length(unique(pChain[,1]))
    acceptance = (paste("Volume =",vol[v],", Total Parameter number =",no.param.par.var[z],": ", nAccepted, "out of ", chainLength, "candidates accepted ( = ",
                        round(100*nAccepted/chainLength), "%)"))
    print(acceptance)
    
    
    # Plotting all parameter time series seperately with a moving average for the overall trend
    ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
    n = 10
    png(file = paste("Parameters_vol_",vol[v],"_par_",no.param.par.var[z], ".png", sep = ""))
    par(mfrow=c(2,3))
    plot(param.final$k,type='p',col="red",main="Utilization coefficient, k",xlab="Days")
    # lines(ma(param.final$k,n),type='l',col="black")
    plot(param.final$Y,type='p',col="chocolate",main="Allocation fraction to Biomass, Y",xlab="Days")
    # lines(ma(param.final$Y,n),type='l',col="black")
    plot(param.final$af,type='p',col="green",main="Allocation fraction to foliage, af",xlab="Days")
    # lines(ma(param.final$af,n),type='l',col="black")
    plot(param.final$as,type='p',col="blue",main="Allocation fraction to stem, as",xlab="Days")
    # lines(ma(param.final$as,n),type='l',col="black")
    plot(param.final$ar,type='p',col="magenta",main="Allocation fraction to root, ar",xlab="Days")
    # lines(ma(param.final$ar,n),type='l',col="black")
    plot(param.final$sf,type='p',col="purple",main="Foliage tunrover rate, sf",xlab="Days")
    # lines(ma(param.final$sf,n),type='l',col="black")
    dev.off()
    
    
    # Plotting all parameter whole iterations for Day 1 only to check the convergance
    png(file = paste("Parameter_iterations_day1_vol_",vol[v],"_par_",no.param.par.var[z], ".png", sep = ""))
    par(mfrow=c(2,3))
    plot(pChain[,1],col="red",main="Utilization coefficient at Day 1",xlab="Iterations",ylab="k")
    plot(pChain[,1+no.param],col="green",main="Alloc frac to Biomass at Day 1",xlab="Iterations",ylab="Y")
    plot(pChain[,1+2*no.param],col="magenta",main="Alloc frac to foliage at Day 1",xlab="Iterations",ylab="af")
    plot(pChain[,1+3*no.param],col="blue",main="Alloc frac to stem at Day 1",xlab="Iterations",ylab="as")
    plot(pChain[,1+4*no.param],col="green",main="Foliage turnover at Day 1",xlab="Iterations",ylab="sf")
    plot(pChain[,1+5*no.param],col="magenta",main="Log-likelihood",xlab="Iterations",ylab="Log-likelihood")
    dev.off()
    
    
    # Store the final mean parameter values
    param.mean[q,c(1:6)] = colMeans(param.final[ , c(1:6)])
    param.mean$volume[q] = vol[v]
    param.mean$no.param[q] = no.param.par.var[z]
    
    
    # Calcualte LogLi, AIC, BIC, Time to find the most accurate model for best balance between model fit and complexity
    logLi <- matrix(0, nrow=length(GPP.data$Date), ncol = 1) # Initialising the logLi
    for (i in 1:length(GPP.data$Date)) {
      if (!is.na(data$Mleaf[i])) {
        logLi[i] = - 0.5*((output.final$Mleaf.modelled[i] - data$Mleaf[i])/data$Mleaf_SD[i])^2 - log(data$Mleaf_SD[i])
      }
      if (!is.na(data$Mstem[i])) {
        logLi[i] = logLi[i] - 0.5*((output.final$Mstem.modelled[i] - data$Mstem[i])/data$Mstem_SD[i])^2 - log(data$Mstem_SD[i])
      }
      if (!is.na(data$Mroot[i])) {
        logLi[i] = logLi[i] - 0.5*((output.final$Mroot.modelled[i] - data$Mroot[i])/data$Mroot_SD[i])^2 - log(data$Mroot_SD[i])
      }
      if (!is.na(data$Sleaf[i])) {
        logLi[i] = logLi[i] - 0.5*((output.final$Sleaf.modelled[i] - data$Sleaf[i])/data$Sleaf_SD[i])^2 - log(data$Sleaf_SD[i])
      }
    }
    aic.bic[q,1] <- sum(logLi)
    
    k1 = 2 # k = 2 for the usual AIC
    npar = no.param*no.var # npar = total number of parameters in the fitted model
    aic.bic[q,2] = -2*aic.bic[q,1] + k1*npar
    
    n = sum(!is.na(data$Sleaf)) + sum(!is.na(data$Mleaf)) + sum(!is.na(data$Mstem)) + sum(!is.na(data$Mroot))
    k2 = log(n) # n being the number of observations for the so-called BIC
    aic.bic[q,3] = -2*aic.bic[q,1] + k2*npar
    
    time$end.time[q] <- Sys.time()
    time$time.taken[q] <- time$end.time[q] - time$start.time[q]
    aic.bic[q,4] = time$time.taken[q]
    aic.bic[q,5] = vol[v]
    aic.bic[q,6] = no.param.par.var[z]
    }
}


# Set working directory for saving figures
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM/Results")


# Plot modelled vs measured data ("Mleaf","Mstem","Mroot","Sleaf") against "volume" and "Total No of param"
meas = as.factor(c("Mleaf","Mstem","Mroot","Sleaf"))
res = as.factor(c("Mleaf.modelled","Mstem.modelled","Mroot.modelled","Sleaf.modelled"))
error = as.factor(c("Mleaf_SD","Mstem_SD","Mroot_SD","Sleaf_SD"))
pd <- position_dodge(3) # move the overlapped errorbars horizontally
for (p in 1:length(meas)) {
  summary.data.Cpool = subset(summary.data,variable==meas[p])
  summary.output.Cpool = subset(summary.output,variable==res[p])
  summary.error.Cpool = subset(summary.error,variable==error[p])
  summary.error.Cpool$parameter = summary.data.Cpool$value
  # p3 = ggplot() +
  #   geom_point(position=pd,data=summary.error.Cpool, aes(x = Date, y = parameter, ymin=parameter-value, ymax=parameter+value, group = volume, shape=factor(volume))) +
  #   geom_errorbar(data=summary.error.Cpool, aes(x = Date, y = parameter, ymin=parameter-value, ymax=parameter+value),colour="grey",width=2,position=pd) +
  #   geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,no.param), shape=factor(volume), colour=factor(no.param))) + 
  #   ylab(as.character(meas[p])) +
  #   ggtitle("C pools - Measured (points) vs Modelled (lines)") +
  #   labs(shape="Soil Volume", colour="Total No of Parameter") +
  #   theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
  # p3
  p3 = ggplot(summary.error.Cpool, aes(x=Date, y=parameter, shape=volume, group=volume)) + 
    geom_errorbar(aes(ymin=parameter-value, ymax=parameter+value), colour="grey", width=3, position=pd) +
    geom_line(position=pd,data = summary.output.Cpool, aes(x = Date, y = value, group = interaction(volume,no.param), linetype=volume, colour=no.param)) + 
    geom_point(position=pd, size=2) +
    ylab(as.character(meas[p])) +
    ggtitle("C pools - Measured (points) vs Modelled (lines)") +
    labs(shape="Soil Volume", linetype="Soil Volume", colour="Total No of Parameter") +
    theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
  p3
  ggsave(p3,filename=paste(meas[p],"_Measured_vs_Modelled.png",sep=""))
}


# Plot modelled Cstorage against "volume" and "Total No of param"
p4 = ggplot() +
  geom_line(data = summary.Cstorage, aes(x = Date, y = Cstorage.modelled, group = interaction(volume,no.param),colour=factor(no.param), linetype=factor(volume))) + 
  ylab("Cstorage (gC)") +
  ggtitle("Modelled Cstorage") +
  labs(linetype="Soil Volume", colour="Total No of Parameter") +
  theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
p4
ggsave(p4,filename=paste("Cstorage_Modelled.png",sep=""))


# Plot individual modelled parameters ("k","Y","af","as","ar","sf") against "volume" and "Total No of param"
var = as.factor(c("k","Y","af","as","ar","sf"))
for (p in 1:length(var)) {
  summary.param.set = subset(summary.param,variable==var[p])
  pd <- position_dodge(0.5) # move the overlapped errorbars horizontally
  p5 = ggplot() +
    geom_point(position=pd,data = summary.param.set, aes(x = Date, y = Parameter,  group = interaction(volume,no.param), colour=factor(no.param), shape=factor(volume))) +
    geom_line(position=pd,data = summary.param.set, aes(x = Date, y = Parameter,  group = interaction(volume,no.param), colour=factor(no.param), linetype=factor(volume))) +
    xlab("Days") +
    ylab(as.character(var[p])) +
    ggtitle("Modelled coefficients") +
    labs(shape="Soil Volume", linetype="Soil Volume", colour="Total No of Parameter") +
    theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
  ggsave(p5,filename=paste(var[p],"_over_time.png",sep=""))
  p5
}


# Plot modelled parameter means ("k","Y","af","as","ar","sf") against "volume" and "Total No of param"
melted.param.mean = melt(param.mean, id.vars=c("no.param","volume"))
pd <- position_dodge(0.1)
p6 = ggplot() +
  geom_point(position=pd, data = melted.param.mean, aes(x = variable, y = value, group = interaction(volume,no.param), colour=factor(no.param), shape=factor(volume))) +
  xlab("Allocation fractions") +
  ylab("Value of the coefficients") +
  ggtitle("Modelled allocation fractions") +
  labs(shape="Soil Volume", colour="Total No of Parameter") +
  theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
p6
ggsave(p6,filename=paste("Modelled_mean_allocation_fractions.png"))


# Plot Model Measures ("logLi","aic","bic","time") against "volume" and "Total No of param"
names(aic.bic) <- c("logLi","aic","bic","time","volume","no.param")
write.csv(aic.bic, file = "/Users/kashifmahmud/WSU/ARC_project/CBM/Results/logli_aic_bic_time.csv", row.names = FALSE)
melted.aic.bic = melt(aic.bic, id.vars=c("no.param","volume"))
pd <- position_dodge(0.1)
p7 = ggplot(data = melted.aic.bic, aes(x = variable, y = value, group = interaction(volume,no.param), colour=factor(no.param), shape=factor(volume))) +
  geom_point(position=pd) +
  xlab("Model Measures") +
  ylab("LogLi, AIC, BIC, Time") +
  ggtitle("LogLi, AIC, BIC, Time for various models") +
  labs(shape="Soil Volume", colour="Total No of Parameter") +
  theme(legend.title = element_text(colour="chocolate", size=10, face="bold"))
p7
ggsave(p7,filename=paste("LogLi_aic_bic_time.png"))




