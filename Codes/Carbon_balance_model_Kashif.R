# Carbon Balance Model 
# R code to import and process data, and then model the carbon pools and fluxes
rm(list=ls())
setwd("/Users/kashifmahmud/WSU/ARC_project/Data_files/Court_data")

# Calculate Cleaf pool
# inputs:
number_of_pots = 7; # number of soil manipulation tests with volumes of 5, 10, 15, 20, 25, 35, Free(1000) Litres
number_of_plots = 7; # number of replicas

# Court's data files - Soil Manipulation test
# plot_summary.xlsx - volume designations by plot and pot
plot.summary = read.csv("plot_summary.csv")
head(plot.summary)
# Find different soil volume pots from the data set
vols = unique(plot.summary$volume, fromLast = TRUE)
vols = vols[order(vols)]

# leaf_data.xlsx - Leaf mass and SLA data (fortnightly measurements)
leaf.data = read.csv("leaf_data.csv")
leaf.data$Date = as.Date(leaf.data$Date, format = "%d/%m/%Y")
head(leaf.data)
class(leaf.data$Date)
# leaf_count.xlsx - Leaf count data (weekly measurements)
leaf.count = read.csv("leaf_count_1.csv")
names(leaf.count)[1:2] <- c("plot", "pot")
#leaf.count = leaf.count[-2, ]
# lc.date = leaf.count[1,3:ncol(leaf.count)]
# lc.date <- parse_date_time(lc.date,"d m y")
lc <- merge(plot.summary,leaf.count,by=c("pot","plot"))
lc.read = read.csv("leaf_count_1.csv", header = F)
lc.final <- as.data.frame(t(lc.read[1,3:ncol(leaf.count)]))
names(lc.final)[1] <- "Date"
lc.final$Date = as.Date(lc.final$Date, format = "%m/%d/%Y")
# library(lubridate)
# lc.final$Date <- parse_date_time(lc.final$D,"m d y")
# keeps <- "Date"
# lc.final = lc.final[ , keeps, drop = FALSE]
dimnames(lc.final)[[1]] <- 1:nrow(lc.final)

################### Count total leaf numbers for various soil manipulation tests (Weekly data)
for(i in 1:length(vols)) { #-- Create objects  'lc.idn.1', 'lc.idn.2', ... 'lc.idn.7' --
  lc.idn = subset(lc,volume==vols[i]) 
  lc.idn[nrow(lc.idn)+1, ] = colMeans(lc.idn, na.rm = TRUE) # R8 = Average of leaf counts
  lc.idn[] <- round(lc.idn)
  nam <- paste("volume_", vols[i], sep = "")
  lc.final[i+1] <- as.data.frame(t(lc.idn[nrow(lc.idn), 4:ncol(lc.idn)]))
  # lc.final[2,1:ncol(lc.final)] = lc.idn[nrow(lc.idn), 4:ncol(lc.idn)]
  names(lc.final)[i+1] <- nam
}

################### Interpolate leaf counts for daily data from weekly measurements
library(reshape2)
library(ggplot2)
lc.melt <- melt(lc.final, id.vars = "Date")
names(lc.melt)[2:3] <- c("volume", "Leaf_Count")

# Perform Cubic Spline interpolation to get daily leaf count
lc.daily = data.frame(seq(lc.melt$Date[1], lc.melt$Date[nrow(lc.melt)], 1))
names(lc.daily)[1] <- "Date"
number_of_day = nrow(lc.daily)
# lc.daily$Date = as.Date(lc.daily$Date, format = "%d/%m/%Y")
for(i in 1:length(vols)) {
  func = splinefun(x=lc.final$Date, y=lc.final[ ,i+1], method="fmm",  ties = mean)
  lc.daily[ ,i+1] = round(func(lc.daily$Date))
}
names(lc.daily)[2:ncol(lc.daily)] <- c("volume_5", "volume_10", "volume_15", "volume_20", "volume_25", "volume_35", "volume_1000")
lc.daily.melt <- melt(lc.daily, id.vars = "Date")
names(lc.daily.melt)[2:3] <- c("volume", "Leaf_Count")
# plot interpolated daily leaf counts
png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/Leaf_count_daily.png")
ggplot() +
  geom_line(data = lc.daily.melt, aes(x = Date, y = Leaf_Count, group = volume, colour=factor(volume))) + 
  xlab("Date") +
  ylab("Daily leaf count") +
  ggtitle("Daily leaf count")
dev.off()


################### Harvested leaf mass and leaf area
harvest_data = read.csv("/Users/kashifmahmud/WSU/ARC_project/Data_files/Court_data/harvest aboveground mass.csv")
harvest_data <- merge(plot.summary,harvest_data,by=c("pot","plot"))
for(i in 1:length(vols)) {
  hd.idn = subset(harvest_data,volume==vols[i]) 
  hd.idn[nrow(hd.idn)+1, ] = colMeans(hd.idn, na.rm = TRUE) # R8 = Average of leaf counts
  hd.idn[nrow(hd.idn)+1, ] = apply(hd.idn, 2, sd) # R9 = Standard deviation of leaf counts
  hd.idn$newleaf_count <- round(hd.idn$newleaf_count)
  hd.idn$leaf_count <- round(hd.idn$leaf_count)
  hd.idn$leaf_mass <- hd.idn$Leafmass.bag - hd.idn$Leaf_bag
  hd.idn$stem_mass <- hd.idn$stemmass.bag - hd.idn$stem_bag
  keeps <- c("volume", "leaf_area", "leaf_count", "leaf_mass", "stem_mass", "Croot", "Froot")
  hd = hd.idn[ , keeps, drop = FALSE]
  dimnames(hd)[[1]] <- c(1:7, "mean", "SD")
  if (i == 1) {
    hd.final <- hd[0,]
    }
  hd.final[i, ] <- hd["mean", ]
}
dimnames(hd.final)[[1]] <- c(1:7)

################### Comparison between leaf area directly from harvest and considering the SLA changes
# Calculate leaf area directly from harvest data
# Leaf area (t) = Leaf area (T) * Leaf count (t) / Leaf count (T); t = time, T = time of harvest
la.daily = data.frame(lc.daily$Date)
for(i in 1:length(vols)) {
  la.daily[ , i+1] = lc.daily[ , i+1] * hd.final$leaf_area[i] / hd.final$leaf_count[i]
}
names(la.daily)[1:ncol(la.daily)] <- c("Date", "volume_5", "volume_10", "volume_15", "volume_20", "volume_25", "volume_35", "volume_1000")

la.daily.melt <- melt(la.daily, id.vars = "Date")
names(la.daily.melt)[2:3] <- c("volume", "Leaf_area")
# plot interpolated daily leaf area from harvest data
png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/Leaf_area_daily.png")
ggplot() +
  geom_line(data = la.daily.melt, aes(x = Date, y = Leaf_area, group = volume, colour=factor(volume))) +
  ylab("Leaf area (cm^2)") +
  xlab("Days") +
  ggtitle("Leaf area from harvest data")
dev.off()

# Calculate leaf area considering the SLA changes
hd.final$sla_harvest = hd.final$leaf_area / hd.final$leaf_mass / 100^2 
keeps <- c("Date", "volume", "sla", "starch_mgperg", "sugars_mgperg")
sla.gas = leaf.data[ , keeps, drop = FALSE]
names(sla.gas)[1:2] <- c("Date", "volume")

for(i in 1:length(vols)) {
  sla.idn = subset(sla.gas,volume==vols[i]) 
  for(j in 1:length(unique(sla.idn$Date))) {
    sla.idn.date = subset(sla.idn, Date == unique(sla.idn$Date)[j])
    sla.idn.date[nrow(sla.idn.date)+1, 2:ncol(sla.idn.date)] = colMeans(sla.idn.date[2:ncol(sla.idn.date)], na.rm = TRUE) # R7 = Average of leaf data
    sla.idn.date[nrow(sla.idn.date)+1, 2:ncol(sla.idn.date)] = apply(sla.idn.date[2:ncol(sla.idn.date)], 2, sd) # R8 = Standard deviation of leaf counts
    sla.idn.date$Date = sla.idn.date[1,1]
    dimnames(sla.idn.date)[[1]] <- c(1:(nrow(sla.idn.date)-2), "Mean", "SD")
    if (i == 1 && j == 1) {
      sla.final <- sla.idn.date[0,]
    }
    sla.final[j+(i-1)*length(unique(sla.idn$Date)), ] <- sla.idn.date["Mean", ]
    sla.final$sla_SD[j+(i-1)*length(unique(sla.idn$Date))] <- sla.idn.date["SD", 3]
    sla.final$starch_SD[j+(i-1)*length(unique(sla.idn$Date))] <- sla.idn.date["SD", 4]
    sla.final$sugars_SD[j+(i-1)*length(unique(sla.idn$Date))] <- sla.idn.date["SD", 5]
  }
}

# Assign numeric volumes in the "volume" variable column
lc.daily.m = lc.daily.melt
lc.daily.m$volume = as.numeric(lc.daily.m$volume)
for(i in length(vols):1) {
  ind = which(lc.daily.m$volume %in% i)
  lc.daily.m$volume[ind] = vols[i]
}
# Take the rows with GEC Date
lc.daily.m = lc.daily.m[lc.daily.m$Date %in% sla.gas$Date, ]

# Merge both data frame to get the leaf counts for GEC Date
leaf.data.gas = merge(sla.final,lc.daily.m,by=c("Date","volume"))
leaf.data.gas = leaf.data.gas[with(leaf.data.gas, order(volume)), ]

# Modify sla values considering changes of sla due to new leaves only
for(i in 1:length(vols)) {
  leaf.data.gas$sla_corrected[1+(i-1)*length(unique(leaf.data.gas$Date))] = leaf.data.gas$sla[1+(i-1)*length(unique(leaf.data.gas$Date))]
  for(j in 2:length(unique(leaf.data.gas$Date))) {
    if (leaf.data.gas$Leaf_Count[j+(i-1)*length(unique(leaf.data.gas$Date))] > leaf.data.gas$Leaf_Count[j-1+(i-1)*length(unique(leaf.data.gas$Date))]) {
      leaf.data.gas$sla_corrected[j+(i-1)*length(unique(leaf.data.gas$Date))] = (leaf.data.gas$sla_corrected[j-1+(i-1)*length(unique(leaf.data.gas$Date))] * 
                                                                                    leaf.data.gas$Leaf_Count[j-1+(i-1)*length(unique(leaf.data.gas$Date))] + ((leaf.data.gas$Leaf_Count[j+(i-1)*length(unique(leaf.data.gas$Date))] - 
                                                                                                                                                                 leaf.data.gas$Leaf_Count[j-1+(i-1)*length(unique(leaf.data.gas$Date))]) * leaf.data.gas$sla[j+(i-1)*length(unique(leaf.data.gas$Date))]) ) / 
        leaf.data.gas$Leaf_Count[j+(i-1)*length(unique(leaf.data.gas$Date))]
    } else {
      leaf.data.gas$sla_corrected[j+(i-1)*length(unique(leaf.data.gas$Date))] = leaf.data.gas$sla[j+(i-1)*length(unique(leaf.data.gas$Date))]
    }
  }
}
leaf.data.gas = leaf.data.gas[with(leaf.data.gas, order(volume)), ]
# Plot the differences in SLA results (Measured SLA and Modified SLA)
ggplot() +
  geom_line(data = leaf.data.gas, aes(x = Date, y = sla, group = volume, colour=factor(volume))) + 
  xlab("Date") +
  ylab("SLA") +
  ggtitle("Measured SLA")
  
ggplot() +
   geom_line(data = leaf.data.gas, aes(x = Date, y = sla_corrected, group = volume, colour=factor(volume))) + 
  xlab("Date") +
  ylab("SLA") +
  ggtitle("Modified SLA")
  
# par(mfrow=c(1,2))
with(subset(leaf.data.gas,volume==35),plot(Date,sla))
with(subset(leaf.data.gas,volume==35),points(Date,sla_corrected,col="red"))

# png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/sla_gas.png")
# ggplot() +
#   geom_line(data = leaf.data.gas, aes(x = Date, y = sla, group = volume, colour=factor(volume)))
# dev.off()
# png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/sla_gas_corrected.png")
# ggplot() +
#   geom_line(data = leaf.data.gas, aes(x = Date, y = sla_corrected, group = volume, colour=factor(volume)))
# dev.off()


# Leaf area calcualation based on modified SLA
######## NO NEED AS PER THE RESULTS AND DISCUSSION WITH BELINDA #############
######## THE SLA CHANGES ONLY REPRESENT THE UPPER CANOPY LEAVES (NEW ONES) #############


################### Calculate leaf mass (leaf carbon pool) directly from harvest data
# Leaf mass (t) = Leaf mass (T) * Leaf count (t) / Leaf count (T); t = time, T = time of harvest
lm.daily = data.frame(lc.daily$Date)
for(i in 1:length(vols)) {
  lm.daily[ , i+1] = lc.daily[ , i+1] * hd.final$leaf_mass[i] / hd.final$leaf_count[i]
  lm.daily[ , i+1] = lm.daily[ , i+1] * 0.5 # Convert gm of DM to gm of C (multiplying by 0.5)
}
names(lm.daily)[1:ncol(lm.daily)] <- c("Date", "volume_5", "volume_10", "volume_15", "volume_20", "volume_25", "volume_35", "volume_1000")

lm.daily.melt <- melt(lm.daily, id.vars = "Date")
names(lm.daily.melt)[2:3] <- c("volume", "leafmass")
# plot interpolated daily leaf mass from harvest data
png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/Leaf_mass_daily.png")
ggplot() +
  geom_line(data = lm.daily.melt, aes(x = Date, y = leafmass, group = volume, colour=factor(volume))) +
  ylab("Leaf mass (gmC)") +
  xlab("Days") +
  ggtitle("Leaf mass from harvest data")
dev.off()


################### Analyse stem height diameter to estimate Stem carbon pool
# Import weekly height diameter data for 3 months (Diameter is in mm; Height is in cm)
height.dia <- read.csv("height_diameter.csv")
library(lubridate)
height.dia$Date <- parse_date_time(height.dia$Date,"d m y")
height.dia$Date = as.Date(height.dia$Date, format = "%d/%m/%Y")
height.dia <- merge(plot.summary,height.dia,by=c("pot","plot"))
height.dia = height.dia[ ,c(-1,-2)]
height.dia = height.dia[with(height.dia, order(volume,Date)), ]
height.dia <- height.dia[c("Date", "volume", "height", "diameter")]

for(i in 1:length(vols)) {
  height.dia.idn = subset(height.dia,volume==vols[i]) 
  for(j in 1:length(unique(height.dia.idn$Date))) {
    height.dia.idn.date = subset(height.dia.idn, Date == unique(height.dia.idn$Date)[j])
    height.dia.idn.date[nrow(height.dia.idn.date)+1, 2:ncol(height.dia.idn.date)] = colMeans(height.dia.idn.date[2:ncol(height.dia.idn.date)], na.rm = TRUE) # R7 = Average of leaf data
    height.dia.idn.date[nrow(height.dia.idn.date)+1, 2:ncol(height.dia.idn.date)] = apply(height.dia.idn.date[2:ncol(height.dia.idn.date)], 2, sd) # R8 = Standard deviation of leaf counts
    height.dia.idn.date$Date = height.dia.idn.date$Date[1]
    dimnames(height.dia.idn.date)[[1]] <- c(1:(nrow(height.dia.idn.date)-2), "Mean", "SD")
    if (i == 1 && j == 1) {
      height.dia.final <- height.dia.idn.date[0,]
    }
    height.dia.final[j+(i-1)*length(unique(height.dia.idn$Date)), ] <- height.dia.idn.date["Mean", ]
    height.dia.final$height_SD[j+(i-1)*length(unique(height.dia.idn$Date))] <- height.dia.idn.date["SD", 3]
    height.dia.final$dia_SD[j+(i-1)*length(unique(height.dia.idn$Date))] <- height.dia.idn.date["SD", 4]
  }
}

################### Linear regression model fitting [log(wood_mass) = b(1) + b(2)*log(dia) + b(3)*log(height)]
# Fit the model with initial data (10) and harvested data for free seedlings (7)
# Import initial seedling data
initial.data <- read.csv("seedling_initial.csv")
initial.data[nrow(initial.data)+1, 2:ncol(initial.data)] = colMeans(initial.data[2:ncol(initial.data)], na.rm = TRUE) # R7 = Average of leaf data
initial.data[nrow(initial.data)+1, 2:ncol(initial.data)] = apply(initial.data[2:ncol(initial.data)], 2, sd) # R8 = Standard deviation of leaf counts
dimnames(initial.data)[[1]] <- c(1:(nrow(initial.data)-2), "Mean", "SD")

# Import harvested seedling data for all different treatments
end.data <- read.csv("seedling_mass.csv")
end.data = end.data[ ,c(-1,-2)]
end.data = end.data[with(end.data, order(volume)), ]
for(i in 1:length(vols)) {
  end.data.idn = subset(end.data,volume==vols[i]) 
  end.data.idn[nrow(end.data.idn)+1, 2:ncol(end.data.idn)] = colMeans(end.data.idn[2:ncol(end.data.idn)], na.rm = TRUE) # R7 = Average of leaf data
  end.data.idn[nrow(end.data.idn)+1, 2:ncol(end.data.idn)] = apply(end.data.idn[2:ncol(end.data.idn)], 2, sd) # R8 = Standard deviation of leaf counts
  end.data.idn$volume = end.data.idn$volume[1]
  dimnames(end.data.idn)[[1]] <- c(1:(nrow(end.data.idn)-2), "Mean", "SD")
  if (i == 1) {
    end.data.final <- end.data.idn[0,]
  }
  end.data.final[i, ] <- end.data.idn["Mean", ]
  end.data.final$courseroot_SD[i] <- end.data.idn["SD", 2]
  end.data.final$fineroot_SD[i] <- end.data.idn["SD", 3]
  end.data.final$leafmass_SD[i] <- end.data.idn["SD", 4]
  end.data.final$stemmass_SD[i] <- end.data.idn["SD", 5]
  end.data.final$totalmass_SD[i] <- end.data.idn["SD", 6]
}

# Linear regression model fitting: log(stem_mass) = b(1) + b(2)*log(dia) + b(3)*log(height)
stemmass = c(initial.data$wood_mass[1:10], hd$stem_mass[1:7])
p = subset(height.dia.idn,Date=="2013-05-21")
height = c(initial.data$height[1:10], p$height)
diameter = c(initial.data$diameter_15[1:10], p$diameter)

model.fit = data.frame(stemmass, height, diameter)
write.csv(model.fit, file = "stemmass_height_diameter_data.csv", row.names = FALSE)
fit <- lm(log(stemmass) ~ log(diameter) + log(height), data=model.fit)
summary(fit) # show results
coefficients(fit) # model coefficients
cat("Linear regression model fitting: log(stem_mass) = ", coefficients(fit)[1], "+", coefficients(fit)[2], 
    "* log(diameter) +", coefficients(fit)[3], "* log(height)")

# Estimate the stemmass from the fitted linear regression equation
eq = function(x,y){exp(coefficients(fit)[1] + coefficients(fit)[2] * log(x)  + coefficients(fit)[3] * log(y))}
x = height.dia.final$diameter
y = height.dia.final$height
z = height.dia.final$stemmass = eq(x,y)
# height.dia.final$stemmass = exp(coefficients(fit)[1] + coefficients(fit)[2] * log(height.dia.final$diameter) 
# + coefficients(fit)[3] * log(height.dia.final$height))

# Plotting observation and modelled data
png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/Height_Stem mass.png")
plot(model.fit$height,model.fit$stemmass,col="red",main="Height vs Stem mass", pch=15)
lines(y,z,type="p",xlab="Height (cm)", ylab="Stem mass (g)",col="green", pch=20)
legend('topleft', c("Measurements", "Modelled data"), lty=1, col=c('red','green'), bty='n', cex=0.75)
dev.off()
png(file = "/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/Diameter_Stem mass.png")
plot(model.fit$diameter,model.fit$stemmass, col="red",main="Diameter vs Stem mass")
lines(x,z,type="p",col="green", xlab="Diameter (cm)", ylab="Stem mass (g)")
legend('topleft', c("Measurements", "Modelled data"), lty=1, col=c('red','green'), bty='n', cex=0.75)
dev.off()
# plot(model.fit$height,model.fit$diameter,col="red",main="Height vs Diameter")
# lines(y,x,type="p",col="green", xlab="Height (cm)", ylab="Diameter (cm)")
# legend('topleft', c("Measurements", "Modelled data"), lty=1, col=c('red','green'), bty='n', cex=0.75)

# #########
# # 3D plots for stemmass with height and diameter
# # install.packages("plotly")
# library("plotly")
# # install.packages("scatterplot3d")
# library("scatterplot3d")
# par(mfrow=c(1,1))
# q <- scatterplot3d(x, y, z, type="l", xlab="Diameter", ylab="Height", zlab="Stem mass", box=FALSE)
# # plot_ly(height.dia.final, x = ~diameter, y = ~height, z = ~stemmass, color = ~time)
# # install.packages("car")
# library(car)
# scatter3d(x = x, y = y, z = z, point.col = "blue", groups = factor(height.dia.final$volume), surface=FALSE, 
#           xlab = "Diameter (cm)", ylab = "Height (cm)", zlab = "Stem mass (g)")
# 
# # install.packages("plot3D")
# library(plot3D)
# scatter3D(diameter, height, stemmass, xlab = "Dia", ylab = "Height", zlab = "Stem mass") 
# par(new=TRUE)
# lines3D(x, y, z, xlab = "Dia", ylab = "Height", zlab = "Stem mass") 




################### Calculate daily GPP, NPP, R.plant(Plant respiration)
##use modelled C net and gross to determine respiration

##read data
# install.packages("rio")
library("rio")
finalmass = import("https://raw.githubusercontent.com/CourtneyCampany/EucPVE/master/calculated%20data/harvest_mass_means.csv")
# write.csv(finalmass, "harvest_mass_means.csv", row.names=FALSE)
# finalmass <- read.csv("/Users/kashifmahmud/WSU/ARC_project/Data_files/Carbon_balance_model_data/harvest_mass_means.csv")

# la_pred <- import("https://raw.githubusercontent.com/CourtneyCampany/EucPVE/master/calculated%20data/LApred_volume.csv")
# la_Kashif = read.csv("/Users/kashifmahmud/WSU/ARC_project/Data_files/Court_data/LA_Kashif.csv", header=F, sep=",")
# la_pred[ , 4] <- la_Kashif[ , 1] # Leaf area calculation from harvested data
# la_pred[ , 5] <- la_Kashif[ , 2] # Leaf area calculation considering changes in SLA over time
# names(la_pred)[3:5] <- c("LA", "LA", "LA") # C3 = LA_court, C4 = LA_harvest, C5 = LA_sla
# la_pred$Date <- as.Date(la_pred$Date)

la.daily.m = la.daily.melt # Leaf area calculated from harvested data
la.daily.m$volume = as.numeric(la.daily.m$volume)
for(i in length(vols):1) {
  ind = which(la.daily.m$volume %in% i)
  la.daily.m$volume[ind] = vols[i]
}
la.daily.m$Leaf_area = la.daily.m$Leaf_area / (100*100) # unit conversion from cm2 to m2

Cday_net <- read.csv("https://raw.githubusercontent.com/CourtneyCampany/EucPVE/master/calculated%20data/Aleaf_model/cday_120_clean.csv")
Cday_net$Date <- as.Date(Cday_net$Date)

Cday_gross <- read.csv("https://raw.githubusercontent.com/CourtneyCampany/EucPVE/master/calculated%20data/Aleaf_model/cday_120_clean_gross.csv")
Cday_gross$Date <- as.Date(Cday_gross$Date)

##Cday with leaf area needs self shading (use slope intercept, from 5-free vol)
sigma <- read.csv("https://raw.githubusercontent.com/CourtneyCampany/EucPVE/master/calculated%20data/M_leafarea_model.csv")

##function to generate total plant daily C gain-----------------------------------------------------------------------------
modelledC_func <- function(leafarea, shading, Cday){ #leafarea dfr, self shading dfr, and modelled C gain (gm2)
  dailyCnet <- merge(leafarea,Cday)
  ##this needs to include self shadeing (M as a linear function of leaf area)
  dailyCnet <- merge(dailyCnet, shading[, c(2,3,5)], by="volume")
  dailyCnet$M <- with(dailyCnet, b*Leaf_area+intercept)
  #calculate total daily C gain with self shading
  dailyCnet$tdcg <- with(dailyCnet, Leaf_area * carbon_day * M)
  print("successfully calculated total daily carbon with modelled Cgain, leaf area and self shading")
  
  dailyCnet$volume <- as.factor(dailyCnet$volume)
  return(dailyCnet)
}

##calculate net and gross C gain, then quantify leaf respiration--------------------------------------------------------------
dailyCnet <- modelledC_func(la.daily.m, sigma, Cday_net)
# dailyCnet <- modelledC_func(la_pred[, c(1:2,3)], sigma, Cday_net)

names(dailyCnet)[8] <- "tdc_net"
names(dailyCnet)[4] <- "cday_net"

dailyCgross <- modelledC_func(la.daily.m, sigma, Cday_gross)
# dailyCgross <- modelledC_func(la_pred[, c(1:2,3)], sigma, Cday_gross)
names(dailyCgross)[8] <- "tdc_gross"
names(dailyCgross)[4] <- "cday_gross"

dailyC <- merge(dailyCgross[, c(1:2, 4, 8)], dailyCnet[,c(1:2,4,8)], by=c("Date", "volume")) # Unit = micromol CO2 m-2 d-1
# Unit conversion to gDM m-2 d-1: (12/44) is for micromol CO2 to gC; 0.5 is for gC to gDM
dailyC[ , c(3:6)] <- dailyC[ , c(3:6)] * (12/44) / 0.5 # Unit = gDM m-2 d-1

ggplot(data = dailyC, aes(x = Date, y = tdc_net, group = volume, colour=factor(volume))) + 
  geom_point(size=1) +
  xlab("Date") +
  ylab("Net total daily C (gDM m-2 d-1)") +
  ggtitle("Net total daily C")
# with(dailyC, plot(Date,tdc_net, col=volume, xlab="date", ylab="Net total daily C (gDM m-2 d-1)"))

dailyC$R.plant <- with(dailyC, tdc_gross-tdc_net) # Unit = gDM m-2 d-1

#Calculate total seedling C gain over experiment (120d) and compare to final harvest mass C
library("doBy")
plantCnet <- doBy::summaryBy(tdc_net+R.plant+tdc_gross ~ volume, FUN=sum, data=dailyC)
par(mfrow=c(1,1))

# Plot both the results for comparison
plantCnet$volume = as.numeric(as.character(plantCnet$volume))
plot(plantCnet$volume,plantCnet$tdc_net.sum,type="p",col="red",xlab="Soil volume (L)", ylab="Net total daily carbon (g)",
     ylim=c(0,max(plantCnet$tdc_net.sum,finalmass$mass)))
lines(finalmass$volume,finalmass$mass,type="p",col="green")
legend('topleft', c("Modelled data", "Measurements from harvest"), lty=1, col=c('red','green'), bty='n', cex=0.75)
# write.csv(finalmass, "finalmass.csv", row.names = FALSE)



#################### Rdark prediction through time using rdarkq10 equation by volume
#import site weather data, rename variables, format date stuff
eucpve_met <- read.csv("https://raw.githubusercontent.com/kashifmahmud/EucPVE/master/calculated%20data/eucpve_met.csv")
eucpve_met <- eucpve_met[ , c(1,4)]
names(eucpve_met)[2] <- "temp"
eucpve_met$Date <- as.Date(eucpve_met$DateTime15)

#need to turn the datetime 15 into hms
#install.packages("lubridate")
library("lubridate")
eucpve_met$DateTime15 <- ymd_hms(eucpve_met$DateTime15)
eucpve_met$time15 <- format(eucpve_met$DateTime15, format='%H:%M:%S')

#subset by Date range of experiment
eucpve_met1 <- subset(eucpve_met[, 2:4], Date  >= "2013-01-21" & Date  <= "2013-05-21")

#Rdark Q10 equations by volume
rdarkq10 <- read.csv("https://raw.githubusercontent.com/kashifmahmud/EucPVE/master/calculated%20data/rdarkq10.csv")

#---------------------------------------------------------------------------------------------
#A_model, use temperature data from site to generate daily total leaf dark respiration to enter into the exponential Rd model
volume <- data.frame(vols)
A_model <- merge(eucpve_met1, volume)
names(A_model)[4] <- "volume"
#need to calculate Rdark through time using rdarkq10 equation by volume
A_model <- merge(A_model, rdarkq10[,1:2], by="volume")
q25_drake <- 1.86

# refit RD (gC m-2 leaf s-1)
A_model$Rd_pred <- with(A_model, rd12.3 * q25_drake^((temp-12.3)/10)) # unit (micromol CO2 m-2 leaf d-1)
A_model <- merge(A_model, hd.final[,c(1,8)], by="volume")
A_model$Rd_daily <- with(A_model, Rd_pred * sla_harvest) # unit conversion to (micromol CO2 g-1 leaf d-1)
with(A_model, plot(temp,Rd_daily, col=volume))

library("doBy")
Rd <- summaryBy(Rd_daily ~ Date+volume, data=A_model, FUN=sum, keep.names=TRUE ) # Sum of all same day Rd
Rd$Rd_daily = Rd$Rd_daily / (24*4) # Average of all 96 data from a day
Rd$Rd_daily = Rd$Rd_daily * (12/44)  # unit conversion from micromol CO2 to gC
# Rd$Rd_daily = Rd$Rd_daily / 0.5 # Unit = gDM m-2 d-1
write.csv(Rd, "Rd_pred.csv", row.names=FALSE)


############## Interpolate daily Cleaf, Cstem, Croot
library(dplyr)
library(zoo)

C.pool = Rd[ ,c(1,2)]
C.pool = C.pool[with(C.pool, order(volume)), ]

x2 <- merge(height.dia.final[ , c(1:2,7)], C.pool, by=c("Date", "volume"), all=TRUE)
for(i in 1:length(vols)) {
  x2.idn = data.frame(subset(x2[ , c(1,3)],volume==vols[i]))
  x2.idn = data.frame(Date = x2.idn$Date) %>%
    full_join(x2.idn, by = "Date") %>%
    mutate(stemmass = na.approx(stemmass))
  C.pool$stemmass[(1+(i-1)*number_of_day) : (i*number_of_day)] = x2.idn$stemmass
}
lm.daily.m = lm.daily.melt # Leaf area calculated from harvested data
lm.daily.m$volume = as.numeric(lm.daily.m$volume)
for(i in length(vols):1) {
  ind = which(lm.daily.m$volume %in% i)
  lm.daily.m$volume[ind] = vols[i]
}
C.pool <- merge(C.pool, lm.daily.m, by=c("Date", "volume"))
C.pool = C.pool[with(C.pool, order(volume)), ]

root_mass_initial = data.frame(Date = rep(as.Date("2013-01-21", each=length(vols))), volume = hd.final$volume, 
                               rootmass = rep(initial.data$root_mass[11], each=length(vols)))
hd.final$rootmass = hd.final$Croot+hd.final$Froot
root_mass_final = data.frame(Date = rep(as.Date("2013-05-21", each=length(vols))), volume = hd.final$volume, rootmass = hd.final$rootmass)
root_mass = rbind(root_mass_initial,root_mass_final)


# Estimate daily Plant respiration, rd.plant = Rd(Cleaf + Cstem + Croot); Rd = Plant respiration (g-1 plant weight) calcualted from 
# Leaf respiration and temperature (15 mins data)
x3 <- merge(root_mass, C.pool, by=c("Date", "volume"), all=TRUE)
for(i in 1:length(vols)) {
  x3.idn = data.frame(subset(x3[ , c(1,3)],volume==vols[i]))
  x3.idn = data.frame(Date = x3.idn$Date) %>%
    full_join(x3.idn, by = "Date") %>%
    mutate(rootmass = na.approx(rootmass))
  C.pool$rootmass[(1+(i-1)*number_of_day) : (i*number_of_day)] = x3.idn$rootmass
}
C.pool$totalCmass = C.pool$leafmass + C.pool$stemmass + C.pool$rootmass
C.pool$rd.plant = C.pool$totalCmass * Rd$Rd_daily


########## Comparison between whole plant respiration calculated from leaf respiration {rd.plant = Rd(Cleaf + Cstem + Croot)} and 
########## modelled Cday (gross and net) from Court's data
ggplot(dailyC, aes(x = Date, y = R.plant, group = volume, colour=factor(volume))) + 
  geom_point(size=0.5) +
  geom_line(data = C.pool, aes(x = Date, y = rd.plant, group = volume, colour=factor(volume))) + 
  xlab("Date") +
  ylab("Plant respiration (gC plant-1 d-1)") +
  # scale_colour_manual(values=c("points","line")) + 
  ggtitle("Lines = From measured leaf resp; Points = From modelled Cday (gross-net)")





################### Find leaf storage (tnc, fortnightly data) for corresponding Dates (from Court's leaf_data file: represents Gas measurement campaign)
keeps <- c("Date", "volume", "starch_mgperg", "sugars_mgperg")
tnc.gas = leaf.data[ , keeps, drop = FALSE]
names(tnc.gas)[1:2] <- c("Date", "volume")
tnc.gas$tnc = tnc.gas$starch_mgperg + tnc.gas$sugars_mgperg
keeps <- c("Date", "volume", "tnc") # Unit = mg g-1leaf
tnc.gas = tnc.gas[ , keeps, drop = FALSE]

for(i in 1:length(vols)) {
  tnc.idn = subset(tnc.gas,volume==vols[i]) 
  for(j in 1:length(unique(tnc.idn$Date))) {
    tnc.idn.date = subset(tnc.idn, Date == unique(tnc.idn$Date)[j])
    tnc.idn.date[nrow(tnc.idn.date)+1, 2:ncol(tnc.idn.date)] = colMeans(tnc.idn.date[2:ncol(tnc.idn.date)], na.rm = TRUE) # R7 = Average of tnc
    tnc.idn.date[nrow(tnc.idn.date)+1, 2:ncol(tnc.idn.date)] = apply(tnc.idn.date[2:ncol(tnc.idn.date)], 2, sd) # R8 = Standard deviation of tnc
    tnc.idn.date$Date = tnc.idn.date[1,1]
    dimnames(tnc.idn.date)[[1]] <- c(1:(nrow(tnc.idn.date)-2), "Mean", "SD")
    if (i == 1 && j == 1) {
      tnc.final <- tnc.idn.date[0,]
    }
    tnc.final[j+(i-1)*length(unique(tnc.idn$Date)), ] <- tnc.idn.date["Mean", ]
    tnc.final$tnc_SD[j+(i-1)*length(unique(tnc.idn$Date))] <- tnc.idn.date["SD", 3]
  }
}
# Unit conversion from (mg g-1leaf) to  (g plant-1)
lm.daily.m$Date = as.Date(lm.daily.m$Date)
lm.daily.gas = lm.daily.m[lm.daily.m$Date %in% as.Date(c(unique(tnc.final$Date))), ]
tnc.final$tnc = tnc.final$tnc * lm.daily.gas$leafmass / 1000 # Unit = g plant-1
tnc.final$tnc_SD = tnc.final$tnc_SD * lm.daily.gas$leafmass / 1000 # Unit = g plant-1

# Plot tnc data
ggplot(data = dailyC, aes(x = Date, y = tdc_net, group = volume, colour=factor(volume))) + 
  geom_point(size=1) +
  xlab("Date") +
  ylab("Total N concentration (g plant-1)") +
  ggtitle("Total N concentration (g plant-1)")



