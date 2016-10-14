# Analysing Plant Storage (TNC) partitioning for Cstorage pool prediction
# Developed by Kashif Mahmud and Belinda Medlyn (October 2016)
# k.mahmud@westernsydney.edu.au

rm(list=ls())
setwd("/Users/kashifmahmud/WSU/ARC_project/CBM/Data_files")

##### TNC partitioning to tree organs (without considering organ biomass)
# data file is downloaded from "https://hiev.uws.edu.au/data_files/99499" -
# GHS30_Eglob-TxCxW_carbohydrates_20110117-20110321_L1.csv
carbohydrates = read.csv("Duan_carbohydrates.csv")
carbohydrates = subset(carbohydrates,CO2 == 400 & Temp == "Amb" & Water == "Well watered")
carbohydrates$tnc = carbohydrates$StarchW + carbohydrates$SolSugW # Unit = mg of tnc per g of dry weight biomass
carbohydrates$tnc = carbohydrates$tnc /10 # Unit = % of dry weight biomass


##### Total TNC calculation considering tree organ biomass partitioning
# data file is downloaded from "https://hiev.uws.edu.au/data_files/99503" -
# GHS30_Eglob-TxCxW_harvest_20110117-20110321_L1.csv
harvest = read.csv("Duan_harvest.csv")
harvest = subset(harvest,CO2 == 400 & Temp == "Amb" & Water == "Well watered")

leaf.tnc = subset(carbohydrates,Organ == "Leaf") # Unit = % of dry weight leafmass
stem.tnc = subset(carbohydrates,Organ == "Stem") # Unit = % of dry weight stemmass
root.tnc = subset(carbohydrates,Organ == "Root") # Unit = % of dry weight rootmass

tnc = data.frame(harvest$LeafDW,leaf.tnc$tnc,harvest$StemDW,stem.tnc$tnc,harvest$RootDW,root.tnc$tnc)
names(tnc) <- c("leafDW","leaf.tnc","stemDW","stem.tnc","rootDW","root.tnc") 
tnc$total.leaf.tnc = tnc$leaf.tnc * tnc$leafDW / 100 # Unit = gC
tnc$total.stem.tnc = tnc$stem.tnc * tnc$stemDW / 100 # Unit = gC
tnc$total.root.tnc = tnc$root.tnc * tnc$rootDW / 100 # Unit = gC

tnc$leaf_to_all = tnc$total.leaf.tnc / (tnc$total.leaf.tnc + tnc$total.stem.tnc + tnc$total.root.tnc) * 100 # Unit = %
tnc$stem_to_all = tnc$total.stem.tnc / (tnc$total.leaf.tnc + tnc$total.stem.tnc + tnc$total.root.tnc) * 100 # Unit = %
tnc$root_to_all = tnc$total.root.tnc / (tnc$total.leaf.tnc + tnc$total.stem.tnc + tnc$total.root.tnc) * 100 # Unit = %

matplot(tnc[ , 10:12],type = c("b"),pch=1,col = 1:3,xlab="sample number",ylab="ratios",main="Storage (TNC) partitioning") #plot
legend("topright", legend = c("Leaf TNC ratio", "Stem TNC ratio", "Root TNC ratio"), col=1:3, pch=0.75) # optional legend

tnc[nrow(tnc)+1, ] = colMeans(tnc, na.rm = TRUE) # R7 = Average of data
tnc[nrow(tnc)+1, ] = apply(tnc, 2, sd) # R8 = Standard deviation of data
dimnames(tnc)[[1]] <- c(1:6, "mean", "SD")
write.csv(tnc, file = "tnc_partitioning.csv", row.names = FALSE)
cat("Total Leaf TNC : Total Stem TNC : Total Root TNC = ", tnc$leaf_to_all[7], ":", tnc$stem_to_all[7], 
    ":", tnc$root_to_all[7])


