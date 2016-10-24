AUTHOR: Kashif Mahmud and Belinda Medlyn
k.mahmud@westernsydney.edu.au
01/11/2016

Carbon Balance Model (CBM)

A simple Carbon Balance Model (CBM) for an individual free growing plant seedling (Eucalyptus Tereticornis) for daily time scale:

The folder "Codes" includes several R scripts to process and run MCMC for CBM: 

1. R script “Carbon_balance_model_Kashif.R” for data processing

2. R script “MCMC_CBM_kashif_v6.R” performing a first test case of MCMC applying in a simple carbon balance model with a very simple synthetic data set (having random gap data). 

3. R script “MCMC_CBM_kashif_v8.R” where MCMC is applied with soil manipulation pot experiment data for free seedling, however need to twick the parameters for soil restricted cases.

4. “MCMC_CBM_kashif_v9.R” script applies MCMC with soil manipulation pot experiment data for free seedling. This version considers only daily parameter set for 5 variables ("k","Y","af","as","sf”).

5. The latest version “MCMC_CBM_kashif_v10.R” put on MCMC with soil manipulation pot experiment data for all treatments (including the free seedling). This version considers either daily/weekly/monthly/just one parameter set etc. for all 5 variables ("k","Y","af","as","sf”). So we can set the parameters for various time frames. 
  

