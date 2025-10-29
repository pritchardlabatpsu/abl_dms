#title: "Predicting IC50s and relative viabilities from net growth rates"
#author: "Scott Leighow"
#date: "2024-10-24"

### Predicting IC50s from net growth data and predicting  relative viability estimates###
## Set up
library(dr4pl)
library(dplyr)
library(pROC)
library(ggplot2)

## Import and parse data
ngr.df = read.csv("merged_table_imatinib_full_kinase_ngr_with_resmuts.csv",
                  header=T,stringsAsFactors=F)

screen.df = ngr.df

## For each variant, fit dose response curve

# Two replicates at three doses
doses = c(rep(screen.df$dose.low[1],2),
          rep(screen.df$dose.medium[1],2),
          rep(screen.df$dose.high[1],2))

# Initiate dr4pl fit columns
screen.df$IC50 = NA
screen.df$Hill.slope = NA
screen.df$convergence = NA
screen.df$RSS = NA

for (i in 1:nrow(screen.df)) {
  
  species.i = screen.df$species[i]
  
  ## Convert net growth rates into expected relative viability estimates for a 3 day IC50 assay
  mean.netgr.nodrug = mean(screen.df$netgr_obs_rep1.il3[i],
                           screen.df$netgr_obs_rep2.il3[i])
  
  netgrs = c(screen.df$netgr_corr_rep1.low[i],
             screen.df$netgr_corr_rep2.low[i],
             screen.df$netgr_corr_rep1.medium[i],
             screen.df$netgr_corr_rep2.medium[i],
             screen.df$netgr_corr_rep1.high[i],
             screen.df$netgr_corr_rep2.high[i])
  
  
  t.assay = 72 # in hours

  rel.viabs = exp(netgrs*t.assay)/exp(0.055*t.assay)
  
  ## Fit dose response
  
  # Force max and min values to be (near) 1 and 0, respectively
  fit = dr4pl(rel.viabs~doses,
              method.init="logistic",
              init.parm = dr4pl_theta(.99,100,-2,0.01),
              upperl=c(1,Inf,Inf,0.02),lowerl=c(0.98,-Inf,-Inf,0))
  
  rsdls = residuals(fit)
  
  plot(fit)
  
  ## Record parameters
  screen.df$IC50[i] = coef(fit)[2]
  screen.df$Hill.slope[i] = coef(fit)[3]
  screen.df$convergence[i] = fit$convergence
  screen.df$RSS[i] = sum(rsdls^2)
  
}

## Go through range of drug doses spanning low to high dose

dose.rnge <- sort(unique(c(seq(10, 10000, by = 10), 444, 760, 916)))

#Calculate rel.viab and alpha for Each Variant and Dose
for (dose in dose.rnge) {
  for (i in 1:nrow(screen.df)) {
    IC50.i <- screen.df$IC50[i]
    hill.i <- screen.df$Hill.slope[i]
    
    pred.resp <- MeanResponse(theta = c(1, IC50.i, hill.i, 0), x = dose)
    screen.df[i, paste0(dose, '.rel.viab')] <- pred.resp
    screen.df[i, paste0(dose, '.alpha')] <- -log(pred.resp) / 72  # t.assay
  }
}


write.csv(screen.df, "all_variants_w_IC50_relative_viability.csv", row.names = FALSE)
