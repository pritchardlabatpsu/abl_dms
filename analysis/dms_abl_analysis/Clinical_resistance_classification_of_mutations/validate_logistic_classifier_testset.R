#title: "Validating the logistic classifier on an independent test set"
#author: "Marta Tomaszkiewicz"
#date: "2024-10-24"


###Validating the logistic classifier on an independent test set: 
#1 Fit concentration–response curves, compute predicted relative viability at each dose, 
#2 Derive the optimal probability cutoff (Youden index)
#3 Derive corresponding viability cutoff to classify variants as resistant vs. non-resistant

## Set up
library(dr4pl)
library(dplyr)
library(pROC)
library(ggplot2)

## Import and parse data
ngr.df = read.csv("merged_table_imatinib_full_kinase_ngr_with_resmuts.csv",
                  header=T,stringsAsFactors=F)

# Known non-resistant mutations
non_resistant_muts = c('A337V', 'P465S', 'V468F', 'I502L', 'V299L')

#Test set of 15 Literature_curated mutations with
resistant_muts <-c("L387M", "V379I", "L384M", "F311L", "E279K", "M388L", 
                   "F317V", "L387F", "S348L", "D325N", "E255D", "E275K", 
                   "F317C", "L248R", "T315V")

ngr.df <- ngr.df %>%
  filter(species %in% c(non_resistant_muts, resistant_muts))

screen.df = ngr.df

screen.df$resmuts = screen.df$species %in% resistant_muts & !(screen.df$species %in% non_resistant_muts)

write.csv(screen.df, "test_dataset.csv", row.names = FALSE)

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

write.csv(screen.df, "IC50_all_corrected.csv", row.names = FALSE)

## Go through range of drug doses spanning low to high dose

dose.rnge <- sort(unique(c(seq(10, 10000, by = 10), 444, 760, 916)))

# Add resistance ground truth as a logical vector, explicitly handling non-resistant cases
screen.df$resmuts = screen.df$species %in% resistant_muts & !(screen.df$species %in% non_resistant_muts)

# Initialize results dataframe
mdl.df <- data.frame(
  dose = dose.rnge,
  rel.viab.AUC = NA,
  optimal.prob.threshold = NA,
  accuracy = NA,
  resistance_threshold = NA
)

# Initialize lists to store GLM fits and predictions
glm_models <- list()
glm_predictions <- list()

# Initialize predictions as empty data frame
predictions <- data.frame()

# Loop over each dose
for (dose in dose.rnge) {
  
  # Compute predicted response and alpha per variant
  for (i in 1:nrow(screen.df)) {
    IC50.i <- screen.df$IC50[i]
    hill.i <- screen.df$Hill.slope[i]
    
    pred.resp <- MeanResponse(theta = c(1, IC50.i, hill.i, 0),
                              x = dose)
    screen.df[i, paste0(dose, ".rel.viab")] <- pred.resp
    screen.df[i, paste0(dose, ".alpha")] <- -log(pred.resp) / t.assay
  }
  
  # Prepare sub-data for logistic regression
  sub.df <- data.frame(
    resmuts = screen.df$resmuts,
    rel.viab = screen.df[[paste0(dose, ".rel.viab")]],
    alpha = screen.df[[paste0(dose, ".alpha")]]
  )
  
  # Fit logistic regression model
  glm.fits <- glm(resmuts ~ rel.viab, data = sub.df, family = binomial)
  glm.probs <- predict(glm.fits, type = "response")
  
  # Store model and predictions
  glm_models[[as.character(dose)]] <- glm.fits
  glm_predictions[[as.character(dose)]] <- glm.probs
  
  # Compute ROC using predicted probabilities
  ROC.rel.viab <- roc(sub.df$resmuts,
                      glm.probs,
                      direction="auto",
                      levels = c(F, T))
  
  plot(ROC.rel.viab, main = paste("ROC Curve - dose", dose),
       col = "blue", lwd = 2)
  text(0.6, 0.2, paste("AUC:", round(auc(ROC.rel.viab), 2)),
       col = "red", cex = 1.2)
  
  # Extract optimal threshold based on Youden index
  roc_thresholds <- data.frame(
    thresholds = ROC.rel.viab$thresholds,
    sensitivity = ROC.rel.viab$sensitivities,
    specificity = ROC.rel.viab$specificities
  )
  roc_thresholds$youden_index <- roc_thresholds$sensitivity + roc_thresholds$specificity - 1
  optimal_threshold <- roc_thresholds$thresholds[which.max(roc_thresholds$youden_index)]
  
  # Classify variants and compute accuracy
  glm.pred <- ifelse(glm.probs > optimal_threshold, 1, 0)
  accuracy <- mean(glm.pred == sub.df$resmuts)
  
  # Convert optimal probability threshold → rel.viab cutoff
  coefs <- coef(glm.fits)
  intercept <- coefs[1]
  slope <- coefs[2]
  logit_val <- log(optimal_threshold / (1 - optimal_threshold))
  relviab_cutoff <- (logit_val - intercept) / slope
  
  # Save metrics
  mdl.df$rel.viab.AUC[mdl.df$dose == dose] <- as.numeric(auc(ROC.rel.viab))
  mdl.df$optimal.prob.threshold[mdl.df$dose == dose] <- optimal_threshold
  mdl.df$accuracy[mdl.df$dose == dose] <- accuracy
  mdl.df$resistance_threshold[mdl.df$dose == dose] <- relviab_cutoff
  
  # Store predictions for this dose
  pred_df <- data.frame(
    dose = dose,
    species = screen.df$species,
    rel.viab = sub.df$rel.viab,
    actual = screen.df$resmuts,
    glm_prob = glm.probs,
    glm_pred = glm.pred,
    correct = glm.pred == sub.df$resmuts,
    threshold_used = optimal_threshold
  )
  
  predictions <- rbind(predictions, pred_df)
}

# Inspect final model summary
print(mdl.df)

write.csv(predictions, "TEST_Predictions_per_dose.csv", row.names = FALSE)
write.csv(mdl.df, "TEST_resistance_threshold_per_dose.csv", row.names = FALSE)
