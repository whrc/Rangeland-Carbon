library(BayesianTools)
library(sensitivity)
library(robustbase)
setwd("C:/Users/username/path-to-the-calibration-data-folder")
data <- read.table('file-with-fPAR-and-covariates-for-GPP-model.csv', header=T, sep = ",") 
data <- na.omit(data) # Apply if there are NAs 
data$GPP_measured <- data$GPP_measured*24*3600/(10^6)*12 # Convert unit to g C/m2 if not done so
data$site_num = as.numeric(as.factor(with(data,paste(Site)))) 
site <- unique(data[c("site_num","Site")])
obs = data['GPP_measured']
sd(obs[,1])
length(unique(data$Site))

# Set and calibrate GPP parameters
refPars <- VSEMgetDefaults()
refPars[12,] <- c(1.001, 1, 1.4)
rownames(refPars)[12] <- "LUEmax"
refPars[13,] <- c(-37.607, -65, -15)
rownames(refPars)[13] <- "MIN_Tmn"
refPars[14,] <- c(13.215, 5, 24)
rownames(refPars)[14] <- "MAX_Tmx"
refPars[15,] <- c(0.348, 0, 1.2)
rownames(refPars)[15] <- "MIN_VPD"
refPars[16,] <- c(1.841, 1.2, 5)
rownames(refPars)[16] <- "MAX_VPD"
refPars[17,] <- c(0.085, 0, 0.2)
rownames(refPars)[17] <- "MIN_SMrz"
refPars[18,] <- c(0.998, 0.85, 1)
rownames(refPars)[18] <- "MAX_SMrz"
refPars[19,] <- c(3.30, 0.1, 5)
rownames(refPars)[19] <- "error-sd"
refPars=refPars[12:19,]
refPars

# Define functions to calculate GPP
calcGPP <- function(param, tsoil_in, sm_in, vpd_in, SW_IN, fPAR){
  LUEmax = param[1]
  MIN_Tmn = param[2]
  MAX_Tmx = param[3]
  MIN_VPD = param[4]
  MAX_VPD = param[5]
  MIN_SMrz = param[6]  
  MAX_SMrz = param[7]  
  tmult = ifelse (tsoil_in < MIN_Tmn, 0.0, 
                  ifelse (tsoil_in > MAX_Tmx, 1.0,
                          (tsoil_in-(MIN_Tmn))/(MAX_Tmx-(MIN_Tmn))))
  smult = ifelse (sm_in < MIN_SMrz, 0.0,
                  ifelse (sm_in > MAX_SMrz, 1.0,
                          (sm_in-(MIN_SMrz))/(MAX_SMrz-(MIN_SMrz))))
  wmult = ifelse (vpd_in < MIN_VPD, 0.0, 
                  ifelse (vpd_in > MAX_VPD, 1.0,
                          (vpd_in-(MIN_VPD))/(MAX_VPD-(MIN_VPD))))
  LUE = LUEmax * tmult * smult * wmult
  GPP = LUE * SW_IN * fPAR * 0.45
  GPP = data.frame(GPP)
  return(GPP)
}
result = calcGPP(refPars$best, data$TS, data$SWC2, data$VPD, data$SW_IN_NLDAS, data$STARFM_fPAR) 

# Likelihopod function
parSel = c(1:7)
likelihood <- function(par, sum = TRUE){
  x = refPars$best
  x[parSel] = par
  predicted <- calcGPP(x, data$TS, data$SWC2, data$VPD, data$SW_IN_NLDAS, data$STARFM_fPAR)
  diff <- c(predicted[,1] - obs[,1])
  llValues <- dnorm(diff, sd = x[8], log = TRUE)  
  if (sum == FALSE) return(llValues)
  else return(sum(llValues))
}
likelihood(refPars$best[parSel])

# Provide prior and setup
prior <- createUniformPrior(lower = refPars$lower[parSel], 
                            upper = refPars$upper[parSel], best = refPars$best[parSel])
bayesianSetup <- createBayesianSetup(likelihood, prior, names = rownames(refPars)[parSel])
settings <- list(iterations = 10000, nrChains = 3) 
plotSensitivity(bayesianSetup, selection = NULL, equalScale = T)

# Run Bayesian MCMC
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)

# Check results
par(mar=c(1,1,1,1))
summary(out)
marginalPlot(out)
gelmanDiagnostics(out)

# Check prior and posterior predictions
createPredictions <- function(par){
  x = refPars$best
  x[parSel] = par
  predicted <- calcGPP(x, data$TS, data$SWC2, data$VPD, data$SW_IN_NLDAS, data$STARFM_fPAR) 
  return(predicted[,1])
}
createError <- function(mean, par){
  return(rnorm(length(mean), mean = mean, sd = refPars$best[8]))
}
plotTimeSeriesResults(sampler = out, model = createPredictions, observed = obs[,1],
                      error = createError, prior = TRUE, main = "Prior predictive")
plotTimeSeriesResults(sampler = out, model = createPredictions, observed = obs[,1],
                      error = createError, main = "Posterior predictive")

# Get predicted vs. observed
post <- data.frame(getSample(out))
post_mean <- colMedians(as.matrix(post))
newresult = calcGPP(post_mean, data$TS, data$SWC2, data$VPD, data$SW_IN_NLDAS, data$STARFM_fPAR)
comparison_before <- data.frame(obs,result,data$Site)
comparison_after <- data.frame(obs,newresult,data$Site)
par(mar=c(2.5,2.5,2.5,2.5))
lr_before <- lm(GPP_measured ~ GPP, data = comparison_before)
summary(lr_before)
plot(comparison_before$GPP_measured,comparison_before$GPP,col=factor(comparison_before$data.Site),
     main = "Observed vs. predicted GPP before calibration", 
     xlab = "Observed GPP (%)",
     ylab = "Predicted GPP (%)",
     xlim = c(0,25), # Set bounds for visualization
     ylim = c(0,25))
text(5,15, paste0("Model R2= ",round(summary(lr_before)$r.squared,digits = 3)),cex = 1.5,col="red") # Change parameters for visualization
legend("topright",legend=levels(factor(comparison_before$data.Site)),pch = 1,
       col= factor(levels(factor(comparison_before$data.Site))))
RMSE <- sqrt(mean((comparison_before$GPP_measured - comparison_before$GPP)^2))
RMSE
mbe <- mean(comparison_before$GPP_measured - comparison_before$GPP)
mbe
lr_after <- lm(GPP_measured ~ GPP, data = comparison_after)
summary(lr_after)
plot(comparison_after$GPP_measured,comparison_after$GPP,col=factor(comparison_after$data.Site),
     main = "Observed vs. predicted GPP after calibration", 
     xlab = "Observed GPP (%)",
     ylab = "Predicted GPP (%)",
     xlim = c(0,25),
     ylim = c(0,25))
text(5,15, paste0("Model R2= ",round(summary(lr_after)$r.squared,digits = 3)),cex = 1.5,col="red")
legend("topright",legend=levels(factor(comparison_after$data.Site)),pch = 1,
       col= factor(levels(factor(comparison_after$data.Site))))
RMSE <- sqrt(mean((comparison_after$GPP_measured - comparison_after$GPP)^2))
RMSE
mbe <- mean(comparison_after$GPP_measured - comparison_after$GPP)
mbe

# Get distribution of model parameters
all_out <- data.frame(rbind(out[[1]]$Z,out[[2]]$Z,out[[3]]$Z))
quantiles_X1 <- quantile(all_out$X1, probs = c(.05, .25, .50, .75, .95))
median_X1 <- median(all_out$X1)
quantiles_X2 <- quantile(all_out$X2, probs = c(.05, .25, .50, .75, .95))
median_X2 <- median(all_out$X2)
quantiles_X3 <- quantile(all_out$X3, probs = c(.05, .25, .50, .75, .95))
median_X3 <- median(all_out$X3)
quantiles_X4 <- quantile(all_out$X4, probs = c(.05, .25, .50, .75, .95))
median_X4 <- median(all_out$X4)
quantiles_X5 <- quantile(all_out$X5, probs = c(.05, .25, .50, .75, .95))
median_X5 <- median(all_out$X5)
quantiles_X6 <- quantile(all_out$X6, probs = c(.05, .25, .50, .75, .95))
median_X6 <- median(all_out$X6)
quantiles_X7 <- quantile(all_out$X7, probs = c(.05, .25, .50, .75, .95))
median_X7 <- median(all_out$X7)
sd_X1 <- sd(all_out$X1)
sd_X2 <- sd(all_out$X2)
sd_X3 <- sd(all_out$X3)
sd_X4 <- sd(all_out$X4)
sd_X5 <- sd(all_out$X5)
sd_X6 <- sd(all_out$X6)
sd_X7 <- sd(all_out$X7)
results_table <- data.frame(
  Variable = c("LUEmax","MIN_Tmn","MAX_Tmx","MIN_VPD","MAX_VPD","MIN_SMrz","MAX_SMrz"),
  Quantile_5th = c(quantiles_X1[1], quantiles_X2[1], quantiles_X3[1], quantiles_X4[1], quantiles_X5[1], quantiles_X6[1], quantiles_X7[1]),
  Quantile_25th = c(quantiles_X1[2], quantiles_X2[2], quantiles_X3[2], quantiles_X4[2], quantiles_X5[2], quantiles_X6[2], quantiles_X7[2]),
  Median = c(median_X1, median_X2, median_X3, median_X4, median_X5, median_X6, median_X7),
  Quantile_75th = c(quantiles_X1[4], quantiles_X2[4], quantiles_X3[4], quantiles_X4[4], quantiles_X5[4], quantiles_X6[4], quantiles_X7[4]),
  Quantile_95th = c(quantiles_X1[5], quantiles_X2[5], quantiles_X3[5], quantiles_X4[5], quantiles_X5[5], quantiles_X6[5], quantiles_X7[5]),
  Standard_Deviation = c(sd_X1, sd_X2, sd_X3, sd_X4, sd_X5, sd_X6, sd_X7)
)
print(results_table)


