library(BayesianTools)
library(sensitivity)
library(robustbase)
setwd("C:/Users/username/path-to-the-calibration-data-folder")
data <- read.table('file-with-fPAR-and-covariates-for-GPP-model.csv', header=T, sep = ",") 
data <- na.omit(data) # Apply if there are NAs 
data$GPP_measured <- data$GPP_measured*24*3600/(10^6)*12 # Convert unit to g C/m2 if not done so
data$site_num = as.numeric(as.factor(with(data,paste(Site)))) 
site <- unique(data[c("site_num","Site")])

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

# Calculate error metrics
RMSE_b <- data.frame()
r2_b <- data.frame()
MBE_b <- data.frame()
RMSE_a <- data.frame()
r2_a <- data.frame()
MBE_a <- data.frame()
posterior <- data.frame()
for (i in 1:length(unique(data$site_num))){
  data_i <- data[which(data$site_num == i),]
  obs = data_i['GPP_measured']
  likelihood <- function(par, sum = TRUE){
    x = refPars$best
    x[parSel] = par
    predicted <- calcGPP(x, data_i$TS, data_i$SWC2, data_i$VPD, data_i$SW_IN_NLDAS, data_i$STARFM_fPAR)
    diff <- c(predicted[,1] - obs[,1])
    llValues <- dnorm(diff, sd = x[8], log = TRUE)  
    if (sum == FALSE) return(llValues)
    else return(sum(llValues))
  }
  bayesianSetup <- createBayesianSetup(likelihood, prior, names = rownames(refPars)[parSel])
  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
  post <- data.frame(getSample(out))
  post_mean <- colMeans(post)
  result = calcGPP(refPars$best, data_i$TS, data_i$SWC2, data_i$VPD, data_i$SW_IN_NLDAS, data_i$STARFM_fPAR)
  newresult = calcGPP(post_mean, data_i$TS, data_i$SWC2, data_i$VPD, data_i$SW_IN_NLDAS, data_i$STARFM_fPAR)
  comparison_before <- data.frame(obs,result)
  comparison_after <- data.frame(obs,newresult)
  RMSE_before <- sqrt(mean((comparison_before$GPP_measured - comparison_before$GPP)^2))
  r2_before <- cor(comparison_before$GPP_measured, comparison_before$GPP)^2
  mbe_before <- mean(comparison_before$GPP_measured - comparison_before$GPP)
  RMSE_after <- sqrt(mean((comparison_after$GPP_measured - comparison_after$GPP)^2))
  r2_after <- cor(comparison_after$GPP_measured, comparison_after$GPP)^2
  mbe_after <- mean(comparison_after$GPP_measured - comparison_after$GPP)
  RMSE_b <- rbind(RMSE_b, RMSE_before)
  r2_b <- rbind(r2_b, r2_before)
  MBE_b <- rbind(MBE_b, mbe_before)
  RMSE_a <- rbind(RMSE_a, RMSE_after)
  r2_a <- rbind(r2_a, r2_after)
  MBE_a <- rbind(MBE_a, mbe_after)
  posterior <- rbind(posterior,post_mean)
}
colnames(posterior) <- c("LUEmax", "MIN_Tmn", "MAX_Tmx", "MIN_VPD", "MAX_VPD", "MIN_SMrz", "MAX_SMrz")
colnames(RMSE_b) <- "RMSE_b"
colnames(r2_b) <- "r2_b"
colnames(MBE_b) <- "MBE_b"
colnames(RMSE_a) <- "RMSE_a"
colnames(r2_a) <- "r2_a"
colnames(MBE_a) <- "MBE_a"
individual_result <- data.frame(site, posterior,RMSE_b,r2_b,MBE_b,RMSE_a,r2_a,MBE_a)


