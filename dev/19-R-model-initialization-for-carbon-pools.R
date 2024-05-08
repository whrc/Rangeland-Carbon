library(dplyr)
library(BayesianTools)
library(sensitivity)
library(robustbase)
setwd("C:/Users/username/path-to-the-calibration-data-folder")
spin <- read.table('file-with-inputs-needed-for-model-spinup', header=T, sep = ",")
spin$Site <- sapply(strsplit(spin$Site, ".csv"),`[`, 1)
spin$site_num = as.numeric(as.factor(with(spin,paste(Site))))
site <- unique(spin[c("site_num","Site")])
site <- site %>% arrange(site_num)
listG <- c("LS1","Wjs","SCg","SdH","xAE","xCP","xDC","xKA","xKZ","xNG","xWD",
           "AR1","AR2","ARb","ARc","Aud","BMM","BRG","Bkg","Ctn","Dia","FPe",
           "Hn2","IB2","KFS","KLS","KM2","KM3","KM4","Kon","RFW","Ro4","Shd",
           "Var","A32","Snd","Tx1","Tx2", "Snf") # Add a list of sites for the specific vegetation type
listS <- c("LS2","xJR","xMB","xNQ","xYE","Jo1","Rls","Hn1","Hn3",
           "Rms","Rwe","Rwf","Rws","SRM","Seg","Ses","Wkg","Wdn","Cop")
listT <- c("Mpj","Ton","xCL","xSJ","Fwf","CZ1")

# Calculate GPP
# Use default parameters before calibration to get initial model spin-up results
refPars <- VSEMgetDefaults()
refPars[12,] <- c(1.030, 1, 1.4)
rownames(refPars)[12] <- "LUEmax"
refPars[13,] <- c(-51.743, -60, -30)
rownames(refPars)[13] <- "MIN_Tmn"
refPars[14,] <- c(17.258, 12, 22)
rownames(refPars)[14] <- "MAX_Tmx"
refPars[15,] <- c(0.439, 0, 1.5)
rownames(refPars)[15] <- "MIN_VPD"
refPars[16,] <- c(1.734, 1.5, 10)
rownames(refPars)[16] <- "MAX_VPD"
refPars[17,] <- c(0.064, 0, 0.2)
rownames(refPars)[17] <- "MIN_SMrz"
refPars[18,] <- c(0.976, 0.85, 1)
rownames(refPars)[18] <- "MAX_SMrz"
refPars[19,] <- c(3.30, 0.1, 5)
rownames(refPars)[19] <- "error-sd"
refPars1=refPars[12:19,]
refPars1 # Perennial-annual
refPars <- VSEMgetDefaults()
refPars[12,] <- c(1.234, 1, 1.6)
rownames(refPars)[12] <- "LUEmax"
refPars[13,] <- c(-9.108, -15, 0)
rownames(refPars)[13] <- "MIN_Tmn"
refPars[14,] <- c(16.409, 10, 25)
rownames(refPars)[14] <- "MAX_Tmx"
refPars[15,] <- c(0.318, 0, 1.2)
rownames(refPars)[15] <- "MIN_VPD"
refPars[16,] <- c(1.776, 1.2, 15)
rownames(refPars)[16] <- "MAX_VPD"
refPars[17,] <- c(0.006, 0, 0.2)
rownames(refPars)[17] <- "MIN_SMrz"
refPars[18,] <- c(0.87, 0.55, 1)
rownames(refPars)[18] <- "MAX_SMrz"
refPars[19,] <- c(1.55, 0.1, 5)
rownames(refPars)[19] <- "error-sd"
refPars2=refPars[12:19,]
refPars2 # Grass-shrub
refPars <- VSEMgetDefaults()
refPars[12,] <- c(1.100, 1, 1.5)
rownames(refPars)[12] <- "LUEmax"
refPars[13,] <- c(-27.192, -30, -15)
rownames(refPars)[13] <- "MIN_Tmn"
refPars[14,] <- c(36.283, 30, 45)
rownames(refPars)[14] <- "MAX_Tmx"
refPars[15,] <- c(0.668, 0, 1.5)
rownames(refPars)[15] <- "MIN_VPD"
refPars[16,] <- c(1.854, 1.5, 5)
rownames(refPars)[16] <- "MAX_VPD"
refPars[17,] <- c(0.109, 0, 0.25)
rownames(refPars)[17] <- "MIN_SMrz"
refPars[18,] <- c(0.924, 0.75, 1)
rownames(refPars)[18] <- "MAX_SMrz"
refPars[19,] <- c(1.92, 0.1, 5)
rownames(refPars)[19] <- "error-sd"
refPars3=refPars[12:19,]
refPars3 # Grass-tree

# Define function to calculate GPP
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

# Get GPP estimates for the spinup files
result_spin1 = calcGPP(refPars1$best, spin$TS, spin$SWC2, spin$VPD, spin$SW_IN_NLDAS, spin$STARFM_fPAR) 
result_spin2 = calcGPP(refPars2$best, spin$TS, spin$SWC2, spin$VPD, spin$SW_IN_NLDAS, spin$STARFM_fPAR) 
result_spin3 = calcGPP(refPars3$best, spin$TS, spin$SWC2, spin$VPD, spin$SW_IN_NLDAS, spin$STARFM_fPAR) 
spin <- data.frame(spin,result_spin1,result_spin2,result_spin3)
spin[which(spin$Site %in% listS),]$GPP <- spin[which(spin$Site %in% listS),]$GPP.1
spin[which(spin$Site %in% listT),]$GPP <- spin[which(spin$Site %in% listT),]$GPP.2

# Change unit to match with the NEE model, if needed
spin$SWC1 <- spin$SWC1 *100
spin$SWC2 <- spin$SWC2 *100
spin$TS <- spin$TS + 273.15


# Set default parameters for spinup and NEE model 
refPars <- VSEMgetDefaults()
refPars[12,] <- c(0.46, 0.4, 0.5) 
rownames(refPars)[12] <- "F_Ra"
refPars[13,] <- c(0.66, 0.5, 0.8) 
rownames(refPars)[13] <- "F_Root"
refPars[14,] <- c(0.0035, 0.0005, 0.08) 
rownames(refPars)[14] <- "k_AGL"
refPars[15,] <- c(0.0027, 0.0005, 0.07) 
rownames(refPars)[15] <- "k_BGL"
refPars[16,] <- c(0.00066, 0.00005, 0.0013) 
rownames(refPars)[16] <- "k_POC"
refPars[17,] <- c(0.000255, 0.00001, 0.0005) 
rownames(refPars)[17] <- "k_HOC"
refPars[18,] <- c(0.279, 0.05, 0.5) 
rownames(refPars)[18] <- "CSE_AGL"
refPars[19,] <- c(0.264, 0.05, 0.5)
rownames(refPars)[19] <- "CSE_BGL"
refPars[20,] <- c(0.491, 0.2, 0.8) 
rownames(refPars)[20] <- "CSE_POC"
refPars[21,] <- c(0.460, 0.15, 0.8) 
rownames(refPars)[21] <- "CSE_HOC"
refPars[22,] <- c(0, 0, 1) 
rownames(refPars)[22] <- "MIN_SMsf"
refPars[23,] <- c(45.1, 40, 60) 
rownames(refPars)[23] <- "MAX_SMsf"
refPars[24,] <- c(240.71, 220, 260) 
rownames(refPars)[24] <- "beta0"
refPars[25,] <- c(66.02, 60, 70) 
rownames(refPars)[25] <- "beta1"
refPars[26,] <- c(227.13, 200, 250) 
rownames(refPars)[26] <- "beta2"
refPars[27,] <- c(30, 0, 40) 
rownames(refPars)[27] <- "MIN_SMsd"
refPars[28,] <- c(70, 40, 100) 
rownames(refPars)[28] <- "MAX_SMsd"
refPars[29,] <- c(0.2, 0.005, 0.3) 
rownames(refPars)[29] <- "h_SMsd"
refPars[30,] <- c(0.05, 0.001, 0.1)  
rownames(refPars)[30] <- "l_SMsd" 
refPars[31,] <- c(0.95, 0.6, 1) 
rownames(refPars)[31] <- "snsd"
refPars[32,] <- c(0.15, 0.008, 0.2) 
rownames(refPars)[32] <- "frt"
refPars[33,] <- c(40, 0, 50) 
rownames(refPars)[33] <- "MIN_SMrd"
refPars[34,] <- c(80, 50, 100) 
rownames(refPars)[34] <- "MAX_SMrd"
refPars[35,] <- c(275, 200, 280) 
rownames(refPars)[35] <- "rdt"
refPars[36,] <- c(0.07, 0.002, 0.1) 
rownames(refPars)[36] <- "h_SMrd"  
refPars[37,] <- c(0.02, 0.001, 0.03)  
rownames(refPars)[37] <- "l_SMrd"  
refPars[38,] <- c(3, 0.1, 5) 
rownames(refPars)[38] <- "error-sd"
refPars=refPars[12:38,]
refPars

# Define initial C pools
C_AGB = 200
C_BGB = 700
C_AGL = 200 
C_BGL = 700 
C_POC = 1000 
C_HOC = 5000 

# Model spinup
spin_year = 2002
len = length(unique(spin$site_num))
resultspinup = data.frame()
calcSOC <- function(param, num_sp, tsoil_sp, sm_sp, sm2_sp, clay_sp, GPP_sp){
  F_Ra = param[1]
  F_Root = param[2]
  k_AGL = param[3]
  k_BGL = param[4]
  k_POC = param[5]
  k_HOC = param[6]
  CSE_AGL = param[7]
  CSE_BGL = param[8]
  CSE_POC = param[9]
  CSE_HOC = param[10]
  MIN_SMsf = param[11]  
  MAX_SMsf = param[12] 
  beta0 = param[13]
  beta1 = param[14]
  beta2 = param[15]
  MIN_SMsd = param[16]
  MAX_SMsd = param[17]
  h_SMsd = param[18]
  l_SMsd = param[19]
  snsd = param[20]
  frt = param[21] 
  MIN_SMrd = param[22]
  MAX_SMrd = param[23]
  rdt = param[24]
  h_SMrd = param[25] 
  l_SMrd = param[26] 
  for (k in 1:len){
    tsoil_sp_k <- tsoil_sp[num_sp == k] 
    sm_sp_k <- sm_sp[num_sp == k] 
    sm2_sp_k <- sm2_sp[num_sp == k] 
    clay_sp_k <- clay_sp[num_sp == k] 
    GPP_sp_k <- GPP_sp[num_sp == k] 
    Tmult_sp = exp(beta0*((1/beta1)-(1/(tsoil_sp_k-beta2))))
    Wmult_sp = ifelse (sm_sp_k < MIN_SMsf, 0,
                       ifelse (sm_sp_k > MAX_SMsf, 1,
                               (sm_sp_k-(MIN_SMsf))/(MAX_SMsf-(MIN_SMsf))*1))
    Smult_sp = 0.7+0.3*(exp(-0.08*clay_sp_k))
    fsdethg_sp = ifelse (sm_sp_k < MIN_SMsd, h_SMsd,
                         ifelse (sm_sp_k > MAX_SMsd, l_SMsd,
                                 h_SMsd-(sm_sp_k-MIN_SMsd)/(MAX_SMsd-MIN_SMsd)*(h_SMsd-l_SMsd)))
    fsdeths_sp = snsd
    fallrt_sp = frt
    rdrg_sp = ifelse (sm2_sp_k < MIN_SMrd, h_SMrd,
                      ifelse (sm2_sp_k > MAX_SMrd, l_SMrd,
                              h_SMrd-(sm2_sp_k-MIN_SMrd)/(MAX_SMrd-MIN_SMrd)*(h_SMrd-l_SMrd))) 
    rdrs_sp = rdt
    NPP_sp = matrix(NA, nrow = 73, ncol = 1)
    AGB_sp = matrix(NA, nrow = 73, ncol = spin_year)
    BGB_sp = matrix(NA, nrow = 73, ncol = spin_year)
    ddt_C_AGL_sp = matrix(NA, nrow = 73, ncol = spin_year)
    ddt_C_BGL_sp = matrix(NA, nrow = 73, ncol = spin_year)
    ddt_C_POC_sp = matrix(NA, nrow = 73, ncol = spin_year)
    ddt_C_HOC_sp = matrix(NA, nrow = 73, ncol = spin_year)
    R_AGL_sp = matrix(NA, nrow = 73, ncol = spin_year)
    R_BGL_sp = matrix(NA, nrow = 73, ncol = spin_year)
    R_POC_sp = matrix(NA, nrow = 73, ncol = spin_year)
    R_HOC_sp = matrix(NA, nrow = 73, ncol = spin_year)
    Rh_sp = matrix(NA, nrow = 73, ncol = spin_year)
    Ra_sp = matrix(NA, nrow = 73, ncol = 1)
    NEE_sp = matrix(NA, nrow = 73, ncol = spin_year)
    for (i in 1:spin_year){
      for (j in 1: 73){
        days = ifelse(j < 73, 5, ifelse(i %% 4 == 0, 6, 5))
        fsdeth_sp = ifelse(j < 46, fsdethg_sp[j], fsdeths_sp) 
        rdr_sp = ifelse(tsoil_sp[j] > rdrs_sp, rdrg_sp[j], 0) 
        GPP_sp_k[j] = GPP_sp_k[j] * days
        NPP_sp[j] = GPP_sp_k[j] * (1-F_Ra)
        AGB_sp[j,i] = ifelse((i==1)&(j==1), C_AGB * (1 - fsdeth_sp*fallrt_sp) + NPP_sp[1] * (1.0 - F_Root), 
                             ifelse(j==1, AGB_sp[73,i-1] * (1 - fsdeth_sp*fallrt_sp) + NPP_sp[1] * (1.0 - F_Root), 
                                    AGB_sp[j-1,i] * (1 - fsdeth_sp*fallrt_sp) + NPP_sp[j] * (1.0 - F_Root))) 
        BGB_sp[j,i] = ifelse((i==1)&(j==1), C_BGB * (1 - rdr_sp) + NPP_sp[1] * F_Root, 
                             ifelse(j==1, BGB_sp[73,i-1] * (1-rdr_sp) + NPP_sp[1] * F_Root,
                                    BGB_sp[j-1,i] * (1 - rdr_sp) + NPP_sp[j] * F_Root))
        ddt_C_AGL_sp[j,i] = ifelse((i==1)&(j==1), AGB_sp[1, 1] * fsdeth_sp * fallrt_sp - (Tmult_sp[1] * Wmult_sp[1] * k_AGL * days * C_AGL),
                                   ifelse(j==1, AGB_sp[73, i-1] * fsdeth_sp * fallrt_sp - (Tmult_sp[j] * Wmult_sp[j] * k_AGL * days * C_AGL),
                                          AGB_sp[j-1,i] * fsdeth_sp * fallrt_sp - (Tmult_sp[j] * Wmult_sp[j] * k_AGL * days * C_AGL)))
        ddt_C_BGL_sp[j,i] = ifelse((i==1)&(j==1), BGB_sp[1, 1] * rdr_sp - (Tmult_sp[1] * Wmult_sp[1] * k_BGL * days * C_BGL), 
                                   ifelse(j==1, BGB_sp[73, i-1] * rdr_sp - (Tmult_sp[j] * Wmult_sp[j] * k_BGL * days * C_BGL),
                                          BGB_sp[j-1,i] * rdr_sp - (Tmult_sp[j] * Wmult_sp[j] * k_BGL * days * C_BGL)))
        ddt_C_POC_sp[j,i] = (Tmult_sp[j] * Wmult_sp[j] * k_AGL * days * CSE_AGL * C_AGL + Tmult_sp[j] * Wmult_sp[j] * k_BGL * days * CSE_BGL * C_BGL) - Tmult_sp[j] * Wmult_sp[j] * Smult_sp[j] * k_POC * days * C_POC
        ddt_C_HOC_sp[j,i] = Tmult_sp[j] * Wmult_sp[j] * Smult_sp[j] * k_POC * days * C_POC * CSE_POC - Tmult_sp[j] * Wmult_sp[j] * Smult_sp[j] * k_HOC * days * C_HOC * (1 - CSE_HOC)
        C_AGL = (ddt_C_AGL_sp[j,i]) + C_AGL
        C_BGL = (ddt_C_BGL_sp[j,i]) + C_BGL
        C_POC = (ddt_C_POC_sp[j,i]) + C_POC
        C_HOC = (ddt_C_HOC_sp[j,i]) + C_HOC
        TSOC = C_POC + C_HOC + 1000
        R_AGL_sp[j,i] = (Tmult_sp[j] * Wmult_sp[j] * k_AGL * days * (1-CSE_AGL)) * C_AGL
        R_BGL_sp[j,i] = (Tmult_sp[j] * Wmult_sp[j] * k_BGL * days * (1-CSE_BGL)) * C_BGL
        R_POC_sp[j,i] = (Tmult_sp[j] * Wmult_sp[j] * Smult_sp[j] * k_POC * days * (1-CSE_POC)) * C_POC
        R_HOC_sp[j,i] = (Tmult_sp[j] * Wmult_sp[j] * Smult_sp[j] * k_HOC * days * (1-CSE_HOC)) * C_HOC
        Rh_sp[j,i] = (-1)*(R_AGL_sp[j,i] + R_BGL_sp[j,i] + R_POC_sp[j,i] + R_HOC_sp[j,i]) 
        Ra_sp[j] = NPP_sp[j] - GPP_sp_k[j] 
        NEE_sp[j,i] = (GPP_sp_k[j] + (Rh_sp[j,i] + Ra_sp[j])) 
        NEE_sp[j,i] = NEE_sp[j,i]/days
        GPP_sp_k[j] = GPP_sp_k[j]/days
        Rh_sp[j,i] = Rh_sp[j,i]/days
        Ra_sp[j] = Ra_sp[j]/days
        NPP_sp[j] = NPP_sp[j]/days
      }
    }
    result_list <- data.frame(AGB_sp[73,2002], BGB_sp[73,2002], C_AGL, C_BGL, C_POC, C_HOC, TSOC) 
    resultspinup <- rbind(resultspinup, result_list)
  }
  return(data.frame(resultspinup[[1]],resultspinup[[2]],resultspinup[[3]],resultspinup[[4]],resultspinup[[5]],resultspinup[[6]],resultspinup[[7]]))
}
result_spinup = calcSOC(refPars$best, spin$site_num, spin$TS, spin$SWC1, spin$SWC2, spin$Clay, spin$GPP) 
colnames(result_spinup) <- c("C_AGB", "C_BGB", "C_AGL", "C_BGL", "C_POC", "C_HOC", "TOC")
C_Pool <- cbind(site,result_spinup)
colnames(C_Pool)[1:2] <- c("Site_num", "Site")
