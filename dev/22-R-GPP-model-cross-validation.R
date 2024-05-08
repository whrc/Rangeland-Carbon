library(dplyr)
library(BayesianTools)
library(robustbase)
setwd("C:/Users/username/path-to-the-calibration-data-folder")

# Import data
i_site=1 # Change number here for different sites
NEE <- read.table('file-with-fPAR-and-covariates-for-NEE-model.csv', header=T, sep = ",") 
NEE <- NEE %>% arrange(Site)
NEE$site_num = as.numeric(as.factor(with(NEE,paste(Site))))
site <- unique(NEE[c("site_num","Site")])
site <- site %>% arrange(site_num)
NEE$NEE_measured <- NEE$NEE_measured*24*3600/(10^6)*12 # Convert unit to g C/m2 if not done so
obs = NEE['NEE_measured']
obs <- na.omit(obs)
date <- NEE[which(!is.na(NEE$NEE_measured)),]['Date']
sitelist <- data.frame(NEE[(!is.na(NEE$NEE_measured)),]$Site)
NEE <- NEE[c("Site","site_num","Date","SW_IN_NLDAS","SWC1","SWC2","VPD","TS","Clay","STARFM_fPAR","NEE_measured")]

# Import model spin-up data
spin <- read.table('file-to-spin-up-dataset', header=T, sep = ",")
spin$Site <- sapply(strsplit(spin$Site, ".csv"),`[`, 1)
spin <- spin[!spin$Site %in% list,]
spin <- merge(spin,site,by= "Site")
spin <- spin[c("Site","site_num","DOY","SW_IN_NLDAS","SWC1","SWC2","VPD","TS","Clay","STARFM_fPAR")]

# Define GPP model parameters
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

# Get GPP estimates for the spin-up and NEE files
result_spin = calcGPP(refPars$best, spin$TS, spin$SWC2, spin$VPD, spin$SW_IN_NLDAS, spin$STARFM_fPAR) 
spin <- data.frame(spin,result_spin)
result_NEE = calcGPP(refPars$best, NEE$TS, NEE$SWC2, NEE$VPD, NEE$SW_IN_NLDAS, NEE$STARFM_fPAR) 
NEE <- data.frame(NEE,result_NEE)

# Change unit to match with NEE modeling
spin$SWC1 <- spin$SWC1 *100
spin$SWC2 <- spin$SWC2 *100
spin$TS <- spin$TS + 273.15
NEE$SWC1 <- NEE$SWC1 *100
NEE$SWC2 <- NEE$SWC2 *100
NEE$TS <- NEE$TS + 273.15

# Set parameters for spin-up and NEE model
refPars <- VSEMgetDefaults()
refPars[12,] <- c(0.447, 0.35, 0.55) 
rownames(refPars)[12] <- "F_Ra"
refPars[13,] <- c(0.629, 0.5, 0.72) 
rownames(refPars)[13] <- "F_Root"
refPars[14,] <- c(0.0184, 0.001, 0.07) 
rownames(refPars)[14] <- "k_AGL"
refPars[15,] <- c(0.0256, 0.0005, 0.05) 
rownames(refPars)[15] <- "k_BGL"
refPars[16,] <- c(0.0059, 0.0001, 0.015) #Multiply by 10 here to allow result exports, divide back in the next step
rownames(refPars)[16] <- "k_POC"
refPars[17,] <- c(0.0027, 0.00005, 0.008) #Multiply by 10 here to allow result exports, divide back in the next step
rownames(refPars)[17] <- "k_HOC"
refPars[18,] <- c(0.397, 0.05, 0.65) 
rownames(refPars)[18] <- "CSE_AGL"
refPars[19,] <- c(0.336, 0.1, 0.7) 
rownames(refPars)[19] <- "CSE_BGL"
refPars[20,] <- c(0.443, 0.2, 0.7)
rownames(refPars)[20] <- "CSE_POC"
refPars[21,] <- c(0.479, 0.2, 0.75)
rownames(refPars)[21] <- "CSE_HOC"
refPars[22,] <- c(10.45, 0, 20) 
rownames(refPars)[22] <- "MIN_SMsf"
refPars[23,] <- c(62.72, 35, 100) 
rownames(refPars)[23] <- "MAX_SMsf"
refPars[24,] <- c(242.48, 225, 265) 
rownames(refPars)[24] <- "beta0"
refPars[25,] <- c(67.87, 55, 70) 
rownames(refPars)[25] <- "beta1"
refPars[26,] <- c(227, 205, 250) 
rownames(refPars)[26] <- "beta2"
refPars[27,] <- c(21.68, 0, 38) 
rownames(refPars)[27] <- "MIN_SMsd"
refPars[28,] <- c(67.61, 38, 100) 
rownames(refPars)[28] <- "MAX_SMsd"
refPars[29,] <- c(0.14, 0.005, 0.32) 
rownames(refPars)[29] <- "h_SMsd"
refPars[30,] <- c(0.0414, 0.005, 0.1) 
rownames(refPars)[30] <- "l_SMsd"
refPars[31,] <- c(0.74, 0.6, 0.99) 
rownames(refPars)[31] <- "snsd"
refPars[32,] <- c(0.146, 0.02, 0.25)
rownames(refPars)[32] <- "frt"
refPars[33,] <- c(29.2, 0, 50) 
rownames(refPars)[33] <- "MIN_SMrd"
refPars[34,] <- c(79.16, 50, 110) 
rownames(refPars)[34] <- "MAX_SMrd"
refPars[35,] <- c(260.8, 230, 290) 
rownames(refPars)[35] <- "rdt"
refPars[36,] <- c(0.048, 0.005, 0.15) 
rownames(refPars)[36] <- "h_SMrd"
refPars[37,] <- c(0.0169, 0.002, 0.03) 
rownames(refPars)[37] <- "l_SMrd"
refPars[38,] <- c(3, 0.1, 5)
rownames(refPars)[38] <- "error-sd"
refPars=refPars[12:38,]
refPars

# Initial C pool setup
C_Pool <- read.table('file-to-initial-C-pool-values-stimated-with-default-parameters.csv', header=T, sep = ",")
C_Pool <- C_Pool[,c("Site","C_AGB","C_BGB","C_AGL","C_BGL","C_POC","C_HOC")]
C_Pool <- merge(C_Pool, site, by = "Site")

# Split sites into calibration and validation
NEE_c <- NEE[which(NEE$site_num != i_site),]
NEE_i <- NEE[which(NEE$site_num == i_site),]
spin_c <- spin[which(spin$site_num != i_site),]
spin_i <- spin[which(spin$site_num == i_site),]
obs_c = na.omit(NEE_c['NEE_measured'])
obs_i = na.omit(NEE_i['NEE_measured'])
site_num_list <- unique(NEE_c$site_num)

# Spinup and NEE model
spin_year = 2002
len = length(unique(NEE_c$site_num))
resultspinup = data.frame()
resultnee = data.frame()
calcSOC <- function(param, num_sp, tsoil_sp, sm_sp, sm2_sp, clay_sp, GPP_sp, 
                    tsoil_in, sm_in, sm2_in, clay_in, GPP_in, num_in, NEE_m){
  len = ifelse(length(unique(num_sp)) == 1,1,len)
  F_Ra = param[1]
  F_Root = param[2]
  k_AGL = param[3]
  k_BGL = param[4]
  k_POC = param[5]/10 
  k_HOC = param[6]/10 
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
    k = ifelse(length(unique(num_sp)) == 1,unique(num_sp),site_num_list[k]) 
    C_AGB = C_Pool[which(C_Pool$site_num == k),c("C_AGB")]
    C_BGB = C_Pool[which(C_Pool$site_num == k),c("C_BGB")]
    C_AGL = C_Pool[which(C_Pool$site_num == k),c("C_AGL")] 
    C_BGL = C_Pool[which(C_Pool$site_num == k),c("C_BGL")]
    C_POC = C_Pool[which(C_Pool$site_num == k),c("C_POC")] 
    C_HOC = C_Pool[which(C_Pool$site_num == k),c("C_HOC")]
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
        days = 5
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
        R_AGL_sp[j,i] = (Tmult_sp[j] * Wmult_sp[j] * k_AGL * days * (1-CSE_AGL)) * C_AGL
        R_BGL_sp[j,i] = (Tmult_sp[j] * Wmult_sp[j] * k_BGL * days * (1-CSE_BGL)) * C_BGL
        R_POC_sp[j,i] = (Tmult_sp[j] * Wmult_sp[j] * Smult_sp[j] * k_POC * days * (1-CSE_POC)) * C_POC
        R_HOC_sp[j,i] = (Tmult_sp[j] * Wmult_sp[j] * Smult_sp[j] * k_HOC * days * (1-CSE_HOC)) * C_HOC
        Rh_sp[j,i] = (-1)*(R_AGL_sp[j,i] + R_BGL_sp[j,i] + R_POC_sp[j,i] + R_HOC_sp[j,i]) 
        Ra_sp[j] = NPP_sp[j] - GPP_sp[j] 
        NEE_sp[j,i] = (GPP_sp[j] + (Rh_sp[j,i] + Ra_sp[j])) 
        NEE_sp[j,i] = NEE_sp[j,i]/days
        GPP_sp_k[j] = GPP_sp_k[j]/days
        Rh_sp[j,i] = Rh_sp[j,i]/days
        Ra_sp[j] = Ra_sp[j]/days
        NPP_sp[j] = NPP_sp[j]/days
      }
    }
    result_list <- data.frame(AGB_sp[73,2002], BGB_sp[73,2002],C_AGL, C_BGL, C_POC, C_HOC)
    C_AGB_k = result_list[1,1]
    C_BGB_k = result_list[1,2]
    C_AGL_k = result_list[1,3]
    C_BGL_k = result_list[1,4]
    C_POC_k = result_list[1,5]
    C_HOC_k = result_list[1,6]
    tsoil_in_k <- tsoil_in[num_in == k] 
    sm_in_k <- sm_in[num_in == k] 
    sm2_in_k <- sm2_in[num_in == k] 
    clay_in_k <- clay_in[num_in == k] 
    GPP_in_k <- GPP_in[num_in == k] 
    NEE_m_k = NEE_m[num_in == k]
    Tmult_in = exp(beta0*((1/beta1)-(1/(tsoil_in_k-beta2))))
    Wmult_in = ifelse (sm_in_k < MIN_SMsf, 0,
                       ifelse (sm_in_k > MAX_SMsf, 1,
                               (sm_in_k-(MIN_SMsf))/(MAX_SMsf-(MIN_SMsf))*1))
    Smult_in = 0.7+0.3*(exp(-0.08*clay_in_k))
    fsdethg_in = ifelse (sm_in_k < MIN_SMsd, h_SMsd,
                         ifelse (sm_in_k > MAX_SMsd, l_SMsd,
                                 h_SMsd-(sm_in_k-MIN_SMsd)/(MAX_SMsd-MIN_SMsd)*(h_SMsd-l_SMsd)))
    fsdeths_in = snsd
    fallrt_in = frt
    rdrg_in = ifelse (sm2_in_k < MIN_SMrd, h_SMrd,
                      ifelse (sm2_in_k > MAX_SMrd, l_SMrd,
                              h_SMrd-(sm2_in_k-MIN_SMrd)/(MAX_SMrd-MIN_SMrd)*(h_SMrd-l_SMrd))) 
    rdrs_in = rdt
    NPP_in = matrix(NA, nrow = 1412, ncol = 1)
    AGB_in = matrix(NA, nrow = 1412, ncol = 1)
    BGB_in = matrix(NA, nrow = 1412, ncol = 1)
    ddt_C_AGL_in = matrix(NA, nrow = 1412, ncol = 1)
    ddt_C_BGL_in = matrix(NA, nrow = 1412, ncol = 1)
    ddt_C_POC_in = matrix(NA, nrow = 1412, ncol = 1)
    ddt_C_HOC_in = matrix(NA, nrow = 1412, ncol = 1)
    TSOC_in = matrix(NA, nrow = 1412, ncol = 1)
    R_AGL_in = matrix(NA, nrow = 1412, ncol = 1)
    R_BGL_in = matrix(NA, nrow = 1412, ncol = 1)
    R_POC_in = matrix(NA, nrow = 1412, ncol = 1)
    R_HOC_in = matrix(NA, nrow = 1412, ncol = 1)
    Rh_in = matrix(NA, nrow = 1412, ncol = 1)
    Ra_in = matrix(NA, nrow = 1412, ncol = 1)
    NEE_in = matrix(NA, nrow = 1412, ncol = 1)
    for (n in 1: 1412){
      days = 5
      fsdeth_in = ifelse((n %% 73 < 46)&(n %% 73 > 0), fsdethg_sp[n %% 73], fsdeths_sp) 
      rdr_in = ifelse(tsoil_in[n] > rdrs_sp, rdrg_in[n], 0) 
      GPP_in_k[n] = GPP_in_k[n] * days
      NPP_in[n] = GPP_in_k[n] * (1-F_Ra)
      AGB_in[n] = ifelse(n == 1, C_AGB_k * (1 - fsdeth_in*fallrt_in) + NPP_in[1] * (1.0 - F_Root), 
                         AGB_in[n-1] * (1 - fsdeth_in*fallrt_in) + NPP_in[n] * (1.0 - F_Root)) 
      BGB_in[n] = ifelse(n == 1, C_BGB_k * (1 - rdr_in) + NPP_in[1] * F_Root, 
                         BGB_in[n-1] * (1 - rdr_in) + NPP_in[n] * F_Root)
      ddt_C_AGL_in[n] = ifelse(n == 1, AGB_in[1] * fsdeth_in * fallrt_in - (Tmult_in[1] * Wmult_in[1] * k_AGL * days * C_AGL_k),
                               AGB_in[n-1] * fsdeth_in * fallrt_in - (Tmult_in[n] * Wmult_in[n] * k_AGL * days * C_AGL_k))
      ddt_C_BGL_in[n] = ifelse(n == 1, BGB_in[1] * rdr_in - (Tmult_in[1] * Wmult_in[1] * k_BGL * days * C_BGL_k), 
                               BGB_in[n-1] * rdr_in - (Tmult_in[n] * Wmult_in[n] * k_BGL * days * C_BGL_k))
      ddt_C_POC_in[n] = (Tmult_in[n] * Wmult_in[n] * k_AGL * days * CSE_AGL * C_AGL_k + Tmult_in[n] * Wmult_in[n] * k_BGL * days * CSE_BGL * C_BGL_k) - Tmult_in[n] * Wmult_in[n] * Smult_in[n] * k_POC * days * C_POC_k
      ddt_C_HOC_in[n] = Tmult_in[n] * Wmult_in[n] * Smult_in[n] * k_POC * days * C_POC_k * CSE_POC - Tmult_in[n] * Wmult_in[n] * Smult_in[n] * k_HOC * days * C_HOC_k * (1 - CSE_HOC)
      C_AGL_k = (ddt_C_AGL_in[n]) + C_AGL_k
      C_BGL_k = (ddt_C_BGL_in[n]) + C_BGL_k
      C_POC_k = (ddt_C_POC_in[n]) + C_POC_k
      C_HOC_k = (ddt_C_HOC_in[n]) + C_HOC_k
      TSOC_in[n] = C_POC_k + C_HOC_k + 1000
      R_AGL_in[n] = (Tmult_in[n] * Wmult_in[n] * k_AGL * days * (1-CSE_AGL)) * C_AGL_k
      R_BGL_in[n] = (Tmult_in[n] * Wmult_in[n] * k_BGL * days * (1-CSE_BGL)) * C_BGL_k
      R_POC_in[n] = (Tmult_in[n] * Wmult_in[n] * Smult_in[n] * k_POC * days * (1-CSE_POC)) * C_POC_k
      R_HOC_in[n] = (Tmult_in[n] * Wmult_in[n] * Smult_in[n] * k_HOC * days * (1-CSE_HOC)) * C_HOC_k
      Rh_in[n] = (-1)*(R_AGL_in[n] + R_BGL_in[n] + R_POC_in[n] + R_HOC_in[n]) 
      Ra_in[n] = NPP_in[n] - GPP_in_k[n] 
      NEE_in[n] = (GPP_in_k[n] + (Rh_in[n] + Ra_in[n])) 
      NEE_in[n] = NEE_in[n]/days
      GPP_in_k[n] = GPP_in_k[n]/days
      Rh_in[n] = Rh_in[n]/days
      Ra_in[n] = Ra_in[n]/days
      NPP_in[n] = NPP_in[n]/days
    }
    NEE_p <- NEE_in[!is.na(NEE_m_k)] 
    out_list <- data.frame(NEE_p)
    resultnee <- rbind(resultnee, out_list)
  }
  return(resultnee)
}
result_before <- calcSOC(refPars$best, spin_c$site_num, spin_c$TS, spin_c$SWC1, spin_c$SWC2, spin_c$Clay, spin_c$GPP, 
                         NEE_c$TS, NEE_c$SWC1, NEE_c$SWC2, NEE_c$Clay, NEE_c$GPP, NEE_c$site_num, NEE_c$NEE_measured) 

# Set parameters
parSel = c(1:26) # All parameters
prior <- createUniformPrior(lower = refPars$lower[parSel], 
                            upper = refPars$upper[parSel], best = refPars$best[parSel])
settings <- list(iterations = 5000, nrChains = 3) 

# Carry out calibration
likelihood <- function(par, sum = TRUE){
  x = refPars$best
  x[parSel] = par
  predicted <- calcSOC(refPars$best, spin_c$site_num, spin_c$TS, spin_c$SWC1, spin_c$SWC2, spin_c$Clay, spin_c$GPP, 
                       NEE_c$TS, NEE_c$SWC1, NEE_c$SWC2, NEE_c$Clay, NEE_c$GPP, NEE_c$site_num, NEE_c$NEE_measured) 
  diff <- c(predicted[,1] - obs_c[,1])
  llValues <- dnorm(diff, sd = x[27], log = TRUE)  
  if (sum == FALSE) return(llValues)
  else return(sum(llValues))
}
bayesianSetup <- createBayesianSetup(likelihood, prior, names = rownames(refPars)[parSel])
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
post <- data.frame(getSample(out))
post_mean <- colMedians(as.matrix(post))
newresult = calcSOC(post_mean, spin_i$site_num, spin_i$TS, spin_i$SWC1, spin_i$SWC2, spin_i$Clay, spin_i$GPP, 
                    NEE_i$TS, NEE_i$SWC1, NEE_i$SWC2, NEE_i$Clay, NEE_i$GPP, NEE_i$site_num, NEE_i$NEE_measured)

# Interpret results
comparison_before <- data.frame(obs_c,result_before)
comparison_before$Site <- NEE_c[which(!is.na(NEE_c$NEE_measured)),]$Site
comparison_before$Date <- NEE_c[which(!is.na(NEE_c$NEE_measured)),]$Date
RMSE <- sqrt(mean((comparison_before$NEE_measured - comparison_before$NEE_p)^2))
R2 <- cor(comparison_before$NEE_measured, comparison_before$NEE_p)^2
MBE <- mean(comparison_before$NEE_measured - comparison_before$NEE_p)
metric_before <- data.frame(RMSE,R2,MBE)
comparison_after <- data.frame(obs_i,newresult)
comparison_after$Site <- unique(NEE_i$Site)
comparison_after$Date <- NEE_i[which(!is.na(NEE_i$NEE_measured)),]$Date
rmse <- sqrt(mean((comparison_after$NEE_measured - comparison_after$NEE_p)^2))
r2 <- cor(comparison_after$NEE_measured, comparison_after$NEE_p)^2
mbe <- mean(comparison_after$NEE_measured - comparison_after$NEE_p)
metric_after <- data.frame(rmse,r2,mbe)
model_params <- data.frame(post_mean)
