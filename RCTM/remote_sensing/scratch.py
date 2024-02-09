#Spinup
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
                                          BGB_sp[j-1,i] * rdr_sp) - (Tmult_sp[j] * Wmult_sp[j] * k_BGL * days * C_BGL))
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
write.csv(C_Pool,"Spinup/C_Pool_STARFM.csv")
