#### Header ####################################################################
## name:             ACOT4_paper_models_and_figures_02
## 
## purpose:          Neatly, succinctly, and reproducibly run analysis of ACOT4
##                   variants in NIHS and UKBB as presented in the paper.
##                   
##                   
## copyright (C) 2024 Jacob William Ilijev Meyjes-Brown
## written by: JWIMB
##
## first written: 14 October 2024 (based on older analyses)
## last modified: 11 November 2024
################################################################################

# 1.0 set up ###################################################################
library(tidyverse)
library(magrittr)
library(lm.beta)
library(ggpubr)
library(sjPlot)
library(readxl)

setwd(dir = "~/ACOT4/")

#  1.1 import data -------------------------------------------------------------
NIHS_2000 <- 
  read_excel("DATA/NI/Working_Data/NIHS2000_phenotype_JMB_working_copy_03.xlsx",
             na = "NA")


UKBB_WES  <- 
  read_csv("DATA/ukbb.wes200k.acot4/UKBB_ACOT4_WES_200k_composite.csv")

#  1.2 create carrier variables ------------------------------------------------
NIHS_2000 <- 
  mutate(NIHS_2000,
         MNV_carrier = case_when(
           ACOT4_simple_geno >= 1 ~ 1,
           ACOT4_simple_geno == 0 ~ 0
         )
  )

NIHS_2000 <- 
  mutate(NIHS_2000,
         WT_carrier = case_when(
           ACOT4_simple_geno <= 1 ~ 1,
           ACOT4_simple_geno == 2 ~ 0
         )
  )

UKBB_WES <- 
  mutate(UKBB_WES,
         MNV_carrier = case_when(
           rs35724886 >= 1 ~ 1,
           rs35724886 == 0 ~ 0
         )
  )

UKBB_WES <- 
  mutate(UKBB_WES,
         WT_carrier = case_when(
           rs35724886 <= 1 ~ 1,
           rs35724886 == 2 ~ 0
         )
  )


#  1.3 functions ---------------------------------------------------------------

readable_logit <- function(x){   # x = a logistic regression model
    #-----------------------------------------#
  SUMMARY_x <- summary(x)
  CONFINT_x <- exp(confint(x))
  CALL_x    <- SUMMARY_x[1]
  TABLE_x   <- coef(SUMMARY_x)
  COEF_x    <- coef(x)
  OR_x      <- exp(COEF_x)
    #-----------------------------------------#
  READABLE  <- round(
    cbind(
      TABLE_x[, 1:2], 
  OR_x, 
  CONFINT_x, 
  TABLE_x[, 3:4]
    ), 
  digits = 4
  )
    #-----------------------------------------#
  print(CALL_x)
  print(READABLE)
    #-----------------------------------------#
}


# 2.0 models ###################################################################
NIHS_2000 %>% 
  select(c( SBP, DBP, EBP,CHOL, high_chol, HDL, lowHDL)) %>%
  summary()  
  
NIHS_2000 %>% 
  select(c( SBP, DBP, EBP,CHOL, high_chol, HDL, lowHDL)) %>%
  lapply(sd, na.rm = T)

NIHS_2000 %>% 
  filter(Gender == "female") %>% 
  select(c( SBP, DBP, EBP,CHOL, high_chol, HDL, lowHDL)) %>%
  summary()  
  
NIHS_2000 %>% 
  filter(Gender == "female") %>% 
  select(c( SBP, DBP, EBP,CHOL, high_chol, HDL, lowHDL)) %>%
  lapply(sd, na.rm = T)

NIHS_2000 %>% 
  filter(Gender == "male") %>% 
  select(c( SBP, DBP, EBP,CHOL, high_chol, HDL, lowHDL)) %>%
  summary()  
  
NIHS_2000 %>% 
  filter(Gender == "male") %>% 
  select(c( SBP, DBP, EBP,CHOL, high_chol, HDL, lowHDL)) %>%
  lapply(sd, na.rm = T)

NIHS_2000 %>% 
  filter(Gender == "male") %>% 
  select(c( SBP, DBP, EBP,CHOL, high_chol, HDL, lowHDL)) %>%
  lapply(sum, na.rm = T)
#------------------------------------------------------------------------------#
UKBB_WES %>% 
  select(c(av_sbp, av_dbp, diagnosis_BP_ext, 
           CHOL, HiChol_or_meds, HDL, diagnosis_HDL_ext)) %>%
  summary()  
  
UKBB_WES %>% 
  select(c(av_sbp, av_dbp, diagnosis_BP_ext, 
           CHOL, HiChol_or_meds, HDL, diagnosis_HDL_ext)) %>%
  lapply(sd, na.rm = T)
  
UKBB_WES %>% 
  filter(Sex == "female") %>% 
  select(c(av_sbp, av_dbp, diagnosis_BP_ext, 
           CHOL, HiChol_or_meds, HDL, diagnosis_HDL_ext)) %>%
  summary()  
  
UKBB_WES %>% 
  filter(Sex == "female") %>% 
  select(c(av_sbp, av_dbp, diagnosis_BP_ext, 
           CHOL, HiChol_or_meds, HDL, diagnosis_HDL_ext)) %>%
  lapply(sd, na.rm = T)
  
  
UKBB_WES %>% 
  filter(Sex == "male") %>% 
  select(c(av_sbp, av_dbp, diagnosis_BP_ext, 
           CHOL, HiChol_or_meds, HDL, diagnosis_HDL_ext)) %>%
  summary()  
  
UKBB_WES %>% 
  filter(Sex == "male") %>% 
  select(c(av_sbp, av_dbp, diagnosis_BP_ext, 
           CHOL, HiChol_or_meds, HDL, diagnosis_HDL_ext)) %>%
  lapply(sd, na.rm = T)
  
# 2.1 blood pressure -----------------------------------------------------------
#  2.1.1 NIHS systolic ---------------------------------------------------------

NIHS_model_spb_lin_a <- 
  lm(data = NIHS_2000,
     formula = SBP ~ 
  ACOT4_simple_geno +
  AGE_2000NIHS + 
  Gender + 
  COREPED_1 + 
  CHOL
  )
summary(NIHS_model_spb_lin_a)
summary(lm.beta(NIHS_model_spb_lin_a))  
confint(lm.beta(NIHS_model_spb_lin_a))  

# p = 0.13

NIHS_model_spb_lin_g <- 
  lm(data = NIHS_2000,
     formula = SBP ~ 
  as.factor(ACOT4_simple_geno) +
  AGE_2000NIHS + 
  Gender + 
  COREPED_1 + 
  CHOL
  )
summary(NIHS_model_spb_lin_g)
summary(lm.beta(NIHS_model_spb_lin_g))
confint(lm.beta(NIHS_model_spb_lin_g))

# relevant output -------------------------------------------------------------#
# Coefficients:
#                               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   70.13539    4.60819  15.220  < 2e-16 ***
# as.factor(ACOT4_simple_geno)1 -1.16317    1.75376  -0.663  0.50744    
# as.factor(ACOT4_simple_geno)2 -7.21139    4.09495  -1.761  0.07876 .  
# ...
# relevant output end ---------------------------------------------------------#

NIHS_model_spb_lin_mnv <- 
  lm(data = NIHS_2000,
     formula = SBP ~ 
  MNV_carrier +
  AGE_2000NIHS + 
  Gender + 
  COREPED_1 + 
  CHOL
  )
summary(NIHS_model_spb_lin_mnv)
summary(lm.beta(NIHS_model_spb_lin_mnv))
confint(lm.beta(NIHS_model_spb_lin_mnv)) %>% round(3)
# p = 0.27


NIHS_model_spb_lin_wt <- 
  lm(data = NIHS_2000,
     formula = SBP ~ 
  WT_carrier +
  AGE_2000NIHS + 
  Gender + 
  COREPED_1 + 
  CHOL
  )
summary(NIHS_model_spb_lin_wt)
summary(lm.beta(NIHS_model_spb_lin_wt))
confint(lm.beta(NIHS_model_spb_lin_wt)) %>% round(3)

# p = 0.09

#  2.1.2 NIHS diastolic --------------------------------------------------------

NIHS_model_dpb_lin_a <- 
  lm(data = NIHS_2000,
     formula = DBP ~
  ACOT4_simple_geno +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  CHOL
  )
summary(NIHS_model_dpb_lin_a)
# p = 0.45

NIHS_model_dpb_lin_g <- 
  lm(data = NIHS_2000,
     formula = DBP ~
  as.factor(ACOT4_simple_geno) +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  CHOL
  )
summary(NIHS_model_dpb_lin_g)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                   51.92013    2.85277  18.200  < 2e-16 ***
# as.factor(ACOT4_simple_geno)1  0.62964    1.08569   0.580   0.5622    
# as.factor(ACOT4_simple_geno)2 -5.08604    2.53504  -2.006   0.0453 *  
# ...
# relevant output end ---------------------------------------------------------#

NIHS_model_dpb_lin_mnv <- 
  lm(data = NIHS_2000,
     formula = DBP ~ 
  MNV_carrier +
  AGE_2000NIHS + 
  Gender + 
  COREPED_1 + 
  CHOL
  )
summary(NIHS_model_dpb_lin_mnv)
# p = 0.99


NIHS_model_dpb_lin_wt <- 
  lm(data = NIHS_2000,
     formula = DBP ~ 
  WT_carrier +
  AGE_2000NIHS + 
  Gender + 
  COREPED_1 + 
  CHOL
  )
summary(NIHS_model_dpb_lin_wt)
summary(lm.beta(NIHS_model_dpb_lin_wt))
confint(lm.beta(NIHS_model_dpb_lin_wt))
# p = 0.034


#  2.1.3 NIHS elevated BP ------------------------------------------------------
NIHS_2000 <- 
  mutate(NIHS_2000,
         EBP = case_when(
           SBP >= 130 | DBP >= 85 ~ 1,
           SBP <  130 & DBP <  85 ~ 0),
           .after = DBP
           )
summary(NIHS_2000$EBP)
summary(as.factor(NIHS_2000$EBP))

NIHS_2000 <- 
  mutate(NIHS_2000,
         Measured_hypertension = case_when(
           SBP >= 140 | DBP >= 90 ~ 1,
           SBP <  140 & DBP <  90 ~ 0
         ), .after = EBP)



NIHS_2000 %>% select(Measured_hypertension) %>% summary()
NIHS_2000 %>% filter(AGE_2000NIHS >= 65) %>% select(Measured_hypertension) %>% summary()
NIHS_2000 %>% filter(AGE_2000NIHS >= 65) %>% select(EBP) %>% summary()

NIHS_2000_EBP_TABLE <- 
  table(NIHS_2000$EBP,
        NIHS_2000$ACOT4_simple_geno)

NIHS_2000_EBP_TABLE
#     0   1   2
# 0 190 105  21
# 1 169  93   4

NIHS_2000_HTN_TABLE <- 
  table(NIHS_2000$Measured_hypertension,
        NIHS_2000$ACOT4_simple_geno)

NIHS_2000_HTN_TABLE
#     0   1   2
# 0 241 133  22
# 1 118  65   3

round(
  prop.table(
    NIHS_2000_EBP_TABLE, margin = 1), 
digits = 3)
#       0     1     2
# 0 0.601 0.332 0.066
# 1 0.635 0.350 0.015

round(
  prop.table(
    NIHS_2000_HTN_TABLE, margin = 1), 
digits = 3)
#       0     1     2
# 0 0.609 0.336 0.056
# 1 0.634 0.349 0.016

fisher.test(NIHS_2000_EBP_TABLE)
# data:  NIHS_2000_EBP_TABLE
# p-value = 0.00766
# alternative hypothesis: two.side

fisher.test(NIHS_2000_HTN_TABLE)
# data:  NIHS_2000_HTN_TABLE
# p-value = 0.08059
# alternative hypothesis: two.sided

NIHS_model_epb_log_a <- 
  glm(data = NIHS_2000,
      formula = EBP ~
  ACOT4_simple_geno +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  CHOL,
family = binomial
  )
summary(NIHS_model_epb_log_a)
readable_logit(NIHS_model_epb_log_a)
# p = 0.08

NIHS_model_epb_log_g <- 
  glm(data = NIHS_2000,
      formula = EBP ~
  as.factor(ACOT4_simple_geno) +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  CHOL,
family = binomial
  )
summary(NIHS_model_epb_log_g)
readable_logit(NIHS_model_epb_log_g)
# relevant output -------------------------------------------------------------#
#                      Estimate Std. Error   OR_x  2.5 % 97.5 % z value Pr(>|z|)
# (Intercept)           -5.0871     0.6145 0.0062 0.0018 0.0199 -8.2780   0.0000
# as.factor(ACOT4...)1  -0.0817     0.2021 0.9216 0.6194 1.3691 -0.4043   0.6860
# as.factor(ACOT4...)2  -1.3751     0.5848 0.2528 0.0696 0.7291 -2.3512   0.0187
# ...
# relevant output end ---------------------------------------------------------#

NIHS_model_ebp_log_mnv <- 
  glm(data = NIHS_2000,
      formula = EBP ~
        MNV_carrier +
        AGE_2000NIHS +
        Gender +
        COREPED_1 +
        CHOL,
      family = binomial
      )
summary(NIHS_model_ebp_log_mnv)
readable_logit(NIHS_model_ebp_log_mnv)

NIHS_model_ebp_log_wt <- 
  glm(data = NIHS_2000,
      formula = EBP ~
        WT_carrier +
        AGE_2000NIHS +
        Gender +
        COREPED_1 +
        CHOL,
      family = binomial
      )
summary(NIHS_model_ebp_log_wt)
readable_logit(NIHS_model_ebp_log_wt)



#  2.1.4 UKBB systolic ---------------------------------------------------------
UKBB_model_sbp_lin_a_CHOL <- 
  lm(data = UKBB_WES,
     formula = av_sbp ~ 
  rs35724886 +
  Age +
  Sex +
  CHOL
  )
summary(UKBB_model_sbp_lin_a_CHOL)
summary(lm.beta(UKBB_model_sbp_lin_a_CHOL))
confint(lm.beta(UKBB_model_sbp_lin_a_CHOL)) %>% round(3)
# p = 0.0041

UKBB_model_sbp_lin_a_HDL <- 
  lm(data = UKBB_WES,
     formula = av_sbp ~ 
  rs35724886 +
  Age +
  Sex +
  HDL
  )
summary(UKBB_model_sbp_lin_a_HDL)
# p = 0.01

UKBB_model_sbp_lin_a_lowHDL <- 
  lm(data = UKBB_WES,
     formula = av_sbp ~ 
  rs35724886 +
  Age +
  Sex +
  diagnosis_HDL_ext
  )
summary(UKBB_model_sbp_lin_a_lowHDL)
# p = 0.02

#------------------------------------------------------------------------------#

UKBB_model_sbp_lin_g_CHOL <- 
  lm(data = UKBB_WES,
     formula = av_sbp ~ 
  as.factor(rs35724886) +
  Age +
  Sex +
  CHOL
  )
summary(UKBB_model_sbp_lin_g_CHOL)
summary(lm.beta(UKBB_model_sbp_lin_g_CHOL))
confint(lm.beta(UKBB_model_sbp_lin_g_CHOL)) %>% round(3)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            83.439834   0.337109 247.516   <2e-16 ***
# as.factor(rs35724886)1 -0.185058   0.085068  -2.175   0.0296 *  
# as.factor(rs35724886)2 -0.420376   0.191644  -2.194   0.0283 *  
# ...
# relevant output end ---------------------------------------------------------#

UKBB_model_sbp_lin_g_HDL <- 
  lm(data = UKBB_WES,
     formula = av_sbp ~ 
  as.factor(rs35724886) +
  Age +
  Sex +
  HDL
  )
summary(UKBB_model_sbp_lin_g_HDL)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            92.411198   0.339876 271.897  < 2e-16 ***
# as.factor(rs35724886)1 -0.158080   0.089286  -1.770   0.0766 .  
# as.factor(rs35724886)2 -0.414965   0.200841  -2.066   0.0388 *  
# ...
# relevant output end ---------------------------------------------------------#

UKBB_model_sbp_lin_g_lowHDL <- 
  lm(data = UKBB_WES,
     formula = av_sbp ~ 
  as.factor(rs35724886) +
  Age +
  Sex +
  diagnosis_HDL_ext
  )
summary(UKBB_model_sbp_lin_g_lowHDL)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            93.239957   0.275236 338.764   <2e-16 ***
# as.factor(rs35724886)1 -0.153524   0.083548  -1.838   0.0661 .  
# as.factor(rs35724886)2 -0.319388   0.188047  -1.698   0.0894 .  
# ...
# relevant output end ---------------------------------------------------------#


#------------------------------------------------------------------------------#

UKBB_model_sbp_lin_mnv_CHOL <- 
  lm(data = UKBB_WES,
     formula = av_sbp ~ 
       MNV_carrier +
       Age +
       Sex +
       CHOL
  )
summary(UKBB_model_sbp_lin_mnv_CHOL)
summary(lm.beta(UKBB_model_sbp_lin_mnv_CHOL))
confint(lm.beta(UKBB_model_sbp_lin_mnv_CHOL)) %>% round(3)

UKBB_model_sbp_lin_wt_CHOL <- 
  lm(data = UKBB_WES,
     formula = av_sbp ~ 
       WT_carrier +
       Age +
       Sex +
       CHOL
  )
summary(UKBB_model_sbp_lin_wt_CHOL)
summary(lm.beta(UKBB_model_sbp_lin_wt_CHOL))
confint(lm.beta(UKBB_model_sbp_lin_wt_CHOL)) %>% round(3)


#  2.1.4 UKBB diastolic --------------------------------------------------------
UKBB_model_dbp_lin_a_CHOL <- 
  lm(data = UKBB_WES,
     formula = av_dbp ~ 
  rs35724886 +
  Age +
  Sex +
  CHOL
  )
summary(UKBB_model_dbp_lin_a_CHOL)
# p = 5.53e-05

UKBB_model_dbp_lin_a_HDL <- 
  lm(data = UKBB_WES,
     formula = av_dbp ~ 
  rs35724886 +
  Age +
  Sex +
  HDL
  )
summary(UKBB_model_dbp_lin_a_HDL)
# p = 0.003

UKBB_model_dbp_lin_a_lowHDL <- 
  lm(data = UKBB_WES,
     formula = av_dbp ~ 
  rs35724886 +
  Age +
  Sex +
  diagnosis_HDL_ext
  )
summary(UKBB_model_dbp_lin_a_lowHDL)
# p = 0.0005

UKBB_model_dbp_lin_g_CHOL <- 
  lm(data = UKBB_WES,
     formula = av_dbp ~ 
  as.factor(rs35724886) +
  Age +
  Sex +
  CHOL
  )
summary(UKBB_model_dbp_lin_g_CHOL)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            70.386570   0.193047 364.608  < 2e-16 ***
# as.factor(rs35724886)1 -0.168384   0.048715  -3.457 0.000547 ***
# as.factor(rs35724886)2 -0.286237   0.109746  -2.608 0.009103 ** 
# ...
# relevant output end ---------------------------------------------------------#

UKBB_model_dbp_lin_g_HDL <- 
  lm(data = UKBB_WES,
     formula = av_dbp ~ 
  as.factor(rs35724886) +
  Age +
  Sex +
  HDL
  )
summary(UKBB_model_dbp_lin_g_HDL)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            79.194085   0.195708 404.655   <2e-16 ***
# as.factor(rs35724886)1 -0.127445   0.051413  -2.479   0.0132 *  
# as.factor(rs35724886)2 -0.225455   0.115648  -1.949   0.0512 .  
# ...
# relevant output end ---------------------------------------------------------#

UKBB_model_dbp_lin_g_lowHDL <- 
  lm(data = UKBB_WES,
     formula = av_dbp ~ 
  as.factor(rs35724886) +
  Age +
  Sex +
  diagnosis_HDL_ext
  )
summary(UKBB_model_dbp_lin_g_lowHDL)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            77.825380   0.158421 491.256  < 2e-16 ***
# as.factor(rs35724886)1 -0.145168   0.048089  -3.019  0.00254 ** 
# as.factor(rs35724886)2 -0.239560   0.108236  -2.213  0.02688 *  
# ...
# relevant output end ---------------------------------------------------------#

#  2.1.4 UKBB elevated BP ------------------------------------------------------
UKBB_model_ebp_s_log_a_CHOL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP ~
  rs35724886 +
  Age +
  Sex +
  CHOL,
family = binomial
  )
summary(UKBB_model_ebp_s_log_a_CHOL)
readable_logit(UKBB_model_ebp_s_log_a_CHOL)
# p = 0.00046

UKBB_model_ebp_s_log_a_HDL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP ~
  rs35724886 +
  Age +
  Sex +
  HDL,
family = binomial
  )
summary(UKBB_model_ebp_s_log_a_HDL)
# p = 0.003

UKBB_model_ebp_s_log_a_lowHDL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP ~
  rs35724886 +
  Age +
  Sex +
  diagnosis_HDL_ext,
family = binomial
  )
summary(UKBB_model_ebp_s_log_a_lowHDL)
# p = 0.002

#------------------------------------------------------------------------------#

UKBB_model_ebp_e_log_a_CHOL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP_ext ~
  rs35724886 +
  Age +
  Sex +
  CHOL,
family = binomial
  )
summary(UKBB_model_ebp_e_log_a_CHOL)
# p = 9.48e-05

UKBB_model_ebp_e_log_a_HDL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP_ext ~
  rs35724886 +
  Age +
  Sex +
  HDL,
family = binomial
  )
summary(UKBB_model_ebp_e_log_a_HDL)
# p = 0.0005

UKBB_model_ebp_e_log_a_lowHDL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP_ext ~
  rs35724886 +
  Age +
  Sex +
  diagnosis_HDL_ext,
family = binomial
  )
summary(UKBB_model_ebp_e_log_a_lowHDL)
# p = 0.0007

#------------------------------------------------------------------------------#

UKBB_model_ebp_s_log_g_CHOL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP ~
  as.factor(rs35724886) +
  Age +
  Sex +
  CHOL,
family = binomial
  )
summary(UKBB_model_ebp_s_log_g_CHOL)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                          Estimate Std. Error  z value Pr(>|z|)    
# (Intercept)            -4.4314640  0.0442663 -100.109  < 2e-16 ***
# as.factor(rs35724886)1 -0.0361699  0.0109972   -3.289  0.00101 ** 
# as.factor(rs35724886)2 -0.0478558  0.0247000   -1.937  0.05269 .  
# ...
# relevant output end ---------------------------------------------------------#

UKBB_model_ebp_s_log_g_HDL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP ~
  as.factor(rs35724886) +
  Age +
  Sex +
  HDL,
family = binomial
  )
summary(UKBB_model_ebp_s_log_g_HDL)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                          Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -3.2389145  0.0431215 -75.111  < 2e-16 ***
# as.factor(rs35724886)1 -0.0308332  0.0114231  -2.699  0.00695 ** 
# as.factor(rs35724886)2 -0.0444509  0.0256302  -1.734  0.08286 .  
# ...
# relevant output end ---------------------------------------------------------#

UKBB_model_ebp_s_log_g_lowHDL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP ~
  as.factor(rs35724886) +
  Age +
  Sex +
  diagnosis_HDL_ext,
family = binomial
  )
summary(UKBB_model_ebp_s_log_g_lowHDL)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                          Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -3.3916928  0.0351756 -96.422  < 2e-16 ***
# as.factor(rs35724886)1 -0.0301720  0.0106811  -2.825  0.00473 ** 
# as.factor(rs35724886)2 -0.0417904  0.0239634  -1.744  0.08117 .  
# ...
# relevant output end ---------------------------------------------------------#


#------------------------------------------------------------------------------#

UKBB_model_ebp_e_log_a_CHOL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP_ext ~
  rs35724886 +
  Age +
  Sex +
  CHOL,
family = binomial
  )
summary(UKBB_model_ebp_e_log_a_CHOL)
readable_logit(UKBB_model_ebp_e_log_a_CHOL)
# p = 9.48e-05

UKBB_model_ebp_e_log_a_HDL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP_ext ~
  rs35724886 +
  Age +
  Sex +
  HDL,
family = binomial
  )
summary(UKBB_model_ebp_e_log_a_HDL)
# p = 0.0005

UKBB_model_ebp_e_log_a_lowHDL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP_ext ~
  rs35724886 +
  Age +
  Sex +
  diagnosis_HDL_ext,
family = binomial
  )
summary(UKBB_model_ebp_e_log_a_lowHDL)
# p = 0.0007

#------------------------------------------------------------------------------#

UKBB_model_ebp_e_log_g_CHOL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP_ext ~
  as.factor(rs35724886) +
  Age +
  Sex +
  CHOL,
family = binomial
  )
summary(UKBB_model_ebp_e_log_g_CHOL)
readable_logit(UKBB_model_ebp_e_log_g_CHOL)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                         Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -4.363744   0.045312 -96.305  < 2e-16 ***
# as.factor(rs35724886)1 -0.042120   0.011522  -3.656 0.000257 ***
# as.factor(rs35724886)2 -0.055972   0.025848  -2.165 0.030353 *  
# ...
# relevant output end ---------------------------------------------------------#

UKBB_model_ebp_e_log_g_HDL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP_ext ~
  as.factor(rs35724886) +
  Age +
  Sex +
  HDL,
family = binomial
  )
summary(UKBB_model_ebp_e_log_g_HDL)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                         Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -3.464293   0.044783 -77.357  < 2e-16 ***
# as.factor(rs35724886)1 -0.037684   0.012042  -3.129  0.00175 ** 
# as.factor(rs35724886)2 -0.056012   0.026998  -2.075  0.03802 *  
# ...
# relevant output end ---------------------------------------------------------#

UKBB_model_ebp_e_log_g_lowHDL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP_ext ~
  as.factor(rs35724886) +
  Age +
  Sex +
  diagnosis_HDL_ext,
family = binomial
  )
summary(UKBB_model_ebp_e_log_g_lowHDL)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                          Estimate Std. Error  z value Pr(>|z|)    
# (Intercept)            -3.9737597  0.0371754 -106.892  < 2e-16 ***
# as.factor(rs35724886)1 -0.0351227  0.0112903   -3.111  0.00187 ** 
# as.factor(rs35724886)2 -0.0492733  0.0252938   -1.948  0.05141 .  
# ...
# relevant output end ---------------------------------------------------------#

#------------------------------------------------------------------------------#

UKBB_model_ebp_e_log_mnv_CHOL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP_ext ~
        MNV_carrier +
        Age +
        Sex +
        CHOL,
      family = binomial
  )
summary(UKBB_model_ebp_e_log_mnv_CHOL)
readable_logit(UKBB_model_ebp_e_log_mnv_CHOL)

UKBB_model_ebp_e_log_wt_CHOL <- 
  glm(data = UKBB_WES,
      formula = diagnosis_BP_ext ~
        WT_carrier +
        Age +
        Sex +
        CHOL,
      family = binomial
  )
summary(UKBB_model_ebp_e_log_wt_CHOL)
readable_logit(UKBB_model_ebp_e_log_wt_CHOL)

# 2.2 HDL cholesterol ----------------------------------------------------------
#  2.2.1 NIHS HDL  -------------------------------------------------------------




NIHS_model_hdl_lin_a <- 
  lm(data = NIHS_2000,
     formula = HDL ~
  ACOT4_simple_geno + 
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  EBP
  )
summary(NIHS_model_hdl_lin_a)
summary(lm.beta(NIHS_model_hdl_lin_a))
confint(lm.beta(NIHS_model_hdl_lin_a))

# p = 0.97


NIHS_model_hdl_lin_g <- 
  lm(data = NIHS_2000,
     formula = HDL ~
  as.factor(ACOT4_simple_geno) + 
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  EBP
  )
summary(NIHS_model_hdl_lin_g)
summary(lm.beta(NIHS_model_hdl_lin_g))
confint(lm.beta(NIHS_model_hdl_lin_g))
# relevant output -------------------------------------------------------------#
# Coefficients:
#                                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                    1.5209129  0.0543481  27.985  < 2e-16 ***
# as.factor(ACOT4_simple_geno)1  0.0200178  0.0302059   0.663   0.5078    
# as.factor(ACOT4_simple_geno)2 -0.0663091  0.0722226  -0.918   0.3589    
# ...
# relevant output end ---------------------------------------------------------#

NIHS_model_hdl_lin_mnv <- 
  lm(data = NIHS_2000,
     formula = HDL ~
  MNV_carrier +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  EBP
  )
summary(NIHS_model_hdl_lin_mnv)
summary(lm.beta(NIHS_model_hdl_lin_mnv))
confint(lm.beta(NIHS_model_hdl_lin_mnv)) %>% round(3)

NIHS_model_hdl_lin_wt <- 
    lm(data = NIHS_2000,
     formula = HDL ~
  WT_carrier +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  EBP
  )
summary(NIHS_model_hdl_lin_wt)
summary(lm.beta(NIHS_model_hdl_lin_wt))
confint(lm.beta(NIHS_model_hdl_lin_wt)) %>% round(3)


#  2.2.2 NIHS low HDL-----------------------------------------------------------

NIHS_2000 <-
  mutate(NIHS_2000,
         lowHDL = case_when(
           Gender == "female" & HDL <  1.3 ~ 1,
Gender ==   "male" & HDL <  1.0 ~ 1,
Gender == "female" & HDL >= 1.3 ~ 0,
Gender ==   "male" & HDL >= 1.0 ~ 0),
.after = LOW_HDL
  )
summary(NIHS_2000$lowHDL)

NIHS_2000_lHDL_TABLE <- 
  table(NIHS_2000$lowHDL,
        NIHS_2000$ACOT4_simple_geno)

NIHS_2000_lHDL_TABLE
#     0   1   2
# 0 275 154  13
# 1  86  45  11

round(
  prop.table(
    NIHS_2000_lHDL_TABLE),
digits = 3)
#       0     1     2
# 0 0.471 0.264 0.022
# 1 0.147 0.077 0.019

fisher.test(NIHS_2000_lHDL_TABLE)
# data:  NIHS_2000_lHDL_TABLE
# p-value = 0.05362
# alternative hypothesis: two.sided

NIHS_model_hdl_log_a <- 
  glm(data = NIHS_2000,
      formula = lowHDL ~
  ACOT4_simple_geno +
  AGE_2000NIHS + 
  Gender +
  COREPED_1 +
  EBP,
family = binomial
  )
readable_logit(NIHS_model_hdl_log_a)
# p = 0.23

NIHS_model_hdl_log_g <- 
  glm(data = NIHS_2000,
      formula = lowHDL ~
  as.factor(ACOT4_simple_geno) +
  AGE_2000NIHS + 
  Gender +
  COREPED_1 +
  EBP,
family = binomial
  )
readable_logit(NIHS_model_hdl_log_g)
# relevant output -------------------------------------------------------------#
#                      Estimate Std. Error   OR_x  2.5 % 97.5 % z value Pr(>|z|)
# (Intercept)           -0.6767     0.3775 0.5083 0.2407 1.0599 -1.7927   0.0730
# as.factor(ACOT4...)1  -0.0735     0.2173 0.9291 0.6034 1.4165 -0.3384   0.7351
# as.factor(ACOT4...)2   1.0554     0.4427 2.8732 1.1884 6.8618  2.3840   0.0171
# ...
# relevant output end ---------------------------------------------------------#

NIHS_model_hdl_log_mnv <- 
  glm(data = NIHS_2000,
      formula = lowHDL ~
  MNV_carrier +
  AGE_2000NIHS + 
  Gender +
  COREPED_1 +
  EBP,
family = binomial
  )
readable_logit(NIHS_model_hdl_log_mnv)
  
NIHS_model_hdl_log_wt <- 
  glm(data = NIHS_2000,
      formula = lowHDL ~
  WT_carrier +
  AGE_2000NIHS + 
  Gender +
  COREPED_1 +
  EBP,
family = binomial
  )
readable_logit(NIHS_model_hdl_log_wt)
  


#  2.2.3 UKBB HDL --------------------------------------------------------------
UKBB_model_hdl_lin_a_dBPe <- 
  lm(data = UKBB_WES,
     formula = HDL ~ 
  rs35724886 +
  Age +
  Sex +
  chol_meds +
  diagnosis_BP_ext
  )
summary(UKBB_model_hdl_lin_a_dBPe)
summary(lm.beta(UKBB_model_hdl_lin_a_dBPe))
confint(lm.beta(UKBB_model_hdl_lin_a_dBPe)) %>% round(4)
# p = 0.15

UKBB_model_hdl_lin_g_dBPe <- 
  lm(data = UKBB_WES,
     formula = HDL ~ 
  as.factor(rs35724886) +
  Age +
  Sex +
    chol_meds +
    diagnosis_BP_ext
  )
summary(UKBB_model_hdl_lin_g_dBPe)
summary(lm.beta(UKBB_model_hdl_lin_g_dBPe)) %>% coef() %>% round(4)
confint(lm.beta(UKBB_model_hdl_lin_g_dBPe)) %>% round(3)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                          Estimate Std. Error  t value Pr(>|t|)    
# (Intercept)             1.4294145  0.0059322  240.960  < 2e-16 ***
# as.factor(rs35724886)1  0.0055996  0.0017965    3.117  0.00183 ** 
# as.factor(rs35724886)2 -0.0050707  0.0040423   -1.254  0.20969    
# ...
# relevant output end ---------------------------------------------------------#

UKBB_model_hdl_lin_mnv_dBPe <- 
  lm(data = UKBB_WES,
     formula = HDL ~ 
  MNV_carrier +
  Age +
  Sex +
    chol_meds +
      diagnosis_BP_ext
  )
summary(UKBB_model_hdl_lin_mnv_dBPe)
summary(lm.beta(UKBB_model_hdl_lin_mnv_dBPe)) %>% coef() %>% round(4)
confint(lm.beta(UKBB_model_hdl_lin_mnv_dBPe)) %>% round(3)


UKBB_model_hdl_lin_wt_dBPe <- 
  lm(data = UKBB_WES,
     formula = HDL ~ 
  WT_carrier +
  Age +
  Sex +
    chol_meds +
    
  diagnosis_BP_ext
  )
summary(UKBB_model_hdl_lin_wt_dBPe)
summary(lm.beta(UKBB_model_hdl_lin_wt_dBPe)) %>% coef() %>% round(4)
confint(lm.beta(UKBB_model_hdl_lin_wt_dBPe)) %>% round(3)



#  2.2.4 UKBB low HDL ----------------------------------------------------------
UKBB_model_hdl_log_a_dBPe <- 
  glm(data = UKBB_WES,
      formula = diagnosis_HDL_ext ~
  rs35724886 +
  Age +
  Sex +
  diagnosis_BP_ext,
family = binomial
  )
summary(UKBB_model_hdl_log_a_dBPe)
readable_logit(UKBB_model_hdl_log_a_dBPe)
# p = 0.0003

UKBB_model_hdl_log_g_dBPe <- 
  glm(data = UKBB_WES,
      formula = diagnosis_HDL_ext ~
  as.factor(rs35724886) +
  Age +
  Sex +
  diagnosis_BP_ext,
family = binomial
  )
summary(UKBB_model_hdl_log_g_dBPe)
readable_logit(UKBB_model_hdl_log_g_dBPe)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                          Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -3.1743621  0.0373180 -85.063  < 2e-16 ***
# as.factor(rs35724886)1 -0.0373985  0.0107348  -3.484 0.000494 ***
# as.factor(rs35724886)2 -0.0472731  0.0242653  -1.948 0.051394 .  
# ...
# relevant output end ---------------------------------------------------------#

UKBB_model_hdl_log_mnv_dBPe <- 
  glm(data = UKBB_WES,
      formula = diagnosis_HDL_ext ~
  MNV_carrier +
  Age +
  Sex +
  diagnosis_BP_ext,
family = binomial
  )
summary(UKBB_model_hdl_log_mnv_dBPe)
readable_logit(UKBB_model_hdl_log_mnv_dBPe)

UKBB_model_hdl_log_wt_dBPe <- 
  glm(data = UKBB_WES,
      formula = diagnosis_HDL_ext ~
  WT_carrier +
  Age +
  Sex +
  diagnosis_BP_ext,
family = binomial
  )
summary(UKBB_model_hdl_log_wt_dBPe)
readable_logit(UKBB_model_hdl_log_wt_dBPe)


# 2.3 LDL cholesterol
#  2.3.1 NIHS LDL --------------------------------------------------------------
NIHS_model_ldl_lin_a <- 
  lm(data = NIHS_2000,
     formula = LDL ~ 
       ACOT4_simple_geno +
       AGE_2000NIHS +
       Gender +
       COREPED_1)
summary(NIHS_model_ldl_lin_a)

NIHS_model_ldl_lin_g <- 
  lm(data = NIHS_2000,
     formula = LDL ~ 
       as.factor(ACOT4_simple_geno) +
       AGE_2000NIHS +
       Gender +
       COREPED_1)
summary(NIHS_model_ldl_lin_g)

NIHS_model_ldl_lin_mnv <- 
  lm(data = NIHS_2000,
     formula = LDL ~ 
       MNV_carrier +
       AGE_2000NIHS +
       Gender +
       COREPED_1)
summary(NIHS_model_ldl_lin_mnv)

NIHS_model_ldl_lin_wt <- 
  lm(data = NIHS_2000,
     formula = LDL ~ 
       WT_carrier +
       AGE_2000NIHS +
       Gender +
       COREPED_1)
summary(NIHS_model_ldl_lin_wt)

#  2.3.2 UKBB LDL ---------------------------------------------------------------
UKBB_model_ldl_lin_a_dBPe <- 
  lm(data = UKBB_WES,
     formula = LDL ~
       rs35724886 +
       Age +
       Sex +
       chol_meds +
       diagnosis_BP_ext)
summary(UKBB_model_ldl_lin_a_dBPe)

UKBB_model_ldl_lin_g_dBPe <- 
  lm(data = UKBB_WES,
     formula = LDL ~
       as.factor(rs35724886) +
       Age +
       Sex +
chol_meds+
              diagnosis_BP_ext)
summary(UKBB_model_ldl_lin_g_dBPe)

UKBB_model_ldl_lin_mnv_dBPe <- 
  lm(data = UKBB_WES,
     formula = LDL ~
       MNV_carrier +
       Age +
       Sex +
chol_meds +
              diagnosis_BP_ext)
summary(UKBB_model_ldl_lin_mnv_dBPe)

UKBB_model_ldl_lin_wt_dBPe <- 
  lm(data = UKBB_WES,
     formula = LDL ~
       WT_carrier +
       Age +
       Sex +
       chol_meds +
       diagnosis_BP_ext)
summary(UKBB_model_ldl_lin_wt_dBPe)


# 2.3 total cholesterol --------------------------------------------------------
#  2.3.1 NIHS total cholesterol ------------------------------------------------
NIHS_model_chol_lin_a <-
  lm(data = NIHS_2000,
     formula = CHOL ~
  ACOT4_simple_geno +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  EBP
  )
summary(NIHS_model_chol_lin_a)
summary(lm.beta(NIHS_model_chol_lin_a))
confint(lm.beta(NIHS_model_chol_lin_a))
#p = 0.044

NIHS_model_chol_lin_g <-
  lm(data = NIHS_2000,
     formula = CHOL ~
  as.factor(ACOT4_simple_geno) +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  EBP
  )
summary(NIHS_model_chol_lin_g)
summary(lm.beta(NIHS_model_chol_lin_g))
confint(lm.beta(NIHS_model_chol_lin_g))
# relevant output -------------------------------------------------------------#
# Coefficients:
#                                Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                    4.494824   0.169324  26.546  < 2e-16 ***
# as.factor(ACOT4_simple_geno)1  0.211774   0.093948   2.254  0.02456 *  
# as.factor(ACOT4_simple_geno)2  0.142896   0.221353   0.646  0.51883    
# ...
# relevant output end ---------------------------------------------------------#

NIHS_model_chol_lin_mnv <- 
  lm(data = NIHS_2000,
     formula = CHOL ~
  MNV_carrier +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  EBP
  )
summary(NIHS_model_chol_lin_mnv)
summary(lm.beta(NIHS_model_chol_lin_mnv))
confint(lm.beta(NIHS_model_chol_lin_mnv)) %>% round(3)
  

NIHS_model_chol_lin_wt <- 
  lm(data = NIHS_2000,
     formula = CHOL ~
  WT_carrier +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  EBP
  )
summary(NIHS_model_chol_lin_wt)
summary(lm.beta(NIHS_model_chol_lin_wt))
confint(lm.beta(NIHS_model_chol_lin_wt)) %>% round(3)
  

#  2.3.2 NIHS high cholesterol -------------------------------------------------
NIHS_2000 <- 
  mutate(NIHS_2000,
         high_chol = case_when(
           CHOL >= 5.5 ~ 1,
CHOL <  5.5 ~ 0),
.after = HIGH_CHOL
  )
summary(NIHS_2000$high_chol)

NIHS_2000_hiCHOL_TABLE <- 
  table(NIHS_2000$high_chol,
        NIHS_2000$ACOT4_simple_geno)

NIHS_2000_hiCHOL_TABLE
#     0   1   2
# 0 182  93  12
# 1 180 108  13

round(
  prop.table(
    NIHS_2000_hiCHOL_TABLE),
digits = 3)


fisher.test(NIHS_2000_hiCHOL_TABLE)
# data:  NIHS_2000_hiCHOL_TABLE
# p-value = 0.6577
# alternative hypothesis: two.sided

NIHS_model_chol_log_a <- 
  glm(data = NIHS_2000,
      formula = high_chol ~
  ACOT4_simple_geno +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  EBP,
family = binomial
  )
summary(NIHS_model_chol_log_a)
readable_logit(NIHS_model_chol_log_a)
# p = 0.29

NIHS_model_chol_log_g <- 
  glm(data = NIHS_2000,
      formula = high_chol ~
  as.factor(ACOT4_simple_geno) +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  EBP,
family = binomial
  )
summary(NIHS_model_chol_log_g)
readable_logit(NIHS_model_chol_log_g)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                      Estimate Std. Error   OR_x  2.5 % 97.5 % z value Pr(>|z|)
# (Intercept)           -1.5723     0.3407 0.2076 0.1053 0.4013 -4.6147   0.0000
# as.factor(ACOT4...)1   0.1544     0.1837 1.1669 0.8143 1.6744  0.8402   0.4008
# as.factor(ACOT4...)2   0.3265     0.4291 1.3861 0.5954 3.2556  0.7608   0.4468
# ... 
# relevant output end ---------------------------------------------------------#

NIHS_model_chol_log_mnv <- 
  glm(data = NIHS_2000,
      formula = high_chol ~
  MNV_carrier +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  EBP,
family = binomial
  )
summary(NIHS_model_chol_log_mnv)
readable_logit(NIHS_model_chol_log_mnv)


NIHS_model_chol_log_wt <- 
  glm(data = NIHS_2000,
      formula = high_chol ~
  WT_carrier +
  AGE_2000NIHS +
  Gender +
  COREPED_1 +
  EBP,
family = binomial
  )
summary(NIHS_model_chol_log_wt)
readable_logit(NIHS_model_chol_log_wt)




#  2.3.3 UKBB total cholesterol ------------------------------------------------
UKBB_model_chol_lin_a_dBPe <- 
  lm(data = UKBB_WES,
     formula = CHOL ~ 
  rs35724886 +
  Age +
  Sex +
  chol_meds +  # high chol --> chol meds --> normal chol ? -JMB
  diagnosis_BP_ext
  )
summary(UKBB_model_chol_lin_a_dBPe)
summary(lm.beta(UKBB_model_chol_lin_a_dBPe)) 
confint(lm.beta(UKBB_model_chol_lin_a_dBPe))%>% round(3)
# p = 0.0006

UKBB_model_chol_lin_g_dBPe <- 
  lm(data = UKBB_WES,
     formula = CHOL ~ 
  as.factor(rs35724886) +
  Age +
  Sex +
  chol_meds +  # high chol --> chol meds --> normal chol ? -JMB
  diagnosis_BP_ext
  )
summary(UKBB_model_chol_lin_g_dBPe)
summary(lm.beta(UKBB_model_chol_lin_g_dBPe))
confint(lm.beta(UKBB_model_chol_lin_g_dBPe)) %>% round(3)

# relevant output -------------------------------------------------------------#
# Coefficients:
#                          Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             5.3806795  0.0182714 294.486  < 2e-16 ***
# as.factor(rs35724886)1  0.0181780  0.0055356   3.284  0.00102 ** 
# as.factor(rs35724886)2  0.0231207  0.0124736   1.854  0.06380 .  
# ...
# relevant output end ---------------------------------------------------------#

UKBB_model_chol_lin_mnv_dBPe <- 
  lm(data = UKBB_WES,
     formula = CHOL ~ 
  MNV_carrier +
  Age +
  Sex +
  chol_meds +  # high chol --> chol meds --> normal chol ? -JMB
  diagnosis_BP_ext
  )
summary(UKBB_model_chol_lin_mnv_dBPe)
summary(lm.beta(UKBB_model_chol_lin_mnv_dBPe))
confint(lm.beta(UKBB_model_chol_lin_mnv_dBPe)) %>% round(3)

UKBB_model_chol_lin_wt_dBPe <- 
  lm(data = UKBB_WES,
     formula = CHOL ~ 
  WT_carrier +
  Age +
  Sex +
  chol_meds +  # high chol --> chol meds --> normal chol ? -JMB
  diagnosis_BP_ext
  )
summary(UKBB_model_chol_lin_wt_dBPe)
summary(lm.beta(UKBB_model_chol_lin_wt_dBPe))
confint(lm.beta(UKBB_model_chol_lin_wt_dBPe)) %>% round(3)

#  2.3.4 UKBB high cholesterol -------------------------------------------------
UKBB_WES <- 
  mutate(UKBB_WES,
         HiChol_or_meds = case_when(
           CHOL >= 5.5 | chol_meds == 1 ~ 1,
CHOL <  5.5 & chol_meds == 0 ~ 0),
.after = chol_ge_5.5mmolL
  )
summary(UKBB_WES$HiChol_or_meds)


UKBB_model_chol_log_a_dBPe <- 
  glm(data = UKBB_WES,
      formula = HiChol_or_meds ~
  rs35724886 +
  Age +
  Sex +
  diagnosis_BP_ext,
family = binomial
  )
summary(UKBB_model_chol_log_a_dBPe)
readable_logit(UKBB_model_chol_log_a_dBPe)
# p = 0.67

UKBB_model_chol_log_g_dBPe <- 
  glm(data = UKBB_WES,
      formula = HiChol_or_meds ~
  as.factor(rs35724886) +
  Age +
  Sex +
  diagnosis_BP_ext,
family = binomial
  )
summary(UKBB_model_chol_log_g_dBPe)
readable_logit(UKBB_model_chol_log_g_dBPe)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                          Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -0.4179418  0.0330953 -12.628   <2e-16 ***
# as.factor(rs35724886)1  0.0199015  0.0100696   1.976   0.0481 *  
# as.factor(rs35724886)2  0.0219644  0.0227060   0.967   0.3334    
# ...
# relevant output end ---------------------------------------------------------#

UKBB_model_chol_log_mnv_dBPe <- 
  glm(data = UKBB_WES,
      formula = HiChol_or_meds ~
  MNV_carrier +
  Age +
  Sex +
  diagnosis_BP_ext,
family = binomial
  )
summary(UKBB_model_chol_log_mnv_dBPe)
readable_logit(UKBB_model_chol_log_mnv_dBPe)

UKBB_model_chol_log_wt_dBPe <- 
  glm(data = UKBB_WES,
      formula = HiChol_or_meds ~
  WT_carrier +
  Age +
  Sex +
  diagnosis_BP_ext,
family = binomial
  )
summary(UKBB_model_chol_log_wt_dBPe)
readable_logit(UKBB_model_chol_log_wt_dBPe)



# 2.4 stroke -------------------------------------------------------------------
#  2.4.1 UKBB stroke -----------------------------------------------------------
#  no stroke data for NIHS 2000

summary(UKBB_WES$stroke)
summary(as.factor(UKBB_WES$stroke))

UKBB_model_stroke_log_a_ <- 
  glm(data = UKBB_WES,
      formula = stroke ~
  rs35724886 +
  Age +
  Sex,
family = binomial
  )
summary(UKBB_model_stroke_log_a_)
readable_logit(UKBB_model_stroke_log_a_)

# p = 0.72

UKBB_model_stroke_log_g_ <- 
  glm(data = UKBB_WES,
      formula = stroke ~
  as.factor(rs35724886) +
  Age +
  Sex,
family = binomial
  )
summary(UKBB_model_stroke_log_g_)
readable_logit(UKBB_model_stroke_log_g_)

UKBB_model_stroke_log_mnv_ <- 
  glm(data = UKBB_WES,
      formula = stroke ~
  MNV_carrier +
  Age +
  Sex,
family = binomial
  )
summary(UKBB_model_stroke_log_mnv_)
readable_logit(UKBB_model_stroke_log_mnv_)

UKBB_model_stroke_log_wt_ <- 
  glm(data = UKBB_WES,
      formula = stroke ~
  WT_carrier +
  Age +
  Sex,
family = binomial
  )
summary(UKBB_model_stroke_log_wt_)
readable_logit(UKBB_model_stroke_log_wt_)





# relevant output -------------------------------------------------------------#
# Coefficients:
#                         Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -9.171477   0.175116 -52.374   <2e-16 ***
# as.factor(rs35724886)1  0.068123   0.040578   1.679   0.0932 .  
# as.factor(rs35724886)2 -0.140715   0.099617  -1.413   0.1578    
# ...
# relevant output end ---------------------------------------------------------#

#------------------------------------------------------------------------------#
UKBB_model_stroke_log_a_hdl_ebp <- 
  glm(data = UKBB_WES,
      formula = stroke ~
  rs35724886 +
  Age +
  Sex +
  diagnosis_HDL_ext +
  diagnosis_BP_ext,
family = binomial
  )
summary(UKBB_model_stroke_log_a_hdl_ebp)
# p = 0.41

UKBB_model_stroke_log_g_hdl_ebp <- 
  glm(data = UKBB_WES,
      formula = stroke ~
  as.factor(rs35724886) +
  Age +
  Sex +
  diagnosis_HDL_ext +
  diagnosis_BP_ext,
family = binomial
  )
summary(UKBB_model_stroke_log_g_hdl_ebp)
# relevant output -------------------------------------------------------------#
# Coefficients:
#                         Estimate Std. Error z value Pr(>|z|)    
# (Intercept)            -8.938775   0.173520 -51.514  < 2e-16 ***
# as.factor(rs35724886)1  0.087408   0.040926   2.136   0.0327 *  
# as.factor(rs35724886)2 -0.118887   0.100336  -1.185   0.2361    
# ...
# relevant output end ---------------------------------------------------------#

UKBB_model_stroke_log_mnv_hdl_ebp <- 
  glm(data = UKBB_WES,
      formula = stroke ~
  MNV_carrier +
  Age +
  Sex +
  diagnosis_HDL_ext +
  diagnosis_BP_ext,
family = binomial
  )
summary(UKBB_model_stroke_log_mnv_hdl_ebp)
readable_logit(UKBB_model_stroke_log_mnv_hdl_ebp)


UKBB_model_stroke_log_wt_hdl_ebp <- 
  glm(data = UKBB_WES,
      formula = stroke ~
  WT_carrier +
  Age +
  Sex +
  diagnosis_HDL_ext +
  diagnosis_BP_ext,
family = binomial
  )
summary(UKBB_model_stroke_log_wt_hdl_ebp)
readable_logit(UKBB_model_stroke_log_wt_hdl_ebp)




# 2.5 medication ---------------------------------------------------------------

bp_meds_table <- table(UKBB_WES$BP_meds,
                       UKBB_WES$rs35724886)
bp_meds_table
#        0      1      2
# 0 100433  52094   7338
# 1  25833  13144   1768
prop.table(bp_meds_table, margin = 1)
#            0          1          2
# 0 0.62823632 0.32586245 0.04590123
# 1 0.63401644 0.32259173 0.04339183
chisq.test(bp_meds_table)
# X-squared = 7.2961, df = 2, p-value = 0.02604

plot(bp_meds_table)

UKBB_model_BP_meds_log_a <- 
  glm(data = UKBB_WES,
      formula = BP_meds ~
        rs35724886 +
        Age +
        Sex,
      family = binomial
)

summary(UKBB_model_BP_meds_log_a)
readable_logit(UKBB_model_BP_meds_log_a)
#             Estimate Std. Error   OR_x  2.5 % 97.5 %   z value Pr(>|z|)
# (Intercept)  -7.1436     0.0513 0.0008 0.0007 0.0009 -139.2558        0
# rs35724886   -0.0482     0.0101 0.9530 0.9343 0.9720   -4.7739 1.81e-06
# Age           0.0962     0.0008 1.1010 1.0991 1.1028  113.9540        0
# Sexmale       0.4071     0.0116 1.5024 1.4687 1.5369   35.1133        0


UKBB_model_BP_meds_log_g <- 
  glm(data = UKBB_WES,
      formula = BP_meds ~
        as.factor(rs35724886) +
        Age +
        Sex,
      family = binomial
)

summary(UKBB_model_BP_meds_log_g)
readable_logit(UKBB_model_BP_meds_log_g)
#                        Estimate Std. Error   OR_x  2.5 % 97.5 %   z value Pr(>|z|)
# (Intercept)             -7.1433     0.0513 0.0008 0.0007 0.0009 -139.2170   0.0000
# as.factor(rs35724886)1  -0.0501     0.0125 0.9511 0.9280 0.9747   -4.0018   0.0001
# as.factor(rs35724886)2  -0.0910     0.0286 0.9130 0.8631 0.9654   -3.1856   0.0014
# Age                      0.0962     0.0008 1.1010 1.0991 1.1028  113.9527   0.0000
# Sexmale                  0.4071     0.0116 1.5024 1.4687 1.5369   35.1139   0.0000


UKBB_model_BP_meds_log_wt <- 
  glm(data = UKBB_WES,
      formula = BP_meds ~
        WT_carrier +
        Age +
        Sex,
      family = binomial
)

summary(UKBB_model_BP_meds_log_wt)
readable_logit(UKBB_model_BP_meds_log_wt)
#             Estimate Std. Error   OR_x  2.5 % 97.5 %   z value Pr(>|z|)
# (Intercept)  -7.2307     0.0580 0.0007 0.0006 0.0008 -124.6539   0.0000
# WT_carrier    0.0738     0.0282 1.0766 1.0189 1.1382    2.6147   0.0089
# Age           0.0961     0.0008 1.1009 1.0991 1.1027  113.9054   0.0000
# Sexmale       0.4071     0.0116 1.5025 1.4687 1.5370   35.1177   0.0000


UKBB_model_BP_meds_log_mnv <- 
  glm(data = UKBB_WES,
      formula = BP_meds ~
        MNV_carrier +
        Age +
        Sex,
      family = binomial
)

summary(UKBB_model_BP_meds_log_mnv)
readable_logit(UKBB_model_BP_meds_log_mnv)

#------------------------------------------------------------------------------#
chol_meds_table <- table(UKBB_WES$chol_meds,
                         UKBB_WES$rs35724886)
chol_meds_table
#        0      1      2
# 0 104955  54243   7628
# 1  21311  10995   1478

UKBB_model_CHOL_meds_log_a <- 
  glm(data = UKBB_WES,
      formula = chol_meds ~
        rs35724886 +
        Age +
        Sex,
      family = binomial
)

summary(UKBB_model_CHOL_meds_log_a)
readable_logit(UKBB_model_CHOL_meds_log_a)
# Estimate Std. Error   OR_x  2.5 % 97.5 %   z value Pr(>|z|)
# (Intercept)  -8.3752     0.0591 0.0002 0.0002 0.0003 -141.6350    0.000
# rs35724886   -0.0338     0.0109 0.9667 0.9462 0.9876   -3.0951    0.002
# Age           0.1095     0.0010 1.1157 1.1136 1.1178  114.1649    0.000
# Sexmale       0.7290     0.0127 2.0731 2.0223 2.1252   57.5376    0.000

UKBB_model_CHOL_meds_log_g <- 
  glm(data = UKBB_WES,
      formula = chol_meds ~
        as.factor(rs35724886) +
        Age +
        Sex,
      family = binomial
)

summary(UKBB_model_CHOL_meds_log_g)
readable_logit(UKBB_model_CHOL_meds_log_g)
#                        Estimate Std. Error   OR_x  2.5 % 97.5 %   z value Pr(>|z|)
# (Intercept)             -8.3752     0.0591 0.0002 0.0002 0.0003 -141.6052   0.0000
# as.factor(rs35724886)1  -0.0339     0.0136 0.9667 0.9413 0.9927   -2.4976   0.0125
# as.factor(rs35724886)2  -0.0674     0.0309 0.9348 0.8796 0.9929   -2.1816   0.0291
# Age                      0.1095     0.0010 1.1157 1.1136 1.1178  114.1619   0.0000
# Sexmale                  0.7290     0.0127 2.0731 2.0223 2.1252   57.5374   0.0000


UKBB_model_CHOL_meds_log_wt <- 
  glm(data = UKBB_WES,
      formula = chol_meds ~
        WT_carrier +
        Age +
        Sex,
      family = binomial
)

summary(UKBB_model_CHOL_meds_log_wt)
readable_logit(UKBB_model_CHOL_meds_log_wt)

UKBB_model_CHOL_meds_log_mnv <- 
  glm(data = UKBB_WES,
      formula = chol_meds ~
        MNV_carrier +
        Age +
        Sex,
      family = binomial
)

summary(UKBB_model_CHOL_meds_log_mnv)
readable_logit(UKBB_model_CHOL_meds_log_mnv)








# 3.0 plots ####################################################################

# 3.1 blood pressure -----------------------------------------------------------

ggplot(data = NIHS_2000, 
       mapping = aes(x = ACOT4_simple_geno, y = SBP, group = ACOT4_simple_geno))+
  geom_jitter(width = 0.33, na.rm = TRUE, color = "darkgrey", alpha = 0.7)+
  geom_boxplot(outliers = FALSE, na.rm = TRUE, alpha = 0.4, fill = "lightblue") +
  scale_x_discrete()+
  xlab("ACOT4 MNV genotype")+
  ylab("Systolic blood pressure (mmHg)")+
  theme_dark()+
  theme(plot.background = element_rect(fill = "lightgrey"),
        panel.background = element_rect(fill = "grey"),
        panel.grid = element_blank())


ggplot(data = UKBB_WES, 
       mapping = aes(x = rs35724886, y = av_sbp, group = rs35724886))+
  geom_jitter(width = 0.33, na.rm = TRUE, color = "darkgrey", alpha = 0.7)+
  geom_boxplot(outliers = FALSE, na.rm = TRUE, alpha = 0.4, fill = "lightgreen") +
  scale_x_discrete()+
  xlab("ACOT4 MNV genotype")+
  ylab("Systolic blood pressure (mmHg)")+
  theme_dark()+
  theme(plot.background = element_rect(fill = "lightgrey"),
        panel.background = element_rect(fill = "grey"),
        panel.grid = element_blank())

column_plot_data <- UKBB_WES %>% 
 arrange(rs35724886,BP_meds)

# column_plot_data <- 
  # column_plot_data %>% mutate(bp_meds_prop =  

sum(UKBB_WES$BP_meds, na.rm = TRUE)/sum(!is.na(UKBB_WES$BP_meds))

bp_meds_counts <- as.data.frame(bp_meds_table)

ggplot(bp_meds_counts,
       aes(x = Var2, y = Freq, fill = Var1))+
  geom_col(position = "dodge")




ggplot(data = UKBB_WES,
       mapping 



# 3.2 hdl ----------------------------------------------------------------------

ggplot(data = NIHS_2000, 
       mapping = aes(x = ACOT4_simple_geno, y = HDL, group = ACOT4_simple_geno))+
  geom_jitter(width = 0.33, na.rm = TRUE, color = "darkgrey", alpha = 0.7)+
  geom_boxplot(outliers = FALSE, na.rm = TRUE, alpha = 0.4, fill = "lightblue") +
  scale_x_discrete()+
  xlab("ACOT4 MNV genotype")+
  ylab("HDL cholesterol (mmol/L)")+
  theme_dark()+
  theme(plot.background = element_rect(fill = "lightgrey"),
        panel.background = element_rect(fill = "grey"),
        panel.grid = element_blank())

ggplot(data = UKBB_WES, 
       mapping = aes(x = rs35724886, y = HDL, group = rs35724886))+
  geom_jitter(width = 0.33, na.rm = TRUE, color = "darkgrey", alpha = 0.7)+
  geom_boxplot(outliers = FALSE, na.rm = TRUE, alpha = 0.4, fill = "lightgreen") +
  scale_x_discrete()+
  xlab("ACOT4 MNV genotype")+
  ylab("HDL cholesterol (mmol/L)")+
  theme_dark()+
  theme(plot.background = element_rect(fill = "lightgrey"),
        panel.background = element_rect(fill = "grey"),
        panel.grid = element_blank())

# 3.3 total cholesterol --------------------------------------------------------

ggplot(data = NIHS_2000, 
       mapping = aes(x = ACOT4_simple_geno, y = CHOL, group = ACOT4_simple_geno))+
  geom_jitter(width = 0.33, na.rm = TRUE, color = "darkgrey", alpha = 0.7)+
  geom_boxplot(outliers = FALSE, na.rm = TRUE, alpha = 0.4, fill = "lightblue") +
  scale_x_discrete()+
  xlab("ACOT4 MNV genotype")+
  ylab("Total cholesterol (mmol/L)")+
  theme_dark()+
  theme(plot.background = element_rect(fill = "lightgrey"),
        panel.background = element_rect(fill = "grey"),
        panel.grid = element_blank())

ggplot(data = UKBB_WES, 
       mapping = aes(x = rs35724886, y = CHOL, group = rs35724886))+
  geom_jitter(width = 0.33, na.rm = TRUE, color = "darkgrey", alpha = 0.7)+
  geom_boxplot(outliers = FALSE, na.rm = TRUE, alpha = 0.4, fill = "lightgreen") +
  scale_x_discrete()+
  xlab("ACOT4 MNV genotype")+
  ylab("Total cholesterol (mmol/L)")+
  theme_dark()+
  theme(plot.background = element_rect(fill = "lightgrey"),
        panel.background = element_rect(fill = "grey"),
        panel.grid = element_blank() )

# 3.4 medication use -----------------------------------------------------------

Medication_data <- UKBB_WES %>% select(c(eid, rs35724886, BP_meds, chol_meds))

Medication_data %<>% pivot_longer(cols = c(BP_meds, chol_meds),
                                 names_to = "meds_type",
                                 values_to = "meds_use")

ggplot(filter(Medication_data, !is.na(rs35724886)),
       aes(x = meds_type))+
  geom_bar(aes(fill = as.factor(meds_use)))+
  facet_grid( ~ rs35724886)

ggplot(filter(Medication_data, !is.na(rs35724886)),
       aes(x = meds_type))+
  geom_bar(aes(fill = as.factor(meds_use)),
           # position = "fill"
           )+
  facet_grid(~ rs35724886, 
             # switch = "",
             labeller = as_labeller(c("0" = "WT", "1" = "Het", "2" = "Hom")))+
  theme_blank()+
  theme(axis.text.x = element_text(angle = -45, hjust = 0))+
  labs(#title = "Medication use according to ACOT4\nMNV genotype in the UK Biobank",
  x = "medication type",
  y = "count",
  fill = "medication use")+
  scale_x_discrete(labels = c("BP", "cholesterol"))+
  scale_y_continuous(breaks = c(0, 50e3, 100e3), labels = c("0", "50k", "100k"))+
  scale_fill_brewer(labels = c("0" = "no", "1" = "yes"),
                    palette = "Set2"  )
  



# Medication_data %>% 
#   filter(rs35724886 == 0) %>% # nrow()  # 126266
#   select(BP_meds) %>% #sum(na.rm = TRUE)#  25833
#   is.na() %>% sum() #                          0
# 
# Medication_data %>% 
#   filter(rs35724886 == 1) %>%  #nrow()   # 65238
#   select(BP_meds) %>% #sum(na.rm = TRUE)#  13144
#   is.na() %>% sum() #                          0
# 
# Medication_data %>% 
#   filter(rs35724886 == 2) %>%  #nrow()   # 9106
#   select(BP_meds) %>% #sum(na.rm = TRUE)#  1768
#   is.na() %>% sum() #                         0
# 
# Medication_data %>% 
#   filter(rs35724886 == 0) %>% # nrow()    # 126266
#   select(chol_meds) %>% #sum(na.rm = TRUE)#  21311
#   is.na() %>% sum() #                            0
# 
# Medication_data %>% 
#   filter(rs35724886 == 1) %>%  #nrow()    # 65238
#   select(chol_meds) %>% #sum(na.rm = TRUE)# 10995
#   is.na() %>% sum() #                          0
# 
# Medication_data %>% 
#   filter(rs35724886 == 2) %>%  #nrow()    #   9106
#   select(chol_meds) %>% #sum(na.rm = TRUE)#   1478
#   is.na() %>% sum() #                            0
# 
# # geno_count = c(126266, 65238, 9106),
# 
# rs35724886_fct = factor(
#   c(rep("wt", times = 126266),
#     rep("ht", times =  65238),
#     rep("hm", times =   9106))
# )
# BP_MEDS_count = factor(
#   c(rep("yes", times = 25833), rep("no", times = (126266 - 25833)),
#     rep("yes", times = 13144), rep("no", times =  (65238 - 13144)),
#     rep("yes", times =  1768), rep("no", times =   (9106 -  1768)))
# )
# 
# 
# CHOL_MEDS_count = factor(
#   c(rep("yes", times = 21311), rep("no", times = (126266 - 21311)),
#     rep("yes", times = 10995), rep("no", times =  (65238 - 10995)),
#     rep("yes", times =  1478), rep("no", times =   (9106 -  1478)))
# )
#                     )
# PLOT_Medication_data <- data.frame(
#   rs35724886_fct = rs35724886_fct,
#   BP_MEDS_count = BP_MEDS_count,
#   CHOL_MEDS_count= CHOL_MEDS_count
# )
# 
# PLOT_Medication_data %>% slice_sample(n = 10) 
# 
# PLOT_Medication_data %<>% pivot_longer(cols = -rs35724886_fct, 
#                                        names_to = "meds_type",
#                                        values_to = "use_meds")




# Footer - tidy up and quit ####################################################
# dev.off()
# rm(list = ls())
# #q()