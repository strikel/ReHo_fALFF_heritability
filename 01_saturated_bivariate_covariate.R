rm(list = ls())
setwd("~/GitHub/qtab/rsfMRI/CerebCortex_ms")

#### Load libraries ####
library(OpenMx)
library(dplyr)
library(stringr)
source("miFunctions.R")
mxOption(NULL, "Default optimizer", "NPSOL") # NPSOL SLSQP CSOLNP
mxOption(NULL, key="Number of Threads", value=4)

Saturated_bivariate_covariate <- function(phenotype, twin.data) {
  print(phenotype)
  covariates <- c("sex_female", "age_years_ses01", "age_years_ses02", "fd_meanZ_ses01", "fd_meanZ_ses02")
  nc <- length(covariates)
  
  # Recode missing covariates to use definition variables
  for (x in covariates) {
    twin.data[, paste0(x, "_01")][is.na(twin.data[, paste0(x, "_01")])] <- -999
    twin.data[, paste0(x, "_02")][is.na(twin.data[, paste0(x, "_02")])] <- -999
  }
  
  # Select variables
  phenotypes <- c(
    paste0(phenotype, "_ses01"),
    paste0(phenotype, "_ses02")
  )
  selVars <- c(paste0(phenotypes, "_01"), paste0(phenotypes, "_02"))
  covVars <- c(paste0(covariates, "_01"), paste0(covariates, "_02"))
  
  # Select Data for Analysis
  useVars <- c(selVars, covVars)
  mzData <- subset(twin.data, zyg == 1 | zyg == 2, useVars)
  dzData <- subset(twin.data, zyg == 3 | zyg == 4 | zyg == 5, useVars)
  
  # Set Starting Values
  nt <- 2 # number of twins/sibs
  nv <- length(phenotypes) # number of variables
  ntv <- nv * nt # number of total variables
  svBe <- 0.01 # start value for regressions
  svMe <- colMeans(twin.data[, selVars], na.rm = T)
  svVa <- diag(cov(twin.data[, selVars], use = "pairwise.complete.obs"))
  
  # # Raw means (SD)
  mean.val.ses01 <- mean(unlist(twin.data[, c(selVars[1], selVars[3])]), na.rm = T)
  mean.val.ses02 <- mean(unlist(twin.data[, c(selVars[2], selVars[4])]), na.rm = T)
  sd.val.ses01 <- sd(unlist(twin.data[, c(selVars[1], selVars[3])]), na.rm = T)
  sd.val.ses02 <- sd(unlist(twin.data[, c(selVars[2], selVars[4])]), na.rm = T)
  
  # Create Labels
  labMeMZ <- labVars("meanMZ", selVars)
  labMeDZ <- labVars("meanDZ", selVars)
  labMeZ <- labVars("meanZ", selVars)
  
  labCvMZ <- labLower("covMZ", ntv)
  labCvDZ <- labLower("covDZ", ntv)
  labCvZ <- labLower("covZ", ntv)
  
  labVaMZ <- labDiag("covMZ", ntv)
  labVaDZ <- labDiag("covDZ", ntv)
  labVaZ <- labDiag("covZ", ntv)
  
  # ------------------------------------------------------------------------------
  # PREPARE MODEL
  # Saturated Model
  
  # Create Matrices for Covariates and Linear Regression Coefficients
  ageLabels <- paste(rep(c("age_years_ses01", "age_years_ses02"), times = nt), rep(1:nt, each = nv), sep = "_0")
  sexLabels <- paste(rep(c("sex_female", "sex_female"), times = nt), rep(1:nt, each = nv), sep = "_0")
  outlierLabels <- paste(rep(c("fd_meanZ_ses01", "fd_meanZ_ses02"), times = nt), rep(1:nt, each = nv), sep = "_0")
  
  defAge <- mxMatrix(type = "Full", nrow = 1, ncol = nv * nt, free = F, labels = paste0(rep("data.", times = nv * nt), ageLabels), name = "defAge")
  defSex <- mxMatrix(type = "Full", nrow = 1, ncol = nv * nt, free = F, labels = paste0(rep("data.", times = nv * nt), sexLabels), name = "defSex")
  defMotion <- mxMatrix(type = "Full", nrow = 1, ncol = nv * nt, free = F, labels = paste0(rep("data.", times = nv * nt), outlierLabels), name = "defMotion")
  
  beta_age_mean_labels <- rep(paste0("beta_age_mean_coef_", 1:nv), times = nt)
  beta_sex_mean_labels <- rep(paste0("beta_sex_mean_coef_", 1:nv), times = nt)
  beta_outlier_mean_labels <- rep(paste0("beta_Motion_mean_coef_", 1:nv), times = nt)
  
  betaAge <- mxMatrix(type = "Full", nrow = 1, ncol = nt * nv, free = TRUE, values = svBe, labels = beta_age_mean_labels, name = "betaAge")
  betaSex <- mxMatrix(type = "Full", nrow = 1, ncol = nt * nv, free = TRUE, values = svBe, labels = beta_sex_mean_labels, name = "betaSex")
  betaMotion <- mxMatrix(type = "Full", nrow = 1, ncol = nt * nv, free = TRUE, values = svBe, labels = beta_outlier_mean_labels, name = "betaMotion")
  
  # Algebra for expected Mean Matrices in MZ & DZ twins
  meanMZ <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, labels = labMeMZ, name = "meanMZ")
  meanDZ <- mxMatrix(type = "Full", nrow = 1, ncol = ntv, free = TRUE, values = svMe, labels = labMeDZ, name = "meanDZ")
  
  expMeanMZ <- mxAlgebra(expression = meanMZ + betaAge * defAge + betaSex * defSex + betaMotion * defMotion, name = "expMeanMZ")
  expMeanDZ <- mxAlgebra(expression = meanDZ + betaAge * defAge + betaSex * defSex + betaMotion * defMotion, name = "expMeanDZ")
  
  # Twin Correlations
  expCorMZ <- mxAlgebra(cov2cor(expCovMZ), name = "expCorMZ")
  expCorDZ <- mxAlgebra(cov2cor(expCovDZ), name = "expCorDZ")
  
  # Create Algebra for expected Variance/Covariance Matrices
  expCovMZ <- mxMatrix(type = "Symm", nrow = ntv, free = T, values = valDiag(svVa, ntv), labels = labCvMZ, name = "expCovMZ")
  expCovDZ <- mxMatrix(type = "Symm", nrow = ntv, free = T, values = valDiag(svVa, ntv), labels = labCvDZ, name = "expCovDZ")
  
  # Data objects for Multiple Groups
  dataMZ <- mxData(observed = mzData, type = "raw")
  dataDZ <- mxData(observed = dzData, type = "raw")
  
  # Objective objects for Multiple Groups
  funML <- mxFitFunctionML()
  objMZ <- mxExpectationNormal(covariance = "expCovMZ", means = "expMeanMZ", dimnames = selVars)
  objDZ <- mxExpectationNormal(covariance = "expCovDZ", means = "expMeanDZ", dimnames = selVars)
  
  # Combine Groups
  modelMZ <- mxModel("MZ", meanMZ, betaAge, betaSex, betaMotion, defAge, defSex, defMotion, expMeanMZ, expCorMZ, expCovMZ, dataMZ, objMZ, funML)
  modelDZ <- mxModel("DZ", meanDZ, betaAge, betaSex, betaMotion, defAge, defSex, defMotion, expMeanDZ, expCorDZ, expCovDZ, dataDZ, objDZ, funML)
  minus2ll <- mxAlgebra(MZ.objective + DZ.objective, name = "minus2sumloglikelihood")
  obj <- mxFitFunctionAlgebra("minus2sumloglikelihood")
  ciCor <- mxCI(c("MZ.expCorMZ[3,1]", "DZ.expCorDZ[3,1]", "MZ.expCorMZ[4,2]", "DZ.expCorDZ[4,2]",
                  "MZ.expCorMZ[3,2]", "MZ.expCorMZ[4,1]", "DZ.expCorDZ[3,2]", "DZ.expCorDZ[4,1]"))
  twinSatModel <- mxModel("twinSat", modelMZ, modelDZ, minus2ll, obj, ciCor)
  
  # Run Saturated Model
  twinSatFit <- mxTryHard(twinSatModel, intervals = F, extraTries = 50)
  twinSatSumm <- summary(twinSatFit)
  
  #### Assumption testing ####
  # Equate means across twin order
  MModel1 <- twinSatModel
  MModel1 <- omxSetParameters(MModel1, label = c(labMeMZ[1], labMeMZ[3]), newlabels = "mMZ_ses01")
  MModel1 <- omxSetParameters(MModel1, label = c(labMeMZ[2], labMeMZ[4]), newlabels = "mMZ_ses02") 
  MModel1 <- omxSetParameters(MModel1, label = c(labMeDZ[1], labMeDZ[3]), newlabels = "mDZ_ses01")
  MModel1 <- omxSetParameters(MModel1, label = c(labMeDZ[2], labMeDZ[4]), newlabels = "mDZ_ses02")
  MModel1 <- omxAssignFirstParameters(MModel1)
  MModel1Fit <- mxTryHard(MModel1, intervals = F, extraTries = 50)
  
  # Equate means between zyg
  MModel2 <- MModel1
  MModel2 <- omxSetParameters(MModel2, label = c("mMZ_ses01", "mDZ_ses01"), newlabels = "Mean_ses01")
  MModel2 <- omxSetParameters(MModel2, label = c("mMZ_ses02", "mDZ_ses02"), newlabels = "Mean_ses02")
  MModel2 <- omxAssignFirstParameters(MModel2)
  MModel2Fit <- mxTryHard(MModel2, intervals = F, extraTries = 50)
  
  # Equate variance across twin order
  VModel1 <- MModel2
  VModel1 <- omxSetParameters(VModel1, label = c(labVaMZ[1], labVaMZ[3]), newlabels = "varMZ_ses01")
  VModel1 <- omxSetParameters(VModel1, label = c(labVaMZ[2], labVaMZ[4]), newlabels = "varMZ_ses02")
  VModel1 <- omxSetParameters(VModel1, label = c(labVaDZ[1], labVaDZ[3]), newlabels = "varDZ_ses01")
  VModel1 <- omxSetParameters(VModel1, label = c(labVaDZ[2], labVaDZ[4]), newlabels = "varDZ_ses02")
  VModel1 <- omxAssignFirstParameters(VModel1)
  VModel1Fit <- mxTryHard(VModel1, intervals = F, extraTries = 50)
  
  # Equate variance between zyg
  VModel2 <- VModel1
  VModel2 <- omxSetParameters(VModel2, label = c("varMZ_ses01", "varDZ_ses01"), newlabels = "varses01")
  VModel2 <- omxSetParameters(VModel2, label = c("varMZ_ses02", "varDZ_ses02"), newlabels = "varses02")
  VModel2 <- omxAssignFirstParameters(VModel2)
  VModel2Fit <- mxTryHard(VModel2, intervals = T, extraTries = 50)
  VModel2Summ <- summary(VModel2Fit)
  
  # Single CTCT correlation
  CTCT <- VModel2
  CTCT <- omxSetParameters(CTCT, label = c("covMZ41", "covMZ32"), newlabels = "covMZ_CTCT")
  CTCT <- omxSetParameters(CTCT, label = c("covDZ41", "covDZ32"), newlabels = "covDZ_CTCT")
  CTCT <- mxModel(CTCT, remove=TRUE, CTCT$intervals[])
  CTCT <- omxAssignFirstParameters(CTCT)
  ciCor.CTCT <- mxCI(c("MZ.expCorMZ[3,2]", "DZ.expCorDZ[3,2]"))
  CTCT <- mxModel(CTCT, ciCor.CTCT)
  CTCTFit <- mxTryHard(CTCT, intervals = T, extraTries = 50)
  CTCTSumm <- summary(CTCTFit)
  
  #### Phenotypic Correlation ses01-ses02
  rVModel <- VModel2
  rVModel <- omxSetParameters(rVModel, label = c("covMZ21", "covDZ21", "covMZ43", "covDZ43"), newlabels = "rV")
  rVModel <- mxModel(rVModel, remove=TRUE, rVModel$intervals[])
  rVModel <- omxAssignFirstParameters(rVModel)
  ciCor.rV <- mxCI("MZ.expCorMZ[2,1]")
  rVModel <- mxModel(rVModel, ciCor.rV)
  rVModelFit <- mxTryHard(rVModel, intervals = T, extraTries = 50)
  rVModelSum <- summary(rVModelFit)
  
  #### Test Covariates ####
  # Drop Sex
  noSex_ses01 <- VModel2
  noSex_ses01 <- omxSetParameters(model = noSex_ses01, labels = "beta_sex_mean_coef_1", free = FALSE, values = 0, name = "noSexModel_ses01")
  noSex_ses01Fit <- mxTryHard(noSex_ses01, extraTries = 100, intervals = F)
  noSex_ses01Summ <- summary(noSex_ses01Fit)
  
  noSex_ses02 <- VModel2
  noSex_ses02 <- omxSetParameters(model = noSex_ses02, labels = "beta_sex_mean_coef_2", free = FALSE, values = 0, name = "noSexModel_ses02")
  noSex_ses02Fit <- mxTryHard(noSex_ses02, extraTries = 100, intervals = F)
  noSex_ses02Summ <- summary(noSex_ses02Fit)
  
  # Drop Age
  noAge_ses01 <- VModel2
  noAge_ses01 <- omxSetParameters(model = noAge_ses01, labels = "beta_age_mean_coef_1", free = FALSE, values = 0, name = "noAgeModel_ses01")
  noAge_ses01Fit <- mxTryHard(noAge_ses01, extraTries = 100, intervals = F)
  noAge_ses01Summ <- summary(noAge_ses01Fit)
  
  noAge_ses02 <- VModel2
  noAge_ses02 <- omxSetParameters(model = noAge_ses02, labels = "beta_age_mean_coef_2", free = FALSE, values = 0, name = "noAgeModel_ses02")
  noAge_ses02Fit <- mxTryHard(noAge_ses02, extraTries = 100, intervals = F)
  noAge_ses02Summ <- summary(noAge_ses02Fit)
  
  # Drop Motion
  noMotion_ses01 <- VModel2
  noMotion_ses01 <- omxSetParameters(model = noMotion_ses01, labels = "beta_Motion_mean_coef_1", free = FALSE, values = 0, name = "noMotionModel_ses01")
  noMotion_ses01Fit <- mxTryHard(noMotion_ses01, extraTries = 100, intervals = F)
  noMotion_ses01Summ <- summary(noMotion_ses01Fit)
  
  noMotion_ses02 <- VModel2
  noMotion_ses02 <- omxSetParameters(model = noMotion_ses02, labels = "beta_Motion_mean_coef_2", free = FALSE, values = 0, name = "noMotionModel_ses02")
  noMotion_ses02Fit <- mxTryHard(noMotion_ses02, extraTries = 100, intervals = F)
  noMotion_ses02Summ <- summary(noMotion_ses02Fit)
  
  # Are ses01/ses02 means different?
  equalMeans <- VModel2
  equalMeans <- omxSetParameters(equalMeans, label = c("Mean_ses01", "Mean_ses02"), newlabels = "Mean")
  equalMeans <- omxAssignFirstParameters(equalMeans)
  equalMeansFit <- mxTryHard(equalMeans, extraTries = 100, intervals = F)
  equalMeansSumm <- summary(equalMeansFit)
  
  # LRT
  Drop_sex_ses01 <- mxCompare(VModel2Fit, noSex_ses01Fit)
  Drop_sex_ses02 <- mxCompare(VModel2Fit, noSex_ses02Fit)
  Drop_age_ses01 <- mxCompare(VModel2Fit, noAge_ses01Fit)
  Drop_age_ses02 <- mxCompare(VModel2Fit, noAge_ses02Fit)
  Drop_Motion_ses01 <- mxCompare(VModel2Fit, noMotion_ses01Fit)
  Drop_Motion_ses02 <- mxCompare(VModel2Fit, noMotion_ses02Fit)
  AssumpLRT1 <- mxCompare(twinSatFit, MModel1Fit)
  AssumpLRT2 <- mxCompare(MModel1Fit, MModel2Fit)
  AssumpLRT3 <- mxCompare(MModel2Fit, VModel1Fit)
  AssumpLRT4 <- mxCompare(VModel1Fit, VModel2Fit)
  SingleCTCT <- mxCompare(VModel2Fit, CTCTFit)
  EquateMeansLRT <- mxCompare(VModel2Fit, equalMeansFit)
  
  Results <- data.frame(
    Phenotype = phenotype[1],
    ObservedStatistics = twinSatSumm$observedStatistics,
    Mean_ses01 = mean.val.ses01,
    Mean_ses02 = mean.val.ses02,
    SD_ses01 = sd.val.ses01,
    SD_ses02 = sd.val.ses02,
    
    rMZ_ses01_95CI = paste0(sprintf("%.2f", round(VModel2Summ$CI$estimate[1], 2)), " (", sprintf("%.2f", round(VModel2Summ$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(VModel2Summ$CI$ubound[1], 2)), ")"),
    rDZ_ses01_95CI = paste0(sprintf("%.2f", round(VModel2Summ$CI$estimate[2], 2)), " (", sprintf("%.2f", round(VModel2Summ$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(VModel2Summ$CI$ubound[2], 2)), ")"),
    rMZ_ses02_95CI = paste0(sprintf("%.2f", round(VModel2Summ$CI$estimate[3], 2)), " (", sprintf("%.2f", round(VModel2Summ$CI$lbound[3], 2)), ", ", sprintf("%.2f", round(VModel2Summ$CI$ubound[3], 2)), ")"),
    rDZ_ses02_95CI = paste0(sprintf("%.2f", round(VModel2Summ$CI$estimate[4], 2)), " (", sprintf("%.2f", round(VModel2Summ$CI$lbound[4], 2)), ", ", sprintf("%.2f", round(VModel2Summ$CI$ubound[4], 2)), ")"),
    
    rMZ_Twin01_ses01_Twin02_ses02_95CI = paste0(sprintf("%.2f", round(VModel2Summ$CI$estimate[5], 2)), " (", sprintf("%.2f", round(VModel2Summ$CI$lbound[5], 2)), ", ", sprintf("%.2f", round(VModel2Summ$CI$ubound[5], 2)), ")"),
    rMZ_Twin02_ses02_Twin01_ses01_95CI = paste0(sprintf("%.2f", round(VModel2Summ$CI$estimate[6], 2)), " (", sprintf("%.2f", round(VModel2Summ$CI$lbound[6], 2)), ", ", sprintf("%.2f", round(VModel2Summ$CI$ubound[6], 2)), ")"),
    rDZ_Twin01_ses01_Twin02_ses02_95CI = paste0(sprintf("%.2f", round(VModel2Summ$CI$estimate[7], 2)), " (", sprintf("%.2f", round(VModel2Summ$CI$lbound[7], 2)), ", ", sprintf("%.2f", round(VModel2Summ$CI$ubound[7], 2)), ")"),
    rDZ_Twin02_ses02_Twin01_ses01_95CI = paste0(sprintf("%.2f", round(VModel2Summ$CI$estimate[8], 2)), " (", sprintf("%.2f", round(VModel2Summ$CI$lbound[8], 2)), ", ", sprintf("%.2f", round(VModel2Summ$CI$ubound[8], 2)), ")"),
    
    rMZ_CTCT_95CI = paste0(sprintf("%.2f", round(CTCTSumm$CI$estimate[1], 2)), " (", sprintf("%.2f", round(CTCTSumm$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(CTCTSumm$CI$ubound[1], 2)), ")"),
    rDZ_CTCT_95CI = paste0(sprintf("%.2f", round(CTCTSumm$CI$estimate[2], 2)), " (", sprintf("%.2f", round(CTCTSumm$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(CTCTSumm$CI$ubound[2], 2)), ")"),
    
    rV_ses01_ses02_95CI = paste0(sprintf("%.2f", round(rVModelSum$CI$estimate[1], 2)), " (", sprintf("%.2f", round(rVModelSum$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(rVModelSum$CI$ubound[1], 2)), ")"),
    
    age_ses01_Estimate = VModel2Summ$parameters$Estimate[3],
    age_ses02_Estimate = VModel2Summ$parameters$Estimate[4],
    sex_ses01_Estimate = VModel2Summ$parameters$Estimate[5],
    sex_ses02_Estimate = VModel2Summ$parameters$Estimate[6],
    Motion_ses01_Estimate = VModel2Summ$parameters$Estimate[7],
    Motion_ses02_Estimate = VModel2Summ$parameters$Estimate[8],
    sex_ses01_Pval = Drop_sex_ses01$p[2],
    sex_ses02_Pval = Drop_sex_ses02$p[2],
    age_ses01_Pval = Drop_age_ses01$p[2],
    age_ses02_Pval = Drop_age_ses02$p[2],
    Motion_ses01_Pval = Drop_Motion_ses01$p[2],
    Motion_ses02_Pval = Drop_Motion_ses02$p[2],
    assump1_pval = AssumpLRT1$p[2],
    assump2_pval = AssumpLRT2$p[2],
    assump3_pval = AssumpLRT3$p[2],
    assump4_pval = AssumpLRT4$p[2],
    ctct_single_pval = SingleCTCT$p[2],
    equateMeans_pval = EquateMeansLRT$p[2],
    stringsAsFactors = FALSE
  )
  return(Results)
}

#### Run models ####
qtab.data <- readRDS("QTAB_wide_waves.RDS")

# Run for single phenotype
Saturated_bivariate_covariate(phenotype = "MFG_R_7_3_reho", twin.data = qtab.data)

# Run for a list of phenotypes
variable_list <- readLines("BN_variables.txt")
reho <- paste0(variable_list, "_reho")
alff <- paste0(variable_list, "_alff")

results.reho <- as_tibble(lapply(reho, Saturated_bivariate_covariate, twin.data = qtab.data) %>% bind_rows())
results.alff <- as_tibble(lapply(alff, Saturated_bivariate_covariate, twin.data = qtab.data) %>% bind_rows())

#### Multiple testing correction ####
results.reho$sex_ses01_Pval_fdr <- p.adjust(results.reho$sex_ses01_Pval, method = "BH")
results.reho$sex_ses02_Pval_fdr <- p.adjust(results.reho$sex_ses02_Pval, method = "BH")

results.reho$age_ses01_Pval_fdr <- p.adjust(results.reho$age_ses01_Pval, method = "BH")
results.reho$age_ses02_Pval_fdr <- p.adjust(results.reho$age_ses02_Pval, method = "BH")

results.reho$Motion_ses01_Pval_fdr <- p.adjust(results.reho$Motion_ses01_Pval, method = "BH")
results.reho$Motion_ses02_Pval_fdr <- p.adjust(results.reho$Motion_ses02_Pval, method = "BH")

results.reho$assump1_pval_fdr <- p.adjust(results.reho$assump1_pval, method = "BH")
results.reho$assump2_pval_fdr <- p.adjust(results.reho$assump2_pval, method = "BH")
results.reho$assump3_pval_fdr <- p.adjust(results.reho$assump3_pval, method = "BH")
results.reho$assump4_pval_fdr <- p.adjust(results.reho$assump4_pval, method = "BH")

results.alff$sex_ses01_Pval_fdr <- p.adjust(results.alff$sex_ses01_Pval, method = "BH")
results.alff$sex_ses02_Pval_fdr <- p.adjust(results.alff$sex_ses02_Pval, method = "BH")

results.alff$age_ses01_Pval_fdr <- p.adjust(results.alff$age_ses01_Pval, method = "BH")
results.alff$age_ses02_Pval_fdr <- p.adjust(results.alff$age_ses02_Pval, method = "BH")

results.alff$Motion_ses01_Pval_fdr <- p.adjust(results.alff$Motion_ses01_Pval, method = "BH")
results.alff$Motion_ses02_Pval_fdr <- p.adjust(results.alff$Motion_ses02_Pval, method = "BH")

results.alff$assump1_pval_fdr <- p.adjust(results.alff$assump1_pval, method = "BH")
results.alff$assump2_pval_fdr <- p.adjust(results.alff$assump2_pval, method = "BH")
results.alff$assump3_pval_fdr <- p.adjust(results.alff$assump3_pval, method = "BH")
results.alff$assump4_pval_fdr <- p.adjust(results.alff$assump4_pval, method = "BH")

write.csv(results.reho, "saturated_bivariate_covariate_results_reho.csv", row.names = F)
write.csv(results.alff, "saturated_bivariate_covariate_results_alff.csv", row.names = F)
