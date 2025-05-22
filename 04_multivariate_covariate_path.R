rm(list = ls())
setwd("~/GitHub/qtab/rsfMRI/CerebCortex_ms")

#### Load libraries ####
library(OpenMx)
library(umx)
library(tidyverse)
library(parallel)
source("miFunctions.R")
mxOption(NULL, "Default optimizer", "NPSOL") # SLSQP CSOLNP NPSOL
mxOption(NULL, key="Number of Threads", value=2)

Chol_path <- function(phenotype1, phenotype2, twin.data){
  phenotypes <- c(paste0(phenotype1, "_ses01"), paste0(phenotype2, "_ses01"), paste0(phenotype1, "_ses02"), paste0(phenotype2, "_ses02"))
  # Scale variables
  twin.data <- umx_scale_wide_twin_data(varsToScale = phenotypes, sep = "_0", data = twin.data)
  # Regress out covariates for speed & to avoid errors with def vars
  twin.data <- umx_residualize(var = phenotypes[1], covs = c("sex_female", "age_years_ses01", "fd_mean_ses01"), suffixes = c("_01", "_02"), data = twin.data)
  twin.data <- umx_residualize(var = phenotypes[3], covs = c("sex_female", "age_years_ses02", "fd_mean_ses02"), suffixes = c("_01", "_02"), data = twin.data)
  twin.data <- umx_residualize(var = phenotypes[2], covs = c("sex_female", "age_years_ses01"), suffixes = c("_01", "_02"), data = twin.data)
  twin.data <- umx_residualize(var = phenotypes[4], covs = c("sex_female", "age_years_ses02"), suffixes = c("_01", "_02"), data = twin.data)
  
  selVars <- c(paste0(phenotypes, "_01"), paste0(phenotypes, "_02"))
  nv <- length(phenotypes)
  nt <- 2
  ntv <- nv * nt
  
  # Select Data for Analysis
  mzData <- subset(twin.data, zyg<=2, selVars)
  dzData <- subset(twin.data, zyg>=3, selVars)
  
  # Set Starting Values 
  svMe <- 0
  svPa <- 0.4
  svPe <- 0.6
  
  # ----------------------------------------------------------------------------------------------------------------------
  # PREPARE MODEL
  # ACE Model
  meanG <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labVars("mean", phenotypes), name="meanG")
  
  # Create Matrices for Path Coefficients
  pathA <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv), labels=labLower("a",nv), name="a" )
  pathC <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv), labels=labLower("c",nv), name="c" )
  pathE <- mxMatrix( type="Lower", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv), labels=labLower("e",nv), name="e" )
  
  # Create Algebra for Variance Components
  covA <- mxAlgebra( expression=a %*% t(a), name="A" )
  covC <- mxAlgebra( expression=c %*% t(c), name="C" )
  covE <- mxAlgebra( expression=e %*% t(e), name="E" )
  
  # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
  covP <- mxAlgebra( expression= A+C+E, name="V" )
  covMZ <- mxAlgebra( expression= A+C, name="cMZ" )
  covDZ <- mxAlgebra( expression= 0.5%x%A+ C, name="cDZ" )
  expCovMZ <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
  expCovDZ <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )
  
  # Create Data Objects for Multiple Groups
  dataMZ <- mxData( observed=mzData, type="raw" )
  dataDZ <- mxData( observed=dzData, type="raw" )
  
  # Create Expectation Objects for Multiple Groups
  expMZ <- mxExpectationNormal( covariance="expCovMZ", means="meanG", dimnames=selVars )
  expDZ <- mxExpectationNormal( covariance="expCovDZ", means="meanG", dimnames=selVars )
  funML <- mxFitFunctionML()
  
  # Create Algebra for Standardization
  matI <- mxMatrix(type="Iden", nrow=nv, ncol=nv, name="I")
  invSD <- mxAlgebra(expression=solve(sqrt(I*V)), name="iSD")
  
  # Calculate Genetic and Environmental Correlations
  corA <- mxAlgebra(expression=solve(sqrt(I*A))%&%A, name ="rA")
  corC <- mxAlgebra(expression=solve(sqrt(I*C))%&%C, name ="rC")
  corE <- mxAlgebra(expression=solve(sqrt(I*E))%&%E, name ="rE")
  corV <- mxAlgebra(expression=solve(sqrt(I*V))%&%V, name ="rV")
  
  # Create Algebra for Standardised Variance Components
  SA <- mxAlgebra(expression=A/V, name="SA")
  SC <- mxAlgebra(expression=C/V, name="SC")
  SE <- mxAlgebra(expression=E/V, name="SE")
  
  # Create Model Objects for Multiple Groups
  pars <- list(meanG, pathA, pathC, pathE, covA, covC, covE, covP, funML)
  modelMZ <- mxModel(pars, covMZ, expCovMZ, dataMZ, expMZ, name="MZ")
  modelDZ <- mxModel(pars, covDZ, expCovDZ, dataDZ, expDZ, name="DZ")
  multi <- mxFitFunctionMultigroup(c("MZ","DZ"))
  
  # Create Confidence Interval Objects
  ciModel <- mxCI(c("rV[2,1]", "rA[2,1]", "rC[2,1]", "rE[2,1]",
                    "rV[4,3]", "rA[4,3]", "rC[4,3]", "rE[4,3]"))
  
  # Build Model with Confidence Intervals
  calc <- list(SA, SC, SE, matI, invSD, corA, corC, corE, corV, ciModel)
  modelACE <- mxModel("ACE_ACE", pars, modelMZ, modelDZ, multi, calc)
  modelACE <- mxAutoStart(modelACE)
  # ----------------------------------------------------------------------------------------------------------------------
  
  # RUN MODEL
  # Run ACE Model
  fitACE <- mxTryHard(modelACE, intervals=FALSE, extraTries = 100)
  fitACE <- omxParallelCI(fitACE)
  sumACE <- summary(fitACE)

  # AE-AE
  modelAE <- mxModel(fitACE, name="AE_AE")
  modelAE <- omxSetParameters( modelAE, labels=labLower("c",nv), free=FALSE, values=0)
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['rC[2,1]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['rC[4,3]']])
  fitAE <- mxTryHard(modelAE, intervals=FALSE, extraTries = 100)
  fitAE <- omxParallelCI(fitAE)
  sumAE <- summary(fitAE)
  
  # CE-CE
  modelCE <- mxModel(fitACE, name="CE_CE")
  modelCE <- omxSetParameters( modelCE, labels=labLower("a",nv), free=FALSE, values=0)
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['rA[2,1]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['rA[4,3]']])
  fitCE <- mxTryHard(modelCE, intervals=FALSE, extraTries = 100)
  fitCE <- omxParallelCI(fitCE)
  sumCE <- summary(fitCE)
  
  # E-E
  modelE <- mxModel(fitCE, name="E_E")
  modelE <- omxSetParameters(modelE, labels=labLower("c",nv), free=FALSE, values=0)
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['rC[2,1]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['rC[4,3]']])
  fitE <- mxTryHard(modelE, intervals=FALSE, extraTries = 100)
  fitE <- omxParallelCI(fitE)
  sumE <- summary(fitE)

  model.comparison <- mxCompare(fitACE, c(fitAE, fitCE, fitE))
  
  Results <- data.frame(
    Phenotype = phenotype1,
    Scale = phenotype2,
    fitACE_AIC = model.comparison$AIC[1],
    fitAE_AIC = model.comparison$AIC[2],
    fitCE_AIC = model.comparison$AIC[3],
    fitE_AIC = model.comparison$AIC[4],

    ## ACE model ##
    ACE_rV_ses01 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[1], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[1], 2)), ")"),
    ACE_rA_ses01 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[2], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[2], 2)), ")"),
    ACE_rC_ses01 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[3], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[3], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[3], 2)), ")"),
    ACE_rE_ses01 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[4], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[4], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[4], 2)), ")"),
    ACE_rV_ses02 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[5], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[5], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[5], 2)), ")"),
    ACE_rA_ses02 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[6], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[6], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[6], 2)), ")"),
    ACE_rC_ses02 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[7], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[7], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[7], 2)), ")"),
    ACE_rE_ses02 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[8], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[8], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[8], 2)), ")"),
    
    ## AE model ##
    AE_rV_ses01 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[1], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[1], 2)), ")"),
    AE_rA_ses01 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[2], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[2], 2)), ")"),
    AE_rE_ses01 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[3], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[3], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[3], 2)), ")"),
    AE_rV_ses02 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[4], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[4], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[4], 2)), ")"),
    AE_rA_ses02 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[5], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[5], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[5], 2)), ")"),
    AE_rE_ses02 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[6], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[6], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[6], 2)), ")"),

    ## CE model ##
    CE_rV_ses01 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[1], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[1], 2)), ")"),
    CE_rC_ses01 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[2], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[2], 2)), ")"),
    CE_rE_ses01 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[3], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[3], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[3], 2)), ")"),
    CE_rV_ses02 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[4], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[4], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[4], 2)), ")"),
    CE_rC_ses02 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[5], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[5], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[5], 2)), ")"),
    CE_rE_ses02 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[6], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[6], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[6], 2)), ")"),

    ## E model ##
    E_rV_ses01 = paste0(sprintf("%.2f", round(sumE$CI$estimate[1], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[1], 2)), ")"),
    E_rE_ses01 = paste0(sprintf("%.2f", round(sumE$CI$estimate[2], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[2], 2)), ")"),
    E_rV_ses02 = paste0(sprintf("%.2f", round(sumE$CI$estimate[3], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[3], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[3], 2)), ")"),
    E_rE_ses02 = paste0(sprintf("%.2f", round(sumE$CI$estimate[4], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[4], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[4], 2)), ")")
      )
  return(Results)
}

qtab.data <- readRDS("QTAB_wide_waves.RDS")

# Run for a single phenotype
Chol_path(phenotype1 = "FuG_L_3_1_alff", phenotype2 = "SCAS_score", twin.data = qtab.data)
 
#### Posthoc tests ####
x1 <- Chol_path(phenotype1 = "ITG_R_7_1_reho", phenotype2 = "SCAS_score", twin.data = qtab.data)
x2 <- Chol_path(phenotype1 = "FuG_L_3_1_reho", phenotype2 = "SCAS_score", twin.data = qtab.data)
x3 <- Chol_path(phenotype1 = "PhG_R_6_1_reho", phenotype2 = "SCAS_score", twin.data = qtab.data)
x4 <- Chol_path(phenotype1 = "PhG_L_6_1_reho", phenotype2 = "SMFQ_score", twin.data = qtab.data)
x5 <- Chol_path(phenotype1 = "PhG_L_6_5_reho", phenotype2 = "SMFQ_score", twin.data = qtab.data)
x6 <- Chol_path(phenotype1 = "ITG_L_7_3_reho", phenotype2 = "SPHERE_anxdep_score", twin.data = qtab.data)
x7 <- Chol_path(phenotype1 = "FuG_L_3_1_reho", phenotype2 = "SPHERE_anxdep_score", twin.data = qtab.data)
x8 <- Chol_path(phenotype1 = "PhG_L_6_1_reho", phenotype2 = "SPHERE_anxdep_score", twin.data = qtab.data)
reho_results <- rbind(x1, x2, x3, x4, x5, x6, x7, x8)
reho_results$smallest_aic <- colnames(reho_results)[3:6][apply(reho_results[, 3:6],1,which.min)]
saveRDS(reho_results, "ReHo_scores_multi.RDS")

x1 <- Chol_path(phenotype1 = "OrG_R_6_2_alff", phenotype2 = "SCAS_score", twin.data = qtab.data)
x2 <- Chol_path(phenotype1 = "ITG_R_7_1_alff", phenotype2 = "SCAS_score", twin.data = qtab.data)
x3 <- Chol_path(phenotype1 = "ITG_L_7_3_alff", phenotype2 = "SCAS_score", twin.data = qtab.data)
x4 <- Chol_path(phenotype1 = "ITG_R_7_4_alff", phenotype2 = "SCAS_score", twin.data = qtab.data)
x5 <- Chol_path(phenotype1 = "FuG_L_3_1_alff", phenotype2 = "SCAS_score", twin.data = qtab.data)
x6 <- Chol_path(phenotype1 = "ITG_L_7_3_alff", phenotype2 = "SMFQ_score", twin.data = qtab.data)
x7 <- Chol_path(phenotype1 = "FuG_L_3_1_alff", phenotype2 = "SMFQ_score", twin.data = qtab.data)
x8 <- Chol_path(phenotype1 = "ITG_L_7_3_alff", phenotype2 = "SPHERE_anxdep_score", twin.data = qtab.data)
x9 <- Chol_path(phenotype1 = "FuG_L_3_1_alff", phenotype2 = "SPHERE_anxdep_score", twin.data = qtab.data)
falff_results <- rbind(x1, x2, x3, x4, x5, x6, x7, x8, x9)
falff_results$smallest_aic <- colnames(falff_results)[3:6][apply(falff_results[, 3:6],1,which.min)]
saveRDS(falff_results, "fALFF_scores_multi.RDS")
