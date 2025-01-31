rm(list = ls())
setwd("~/GitHub/ReHo_fALFF_heritability")

#### Load libraries ####
library(OpenMx)
library(tidyverse)
source("miFunctions.R")
mxOption(NULL, "Default optimizer", "NPSOL") # SLSQP CSOLNP NPSOL
mxOption(NULL, key="Number of Threads", value=8)

Chol_path <- function(phenotype, twin.data){
  print(phenotype)
  phenotypes <- c(paste0(phenotype, "_ses01"), paste0(phenotype, "_ses02"))
  covariates_ses01 = c("sex_female", "age_years_ses01", "fd_mean_ses01")
  covariates_ses02 = c("sex_female", "age_years_ses02", "fd_mean_ses02")
  # OpenMx does not tolerate missing values for definition variables.
  # Recode any missing definition variables as -999
  # BUT!!! Make sure there are not any cases of missing definition variables
  # with a phenotype present
  y <- phenotypes[1]
  for (x in covariates_ses01) {
    twin01.missing <- twin.data[, paste0(y, "_01")][is.na(twin.data[, paste0(x, "_01")])]
    twin02.missing <- twin.data[, paste0(y, "_02")][is.na(twin.data[, paste0(x, "_02")])]
    stopifnot(is.na(twin01.missing))
    stopifnot(is.na(twin02.missing))
    twin.data[, paste0(x, "_01")][is.na(twin.data[, paste0(x, "_01")])] <- -999
    twin.data[, paste0(x, "_02")][is.na(twin.data[, paste0(x, "_02")])] <- -999
  }
  
  y <- phenotypes[2]
  for (x in covariates_ses02) {
    twin01.missing <- twin.data[, paste0(y, "_01")][is.na(twin.data[, paste0(x, "_01")])]
    twin02.missing <- twin.data[, paste0(y, "_02")][is.na(twin.data[, paste0(x, "_02")])]
    stopifnot(is.na(twin01.missing))
    stopifnot(is.na(twin02.missing))
    twin.data[, paste0(x, "_01")][is.na(twin.data[, paste0(x, "_01")])] <- -999
    twin.data[, paste0(x, "_02")][is.na(twin.data[, paste0(x, "_02")])] <- -999
  }
  
  covariates = c("sex_female", "age_years_ses01", "fd_mean_ses01", "age_years_ses02", "fd_mean_ses02")
  nc <- length(covariates_ses01)
  
  # Select variables
  phenotypes <- c(paste0(phenotype, "_ses01"), paste0(phenotype, "_ses02"))
  selVars <- c(paste0(phenotypes, "_01"), paste0(phenotypes, "_02"))
  covVars <- c(paste0(covariates, "_01"), paste0(covariates, "_02"))
  nv <- length(phenotypes)
  nt <- 2
  ntv <- nv * nt
  
  # Covariate labels
  ageLabels <- paste(rep(c("age_years_ses01", "age_years_ses02"), times = nt), rep(1:nt, each = nv), sep = "_0")
  sexLabels <- paste(rep(c("sex_female", "sex_female"), times = nt), rep(1:nt, each = nv), sep = "_0")
  motionLabels <- paste(rep(c("fd_mean_ses01", "fd_mean_ses02"), times = nt), rep(1:nt, each = nv), sep = "_0")
  
  # Beta labels
  beta_age_mean_labels <- rep(paste0("beta_age_mean_coef_", 1:nv), times = nt)
  beta_sex_mean_labels <- rep(paste0("beta_sex_mean_coef_", 1:nv), times = nt)
  beta_motion_mean_labels <- rep(paste0("beta_motion_mean_coef_", 1:nv), times = nt)
  
  # Select Data for Analysis
  useVars <- c(selVars, covVars)
  mzData <- subset(twin.data, zyg<=2, useVars)
  dzData <- subset(twin.data, zyg>=3, useVars)
  
  # Set Starting Values 
  svMe <- as.numeric(colMeans(twin.data[, selVars], na.rm = T)[1:nv])
  svPa <- var(twin.data[, selVars], na.rm = T)[1]/2
  svPe <- var(twin.data[, selVars], na.rm = T)[1]
  svBe <- 0.01 # start value for regressions
  
  # ----------------------------------------------------------------------------------------------------------------------
  # PREPARE MODEL
  # ACE Model
  # Intercept Coefficient for Means
  mean <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labVars("mean", phenotypes), name="InterceptMean")
  
  # Create Betas for Covariates
  betaAge <- mxMatrix(type = "Full", nrow = 1, ncol = nt * nv, free = TRUE, values = svBe, labels = beta_age_mean_labels, name = "betaAge")
  betaSex <- mxMatrix(type = "Full", nrow = 1, ncol = nt * nv, free = TRUE, values = svBe, labels = beta_sex_mean_labels, name = "betaSex")
  betaMotion <- mxMatrix(type = "Full", nrow = 1, ncol = nt * nv, free = TRUE, values = svBe, labels = beta_motion_mean_labels, name = "betaMotion")
  
  # Create Definition Variables for Covariates
  defAge <- mxMatrix(type = "Full", nrow = 1, ncol = nv * nt, free = F, labels = paste0(rep("data.", times = nv * nt), ageLabels), name = "defAge")
  defSex <- mxMatrix(type = "Full", nrow = 1, ncol = nv * nt, free = F, labels = paste0(rep("data.", times = nv * nt), sexLabels), name = "defSex")
  defMotion <- mxMatrix(type = "Full", nrow = 1, ncol = nv * nt, free = F, labels = paste0(rep("data.", times = nv * nt), motionLabels), name = "defMotion")
  
  # Expected Means
  expMean <- mxAlgebra(expression = InterceptMean + betaAge * defAge + betaSex * defSex + betaMotion * defMotion, name = "expMean")
  
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
  expMZ <- mxExpectationNormal( covariance="expCovMZ", means="expMean", dimnames=selVars )
  expDZ <- mxExpectationNormal( covariance="expCovDZ", means="expMean", dimnames=selVars )
  funML <- mxFitFunctionML()
  
  # Create Algebra for Standardization
  matI <- mxMatrix(type="Iden", nrow=nv, ncol=nv, name="I")
  invSD <- mxAlgebra(expression=solve(sqrt(I*V)), name="iSD")
  
  # Create Algebra for Standardized Path Coefficients
  pathA_std <- mxAlgebra(expression=iSD %*% a, name="a_std")
  pathC_std <- mxAlgebra(expression=iSD %*% c, name="c_std")
  pathE_std <- mxAlgebra(expression=iSD %*% e, name="e_std")
  
  # Create Algebra for Standardized Path Coefficients^2
  pathA_std2 <- mxAlgebra(a_std%^%2, name="a_std2")
  pathC_std2 <- mxAlgebra(c_std%^%2, name="c_std2")
  pathE_std2 <- mxAlgebra(e_std%^%2, name="e_std2")

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
  pars <- list(pathA, pathC, pathE, covA, covC, covE, covP, funML)
  modelMZ <- mxModel(pars, mean, expMean, betaAge, betaSex, betaMotion, defAge, defSex, defMotion, covMZ, expCovMZ, dataMZ, expMZ, name="MZ")
  modelDZ <- mxModel(pars, mean, expMean, betaAge, betaSex, betaMotion, defAge, defSex, defMotion, covDZ, expCovDZ, dataDZ, expDZ, name="DZ")
  multi <- mxFitFunctionMultigroup(c("MZ","DZ"))
  
  # Create Confidence Interval Objects
  ciModel <- mxCI(c("A[1,1]", "A[2,2]",
                    "C[1,1]", "C[2,2]",
                    "E[1,1]", "E[2,2]",
                    "V[1,1]", "V[2,2]",
                    "SA[1,1]", "SA[2,2]",
                    "SC[1,1]", "SC[2,2]",
                    "SE[1,1]", "SE[2,2]",
                    "a[1,1]", "a[2,1]", "a[2,2]",
                    "c[1,1]", "c[2,1]", "c[2,2]",
                    "e[1,1]", "e[2,1]", "e[2,2]",
                    "a_std[1,1]", "a_std[2,1]", "a_std[2,2]",
                    "c_std[1,1]", "c_std[2,1]", "c_std[2,2]",
                    "e_std[1,1]", "e_std[2,1]", "e_std[2,2]",
                    "a_std2[1,1]", "a_std2[2,1]", "a_std2[2,2]",
                    "c_std2[1,1]", "c_std2[2,1]", "c_std2[2,2]",
                    "e_std2[1,1]", "e_std2[2,1]", "e_std2[2,2]",
                    "rV[2,1]", "rA[2,1]", "rC[2,1]", "rE[2,1]"))
  
  # Build Model with Confidence Intervals
  calc <- list( pathA_std, pathC_std, pathE_std, pathA_std2, pathC_std2, pathE_std2, SA, SC, SE, matI, invSD, corA, corC, corE, corV, ciModel)
  modelACE <- mxModel("ACE_ACE", pars, modelMZ, modelDZ, multi, calc)
  # ----------------------------------------------------------------------------------------------------------------------
  
  # RUN MODEL
  # Run ACE Model
  fitACE <- mxTryHard(modelACE, intervals=FALSE, extraTries = 100)
  fitACE <- omxParallelCI(fitACE)
  sumACE <- summary(fitACE)
  
  # AE-AE
  modelAE <- mxModel(fitACE, name="AE_AE")
  modelAE <- omxSetParameters( modelAE, labels=c("c11", "c21", "c22"), free=FALSE, values=0)
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['C[1,1]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['C[2,2]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['SC[1,1]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['SC[2,2]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['c[1,1]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['c[2,1]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['c[2,2]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['c_std[1,1]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['c_std[2,1]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['c_std[2,2]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['c_std2[1,1]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['c_std2[2,1]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['c_std2[2,2]']])
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['rC[2,1]']])
  fitAE <- mxTryHard(modelAE, intervals=FALSE, extraTries = 100)
  fitAE <- omxParallelCI(fitAE)
  sumAE <- summary(fitAE)
  
  # drop a22
  modelAE2 <- mxModel(fitAE, name="AE_AE")
  modelAE2 <- omxSetParameters( modelAE2, labels="a22", free=FALSE, values = 0)
  fitAE2 <- mxTryHard(modelAE2, intervals=FALSE, extraTries = 100)
  
  # a21=a11 (a22=0)
  modelAE3 <- mxModel(fitAE2, name="AE_AE")
  modelAE3 <- omxSetParameters( modelAE3, labels="a21", newlabels = "a11", free=TRUE)
  modelAE3 <- omxAssignFirstParameters(modelAE3)
  fitAE3 <- mxTryHard(modelAE3, intervals=FALSE, extraTries = 100)
  
  # a21=a11 (a22 free)
  modelAE4 <- mxModel(fitAE, name="AE_AE")
  modelAE4 <- omxSetParameters( modelAE4, labels="a21", newlabels = "a11", free=TRUE)
  modelAE4 <- omxAssignFirstParameters(modelAE4)
  fitAE4 <- mxTryHard(modelAE4, intervals=FALSE, extraTries = 100)
  
  drop_a22 <- mxCompare(fitAE, fitAE2)
  equate_a_wo_a22 <- mxCompare(fitAE2, fitAE3)
  equate_a_w_a22 <- mxCompare(fitAE, fitAE4)
  
  # CE-CE
  modelCE <- mxModel(fitACE, name="CE_CE")
  modelCE <- omxSetParameters( modelCE, labels=c("a11", "a21", "a22"), free=FALSE, values=0)
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['A[1,1]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['A[2,2]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['SA[1,1]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['SA[2,2]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['a[1,1]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['a[2,1]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['a[2,2]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['a_std[1,1]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['a_std[2,1]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['a_std[2,2]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['a_std2[1,1]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['a_std2[2,1]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['a_std2[2,2]']])
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['rA[2,1]']])
  fitCE <- mxTryHard(modelCE, intervals=FALSE, extraTries = 100)
  fitCE <- omxParallelCI(fitCE)
  sumCE <- summary(fitCE)
  
  # drop c22
  modelCE2 <- mxModel(fitCE, name="CE_CE")
  modelCE2 <- omxSetParameters( modelCE2, labels="c22", free=FALSE, values = 0)
  fitCE2 <- mxTryHard(modelCE2, intervals=FALSE, extraTries = 100)
  
  # c21=c11 (c22=0)
  modelCE3 <- mxModel(fitCE2, name="CE_CE")
  modelCE3 <- omxSetParameters( modelCE3, labels="c21", newlabels = "c11", free=TRUE)
  modelCE3 <- omxAssignFirstParameters(modelCE3)
  fitCE3 <- mxTryHard(modelCE3, intervals=FALSE, extraTries = 100)
  
  # c21=c11 (c22 free)
  modelCE4 <- mxModel(fitCE, name="CE_CE")
  modelCE4 <- omxSetParameters( modelCE4, labels="c21", newlabels = "c11", free=TRUE)
  modelCE4 <- omxAssignFirstParameters(modelCE4)
  fitCE4 <- mxTryHard(modelCE4, intervals=FALSE, extraTries = 100)
  
  drop_c22 <- mxCompare(fitCE, fitCE2)
  equate_c_wo_c22 <- mxCompare(fitCE2, fitCE3)
  equate_c_w_c22 <- mxCompare(fitCE, fitCE4)
  
  # E-E
  modelE <- mxModel(fitCE, name="E_E")
  modelE <- omxSetParameters(modelE, labels=c("c11", "c21", "c22"), free=FALSE, values=0)
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['C[1,1]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['C[2,2]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['SC[1,1]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['SC[2,2]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['c[1,1]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['c[2,1]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['c[2,2]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['c_std[1,1]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['c_std[2,1]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['c_std[2,2]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['c_std2[1,1]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['c_std2[2,1]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['c_std2[2,2]']])
  modelE <- mxModel(modelE, remove=TRUE, modelE$intervals[['rC[2,1]']])
  fitE <- mxTryHard(modelE, intervals=FALSE, extraTries = 100)
  fitE <- omxParallelCI(fitE)
  sumE <- summary(fitE)

  model.comparison <- mxCompare(fitACE, c(fitAE, fitCE, fitE))
  
  Results <- data.frame(
    Phenotype = phenotype,
    fitACE_AIC = model.comparison$AIC[1],
    fitAE_AIC = model.comparison$AIC[2],
    fitCE_AIC = model.comparison$AIC[3],
    fitE_AIC = model.comparison$AIC[4],

    ## ACE model ##
    ACE_A1 = paste0(sprintf("%.3f", round(sumACE$CI$estimate[1], 3)), " (", sprintf("%.3f", round(sumACE$CI$lbound[1], 3)), ", ", sprintf("%.3f", round(sumACE$CI$ubound[1], 3)), ")"),
    ACE_A2 = paste0(sprintf("%.3f", round(sumACE$CI$estimate[2], 3)), " (", sprintf("%.3f", round(sumACE$CI$lbound[2], 3)), ", ", sprintf("%.3f", round(sumACE$CI$ubound[2], 3)), ")"),
    ACE_C1 = paste0(sprintf("%.3f", round(sumACE$CI$estimate[3], 3)), " (", sprintf("%.3f", round(sumACE$CI$lbound[3], 3)), ", ", sprintf("%.3f", round(sumACE$CI$ubound[3], 3)), ")"),
    ACE_C2 = paste0(sprintf("%.3f", round(sumACE$CI$estimate[4], 3)), " (", sprintf("%.3f", round(sumACE$CI$lbound[4], 3)), ", ", sprintf("%.3f", round(sumACE$CI$ubound[4], 3)), ")"),
    ACE_E1 = paste0(sprintf("%.3f", round(sumACE$CI$estimate[5], 3)), " (", sprintf("%.3f", round(sumACE$CI$lbound[5], 3)), ", ", sprintf("%.3f", round(sumACE$CI$ubound[5], 3)), ")"),
    ACE_E2 = paste0(sprintf("%.3f", round(sumACE$CI$estimate[6], 3)), " (", sprintf("%.3f", round(sumACE$CI$lbound[6], 3)), ", ", sprintf("%.3f", round(sumACE$CI$ubound[6], 3)), ")"),
    ACE_V1 = paste0(sprintf("%.3f", round(sumACE$CI$estimate[7], 3)), " (", sprintf("%.3f", round(sumACE$CI$lbound[7], 3)), ", ", sprintf("%.3f", round(sumACE$CI$ubound[7], 3)), ")"),
    ACE_V2 = paste0(sprintf("%.3f", round(sumACE$CI$estimate[8], 3)), " (", sprintf("%.3f", round(sumACE$CI$lbound[8], 3)), ", ", sprintf("%.3f", round(sumACE$CI$ubound[8], 3)), ")"),
    ACE_SA1 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[9], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[9], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[9], 2)), ")"),
    ACE_SA2 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[10], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[10], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[10], 2)), ")"),
    ACE_SC1 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[11], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[11], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[11], 2)), ")"),
    ACE_SC2 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[12], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[12], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[12], 2)), ")"),
    ACE_SE1 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[13], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[13], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[13], 2)), ")"),
    ACE_SE2 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[14], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[14], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[14], 2)), ")"),
    
    ACE_a11 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[15], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[15], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[15], 2)), ")"),
    ACE_a21 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[16], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[16], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[16], 2)), ")"),
    ACE_a22 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[17], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[17], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[17], 2)), ")"),
    ACE_c11 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[18], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[18], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[18], 2)), ")"),
    ACE_c21 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[19], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[19], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[19], 2)), ")"),
    ACE_c22 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[20], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[20], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[20], 2)), ")"),
    ACE_e11 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[21], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[21], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[21], 2)), ")"),
    ACE_e21 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[22], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[22], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[22], 2)), ")"),
    ACE_e22 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[23], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[23], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[23], 2)), ")"),
    
    ACE_a11_std = paste0(sprintf("%.2f", round(sumACE$CI$estimate[24], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[24], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[24], 2)), ")"),
    ACE_a21_std = paste0(sprintf("%.2f", round(sumACE$CI$estimate[25], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[25], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[25], 2)), ")"),
    ACE_a22_std = paste0(sprintf("%.2f", round(sumACE$CI$estimate[26], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[26], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[26], 2)), ")"),
    ACE_c11_std = paste0(sprintf("%.2f", round(sumACE$CI$estimate[27], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[27], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[27], 2)), ")"),
    ACE_c21_std = paste0(sprintf("%.2f", round(sumACE$CI$estimate[28], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[28], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[28], 2)), ")"),
    ACE_c22_std = paste0(sprintf("%.2f", round(sumACE$CI$estimate[29], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[29], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[29], 2)), ")"),
    ACE_e11_std = paste0(sprintf("%.2f", round(sumACE$CI$estimate[30], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[30], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[30], 2)), ")"),
    ACE_e21_std = paste0(sprintf("%.2f", round(sumACE$CI$estimate[31], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[31], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[31], 2)), ")"),
    ACE_e22_std = paste0(sprintf("%.2f", round(sumACE$CI$estimate[32], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[32], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[32], 2)), ")"),
    
    ACE_a11_std2 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[33], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[33], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[33], 2)), ")"),
    ACE_a21_std2 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[34], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[34], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[34], 2)), ")"),
    ACE_a22_std2 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[35], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[35], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[35], 2)), ")"),
    ACE_c11_std2 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[36], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[36], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[36], 2)), ")"),
    ACE_c21_std2 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[37], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[37], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[37], 2)), ")"),
    ACE_c22_std2 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[38], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[38], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[38], 2)), ")"),
    ACE_e11_std2 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[39], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[39], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[39], 2)), ")"),
    ACE_e21_std2 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[40], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[40], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[40], 2)), ")"),
    ACE_e22_std2 = paste0(sprintf("%.2f", round(sumACE$CI$estimate[41], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[41], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[41], 2)), ")"),
    
    ACE_rV = paste0(sprintf("%.2f", round(sumACE$CI$estimate[42], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[42], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[42], 2)), ")"),
    ACE_rA = paste0(sprintf("%.2f", round(sumACE$CI$estimate[43], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[43], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[43], 2)), ")"),
    ACE_rC = paste0(sprintf("%.2f", round(sumACE$CI$estimate[44], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[44], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[44], 2)), ")"),
    ACE_rE = paste0(sprintf("%.2f", round(sumACE$CI$estimate[45], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[45], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[45], 2)), ")"),
    
    ## AE model ##
    AE_A1 = paste0(sprintf("%.3f", round(sumAE$CI$estimate[1], 3)), " (", sprintf("%.3f", round(sumAE$CI$lbound[1], 3)), ", ", sprintf("%.3f", round(sumAE$CI$ubound[1], 3)), ")"),
    AE_A2 = paste0(sprintf("%.3f", round(sumAE$CI$estimate[2], 3)), " (", sprintf("%.3f", round(sumAE$CI$lbound[2], 3)), ", ", sprintf("%.3f", round(sumAE$CI$ubound[2], 3)), ")"),
    AE_E1 = paste0(sprintf("%.3f", round(sumAE$CI$estimate[3], 3)), " (", sprintf("%.3f", round(sumAE$CI$lbound[3], 3)), ", ", sprintf("%.3f", round(sumAE$CI$ubound[3], 3)), ")"),
    AE_E2 = paste0(sprintf("%.3f", round(sumAE$CI$estimate[4], 3)), " (", sprintf("%.3f", round(sumAE$CI$lbound[4], 3)), ", ", sprintf("%.3f", round(sumAE$CI$ubound[4], 3)), ")"),
    AE_V1 = paste0(sprintf("%.3f", round(sumAE$CI$estimate[5], 3)), " (", sprintf("%.3f", round(sumAE$CI$lbound[5], 3)), ", ", sprintf("%.3f", round(sumAE$CI$ubound[5], 3)), ")"),
    AE_V2 = paste0(sprintf("%.3f", round(sumAE$CI$estimate[6], 3)), " (", sprintf("%.3f", round(sumAE$CI$lbound[6], 3)), ", ", sprintf("%.3f", round(sumAE$CI$ubound[6], 3)), ")"),
    AE_SA1 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[7], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[7], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[7], 2)), ")"),
    AE_SA2 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[8], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[8], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[8], 2)), ")"),
    AE_SE1 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[9], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[9], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[9], 2)), ")"),
    AE_SE2 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[10], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[10], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[10], 2)), ")"),
    
    AE_a11 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[11], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[11], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[11], 2)), ")"),
    AE_a21 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[12], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[12], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[12], 2)), ")"),
    AE_a22 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[13], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[13], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[13], 2)), ")"),
    AE_e11 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[14], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[14], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[14], 2)), ")"),
    AE_e21 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[15], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[15], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[15], 2)), ")"),
    AE_e22 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[16], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[16], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[16], 2)), ")"),
    
    AE_a11_std = paste0(sprintf("%.2f", round(sumAE$CI$estimate[17], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[17], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[17], 2)), ")"),
    AE_a21_std = paste0(sprintf("%.2f", round(sumAE$CI$estimate[18], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[18], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[18], 2)), ")"),
    AE_a22_std = paste0(sprintf("%.2f", round(sumAE$CI$estimate[19], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[19], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[19], 2)), ")"),
    AE_e11_std = paste0(sprintf("%.2f", round(sumAE$CI$estimate[20], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[20], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[20], 2)), ")"),
    AE_e21_std = paste0(sprintf("%.2f", round(sumAE$CI$estimate[21], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[21], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[21], 2)), ")"),
    AE_e22_std = paste0(sprintf("%.2f", round(sumAE$CI$estimate[22], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[22], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[22], 2)), ")"),
    
    AE_a11_std2 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[23], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[23], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[23], 2)), ")"),
    AE_a21_std2 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[24], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[24], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[24], 2)), ")"),
    AE_a22_std2 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[25], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[25], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[25], 2)), ")"),
    AE_e11_std2 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[26], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[26], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[26], 2)), ")"),
    AE_e21_std2 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[27], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[27], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[27], 2)), ")"),
    AE_e22_std2 = paste0(sprintf("%.2f", round(sumAE$CI$estimate[28], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[28], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[28], 2)), ")"),
    
    AE_rV = paste0(sprintf("%.2f", round(sumAE$CI$estimate[29], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[29], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[29], 2)), ")"),
    AE_rA = paste0(sprintf("%.2f", round(sumAE$CI$estimate[30], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[30], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[30], 2)), ")"),
    AE_rE = paste0(sprintf("%.2f", round(sumAE$CI$estimate[31], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[31], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[31], 2)), ")"),
  
    AE_drop_a22 = drop_a22$p[2],
    AE_equate_a_wo_a22 = equate_a_wo_a22$p[2],
    AE_equate_a_w_a22 = equate_a_w_a22$p[2],
    
    ## CE model ##
    CE_C1 = paste0(sprintf("%.3f", round(sumCE$CI$estimate[1], 3)), " (", sprintf("%.3f", round(sumCE$CI$lbound[1], 3)), ", ", sprintf("%.3f", round(sumCE$CI$ubound[1], 3)), ")"),
    CE_C2 = paste0(sprintf("%.3f", round(sumCE$CI$estimate[2], 3)), " (", sprintf("%.3f", round(sumCE$CI$lbound[2], 3)), ", ", sprintf("%.3f", round(sumCE$CI$ubound[2], 3)), ")"),
    CE_E1 = paste0(sprintf("%.3f", round(sumCE$CI$estimate[3], 3)), " (", sprintf("%.3f", round(sumCE$CI$lbound[3], 3)), ", ", sprintf("%.3f", round(sumCE$CI$ubound[3], 3)), ")"),
    CE_E2 = paste0(sprintf("%.3f", round(sumCE$CI$estimate[4], 3)), " (", sprintf("%.3f", round(sumCE$CI$lbound[4], 3)), ", ", sprintf("%.3f", round(sumCE$CI$ubound[4], 3)), ")"),
    CE_V1 = paste0(sprintf("%.3f", round(sumCE$CI$estimate[5], 3)), " (", sprintf("%.3f", round(sumCE$CI$lbound[5], 3)), ", ", sprintf("%.3f", round(sumCE$CI$ubound[5], 3)), ")"),
    CE_V2 = paste0(sprintf("%.3f", round(sumCE$CI$estimate[6], 3)), " (", sprintf("%.3f", round(sumCE$CI$lbound[6], 3)), ", ", sprintf("%.3f", round(sumCE$CI$ubound[6], 3)), ")"),
    CE_SC1 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[7], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[7], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[7], 2)), ")"),
    CE_SC2 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[8], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[8], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[8], 2)), ")"),
    CE_SE1 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[9], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[9], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[9], 2)), ")"),
    CE_SE2 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[10], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[10], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[10], 2)), ")"),
    
    CE_c11 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[11], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[11], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[11], 2)), ")"),
    CE_c21 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[12], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[12], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[12], 2)), ")"),
    CE_c22 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[13], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[13], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[13], 2)), ")"),
    CE_e11 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[14], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[14], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[14], 2)), ")"),
    CE_e21 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[15], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[15], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[15], 2)), ")"),
    CE_e22 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[16], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[16], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[16], 2)), ")"),
    
    CE_c11_std = paste0(sprintf("%.2f", round(sumCE$CI$estimate[17], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[17], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[17], 2)), ")"),
    CE_c21_std = paste0(sprintf("%.2f", round(sumCE$CI$estimate[18], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[18], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[18], 2)), ")"),
    CE_c22_std = paste0(sprintf("%.2f", round(sumCE$CI$estimate[19], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[19], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[19], 2)), ")"),
    CE_e11_std = paste0(sprintf("%.2f", round(sumCE$CI$estimate[20], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[20], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[20], 2)), ")"),
    CE_e21_std = paste0(sprintf("%.2f", round(sumCE$CI$estimate[21], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[21], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[21], 2)), ")"),
    CE_e22_std = paste0(sprintf("%.2f", round(sumCE$CI$estimate[22], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[22], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[22], 2)), ")"),
    
    CE_c11_std2 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[23], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[23], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[23], 2)), ")"),
    CE_c21_std2 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[24], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[24], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[24], 2)), ")"),
    CE_c22_std2 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[25], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[25], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[25], 2)), ")"),
    CE_e11_std2 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[26], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[26], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[26], 2)), ")"),
    CE_e21_std2 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[27], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[27], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[27], 2)), ")"),
    CE_e22_std2 = paste0(sprintf("%.2f", round(sumCE$CI$estimate[28], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[28], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[28], 2)), ")"),
    
    CE_rV = paste0(sprintf("%.2f", round(sumCE$CI$estimate[29], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[29], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[29], 2)), ")"),
    CE_rC = paste0(sprintf("%.2f", round(sumCE$CI$estimate[30], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[30], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[30], 2)), ")"),
    CE_rE = paste0(sprintf("%.2f", round(sumCE$CI$estimate[31], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[31], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[31], 2)), ")"),
    
    CE_drop_c22 = drop_c22$p[2],
    CE_equate_c_wo_c22 = equate_c_wo_c22$p[2],
    CE_equate_c_w_c22 = equate_c_w_c22$p[2],
    
    ## E model ##
    E_E1 = paste0(sprintf("%.3f", round(sumE$CI$estimate[1], 3)), " (", sprintf("%.3f", round(sumE$CI$lbound[1], 3)), ", ", sprintf("%.3f", round(sumE$CI$ubound[1], 3)), ")"),
    E_E2 = paste0(sprintf("%.3f", round(sumE$CI$estimate[2], 3)), " (", sprintf("%.3f", round(sumE$CI$lbound[2], 3)), ", ", sprintf("%.3f", round(sumE$CI$ubound[2], 3)), ")"),
    E_V1 = paste0(sprintf("%.3f", round(sumE$CI$estimate[3], 3)), " (", sprintf("%.3f", round(sumE$CI$lbound[3], 3)), ", ", sprintf("%.3f", round(sumE$CI$ubound[3], 3)), ")"),
    E_V2 = paste0(sprintf("%.3f", round(sumE$CI$estimate[4], 3)), " (", sprintf("%.3f", round(sumE$CI$lbound[4], 3)), ", ", sprintf("%.3f", round(sumE$CI$ubound[4], 3)), ")"),
    E_SE1 = paste0(sprintf("%.2f", round(sumE$CI$estimate[5], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[5], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[5], 2)), ")"),
    E_SE2 = paste0(sprintf("%.2f", round(sumE$CI$estimate[6], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[6], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[6], 2)), ")"),
    
    E_e11 = paste0(sprintf("%.2f", round(sumE$CI$estimate[7], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[7], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[7], 2)), ")"),
    E_e21 = paste0(sprintf("%.2f", round(sumE$CI$estimate[8], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[8], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[8], 2)), ")"),
    E_e22 = paste0(sprintf("%.2f", round(sumE$CI$estimate[9], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[9], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[9], 2)), ")"),
    
    E_e11_std = paste0(sprintf("%.2f", round(sumE$CI$estimate[10], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[10], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[10], 2)), ")"),
    E_e21_std = paste0(sprintf("%.2f", round(sumE$CI$estimate[11], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[11], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[11], 2)), ")"),
    E_e22_std = paste0(sprintf("%.2f", round(sumE$CI$estimate[12], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[12], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[12], 2)), ")"),
    
    E_e11_std2 = paste0(sprintf("%.2f", round(sumE$CI$estimate[13], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[13], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[13], 2)), ")"),
    E_e21_std2 = paste0(sprintf("%.2f", round(sumE$CI$estimate[14], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[14], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[14], 2)), ")"),
    E_e22_std2 = paste0(sprintf("%.2f", round(sumE$CI$estimate[15], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[15], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[15], 2)), ")"),
    
    E_rV = paste0(sprintf("%.2f", round(sumE$CI$estimate[16], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[16], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[16], 2)), ")"),
    E_rE = paste0(sprintf("%.2f", round(sumE$CI$estimate[17], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[17], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[17], 2)), ")")
    
  )
  write.csv(Results, paste0("bivariate_covariate_", phenotype, "_result.csv"), row.names = F)
}

qtab.data <- readRDS("QTAB_wide_waves.RDS")

# Run for a single phenotype
Chol_path(phenotype = "OrG_R_6_4_reho", twin.data = qtab.data)

# Run for list of phenotypes
variable_list <- readLines("BN_variables.txt")[1:210]
reho <- paste0(variable_list, "_reho")
alff <- paste0(variable_list, "_alff")

lapply(reho, Chol_path, twin.data = qtab.data)
lapply(alff, Chol_path, twin.data = qtab.data)
