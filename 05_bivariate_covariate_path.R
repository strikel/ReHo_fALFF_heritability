rm(list = ls())
setwd("~/GitHub/qtab/rsfMRI/CerebCortex_ms")

#### Load libraries ####
library(OpenMx)
library(tidyverse)
library(umx)
source("miFunctions.R")
mxOption(NULL, "Default optimizer", "SLSQP") # SLSQP CSOLNP NPSOL
mxOption(NULL, key="Number of Threads", value=2)

Chol_path <- function(phenotype, twin.data){
  print(phenotype)
  phenotypes <- c(paste0(phenotype, "_reho_ses01"), paste0(phenotype, "_alff_ses01"))
  #phenotypes <- c(paste0(phenotype, "_alff_ses01"), paste0(phenotype, "_reho_ses01"))
  # Scale variables
  twin.data <- umx_scale_wide_twin_data(varsToScale = phenotypes, sep = "_0", data = twin.data)
  covariates = c("sex_female", "age_years_ses01", "fd_mean_ses01")
  # OpenMx does not tolerate missing values for definition variables.
  # Recode any missing definition variables as -999
  # BUT!!! Make sure there are not any cases of missing definition variables
  # with a phenotype present
  for (y in phenotypes){
    for (x in covariates) {
      twin01.missing <- twin.data[, paste0(y, "_01")][is.na(twin.data[, paste0(x, "_01")])]
      twin02.missing <- twin.data[, paste0(y, "_02")][is.na(twin.data[, paste0(x, "_02")])]
      stopifnot(is.na(twin01.missing))
      stopifnot(is.na(twin02.missing))
      twin.data[, paste0(x, "_01")][is.na(twin.data[, paste0(x, "_01")])] <- -999
      twin.data[, paste0(x, "_02")][is.na(twin.data[, paste0(x, "_02")])] <- -999
    }
  }
  nc <- length(covariates)
  
  # Select variables
  selVars <- c(paste0(phenotypes, "_01"), paste0(phenotypes, "_02"))
  covVars <- c(paste0(covariates, "_01"), paste0(covariates, "_02"))
  nv <- length(phenotypes)
  nt <- 2
  ntv <- nv * nt
  
  # Covariate labels
  ageLabels <- paste(rep(c("age_years_ses01"), times = nt), rep(1:nt, each = nv), sep = "_0")
  sexLabels <- paste(rep(c("sex_female"), times = nt), rep(1:nt, each = nv), sep = "_0")
  motionLabels <- paste(rep(c("fd_mean_ses01"), times = nt), rep(1:nt, each = nv), sep = "_0")
  
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
  ciModel <- mxCI(c("rV[2,1]", "rA[2,1]", "rC[2,1]", "rE[2,1]"))
  
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
  modelAE <- mxModel(modelAE, remove=TRUE, modelAE$intervals[['rC[2,1]']])
  fitAE <- mxTryHard(modelAE, intervals=FALSE, extraTries = 100)
  fitAE <- omxParallelCI(fitAE)
  sumAE <- summary(fitAE)
  
  # CE-CE
  modelCE <- mxModel(fitACE, name="CE_CE")
  modelCE <- omxSetParameters( modelCE, labels=c("a11", "a21", "a22"), free=FALSE, values=0)
  modelCE <- mxModel(modelCE, remove=TRUE, modelCE$intervals[['rA[2,1]']])
  fitCE <- mxTryHard(modelCE, intervals=FALSE, extraTries = 100)
  fitCE <- omxParallelCI(fitCE)
  sumCE <- summary(fitCE)
  
  # E-E
  modelE <- mxModel(fitCE, name="E_E")
  modelE <- omxSetParameters(modelE, labels=c("c11", "c21", "c22"), free=FALSE, values=0)
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
    ACE_rV = paste0(sprintf("%.2f", round(sumACE$CI$estimate[1], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[1], 2)), ")"),
    ACE_rA = paste0(sprintf("%.2f", round(sumACE$CI$estimate[2], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[2], 2)), ")"),
    ACE_rC = paste0(sprintf("%.2f", round(sumACE$CI$estimate[3], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[3], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[3], 2)), ")"),
    ACE_rE = paste0(sprintf("%.2f", round(sumACE$CI$estimate[4], 2)), " (", sprintf("%.2f", round(sumACE$CI$lbound[4], 2)), ", ", sprintf("%.2f", round(sumACE$CI$ubound[4], 2)), ")"),
    
    ## AE model ##
    AE_rV = paste0(sprintf("%.2f", round(sumAE$CI$estimate[1], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[1], 2)), ")"),
    AE_rA = paste0(sprintf("%.2f", round(sumAE$CI$estimate[2], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[2], 2)), ")"),
    AE_rE = paste0(sprintf("%.2f", round(sumAE$CI$estimate[3], 2)), " (", sprintf("%.2f", round(sumAE$CI$lbound[3], 2)), ", ", sprintf("%.2f", round(sumAE$CI$ubound[3], 2)), ")"),
    
    ## CE model ##
    CE_rV = paste0(sprintf("%.2f", round(sumCE$CI$estimate[1], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[1], 2)), ")"),
    CE_rC = paste0(sprintf("%.2f", round(sumCE$CI$estimate[2], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[2], 2)), ")"),
    CE_rE = paste0(sprintf("%.2f", round(sumCE$CI$estimate[3], 2)), " (", sprintf("%.2f", round(sumCE$CI$lbound[3], 2)), ", ", sprintf("%.2f", round(sumCE$CI$ubound[3], 2)), ")"),
    
    ## E model ##
    E_rV = paste0(sprintf("%.2f", round(sumE$CI$estimate[1], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[1], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[1], 2)), ")"),
    E_rE = paste0(sprintf("%.2f", round(sumE$CI$estimate[2], 2)), " (", sprintf("%.2f", round(sumE$CI$lbound[2], 2)), ", ", sprintf("%.2f", round(sumE$CI$ubound[2], 2)), ")")
  )
  return(Results)
}

qtab.data <- readRDS("QTAB_wide_waves.RDS")

# Run for a single phenotype
Chol_path(phenotype = "LOcC_L_4_4", twin.data = qtab.data)

# Run for list of phenotypes
variable_list <- readLines("BN_variables.txt")

results <- as_tibble(lapply(variable_list, Chol_path, twin.data = qtab.data) %>% bind_rows())

saveRDS(results, "ReHo_fALFF_ses01_correlations.RDS")
