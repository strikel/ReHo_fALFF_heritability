library(tidyverse)
library(lmerTest)
rm(list = ls())

setwd("~/GitHub/qtab/rsfMRI/CerebCortex_ms")

# SMFQ_score gives the boundary (singular fit) for many ReHo/ALFF ROIs
# Including it as a fixed effect seems to leave no variance for the participant_id random effect
# Results are identical w & w/out participant_id in these instances - only the singular fit warning message is different

qtab.data <- readRDS("QTAB_long_waves.RDS")
qtab.data <- 
  qtab.data %>% filter(!is.na(age_months))
qtab.data <- 
  qtab.data %>% filter(!is.na(fd_meanZ))
qtab.data$handedness_right <- ifelse(qtab.data$handedness=="Right", 1, 0)
qtab.data$handedness_left <- ifelse(qtab.data$handedness=="Left", 1, 0)
qtab.data$participant_id <- as.factor(qtab.data$participant_id_OpenNeuro)
qtab.data$Pair <- as.factor(qtab.data$Pair)
qtab.data$M <- as.factor(qtab.data$M)
qtab.data$age_years <- qtab.data$age_months/12
qtab.data$session_bin <- ifelse(qtab.data$session==1, 0, 1)

# Can centre variables to help interpretation of intercept
qtab.data$age_years_cent <- qtab.data$age_years - mean(qtab.data$age_years, na.rm = T)
qtab.data$SEI_cent <- qtab.data$SEI - mean(qtab.data$SEI, na.rm = T)

m1 <- lmer(SMFQ_score ~ 1 + age_years_cent + sex_female + fd_mean + session_bin + FuG_L_3_1_alff + (1 | Pair) + (1 | participant_id), data = qtab.data)
effectsize::standardize_parameters(m1, method = "basic")

#### lmer ####
for (iii in c("reho","alff")){
  for (ii in c("CrystallizedComposite_uncorr_stan","SCAS_score", "SPHERE_anxdep_score", "SMFQ_score")){
    mri.rois <- readLines(paste0("BN_variables.txt"))[1:210]
    results.df <- data.frame(matrix(nrow = length(mri.rois), ncol = 12))
    colnames(results.df) <- c("ROI", "error_message",
                              "age_beta", "age_pvalue", 
                              "sex_beta", "sex_pvalue",
                              "fd_beta", "fd_pvalue",
                              "session_beta", "session_pvalue",
                              "ROI_beta", "ROI_pvalue")
    reho.rois <- paste0(mri.rois, "_", iii)
    results.df$ROI <- reho.rois
    
    for (i in reho.rois){
      mod.tmp <- lmer(get(ii) ~ 1 + age_years_cent + sex_female + fd_mean + session_bin + get(i) + (1 | Pair) + (1 | participant_id), data = qtab.data)
      mod.summ <- summary(mod.tmp)
      results.df[which(results.df$ROI==i), c("age_beta", "age_pvalue")] <- c(mod.summ$coefficients[2,1], mod.summ$coefficients[2,5])
      results.df[which(results.df$ROI==i), c("sex_beta", "sex_pvalue")] <- c(mod.summ$coefficients[3,1], mod.summ$coefficients[3,5])
      results.df[which(results.df$ROI==i), c("fd_beta", "fd_pvalue")] <- c(mod.summ$coefficients[4,1], mod.summ$coefficients[4,5])
      results.df[which(results.df$ROI==i), c("session_beta", "session_pvalue")] <- c(mod.summ$coefficients[5,1], mod.summ$coefficients[5,5])
      results.df[which(results.df$ROI==i), c("ROI_beta", "ROI_pvalue")] <- c(mod.summ$coefficients[6,1], mod.summ$coefficients[6,5])
      error.message <- mod.tmp@optinfo$conv$lme4$messages
      if (!is.null(error.message)){
        results.df[which(results.df$ROI==i), "error_message"] <- error.message
      }
    }
    results.df$ROI_pvalue_fdr<- p.adjust(results.df$ROI_pvalue, method = "BH")
    write.table(results.df, paste0("01_lmer_", iii, "_", ii, ".txt"), sep = '\t', quote = F, col.names = T, row.names = F)
  }
}
