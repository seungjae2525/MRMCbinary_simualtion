################################################################################
################################################################################
## Load packages
library(stringr)
library(openxlsx)
library(ggplot2)
library(patchwork)

## Automatically set working directory to the location of this script
# setwd("~/")
args <- commandArgs(trailingOnly = FALSE)
print("Command arguments:")
print(args)
script_path <- NULL
file_arg_index <- grep("--file=", args)
if (length(file_arg_index) > 0) {
  filepath <- sub("^--file=", "", args[file_arg_index])
  script_path <- normalizePath(filepath)
}
if (is.null(script_path)) {
  if (!is.null(sys.frames()) && !is.null(sys.frames()[[1]]$ofile)) {
    script_path <- normalizePath(sys.frames()[[1]]$ofile)
  }
}
if (is.null(script_path) && requireNamespace("rstudioapi", quietly = TRUE) && 
    rstudioapi::isAvailable()) {
  script_path <- rstudioapi::getActiveDocumentContext()$path
}
if (!is.null(script_path) && file.exists(script_path)) {
  setwd(dirname(script_path))
  print(paste("Working directory set to:", getwd()))
} else {
  stop("Unable to set working directory automatically.")
}

## Load all simulation results
rm(list = ls())
load(file = "Simulation in Sec 2.3.RData")

##
re_function <- function(results_all = results.sim, n.par = 1, xxx="fit_CLR") {
  
  xxx2 <- xxx
  xxx <- strsplit(xxx, "_s")[[1]]
  
  model <- ifelse(grepl("CLR", xxx), "CLR",
                  ifelse(grepl("GLMM", xxx), "GLMM", 
                         ifelse(grepl("GEE", xxx), "GEE", "GLM")))
  
  results <- results_all[[n.par]]
  I <- par.comb[n.par, 1]
  n_J <- par.comb[n.par, 2]
  fit_results <- purrr::map(results, ~ .x[[xxx]])
  ##############################################################################
  ##
  ew_message <- unlist(purrr::map(results, ~ .x[[paste0(xxx,"_message")]]))
  ew_message <- ifelse(grepl(pattern = "Model failed to converge", ew_message), 
                       "WWWW Model failed to converge with max|grad|", ew_message)
  
  ## Error
  e_which <- grepl("EEEE", ew_message)
  ## GEE negative variance error
  if (model == "GEE") {
    gee_e_which <- unlist(purrr::map(results, ~ .x[[paste0(xxx,"_e_which")]]))
    e_which <- e_which | gee_e_which
  } 
  
  ## Additional Error: |Estimate| > 100 --> Error!
  strange_se <- sapply(seq_along(fit_results), function(i) {
    summary_x <- fit_results[[i]]
    if (sum(is.na(summary_x)) > 0) {
      est <- NA
      se <- NA
    } else if (model == "CLR") {
      est <- summary_x[, 1]
      se <- summary_x[, 3]
    } else if (model %in% c("GLM", "GEE", "GLMM")) {
      est <- summary_x[, 1]
      se <- summary_x[, 2]
    } else {
      stop("")
    }
    temp1 <- any(is.nan(est)) | any(is.nan(se)) | any(is.na(est)) | any(is.na(se))
    temp2 <- ifelse(temp1 == TRUE, FALSE, any(abs(est) > 100))
    c(temp1, temp2)
  }) 
  nan_which <- strange_se[1,]
  big_est_which <- strange_se[2,]
  
  ## Total error
  e_which <- e_which | nan_which | big_est_which
  
  ## Total error message
  e_message <- ew_message[e_which]
  if (length(e_message) == 0) {
    e_message <- "No error"
  } else {
    e_message <- table(e_message)
  }
  
  ## Warning
  w_which <- grepl("WWWW", ew_message)
  
  ## Total warning message
  w_message <- ew_message[w_which]
  if (length(w_message) == 0) {
    w_message <- "No warning"
  } else {
    w_message <- table(w_message)
  }
  
  ## Both error and warning
  ew_which <- e_which | w_which
  else_which <- grepl("No warnings or errors", ew_message)
  
  ##############################################################################
  ## For both errors and warnings
  fit_results_ew <- lapply(seq_along(fit_results), function(i) {
    if (!ew_which[i]) fit_results[[i]] else NULL
  })
  fit_results_ew_remove <- fit_results_ew
  fit_results_ew_remove[sapply(fit_results_ew_remove, is.null)] <- NULL
  
  x <- fit_results_ew_remove
  
  ##############################################################################
  if (model == "CLR") {
    results_temp <- sapply(1:length(x), function(i) {
      summary_x <- x[[i]]
      
      est <- summary_x[, 1]
      se <- summary_x[, 3]
      
      lb <- est - qnorm(0.975) * se
      ub <- est + qnorm(0.975) * se
      
      ll1 <- cbind(est, se, lb, ub)
      ll1 <- data.frame(ll1)
      names(ll1) <- c("coef", "se (coef)", "95CI (Lower)", "95CI (Upper)")
      list(dplyr::lst(ll1))
    })
  } else if (model == "GLM") {
    if (xxx2 == "fit_GLM_s") {
      results_temp <- sapply(1:length(x), function(i) {
        summary_x <- x[[i]]
        
        est <- summary_x[, 1] / (phi * lambda)
        se <- summary_x[, 2] / (phi * lambda)
        
        lb <- est - qnorm(0.975) * se
        ub <- est + qnorm(0.975) * se
        
        ll1 <- cbind(est, se, lb, ub)[-1,]
        ll1 <- data.frame(ll1)
        names(ll1) <- c("coef", "se (coef)", "95CI (Lower)", "95CI (Upper)")
        list(dplyr::lst(ll1))
      })
    } else {
      results_temp <- sapply(1:length(x), function(i) {
        summary_x <- x[[i]]
        
        est <- summary_x[, 1]
        se <- summary_x[, 2]
        
        lb <- est - qnorm(0.975) * se
        ub <- est + qnorm(0.975) * se
        
        ll1 <- cbind(est, se, lb, ub)[-1,]
        ll1 <- data.frame(ll1)
        names(ll1) <- c("coef", "se (coef)", "95CI (Lower)", "95CI (Upper)")
        list(dplyr::lst(ll1))
      }) 
    }
  } else if (model == "GEE") {
    if (xxx2 %in% c("fit_GEE1_s", "fit_GEE2_s")) {
      results_temp <- sapply(1:length(x), function(i) {
        summary_x <- x[[i]]
        
        est <- summary_x[, 1] / (phi * lambda)
        se <- summary_x[, 2] / (phi * lambda)
        
        lb <- est - qnorm(0.975) * se
        ub <- est + qnorm(0.975) * se
        
        ll1 <- cbind(est, se, lb, ub)[-1,]
        ll1 <- data.frame(ll1)
        names(ll1) <- c("coef", "se (coef)", "95CI (Lower)", "95CI (Upper)")
        list(dplyr::lst(ll1))
      })
    } else {
      results_temp <- sapply(1:length(x), function(i) {
        summary_x <- x[[i]]
        
        est <- summary_x[, 1]
        se <- summary_x[, 2]
        
        lb <- est - qnorm(0.975) * se
        ub <- est + qnorm(0.975) * se
        
        ll1 <- cbind(est, se, lb, ub)[-1,]
        ll1 <- data.frame(ll1)
        names(ll1) <- c("coef", "se (coef)", "95CI (Lower)", "95CI (Upper)")
        list(dplyr::lst(ll1))
      })
    }
  } else if (model == "GLMM") {
    results_temp <- sapply(1:length(x), function(i) {
      summary_x <- x[[i]]
      
      est <- summary_x[, "Estimate"]
      se <- summary_x[, "Std. Error"]
      
      lb <- est - qnorm(0.975) * se
      ub <- est + qnorm(0.975) * se
      
      ll1 <- data.frame(cbind(est[-1], se[-1], lb[-1], ub[-1]))
      names(ll1) <- c("coef", "se (coef)", "95CI (Lower)", "95CI (Upper)")
      list(dplyr::lst(ll1))
    })
  }
  
  ## Extract ll1 elements
  re1 <- purrr::map(results_temp, ~ .x$ll1)
  
  ## Get the number of rows in the first element of re1 to determine
  n_rows <- nrow(re1[[1]])
  
  ## Populate the mats list
  mats <- lapply(1:n_rows, function(j) {
    do.call(rbind, lapply(re1, function(x) x[j, ]))
  })
  
  ###
  ## Initialize vectors to store results
  est <- bias <- se_mean <- se_median <- lower_ci <- upper_ci <- numeric(n_rows)
  coverage <- b_analysis <- warn <- numeric(n_rows)
  for (ii in 1:n_rows) {
    ## Modality effects
    all_effect <- cbind(beta_2k, beta_3k)
    mat_new <- mats[[ii]]
    
    ## Calculate estimates, standard errors, confidence intervals, and coverage
    est[ii] <- mean(mat_new[, "coef"])
    bias[ii] <- mean(mat_new[, "coef"] - all_effect[,ii])
    se_mean[ii] <- mean(mat_new[, "se (coef)"])
    se_median[ii] <- median(mat_new[, "se (coef)"])
    lower_ci[ii] <- median(mat_new[, "95CI (Lower)"])
    upper_ci[ii] <- median(mat_new[, "95CI (Upper)"])
    coverage[ii] <- mean(mat_new[, "95CI (Lower)"] <= all_effect[,ii] & 
                           all_effect[,ii] <= mat_new[, "95CI (Upper)"])
    # warn[ii] <- length(warn_which)
    b_analysis[ii] <- nrow(all_effect)
  }
  
  ## Create a data frame with the results
  datt1 <- data.frame(est = est,
                      bias = bias,
                      se_mean = se_mean,
                      se_median = se_median,
                      lower_ci = lower_ci,
                      upper_ci = upper_ci,
                      coverage = coverage,
                      error_gee = ifelse(model == "gee", sum(gee_e_which), 0), # GEE Error
                      error_nan = sum(nan_which), # NaN estimate Error
                      error_est = sum(big_est_which), # Large estimate Error
                      error_total = sum(e_which), # Total Error 
                      warning_total = sum(w_which), # Total Warning
                      error_warning_total = sum(ew_which), # Total Error + Warning
                      else_ew = sum(else_which), # Error or Warning in model output message
                      b_analysis = length(re1)) # Number used in analysis
  return(list(datt1 = datt1, e_message = e_message, w_message = w_message))
}

##
all_method <- c("fit_CLR", 
                "fit_GLM", "fit_GLM_s",  
                "fit_GEE1", "fit_GEE1_s",  
                "fit_GEE2", "fit_GEE2_s",
                "fit_glmer1", "fit_glmer2",
                "fit_glmmTMB1", "fit_glmmTMB2")
all_results <- c()
for (iii in 1:nrow(par.comb)) {
  all_results_temp <- c()
  for (kkk in 1:length(all_method)) {
    Re_temp <- re_function(results_all = results.sim, n.par = iii, xxx = all_method[kkk])
    Re_temp <- Re_temp$datt1
    Re_temp$xxx <- all_method[kkk]
    Re_temp$beta <- c("beta_2k", "beta_3k")
    Re_temp$I <- par.comb[iii, 1]
    Re_temp$n_J <- par.comb[iii, 2]
    all_results_temp <- rbind(all_results_temp, Re_temp)
    rm(Re_temp)
  }
  all_results <- rbind(all_results, all_results_temp)
  rm(all_results_temp)
}


###
reprot_f <- function(dat = all_results, B = 1000) {
  
  n_K_1 <- length(unique(all_results$beta))
  n_K <- n_K_1+1
  
  ##
  Parameter <- ifelse(all_results$beta == "beta_2k", 
                      "$\\beta_{2, z}^{(1)}$", "$\\beta_{3, z}^{(1)}$")
  
  ##
  Method <- stringr::str_split(all_results$xxx, "_", simplify = T)
  Method <- ifelse(Method[,3] == "", Method[,2],
                   paste0(Method[,2], "_", Method[,3]))
  
  Method <- gsub("CLR", "CLR", Method)
  Method <- gsub("GLMM", "GLMM", Method)
  Method <- gsub("GLM", "GLM", Method)
  Method <- gsub("GEE", "GEE", Method)
  
  n_I <- all_results$I
  n_J <- all_results$n_J
  
  Estimate <- sprintf('%.3f', round(all_results$est, 3))
  Bias <- sprintf('%.3f', round(all_results$bias, 3))
  SE <- sprintf('%.3f', round(all_results$se_mean, 3))
  CI <- paste0("[", sprintf('%.3f', round(all_results$lower_ci, 3)), ", ",
               sprintf('%.3f', round(all_results$upper_ci, 3)), "]")
  Coverage <- sprintf('%.3f', round(all_results$coverage, 3))
  Removed_temp <- B - all_results$b_analysis
  Removed <- paste0(Removed_temp, " (", 
                    sprintf('%.2f', round(100*Removed_temp/B, 2)), "%)")
  Removed <- ifelse(Removed == "0 (0.00%)", "0 (0%)", Removed)
  
  return_dat <- data.frame(Method, "Removed cases" = Removed, 
                           n_I, n_J, 
                           Parameter, Estimate, Bias, SE, 
                           "Confidence interval" = CI, Coverage, 
                           check.names = FALSE)
  return(return_dat)
}

report_table <- reprot_f(dat = all_results, B = 1000)

##
report_table$Method <- factor(report_table$Method, levels = unique(report_table$Method))
report_table <- report_table[order(report_table$n_I, report_table$Method), ]

## xlsx save
wb <- createWorkbook()
addWorksheet(wb, "Sheet1")
writeData(wb, "Sheet1", report_table, startRow = 1, startCol = 1, rowNames = FALSE)
saveWorkbook(wb, file = "Simulation in Sec 2.3.xlsx", overwrite = TRUE)


################################################################################
################################################################################
## Results for n = 100
temp_name <- "Simulation in Sec 2.3"

all_results_n100 <- all_results[all_results$I == 100, ]
plot_modality_n100 <- data.frame(
  Method = rep(ifelse(stringr::str_split(all_results$xxx, "_", simplify = T)[,3] == "", 
                      stringr::str_split(all_results$xxx, "_", simplify = T)[,2],
                      paste0(stringr::str_split(all_results$xxx, "_", simplify = T)[,2], "*")), 4),
  Metric = c(rep("Absolute bias", nrow(all_results_n100)), 
             rep("Standard error", nrow(all_results_n100)),
             rep("Coverage", nrow(all_results_n100)),
             rep("Proportion of failed", nrow(all_results_n100))),
  Target = rep(all_results_n100$n_J, 4), 
  Parameter = rep(all_results_n100$beta, 4),
  Value = c(abs(as.numeric(all_results_n100$bias)), 
            as.numeric(all_results_n100$se_mean), 
            as.numeric(all_results_n100$coverage), 
            as.numeric((1000 - all_results_n100$b_analysis)/1000)))

#
plot_modality_n100$Method <- gsub("CLR", "CLR", plot_modality_n100$Method)
plot_modality_n100$Method <- gsub("GLMM", "GLMM", plot_modality_n100$Method)
plot_modality_n100$Method <- gsub("GLM", "GLM", plot_modality_n100$Method)
plot_modality_n100$Method <- gsub("GEE", "GEE", plot_modality_n100$Method)

#
plot_modality_n100$Target <- ifelse(plot_modality_n100$Target == 3, "3 Readers",
                                    ifelse(plot_modality_n100$Target == 5, "5 Readers", "10 Readers"))

#
plot_modality_n100$Target <- factor(plot_modality_n100$Target, 
                                    levels = c("3 Readers", "5 Readers", "10 Readers"))

#
plot_modality_n100$Method <- factor(plot_modality_n100$Method, 
                                    levels = c("CLR", 
                                               "GLM", "GLM*", 
                                               "GEE1", "GEE1*", 
                                               "GEE2", "GEE2*", 
                                               "glmer1", "glmer2", "glmmTMB1", "glmmTMB2"))


################################################################################
## Plot for beta2
library(ggplot2)
library(patchwork)
plot_modality_n100_beta2 <- plot_modality_n100[plot_modality_n100$Parameter == "beta_2k", ]
plot_A <- ggplot(subset(plot_modality_n100_beta2, Metric == "Absolute bias"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 0.7, 0.05),
                     limits = c(0, 0.7)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Absolute bias") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

plot_B <- ggplot(subset(plot_modality_n100_beta2, Metric == "Standard error"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 0.6, 0.05),
                     limits = c(0, 0.6)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Standard error") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

plot_C <- ggplot(subset(plot_modality_n100_beta2, Metric == "Coverage"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0.00, 1.0, 0.05),
                     limits = c(0.00, 1)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Coverage") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

plot_D <- ggplot(subset(plot_modality_n100_beta2, Metric == "Proportion of failed"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 0.10, 0.02),
                     limits = c(0, 0.10)) +
  geom_hline(aes(yintercept = 0.0), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Proportion of failed") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

## 14.70 X 8 pdf save: 2X2 PLOT
pdf(paste0(temp_name, "_n100_b2.pdf"), width = 14.70, height = 8, family = "Times")
(plot_A + plot_B) / 
  (plot_C + plot_D) + 
  plot_layout(heights = c(1, 1)) & 
  theme(plot.margin = margin(t = 2, r = 3, b = -10, l = 3))
dev.off()


################################################################################
## Plot for beta3
plot_modality_n100_beta3 <- plot_modality_n100[plot_modality_n100$Parameter == "beta_3k", ]
plot_A <- ggplot(subset(plot_modality_n100_beta3, Metric == "Absolute bias"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 0.70, 0.05),
                     limits = c(0, 0.70)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Absolute bias") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

plot_B <- ggplot(subset(plot_modality_n100_beta3, Metric == "Standard error"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 0.65, 0.05),
                     limits = c(0, 0.65)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Standard error") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

plot_C <- ggplot(subset(plot_modality_n100_beta3, Metric == "Coverage"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0.00, 1.0, 0.05),
                     limits = c(0.00, 1)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Coverage") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

plot_D <- ggplot(subset(plot_modality_n100_beta3, Metric == "Proportion of failed"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 0.10, 0.02),
                     limits = c(0, 0.10)) +
  geom_hline(aes(yintercept = 0.0), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Proportion of failed") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

## 14.70 X 8 pdf save: 2X2 PLOT
pdf(paste0(temp_name, "_n100_b3.pdf"), width = 14.70, height = 8, family = "Times")
(plot_A + plot_B) / 
  (plot_C + plot_D) + 
  plot_layout(heights = c(1, 1)) & 
  theme(plot.margin = margin(t = 2, r = 3, b = -10, l = 3))
dev.off()


################################################################################
################################################################################
## Results for n = 500
all_results_n500 <- all_results[all_results$I == 500, ]
plot_modality_n500 <- data.frame(
  Method = rep(ifelse(stringr::str_split(all_results$xxx, "_", simplify = T)[,3] == "", 
                      stringr::str_split(all_results$xxx, "_", simplify = T)[,2],
                      paste0(stringr::str_split(all_results$xxx, "_", simplify = T)[,2], "*")), 4),
  Metric = c(rep("Absolute bias", nrow(all_results_n500)), 
             rep("Standard error", nrow(all_results_n500)),
             rep("Coverage", nrow(all_results_n500)),
             rep("Proportion of failed", nrow(all_results_n500))),
  Target = rep(all_results_n500$n_J, 4), 
  Parameter = rep(all_results_n500$beta, 4),
  Value = c(abs(as.numeric(all_results_n500$bias)), 
            as.numeric(all_results_n500$se_mean), 
            as.numeric(all_results_n500$coverage), 
            as.numeric((1000 - all_results_n500$b_analysis)/1000)))
#
plot_modality_n500$Method <- gsub("CLR", "CLR", plot_modality_n500$Method)
plot_modality_n500$Method <- gsub("GLMM", "GLMM", plot_modality_n500$Method)
plot_modality_n500$Method <- gsub("GLM", "GLM", plot_modality_n500$Method)
plot_modality_n500$Method <- gsub("GEE", "GEE", plot_modality_n500$Method)

#
plot_modality_n500$Target <- ifelse(plot_modality_n500$Target == 3, "3 Readers",
                                    ifelse(plot_modality_n500$Target == 5, "5 Readers", "10 Readers"))

#
plot_modality_n500$Target <- factor(plot_modality_n500$Target, 
                                    levels = c("3 Readers", "5 Readers", "10 Readers"))

#
plot_modality_n500$Method <- factor(plot_modality_n500$Method, 
                                    levels = c("CLR", 
                                               "GLM", "GLM*", 
                                               "GEE1", "GEE1*", 
                                               "GEE2", "GEE2*", 
                                               "glmer1", "glmer2", "glmmTMB1", "glmmTMB2"))


## Plot for beta2
plot_modality_n500_beta2 <- plot_modality_n500[plot_modality_n500$Parameter == "beta_2k", ]
plot_A <- ggplot(subset(plot_modality_n500_beta2, Metric == "Absolute bias"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 0.70, 0.05),
                     limits = c(0, 0.70)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Absolute bias") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

plot_B <- ggplot(subset(plot_modality_n500_beta2, Metric == "Standard error"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 0.27, 0.03),
                     limits = c(0, 0.27)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Standard error") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

plot_C <- ggplot(subset(plot_modality_n500_beta2, Metric == "Coverage"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0.00, 1.0, 0.05),
                     limits = c(0.00, 1)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Coverage") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

plot_D <- ggplot(subset(plot_modality_n500_beta2, Metric == "Proportion of failed"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 0.12, 0.02),
                     limits = c(0, 0.12)) +
  geom_hline(aes(yintercept = 0.0), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Proportion of failed") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

## 14.70 X 8 pdf save: 2X2 PLOT
pdf(paste0(temp_name, "_n500_b2.pdf"), width = 14.70, height = 8, family = "Times")
(plot_A + plot_B) / 
  (plot_C + plot_D) + 
  plot_layout(heights = c(1, 1)) & 
  theme(plot.margin = margin(t = 2, r = 3, b = -10, l = 3))
dev.off()


################################################################################
## Plot for beta3
plot_modality_n500_beta3 <- plot_modality_n500[plot_modality_n500$Parameter == "beta_3k", ]
plot_A <- ggplot(subset(plot_modality_n500_beta3, Metric == "Absolute bias"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 0.70, 0.05),
                     limits = c(0, 0.70)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Absolute bias") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

plot_B <- ggplot(subset(plot_modality_n500_beta3, Metric == "Standard error"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 0.30, 0.03),
                     limits = c(0, 0.30)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Standard error") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

plot_C <- ggplot(subset(plot_modality_n500_beta3, Metric == "Coverage"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0.0, 1.0, 0.05),
                     limits = c(0.0, 1)) +
  geom_hline(aes(yintercept = 0.95), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Coverage") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

plot_D <- ggplot(subset(plot_modality_n500_beta3, Metric == "Proportion of failed"), 
                 aes(x = "", y = Value)) +
  facet_grid(. ~ Method) + 
  geom_point(size = 1.6, position = position_dodge(width = 0.9), aes(color = Target)) + 
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("3 Readers", "5 Readers", "10 Readers")) +
  scale_y_continuous(expand = c(0.02, 0), breaks = seq(0, 0.12, 0.02),
                     limits = c(0, 0.12)) +
  geom_hline(aes(yintercept = 0.0), linetype = "dashed", color = "grey20", linewidth=0.5) +
  theme_bw() +
  ggtitle("Proportion of failed") + 
  theme(plot.title = element_text(face="bold", size=15, family="serif"),
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.spacing = unit(0, "lines"), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  theme(strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face="bold", size=8, family="serif")) + 
  theme(axis.text = element_text(face="bold", size=10, family="serif"),
        axis.ticks.x = element_blank()) + 
  theme(legend.key = element_blank(),
        legend.position = "none") +
  labs(x = NULL,
       y = NULL,
       shape = "Number of readers",
       color = "Number of readers")

## 14.70 X 8 pdf save: 2X2 PLOT
pdf(paste0(temp_name, "_n500_b3.pdf"), width = 14.70, height = 8, family = "Times")
(plot_A + plot_B) / 
  (plot_C + plot_D) + 
  plot_layout(heights = c(1, 1)) & 
  theme(plot.margin = margin(t = 2, r = 3, b = -10, l = 3))
dev.off()



