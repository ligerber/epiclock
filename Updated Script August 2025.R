#################################################################################
## This script uses linear mixed models to investigate the effect of social    ##
## bonds on epigenetic age and performs model selection                        ##
#################################################################################


# Load required packages and import data ---------------------------------------

library(ggpubr) #v0.6.0 - used for creating publication-ready plots and data visualisations
library(ggplot2) #v3.5.1 - used for creating various types of plots and data visualisations
library(glmmTMB) #v1.1.10 - used for fitting generalized linear mixed models with complex random effects structures
library(ggeffects) #v1.7.2 - used for visualising effects in statistical models, particularly useful for mixed models
library(gridExtra) #v2.3 - used for arranging multiple plots on a single page
library(car) #v3.1-3 - used for various statistical functions, including regression diagnostics
library(lme4) #v1.1-35.5 - used for fitting linear mixed-effects models
library(lmerTest) #v3.1-3 - used for obtaining p-values for linear mixed-effects models fitted with lme4
library(mediation) #v4.5.0 - used for conducting mediation analysis
library(parameters) #v0.25.0 - used for extracting and computing parameters from various statistical models, including confidence intervals
library(patchwork) #v1.3.0 - used for combining multiple ggplot2 plots into a single figure

# Import data
setwd("C:/Users/GER094/OneDrive - CSIRO/Legacy Publication")
fittedValues <- read.csv("Data/SocialAndBiolAgeData365_MappedTAdu_BowtieQ10_noPUC19.csv") 

# Data Exploration -------------------------------------------------------------

# Calculate AgeAccel
fittedValues$AgeAccelLOIOCV <- residuals(
  lm(fittedValues$transAge_tidymodels_LOIOCV~fittedValues$AgeDB)
)

# Prepare by removing rows with NA's for any of the variables of interest and all female data
dat <- droplevels(fittedValues[!is.na(fittedValues$transAge_tidymodels_LOIOCV) &
                                 fittedValues$Sex!='Female' &
                                 !is.na(fittedValues$NormStrength)&
                                 !is.na(fittedValues$YearSampled), ])

# Filter data for additional analyses where ages are known with higher accuracy
dat_HalfYear <- dat[(dat$Confidence<=183),]

# Filter data for additional analyses on adult males only
dat_Over14Years <- dat[(dat$AgeDB>=14),]

# Set YearSampled, DolphinID and Sex as factors (not characters)
dat$YearSampled <- factor(dat$YearSampled)
dat$DolphinID <- factor(dat$DolphinID)
dat$Sex <- factor(dat$Sex)

# Distribution outcome
ggplot(dat, aes(transAge_tidymodels_LOIOCV)) +
  geom_histogram(fill= "white", col= "black", binwidth= 2)

# Let's have a quick look
nrow(dat)              # The number of observations we can use
nlevels(dat$DolphinID) # This will be the sample size in each of the permutation models below
nlevels(dat$YearSampled) # This will be the maximum n years in any of the permutation models
length(which(dat$Sex=="Female")) # To verify that females were excluded
length(which(dat$Sex=="Male"))
table(dat$DolphinID)

sort(table(dat$DolphinID))

# Model selection --------------------------------------------------------------

# Step 1: Determine the optimal random effects structure
m_null <- lm(transAge_tidymodels_LOIOCV ~ 1, data = dat)
m_re1 <- lmer(transAge_tidymodels_LOIOCV ~ 1 + (1|DolphinID), data = dat, REML = TRUE)
m_re2 <- lmer(transAge_tidymodels_LOIOCV ~ 1 + (1|YearSampled), data = dat, REML = TRUE)
m_re3 <- lmer(transAge_tidymodels_LOIOCV ~ 1 + (1|DolphinID) + (1|YearSampled), data = dat, REML = TRUE)

# Compare random effects models using AIC and BIC
re_models <- list(m_null, m_re1, m_re2, m_re3)
re_model_names <- c("Null", "DolphinID", "YearSampled", "DolphinID + YearSampled")
re_aic_values <- sapply(re_models, AIC)
re_bic_values <- sapply(re_models, BIC)

re_results <- data.frame(
  Model = re_model_names,
  AIC = re_aic_values,
  BIC = re_bic_values,
  dAIC = re_aic_values - min(re_aic_values),
  dBIC = re_bic_values - min(re_bic_values)
)

# Print results for random effects models
cat("Random effects model comparison results:\n")
print(re_results)

# Select the best random effects structure based on AIC
best_re_model <- re_models[[which.min(re_aic_values)]]
cat("\nBest random effects model according to AIC:", re_model_names[which.min(re_aic_values)], "\n")

# Step 2: Fixed effects selection
m0 <- update(best_re_model, . ~ ., REML = FALSE)
m1 <- update(m0, transAge_tidymodels_LOIOCV ~ AgeDB + (1|DolphinID) + (1|YearSampled))
m2 <- update(m1, transAge_tidymodels_LOIOCV ~ AgeDB + NormStrength + (1|DolphinID) + (1|YearSampled))
m3 <- update(m2, transAge_tidymodels_LOIOCV ~ AgeDB + NormStrength + GroupSizeMale + (1|DolphinID) + (1|YearSampled))
m4 <- update(m3, transAge_tidymodels_LOIOCV ~ AgeDB + NormStrength + GroupSizeMale + CV + (1|DolphinID) + (1|YearSampled))
m5 <- update(m4, transAge_tidymodels_LOIOCV ~ AgeDB + NormStrength * GroupSizeMale + CV + (1|DolphinID) + (1|YearSampled))
m6 <- update(m0, transAge_tidymodels_LOIOCV ~ AgeDB + GroupSizeMale + (1|DolphinID) + (1|YearSampled))
m7 <- update(m0, transAge_tidymodels_LOIOCV ~ AgeDB + GroupSizeMale*NormStrength + (1|DolphinID) + (1|YearSampled))

# Compare fixed effects models
models <- list(m0, m1, m2, m3, m4, m5, m6, m7)
model_names <- c("Null", 
                 "AgeDB", 
                 "AgeDB + NormStrength", 
                 "AgeDB + NormStrength + GroupSizeMale", 
                 "AgeDB + NormStrength + GroupSizeMale + CV", 
                 "AgeDB + NormStrength * GroupSizeMale + CV",
                 "AgeDB + GroupSizeMale",
                 "AgeDB + GroupSizeMale * NormStrength")

# Calculate AIC, BIC, and log-likelihood
aic_values <- sapply(models, AIC)
bic_values <- sapply(models, BIC)
logLik_values <- sapply(models, logLik)

# Create a data frame with the results
results <- data.frame(
  Model = model_names,
  AIC = aic_values,
  BIC = bic_values,
  LogLik = logLik_values,
  dAIC = aic_values - min(aic_values),
  dBIC = bic_values - min(bic_values)
)

# Sort by AIC
results_aic <- results[order(results$AIC),]
results_bic <- results[order(results$BIC),]

# Print results
cat("Model comparison results (sorted by AIC):\n")
print(results_aic)
cat("\nModel comparison results (sorted by BIC):\n")
print(results_bic)

# Identify the best model
best_aic <- results_aic$Model[1]
best_bic <- results_bic$Model[1]

cat("\nBest model according to AIC:", best_aic)
cat("\nBest model according to BIC:", best_bic)

# Summary of the best model (using AIC as the criterion)
best_model <- models[[which(model_names == best_aic)]]
summary(best_model)

# UPDATED PLOTTING SECTION - Partial Residual Plots to Address Reviewer Concern ----

# Function to calculate partial residuals for mixed models
calculate_partial_residuals <- function(model, variable) {
  # Get model frame and response
  mf <- model.frame(model)
  response <- model.response(mf)
  
  # Get fixed effects coefficients
  fixed_coefs <- fixef(model)
  
  # Get model matrix for fixed effects
  X <- model.matrix(model, type = "fixed")
  
  # Find column index for the variable of interest
  var_col <- which(colnames(X) == variable)
  
  if (length(var_col) == 0) {
    stop(paste("Variable", variable, "not found in model"))
  }
  
  # Calculate fitted values without the variable of interest
  X_without_var <- X
  X_without_var[, var_col] <- 0
  fitted_without_var <- X_without_var %*% fixed_coefs
  
  # Calculate partial residuals
  partial_residuals <- response - fitted_without_var
  
  # Get the predictor values
  predictor_values <- X[, var_col]
  
  return(data.frame(
    x = predictor_values,
    y = as.numeric(partial_residuals)
  ))
}

# Function to create ggplot2 partial residual plots for mixed models
create_partial_plot_ggplot <- function(model, effect, label) {
  # Calculate partial residuals
  partial_data <- calculate_partial_residuals(model, effect)
  
  # Map original variable names to descriptive names
  effect_labels <- c(
    "AgeDB" = "Chronological Age",
    "NormStrength" = "Social Bond Strength", 
    "GroupSizeMale" = "Group Size"
  )
  
  effect_label <- effect_labels[effect]
  
  ggplot(partial_data, aes(x = x, y = y)) +
    geom_point(color = "black", alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "skyblue4", fill = "lightblue", alpha = 0.2) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.tag = element_text(size = 20, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      strip.text = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      panel.grid.major = element_line(size = 0.5),
      panel.grid.minor = element_line(size = 0.3)
    ) +
    labs(
      title = paste("Partial Effects for", effect_label),
      x = effect_label,
      y = "Partial Residuals (Epigenetic Age)",
      tag = label
    )
}

# Create plots for each fixed effect
fixed_effects <- c("AgeDB", "NormStrength", "GroupSizeMale")
plots <- list()
labels <- LETTERS[1:length(fixed_effects)]

for (i in seq_along(fixed_effects)) {
  plots[[i]] <- create_partial_plot_ggplot(best_model, fixed_effects[i], labels[i])
}

# Combine plots into a single panel
combined_plot <- plots[[1]] + plots[[2]] + plots[[3]] +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Partial Residual Plots for Fixed Effects",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))
  )


# Save the combined plot
ggsave("partial_residual_plots_combined.png", combined_plot, width = 18, height = 6, dpi = 300)


# Effect sizes / standardised effects -------------------------------------

# Create standardized versions of predictors
dat$AgeDB_std <- scale(dat$AgeDB)
dat$NormStrength_std <- scale(dat$NormStrength)
dat$GroupSizeMale_std <- scale(dat$GroupSizeMale)

# Fit model with standardized predictors only (outcome remains unstandardized)
std_model <- lmer(transAge_tidymodels_LOIOCV ~ AgeDB_std + NormStrength_std + GroupSizeMale_std +
                    (1 | DolphinID) + (1 | YearSampled), data = dat)


sd(dat$NormStrength, na.rm=T)
#0.1718419

sd(dat$GroupSizeMale)
# 1.4442

summary(std_model)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: transAge_tidymodels_LOIOCV ~ AgeDB_std + NormStrength_std + GroupSizeMale_std +      (1 | DolphinID) + (1 | YearSampled)
#    Data: dat
# 
# REML criterion at convergence: 234.1
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -1.66319 -0.38664 -0.01388  0.31193  1.78827 
# 
# Random effects:
#  Groups      Name        Variance Std.Dev.
#  DolphinID   (Intercept) 2.469    1.571   
#  YearSampled (Intercept) 6.400    2.530   
#  Residual                2.660    1.631   
# Number of obs: 50, groups:  DolphinID, 38; YearSampled, 13
# 
# Fixed effects:
#                   Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)        13.9990     0.8477  9.5580  16.514 2.40e-08 ***
# AgeDB_std           4.7589     0.6230 41.9683   7.638 1.81e-09 ***
# NormStrength_std   -1.6920     0.6631 38.0291  -2.552   0.0149 *  
# GroupSizeMale_std   1.3325     0.5268 42.5442   2.529   0.0152 *  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#             (Intr) AgDB_s NrmSt_
# AgeDB_std    0.155              
# NrmStrngth_ -0.044 -0.717       
# GrpSzMl_std  0.032  0.187 -0.505

model_parameters(std_model, ci = 0.95)
# # Fixed Effects 
# 
# Parameter         | Coefficient |   SE |         95% CI | t(43) |      p
# ------------------------------------------------------------------------
#   (Intercept)       |       14.00 | 0.85 | [12.29, 15.71] | 16.51 | < .001
# AgeDB std         |        4.76 | 0.62 | [ 3.50,  6.02] |  7.64 | < .001
# NormStrength std  |       -1.69 | 0.66 | [-3.03, -0.35] | -2.55 | 0.014 
# GroupSizeMale std |        1.33 | 0.53 | [ 0.27,  2.40] |  2.53 | 0.015 
# 
# # Random Effects 
# 
# Parameter                   | Coefficient |   SE |       95% CI
# ---------------------------------------------------------------
#   SD (Intercept: DolphinID)   |        1.57 | 0.60 | [0.74, 3.33]
# SD (Intercept: YearSampled) |        2.53 | 0.70 | [1.47, 4.36]
# SD (Residual)               |        1.63 | 0.49 | [0.90, 2.95]
# 
# Uncertainty intervals (equal-tailed) and p-values (two-tailed) computed using a Wald t-distribution
# approximation.

# Confidence intervals for non-standardised effects


# Get confidence intervals
model_parameters(best_model, ci = 0.95)
# 
# 
# # Fixed Effects 
# 
# Parameter     | Coefficient |   SE |          95% CI | t(43) |      p
# ---------------------------------------------------------------------
#   (Intercept)   |        1.75 | 1.72 | [ -1.73,  5.23] |  1.01 | 0.317 
# AgeDB         |        0.74 | 0.09 | [  0.55,  0.93] |  7.98 | < .001
# NormStrength  |       -9.75 | 3.72 | [-17.25, -2.26] | -2.63 | 0.012 
# GroupSizeMale |        0.90 | 0.35 | [  0.20,  1.61] |  2.59 | 0.013 
# 
# # Random Effects 
# 
# Parameter                   | Coefficient |   SE |       95% CI
# ---------------------------------------------------------------
#   SD (Intercept: DolphinID)   |        1.46 | 0.56 | [0.69, 3.12]
# SD (Intercept: YearSampled) |        2.32 | 0.62 | [1.38, 3.90]
# SD (Residual)               |        1.63 | 0.45 | [0.94, 2.81]
# 
# Uncertainty intervals (equal-tailed) and p-values (two-tailed) computed using a Wald t-distribution approximation.


# Mediation analysis with simple linear models ---------------------------------

# Total effect model
fit.totaleffect <- lm(transAge_tidymodels_LOIOCV ~ GroupSizeMale + AgeDB, data = dat)

# Mediator model
fit.mediator <- lm(NormStrength ~ GroupSizeMale + AgeDB, data = dat)

# Full model (DV)
fit.dv <- lm(transAge_tidymodels_LOIOCV ~ NormStrength + GroupSizeMale + AgeDB, data = dat)

# Mediation analysis
set.seed(123)
med.out <- mediate(fit.mediator, fit.dv, treat = 'GroupSizeMale', mediator = 'NormStrength', boot = TRUE, sims = 1000)

# Summary of mediation analysis
med_summary <- summary(med.out)
print(med_summary)

cat("\nInterpretation of Mediation Analysis:\n")
cat("Average Causal Mediation Effect (ACME):", med_summary$d.avg, "\n")
cat("Average Direct Effect (ADE):", med_summary$z.avg, "\n")
cat("Total Effect:", med_summary$tau.coef, "\n")
cat("Proportion Mediated:", med_summary$n.avg, "\n")

if (med_summary$d.avg.p < 0.05) {
  cat("The indirect effect is statistically significant, suggesting mediation.\n")
} else {
  cat("The indirect effect is not statistically significant, suggesting no strong evidence for mediation.\n")
}

# med_summary is 0.052. Thus, the mediation is not significant.


# Function to add correlation coefficient to plots
add_correlation <- function(x, y, ...) {
  cor_coef <- round(cor(x, y, use = "complete.obs"), 2)
  label <- paste("r =", cor_coef)
  ggplot2::annotate("text", x = Inf, y = Inf, label = label, hjust = 1, vjust = 1, ...)
}

# 1. Epigenetic Age vs Social Bond Strength
p1 <- ggplot(data = dat, aes(x = NormStrength, y = transAge_tidymodels_LOIOCV)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Social Bond Strength", y = "Epigenetic Age",
       title = "Epigenetic Age vs Social Bond Strength") +
  theme_minimal() +
  add_correlation(dat$NormStrength, dat$transAge_tidymodels_LOIOCV)

# 2. Epigenetic Age vs Group Size
p2 <- ggplot(data = dat, aes(x = GroupSizeMale, y = transAge_tidymodels_LOIOCV)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Group Size", y = "Epigenetic Age",
       title = "Epigenetic Age vs Group Size") +
  theme_minimal() +
  add_correlation(dat$GroupSizeMale, dat$transAge_tidymodels_LOIOCV)

# 3. Social Bond Strength vs Group Size
p3 <- ggplot(data = dat, aes(x = GroupSizeMale, y = NormStrength)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(x = "Group Size", y = "Social Bond Strength",
       title = "Social Bond Strength vs Group Size") +
  theme_minimal() +
  add_correlation(dat$GroupSizeMale, dat$NormStrength)


# Arrange all plots in a grid
combined_plot <- ggarrange(p1, p2, p3, ncol = 2, nrow = 2)

# Add overall title
final_plot <- annotate_figure(combined_plot, 
                              top = text_grob("Relationships between Epigenetic Age, Social Bond Strength, and Group Size", 
                                              face = "bold", size = 14))

# Display the plot
print(final_plot)

# Save the plot
ggsave("pairwise_relationships.png", final_plot, width = 12, height = 10, dpi = 300)

# Additional visualizations for pairwise relationships

p1 <- ggplot(dat, aes(x = GroupSizeMale, y = transAge_tidymodels_LOIOCV)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(x = "Group Size", y = "Epigenetic Age")

p2 <- ggplot(dat, aes(x = NormStrength, y = transAge_tidymodels_LOIOCV)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(x = "Social Bond Strength", y = "Epigenetic Age")

p3 <- ggplot(dat, aes(x = GroupSizeMale, y = NormStrength)) +
  geom_point() + geom_smooth(method = "lm") +
  labs(x = "Group Size", y = "Social Bond Strength")

grid.arrange(p1, p2, p3, ncol = 2)
ggsave("pairwise_relationships_table2fallacy.png", width = 12, height = 8, dpi = 300)

cat("\nThese analyses help address the Table 2 fallacy by showing how coefficients change across different model specifications and visualizing the pairwise relationships between variables.\n")

# Plot of epigenetic versus chronological ages ---------------------------------


# choose column containing the epigenetic age estimates to be displayed

#error LOIOCV
errorLOIOCV <- fittedValues$AgeDB - fittedValues$transAge_tidymodels_LOIOCV
MAE_LOIOCV <- sqrt(median(errorLOIOCV^2))
MAE_LOIOCV

# pearson correlation LOIOCV
r_LOIOCV <- cor(fittedValues$AgeDB, fittedValues$transAge_tidymodels_LOIOCV, method = "pearson")
r_LOIOCV


#plot LOIOCV or epigenetic ages generated with MammalMethyl package
p = ggplot(fittedValues, aes(x=AgeDB, y = transAge_tidymodels_LOIOCV, color=)) +
  ylim(0,40) + xlim(0,40) +
  ggtitle("LOIOCV") +
  geom_point(shape=1, size=2) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, colour="deepskyblue4") +
  geom_abline(intercept = 0, slope = 1, linetype = 3)+
  xlab("Chronological Age") +
  ylab("Epigenetic Age") +
  theme_bw()

theme(axis.title=element_blank(), axis.text = element_text(size=14, colour = "black"),);p

p + annotate("text", x = 0, y = 0, 
             label = paste("R =", round(r_LOIOCV, 2), "\nMAE =", round(MAE_LOIOCV, 2)), 
             hjust = 0, vjust = -10.5
)

ggsave("Plot_LOIOCV.tiff", units="cm", width=15, height=12, dpi=300)


#####
# Plot epigenetic clock for all individuals where chronological age is known with an accuracy of <=365
####

#error LOIOCV
errorLOIOCV <- fittedValues$transAge_tidymodels_LOIOCV - fittedValues$AgeDB
MAE_LOIOCV <- sqrt(median(errorLOIOCV^2))
MAE_LOIOCV

# pearson correlation LOIOCV
r_LOIOCV <- cor(fittedValues$AgeDB, fittedValues$transAge_tidymodels_LOIOCV, method = "pearson")
r_LOIOCV

fittedValues$inAnalyses <- factor(fittedValues$inAnalyses)
#plot LOIOCV
Plot = ggplot(fittedValues, aes(x=AgeDB, y = transAge_tidymodels_LOIOCV, col = inAnalyses)) +
  ylim(0,40) + xlim(0,40) +
  #geom_point(aes(colour = factor(inAnalyses)), size = 2) +
  geom_point(size=2, (aes(shape = inAnalyses))) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE, colour="deepskyblue4") +
  geom_abline(intercept = 0, slope = 1, linetype = 3)+
  xlab("Chronological Age") +
  ylab("Epigenetic Age") +
  theme_bw()
Plot + scale_color_manual(values = c("black", "black")) + scale_shape_manual(values = c(21,16)) + theme(legend.position = "none")

#ggsave("Plot_EpigeneticAgeAllSamples.tiff", units="cm", width=15, height=12, dpi=300)

# to run model the with estimates gained from various clocks, replace 'transAge_tidymodels_LOIOCV' with one of the following:
# BottlenoseBarratcloughGeneral, BottlenoseBarratcloghSkin, BottlenoseRobeckGeneral, BottlenoseRobeckSkin, OdontoceteGeneral, OdontoceteSkin, Universal3

best_mod <- lmer(transAge_tidymodels_LOIOCV ~ AgeDB + NormStrength + GroupSizeMale + (1|DolphinID) + (1|YearSampled), data = dat_Over14Years)
summary(best_mod)
