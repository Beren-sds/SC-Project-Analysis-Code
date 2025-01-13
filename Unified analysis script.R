# Unified scripts

# --------------------------------------------------------------
# Part 1: Data Loading & Preparation, EDA
# --------------------------------------------------------------

# Relevant libraries
library(dplyr)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(corrplot)
library(car)
library(Matrix)
library(lme4)
library(lmerTest)
library(MASS)
library(effects)
library(emmeans)
library(lattice)
library(nlme)
library(mice)
library(multcomp)
library(kableExtra)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(patchwork)
library(ggthemes)
library(MuMIn)
library(robustlmm)
library(performance)
library(see)
library(ggeffects)
library(effectsize)
library(mgcv)
library(viridis)
library(GGally)
library(broom.mixed)
library(DHARMa)
library(sjPlot)
library(forcats)        # Factor manipulations
library(hrbrthemes)     # Plotting themes
library(ggdist)         # Uncertainty visualisation
library(latex2exp)      # LaTeX labels
library(brms)           # Bayesian fitting
library(bayesplot)      # Plot MCMC output
library(tidybayes)      # Manipulate & visualize posterior draws
library(gtsummary)      # Nice summary tables

# --------------------------------------------------------------
# 1. Data Loading and Preparation
# --------------------------------------------------------------
data <- read.csv("jagger_data.csv")

# Factors & transformations (unified):
data$diel <- factor(data$diel)
data$noct <- factor(data$noct)
data$rep  <- factor(data$rep)

# Multiply dm by 100 for more convenient numeric scale
data$dm_scaled <- data$dm * 100

# Create composite treatment factor (used in frequentist script)
data$treatment <- paste(data$diel, "h day,", data$noct, "h dark")
data$treatment <- factor(data$treatment,
                         levels = c("20 h day, 2 h dark",
                                    "22 h day, 2 h dark",
                                    "24 h day, 2 h dark",
                                    "26 h day, 2 h dark",
                                    "28 h day, 2 h dark",
                                    "20 h day, 4 h dark",
                                    "22 h day, 4 h dark",
                                    "24 h day, 4 h dark",
                                    "26 h day, 4 h dark",
                                    "28 h day, 4 h dark",
                                    "22 h day, 6 h dark",
                                    "24 h day, 6 h dark",
                                    "26 h day, 6 h dark",
                                    "28 h day, 6 h dark"))

# --------------------------------------------------------------
# 2. EDA
# --------------------------------------------------------------
# 2.1 Raw Distribution
p1 <- ggplot(data, aes(x = diel, y = dm_scaled, fill = noct)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
  theme_minimal(base_size = 12) +
  scale_fill_viridis_d(option = "D") +
  labs(x = "Day Length (hours)",
       y = "Dry Matter Content (%)",
       fill = "Dark Period (hours)",
       title = "Distribution of Measurements",
       subtitle = "Box plots show quartiles, points show individual observations")

p2 <- ggplot(data, aes(x = diel, y = dm_scaled, fill = noct)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.3, outlier.shape = NA) +
  theme_minimal(base_size = 12) +
  scale_fill_viridis_d(option = "D") +
  labs(x = "Day Length (hours)",
       y = "Dry Matter Content (%)",
       fill = "Dark Period (hours)",
       title = "Distribution Patterns",
       subtitle = "Violin plots show full data distribution")

# --------------------------------------------------------------
# 3. Treatment Effect & Means
# --------------------------------------------------------------
# 3.1 Group means for (diel, noct) combos (used in plots)
interaction_means <- data %>%
  group_by(diel, noct) %>%
  summarise(
    mean_dm = mean(dm_scaled),
    sd = sd(dm_scaled),
    n = n(),
    se = sd/sqrt(n),
    .groups = "drop"
  )

# 3.2 Interaction Plots
p3 <- ggplot(interaction_means,
             aes(x = diel, y = mean_dm, color = noct, group = noct)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top") +
  scale_color_viridis_d(option = "D") +
  labs(x = "Day Length (hours)",
       y = "Dry Matter Content (%)",
       color = "Dark Period (hours)",
       title = "Interaction Pattern",
       subtitle = "Average effects across treatment combinations")

p4 <- ggplot(interaction_means, aes(x = diel, y = mean_dm)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_dm - se, ymax = mean_dm + se), width = 0.2) +
  facet_wrap(~noct, labeller = label_both) +
  theme_minimal(base_size = 12) +
  labs(x = "Day Length (hours)",
       y = "Dry Matter Content (%)",
       title = "Dark Period-Specific Patterns",
       subtitle = "Separate trends for each dark period length")


# --------------------------------------------------------------
# Part 2: Frequentist Analysis (Linear Mixed Models, Contrasts)
# --------------------------------------------------------------

# 4. Mixed Model
model_oneway <- lmer(dm_scaled ~ diel * noct + (1|rep), data = data, REML = TRUE)

# 5. EMMs (Estimated Marginal Means)
emms_main <- emmeans(model_oneway, ~ diel * noct)
emms_df   <- as.data.frame(emms_main)

# 5.1 Partial averaging for diel effect
diel_subset <- emms_df %>%
  filter(noct %in% c("2","4")) %>%
  group_by(diel) %>%
  summarise(emmean = mean(emmean, na.rm = TRUE))

p5_partial <- ggplot(diel_subset, aes(x = as.numeric(as.character(diel)), y = emmean)) +
  geom_line(linewidth = 1, group = 1) +
  geom_point(size = 3) +
  theme_minimal(base_size = 12) +
  labs(x = "Day Length (hours)",
       y = "Dry Matter Content (%)",
       title = "Day Length Main Effect (Partial Averaging)",
       subtitle = "Averaged over noct=2 and 4 only")

# 5.2 Partial averaging for noct effect
noct_subset <- emms_df %>%
  filter(diel %in% c("22","24","26","28")) %>% 
  group_by(noct) %>%
  summarise(emmean = mean(emmean, na.rm = TRUE))

p6_partial <- ggplot(noct_subset, aes(x = as.numeric(as.character(noct)), y = emmean)) +
  geom_line(linewidth = 1, group = 1) +
  geom_point(size = 3) +
  theme_minimal(base_size = 12) +
  labs(x = "Dark Period (hours)",
       y = "Dry Matter Content (%)",
       title = "Dark Period Main Effect (Partial Averaging)",
       subtitle = "Averaged over diel=22,24,26,28 only")

# --------------------------------------------------------------
# 6. Custom Contrast Testing
# --------------------------------------------------------------
# Use the single-factor model (treatment) to define and test custom contrasts
model_oneway <- lmer(dm_scaled ~ treatment + (1|rep), data = data, REML = TRUE)
emms <- emmeans(model_oneway, ~ treatment)

# Contrast matrix (7 contrasts for 14-level factor)
my_contrasts <- list(
  "Dial20 vs Dial22: 2h vs 4h" = c( 1, -1,  0,  0,  0,
                                   -1,  1,  0,  0,  0,
                                    0,  0,  0,  0),
  "Dial20 vs Dial24: 2h vs 4h" = c( 1,  0, -1,  0,  0,
                                   -1,  0,  1,  0,  0,
                                    0,  0,  0,  0),
  "Dial20 vs Dial26: 2h vs 4h" = c( 1,  0,  0, -1,  0,
                                   -1,  0,  0,  1,  0,
                                    0,  0,  0,  0),
  "Dial20 vs Dial28: 2h vs 4h" = c( 1,  0,  0,  0, -1,
                                   -1,  0,  0,  0,  1,
                                    0,  0,  0,  0),
  "Dial22 vs Dial24: 2h vs 6h" = c( 0,  1, -1,  0,  0,
                                    0,  0,  0,  0,  0,
                                   -1,  1,  0,  0),
  "Dial22 vs Dial26: 2h vs 6h" = c( 0,  1,  0, -1,  0,
                                    0,  0,  0,  0,  0,
                                   -1,  0,  1,  0),
  "Dial22 vs Dial28: 2h vs 6h" = c( 0,  1,  0,  0, -1,
                                    0,  0,  0,  0,  0,
                                   -1,  0,  0,  1)
)

contrast_results <- contrast(emms, method = my_contrasts, adjust = "bonferroni")
contrast_df <- as.data.frame(contrast_results)

Day <- c("20 vs 22", "20 vs 24", "20 vs 26", "20 vs 28",
         "22 vs 24", "22 vs 26", "22 vs 28")
Type <- c("2h vs 4h", "2h vs 4h", "2h vs 4h", "2h vs 4h",
          "2h vs 6h", "2h vs 6h", "2h vs 6h")

contrast_data <- data.frame(
  Contrast = rownames(contrast_df),
  Day = factor(Day),
  Type = factor(Type),
  Estimate = contrast_df$estimate,
  SE = contrast_df$SE,
  p_value = contrast_df$p.value
)

contrast_data$Significance <- ifelse(contrast_data$p_value < 0.05, "*", "ns")

# 7. Visualize Contrasts
p7 <- ggplot(contrast_data, 
             aes(y = reorder(paste(Day, Type), Estimate),
                 x = Estimate, color = Type)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = Estimate - 1.96*SE,
                     xmax = Estimate + 1.96*SE),
                 height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  scale_color_brewer(palette = "Set2") +
  labs(y = "Contrast",
       x = "Estimated Difference",
       title = "Forest Plot of Contrast Estimates",
       subtitle = "Red line indicates no effect")

p8 <- ggplot(contrast_data, aes(x = Day, y = Estimate, color = Type)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = Estimate - 1.96*SE, ymax = Estimate + 1.96*SE),
                width = 0.2,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  scale_color_brewer(palette = "Set2") +
  labs(x = "Day Length (hours)",
       y = "Contrast Estimate",
       color = "Comparison Type",
       title = "Effect Sizes by Day Length",
       subtitle = "Comparing dark period differences")

# --------------------------------------------------------------
# 8. Diagnostics
# --------------------------------------------------------------
resids <- residuals(model_oneway)
fitted_vals <- fitted(model_oneway)

p9 <- ggplot(data.frame(resids = resids), aes(sample = resids)) +
  stat_qq() +
  stat_qq_line() +
  theme_minimal(base_size = 12) +
  labs(title = "Normal Q-Q Plot",
       subtitle = "Check for normality of residuals")

p10 <- ggplot(data.frame(fitted = fitted_vals, resids = resids),
              aes(x = fitted, y = resids)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "loess", se = FALSE, color = "blue", alpha = 0.5) +
  theme_minimal(base_size = 12) +
  labs(x = "Fitted Values",
       y = "Residuals",
       title = "Residuals vs Fitted",
       subtitle = "Check for homoscedasticity")

# --------------------------------------------------------------
# 9. Results and Summaries
# --------------------------------------------------------------
cat("One-way ANOVA Results\n")
anova_table <- anova(model_oneway, type = "III")

print(kable(anova_table,
            caption = "Type III Analysis of Variance",
            digits = 3) %>%
        kable_styling(bootstrap_options = c("striped", "hover")))

means_table <- interaction_means %>%
  mutate(across(where(is.numeric), round, 3)) %>%
  rename("Day Length" = diel,
         "Dark Period" = noct,
         "Mean" = mean_dm,
         "SD" = sd,
         "N" = n,
         "SE" = se)

print(kable(means_table,
            caption = "Descriptive Statistics by Treatment Combination",
            digits = 3) %>%
        kable_styling(bootstrap_options = c("striped", "hover")))

contrast_table <- data.frame(
  Contrast = paste(contrast_data$Day, contrast_data$Type),
  Estimate = contrast_data$Estimate,
  SE = contrast_data$SE,
  "Z value" = contrast_df$t.ratio,
  "P value" = contrast_data$p_value,
  Significance = ifelse(contrast_data$p_value < 0.05, "*", "ns")
) %>%
  mutate(across(where(is.numeric), round, 3))

print(kable(contrast_table,
            caption = "Custom Contrast Analysis Results",
            digits = 3) %>%
        kable_styling(bootstrap_options = c("striped", "hover")))

# Plot objects from frequentist analysis:
p1; p2; p3; p4; p5_partial; p6_partial; p7; p8; p9; p10

# --------------------------------------------------------------
# 10. Variance Components & Random Effects
# --------------------------------------------------------------
vc <- VarCorr(model_oneway)
var_components <- data.frame(
  Component = c("Random Effects (rep)", "Residual"),
  Variance = c(attr(vc$rep, "stddev")^2, attr(vc, "sc")^2),
  Std_Dev = c(attr(vc$rep, "stddev"), attr(vc, "sc")),
  Pct_Total = NA
)
var_components$Pct_Total <- var_components$Variance / sum(var_components$Variance) * 100

variance_table <- kable(var_components,
                        col.names = c("Component", "Variance", "Std. Deviation", "% of Total"),
                        caption = "Variance Components Analysis",
                        digits = 3) %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
variance_table

re <- ranef(model_oneway, condVar = TRUE)
re_df <- data.frame(
  Level = factor(rownames(re$rep)),
  Effect = re$rep[,1],
  SE = sqrt(attr(re$rep, "postVar")[1,,])
)

if(nrow(re_df) > 0) {
  random_effects_plot <- ggplot(re_df, aes(x = reorder(Level, Effect), y = Effect)) +
    geom_point() +
    geom_errorbar(aes(ymin = Effect - 2*SE, ymax = Effect + 2*SE), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Random Effects by Replication",
         subtitle = "Estimates ± 2 Standard Errors",
         x = "Replication Level",
         y = "Random Effect Estimate")
  random_effects_plot
} else {
  random_effects_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No random effects data") +
    theme_void()
  warning("Random effects extraction resulted in empty data frame")
}

treatment_effects <- emmeans(model_oneway, ~ treatment)
pairs_results <- pairs(treatment_effects)
effect_sizes <- as.data.frame(pairs_results)

# Compute t.ratio, p.value, CIs
alpha <- 0.05
effect_sizes <- effect_sizes %>%
  mutate(
    t.ratio  = estimate / SE,
    p.value  = 2 * (1 - pt(abs(t.ratio), df)),
    lower    = estimate - qt(1 - alpha/2, df) * SE,
    upper    = estimate + qt(1 - alpha/2, df) * SE
  ) %>%
  arrange(desc(abs(estimate)))

effect_size_plot <- ggplot(effect_sizes %>% filter(abs(estimate) > mean(abs(estimate))), 
                           aes(x = estimate, y = reorder(contrast, estimate))) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red', alpha = 0.5) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8),
        plot.title  = element_text(size = 11)) +
  labs(title = "Treatment Effect Sizes with 95% Confidence Intervals",
       subtitle = "Effects larger than mean absolute effect size",
       x = "Estimated Difference in Dry Matter Content (%)",
       y = "Treatment Contrast")
effect_size_plot

effect_summary <- effect_sizes %>%
  dplyr::select(contrast, estimate, SE, df, t.ratio, p.value, lower, upper) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

effect_size_table <- kable(effect_summary,
                           caption = "Treatment Contrast Analysis Results",
                           col.names = c("Contrast", "Estimate", "SE", "df", 
                                         "t-ratio", "p-value", "Lower CI", "Upper CI")) %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
effect_size_table

dharma_residuals <- simulateResiduals(model_oneway)
plot(dharma_residuals)

# QQ plot with confidence intervals
compute_qq_bands <- function(residuals, conf = 0.95) {
  n <- length(residuals)
  probs <- ppoints(n)
  theoretical <- qnorm(probs)
  
  ordered <- sort(residuals)
  se <- sqrt(probs * (1 - probs) / n) / dnorm(theoretical)
  
  alpha <- 1 - conf
  z <- qnorm(1 - alpha/2)
  upper <- ordered + z * se * sd(residuals)
  lower <- ordered - z * se * sd(residuals)
  
  data.frame(
    theoretical = theoretical,
    observed    = ordered,
    upper       = upper,
    lower       = lower
  )
}

qq_data <- compute_qq_bands(resid(model_oneway))
qq_plot_advanced <- ggplot(qq_data, aes(x = theoretical, y = observed)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  theme_minimal() +
  labs(title = "Normal Q-Q Plot with 95% Confidence Bands",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")
qq_plot_advanced

performance_metrics <- model_performance(model_oneway)
performance_table <- kable(as.data.frame(performance_metrics),
                           caption = "Model Performance Metrics",
                           digits = 3) %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
performance_table

r2_decomposition <- r2(model_oneway)
r2_table <- kable(as.data.frame(r2_decomposition),
                  caption = "R-squared Decomposition",
                  digits = 3) %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
r2_table


# --------------------------------------------------------------
# Part 3: Bayesian Approach
# --------------------------------------------------------------

# Model formula with random intercept
# Using dm_scaled rather than dm_100 for consistency
bayes_formula <- bf(dm_scaled ~ 1 + diel + noct + diel:noct + (1|rep))

# possible priors
get_prior(formula = bayes_formula, data = data)

# Example prior specification (adjust as needed)
#   Intercept ~ Normal(7.5, 0.5)
#   Coefficients ~ Normal(0, 1)
#   Random Effect SD ~ Student_t(3, 0, 1)
#   Residual Error (sigma) ~ Exponential(1/0.67)
prior_data <- c(
  set_prior("normal(0, 1)", class = "b"),
  set_prior("normal(7.5, 0.5)", class = "Intercept"),
  set_prior("student_t(3, 0, 1)", class = "sd", group = "rep", lb = 0),
  set_prior("exponential(1 / 0.67)", class = "sigma", lb = 0)
)

# Bayesian model
bayesian_mixed_effects_anova <- brm(
  formula   = bayes_formula,
  data      = data,
  seed      = 18,
  prior     = prior_data,
  save_pars = save_pars(all = TRUE),
  control   = list(adapt_delta = 0.98),
  cores     = 4
)

# Summaries
summary(bayesian_mixed_effects_anova)
tbl_regression(bayesian_mixed_effects_anova, intercept = TRUE)

# Bayesian R-squared
bayes_R2(bayesian_mixed_effects_anova, re.form = NULL, re_formula = NULL, summary = TRUE)
post_r2 <- bayes_R2(bayesian_mixed_effects_anova, summary = FALSE)
ggplot(data = data.frame(r2 = post_r2[, 1])) +
  stat_histinterval(aes(x = r2)) +
  xlab(TeX("$R^2$")) +
  ylab("density")

# Posterior predictive checks
pp_check(bayesian_mixed_effects_anova, type = "scatter_avg")
pp_check(bayesian_mixed_effects_anova, type = "scatter_avg_grouped", group = "rep")

# MCMC Diagnostics
mcmc_trace(bayesian_mixed_effects_anova,
           pars = c("b_Intercept", "b_diel22", "b_diel24", "b_diel26",
                    "b_diel28", "b_noct4", "b_noct6", 
                    "b_diel22:noct4", "b_diel24:noct4", 
                    "b_diel26:noct4", "b_diel28:noct4", 
                    "b_diel22:noct6", "b_diel24:noct6",
                    "b_diel26:noct6", "b_diel28:noct6"))

mcmc_trace(bayesian_mixed_effects_anova,
           pars = c("sd_rep__Intercept", "sigma", "Intercept", 
                    "r_rep[1,Intercept]", "r_rep[2,Intercept]", 
                    "r_rep[3,Intercept]", "r_rep[4,Intercept]"))

mcmc_acf(bayesian_mixed_effects_anova,
         pars = c("b_Intercept", "b_diel22", "b_diel24", "b_diel26",
                  "b_diel28", "b_noct4", "b_noct6",  "b_diel22:noct4",
                  "b_diel24:noct4", "b_diel26:noct4", "b_diel28:noct4",
                  "b_diel22:noct6", "b_diel24:noct6","b_diel26:noct6",
                  "b_diel28:noct6"))

# Residual diagnostics
diagnostics_df <- data
diagnostics_df$fitted    <- fitted(bayesian_mixed_effects_anova)[, 1]
diagnostics_df$residuals <- residuals(bayesian_mixed_effects_anova)[, 1]

ggplot(diagnostics_df) +
  geom_point(aes(x = fitted, y = residuals)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red")

ggplot(diagnostics_df) +
  geom_point(aes(x = residuals, y = diel, colour = diel)) +
  theme(legend.position = "bottom")  +
  scale_colour_viridis_d(option = "turbo")

ggplot(diagnostics_df) +
  geom_point(aes(x = residuals, y = noct, colour = noct)) +
  theme(legend.position = "bottom")  +
  scale_colour_viridis_d(option = "turbo")

ggplot(diagnostics_df, aes(sample = residuals)) +
  geom_qq() +
  geom_qq_line(colour = "red")

# Posterior predictions for a missing combo (e.g., diel=20, noct=6)
newdata <- data.frame(
  diel = c(20, 20, 20, 20),
  noct = c(6, 6, 6, 6),
  rep  = c(1, 2, 3, 4)
)
pp <- posterior_predict(bayesian_mixed_effects_anova, newdata = newdata)
mean(pp); sd(pp)

# Marginal Effects
marginal <- conditional_effects(bayesian_mixed_effects_anova, effects = "diel:noct")
marginal_plot <- plot(marginal, plot = FALSE)[[1]] + 
  labs(title = "Marginal effect estimates",
       x = "Day length (h)",
       y = "Dry matter content (%)",
       color = 'Dark period (h)') +
  theme_minimal()
marginal_plot

# Random Effects
ramdom_effects <- as.data.frame(ranef(bayesian_mixed_effects_anova)$rep)
random_effects <- cbind(rep = rownames(ramdom_effects), ramdom_effects)

ggplot(random_effects, aes(y = rep, x = Estimate.Intercept)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = Q2.5.Intercept, xmax = Q97.5.Intercept), width = 0.2) +
  labs(title = "Random effect estimates",
       x = "Random effect",
       y = "Replicate") +
  theme_minimal() +
  geom_vline(xintercept = 0, linetype = "dashed")

# Pairwise Comparisons for diel
bayesian_emm_diel <- emmeans(bayesian_mixed_effects_anova, ~ diel, re_formula = NA)
pairs(bayesian_emm_diel)

# Pairwise Comparisons for noct
bayesian_emm_noct <- emmeans(bayesian_mixed_effects_anova, ~ noct, re_formula = NA)
pairs(bayesian_emm_noct)

