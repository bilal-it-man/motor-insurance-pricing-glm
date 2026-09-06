# 4_cross_validation.R
# Three checks:
#  1. 5-fold CV, per-fold Gini (as before)
#  2. Out-of-fold stacking: every row held out exactly once across the 5
#     folds, pooled into one ~678K-row evaluation instead of five ~135K ones
#  3. Capped-actual evaluation: tests the model against the same 99th-
#     percentile-capped target it was actually built to predict, not the
#     raw uncapped cost that includes claims nothing can predict

library(dplyr)
source("R/functions.R")

dat <- bind_rows(readRDS("data/train.rds"), readRDS("data/test.rds"))

freq_formula <- ClaimNb ~ DrivAgeBand + VehAgeBand + VehPower + VehGas + Area + BonusMalus
sev_formula  <- AvgClaim_capped ~ VehPower + Area

fit_and_evaluate <- function(train_fold, test_fold) {
  freq_model <- glm(freq_formula, family = poisson(link = "log"),
                    offset = log(Exposure), data = train_fold)
  disp <- sum(residuals(freq_model, type = "pearson")^2) / freq_model$df.residual
  if (disp > 1.2) {
    freq_model <- glm(freq_formula, family = quasipoisson(link = "log"),
                      offset = log(Exposure), data = train_fold)
  }
  
  sev_train <- train_fold %>%
    filter(ClaimNb > 0, ClaimAmount > 0) %>%
    mutate(AvgClaim = ClaimAmount / ClaimNb)
  threshold <- quantile(sev_train$AvgClaim, 0.99)
  sev_train <- sev_train %>% mutate(AvgClaim_capped = pmin(AvgClaim, threshold))
  
  sev_model <- glm(sev_formula, family = Gamma(link = "log"),
                   weights = ClaimNb, data = sev_train)
  excess_loading <- mean(pmax(sev_train$AvgClaim - threshold, 0))
  
  test_fold$pred_freq       <- predict(freq_model, newdata = test_fold, type = "response",
                                       newoffset = log(rep(1, nrow(test_fold))))
  test_fold$pred_sev_capped <- predict(sev_model, newdata = test_fold, type = "response")
  test_fold$pure_prem       <- test_fold$pred_freq * (test_fold$pred_sev_capped + excess_loading)
  test_fold$attritional_prem <- test_fold$pred_freq * test_fold$pred_sev_capped
  
  actual_rate      <- test_fold$ClaimAmount / test_fold$Exposure
  actual_freq_rate <- test_fold$ClaimNb / test_fold$Exposure
  
  list(
    model_gini      = gini_curve(actual_rate, test_fold$pure_prem, test_fold$Exposure),
    freq_only_gini  = gini_curve(actual_freq_rate, test_fold$pred_freq, test_fold$Exposure),
    predictions     = test_fold %>%
      select(ClaimNb, ClaimAmount, Exposure, pred_freq, pure_prem, attritional_prem)
  )
}

# ---- K-fold CV ----
set.seed(42)
k <- 5
folds <- sample(rep(1:k, length.out = nrow(dat)))

fold_ginis <- list()
oof_predictions <- list()

for (i in 1:k) {
  cat("Fold", i, "of", k, "...\n")
  res <- fit_and_evaluate(dat[folds != i, ], dat[folds == i, ])
  fold_ginis[[i]]      <- data.frame(fold = i, model_gini = res$model_gini,
                                     freq_only_gini = res$freq_only_gini)
  oof_predictions[[i]] <- res$predictions
}

cv_results <- bind_rows(fold_ginis)
print(cv_results)
cat(sprintf("\nModel Gini across folds: mean %.4f, sd %.4f\n",
            mean(cv_results$model_gini), sd(cv_results$model_gini)))
cat(sprintf("Frequency-only Gini across folds: mean %.4f, sd %.4f\n",
            mean(cv_results$freq_only_gini), sd(cv_results$freq_only_gini)))

write.csv(cv_results, "output/cv_results.csv", row.names = FALSE)

# ---- Out-of-fold stacked evaluation ----
# Every row appears exactly once, held out - one ~678K-row evaluation
# instead of five ~135K-row ones.
oof <- bind_rows(oof_predictions)

oof_actual_rate      <- oof$ClaimAmount / oof$Exposure
oof_actual_freq_rate  <- oof$ClaimNb / oof$Exposure

oof_model_gini <- gini_curve(oof_actual_rate, oof$pure_prem, oof$Exposure)
oof_freq_gini  <- gini_curve(oof_actual_freq_rate, oof$pred_freq, oof$Exposure)

cat(sprintf("\nOut-of-fold pooled (n=%d) - model Gini: %.4f, frequency-only Gini: %.4f\n",
            nrow(oof), oof_model_gini, oof_freq_gini))

# ---- Capped-actual evaluation ----
# Tests the model against the target it was actually built to predict:
# attritional cost (capped at the 99th percentile), evaluated with
# attritional_prem, which excludes the flat excess loading on both sides.
claims_all <- dat %>% filter(ClaimNb > 0, ClaimAmount > 0) %>%
  mutate(AvgClaim = ClaimAmount / ClaimNb)
global_cap <- quantile(claims_all$AvgClaim, 0.99)
cat(sprintf("Global cap for capped-actual evaluation: %.2f\n", global_cap))

oof <- oof %>%
  mutate(ClaimAmount_capped = ifelse(ClaimNb > 0,
                                     pmin(ClaimAmount / ClaimNb, global_cap) * ClaimNb, 0))
oof_capped_actual_rate <- oof$ClaimAmount_capped / oof$Exposure

oof_attritional_gini <- gini_curve(oof_capped_actual_rate, oof$attritional_prem, oof$Exposure)

set.seed(42)
n_reps <- 200
attritional_baseline_ginis <- replicate(n_reps, {
  shuffled <- sample(nrow(oof))
  gini_curve(oof_capped_actual_rate[shuffled], rep(1, nrow(oof)), oof$Exposure[shuffled])
})
attritional_baseline_mean <- mean(attritional_baseline_ginis)
attritional_baseline_sd    <- sd(attritional_baseline_ginis)

cat(sprintf("Attritional (capped-actual) Gini: %.4f\n", oof_attritional_gini))
cat(sprintf("Attritional flat baseline across %d tie-breaks: mean %.4f, sd %.4f\n",
            n_reps, attritional_baseline_mean, attritional_baseline_sd))

# ---- Bootstrap CI on the out-of-fold attritional result ----
set.seed(42)
n_boot <- 2000
boot_attritional_ginis <- replicate(n_boot, {
  idx <- sample(nrow(oof), replace = TRUE)
  gini_curve(oof_capped_actual_rate[idx], oof$attritional_prem[idx], oof$Exposure[idx])
})
ci_attritional <- quantile(boot_attritional_ginis, c(0.025, 0.975))
cat(sprintf("Bootstrap 95%% CI for attritional Gini (%d resamples, n=%d): [%.4f, %.4f]\n",
            n_boot, nrow(oof), ci_attritional[1], ci_attritional[2]))

write.csv(
  data.frame(metric = c("cv_model_gini_mean", "cv_model_gini_sd",
                        "cv_freq_gini_mean", "cv_freq_gini_sd",
                        "oof_model_gini", "oof_freq_gini",
                        "oof_attritional_gini", "attritional_baseline_mean",
                        "attritional_baseline_sd", "attritional_ci_lower",
                        "attritional_ci_upper", "global_cap"),
             value  = c(mean(cv_results$model_gini), sd(cv_results$model_gini),
                        mean(cv_results$freq_only_gini), sd(cv_results$freq_only_gini),
                        oof_model_gini, oof_freq_gini,
                        oof_attritional_gini, attritional_baseline_mean,
                        attritional_baseline_sd, ci_attritional[1],
                        ci_attritional[2], global_cap)),
  "output/cv_summary.csv", row.names = FALSE
)