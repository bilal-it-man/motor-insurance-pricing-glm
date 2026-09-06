# 4_cross_validation.R
# 5-fold CV to check the Gini result holds across splits, plus a bootstrap
# CI on the original split.

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

  test_fold$pred_freq <- predict(freq_model, newdata = test_fold, type = "response",
                                 newoffset = log(rep(1, nrow(test_fold))))
  test_fold$pred_sev  <- predict(sev_model, newdata = test_fold, type = "response") + excess_loading
  test_fold$pure_prem <- test_fold$pred_freq * test_fold$pred_sev

  actual_rate      <- test_fold$ClaimAmount / test_fold$Exposure
  actual_freq_rate <- test_fold$ClaimNb / test_fold$Exposure

  list(
    model_gini     = gini_curve(actual_rate, test_fold$pure_prem, test_fold$Exposure),
    freq_only_gini = gini_curve(actual_freq_rate, test_fold$pred_freq, test_fold$Exposure)
  )
}

# ---- K-fold CV ----
set.seed(42)
k <- 5
folds <- sample(rep(1:k, length.out = nrow(dat)))

cv_results <- lapply(1:k, function(i) {
  cat("Fold", i, "of", k, "...\n")
  res <- fit_and_evaluate(dat[folds != i, ], dat[folds == i, ])
  data.frame(fold = i, model_gini = res$model_gini, freq_only_gini = res$freq_only_gini)
}) %>% bind_rows()

print(cv_results)
cat(sprintf("\nModel Gini across folds: mean %.4f, sd %.4f\n",
            mean(cv_results$model_gini), sd(cv_results$model_gini)))
cat(sprintf("Frequency-only Gini across folds: mean %.4f, sd %.4f\n",
            mean(cv_results$freq_only_gini), sd(cv_results$freq_only_gini)))

write.csv(cv_results, "output/cv_results.csv", row.names = FALSE)

# ---- Bootstrap CI on the original split ----
# Refit once, predict once, then resample the predictions (tests evaluation
# stability, not model-fitting stability - the CV above covers that).
train <- readRDS("data/train.rds")
test  <- readRDS("data/test.rds")

single_split <- fit_and_evaluate(train, test)
cat(sprintf("\nOriginal split - model Gini: %.4f, frequency-only Gini: %.4f\n",
            single_split$model_gini, single_split$freq_only_gini))

freq_model <- glm(freq_formula, family = quasipoisson(link = "log"),
                  offset = log(Exposure), data = train)
sev_train <- train %>%
  filter(ClaimNb > 0, ClaimAmount > 0) %>%
  mutate(AvgClaim = ClaimAmount / ClaimNb)
threshold <- quantile(sev_train$AvgClaim, 0.99)
sev_train <- sev_train %>% mutate(AvgClaim_capped = pmin(AvgClaim, threshold))
sev_model <- glm(sev_formula, family = Gamma(link = "log"), weights = ClaimNb, data = sev_train)
excess_loading <- mean(pmax(sev_train$AvgClaim - threshold, 0))

test$pred_freq <- predict(freq_model, newdata = test, type = "response",
                          newoffset = log(rep(1, nrow(test))))
test$pred_sev  <- predict(sev_model, newdata = test, type = "response") + excess_loading
test$pure_prem <- test$pred_freq * test$pred_sev

set.seed(42)
n_boot <- 2000
boot_ginis <- replicate(n_boot, {
  idx <- sample(nrow(test), replace = TRUE)
  gini_curve(test$ClaimAmount[idx] / test$Exposure[idx], test$pure_prem[idx], test$Exposure[idx])
})

ci <- quantile(boot_ginis, c(0.025, 0.975))
cat(sprintf("Bootstrap 95%% CI for model Gini (%d resamples): [%.4f, %.4f]\n",
            n_boot, ci[1], ci[2]))

write.csv(
  data.frame(metric = c("cv_model_gini_mean", "cv_model_gini_sd",
                        "cv_freq_gini_mean", "cv_freq_gini_sd",
                        "boot_ci_lower", "boot_ci_upper"),
             value  = c(mean(cv_results$model_gini), sd(cv_results$model_gini),
                        mean(cv_results$freq_only_gini), sd(cv_results$freq_only_gini),
                        ci[1], ci[2])),
  "output/cv_summary.csv", row.names = FALSE
)
