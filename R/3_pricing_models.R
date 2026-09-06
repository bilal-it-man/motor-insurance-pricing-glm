# 3_pricing_models.R
# Frequency (Poisson) + severity (Gamma) GLMs -> pure premium, evaluated
# against a flat baseline.

library(dplyr)

source("R/functions.R")

dir.create("output", showWarnings = FALSE)

train <- readRDS("data/train.rds")
test  <- readRDS("data/test.rds")

# ---- Frequency ----
# Region and density were tested and dropped: negligible added power once
# area and bonus-malus are in the model.
freq_model <- glm(
  ClaimNb ~ DrivAgeBand + VehAgeBand + VehPower + VehGas + Area + BonusMalus,
  family = poisson(link = "log"),
  offset = log(Exposure),
  data   = train
)
summary(freq_model)

disp <- sum(residuals(freq_model, type = "pearson")^2) / freq_model$df.residual
cat("Frequency dispersion:", round(disp, 3), "\n")

if (disp > 1.2) {
  cat("Dispersion above 1.2: refitting as quasi-Poisson for correct standard errors.\n")
  freq_model_final <- glm(
    ClaimNb ~ DrivAgeBand + VehAgeBand + VehPower + VehGas + Area + BonusMalus,
    family = quasipoisson(link = "log"),
    offset = log(Exposure),
    data   = train
  )
} else {
  freq_model_final <- freq_model
}

freq_relativities <- exp(coef(freq_model_final))
print(round(freq_relativities, 3))

# ---- Severity ----
# Vehicle age, fuel type, driver age and region tested and dropped: not
# significant. VehPower and Area retained.
#
# Large-loss cap: a small number of claims are big enough to swamp the GLM's
# factor effects rather than inform them. The severity model is fit on
# claims capped at the 99th percentile (the "attritional" layer); the
# capped-out amount is priced separately as a flat per-claim loading rather
# than asking VehPower/Area to explain it, since they don't.
sev_train <- train %>%
  filter(ClaimNb > 0, ClaimAmount > 0) %>%
  mutate(AvgClaim = ClaimAmount / ClaimNb)

large_loss_threshold <- quantile(sev_train$AvgClaim, 0.99)
cat("Large-loss threshold (99th percentile):", round(large_loss_threshold, 2), "\n")

sev_train <- sev_train %>%
  mutate(AvgClaim_capped = pmin(AvgClaim, large_loss_threshold))

sev_model <- glm(
  AvgClaim_capped ~ VehPower + Area,
  family  = Gamma(link = "log"),
  weights = ClaimNb,
  data    = sev_train
)
summary(sev_model)
cat("Severity model fitted on", nrow(sev_train), "claim records\n")

excess_loading <- mean(pmax(sev_train$AvgClaim - large_loss_threshold, 0))
cat("Flat excess loading per claim (above threshold):", round(excess_loading, 2), "\n")

# ---- Pure premium (frequency x severity; excludes expenses/profit/reinsurance) ----
test$pred_freq <- predict(freq_model_final, newdata = test, type = "response",
                          newoffset = log(rep(1, nrow(test))))
test$pred_sev  <- predict(sev_model, newdata = test, type = "response") + excess_loading
test$pure_prem <- test$pred_freq * test$pred_sev
summary(test$pure_prem)

baseline <- sum(train$ClaimAmount) / sum(train$Exposure)
cat("Baseline flat pure premium:", round(baseline, 2), "\n")

# ---- Lift chart ----
# Deciles split by cumulative EXPOSURE, not row count, so no decile is
# dominated by a handful of policies.
test <- test %>%
  arrange(pure_prem) %>%
  mutate(
    cum_exposure = cumsum(Exposure),
    decile       = pmin(10, ceiling(cum_exposure / sum(Exposure) * 10))
  )

lift <- test %>%
  group_by(decile) %>%
  summarise(
    exposure  = sum(Exposure),
    predicted = sum(pure_prem * Exposure) / sum(Exposure),
    actual    = sum(ClaimAmount) / sum(Exposure),
    .groups   = "drop"
  )
print(lift)

p_lift <- lift %>%
  tidyr::pivot_longer(cols = c(predicted, actual), names_to = "series", values_to = "value") %>%
  ggplot2::ggplot(ggplot2::aes(x = decile, y = value, colour = series)) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_point(size = 2) +
  ggplot2::labs(title = "Lift chart: predicted vs actual pure premium by decile",
                x = "Predicted premium decile (1 = lowest risk)",
                y = "Average cost per policy-year", colour = NULL) +
  ggplot2::theme_minimal()
ggplot2::ggsave("figures/lift_chart.png", p_lift, width = 7, height = 4, dpi = 150)

lift_ratio <- lift$actual[lift$decile == 10] / lift$actual[lift$decile == 1]
cat(sprintf("Decile 10 vs decile 1 actual cost ratio: %.2fx\n", lift_ratio))

# ---- Gini ----
# gini_curve() is in R/functions.R. The model's ranking by pure_prem is
# continuous with effectively no ties, so its Gini is stable on its own.
# The flat baseline is a genuine tie (every policy gets the same premium),
# so its Gini depends entirely on how ties get broken. With this dataset's
# heavy-tailed severity, one random tie-break can land the big claims
# anywhere and swing the "no information" Gini noticeably away from zero
# just by chance. Averaging over many random tie-breaks gives a stable
# estimate of what a genuinely uninformative model looks like.

actual_rate <- test$ClaimAmount / test$Exposure
model_gini  <- gini_curve(actual_rate, test$pure_prem, test$Exposure)

# Frequency alone, evaluated against claim count rather than claim cost.
# This isolates whether the part of the model with a clear factor effect
# (frequency) actually separates risk, independent of severity noise.
actual_freq_rate <- test$ClaimNb / test$Exposure
freq_only_gini   <- gini_curve(actual_freq_rate, test$pred_freq, test$Exposure)

set.seed(42)
n_reps <- 200
baseline_ginis <- replicate(n_reps, {
  shuffled <- sample(nrow(test))
  gini_curve(actual_rate[shuffled], rep(1, nrow(test)), test$Exposure[shuffled])
})
baseline_gini_mean <- mean(baseline_ginis)
baseline_gini_sd    <- sd(baseline_ginis)

cat(sprintf("Frequency-only Gini (claim count vs predicted frequency): %.4f\n", freq_only_gini))
cat(sprintf("Model Gini (pure premium): %.4f\n", model_gini))
cat(sprintf("Baseline (flat) Gini across %d random tie-breaks: mean %.4f, sd %.4f\n",
            n_reps, baseline_gini_mean, baseline_gini_sd))

saveRDS(freq_model_final, "output/freq_model.rds")
saveRDS(sev_model,        "output/sev_model.rds")
write.csv(lift, "output/lift_table.csv", row.names = FALSE)
write.csv(
  data.frame(metric = c("baseline_pure_premium", "lift_decile10_vs_decile1", "freq_only_gini",
                        "model_gini", "baseline_gini_mean", "baseline_gini_sd",
                        "large_loss_threshold", "excess_loading"),
             value  = c(baseline, lift_ratio, freq_only_gini, model_gini,
                        baseline_gini_mean, baseline_gini_sd,
                        large_loss_threshold, excess_loading)),
  "output/summary_metrics.csv", row.names = FALSE
)
