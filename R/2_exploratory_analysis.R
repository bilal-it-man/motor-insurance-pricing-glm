# 2_exploratory_analysis.R
# Frequency by rating factor + severity distribution, saved to figures/.

library(dplyr)
library(ggplot2)

dir.create("figures", showWarnings = FALSE)

train <- readRDS("data/train.rds")

# Exposure-weighted frequency by factor level, dropping thin levels (<100 exposure).
freq_by <- function(df, var) {
  df %>%
    group_by(level = .data[[var]]) %>%
    summarise(exposure = sum(Exposure), claims = sum(ClaimNb),
              frequency = claims / exposure, .groups = "drop") %>%
    filter(exposure > 100)
}

plot_freq_by <- function(df, var, title, filename) {
  p <- freq_by(df, var) %>%
    ggplot(aes(x = level, y = frequency)) +
    geom_col(fill = "grey30") +
    labs(title = title, x = var, y = "Claims per policy-year") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path("figures", filename), p, width = 7, height = 4, dpi = 150)
  p
}

rating_factors <- list(
  DrivAgeBand = "Observed claim frequency by driver age band",
  VehAgeBand  = "Observed claim frequency by vehicle age band",
  VehPower    = "Observed claim frequency by vehicle power",
  VehGas      = "Observed claim frequency by fuel type",
  Area        = "Observed claim frequency by area density code",
  Region      = "Observed claim frequency by region"
)

for (var in names(rating_factors)) {
  plot_freq_by(train, var, rating_factors[[var]], paste0("freq_by_", tolower(var), ".png"))
}

train_binned <- train %>%
  mutate(
    BonusMalusBand = cut(BonusMalus, breaks = c(0, 60, 80, 100, 120, 150, Inf), right = FALSE),
    DensityBand    = cut(Density, breaks = quantile(Density, probs = seq(0, 1, 0.2)), include.lowest = TRUE)
  )

plot_freq_by(train_binned, "BonusMalusBand", "Observed claim frequency by bonus-malus band", "freq_by_bonusmalus.png")
plot_freq_by(train_binned, "DensityBand", "Observed claim frequency by population density quintile", "freq_by_density.png")

claims_only <- train %>% filter(ClaimNb > 0, ClaimAmount > 0)

p_sev <- ggplot(claims_only, aes(x = ClaimAmount)) +
  geom_histogram(bins = 60, fill = "grey30") +
  scale_x_log10() +
  labs(title = "Distribution of claim amounts (log scale)", x = "Claim amount", y = "Count") +
  theme_minimal()
ggsave("figures/severity_distribution.png", p_sev, width = 7, height = 4, dpi = 150)

cat(sprintf("Policies with at least one claim: %.2f%%\n", mean(train$ClaimNb > 0) * 100))
cat(sprintf("Overall frequency: %.4f claims per policy-year\n", sum(train$ClaimNb) / sum(train$Exposure)))
cat(sprintf("Claims used for severity: %d\n", nrow(claims_only)))
