# 1_data_prep.R
# Loads freMTPL2freq/freMTPL2sev from CASdatasets, cleans, joins, splits.

if (!requireNamespace("xts", quietly = TRUE)) install.packages("xts")
if (!requireNamespace("zoo", quietly = TRUE)) install.packages("zoo")
if (!requireNamespace("CASdatasets", quietly = TRUE)) {
  install.packages("CASdatasets", repos = "https://cas.uqam.ca/pub/", type = "source")
  # mirrors if that's down:
  # repos = "https://dutangc.perso.math.cnrs.fr/RRepository/pub/"
  # repos = "http://dutangc.free.fr/pub/RRepos/pub/"
}

library(CASdatasets)
library(dplyr)
library(ggplot2)

dir.create("data", showWarnings = FALSE)

data(freMTPL2freq)
data(freMTPL2sev)

freq <- freMTPL2freq
sev  <- freMTPL2sev

# Exposure >1 and ClaimNb outliers are data artefacts; capped not dropped.
# The cluster of ClaimAmount == 1200 in freMTPL2sev is a known dataset
# artefact but kept, since genuine 1200 claims can't be distinguished from it.
freq_clean <- freq %>%
  mutate(
    Exposure = pmin(Exposure, 1),
    ClaimNb  = pmin(as.numeric(ClaimNb), 4)
  ) %>%
  filter(Exposure > 0)

sev_agg <- sev %>%
  group_by(IDpol) %>%
  summarise(ClaimAmount = sum(ClaimAmount), .groups = "drop")

dat <- freq_clean %>%
  left_join(sev_agg, by = "IDpol") %>%
  mutate(ClaimAmount = ifelse(is.na(ClaimAmount), 0, ClaimAmount))

dat <- dat %>%
  mutate(
    DrivAgeBand = cut(DrivAge, breaks = c(17, 21, 26, 31, 41, 51, 71, Inf), right = FALSE),
    VehAgeBand  = cut(VehAge,  breaks = c(0, 1, 5, 10, 15, Inf), right = FALSE),
    VehPower    = as.factor(VehPower),
    VehBrand    = as.factor(VehBrand),
    VehGas      = as.factor(VehGas),
    Region      = as.factor(Region),
    Area        = as.factor(Area)
  )

set.seed(42)
n <- nrow(dat)
train_idx <- sample(seq_len(n), size = floor(0.8 * n))
train <- dat[train_idx, ]
test  <- dat[-train_idx, ]

cat("Rows:", n, "| train:", nrow(train), "| test:", nrow(test), "\n")
cat("Total exposure (train):", round(sum(train$Exposure), 1), "policy-years\n")
cat("Overall claim frequency (train):",
    round(sum(train$ClaimNb) / sum(train$Exposure), 4), "claims per policy-year\n")

saveRDS(train, "data/train.rds")
saveRDS(test,  "data/test.rds")
