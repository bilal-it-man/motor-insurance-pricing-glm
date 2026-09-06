# functions.R
# Shared helpers, sourced by 3_pricing_models.R and 4_cross_validation.R.

# Gini coefficient from a Lorenz curve: policies ranked by `predicted`,
# cumulative exposure on the x-axis, cumulative actual cost/frequency on
# the y-axis.
gini_curve <- function(actual, predicted, exposure) {
  ord <- order(predicted)
  actual_o   <- actual[ord]
  exposure_o <- exposure[ord]
  cum_exposure <- cumsum(exposure_o) / sum(exposure_o)
  cum_loss     <- cumsum(actual_o * exposure_o) / sum(actual_o * exposure_o)
  n <- length(cum_exposure)
  area <- sum((cum_exposure[-1] - cum_exposure[-n]) * (cum_loss[-1] + cum_loss[-n]) / 2)
  1 - 2 * area
}
