# 03_har.R
# HAR-RV (Corsi, 2009): a linear regression of tomorrow's realized variance on
# realized variance averaged over three horizons.
#   RV_t+1 = b0 + b_d * RV_t + b_w * mean(RV_t-4..t) + b_m * mean(RV_t-21..t) + e
# Idea: traders act at different horizons (daily, weekly, monthly), and each group
# reacts to volatility at its own horizon. The sum is a simple long-memory pattern.
#
# Part A: full-sample fit, to read the coefficients.
# Part B: rolling out-of-sample forecasts (used in 04_tree.R).

source(here::here("R", "functions.R"))

daily <- readRDS(here("data", "processed", "daily.rds"))
feat  <- make_features(daily)

# ---- Part A: full-sample fit ---------------------------------------------------------
har <- lm(target ~ rv_d + rv_w + rv_m, data = feat)
har_summary <- summary(har)

# Note: these OLS standard errors ignore autocorrelation in the residuals, so read the
# size of the coefficients and the R-squared, not the t-statistics. (Newey-West errors,
# e.g. sandwich::NeweyWest(), would be the natural upgrade.)
write.csv(har_summary$coefficients, here("output", "tables", "har_coefficients.csv"))
cat("\nHAR-RV coefficients:\n")
print(round(har_summary$coefficients, 4))

b <- coef(har)
summary_tab <- data.frame(
  quantity = c("sum of slopes (b_d + b_w + b_m)", "R-squared (in sample)"),
  value    = c(sum(b[c("rv_d", "rv_w", "rv_m")]), har_summary$r.squared)
)
write.csv(summary_tab, here("output", "tables", "har_summary.csv"), row.names = FALSE)
cat("\nDerived quantities:\n")
print(summary_tab, row.names = FALSE, digits = 4)

# ---- Part B: rolling out-of-sample forecasts -------------------------------------------
fit_har <- function(train) lm(target ~ rv_d + rv_w + rv_m, data = train)
out <- rolling_forecast(feat, fit_har)
saveRDS(out, here("data", "processed", "forecasts_har.rds"))

message("Done: ", nrow(out), " out-of-sample forecasts. Next: run 04_tree.R")
