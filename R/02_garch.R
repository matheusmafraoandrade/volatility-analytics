# 02_garch.R
# GARCH(1,1) estimated on DAILY RETURNS only (no intraday information):
#   r_t        = mu + e_t,   e_t = sigma_t * z_t,   z_t ~ N(0, 1)
#   sigma2_t+1 = omega + alpha * e_t^2 + beta * sigma2_t
#
# Part A: fit on the full sample to read and interpret the parameters.
# Part B: rolling out-of-sample forecasts of next-day variance (used in 04_tree.R).

source(here::here("R", "functions.R"))
library(rugarch)

daily <- readRDS(here("data", "processed", "daily.rds"))
feat  <- make_features(daily)  # used only so GARCH forecasts the same days as HAR and tree

spec <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "norm"
)

# ---- Part A: full-sample fit ---------------------------------------------------------
fit <- ugarchfit(spec, data = daily$ret)
write.csv(fit@fit$robust.matcoef, here("output", "tables", "garch_parameters.csv"))
cat("\nGARCH(1,1) parameters (robust standard errors):\n")
print(round(fit@fit$robust.matcoef, 4))

p <- coef(fit)
persistence <- unname(p["alpha1"] + p["beta1"])
uncond_var  <- unname(p["omega"]) / (1 - persistence)
summary_tab <- data.frame(
  quantity = c("persistence (alpha + beta)",
               "shock half-life (days)",
               "unconditional variance (%^2)",
               "unconditional annualized volatility (%)"),
  value = c(persistence,
            log(0.5) / log(persistence),
            uncond_var,
            sqrt(ANN_FACTOR * uncond_var))
)
write.csv(summary_tab, here("output", "tables", "garch_summary.csv"), row.names = FALSE)
cat("\nDerived quantities:\n")
print(summary_tab, row.names = FALSE, digits = 4)

# Fitted GARCH volatility vs. realized volatility: GARCH is smooth, RV is noisy but unbiased.
png(here("output", "figures", "garch_in_sample.png"), width = 900, height = 450)
plot(daily$date, sqrt(ANN_FACTOR * daily$rv), type = "l", col = "grey70",
     xlab = "", ylab = "Annualized volatility (%)",
     main = "In-sample GARCH volatility vs. realized volatility")
lines(daily$date, sqrt(ANN_FACTOR) * as.numeric(sigma(fit)), col = "red", lwd = 1.5)
legend("topright", c("Realized (5-min)", "GARCH(1,1)"), col = c("grey70", "red"), lty = 1, bty = "n")
dev.off()

# ---- Part B: rolling out-of-sample forecasts -------------------------------------------
# Every REFIT_EVERY days we re-estimate GARCH on the last WINDOW daily returns.
# In between, the parameters stay fixed and we update the variance recursion by hand
# with each new return: sigma2_next = omega + alpha * (r - mu)^2 + beta * sigma2.
n        <- nrow(feat)
forecast <- rep(NA_real_, n)

for (j in seq(WINDOW + 1, n, by = REFIT_EVERY)) {
  s     <- feat$day_index[j]  # forecast origin: end of day s
  fit_j <- ugarchfit(spec, data = daily$ret[(s - WINDOW + 1):s], solver = "hybrid")
  p     <- coef(fit_j)
  sigma2 <- unname(tail(as.numeric(sigma(fit_j)), 1))^2  # model variance of day s

  for (k in j:min(j + REFIT_EVERY - 1, n)) {
    s_k    <- feat$day_index[k]
    sigma2 <- unname(p["omega"] + p["alpha1"] * (daily$ret[s_k] - p["mu"])^2 +
                       p["beta1"] * sigma2)
    forecast[k] <- sigma2  # forecast of the variance of day s_k + 1
  }
}

out <- data.frame(date = feat$target_date, forecast = forecast, realized = feat$target)
out <- out[!is.na(out$forecast), ]
saveRDS(out, here("data", "processed", "forecasts_garch.rds"))

message("Done: ", nrow(out), " out-of-sample forecasts. Next: run 03_har.R")
