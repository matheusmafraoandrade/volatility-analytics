# 01_data_prep.R
# Goal: turn raw 5-minute prices into two daily series
#   ret = daily close-to-close return (%)   -> the only input GARCH uses
#   rv  = daily realized variance (%^2)     -> input of HAR and tree, and the
#                                              "truth" that all models try to forecast
# It also shows WHY intraday data is more informative than daily returns
# (the daily squared return is a very noisy measure of that day's variance).
#
# Input : one CSV with 5-minute EUR/USD candles in data/raw/
# Output: data/processed/daily.rds, plus tables and figures in output/

source(here::here("R", "functions.R"))

# ---- 1. Clean the raw data and build the daily series --------------------------
bars  <- clean_raw(read_raw())
daily <- make_daily(bars)
saveRDS(daily, here("data", "processed", "daily.rds"))

# Sanity check for bad ticks: the largest 5-minute moves should look plausible.
five_min_ret <- 100 * abs(unlist(lapply(split(bars$close, as.Date(bars$time)),
                                        function(x) diff(log(x))), use.names = FALSE))
cat("\nLargest 5-minute absolute returns within a day (%):\n")
print(round(sort(five_min_ret, decreasing = TRUE)[1:5], 3))

# ---- 2. Descriptive statistics ---------------------------------------------------
stats <- data.frame(
  statistic = c("number of days", "first day", "last day",
                "mean daily return (%)", "sd of daily return (%)",
                "mean realized variance (%^2)",
                "annualized volatility from RV (%)"),
  value = c(nrow(daily), as.character(min(daily$date)), as.character(max(daily$date)),
            round(mean(daily$ret), 4), round(sd(daily$ret), 4),
            round(mean(daily$rv), 4),
            round(sqrt(ANN_FACTOR * mean(daily$rv)), 2))
)
write.csv(stats, here("output", "tables", "descriptive_stats.csv"), row.names = FALSE)
cat("\nDescriptive statistics:\n")
print(stats, row.names = FALSE)

# ---- 3. Limitation of daily data: squared return vs. realized variance -------------
# Both are measures of the same unobserved quantity (the day's variance), but the
# squared return comes from ONE observation per day, the realized variance from ~288.
daily$ret2 <- daily$ret^2
acf1 <- function(x) acf(x, plot = FALSE)$acf[2]

proxy <- data.frame(
  measure       = c("squared daily return", "realized variance (5-min)"),
  mean          = c(mean(daily$ret2), mean(daily$rv)),
  sd            = c(sd(daily$ret2), sd(daily$rv)),
  lag1_autocorr = c(acf1(daily$ret2), acf1(daily$rv))
)
proxy$correlation_with_other <- cor(daily$ret2, daily$rv)
write.csv(proxy, here("output", "tables", "proxy_comparison.csv"), row.names = FALSE)
cat("\nSquared return vs. realized variance:\n")
print(proxy, row.names = FALSE, digits = 3)

# ---- 4. Figures --------------------------------------------------------------------
png(here("output", "figures", "realized_volatility.png"), width = 900, height = 450)
plot(daily$date, sqrt(ANN_FACTOR * daily$rv), type = "l",
     xlab = "", ylab = "Annualized realized volatility (%)",
     main = "EUR/USD realized volatility from 5-minute returns")
dev.off()

png(here("output", "figures", "return_vs_realized_vol.png"), width = 800, height = 600)
plot(sqrt(daily$rv), abs(daily$ret), pch = 16, col = rgb(0, 0, 0, 0.3),
     xlab = "Realized volatility from 5-minute data (% per day)",
     ylab = "Absolute daily return (% per day)",
     main = "Same-day volatility: two measures")
abline(lm(abs(ret) ~ sqrt(rv), data = daily), col = "red", lwd = 2)
dev.off()

message("Done. Next: run 02_garch.R")
