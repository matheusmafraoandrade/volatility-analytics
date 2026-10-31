# Shared settings and helper functions for the volatility forecasting project
# (HAR vs. GARCH vs. decision tree). Sourced at the top of every script;
# it is not meant to be run on its own.

library(here)

# ---- Settings ---------------------------------------------------------------
WINDOW      <- 500  # training window, in trading days (about 3 years)
REFIT_EVERY <- 20   # re-estimate each model every 20 trading days (about a month)
MIN_BARS    <- 200  # a full FX day has 288 five-minute bars; shorter days (holidays) are dropped
MAX_FLAT    <- 0.5  # also drop days where over half of the 5-minute returns are exactly zero (stale quotes)
ANN_FACTOR  <- 252  # trading days per year; only used to annualize volatility for reading/plots

for (folder in c("data/processed", "output/figures", "output/tables")) {
  dir.create(here(folder), recursive = TRUE, showWarnings = FALSE)
}

# ---- Data -------------------------------------------------------------------

# Read the raw 5-minute CSV. If no file is given, use the only CSV in data/raw/.
read_raw <- function(file = NULL) {
  if (is.null(file)) {
    files <- list.files(here("data", "raw"), pattern = "\\.csv$", full.names = TRUE)
    if (length(files) != 1) {
      stop("Put exactly one CSV file in data/raw/ (found ", length(files), ").")
    }
    file <- files
  }
  raw <- read.csv(file, stringsAsFactors = FALSE)
  names(raw) <- tolower(names(raw))
  time_col <- grep("time|date", names(raw))[1]

  if (is.numeric(raw[[time_col]])) {
    # Unix timestamp in milliseconds (format produced by dukascopy-node)
    time <- as.POSIXct(raw[[time_col]] / 1000, origin = "1970-01-01", tz = "UTC")
  } else {
    # Text timestamp such as "02.01.2020 00:00:00.000" (manual download from the website)
    time <- as.POSIXct(raw[[time_col]], format = "%d.%m.%Y %H:%M:%OS", tz = "UTC")
  }
  data.frame(time = time, close = raw$close)
}

# Basic cleaning: drop missing/invalid prices, duplicated timestamps and weekends.
clean_raw <- function(raw) {
  n_start <- nrow(raw)
  x <- raw[!is.na(raw$time) & !is.na(raw$close) & raw$close > 0, ]
  x <- x[!duplicated(x$time), ]
  x <- x[order(x$time), ]
  weekday <- as.integer(format(x$time, "%u"))  # 1 = Monday ... 7 = Sunday
  x <- x[weekday <= 5, ]
  message("Cleaning: ", n_start, " -> ", nrow(x), " five-minute bars")
  x
}

# Collapse 5-minute bars into one row per day:
#   ret = close-to-close daily return, in %          (what GARCH uses)
#   rv  = realized variance, in %^2: the sum of squared 5-minute returns within the day
make_daily <- function(bars) {
  by_day <- split(bars, as.Date(bars$time))
  daily <- data.frame(
    date   = as.Date(names(by_day)),
    n_bars = sapply(by_day, nrow),
    close  = sapply(by_day, function(x) x$close[nrow(x)]),
    rv     = sapply(by_day, function(x) sum(diff(log(x$close))^2)) * 100^2,
    share_flat = sapply(by_day, function(x) mean(diff(log(x$close)) == 0))
  )
  daily$ret <- c(NA, diff(log(daily$close))) * 100

  # Keep only days with real trading: enough bars, not mostly stale quotes, rv > 0.
  # (A day with rv = 0 would make the QLIKE loss infinite.)
  keep <- daily$n_bars >= MIN_BARS & daily$share_flat <= MAX_FLAT &
    daily$rv > 0 & !is.na(daily$ret)
  if (any(!keep)) {
    cat("\nDays dropped from the daily series:\n")
    print(daily[!keep, c("date", "n_bars", "share_flat", "rv")], row.names = FALSE)
  }
  message("Daily series: ", sum(keep), " days (", sum(!keep), " dropped)")
  daily <- daily[keep, ]
  rownames(daily) <- NULL
  daily[, c("date", "ret", "rv", "n_bars")]
}

# ---- Features and rolling forecasts -------------------------------------------

# One row per forecast origin (end of day t). The target is the realized variance
# of day t+1. The three regressors are the HAR ones: yesterday's RV, the 5-day
# average and the 22-day average (all including day t).
make_features <- function(daily) {
  n  <- nrow(daily)
  rv <- daily$rv
  mean_last <- function(k) {
    sapply(seq_len(n), function(i) if (i < k) NA else mean(rv[(i - k + 1):i]))
  }
  feat <- data.frame(
    day_index   = seq_len(n),
    rv_d        = rv,
    rv_w        = mean_last(5),
    rv_m        = mean_last(22),
    target      = c(rv[-1], NA),
    target_date = c(daily$date[-1], NA)
  )
  feat <- feat[complete.cases(feat), ]
  rownames(feat) <- NULL
  feat
}

# Rolling-window forecasts for models that use the features above (HAR, tree).
# fit_fun(train) must return a model that works with predict(model, newdata).
# Every REFIT_EVERY days the model is re-estimated on the last WINDOW observations
# whose target is already known, then used to forecast the next REFIT_EVERY days.
rolling_forecast <- function(feat, fit_fun) {
  n <- nrow(feat)
  if (n <= WINDOW) stop("Not enough data: need more than ", WINDOW, " days, have ", n)
  forecast <- rep(NA_real_, n)
  for (j in seq(WINDOW + 1, n, by = REFIT_EVERY)) {
    train <- feat[(j - WINDOW):(j - 1), ]  # targets up to day s_j: no look-ahead
    model <- fit_fun(train)
    rows  <- j:min(j + REFIT_EVERY - 1, n)
    forecast[rows] <- predict(model, newdata = feat[rows, ])
  }
  out <- data.frame(date = feat$target_date, forecast = forecast, realized = feat$target)
  out[!is.na(out$forecast), ]
}

# ---- Loss functions and test ------------------------------------------------------

# Loss of each single forecast (vectors), so they can be averaged or tested.
se_loss <- function(forecast, realized) (realized - forecast)^2

qlike_loss <- function(forecast, realized) {
  forecast <- pmax(forecast, 1e-6)  # QLIKE needs strictly positive forecasts
  realized / forecast - log(realized / forecast) - 1
}

# Diebold-Mariano test of equal predictive accuracy (one-step-ahead forecasts).
# Negative statistic: model A has the lower average loss (is better than B).
# Uses a Newey-West long-run variance because daily losses are autocorrelated.
dm_test <- function(loss_a, loss_b) {
  d    <- loss_a - loss_b
  n    <- length(d)
  lag  <- floor(n^(1 / 3))
  dbar <- mean(d)
  gamma <- function(k) sum((d[(k + 1):n] - dbar) * (d[1:(n - k)] - dbar)) / n
  weights <- 1 - seq_len(lag) / (lag + 1)
  lrv  <- gamma(0) + 2 * sum(weights * sapply(seq_len(lag), gamma))
  stat <- dbar / sqrt(lrv / n)
  c(statistic = stat, p_value = 2 * pnorm(-abs(stat)))
}
