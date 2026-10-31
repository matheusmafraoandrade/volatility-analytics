library(rugarch)
library(ggplot2)

# 1. Load and prepare data
data <- read.csv("data/spy_intraday_1min.csv")
data$datetime <- as.POSIXct(data$datetime, tz = "UTC")
data <- data[order(data$datetime), ]

data$return <- c(NA, diff(log(data$close)))
data <- na.omit(data)
data$return_pct <- 100 * data$return

# 2. Split data chronologically
split_index <- floor(0.7 * nrow(data))
test <- data[(split_index + 1):nrow(data), ]

# 3. Specify GARCH(1,1)
spec <- ugarchspec(
  variance.model = list(
    model = "sGARCH",
    garchOrder = c(1, 1)
  ),
  mean.model = list(
    armaOrder = c(0, 0),
    include.mean = TRUE
  ),
  distribution.model = "std"
)

# 4. Estimate initial model
fit <- ugarchfit(
  spec = spec,
  data = data$return_pct[1:split_index],
  solver = "hybrid"
)

# 5. Initialize rolling backtest
refit_every <- 100
predictions <- numeric(nrow(test))

# Initial conditional variance and parameters
pars <- coef(fit)
omega <- pars["omega"]
alpha <- pars["alpha1"]
beta <- pars["beta1"]
mu <- pars["mu"]

last_residual <- tail(residuals(fit), 1)
last_variance <- tail(sigma(fit)^2, 1)

for (i in seq_len(nrow(test))) {
  
  # Re-estimate parameters every 100 observations
  if (i > 1 && (i - 1) %% refit_every == 0) {
    
    end_index <- split_index + i - 1
    
    fit <- ugarchfit(
      spec = spec,
      data = data$return_pct[1:end_index],
      solver = "hybrid"
    )
    
    pars <- coef(fit)
    omega <- pars["omega"]
    alpha <- pars["alpha1"]
    beta <- pars["beta1"]
    mu <- pars["mu"]
    
    last_residual <- tail(residuals(fit), 1)
    last_variance <- tail(sigma(fit)^2, 1)
  }
  
  # Forecast next-minute conditional variance
  next_variance <- omega +
    alpha * as.numeric(last_residual)^2 +
    beta * as.numeric(last_variance)
  
  predictions[i] <- sqrt(max(next_variance, 0))
  
  # Update the state using the observed return
  observed_return <- test$return_pct[i]
  last_residual <- observed_return - mu
  last_variance <- next_variance
}

# 6. Store forecasts
test$forecast_volatility <- predictions
test$realized_proxy <- abs(test$return_pct)

# 7. Evaluate forecasts
evaluation <- data.frame(
  MAE = mean(abs(
    test$realized_proxy - test$forecast_volatility
  )),
  RMSE = sqrt(mean(
    (test$realized_proxy - test$forecast_volatility)^2
  ))
)

print(evaluation)

# 8. Plot results
ggplot(test, aes(x = datetime)) +
  geom_line(aes(y = realized_proxy, color = "Absolute return")) +
  geom_line(aes(y = forecast_volatility, color = "GARCH forecast")) +
  labs(
    title = "GARCH(1,1) Out-of-Sample Backtest",
    x = "Time",
    y = "Value (% per minute)",
    color = "Series"
  ) +
  theme_minimal()

# 9. Save results
write.csv(
  test,
  "output/tables/spy_garch_backtest.csv",
  row.names = FALSE
)

write.csv(
  evaluation,
  "output/tables/spy_garch_metrics.csv",
  row.names = FALSE
)
