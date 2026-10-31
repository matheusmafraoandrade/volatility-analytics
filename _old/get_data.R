library(quantmod)

# SPY: S&P 500 ETF
data <- getSymbols(
  "SPY",
  src = "yahoo",
  periodicity = "1 minutes",
  auto.assign = FALSE
)

# Convert to data frame
data <- data.frame(
  datetime = index(data),
  coredata(data)
)

names(data) <- c(
  "datetime", "open", "high",
  "low", "close", "volume"
)

# View
head(data)
summary(data)

# Save
write.csv(
  data,
  "data/spy_intraday_1min.csv",
  row.names = FALSE
)
