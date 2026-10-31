# 04_tree.R
# Regression tree on the SAME three regressors as HAR (RV_t, 5-day and 22-day averages).
# A tree splits days into groups using thresholds (e.g. "if the 22-day average RV is
# above x...") and predicts the average next-day RV of each group. Compared with HAR
# the only difference is the functional form: linear vs. piecewise constant (nonlinear).
#
# Part A: full-sample tree, to read the splits.
# Part B: rolling out-of-sample forecasts.
# Part C: comparison of GARCH, HAR and tree (needs 02_garch.R and 03_har.R done first).

source(here::here("R", "functions.R"))
library(rpart)

set.seed(1)  # rpart uses random folds for its internal cross-validation

daily <- readRDS(here("data", "processed", "daily.rds"))
feat  <- make_features(daily)

# Grow a deep tree, then prune it back to the size with the lowest cross-validated
# error. (The random CV folds ignore time order; fine for choosing the tree size.)
fit_tree <- function(train) {
  big <- rpart(target ~ rv_d + rv_w + rv_m, data = train, method = "anova",
               control = rpart.control(minbucket = 20, cp = 0.001))
  best_cp <- big$cptable[which.min(big$cptable[, "xerror"]), "CP"]
  prune(big, cp = best_cp)
}

# ---- Part A: full-sample tree --------------------------------------------------------
tree <- fit_tree(feat)
n_leaves <- sum(tree$frame$var == "<leaf>")
cat("\nFull-sample tree: ", n_leaves, " leaves\n", sep = "")

if (n_leaves > 1) {
  png(here("output", "figures", "tree_full_sample.png"), width = 900, height = 600)
  plot(tree, uniform = TRUE, margin = 0.1, main = "Regression tree for next-day realized variance")
  text(tree, use.n = TRUE, cex = 0.8)
  dev.off()

  importance <- data.frame(variable = names(tree$variable.importance),
                           importance = unname(tree$variable.importance))
  write.csv(importance, here("output", "tables", "tree_importance.csv"), row.names = FALSE)
  cat("\nVariable importance:\n")
  print(importance, row.names = FALSE, digits = 3)
} else {
  cat("The pruned tree is just the root: no split improves the forecast.\n")
}

# ---- Part B: rolling out-of-sample forecasts -------------------------------------------
out <- rolling_forecast(feat, fit_tree)
saveRDS(out, here("data", "processed", "forecasts_tree.rds"))

# ---- Part C: compare the three models ------------------------------------------------
files <- c(GARCH = "forecasts_garch.rds", HAR = "forecasts_har.rds", Tree = "forecasts_tree.rds")
missing_files <- files[!file.exists(here("data", "processed", files))]
if (length(missing_files) > 0) {
  stop("Missing forecast files: ", paste(missing_files, collapse = ", "),
       ". Run 02_garch.R and 03_har.R first.")
}
fc <- lapply(files, function(f) readRDS(here("data", "processed", f)))
stopifnot(identical(fc$GARCH$date, fc$HAR$date), identical(fc$GARCH$date, fc$Tree$date))

if (any(fc$GARCH$realized <= 0)) {
  stop("Some realized variances are zero: QLIKE would be infinite. Re-run 01_data_prep.R.")
}

models    <- names(files)
forecasts <- data.frame(date = fc$GARCH$date, realized = fc$GARCH$realized,
                        GARCH = fc$GARCH$forecast, HAR = fc$HAR$forecast, Tree = fc$Tree$forecast)

# 1. Average losses (lower is better).
loss_se    <- sapply(models, function(m) se_loss(forecasts[[m]], forecasts$realized))
loss_qlike <- sapply(models, function(m) qlike_loss(forecasts[[m]], forecasts$realized))
results <- data.frame(model = models, MSE = colMeans(loss_se), QLIKE = colMeans(loss_qlike),
                      row.names = NULL)
write.csv(results, here("output", "tables", "loss_comparison.csv"), row.names = FALSE)
cat("\nOut-of-sample forecasts: ", nrow(forecasts), " days (",
    format(min(forecasts$date)), " to ", format(max(forecasts$date)), ")\n", sep = "")
cat("\nAverage losses (lower is better):\n")
print(results, row.names = FALSE, digits = 4)

# 2. Diebold-Mariano tests. Negative statistic: model A is better than model B.
pairs <- list(c("GARCH", "HAR"), c("GARCH", "Tree"), c("HAR", "Tree"))
dm <- do.call(rbind, lapply(pairs, function(pr) {
  do.call(rbind, lapply(c("MSE", "QLIKE"), function(loss) {
    L <- if (loss == "MSE") loss_se else loss_qlike
    t <- dm_test(L[, pr[1]], L[, pr[2]])
    data.frame(model_A = pr[1], model_B = pr[2], loss = loss,
               statistic = t[["statistic"]], p_value = t[["p_value"]])
  }))
}))
write.csv(dm, here("output", "tables", "diebold_mariano.csv"), row.names = FALSE)
cat("\nDiebold-Mariano tests (negative statistic: A better than B):\n")
print(dm, row.names = FALSE, digits = 3)

# 3. Mincer-Zarnowitz regressions: realized = a + b * forecast + error.
# A good forecast has a close to 0 and b close to 1 (unbiased); R-squared measures fit.
mz <- do.call(rbind, lapply(models, function(m) {
  f <- lm(forecasts$realized ~ forecasts[[m]])
  data.frame(model = m, intercept = coef(f)[[1]], slope = coef(f)[[2]],
             r_squared = summary(f)$r.squared)
}))
write.csv(mz, here("output", "tables", "mincer_zarnowitz.csv"), row.names = FALSE)
cat("\nMincer-Zarnowitz regressions:\n")
print(mz, row.names = FALSE, digits = 3)

# 4. Figure: last 250 forecast days, in annualized volatility.
last <- tail(forecasts, 250)
ann  <- function(x) sqrt(ANN_FACTOR * pmax(x, 0))
png(here("output", "figures", "forecasts_last_year.png"), width = 900, height = 450)
plot(last$date, ann(last$realized), type = "l", col = "grey70", lwd = 2,
     xlab = "", ylab = "Annualized volatility (%)",
     main = "Next-day volatility: realized vs. forecasts (last 250 days)")
lines(last$date, ann(last$GARCH), col = "red")
lines(last$date, ann(last$HAR),   col = "blue")
lines(last$date, ann(last$Tree),  col = "orange")
legend("topright", c("Realized", models), col = c("grey70", "red", "blue", "orange"),
       lty = 1, bty = "n")
dev.off()

message("Done. All tables are in output/tables/ and figures in output/figures/")
