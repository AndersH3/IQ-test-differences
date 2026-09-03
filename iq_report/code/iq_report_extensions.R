#!/usr/bin/env Rscript
# Supplement to iq_difference_simulation.R; base R only.
# Run from the report directory: Rscript code/iq_report_extensions.R
# Outputs: tables/*.csv and *.tex; figures/*.pdf; results/extensions.txt.
# These are model calculations, not empirical IQ-test observations.
options(width = 105)
for (d in c("tables", "figures", "results")) dir.create(d, showWarnings = FALSE)
sigma <- 15
rho <- 0.6
mu <- 100
k <- 5
s <- sigma * sqrt(1 - rho)
tau <- sqrt(2) * s
a <- function(k) integrate(function(z) k * z * dnorm(z) *
  exp((k - 1) * pnorm(z, log.p = TRUE)), -Inf, Inf,
  rel.tol = 1e-11, subdivisions = 2000)$value
a5 <- 5 / (4 * sqrt(pi)) * (1 + 6 / pi * asin(1 / 3))
# Finch (2016), section 1: second raw moment of the iid normal range.
nu5 <- 2 * (1 + 5 * sqrt(3) / (2 * pi) + 15 / pi^2 * acos(2 / 3) -
  5 * sqrt(3) / (2 * pi^2) * acos(1 / 4))
range_population_sd <- s * sqrt(nu5 - (2 * a5)^2)
folded_mean <- function(delta, tau) {
  2 * tau * dnorm(delta / tau) + delta * (2 * pnorm(delta / tau) - 1)
}
# Stable probability of a standard normal lying in [z,z+w].
normal_interval <- function(z, w) {
  pmax(0, ifelse(z > 0,
    pnorm(z, lower.tail = FALSE) - pnorm(z + w, lower.tail = FALSE),
    pnorm(z + w) - pnorm(z)))
}
range_cdf <- function(w, k = 5) {
  if (w <= 0) return(0)
  integrate(function(z) k * dnorm(z) * normal_interval(z, w)^(k - 1),
    -Inf, Inf, rel.tol = 1e-10, subdivisions = 2000)$value
}
save_table <- function(x, name, digits = 3) {
  write.csv(x, paste0("tables/", name, ".csv"), row.names = FALSE)
  rows <- vapply(seq_len(nrow(x)), function(i) {
    vals <- vapply(x, function(v) if (is.numeric(v))
      formatC(v[i], format = "f", digits = digits) else as.character(v[i]), "")
    paste0(paste(vals, collapse = " & "), " \\\\")
  }, "")
  writeLines(rows, paste0("tables/", name, ".tex"))
}
threshold <- c(5, 10, 15, 20, 30, 40)
probabilities <- data.frame(threshold = threshold,
  abs_exceeds_pct = 200 * pnorm(threshold / tau, lower.tail = FALSE),
  range_exceeds_pct = 100 * ptukey(threshold / s, k, Inf, lower.tail = FALSE))
save_table(probabilities, "probabilities", 2)
probs <- c(.025, .5, .95, .975)
quantiles <- data.frame(probability = probs,
  absolute_difference = tau * qnorm((1 + probs) / 2),
  range = s * qtukey(probs, k, Inf))
save_table(quantiles, "quantiles", 3)
rhos <- c(0, .3, .5, .6, .7, .8, .9, 1)
save_table(data.frame(rho = rhos,
  expected_abs = 2 * sigma * sqrt((1 - rhos) / pi),
  expected_range = 2 * sigma * sqrt(1 - rhos) * a5), "sensitivity", 3)
ks <- c(2, 3, 4, 5, 10, 20, 50)
save_table(data.frame(k = ks, expected_normal_max = vapply(ks, a, 0.0),
  expected_range = 2 * s * vapply(ks, a, 0.0)), "test_counts", 3)
xs <- c(85, 100, 115, 130, 145, 160)
delta <- (1 - rho) * (xs - mu)
save_table(data.frame(first_score = xs, expected_second = mu + rho * (xs - mu),
  expected_signed_gap = delta,
  expected_abs_gap = folded_mean(delta, sigma * sqrt(1 - rho^2))),
  "conditional", 3)
shrinks <- c(sqrt(.4), .4, mean(sqrt(1 - c(.3, .9))))
models <- data.frame(model = c("Joint Gaussian", "Shared/independent mixture",
  "Smooth Gaussian mixture"), expected_abs = 2 * sigma * shrinks / sqrt(pi),
  expected_range = 2 * sigma * shrinks * a5)
save_table(models, "models", 4)
# Independent numerical path for the normal range distribution and quantiles.
cdf_checks <- data.frame(w = c(.5, 1, 2, 3, 4, 6))
cdf_checks$integral <- vapply(cdf_checks$w, range_cdf, 0.0)
cdf_checks$ptukey <- ptukey(cdf_checks$w, 5, Inf)
cdf_checks$error <- cdf_checks$integral - cdf_checks$ptukey
q_integral <- vapply(probs, function(p)
  uniroot(function(w) range_cdf(w) - p, c(.001, 12), tol = 1e-9)$root, 0.0)
stopifnot(max(abs(cdf_checks$error)) < 1e-7,
  max(abs(q_integral - qtukey(probs, 5, Inf))) < 2e-4,
  abs(a(5) - a5) < 1e-10,
  abs(range_cdf(2, 2) - (2 * pnorm(sqrt(2)) - 1)) < 1e-9)

# A separate seeded experiment with continuous, positive-density components.
# Equal probabilities of equicorrelation .3 and .9 give marginal rho=.6.
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(20260904)
n <- 1000000
batch <- 100000
acc_n <- 0
acc_mean <- c(abs = 0, range = 0)
acc_M2 <- c(abs = 0, range = 0)
for (b in seq_len(n / batch)) {
  r <- ifelse(runif(batch) < .5, .3, .9)
  common <- rnorm(batch)
  e <- matrix(rnorm(batch * 5), batch, 5)
  x <- sigma * (sqrt(r) * common + sqrt(1 - r) * e)
  lo <- hi <- x[, 1]
  for (j in 2:5) { lo <- pmin(lo, x[, j]); hi <- pmax(hi, x[, j]) }
  y <- cbind(abs = abs(x[, 1] - x[, 2]), range = hi - lo)
  mb <- colMeans(y)
  m2b <- colSums(sweep(y, 2, mb)^2)
  nt <- acc_n + batch
  change <- mb - acc_mean
  acc_M2 <- acc_M2 + m2b + change^2 * acc_n * batch / nt
  acc_mean <- acc_mean + change * batch / nt
  acc_n <- nt
}
smooth_check <- data.frame(quantity = names(acc_mean),
  theory = unlist(models[3, 2:3], use.names = FALSE),
  simulated = unname(acc_mean), mc_se = sqrt(acc_M2 / (n - 1) / n))
smooth_check$z_mc <- (smooth_check$simulated - smooth_check$theory) / smooth_check$mc_se
save_table(smooth_check, "smooth_simulation", 5)

# Full admissible equicorrelation range: verify projection construction.
H <- diag(5) - matrix(1 / 5, 5, 5)
for (r in c(-.25, -.1, 0, .6, 1)) {
  L <- cbind(sqrt(1 - r) * H, rep(sqrt((1 + 4 * r) / 5), 5))
  target <- (1 - r) * diag(5) + r * matrix(1, 5, 5)
  stopifnot(max(abs(tcrossprod(L) - target)) < 1e-12)
}
# Average correlation .6, but four pairs at 1 and six pairs at 1/3.
average_rho <- (4 + 6 / 3) / 10
block_range <- 2 * sigma * sqrt((1 - 1 / 3) / pi)

# Reproducible vector figures with the exact plotted values in CSV.
ink <- "#183B56"
teal <- "#007F82"
rust <- "#B75D36"
rho_grid <- seq(0, 1, length.out = 201)
rho_curve <- data.frame(rho = rho_grid,
  absolute = 2 * sigma * sqrt((1 - rho_grid) / pi),
  range = 2 * sigma * sqrt(1 - rho_grid) * a5)
write.csv(rho_curve, "tables/figure_correlation.csv", row.names = FALSE)
cairo_pdf("figures/correlation.pdf", width = 7.1, height = 3.7, family = "DejaVu Sans")
par(mar = c(4.2, 4.4, 1.1, .7), las = 1)
matplot(rho_curve$rho, rho_curve[, 2:3], type = "l", lty = 1, lwd = 2.2,
  col = c(teal, ink), xlab = "Common pairwise correlation",
  ylab = "Expected gap (IQ points)", ylim = c(0, 36))
abline(v = .6, col = "#AAAAAA", lty = 3)
points(rep(.6, 2), c(2 * s / sqrt(pi), 2 * s * a5), pch = 19, col = c(teal, ink))
legend("topright", c("Five-test range", "Two-test absolute difference"),
  col = c(ink, teal), lwd = 2.2, bty = "n", cex = .88)
dev.off()
q <- seq(0, 60, by = .15)
distributions <- data.frame(gap = q, abs_cdf = 2 * pnorm(q / tau) - 1,
  range_cdf = ptukey(q / s, 5, Inf))
write.csv(distributions, "tables/figure_distributions.csv", row.names = FALSE)
cairo_pdf("figures/distributions.pdf", width = 7.1, height = 3.7, family = "DejaVu Sans")
par(mar = c(4.2, 4.4, 1.1, .7), las = 1)
matplot(q, distributions[, 2:3], type = "l", lty = 1, lwd = 2.2,
  col = c(teal, ink), xlab = "Gap (IQ points)",
  ylab = "Probability gap is at most this value", ylim = c(0, 1))
abline(h = c(.025, .5, .975), col = "#DDDDDD", lty = 3)
legend("bottomright", c("Two-test absolute difference", "Five-test range"),
  col = c(teal, ink), lwd = 2.2, bty = "n", cex = .88)
dev.off()
xg <- seq(60, 160, by = .5)
d <- (1 - rho) * (xg - mu)
conditional_curve <- data.frame(first_score = xg,
  signed_gap = d, absolute_gap = folded_mean(d, sigma * sqrt(1 - rho^2)))
write.csv(conditional_curve, "tables/figure_conditional.csv", row.names = FALSE)
cairo_pdf("figures/conditional.pdf", width = 7.1, height = 3.7, family = "DejaVu Sans")
par(mar = c(4.2, 4.4, 1.1, .7), las = 1)
matplot(xg, conditional_curve[, 2:3], type = "l", lty = c(2, 1), lwd = 2.2,
  col = c(rust, ink), xlab = "Observed first-test score",
  ylab = "Conditional expected gap (IQ points)")
abline(h = 0, col = "#DDDDDD")
legend("topleft", c("Absolute gap", "Signed gap: first minus second"),
  col = c(ink, rust), lty = c(1, 2), lwd = 2.2, bty = "n", cex = .85)
dev.off()
log <- capture.output({
  cat("REPORT EXTENSIONS: 2026-09-03\nNo empirical IQ data used.\n\n")
  cat("Range CDF quadrature checks:\n"); print(cdf_checks, digits = 12)
  cat("Range quantiles from independent CDF inversion (IQ points):\n")
  print(data.frame(probability = probs, quantile = s * q_integral), digits = 12)
  cat("\nSmooth mixture: n=1000000, seed=20260904\n")
  print(smooth_check, digits = 10, row.names = FALSE)
  cat("\nAverage-correlation block model: average rho =", average_rho,
    "; expected range =", block_range, "\n")
  cat("\nProjection covariance checks passed at rho=-.25,-.1,0,.6,1.\n")
  cat("Exact iid normal range second moment =", format(nu5, digits = 12), "\n")
  cat("Gaussian five-test range population SD =",
    format(range_population_sd, digits = 12), "\n")
  cat("Theoretical range MC SE at n=5000000 =",
    format(range_population_sd / sqrt(5000000), digits = 12), "\n")
  cat("\nAll deterministic checks passed.\n\n"); print(sessionInfo())
})
writeLines(log, "results/extensions.txt")
cat(paste(log, collapse = "\n"), "\n")
