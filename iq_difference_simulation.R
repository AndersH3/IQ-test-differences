#!/usr/bin/env Rscript
# Expected IQ-score differences and ranges, with reproducible Monte Carlo checks.
# Base R only; no package installation is needed.
#
# Run:
#   Rscript iq_difference_simulation.R
#   Rscript iq_difference_simulation.R --n=1000000 --rho=0.6 --sd=15
#   Rscript iq_difference_simulation.R --k=10 --out=iq_10_test_results.txt
#
# Each row represents an independent person taking the same k tests.
# Gaussian model: X_j = mu + sigma*(sqrt(rho)*Z0 + sqrt(1-rho)*Zj).
# Z0,Z1,...,Zk are mutually independent standard normal variables.
# Thus every test has mean mu, SD sigma, and each pair has correlation rho.
#
# E|X1-X2| = 2*sigma*sqrt((1-rho)/pi).
# E(range_k) = 2*sigma*sqrt(1-rho)*k*integral z*phi(z)*Phi(z)^(k-1) dz.
# For k=5, E(max Zj) = 5/(4*sqrt(pi))*(1+6/pi*asin(1/3)).
#
# A second model demonstrates that even normal marginals and the entire
# correlation matrix do NOT determine these expectations without joint
# normality. With probability rho all scores share one normal draw;
# otherwise they are independent. Marginals/correlations are unchanged,
# but expected gaps/ranges scale by (1-rho), rather than sqrt(1-rho).
# This is a mathematical counterexample, not a proposed model of IQ tests.
#
# Numerical references (accessed 2026-09-03):
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Normal.html
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/integrate.html
# https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Tukey.html
# https://arxiv.org/abs/2004.04682 (normal orthant formula used in closed form)

options(width = 120)

parse_args <- function(args) {
  opt <- list(n = 5000000, rho = 0.6, sd = 15, mu = 100,
              k = 5, seed = 20260903, batch = 100000,
              out = "iq_simulation_results.txt")
  if ("--help" %in% args) {
    cat("Usage: Rscript iq_difference_simulation.R [--name=value ...]\n",
        "Options: --n=5000000 --rho=0.6 --sd=15 --mu=100 --k=5\n",
        "         --seed=20260903 --batch=100000 --out=iq_simulation_results.txt\n",
        "rho must be in [0,1]; the simulation uses a nonnegative common factor.\n",
        "The program refuses to overwrite an existing output file.\n", sep = "")
    quit(status = 0)
  }
  for (arg in args) {
    if (!grepl("^--[^=]+=.+$", arg)) stop("Use --name=value: ", arg)
    name <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=", "", arg)
    if (!name %in% names(opt)) stop("Unknown option: ", name)
    opt[[name]] <- if (name == "out") value else suppressWarnings(as.numeric(value))
  }
  num <- unlist(opt[setdiff(names(opt), "out")])
  if (any(!is.finite(num))) stop("All numerical options must be finite.")
  for (key in c("n", "k", "batch", "seed")) {
    if (opt[[key]] != floor(opt[[key]])) stop(key, " must be an integer.")
  }
  if (opt$n < 2 || opt$n > 1e9) stop("n must be between 2 and 1e9.")
  if (opt$k < 2 || opt$k > 1000) stop("k must be between 2 and 1000.")
  if (opt$batch < 1) stop("batch must be positive.")
  if (opt$seed < 0 || opt$seed > .Machine$integer.max) stop("Invalid seed.")
  if (opt$rho < 0 || opt$rho > 1) stop("rho must be in [0,1].")
  if (opt$sd <= 0) stop("sd must be positive.")
  if (file.exists(opt$out)) stop("Output already exists; choose another --out: ", opt$out)
  if (!dir.exists(dirname(opt$out))) stop("The output directory does not exist.")
  opt
}

normal_max_mean <- function(k) {
  # Use infinite integration limits. Phi(z)^(k-1) is computed on the log scale.
  integrate(function(z) k * z * dnorm(z) *
              exp((k - 1) * pnorm(z, log.p = TRUE)),
            lower = -Inf, upper = Inf, rel.tol = 1e-11,
            subdivisions = 2000L, stop.on.error = TRUE)
}

empty_moments <- function(p, labels = NULL) {
  list(n = 0, mean = setNames(numeric(p), labels),
       M2 = setNames(numeric(p), labels))
}

merge_moments <- function(acc, x) {
  # Merge batch means/centered sums of squares, avoiding subtraction of
  # almost equal accumulated raw moments.
  nb <- nrow(x)
  mb <- colMeans(x)
  centered <- sweep(x, 2, mb, "-")
  m2b <- colSums(centered * centered)
  nt <- acc$n + nb
  delta <- mb - acc$mean
  acc$M2 <- acc$M2 + m2b + delta^2 * acc$n * nb / nt
  acc$mean <- acc$mean + delta * nb / nt
  acc$n <- nt
  acc
}

make_theory <- function(opt, max_mean, model) {
  shrink <- if (model == "gaussian") sqrt(1 - opt$rho) else 1 - opt$rho
  c(signed_difference = 0,
    absolute_difference = 2 * opt$sd * shrink / sqrt(pi),
    squared_difference = 2 * opt$sd^2 * (1 - opt$rho),
    test_range = 2 * opt$sd * shrink * max_mean,
    highest_score = opt$mu + opt$sd * shrink * max_mean,
    lowest_score = opt$mu - opt$sd * shrink * max_mean)
}

simulate_model <- function(opt, max_mean, model) {
  theory <- make_theory(opt, max_mean, model)
  acc <- empty_moments(length(theory), names(theory))
  score_n <- 0
  score_mean <- numeric(opt$k)
  score_M2 <- matrix(0, opt$k, opt$k)
  chunks <- ceiling(opt$n / opt$batch)
  for (b in seq_len(chunks)) {
    nb <- min(opt$batch, opt$n - acc$n)
    e <- matrix(rnorm(nb * opt$k), nrow = nb, ncol = opt$k)
    common <- rnorm(nb)
    if (model == "gaussian") {
      x <- opt$mu + opt$sd * (sqrt(opt$rho) * common + sqrt(1 - opt$rho) * e)
    } else {
      x <- opt$mu + opt$sd * e
      shared <- runif(nb) < opt$rho
      if (any(shared)) {
        shared_score <- opt$mu + opt$sd * common[shared]
        x[shared, ] <- matrix(shared_score, nrow = sum(shared), ncol = opt$k)
      }
    }
    lo <- x[, 1]
    hi <- lo
    for (j in 2:opt$k) {
      lo <- pmin(lo, x[, j])
      hi <- pmax(hi, x[, j])
    }
    difference <- x[, 1] - x[, 2]
    statistics <- cbind(difference, abs(difference), difference^2, hi - lo, hi, lo)
    colnames(statistics) <- names(theory)
    acc <- merge_moments(acc, statistics)

    mb <- colMeans(x)
    centered <- sweep(x, 2, mb, "-")
    delta <- mb - score_mean
    nt <- score_n + nb
    score_M2 <- score_M2 + crossprod(centered) + tcrossprod(delta) * score_n * nb / nt
    score_mean <- score_mean + delta * nb / nt
    score_n <- nt
    if (b == 1 || b %% 10 == 0 || b == chunks) {
      message(sprintf("%s: %s / %s simulated people", model,
                      format(acc$n, big.mark = ",", scientific = FALSE),
                      format(opt$n, big.mark = ",", scientific = FALSE)))
    }
  }
  mc_sd <- sqrt(acc$M2 / (acc$n - 1))
  mc_se <- mc_sd / sqrt(acc$n)
  covariance <- score_M2 / (score_n - 1)
  table <- data.frame(
    quantity = names(theory), theory = unname(theory),
    simulation = unname(acc$mean), mc_se = unname(mc_se),
    mc_95_low = unname(acc$mean - qnorm(0.975) * mc_se),
    mc_95_high = unname(acc$mean + qnorm(0.975) * mc_se),
    error_in_mc_se = unname(ifelse(mc_se > 0, (acc$mean - theory) / mc_se, NA_real_))
  )
  list(table = table, means = score_mean, sds = sqrt(diag(covariance)),
       correlation = cov2cor(covariance), moments = acc)
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(opt$seed)

# Deterministic checks: closed forms, numerical quadrature, and covariance.
max_integral <- normal_max_mean(opt$k)
max_mean <- max_integral$value
max5_exact <- 5 / (4 * sqrt(pi)) * (1 + 6 / pi * asin(1 / 3))
max5_integral <- normal_max_mean(5)$value
max2_integral <- normal_max_mean(2)$value
stopifnot(abs(max5_integral - max5_exact) < 1e-9,
          abs(max2_integral - 1 / sqrt(pi)) < 1e-9,
          abs(2 * 15 * sqrt(0.4 / pi) - 10.7047446969166) < 1e-9)
target_cov <- opt$sd^2 * ((1 - opt$rho) * diag(opt$k) +
                         opt$rho * matrix(1, opt$k, opt$k))
factor_loadings <- opt$sd * cbind(rep(sqrt(opt$rho), opt$k),
                                sqrt(1 - opt$rho) * diag(opt$k))
stopifnot(max(abs(tcrossprod(factor_loadings) - target_cov)) < 1e-10 * opt$sd^2)

gaussian <- simulate_model(opt, max_mean, "gaussian")
mixture <- simulate_model(opt, max_mean, "mixture")

# Individual/model quantiles are different from Monte Carlo uncertainty.
probs <- c(0.025, 0.5, 0.975)
residual_sd <- opt$sd * sqrt(1 - opt$rho)
individual_quantiles <- data.frame(
  probability = probs,
  absolute_difference = sqrt(2) * residual_sd * qnorm((1 + probs) / 2),
  test_range = residual_sd * qtukey(probs, nmeans = opt$k, df = Inf)
)

check_message <- function(result) {
  # A 95% Monte Carlo interval need not cover its target on every run.
  # This is a diagnostic, not a rule to keep simulating until it passes.
  z <- result$table$error_in_mc_se
  if (any(abs(z[is.finite(z)]) > 6)) {
    "WARNING: at least one mean differs from theory by more than 6 MC SE; investigate."
  } else {
    "Monte Carlo check: every nondegenerate mean is within 6 MC SE of theory."
  }
}

report <- capture.output({
  cat("IQ TEST DIFFERENCES AND RANGES\n",
      "==============================\n\n", sep = "")
  cat(sprintf("People per model: %.0f; tests per person: %d; seed: %d; batch: %.0f\n",
              opt$n, opt$k, opt$seed, opt$batch))
  cat(sprintf("Marginal mean: %g; marginal SD: %g; every pair correlation: %g\n\n",
              opt$mu, opt$sd, opt$rho))
  cat("MODEL ASSUMPTIONS\n",
      "The main results assume an equicorrelated JOINT normal distribution.\n",
      "Each row is one person; independent rows define the averaging population.\n",
      "Scores are continuous and unrounded. Correlation is a population parameter,\n",
      "not the exact empirical correlation in the simulated finite sample.\n\n", sep = "")
  cat("DETERMINISTIC NUMERICAL CHECKS\n")
  cat(sprintf("E(max of %d independent standard normals): %.12f\n", opt$k, max_mean))
  cat(sprintf("Quadrature estimated absolute error: %.3g\n", max_integral$abs.error))
  cat(sprintf("E(max of 5), exact closed form: %.12f\n", max5_exact))
  cat(sprintf("Closed form minus quadrature: %.3g\n", max5_exact - max5_integral))
  cat("All deterministic checks passed.\n\n")
  cat("JOINT NORMAL MODEL\n")
  print(gaussian$table, row.names = FALSE, digits = 9)
  cat(check_message(gaussian), "\n\n")
  cat("Empirical test means:\n"); print(gaussian$means, digits = 8)
  cat("Empirical test SDs:\n"); print(gaussian$sds, digits = 8)
  cat("Empirical correlation matrix:\n"); print(gaussian$correlation, digits = 6)
  cat("\nINDIVIDUAL QUANTILES UNDER THE JOINT NORMAL MODEL\n")
  print(individual_quantiles, row.names = FALSE, digits = 8)
  cat("These are distribution quantiles, NOT confidence intervals for the mean.\n\n")
  cat("COUNTEREXAMPLE: NORMAL MARGINALS, SAME CORRELATIONS, NON-NORMAL JOINT LAW\n",
      "With probability rho all k scores equal one shared normal draw.\n",
      "Otherwise all k scores are independent normal draws.\n",
      "This is a non-identifiability demonstration, not a realistic test model.\n", sep = "")
  print(mixture$table, row.names = FALSE, digits = 9)
  cat(check_message(mixture), "\n\n")
  cat("Empirical counterexample test means:\n"); print(mixture$means, digits = 8)
  cat("Empirical counterexample test SDs:\n"); print(mixture$sds, digits = 8)
  cat("Empirical counterexample correlation matrix:\n"); print(mixture$correlation, digits = 6)
  cat("\nINTERPRETATION\n",
      "mc_se and mc_95_* quantify simulation error, not uncertainty in rho or\n",
      "the assumptions, and not the spread of an individual person's scores.\n",
      "squared_difference has units of IQ points squared; other means use IQ points.\n",
      "The 10 pairwise gaps within a row of 5 tests are not independent replicates.\n",
      "We use one specified pair and one range per independent person.\n",
      "The simulation supports implementation/numerics; the mathematical derivation\n",
      "proves the expectations under the stated model.\n\n", sep = "")
  cat("REPRODUCIBILITY\n"); print(sessionInfo())
})
writeLines(report, opt$out, useBytes = TRUE)
cat(paste(report, collapse = "\n"), "\n", sep = "")
message("Saved: ", normalizePath(opt$out, mustWork = TRUE))
