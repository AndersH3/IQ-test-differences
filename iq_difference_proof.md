# Expected differences between correlated IQ-test scores

## Result and assumptions

Let the scores of a randomly selected person on the same five tests be
\((X_1,\ldots,X_5)\). Assume:

1. The score vector is **jointly multivariate normal**.
2. Every test has the same population mean \(\mu\) and standard deviation \(\sigma\).
3. Every distinct pair has the same **population Pearson correlation** \(\rho=0.6\).
4. Scores are continuous and unrounded, with no ceiling/floor or practice effects in this model.

For the numerical answers, use \(\sigma=15\), and \(\mu=100\) when reporting absolute score levels. The common mean cancels from differences and ranges.

| Quantity | Exact model expectation / value | IQ points |
|---|---|---:|
| Signed difference, \(E(X_1-X_2)\) | \(0\) | 0 |
| Absolute difference, \(E\lvert X_1-X_2\rvert\) | \(2\sigma\sqrt{(1-\rho)/\pi}\) | 10.704744696917 |
| Root-mean-square difference, \(\sqrt{E[(X_1-X_2)^2]}\) | \(\sigma\sqrt{2(1-\rho)}\) | 13.416407864999 |
| Five-test range, \(E(\max_j X_j-\min_j X_j)\) | \(2\sigma\sqrt{1-\rho}\,a_5\) | 22.065699447377 |

Here

\[
a_5=\frac{5}{4\sqrt\pi}\left[1+\frac6\pi\arcsin\left(\frac13\right)\right]
=1.162964473640520\ldots.
\]

Thus the usual interpretation of “average difference” gives **10.70 IQ points for two tests**, and the expected maximum-minus-minimum over five gives **22.07 points**.

Correlation by itself does not imply these two answers. Section 5 gives a counterexample with the very same normal marginal distributions and correlation matrix.

## 1. A construction with exactly the required joint distribution

For \(0\leq\rho\leq1\), let \(Z_0,Z_1,\ldots,Z_k\) be mutually independent \(N(0,1)\) variables, and define

\[
X_j=\mu+\sigma\sqrt\rho\,Z_0+\sigma\sqrt{1-\rho}\,Z_j,
\qquad j=1,\ldots,k.
\]

Linear combinations of independent normal variables are jointly normal. Direct calculation gives

\[
E(X_j)=\mu,\qquad
\operatorname{Var}(X_j)=\sigma^2[\rho+(1-\rho)]=\sigma^2,
\]

and, for \(i\ne j\),

\[
\operatorname{Cov}(X_i,X_j)=\sigma^2\rho,
\qquad \operatorname{Corr}(X_i,X_j)=\rho.
\]

Consequently this construction has covariance matrix

\[
\Sigma=\sigma^2\big[(1-\rho)I_k+\rho\mathbf1\mathbf1^T\big].
\]

A multivariate normal distribution is determined by its mean vector and covariance matrix (for example, through its characteristic function). Therefore the construction is not just one possible Gaussian example: it represents **every** jointly normal score vector with the stated mean and covariance. See [Taboga's proof of closure under linear transformations](https://www.statlect.com/probability-distributions/normal-distribution-linear-combinations).

The common term is a mathematical representation of shared variation. Calling it a person's “true IQ,” and calling every remaining component measurement error, would require additional psychometric assumptions not used here.

## 2. Two-test difference: proof

Write \(D=X_1-X_2\). Equal means imply \(E(D)=0\), irrespective of normality. The covariance identity gives

\[
\begin{aligned}
\operatorname{Var}(D)
&=\operatorname{Var}(X_1)+\operatorname{Var}(X_2)
 -2\operatorname{Cov}(X_1,X_2)\\
&=2\sigma^2(1-\rho).
\end{aligned}
\]

Under joint normality, the linear contrast \(D\) is itself normal:

\[
D\sim N(0,\tau^2),\qquad \tau=\sigma\sqrt{2(1-\rho)}.
\]

Integrating its density,

\[
\begin{aligned}
E|D|
&=\frac{2}{\tau\sqrt{2\pi}}
  \int_0^\infty d\,e^{-d^2/(2\tau^2)}\,\mathrm d d\\
&=\frac{2\tau}{\sqrt{2\pi}}
=\tau\sqrt{\frac2\pi}
=2\sigma\sqrt{\frac{1-\rho}{\pi}}.
\end{aligned}
\]

The integral equals \(\tau^2\), as follows by the substitution \(u=d^2/(2\tau^2)\). Substituting \(\sigma=15\), \(\rho=0.6\) yields

\[
\boxed{E|X_1-X_2|=30\sqrt{0.4/\pi}=10.704744696917\ldots.}
\]

In contrast,

\[
\sqrt{E(D^2)}=\sqrt{180}=13.416407864999\ldots.
\]

The latter is the root-mean-square difference, not the mean absolute difference. With equal means and the specified variances/correlation, it is determined even without normality. Without joint normality, Cauchy–Schwarz gives \(E|D|\leq\sqrt{E(D^2)}\), but does not force the Gaussian value.

## 3. Range of five tests: proof

Set \(s=\sigma\sqrt{1-\rho}\), and let

\[
M_k=\max_{1\leq j\leq k} Z_j,
\qquad m_k=\min_{1\leq j\leq k} Z_j.
\]

Because the same quantity \(\mu+\sigma\sqrt\rho Z_0\) is added to every score, it cancels **pointwise**, before taking expectations:

\[
\begin{aligned}
R_k&=\max_j X_j-\min_j X_j\\
&=s\big(M_k-m_k\big).
\end{aligned}
\]

Joint sign symmetry of the independent standard normals implies \(m_k\overset d=-M_k\). These variables are integrable, since \(\max_j|Z_j|\leq\sum_j|Z_j|\) and \(E|Z_j|<\infty\). Hence

\[
E(R_k)=s\{E(M_k)-E(m_k)\}=2s\,E(M_k).
\]

Let \(\phi\) and \(\Phi\) denote the standard normal density and distribution function. Independence gives

\[
P(M_k\leq z)=\prod_{j=1}^kP(Z_j\leq z)=\Phi(z)^k.
\]

Differentiating yields \(f_{M_k}(z)=k\phi(z)\Phi(z)^{k-1}\), so

\[
\boxed{
E(R_k)=2\sigma\sqrt{1-\rho}\,k
\int_{-\infty}^{\infty}z\phi(z)\Phi(z)^{k-1}\,\mathrm dz.
}
\]

This is an exact expression, not a quantile approximation to an order statistic. For \(k=5\), the integral gives

\[
E(M_5)=1.162964473640520\ldots,
\]

and consequently

\[
\boxed{
E(R_5)=30\sqrt{0.4}\,(1.162964473640520\ldots)
=22.065699447377\ldots.
}
\]

The value of \(E(M_5)\) also agrees with the five-observation entry in the [University of Washington table of normal order-statistic expectations](https://faculty.washington.edu/heagerty/Books/Biostatistics/TABLES/NormalOrder/index.html). The derivation here does not depend on that table.

## 4. Exact closed form for five: proof

Using \(\phi'(z)=-z\phi(z)\), integration by parts gives

\[
\begin{aligned}
a_5:=E(M_5)
&=5\int_{-\infty}^{\infty}z\phi(z)\Phi(z)^4\,\mathrm dz\\
&=20\int_{-\infty}^{\infty}\phi(z)^2\Phi(z)^3\,\mathrm dz.
\end{aligned}
\]

The boundary term vanishes because \(0\leq\Phi(z)^4\leq1\) and \(\phi(z)\to0\) at either infinity.

For \(Y\sim N(0,1/2)\), its density is \(e^{-y^2}/\sqrt\pi\), and

\[
\phi(y)^2=\frac1{2\sqrt\pi}f_Y(y).
\]

Thus

\[
a_5=\frac{10}{\sqrt\pi}E\{\Phi(Y)^3\}.
\]

Introduce independent standard normals \(U_1,U_2,U_3\), independent of \(Y\). Conditioning on \(Y\) shows

\[
E\{\Phi(Y)^3\}=P(U_1\leq Y,U_2\leq Y,U_3\leq Y).
\]

The variables \(W_j=Y-U_j\) are centered jointly normal. Each has variance \(3/2\), and every pair has covariance \(1/2\); hence their correlations are \(1/3\).

For three centered jointly normal variables with correlations \(r_{12},r_{13},r_{23}\),

\[
P(W_1>0,W_2>0,W_3>0)
=\frac18+\frac1{4\pi}
  \{\arcsin r_{12}+\arcsin r_{13}+\arcsin r_{23}\}.
\]

For completeness, this formula follows from the following elementary argument:

* For a centered standardized normal pair of correlation \(r\), write it as \(A\) and \(rA+\sqrt{1-r^2}B\), with \(A,B\) independent standard normals. The direction of \((A,B)\) is uniform on the circle. The two positivity half-planes intersect in a wedge of angle \(\pi-\arccos r=\pi/2+\arcsin r\). Thus their joint positivity probability is \(1/4+\arcsin r/(2\pi)\).
* If \(S_j=\operatorname{sign}(W_j)\), then \(E(S_iS_j)=2\arcsin r_{ij}/\pi\), using the preceding bivariate result and central symmetry.
* Expand \(\prod_{j=1}^3(1+S_j)/8\), which is the indicator that all three variables are positive, outside a null set. Every odd sign-product expectation is zero by central symmetry. The remaining terms yield the formula above.

This orthant formula is also stated in [Pinasco, Smucler and Zalduendo (2020), introduction](https://arxiv.org/abs/2004.04682).

Putting \(r_{12}=r_{13}=r_{23}=1/3\) gives

\[
\begin{aligned}
a_5
&=\frac{10}{\sqrt\pi}
 \left[\frac18+\frac3{4\pi}\arcsin\left(\frac13\right)\right]\\
&=\boxed{\frac5{4\sqrt\pi}
 \left[1+\frac6\pi\arcsin\left(\frac13\right)\right]}.
\end{aligned}
\]

Combining this with Section 3 yields the fully closed-form range expectation

\[
\boxed{
E(R_5)=\frac{5\sigma\sqrt{1-\rho}}{2\sqrt\pi}
\left[1+\frac6\pi\arcsin\left(\frac13\right)\right].
}
\]

## 5. Why even normal marginals plus correlation 0.6 are insufficient

Let \(B\) be Bernoulli with \(P(B=1)=\rho\), independent of mutually independent standard normals \(Z_0,Z_1,\ldots,Z_k\). Define

\[
\widetilde X_j=\begin{cases}
\mu+\sigma Z_0,&B=1,\\
\mu+\sigma Z_j,&B=0.
\end{cases}
\]

Both conditional marginal distributions are \(N(\mu,\sigma^2)\), so every unconditional marginal is exactly \(N(\mu,\sigma^2)\). Also, for \(i\ne j\),

\[
E\big[(\widetilde X_i-\mu)(\widetilde X_j-\mu)\big]
=\rho\sigma^2+(1-\rho)\,0=\rho\sigma^2.
\]

Therefore the full correlation matrix is identical to that of the Gaussian construction in Section 1. But the joint law is not multivariate normal when \(0<\rho<1\): it puts positive probability on all scores being exactly equal.

When \(B=1\), both the absolute difference and the range are zero. When \(B=0\), the scores are independent. Consequently

\[
E|\widetilde X_1-\widetilde X_2|
=(1-\rho)\frac{2\sigma}{\sqrt\pi},
\]

and

\[
E\big(\max_j\widetilde X_j-\min_j\widetilde X_j\big)
=(1-\rho)\,2\sigma a_k.
\]

At \(\rho=0.6,\sigma=15,k=5\), these equal **6.770275002573** and **13.955573683686** IQ points, rather than 10.704745 and 22.065699.

This is a constructive proof of non-identifiability, not a claim that the counterexample is a realistic model of tests. In particular, “every score is normally distributed” is weaker than “the scores are jointly normal.”

## 6. Simulation design and interpretation

The companion file `iq_difference_simulation.R` uses only base R and defaults to **5,000,000 independent people per model**, each taking five tests. It implements both the Gaussian model and the counterexample. It:

* Uses the fixed seed 20260903 with an explicitly specified random-number generator.
* Simulates in batches of 100,000 people to keep memory use moderate.
* Numerically integrates the maximum-density formula and checks it against the exact five-test formula and the two-test special case.
* Checks the common-factor covariance matrix algebraically.
* Reports empirical test means, standard deviations, and all pair correlations.
* Reports simulation estimates, Monte Carlo standard errors, and approximate 95% Monte Carlo intervals for each expected value.
* Computes one specified pair's difference and one full range per person; it does not incorrectly treat all within-person pair differences as independent observations.
* Reports theoretical individual quantiles using `qnorm()` and `qtukey(..., df=Inf)`.
* Refuses to overwrite an existing results file.

Run:

```bash
Rscript iq_difference_simulation.R
```

Or choose a different run size and output filename:

```bash
Rscript iq_difference_simulation.R --n=1000000 --out=iq_results_1m.txt
```

### Actual verification run

The default script was executed in R 4.3.3 with 5,000,000 people per model. The saved output is `iq_simulation_results.txt`.

| Gaussian-model quantity | Theory | R simulation | Monte Carlo SE | Approximate 95% Monte Carlo interval |
|---|---:|---:|---:|---|
| Signed difference | 0 | 0.011730825 | 0.005998945 | -0.000026890 to 0.023488541 |
| Absolute difference | 10.704744697 | 10.701215470 | 0.003617203 | 10.694125882 to 10.708305057 |
| Five-test range | 22.065699447 | 22.063504833 | 0.003665385 | 22.056320810 to 22.070688857 |

Both requested nonnegative expectations were within one Monte Carlo SE of their theoretical values. The ten empirical off-diagonal correlations ranged from approximately 0.599650 to 0.600389. The exact expression for the standard-normal maximum and R's numerical integral both gave 1.162964473641; the integral's reported absolute-error estimate was approximately \(4.21\times10^{-12}\).

The counterexample run produced an absolute-difference mean of 6.767670308, compared with theory 6.770275003, and a five-test range mean of 13.947077098, compared with theory 13.955573684. Its Monte Carlo SEs were 0.005179889 and 0.008473844, respectively.

Additional small runs verified the two-test special case, correlation zero, correlation one (identical scores), a ten-test run and a final batch shorter than the requested batch size. An independent NumPy simulation and SciPy quadrature also supported the derivation; the tables above are from the actual R run, not the Python check.

Numerical integration uses [R's documented `integrate()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/integrate.html). Normal simulation uses [R's `rnorm()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Normal.html). The studentized range tends to the unstudentized standard-normal range as the degrees of freedom tend to infinity; see [R's range distribution documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Tukey.html).

The expectation formulas are proved analytically. Monte Carlo supports the numerical implementation; it is not the proof. The strong law of large numbers applies because the independent-person statistics have finite expectations, so their sample averages converge almost surely to these model expectations. Finite second moments also justify the usual asymptotic Monte Carlo standard-error intervals.

With SD 15, correlation 0.6 and five tests, the Gaussian model's individual range distribution has approximate 2.5th, 50th and 97.5th percentiles **8.06069, 21.41067 and 39.81649** points. This illustrates why an expected range of 22.07 does not mean each person's range will be close to 22.07. These population-distribution quantiles are not confidence intervals for the expected range.

## 7. Scope and limitations

* These are unconditional population expectations for a fixed set of tests under the model. They are not predictions conditional on an already observed first score or on selection for a high/low score.
* A fixed first score generally introduces regression-to-the-mean conditioning. For the bivariate Gaussian model, \(E(X_2\mid X_1=x)=\mu+\rho(x-\mu)\), so \(E(X_1-X_2\mid X_1=x)=(1-\rho)(x-\mu)\), not necessarily zero.
* The assumed correlation is the correlation of test scores across the modeled population, not automatically test reliability or the correlation within a selected ability group.
* All ten distinct pairs among five tests are assumed to have correlation 0.6. An average pairwise correlation of 0.6, or only the correlations with one reference test, is not enough to specify the five-test Gaussian model.
* Different score means, different standard deviations, non-normal joint dependence, ceilings, rounding, practice effects, or nonuniform pair correlations require a modified model.
* The answers scale linearly with the common SD: for another SD \(s\), multiply the SD-15 answers by \(s/15\).
* The common-factor representation in this document covers nonnegative correlations. Negative equicorrelations are outside the script's supported inputs.

## AI-use note

This analysis, proof exposition and R implementation were produced by ChatGPT in response to the user's question. Public mathematical and R references were checked; the numerical results are generated under explicit model assumptions, not estimated from observed IQ-test data. The accompanying run log identifies the R environment and actual simulation output.
