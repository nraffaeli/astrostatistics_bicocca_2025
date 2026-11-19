# Statistical Inference

These notes summarize the main concepts from Lectures L05–L09 of
[*Astrostatistics Bicocca 2025*](https://github.com/dgerosa/astrostatistics_bicocca_2025).
The organization is by topic rather than lecture, so that frequentist and
Bayesian approaches can be compared directly.

## 1. Statistical inference

### 1.1 What is statistical inference?

Statistical inference is the process of drawing conclusions about a population
or a data-generating process from a sample of observations.

- A **population** is the distribution of possible measurements or objects of
  interest. It can also be the distribution of repeated measurements of one
  object.
- A **sample** is the finite set of observations that is actually collected.
- A **statistic** is any function of the sample, such as the sample mean.
- An **estimator** is a rule that uses the data to estimate a population
  quantity or model parameter.
- A **point estimator** returns one value.
- An **interval estimator** returns a range of plausible values.

Typical inference questions are:

1. What is the best estimate of a parameter? -> **point estimator**
2. How uncertain is that estimate? -> **interval estimator**
3. Are the data compatible with a specified model or hypothesis? -> **hypotesis testing**

Measurements have uncertainties, not simply errors. A good analysis models
the uncertainty and propagates it to the final inference.

### 1.2 Frequentist and Bayesian paradigms

The two paradigms use probability differently:

| Aspect | Frequentist inference | Bayesian inference |
|---|---|---|
| Meaning of probability | Long-run relative frequency over repeated experiments | Degree of belief or uncertainty |
| Parameters | Fixed but unknown quantities | Unknown quantities represented by probability distributions |
| Data | Treated as random before observation and fixed after observation | Observed data update uncertainty about parameters |
| Main result | Estimator and its sampling distribution | Posterior probability distribution |
| Uncertainty interval | Confidence interval/region | Credible interval/region |
| Model comparison | Goodness-of-fit, likelihood methods, and tests | Posterior odds, evidence, and Bayes factors |

### 1.3 General properties of an inference method

Useful properties for evaluating an estimator or inference procedure include:

- **Bias:** whether its expected value differs from the quantity being
  estimated.
- **Precision:** how concentrated repeated estimates are.
- **Efficiency:** how much information is extracted relative to competing
  estimators.
- **Consistency:** whether the estimator approaches the true value as the
  sample size increases.
- **Accuracy:** how close an estimate is to the target, combining systematic
  and random effects.
- **Robustness:** whether conclusions remain reasonably reliable when model
  assumptions are only approximately true.

Robustness matters when:

- outliers are present;
- the data are not exactly Gaussian;
- measurement uncertainties are misspecified; or
- the assumed model is an approximation to reality.

Robustness is a property of a method, estimator, likelihood, or model, not a
third statistical paradigm. Both frequentist and Bayesian analyses can be
robust or non-robust.

Some common estimators are:

| Estimator | Definition or principle | Main properties |
|---|---|---|
| **Sample mean** | $\displaystyle \bar{x}=\frac{1}{N}\sum_{i=1}^{N}x_i$ | Efficient for Gaussian data, but sensitive to outliers |
| **Weighted mean** | $\displaystyle \bar{x}_w=\frac{\sum_i w_i x_i}{\sum_i w_i}$, usually $w_i=1/\sigma_i^2$ | Uses measurements with different uncertainties; optimal for independent Gaussian errors when the weights are correct |
| **Sample median** | Middle-ordered value, or the 50th percentile | More robust to outliers than the mean, but generally less efficient for Gaussian data |
| **Sample variance** | $\displaystyle s^2=\frac{1}{N-1}\sum_{i=1}^{N}(x_i-\bar{x})^2$ | Unbiased for the population variance under independent identical sampling |
| **Maximum-likelihood estimator (MLE)** | $\displaystyle \hat{\theta}_{\mathrm{ML}}=\arg\max_\theta L(\theta)$ | Uses the assumed likelihood; often consistent and asymptotically efficient |
| **Method-of-moments estimator** | Equates sample moments to theoretical moments | Simple and often useful for initialization, but not generally efficient |

No estimator is uniformly best: the choice depends on the data, the assumed
model, the uncertainty structure, and the relative importance of robustness
and efficiency.

**NOTE:** for computing many of those we start from sample !!!

### 1.4 Why model the sample?

A model compresses the main structure of the data into a form we can interpret,
compare, and use for prediction. It also makes assumptions explicit and allows
uncertainties to be propagated to scientific conclusions.

Understanding data means formulating hypotheses and/or building models that explain them.

## 2. Non-parametric methods

### 2.1 Definition

Non-parametric methods model data without imposing a specific finite-parameter
family. They are more data-driven and flexible, but usually require more data
and careful smoothing or regularization.

Non-parametric does not mean assumption-free.

### 2.2 Histograms

A histogram estimates a distribution by dividing the range of the data into
bins and counting observations in each bin. It is useful for visualizing
location, spread, skewness, and multiple modes.

Its appearance depends strongly on bin width and origin: wide bins can hide
structure, while narrow bins can create noisy or apparently multimodal
structure. Empty or low-count bins may produce misleading gaps. The bin centre
is only a representation of a bin, not the measured value, and the resulting
step-like density is an approximation. A histogram is therefore a useful first
view, not a unique estimate of the distribution.

#### Equal-occupancy bins

Variable-width bins can be chosen to contain approximately the same number of
observations. They reduce nearly empty bins, but the result still depends on
the chosen occupancy; bins become wider in sparse regions and narrower in
dense regions.

#### Rules of thumb for the bin width

There is no universally correct bin width. Two useful rules are:

| | | |
|---|---|---|
|**Scott's rule** | $\Delta_b = \frac{3.5\sigma}{N^{1/3}}$ | It works best for approximately Gaussian data and is sensitive to outliers.|
|**Freedman–Diaconis rule** | $ \Delta_b = \frac{2(q_{75}-q_{25})}{N^{1/3}}$ | It is more robust for non-Gaussian or heavy-tailed data.|
-----

These rules are starting points, not prescriptions. Inspecting more than one
reasonable bin width is good practice, and the final bin must include the
maximum data value.

#### Normalized heights and bin errors

For bin $k$ with width $\Delta_{b,k}$, count $n_k$, and total sample size $N$,
the normalized density height is

$$
f_k = \frac{n_k}{\Delta_{b,k}N}.
$$

The approximate Poisson counting error is

$$
\sigma_k = \frac{\sqrt{n_k}}{\Delta_{b,k}N}.
$$

For complete coverage of the sample,

$$
\sum_k f_k\Delta_{b,k}=1
$$

and low-count bins have large relative uncertainties. For empty or very
low-count bins, exact Poisson limits are preferable to the simple
$\sqrt{n_k}$ approximation.

### 2.3 Kernel Density Estimation

Kernel Density Estimation (KDE) replaces each observation with a smooth kernel:

$$
\hat f_h(x)
=\frac{1}{Nh}\sum_{i=1}^{N}
K\left(\frac{x-x_i}{h}\right),
$$

where $K$ is the kernel and $h$ is the bandwidth.

- Small $h$: detailed but noisy estimate, with possible overfitting.
- Large $h$: smooth estimate, with possible loss of structure.

The bandwidth is therefore a bias–variance trade-off. KDE builds a smooth
estimate of the distribution directly from the sample, which can then be used
for inference. Unlike a single summary statistic such as the mean or median,
it retains information about the distribution's shape.

### 2.5 Frequentist and Bayesian perspective

The lectures present histograms, KDE, and the KS test mainly as frequentist or
descriptive tools. Their Bayesian counterparts would introduce a prior over
possible distributions or density functions and infer a posterior over those
distributions. The general comparison remains:

- frequentist methods describe sampling behavior of an estimator or test;
- Bayesian methods describe posterior uncertainty over the distribution itself.

## 3. Parametric models

### 3.1 Definition

A **parametric model** describes the data using a specified family of
probability distributions with a finite-dimensional parameter vector
$\boldsymbol{\theta}$:

$$
p(x \mid \boldsymbol{\theta}, M).
$$

The model M specifies assumptions such as:

- the functional relationship between variables;
- the probability distribution of residuals;
- the measurement uncertainties;
- independence or dependence between observations; and
- the parameters that are known or unknown.

Parametric inference is powerful because it is compact and interpretable, but
its conclusions depend on the adequacy of the assumed model. Comparing two different model is therefore crucial.

### 3.2 Likelihood

For observed data $D=\{x_i\}$, the **likelihood** is

$$
L(\boldsymbol{\theta}) = p(D \mid \boldsymbol{\theta},M).
$$

It is the probability of the data viewed as a function of the parameters.
The likelihood is not a probability distribution for the parameters:
parameters are the argument of the function, not the random outcome.

If observations are conditionally independent,

$$
L(\boldsymbol{\theta})
= \prod_{i=1}^{N} p(x_i\mid\boldsymbol{\theta},M),
$$

and it is usually more convenient to work with the log-likelihood:

$$
\ln L(\boldsymbol{\theta})
= \sum_{i=1}^{N}\ln p(x_i\mid\boldsymbol{\theta},M).
$$

Products become sums, and numerical underflow is avoided.

### 3.3 Fitting a parametric model: frequentist approach

The **maximum-likelihood estimator (MLE)** is

$$
\hat{\boldsymbol{\theta}}_{\mathrm{ML}}
= \underset{\boldsymbol{\theta}}{\operatorname{arg\,max}}\,
L(\boldsymbol{\theta}),
$$

or equivalently the value that maximizes $\ln L$.

---

For Gaussian residuals,

$$
y_i = M(x_i;\boldsymbol{\theta}) + r_i,
\qquad r_i\sim\mathcal{N}(0,\sigma_i),
$$

the log-likelihood is, up to a constant,

$$
\ln L(\boldsymbol{\theta})
= -\frac{1}{2}\sum_i
\left[
\frac{(y_i-M(x_i;\boldsymbol{\theta}))^2}{\sigma_i^2}
\right]
$$

when the $\sigma_i$ values are known. Therefore maximizing the likelihood is
equivalent to minimizing a weighted sum of squared residuals, or
$\chi^2$:

$$
\chi^2(\boldsymbol{\theta})
= \sum_i
\frac{(y_i-M(x_i;\boldsymbol{\theta}))^2}{\sigma_i^2}.
$$

For a straight-line model,

$$
M(x_i;a,b)=a x_i+b,
$$

the fitted parameters are the slope $a$ and intercept $b$. The same
likelihood principle applies to nonlinear models.

---

Important properties often associated with MLEs, particularly in large samples,
are consistency, approximate normality, and efficiency under suitable regularity
conditions. These properties are not automatic when the model is wrong or the
sample is small.

#### Large-sample properties of MLEs

If the model is correct, identifiable, and sufficiently regular, the MLE
$\hat{\boldsymbol{\theta}}_{\mathrm{ML}}$ has three key large-sample properties:

1. **Consistency**

   It converges in probability to the true parameter
   $\boldsymbol{\theta}_0$:

   $$
   \hat{\boldsymbol{\theta}}_{\mathrm{ML}}
   \xrightarrow[N\to\infty]{P}
   \boldsymbol{\theta}_0.
   $$

2. **Asymptotic normality**

   For large $N$, its distribution is approximately Gaussian:

   $$
   \hat{\boldsymbol{\theta}}_{\mathrm{ML}}
   \approx
   \mathcal{N}\left(
   \boldsymbol{\theta}_0,
   F^{-1}
   \right).
   $$

3. **Asymptotic efficiency**

   Its covariance approaches the Cramér–Rao lower bound:

   $$
   \operatorname{Cov}(\hat{\boldsymbol{\theta}}_{\mathrm{ML}})
   \simeq
   F^{-1}.
   $$

   Thus, under these assumptions, the MLE is asymptotically efficient. These
   results need not hold for small samples or a misspecified model.

#### Uncertainty around an MLE

The MLE is only a point estimate. Its uncertainty can be obtained from the
sampling distribution of the estimator or, in regular large-sample problems,
from the local curvature of the log-likelihood around its maximum. A sharply
peaked likelihood indicates that the data constrain the parameter strongly; a
wide or nearly flat likelihood indicates weak information.

The Gaussian approximation is useful when the likelihood is approximately
quadratic near its maximum. It can be misleading for small samples, parameter
boundaries, strong parameter correlations, or non-Gaussian likelihoods. In
those cases the full likelihood or a numerical confidence construction should
be retained rather than summarized by a symmetric error bar.

#### Non-Gaussian likelihoods

The Gaussian likelihood is convenient but is not mandatory. A different
likelihood can describe asymmetric errors, heavy tails, or a different
measurement process. Changing the likelihood changes both the fitted
parameters and their uncertainty. This is one of the places where the
robustness definition from Chapter 1 is important: a model with a heavy-tailed
likelihood can reduce the influence of observations far from the main trend.

The Student's $t$ distribution is a useful heavy-tailed alternative to a
Gaussian likelihood. Its tails allow a larger probability for observations
far from the centre. As its degrees of freedom increase, it approaches a
Gaussian distribution. Thus, the choice of likelihood controls a continuous
trade-off between sensitivity to deviations and robustness to extreme values.

### 3.4 Fisher matrix and parameter errors

The Fisher matrix is the negative curvature of the log-likelihood at the MLE:

$$
F_{jk}
=-\left.
\frac{\partial^2\ln L}
{\partial\theta_j\,\partial\theta_k}
\right|_{\boldsymbol{\theta}=\hat{\boldsymbol{\theta}}_{\mathrm{ML}}}.
$$

For independent data, information from separate observations adds. The inverse
Fisher matrix approximates the parameter covariance:

$$
\operatorname{Cov}(\hat{\boldsymbol{\theta}}_{\mathrm{ML}})
\simeq F^{-1}.
$$

The marginal standard error of parameter $\theta_j$ is

$$
\sigma_j=\sqrt{(F^{-1})_{jj}},
$$

while the off-diagonal entries $(F^{-1})_{jk}$ describe parameter
covariances and correlations. The intrinsic spread $\sigma$ is not
automatically the uncertainty on a fitted parameter such as $\mu$.


### 3.5 Fitting a parametric model: Bayesian approach

Bayesian inference combines the likelihood with a prior distribution:

$$
p(\boldsymbol{\theta}\mid D,M)
=
\frac{p(D\mid\boldsymbol{\theta},M)\,
p(\boldsymbol{\theta}\mid M)}
{p(D\mid M)}.
$$

The terms are:

- $p(\boldsymbol{\theta}\mid D,M)$: posterior distribution;
- $p(D\mid\boldsymbol{\theta},M)$: likelihood;
- $p(\boldsymbol{\theta}\mid M)$: prior distribution;
- $p(D\mid M)$: evidence or marginal likelihood.

The evidence is the normalization:

$$
p(D\mid M)
=\int p(D\mid\boldsymbol{\theta},M)
p(\boldsymbol{\theta}\mid M)\,d\boldsymbol{\theta}.
$$

The posterior can be explored on a grid or with numerical sampling methods.
The **maximum a posteriori (MAP)** estimate is

$$
\hat{\boldsymbol{\theta}}_{\mathrm{MAP}}
= \underset{\boldsymbol{\theta}}{\operatorname{arg\,max}}\,
p(\boldsymbol{\theta}\mid D,M).
$$

With a uniform prior, maximizing the posterior is equivalent to maximizing the
likelihood. In general, however, the prior changes both the location and the
uncertainty of the result.

Bayesian inference does not require reducing the result to one best-fit value:
the full posterior can be used to calculate parameter estimates, uncertainties,
correlations, predictions, and probabilities of scientifically relevant
quantities. Working with the full posterior can be computationally difficult,
especially in high dimensions, but modern numerical methods make this possible.

#### Choosing priors

The prior should represent information available before using the data. The
lectures discuss three guiding ideas:

- **Indifference:** assign equal weight when the parameterization makes equal
  alternatives genuinely comparable.
- **Invariance or consistency:** conclusions should not depend arbitrarily on
  how the same parameter is represented or transformed.
- **Maximum entropy:** when only limited information is known, choose the
  distribution that adds the fewest unsupported assumptions while satisfying
  the known constraints.

Priors are not automatically harmless. A prior can be influential when the
data are weak, when the parameter space is bounded, or when model comparison
depends on the volume of prior parameter space. Prior sensitivity should
therefore be checked when scientific conclusions depend strongly on it.

#### Conjugate and hierarchical models

A **conjugate prior** is selected so that the posterior belongs to the same
distributional family as the prior. Conjugacy can provide analytical
posteriors and useful intuition, although it is a mathematical convenience
rather than a requirement of Bayesian inference.

In **hierarchical Bayes**, parameters themselves are described by higher-level
parameters. This allows information to be shared among related observations or
objects while retaining uncertainty in the population-level quantities. The
result is a joint posterior over all levels of the model.

### 3.6 Model comparison

Model fitting finds the best parameters within one model. **Model comparison**
asks which of several models is adequate or preferred. A good fit does not
prove that a model is true: *garbage in, garbage out*.

#### Frequentist comparison

Cross-validation offers a predictive alternative by comparing performance on
data not used for fitting.

For a fitted model, evaluate the likelihood or, for Gaussian errors, the
chi-squared statistic at the best-fit parameters:

$$
\chi^2_{\min}
=\sum_i
\frac{(y_i-M(x_i;\hat{\boldsymbol{\theta}}))^2}{\sigma_i^2}.
$$

With $N$ data points and $k$ fitted parameters,

$$
\nu=N-k,
\qquad
\chi^2_{\mathrm{red}}=\frac{\chi^2_{\min}}{\nu}.
$$

When the model and uncertainties are appropriate, $\chi^2_{\mathrm{red}}$ is
typically of order one. It is a diagnostic, not proof of model validity.

When models have different numbers of parameters, penalized criteria provide
useful heuristics:

$$
\mathrm{AIC}=2k-2\ln L_{\max},
\qquad
\mathrm{BIC}=k\ln N-2\ln L_{\max}.
$$

Lower values are preferred. AIC and BIC balance fit quality against model
complexity; they are not probabilities that a model is true.

#### Bayesian comparison

Bayesian comparison uses the evidence:

$$
p(D\mid M)
=\int p(D\mid\boldsymbol{\theta},M)
p(\boldsymbol{\theta}\mid M)\,d\boldsymbol{\theta}.
$$

For two models, the **Bayes factor** is

$$
B_{12}
=\frac{p(D\mid M_1)}{p(D\mid M_2)}.
$$

Including prior model probabilities gives posterior odds:

$$
\frac{p(M_1\mid D)}{p(M_2\mid D)}
=B_{12}\,
\frac{p(M_1)}{p(M_2)}.
$$

The evidence averages the likelihood over the prior parameter space. Therefore,
a flexible model is penalized when it fits well only in a small part of a broad
parameter range: this is the Bayesian Occam effect. Bayes factors and posterior
odds depend on the specified models and their priors.

#### Frequentist and Bayesian comparison

| Question | Frequentist approach | Bayesian approach |
|---|---|---|
| Fit quality | $\chi^2$, likelihood, goodness-of-fit | Evidence and posterior |
| Complexity | AIC/BIC or a selected test | Integrated into the evidence |
| Model preference | Lower criterion or better test result | Bayes factor and posterior odds |
| Prior information | Not usually explicit | Explicitly included |

## 4. Uncertainty quantification

### 4.1 Frequentist uncertainty

In the frequentist framework, an estimator has a sampling distribution: the
distribution of estimates obtained from repeated samples under the same
data-generating process.

A **confidence interval** is constructed so that, over repeated experiments,
the procedure contains the fixed true parameter at a specified frequency.
For example, a 95% confidence procedure covers the true value in 95% of
repeated applications under its assumptions.

This does not mean that, after one interval is calculated, there is a 95%
probability that the fixed parameter lies inside that particular interval.

### 4.2 Bayesian uncertainty

A Bayesian posterior is a probability distribution for the unknown parameter
conditional on the observed data and model $p(\theta\mid D)$. A **credible interval** or **credible region** contains a specified amount of
posterior probability. For example:

$$
\int_{\mathcal{C}_{0.95}}p(\theta\mid D)\,d\theta=0.95.
$$

Credible regions can be central or highest-posterior-density regions and can be
asymmetric when the posterior is non-Gaussian.

### 4.3 Nuisance parameters and marginalization

A nuisance parameter is needed by the model but is not the main scientific
quantity of interest. Bayesian inference handles it by marginalization:

$$
p(\theta_1\mid D)
=\int p(\theta_1,\theta_2\mid D)\,d\theta_2.
$$

This retains uncertainty in the nuisance parameter instead of replacing it
with an unjustified fixed value.

### 4.4 Direct comparison

| Feature | Confidence interval | Bayesian credible interval |
|---|---|---|
| Probability statement | About repeated applications of a procedure | About posterior probability for the parameter |
| Requires | Sampling distribution and assumptions about repeated data | Likelihood, prior, and observed data |
| Parameter status | Fixed unknown value | Quantity described probabilistically |
| Shape | Often based on approximations or test inversion | Directly follows the posterior and may be asymmetric |

The two intervals can have similar numerical values in regular, large-sample
problems, but they have different interpretations.

## 5. Hypothesis testing

### 5.1 Frequentist hypothesis testing

Hypothesis testing asks whether data are compatible with a **null hypothesis**
$H_0$. A test statistic $S$ measures a feature of the data, and its null
distribution is used to calculate the **p-value**:

$$
p
=P(S\geq S_{\mathrm{obs}}\mid H_0)
$$

for an upper-tail test. The tail definition depends on the statistic and the
alternative hypothesis.

For one point $x_i$, the upper-tail probability under a null density $h_0(x)$
is

$$
p_i
=\int_{x_i}^{\infty}h_0(x)\,dx
=1-H_0(x_i),
$$

where $H_0$ is the null CDF.

#### Rejection of the null hypothesis

Choose a significance level $\alpha$ before examining the result:

$$
p\leq\alpha
\quad\Longrightarrow\quad
\text{reject }H_0.
$$

Otherwise, *do not reject* $H_0$; this is not the same as proving that $H_0$
is true. The two possible decision errors are:

- **False positive (Type I):** reject a true $H_0$.
- **False negative (Type II):** do not reject a false $H_0$.

#### Parametric and non-parametric tests

The same logic applies to both. Parametric tests specify a finite-parameter
null distribution; non-parametric tests make fewer distributional assumptions.
The KS test is an example of the latter.

#### KS test

The Kolmogorov–Smirnov (KS) test compares cumulative distribution functions. For a dataset, the KS statistic is the maximum CDF discrepancy:

$$
D_N=\sup_x\left|F_N(x)-F_0(x)\right|.
$$

For two samples:

$$
D_{N,M}=\sup_x\left|F_N(x)-G_M(x)\right|.
$$

Under a continuous, fully specified $H_0$ and large $N$,
$\lambda=\sqrt{N}D_N$ follows the Kolmogorov distribution:

$$
P(\sqrt{N}D_N\leq\lambda)
=1-2\sum_{k=1}^{\infty}(-1)^{k-1}
\exp(-2k^2\lambda^2).
$$

For two samples, the corresponding large-sample scaling uses

$$
N_{\mathrm{eff}}=\frac{NM}{N+M},
\qquad
\lambda=\sqrt{N_{\mathrm{eff}}}\,D_{N,M}.
$$

The exact finite-sample null distribution depends on the sample sizes. The
usual KS calibration also assumes continuous data; ties or estimated null
parameters require an adjusted calibration or simulation.

### 5.2 Bayesian hypothesis testing

Bayesian testing treats hypotheses as competing models. For models $M_1$ and
$M_2$, comparison is based on posterior odds:

$$
\frac{p(M_1\mid D)}{p(M_2\mid D)}
=
\frac{p(D\mid M_1)}{p(D\mid M_2)}
\frac{p(M_1)}{p(M_2)}.
$$

The ratio of evidences is the Bayes factor. It measures how the data update
the relative support for the two models.

Unlike a p-value, a Bayes factor directly compares specified alternatives. It
still depends on the model definitions and, through the evidence and prior
model odds, on prior choices.

Bayesian hypothesis testing is therefore comparative: it evaluates the
specified models against one another. A Bayes factor is not an absolute
probability that a model is true unless it is combined with prior model
probabilities, and the result can change if relevant alternatives are omitted.

### 5.3 Frequentist and Bayesian testing

| Question | Frequentist answer | Bayesian answer |
|---|---|---|
| How surprising are the data under $H_0$? | p-value | Likelihood/posterior predictive probability |
| Should a null be rejected? | Compare p-value with $\alpha$ | Compare posterior model odds or Bayes factor |
| Can hypotheses receive probabilities? | Not in the standard framework | Yes, through posterior model probabilities |

The two approaches answer related but different questions and should not be
interpreted as interchangeable terminology.

## 6. Formula sheet

### Common probability distributions

**Uniform distribution**

$$
U(x\mid a,b)=
\begin{cases}
\dfrac{1}{b-a}, & a\le x\le b,\\
0, & \text{otherwise},
\end{cases}
\qquad
\mathbb{E}[X]=\frac{a+b}{2},
\qquad
\operatorname{Var}(X)=\frac{(b-a)^2}{12}.
$$

**Gaussian (normal) distribution**

$$
\mathcal{N}(x\mid\mu,\sigma)
=\frac{1}{\sqrt{2\pi}\sigma}
\exp\!\left[-\frac{(x-\mu)^2}{2\sigma^2}\right],
\qquad
\mathbb{E}[X]=\mu,
\qquad
\operatorname{Var}(X)=\sigma^2.
$$

**Chi-squared distribution** with $\nu$ degrees of freedom:

$$
\chi^2(x\mid\nu)
=\frac{x^{\nu/2-1}e^{-x/2}}
{2^{\nu/2}\Gamma(\nu/2)},
\qquad x\ge 0,
\qquad
\mathbb{E}[X]=\nu,
\qquad
\operatorname{Var}(X)=2\nu.
$$

For large $\nu$:

$$
\frac{\chi^2_\nu-\nu}{\sqrt{2\nu}}
\xrightarrow[\nu\to\infty]{d}\mathcal{N}(0,1),
\qquad
\frac{\chi^2_\nu}{\nu}
\approx\mathcal{N}\!\left(1,\frac{2}{\nu}\right).
$$

**Poisson distribution** for a count $k=0,1,2,\ldots$ with expected rate
$\mu$:

$$
P(K=k\mid\mu)
=\frac{\mu^k e^{-\mu}}{k!},
\qquad
\mathbb{E}[K]=\mu,
\qquad
\operatorname{Var}(K)=\mu.
$$

### Likelihood

$$
L(\theta)=p(D\mid\theta,M).
$$

### Bayes’ theorem

$$
p(\theta\mid D,M)
=\frac{p(D\mid\theta,M)p(\theta\mid M)}
{p(D\mid M)}.
$$

### Evidence

$$
p(D\mid M)=\int p(D\mid\theta,M)p(\theta\mid M)\,d\theta.
$$

### Bayes factor

$$
B_{12}=\frac{p(D\mid M_1)}{p(D\mid M_2)}.
$$

### Posterior odds

$$
\frac{p(M_1\mid D)}{p(M_2\mid D)}
=B_{12}\frac{p(M_1)}{p(M_2)}.
$$

## 7. Final comparison

| Topic | Frequentist inference | Bayesian inference |
|---|---|---|
| Probability | Long-run frequency | Degree of belief |
| Parameter | Fixed unknown value | Probability-distributed unknown |
| Parametric fitting | MLE, least squares, likelihood | Posterior inference, MAP, marginalization |
| Main uncertainty object | Sampling distribution | Posterior distribution |
| Interval | Confidence interval | Credible interval |
| Model comparison | Goodness of fit, likelihood, tests | Evidence, Bayes factor, posterior odds |
| Hypothesis testing | p-values and significance levels | Comparison of posterior model probabilities |
| Prior information | Usually external to the formal procedure | Explicitly encoded in priors |
| Robustness | Can be improved with suitable estimators and models | Can be improved with suitable likelihoods and priors |
| Main risk | Misinterpreting p-values or relying on asymptotic assumptions | Sensitivity to priors or poorly specified models |

The central lesson is that neither paradigm removes the need for scientific
judgment. Both require a defensible data model, an appropriate treatment of
uncertainty, and careful interpretation of the resulting numbers.

## Sources

- [L05 — Classical/Frequentist Statistical Inference I](</Users/ffara/Desktop/LML/astrostatistics_bicocca_2025/lectures/L05_frequentist.ipynb>)
- [L06 — Classical/Frequentist Statistical Inference II](</Users/ffara/Desktop/LML/astrostatistics_bicocca_2025/lectures/L06_frequentist.ipynb>)
- [L07 — Classical/Frequentist Statistical Inference III](</Users/ffara/Desktop/LML/astrostatistics_bicocca_2025/lectures/L07_frequentist.ipynb>)
- [L08 — Bayesian Statistical Inference I](</Users/ffara/Desktop/LML/astrostatistics_bicocca_2025/lectures/L08_bayesian.ipynb>)
- [L09 — Bayesian Statistical Inference II](</Users/ffara/Desktop/LML/astrostatistics_bicocca_2025/lectures/L09_bayesian.ipynb>)
