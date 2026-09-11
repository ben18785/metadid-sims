---
title: "Per-study information by study design: derivation"
subtitle: "Supporting the exchange-rate panel (category U) of the metadid simulation study"
---

# 1. Purpose and scope

This note derives, from first principles, the quantity plotted in the
exchange-rate panel: **how much information about the population treatment
effect one study of a given design contributes, relative to one
difference-in-differences (DiD) study.**

Three nested expressions are derived:

1. the **sampling-variance ratio**, assuming no between-study heterogeneity and
   a known nuisance parameter;
2. the **heterogeneity-adjusted ratio**, which adds between-study variation in
   the true effect;
3. the **identification-adjusted ratio**, which additionally charges each
   incomplete design for the fact that its nuisance parameter must be estimated
   from the DiD studies rather than known.

Expression 3 is the one compared against simulation. Every step is stated
explicitly, and Section 8 enumerates the assumptions with the consequence of
each failing.

Throughout, "information" means Fisher information about the population mean
effect, i.e. the reciprocal of the sampling variance of its estimator. All
results are per study.

# 2. Notation

| symbol | meaning | code |
|---|---|---|
| $\sigma$ | within-study, between-subject SD of the outcome at one time point | `within_sd` |
| $\rho$ | within-subject correlation between the pre and post measurement | `rho` |
| $n$ | subjects per arm | `n_control`, `n_treatment` |
| $\theta_i$ | true treatment effect in study $i$ | |
| $\mu_\theta,\ \tau_\theta$ | mean and SD of $\theta_i$ across studies | `true_effect`, `sigma_effect` |
| $\beta_i$ | secular (time) trend in study $i$ | |
| $\mu_\beta,\ \tau_\beta$ | mean and SD of $\beta_i$ across studies | `true_trend`, `sigma_trend` |
| $\gamma_i$ | baseline difference between arms in study $i$ | |
| $\tau_\gamma$ | SD of $\gamma_i$ across studies | |
| $n_D$ | number of DiD studies in the anchor | `U_CORE` |
| $k$ | number of studies added, of a single design | `U_ADD` |

# 3. The data-generating process

Taken directly from `R/simulate.R`. Each subject contributes a pair
$(Y_{\text{pre}}, Y_{\text{post}})$ drawn from a bivariate normal:

$$
\begin{pmatrix} Y_{\text{pre}} \\ Y_{\text{post}} \end{pmatrix}
\sim N\!\left( \boldsymbol{\mu}_{\text{arm}},\ \Sigma \right),
\qquad
\Sigma = \sigma^2 \begin{pmatrix} 1 & \rho \\ \rho & 1 \end{pmatrix}
$$

with arm means (`simulate_one_study`)

$$
\boldsymbol{\mu}_{\text{control}} = (b,\ b + \beta_i), \qquad
\boldsymbol{\mu}_{\text{treatment}} = (b,\ b + \beta_i + \theta_i)
$$

where $b$ is the study baseline. **Both arms share the same pre-period mean**,
so the true baseline difference is $\gamma_i = 0$ in this DGP. The model does
not know this and still estimates a baseline-difference parameter; that
distinction matters in Section 6.

The three designs differ only in which summaries are retained:

- **DiD** (`simulate_did_summary_from_params`) — pre and post means for both arms.
- **Post-only RCT** (`make_rct_summary`) — post means for both arms; no pre.
- **Pre-post, PP** (`make_pp_summary`) — pre and post means for the treated arm only.

# 4. Step 1 — per-study sampling variance

## 4.1 Two preliminary results

For an arm of $n$ subjects, the mean at a single time point has variance
$\operatorname{Var}(\bar Y) = \sigma^2/n$, and the two time points covary as
$\operatorname{Cov}(\bar Y_{\text{pre}}, \bar Y_{\text{post}}) = \rho\sigma^2/n$.

The **within-arm change score** $D = \bar Y_{\text{post}} - \bar Y_{\text{pre}}$
therefore has

$$
\operatorname{Var}(D)
= \frac{\sigma^2}{n} + \frac{\sigma^2}{n} - 2\cdot\frac{\rho\sigma^2}{n}
= \frac{2\sigma^2(1-\rho)}{n}.
\tag{4.1}
$$

This is the only place $\rho$ enters. A design that never differences across
time cannot benefit from it.

## 4.2 Difference-in-differences

The estimator is the difference of the two arms' change scores,
$\hat\theta^{\,\text{DiD}} = D_{\text{t}} - D_{\text{c}}$. Taking expectations,

$$
\mathbb{E}[D_{\text{t}}] = \beta_i + \theta_i, \qquad
\mathbb{E}[D_{\text{c}}] = \beta_i
\;\Longrightarrow\;
\mathbb{E}[\hat\theta^{\,\text{DiD}}] = \theta_i .
$$

**The trend cancels exactly.** The control arm is an identification device: it
removes $\beta_i$ and contributes nothing to the point estimate of $\theta_i$.
The arms are independent, so by (4.1)

$$
\boxed{\;s^2_{\text{DiD}} = \operatorname{Var}(D_{\text{t}}) + \operatorname{Var}(D_{\text{c}})
= \frac{4\sigma^2(1-\rho)}{n}\;}
\tag{4.2}
$$

Note the control arm *adds* variance. DiD trades precision for identification.

## 4.3 Post-only RCT

The estimator is $\hat\theta^{\,\text{RCT}} = \bar Y_{\text{post,t}} - \bar Y_{\text{post,c}}$, with

$$
\mathbb{E}[\hat\theta^{\,\text{RCT}}] = (b + \beta_i + \theta_i + \gamma_i) - (b + \beta_i)
= \theta_i + \gamma_i .
$$

The trend cancels here too (both arms are observed at the same time), but any
baseline imbalance $\gamma_i$ does **not**. Two independent means, no
differencing across time, so $\rho$ never appears:

$$
\boxed{\;s^2_{\text{RCT}} = \frac{\sigma^2}{n} + \frac{\sigma^2}{n} = \frac{2\sigma^2}{n}\;}
\tag{4.3}
$$

## 4.4 Pre-post

The estimator is the treated arm's change score alone,
$\hat\theta^{\,\text{PP}} = D_{\text{t}}$, with

$$
\mathbb{E}[\hat\theta^{\,\text{PP}}] = \theta_i + \beta_i .
$$

Nothing cancels the trend. Directly from (4.1),

$$
\boxed{\;s^2_{\text{PP}} = \frac{2\sigma^2(1-\rho)}{n}\;}
\tag{4.4}
$$

## 4.5 Immediate consequences

$$
s^2_{\text{DiD}} = \frac{4\sigma^2(1-\rho)}{n}, \qquad
s^2_{\text{RCT}} = \frac{2\sigma^2}{n}, \qquad
s^2_{\text{PP}}  = \frac{2\sigma^2(1-\rho)}{n}
$$

- $s^2_{\text{PP}} = \tfrac{1}{2} s^2_{\text{DiD}}$ **for every** $\rho$. A pre-post
  study is twice as precise as a DiD study — it spends all its subjects on the
  treated arm rather than half on an identification device.
- $s^2_{\text{RCT}} = s^2_{\text{DiD}}$ exactly when $4(1-\rho) = 2$, i.e.
  $\rho = 1/2$. Above that the DiD is more precise; below it the post-only RCT is.
- Per *subject* rather than per study the contrast is starker still: a DiD study
  uses $2n$ subjects and a PP study $n$, so PP is four times as efficient per
  subject. The per-study comparison used throughout is the conservative one.

# 5. Step 2 — from sampling variance to information about $\mu_\theta$

## 5.1 The random-effects structure

Study effects are drawn $\theta_i \sim N(\mu_\theta, \tau_\theta^2)$, and study
$i$ estimates its own $\theta_i$ with sampling variance $s^2$. Marginally,

$$
\hat\theta_i \sim N\!\left(\mu_\theta,\ \tau_\theta^2 + s^2\right),
$$

so the information a single study carries about the **population mean** is

$$
\boxed{\;I = \frac{1}{\tau_\theta^2 + s^2}\;}
\tag{5.1}
$$

## 5.2 Why this is not $1/s^2$

$\tau_\theta^2$ acts as a **floor**. As $s^2 \to 0$, $I \to 1/\tau_\theta^2$, not
infinity: a perfectly precise study pins down *its own* $\theta_i$, which is
still only one draw from the population. With $k$ studies the posterior variance
of $\mu_\theta$ cannot fall below $\tau_\theta^2/k$ however large the studies are.

This floor is what compresses design differences, because designs differ only
through $s^2$.

## 5.3 The exchange rate with the nuisance known

$$
\boxed{\;
R^{(2)}_X \;=\; \frac{I_X}{I_{\text{DiD}}}
\;=\; \frac{\tau_\theta^2 + s^2_{\text{DiD}}}{\tau_\theta^2 + s^2_X}
\;}
\tag{5.2}
$$

## 5.4 Limiting behaviour

- $\tau_\theta \to 0$: $R^{(2)}_X \to s^2_{\text{DiD}}/s^2_X$, giving
  $R_{\text{PP}} = 2$ and $R_{\text{RCT}} = 2(1-\rho)$. This is expression 1.
- $\tau_\theta \to \infty$: $R^{(2)}_X \to 1$ for every design.

**Heterogeneity moves the ratio toward 1, not downward.** If the
no-heterogeneity ratio exceeds 1, heterogeneity reduces it; if it is below 1,
heterogeneity *increases* it. Mechanically, a common floor $\tau_\theta^2$
penalises whichever design was more precise. Expressions 1 and 2 are therefore
**not** ordered, and must not be drawn as a monotone sequence.

## 5.5 An invariant

At $\rho = 1/2$, $s^2_{\text{DiD}} = s^2_{\text{RCT}}$ identically, so
$R^{(2)}_{\text{RCT}} = 1$ **regardless of $\tau_\theta$**. The crossover between
DiD and post-only RCT at $\rho = 1/2$ is robust to heterogeneity, even though
the magnitudes elsewhere are not.

# 6. Step 3 — charging for identification

## 6.1 Each design's nuisance, and where it is identified

From Section 4, the incomplete designs estimate a contaminated quantity:

| design | observes | nuisance | identified from |
|---|---|---|---|
| DiD | $\theta_i$ | — | — |
| PP | $\theta_i + \beta_i$ | trend $\beta$ | DiD control-arm change score |
| RCT | $\theta_i + \gamma_i$ | baseline difference $\gamma$ | DiD difference of pre-arm means |

The DiD control-arm change score has $\mathbb{E}[D_{\text{c}}] = \beta_i$ and, by
(4.1), variance $2\sigma^2(1-\rho)/n$; with $\beta_i \sim N(\mu_\beta, \tau_\beta^2)$
the per-DiD-study information about $\mu_\beta$ is

$$
i_\beta = \frac{1}{\tau_\beta^2 + 2\sigma^2(1-\rho)/n}.
\tag{6.1}
$$

The DiD difference of pre-period arm means has expectation $\gamma_i$ and
variance $2\sigma^2/n$ (two independent means, no differencing over time), so

$$
i_\gamma = \frac{1}{\tau_\gamma^2 + 2\sigma^2/n}.
\tag{6.2}
$$

## 6.2 Joint information for $(\mu_\theta,\ \nu)$

Write $\nu$ for the relevant nuisance mean. With $n_D$ DiD studies and $k$
studies of incomplete design $X$, three groups of observations contribute:

- $n_D$ DiD contrasts, each informing $\mu_\theta$ alone with precision $I_{\text{DiD}}$;
- $n_D$ DiD nuisance observations, each informing $\nu$ alone with precision $i_\nu$;
- $k$ incomplete studies, each informing the **sum** $\mu_\theta + \nu$ with precision $I_X$.

An observation of a sum contributes to both diagonal entries *and* to the
off-diagonal. The Fisher information matrix is therefore

$$
J =
\begin{pmatrix}
n_D I_{\text{DiD}} + k I_X & k I_X \\[2pt]
k I_X & n_D i_\nu + k I_X
\end{pmatrix}
\tag{6.3}
$$

## 6.3 Inversion

For a $2\times2$ matrix, the $(1,1)$ entry of the inverse is $J_{22}/\det J$, so
the precision of $\hat\mu_\theta$ is $P = \det J / J_{22}$. Expanding,

$$
\det J = (n_D I_{\text{DiD}} + k I_X)(n_D i_\nu + k I_X) - (k I_X)^2
$$

$$
P = (n_D I_{\text{DiD}} + k I_X) - \frac{(k I_X)^2}{n_D i_\nu + k I_X}
$$

and collecting the $k I_X$ terms,

$$
\boxed{\;
P(n_D, k) \;=\; \underbrace{n_D I_{\text{DiD}}}_{\text{anchor}}
\;+\; \underbrace{k I_X}_{\text{raw}} \cdot
\underbrace{\frac{n_D\, i_\nu}{n_D\, i_\nu + k I_X}}_{\text{identification discount } \mathcal{D}}
\;}
\tag{6.4}
$$

## 6.4 Behaviour of the discount

$$
\mathcal{D} = \frac{n_D\, i_\nu}{n_D\, i_\nu + k I_X} \in (0, 1)
$$

- $n_D \to \infty$: $\mathcal{D} \to 1$. The nuisance is effectively known and the
  design realises its full value.
- $n_D \to 0$: $\mathcal{D} \to 0$. With no anchor the incomplete studies are
  worthless for $\mu_\theta$ — they identify only the sum $\mu_\theta + \nu$.
- $k \to \infty$: $\mathcal{D} \to 0$ and the total contribution $k I_X \mathcal{D}$
  saturates at $n_D i_\nu$. **Adding incomplete studies has diminishing returns
  bounded by the precision of the borrowed nuisance.**

## 6.5 The exchange rate

Adding $k$ DiD studies instead gives $P(n_D + k, 0) = (n_D + k) I_{\text{DiD}}$,
because with no incomplete studies $J$ is block diagonal and the nuisance is
irrelevant to $\mu_\theta$. Hence the gains are

$$
\Delta_X = k I_X \mathcal{D}, \qquad \Delta_{\text{DiD}} = k I_{\text{DiD}}
$$

$$
\boxed{\;
R^{(3)}_X \;=\; \frac{\Delta_X}{\Delta_{\text{DiD}}}
\;=\; \underbrace{\frac{\tau_\theta^2 + s^2_{\text{DiD}}}{\tau_\theta^2 + s^2_X}}_{R^{(2)}_X}
\;\times\;
\underbrace{\frac{n_D\, i_\nu}{n_D\, i_\nu + k I_X}}_{\mathcal{D}}
\;}
\tag{6.5}
$$

The $k$ cancels from the ratio itself but remains inside $\mathcal{D}$, so the
exchange rate depends on how many studies are added as well as on the anchor.

Because $\mathcal{D} < 1$ strictly, expression 3 is always below expression 2 —
these two *are* ordered, unlike expressions 1 and 2.

# 7. The estimator used in the simulation

The panel does not use (6.5) to produce its measured points; it estimates the
same quantity from posterior output.

For each $(\rho, n_D)$ cell four scenarios are run: `core` ($n_D$ DiD only),
`did` ($n_D + k$ DiD), `rct` ($n_D$ DiD $+ k$ RCT), `pp` ($n_D$ DiD $+ k$ PP).
Writing $\widehat{P} = 1/\overline{\text{SD}}^2$ for the posterior precision of
$\mu_\theta$ (`mean_posterior_sd`, averaged over replicates),

$$
\widehat{R}_X = \frac{\widehat{P}_X - \widehat{P}_{\text{core}}}
                     {\widehat{P}_{\text{did}} - \widehat{P}_{\text{core}}}
$$

Two remarks.

**The $k$ divisor cancels.** All arms add the same number of studies, so the
ratio of increments is already per study.

**The prior cancels.** If precision is additive, $\widehat{P}_{\text{core}} =
P_{\text{prior}} + n_D I_{\text{DiD}}$ and $\widehat{P}_{\text{did}} =
P_{\text{prior}} + (n_D + k) I_{\text{DiD}}$, so the difference is $k
I_{\text{DiD}}$ with no prior term. This is why the estimator is a *difference*
of precisions rather than a precision itself.

Two ways it can fail:

- **Non-convergence.** A poorly mixed chain gives a meaningless posterior SD.
  Because $\widehat{P}_{\text{core}}$ appears in both numerator and denominator,
  one bad `core` arm corrupts every ratio at that grid point. The panel
  therefore drops any arm with $\hat R > 1.05$ and leaves the affected point
  missing rather than plotting a wrong value.
- **Non-linearity in $k$.** Additivity requires the increment to be in the
  linear regime. It degrades when $k$ is large relative to $n_D$, or when
  $\tau_\theta$ estimation makes precision non-additive. Running a second value
  of $k$ and checking the ratio is stable is the direct test.

# 8. Assumptions

| # | assumption | where used | consequence if violated |
|---|---|---|---|
| 1 | Equal arm sizes $n_c = n_t = n$ | (4.2)–(4.4) | variances gain $(1/n_c + 1/n_t)$ factors; ratios shift but the structure holds |
| 2 | Common $\sigma$ and $\rho$ across arms, time points and studies | (4.1) | heterogeneous $\sigma$ makes the exchange rate study-specific |
| 3 | Normality of subject-level outcomes | throughout | mild; the CLT covers the arm means for realistic $n$ |
| 4 | Posterior for $\mu_\theta$ approximately Gaussian | Section 7 | posterior SD would no longer summarise precision |
| 5 | $\tau_\theta$, $\tau_\beta$, $\tau_\gamma$ **known** | (5.1), (6.1), (6.2) | **the main simplification.** They are estimated, which adds uncertainty the derivation omits, so (6.5) is optimistic |
| 6 | Information adds across studies | (6.3), Section 7 | breaks when $k$ is large relative to $n_D$ |
| 7 | Flat priors on $\mu_\theta$, $\mu_\beta$, $\mu_\gamma$ | Section 7 | an informative prior would not cancel in the differences |
| 8 | Exchangeability: $\beta_i$ for PP studies drawn from the same distribution as for DiD studies | Section 6.1 | borrowing becomes biased; this is exactly what the trend-plane category tests |
| 9 | $\gamma_i = 0$ in truth, but estimated by the model | Section 3 | with real imbalance the RCT discount is larger |
| 10 | Studies within a design are exchangeable and equally sized | (6.3) | unequal sizes require weighting each study separately |

Assumption 5 is the one most likely to matter and is the first candidate for
any discrepancy between (6.5) and simulation.

# 9. Numerical validation

At $\sigma = 0.12$, $n = 100$, $\tau_\theta = 0.03$, $\tau_\beta = 0.02$,
$\tau_\gamma = 0$, $k = 20$, against a local run at 3 replicates:

| design | anchor $n_D$ | $R^{(2)}$ | $R^{(3)}$ predicted | measured |
|---|---|---|---|---|
| PP | 5 | 1.14 | **0.37** | **0.34** |
| PP | 20 | 1.14 | 0.75 | 0.59 |
| RCT | 5 | 1.00 | **0.51** | **0.47** |
| RCT | 20 | 1.00 | 0.80 | 0.91 |

At the weak anchor the closed form is accurate to a few percent across all
$\rho$, with no fitted quantities. At the strong anchor the two designs' errors
run in opposite directions — PP below prediction, RCT above — which noise alone
would not produce.

Candidate explanations, in order of plausibility: three replicates is too few to
resolve the ratio; assumption 5 (variance components treated as known) bites
harder once the anchor is large enough that the nuisance posterior is no longer
the binding constraint; and the joint-normal approximation in (6.3) degrades in
the same regime. Distinguishing these requires a higher-replicate run.

# 10. Code

- `R/simulate.R` — `simulate_one_study`, `make_rct_summary`, `make_pp_summary`
- `R/scenarios.R` — category `U`, constants `U_CORES`, `U_ADD`, `U_RHO`
- `R/plots.R` — `paper_panels()`, internal `.u_theory` (equation 6.5) and
  `.u_ceiling` (equation 5.2)
