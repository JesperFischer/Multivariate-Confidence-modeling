# Posterior Predictive Functions - Guide

This document explains the posterior predictive functions for generating and visualizing model predictions.

## Two Main Types of Predictions

When making predictions from hierarchical models, we can show:

1. **Credibility Intervals on the Mean**: Where the population average response lies
2. **Population Marginal Intervals**: Where individual subjects' average responses would fall

Both account for estimation uncertainty, but the second also includes between-subject heterogeneity.

---

## Primary Functions (Recommended for Papers)

### Group-Level Functions

| Function | Estimation Uncertainty | Subject Heterogeneity | What It Shows |
|----------|----------------------|---------------------|---------------|
| `Get_predictive_group()` | ✓ Across draws | ✗ N/A (population mean) | **Credibility interval on population mean** |
| `Get_predictive_new_subject()` | ✓ Across draws | ✓ Samples new subjects | **Population marginal predictions** |

---

## Detailed Function Descriptions

### 1. `Get_predictive_group()` ⭐ **Primary Function**
**Credibility intervals on the population mean**

**What it does:**
- Uses group mean parameters (`gm`) from each posterior draw
- Computes **expected values** (no sampling):
  - Binary: probability itself
  - RT: `exp(mu + sigma^2/2)` (lognormal expected value)
  - Confidence: `inv_logit(conf_mu)` (ordered beta mean)

**Uncertainty shown:** Pure estimation uncertainty about group mean parameters μ_θ

**Mathematical notation:** E[Y | X, μ_θ]

**Best for:** Main results - showing the population-level effect with credibility intervals

**Interpretation:** "We're 90% confident the true population mean response is in this interval"

**Matches standard practice:** ✓ Yes - equivalent to `marginaleffects::predictions()`, `emmeans::emmeans()`

**Shrinks with more subjects?** Yes (fastest) - more subjects → better estimation of group means

---

### 2. `Get_predictive_new_subject()` ⭐ **For Showing  Between subject varability **
**Population marginal predictions (integrating over random effects)**

**What it does:**
- For each posterior draw:
  - Extracts group means μ_θ, Between subject variances, and correlation matrix L
  - Samples n_subjects new hypothetical subjects from MVN(μ_θ, Σ), where Σ is the variance covariance matrix
  - Computes expected values for each hypothetical subject
- Pools across hypothetical subjects and posterior draws

**Uncertainty shown:** 
- Parameter uncertainty: How well do we know μ_θ and Σ?
- Population heterogeneity: How variable are individuals in the population?

**Mathematical notation:** E_{θ_s ~ MVN(μ_θ, Σ)}[E[Y | X, θ_s]]

**Best for:** Showing "where would a new subject from this population fall?"

**Interpretation:** "90% of new subjects from this population would have their mean response in this interval"

**Shrinks with more subjects?** 
- More real subjects in study → narrower intervals (better estimation of μ_θ, Σ)
- More hypothetical subjects in function (n_subjects parameter) → stabilizes around 100-200

**Note:** Use `n_subjects = 100` - sufficient for stable Monte Carlo approximation

---

## Comparison Summary

| Function | Mathematical Notation | Interval Width | What It Shows |
|----------|----------------------|----------------|---------------|
| `Get_predictive_group()` | E[Y \| X, μ_θ] | **Narrower** | "Where is the population mean?" |
| `Get_predictive_new_subject()` | E_s[E[Y \| X, θ_s]] | **Wider** | "Where do individuals fall on average?" |

**Key insight:** The marginal function adds between-subject variability on top of parameter uncertainty

**Variance decomposition:**
```
Var[marginal] = Var(μ_θ | data) + E[Var(θ_s | μ_θ, Σ)]
                    ↑                      ↑
            estimation uncertainty    population heterogeneity
```

---

## Matching to the Emperical Data.

When computing the plots from `group_predictive()`, the data uses:

```r
q5 = mean(RT) - 2 * (sd(RT) / sqrt(n()))
```

This is the **Standard Error of the Mean (SEM)**, which matches:
- ✓ `Get_predictive_group()` - estimation uncertainty only

**Therefore, for comparing model predictions to group-level data, use `Get_predictive_group()`.**


---

### Detailed Explanations:

#### **`Get_predictive_group()` - Marginal over estimation uncertainty**
```
E_{p(μ_θ|data)}[E[Y | X, μ_θ]]
```
- **Marginal over**: The posterior distribution of μ_θ
- **Interpretation**: "The population-level effect of X on Y, accounting for uncertainty about what the population parameters are"
- **Statistical term**: Marginal population-level expectation
- **Equivalent to**: Fixed effects predictions in mixed models

---

#### **`Get_predictive_new_subject()` - Marginal over estimation + population heterogeneity**
```
∫∫ E[Y | X, θ_s] p(θ_s | μ_θ, Σ) p(μ_θ, Σ | data) dθ_s d(μ_θ,Σ)
```
- **Marginal over**: Posterior of (μ_θ, Σ) AND the distribution of subject parameters θ_s
- **Interpretation**: "What responses would we expect from new individuals sampled from this population?"
- **Statistical term**: Marginal population prediction
- **Equivalent to**: Predictions with random effects integrated out

---

