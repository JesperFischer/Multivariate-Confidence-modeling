# Posterior Predictive Functions - Guide

This document explains the different posterior predictive functions and how they handle uncertainty when generating predictions.

## Three Sources of Uncertainty

When making predictions from hierarchical models, there are three distinct sources of uncertainty:

1. **Estimation Uncertainty**: Uncertainty about parameter values (captured by the posterior distribution)
2. **Response Variability**: Stochastic variability in individual responses from the data-generating distributions
3. **Subject-Level Variability**: Between-subject heterogeneity (how much individuals differ from each other)

---

## Overview of Functions

### Group-Level Functions

| Function | Estimation Uncertainty | Response Variability | Subject-Level Variability |
|----------|----------------------|---------------------|--------------------------|
| `Get_predictive_group()` | ✓ Across draws | ✗ Uses expected values | ✗ N/A (group-level) |
| `Get_predictive_group_responses()` | ✓ Across draws | ✓ Samples responses | ✗ N/A (group-level) |

### Subject-Level Functions

| Function | Estimation Uncertainty | Response Variability | Subject-Level Variability |
|----------|----------------------|---------------------|--------------------------|
| `Get_predictive()` | ✓ Across draws | ✓ Samples responses | ✓ Per-subject predictions |
| `Get_predictive_subject_average()` | ✗ Uses median | ✗ Uses expected values | ✓ Across subjects |
| `Get_predictive_subject_average_responses()` | ✗ Uses median | ✓ Samples responses | ✓ Across subjects |

---

## Detailed Function Descriptions

### 1. `Get_predictive_group()`
**Level 1: Group-level expected values**

**What it does:**
- Uses group mean parameters (`gm`) from each posterior draw
- Computes **expected values** (no sampling):
  - Binary: probability itself
  - RT: `exp(mu + sigma^2/2)` (lognormal expected value)
  - Confidence: `inv_logit(conf_mu)` (ordered beta mean)

**Uncertainty shown:** Pure estimation uncertainty about group mean parameters

**Best for:** Showing credibility intervals on the population mean response

**Interpretation:** "We're 90% confident the true population mean response is in this interval"

**Shrinks with more subjects?** Yes (fastest) - more subjects → better estimation of group means

---

### 2. `Get_predictive_group_responses()`
**Level 2a: Group-level sampled responses**

**What it does:**
- Uses group mean parameters (`gm`) from each posterior draw
- **Samples responses** from distributions:
  - Binary: `rbinom()`
  - RT: `rlnorm()`
  - Confidence: `rordbeta()`

**Uncertainty shown:** Estimation uncertainty + distributional variance

**Best for:** Showing what a hypothetical "average" subject's responses would look like

**Interpretation:** "A typical subject's response would fall in this interval 90% of the time"

**Shrinks with more subjects?** Yes (moderate) - parameter estimates improve but response variance remains

---

### 3. `Get_predictive()`
**Level 2b: Subject-level predictions with full posterior**

**What it does:**
- For each subject, for each posterior draw:
  - Uses that subject's parameters from that draw
  - Samples responses from distributions

**Uncertainty shown:** Subject-specific estimation uncertainty + response variability

**Best for:** Posterior predictive checks for individual subjects

**Interpretation:** "This is what we'd expect this specific subject's responses to look like"

**Shrinks with more subjects?** No** - each subject analyzed separately

---

### 4. `Get_predictive_subject_average()`
**Level 3: Between-subject variability (expected values)**

**What it does:**
- For each subject:
  - Uses posterior **median** parameters (point estimate)
  - Computes **expected values** (no sampling)

**Uncertainty shown:** Pure between-subject heterogeneity

**Best for:** Showing how much individuals differ in their expected behavior

**Interpretation:** "90% of individuals have their mean response in this range"

**Shrinks with more subjects?** No* - estimates population variability, which is stable

---

### 5. `Get_predictive_subject_average_responses()`
**Subject heterogeneity + response variability**

**What it does:**
- For each subject:
  - Uses posterior **median** parameters (point estimate)
  - **Samples responses** from distributions (one sample)

**Uncertainty shown:** Between-subject variability + response variability

**Best for:** Showing the full range of individual responses

**Interpretation:** "90% of individual observations fall in this range"

**Shrinks with more subjects?** No* - captures population variability + noise

---

## Using with `group_predictive()`

The `group_predictive()` function pools predictions across all dimensions by doing `group_by(X)` and computing quantiles.

### What Gets Pooled:

| Input Function | Dimensions in Output | What `group_predictive()` Pools | Credibility Intervals Represent |
|----------------|---------------------|-------------------------------|--------------------------------|
| `Get_predictive_group()` | `draw × X` | Across draws | Uncertainty about population mean |
| `Get_predictive_group_responses()` | `draw × X` | Across draws (with sampling) | Prediction interval for avg subject |
| `Get_predictive()` | `subject × draw × X` | Across subjects AND draws | Total uncertainty |
| `Get_predictive_subject_average()` | `subject × X` | Across subjects | Between-subject heterogeneity |
| `Get_predictive_subject_average_responses()` | `subject × X` | Across subjects (with sampling) | Subject heterogeneity + noise |

---

## Interpretation Summary

| Function | Mathematical Notation | CrI Width | Shrinks with N? | Best Description |
|----------|----------------------|-----------|----------------|------------------|
| `Get_predictive_group()` | E[Y \| X, μ_θ] | **Narrowest** | Yes (fastest) | "Uncertainty about the population mean" |
| `Get_predictive_group_responses()` | Y ~ f(X, μ_θ) | Moderate | Yes (moderate) | "Typical subject's response range" |
| `Get_predictive_subject_average()` | E[Y \| X, θ̄_s] | Wide | No | "How much individuals differ" |
| `Get_predictive()` | Y_s ~ f(X, θ_s) | **Widest** | No | "Full range of observations" |
| `Get_predictive_subject_average_responses()` | Y_s ~ f(X, θ̄_s) | Wide+ | No | "Subject differences + noise" |

\* Converges to true population SD, but doesn't shrink to zero  
\** Parameter estimates improve slightly with more data

---

## Matching to Data Analysis

When computing empirical summaries in `group_predictive()`, the data uses:

```r
q5 = mean(RT) - 2 * (sd(RT) / sqrt(n()))
```

This is the **Standard Error of the Mean (SEM)**, which matches:
- ✓ `Get_predictive_group()` - estimation uncertainty only

**Therefore, for comparing model predictions to group-level data, use `Get_predictive_group()`.**

---

## Quick Reference: Which Function Should I Use?

**For comparing to group-level empirical data:**
→ `Get_predictive_group()` - matches SEM approach

**For showing what an average subject does:**
→ `Get_predictive_group_responses()` - includes response variability

**For posterior predictive checks by subject:**
→ `Get_predictive()` - full subject-level predictions

**For showing individual differences:**
→ `Get_predictive_subject_average()` - between-subject variability

**For showing the full range of individual observations:**
→ `Get_predictive_subject_average_responses()` - subjects + noise

---

## Suggested Plot Titles and Captions

### For `Get_predictive_group()`
**Mathematical notation:** E[Y | X, μ_θ]

Where μ_θ represents the group mean parameters.

**Title:** "Population mean predictions: E[Y | X, μ_θ]"

**Caption:** 
> Lines show the posterior mean of the expected response value; shaded regions show [95, 90, 80]% credibility intervals reflecting uncertainty in the group-level parameters μ_θ. For binary: P(Y=1|X); RT: E[log-normal]; Confidence: E[ordered-beta].

---

### For `Get_predictive_group_responses()`
**Mathematical notation:** Y ~ f(X, μ_θ)

Where μ_θ represents the group mean parameters.

**Title:** "Population response predictions: Y ~ f(X, μ_θ)"

**Caption:**
> Lines show the posterior mean of sampled responses; shaded regions show [95, 90, 80]% prediction intervals reflecting both uncertainty in group-level parameters μ_θ and distributional variance. Responses sampled from: binary ~ Bernoulli(μ_θ(X)); RT ~ log-normal(μ(X), σ); Confidence ~ ordered-beta(μ(X), φ).

---

### For `Get_predictive_subject_average()`
**Mathematical notation:** E[Y | X, θ̄_s]

Where θ̄_s represents subject-specific posterior median parameters (bar indicates point estimate).

**Title:** "Between-subject variability: E[Y | X, θ̄_s]"

**Caption:**
> Lines show the median across subjects; shaded regions show [95, 90, 80]% intervals capturing between-subject heterogeneity. Each subject's expected response computed at their posterior median parameters θ̄_s. Shows how much individuals differ in their mean behavior.

---

### For `Get_predictive()` (and `Get_predictive_subject_average_responses()`)
**Mathematical notation:** 
- `Get_predictive()`: Quantiles of {Y_s ~ f(X, θ_s)}
- `Get_predictive_subject_average_responses()`: Quantiles of {Y_s ~ f(X, θ̄_s)}

Where θ_s represents subject-specific parameters (full posterior for `Get_predictive()`) and θ̄_s represents posterior median parameters (for `Get_predictive_subject_average_responses()`). When plotted with `group_predictive()`, the credibility intervals show quantiles across the subject distribution.

**Title (for `Get_predictive_subject_average_responses()`):** "Individual response variability: Quantiles of Y_s ~ f(X, θ̄_s)"

**Caption:**
> Lines show the median across subjects; shaded regions show [95, 90, 80]% intervals capturing between-subject heterogeneity and within-subject response variability. Each subject's response sampled from their distribution at posterior median parameters θ̄_s. Intervals show where 50%, 80%, 90%, 95% of subjects' responses fall.

---

### R Code Examples

```r
# 1. Group-level expected values
group_predictive(pred_pure_group, df) +
  labs(
    title = expression(paste("Population mean predictions: E[Y | X, ", mu[theta], "]")),
    subtitle = "Expected value of responses given group-level parameters",
    caption = "Lines show the posterior mean of the expected response value; shaded regions show [95, 90, 80]% 
credibility intervals reflecting uncertainty in the group-level parameters μ_θ."
  )

# 2. Group-level sampled responses
group_predictive(pred_pure_group_resp, df) +
  labs(
    title = expression(paste("Population response predictions: Y ~ f(X, ", mu[theta], ")")),
    subtitle = "Sampled responses from group-level parameters",
    caption = "Lines show the posterior mean of sampled responses; shaded regions show [95, 90, 80]% 
prediction intervals reflecting both uncertainty in group-level parameters μ_θ and distributional variance."
  )

# 3. Subject-level expected values (between-subject variability)
group_predictive(pred_pure_group_subject, df) +
  labs(
    title = expression(paste("Between-subject variability: Quantiles of E[Y | X, ", bar(theta)[s], "]")),
    subtitle = "Distribution of expected responses across individuals",
    caption = "Lines show the median across subjects; shaded regions show [95, 90, 80]% intervals 
capturing between-subject heterogeneity. Each subject's expected response computed at posterior median θ̄_s."
  )

# 4. Subject-level sampled responses  
group_predictive(pred_pure_group_subject_resp, df) +
  labs(
    title = expression(paste("Individual response variability: Quantiles of {", Y[s], " ~ f(X, ", bar(theta)[s], ")}")),
    subtitle = "Distribution of sampled responses across individuals",
    caption = "Lines show the median across subjects; shaded regions show [95, 90, 80]% intervals 
capturing between-subject heterogeneity and within-subject response variability. 
Intervals show where most subjects' responses fall."
  )
```

---

## Conditional vs Marginal Effects Perspective

Understanding these functions in terms of **conditional** and **marginal** effects helps clarify what you're marginalizing over:

### Conditional Effects
**Conditioning on specific parameter values** - what happens for a particular set of parameters.

### Marginal Effects  
**Marginalizing (averaging) over uncertainty** - integrating out sources of variability.

---

### Function Breakdown by Conditional/Marginal:

| Function | Conditional On | Marginalized Over | Effect Type |
|----------|----------------|-------------------|-------------|
| `Get_predictive_group()` | Group parameters μ_θ | Posterior draws of μ_θ | **Marginal effect of X on E[Y]** at population level |
| `Get_predictive_group_responses()` | Group parameters μ_θ | Posterior draws of μ_θ + response noise | **Marginal predictive distribution** for typical subject |
| `Get_predictive_subject_average()` | Subject parameters θ̄_s (point estimate) | Subjects | **Marginal effect across population** (no param uncertainty) |
| `Get_predictive_subject_average_responses()` | Subject parameters θ̄_s (point estimate) | Subjects + response noise | **Marginal predictive distribution** across subjects |
| `Get_predictive()` | Subject-specific parameters θ_s | Posterior draws + subjects + response noise | **Fully marginal prediction** (all uncertainty) |

---

### Detailed Explanations:

#### **`Get_predictive_group()` - Conditional on population, marginal over estimation**
```
E_p(μ_θ|data)[E[Y|X, μ_θ]]
```
- **Conditional**: Y given X and population parameters μ_θ
- **Marginal**: Over the posterior distribution of μ_θ
- **Interpretation**: "The average effect of X on Y at the population level, accounting for uncertainty about what the population parameters are"
- **Statistical term**: Marginal population-level expectation

---

#### **`Get_predictive_group_responses()` - Conditional on population, marginal over estimation + noise**
```
E_p(μ_θ|data)[Y ~ f(X, μ_θ)]
```
- **Conditional**: On population parameters μ_θ
- **Marginal**: Over posterior of μ_θ AND response distribution
- **Interpretation**: "What responses would we see from a typical subject, accounting for both parameter uncertainty and natural response variability?"
- **Statistical term**: Posterior predictive distribution at population level

---

#### **`Get_predictive_subject_average()` - Conditional on point estimates, marginal over subjects**
```
Quantiles_s[E[Y|X, θ̄_s]]
```
- **Conditional**: On each subject's median parameters θ̄_s (point estimates)
- **Marginal**: Over the distribution of subjects
- **Interpretation**: "How does the effect of X on Y vary across individuals in the population?"
- **Statistical term**: Marginal distribution of subject-level conditional expectations
- **Note**: This shows **population heterogeneity** in the effect of X, not parameter uncertainty

---

#### **`Get_predictive_subject_average_responses()` - Conditional on point estimates, marginal over subjects + noise**
```
Quantiles_s[Y_s ~ f(X, θ̄_s)]
```
- **Conditional**: On each subject's median parameters θ̄_s
- **Marginal**: Over subjects AND response noise
- **Interpretation**: "What is the distribution of individual observations across the population?"
- **Statistical term**: Marginal distribution of individual-level predictions

---

#### **`Get_predictive()` - Fully marginal (when pooled with group_predictive)**
```
E_p(θ_s|data),s[Y_s ~ f(X, θ_s)]
```
- **Conditional**: Only on X
- **Marginal**: Over posterior of θ_s, over subjects, over response noise
- **Interpretation**: "What is our total uncertainty about a response at stimulus level X, considering everything we don't know?"
- **Statistical term**: Fully marginal posterior predictive distribution

---

### Key Insight: Levels of Marginalization

Think of it as progressively averaging out uncertainty:

1. **Most Conditional** (`Get_predictive_group()`): 
   - Condition on: Population structure
   - Marginalize: Only parameter estimation uncertainty
   - Shows: "The" population effect with its uncertainty

2. **Add Response Noise** (`Get_predictive_group_responses()`):
   - Marginalize: Parameter uncertainty + distributional noise
   - Shows: What one typical observation looks like

3. **Add Subject Heterogeneity** (`Get_predictive_subject_average()`):
   - Marginalize: Subject-level variation (no parameter uncertainty)
   - Shows: How individuals differ systematically

4. **Add Everything** (`Get_predictive_subject_average_responses()` or `Get_predictive()`):
   - Marginalize: All sources of uncertainty
   - Shows: Total predictive distribution

---

### In Terms of Mixed Models Terminology:

| Function | Fixed Effects | Random Effects | Residual Error |
|----------|---------------|----------------|----------------|
| `Get_predictive_group()` | ✓ Marginal | N/A | ✗ |
| `Get_predictive_group_responses()` | ✓ Marginal | N/A | ✓ Marginal |
| `Get_predictive_subject_average()` | ✓ Conditional* | ✓ Marginal | ✗ |
| `Get_predictive_subject_average_responses()` | ✓ Conditional* | ✓ Marginal | ✓ Marginal |
| `Get_predictive()` | ✓ Marginal | ✓ Marginal | ✓ Marginal |

\* Uses point estimates (median) so not truly marginal over fixed effects uncertainty

---

### Practical Example:

**Research Question**: "How does stimulus intensity (X) affect reaction time?"

- **`Get_predictive_group()`**: "On average, across the population, a 1-unit increase in X changes RT by β ± SE"
  - Population-level (marginal) effect

- **`Get_predictive_subject_average()`**: "Individual A's RT changes by β_A per unit X, Individual B changes by β_B, etc."
  - Subject-specific (conditional) effects, showing variability

- **`Get_predictive_subject_average_responses()`**: "If we sample one RT from each person, this is the distribution we'd see"
  - Marginal distribution of observations

- **`Get_predictive()`**: "Accounting for all uncertainty, this is the full range of possible RTs we might observe"
  - Fully marginal prediction

---

## Notes

- All functions support both `model = "pure"` and `model = "full"` 
- The "full" model includes metacognitive uncertainty (`meta_un`) and RT uncertainty (`rt_un`) parameters
- Constants (`rt_ndt`, `c0`, `c11`) are averaged across subjects in group-level functions
- Subject-level functions use each subject's own constant values
