# Compound Poisson–Gamma Simulation Shiny Application

This repository contains an R Shiny application that simulates the compound Poisson process defined as:

S(t) = X₁ + X₂ + ... + Xₙ(t)

where:

- N(t) is a Poisson process with rate λ
- Xᵢ are i.i.d. exponential random variables with rate μ
- N(t) and Xᵢ are independent

The application enables interactive experimentation with parameter values, visual analysis of distributions, numerical summarization, and downloadable outputs.

---

## Mathematical Derivation

### 1. Distribution of N(t)
If interarrival times are exponentially distributed with rate λ, the resulting counting process follows:

N(t) ~ Poisson(λt)

Thus:

P(N(t)=n) = exp(-λt) · (λt)ⁿ / n!

### 2. Conditional Distribution of S(t)
Given N(t)=n events occur up to time t:

S(t) = X₁ + X₂ + ... + Xn

Since each Xᵢ ~ Exponential(μ), their sum is Gamma distributed:

S(t) | N(t)=n ~ Gamma(shape = n, rate = μ)

Its density is:

f(s|n) = μⁿ sⁿ⁻¹ exp(-μs) / (n-1)!   for s>0

### 3. Unconditional Distribution of S(t)
Since N(t) is random, S(t) is a **Poisson mixture of Gamma distributions**.

P(S(t)=0) = P(N(t)=0) = exp(-λt)

For s>0:

f(s) = Σ (exp(-λt) (λt)ⁿ/n!) [μⁿ sⁿ⁻¹ exp(-μs)/(n-1)!], summed over n≥1

Closed form exists but numerical simulation provides more intuitive insight.

### 4. Moments

E[N(t)] = λt  
E[S(t)] = (λt)/μ  
Var(S(t)) = (2λt)/μ²  

These theoretical values are shown in the application alongside empirical simulation results.

---

## Features of the Shiny Application

- Adjustable input parameters:
  - Poisson arrival rate λ
  - Exponential rate μ
  - Simulation count
  - Time horizon t
  - Optional random seed for reproducibility
- Drag-based collapsible control panel for clean visualization layout
- Histogram with density curve overlay
- Summary statistics including:
  - Sample mean, variance, quantiles
  - Theoretical E[N(t)] and E[S(t)]
- Dynamic interpretive insights based on parameter magnitude
- Export options:
  - Simulation data as CSV
  - Plot as PNG

---

## Insights from Simulation

| Parameter Change | Observed Effect |
|------------------|-----------------|
| Increasing λ      | More arrival events → larger cumulative value of S(t) |
| Increasing μ      | Faster decay of exponential jumps → decreases S(t) |
| Increasing t      | Both mean and spread of S(t) increase |
| Small λt          | Right-skewed, often with mass near zero |
| Large λt          | Distribution becomes smoother, closer to normality (CLT behaviour) |

Key takeaway:  
S(t) behaves like the total accumulated jump size over time. Larger λ or t increases accumulation frequency; μ determines individual jump size.

---


---

## Installation Requirements

Install required packages:

```r
install.packages(c("shiny", "shinyjqui", "shinyWidgets", "ggplot2"))

