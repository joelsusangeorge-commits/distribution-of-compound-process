# Compound Poisson Process Simulation Using R Shiny

This repository contains an R Shiny application that simulates and visualizes a Compound Poisson–Exponential process of the form:

S(t) = X₁ + X₂ + ... + Xₙ , where N(t) denotes the number of events up to time t.

The model assumptions are:

- N(t) follows a Poisson process with rate λ
- Xᵢ are independent exponential random variables with rate μ
- N(t) and Xᵢ are independent processes

The purpose of the application is to study the distribution of S(t), generate simulated realizations, visualize histograms, and explore the sensitivity of the model with respect to λ, μ, and time. The Shiny interface allows live experimentation and comparison across time scales.

---

## Mathematical Foundation

### Distribution of Event Count N(t)
If interarrival times are exponential with rate λ, then:

N(t) follows a Poisson distribution with mean λt.

### Conditional Distribution of S(t)
Given N(t) = n:

S(t) represents the sum of n exponential(μ) random variables, which implies it follows a Gamma distribution with:

- shape parameter = n  
- rate parameter = μ

### Unconditional Behavior of S(t)

- P(S(t) = 0) = exp(-λt)
- For S(t) > 0, its density is a mixture of Gamma densities weighted by Poisson probabilities

### Mean and Variance

E[S(t)] = (λt) / μ  
Var(S(t)) = (2λt) / μ²

These results provide theoretical support for simulation outputs and help interpret plots generated in the application.

---

## Project Goals

- Analyse the compound Poisson–Exponential system analytically
- Simulate realizations of S(t) and compare with theory
- Visualize histograms for different values of t
- Construct an interactive visualization tool using Shiny
- Observe how λ and μ influence the distribution of S(t)
  
Histograms are generated and compared for:

- t = 10  
- t = 100  
- t = 1000  
- t = 10000

In addition, users can generate histograms at custom time points.

---

## Application Features

- Histogram visualization of simulated aggregate values S(t)
- Predefined time plots + custom user-defined time input
- Real-time simulation with adjustable parameters:
  - Rate of arrivals (λ)
  - Rate of jumps (μ)
  - Total simulation count
  - Time horizon
- Sample path visualization for selected parameter sets
- Parameter sensitivity analysis through interactive UI controls
- Lightweight, responsive layout suitable for web deployment

---

## How to Run

### Install Required Package

```r
install.packages("shiny")
