# Bayesian Zero-Inflated Generalized Poisson Model

This project studies a Bayesian zero-inflated generalized Poisson (ZIGP) regression model for count data with excess zeros and non-Poisson dispersion. The analysis implements the model estimation, model comparison, and case influence diagnostics in R using Markov Chain Monte Carlo methods.

The project applies the method to the `Couple` dataset, which is used to analyze unwanted pursuit behavior in the context of couple separation trajectories. The work reproduces and explains methodology from the referenced statistical literature where the original analysis code was not publicly available.

## Project Focus

- Model count outcomes with excess zeros using a ZIGP regression structure
- Estimate posterior distributions for regression and dispersion parameters with MCMC
- Compare models using Bayesian criteria such as DIC and CPO
- Evaluate case influence with Kullback-Leibler divergence diagnostics
- Present the statistical methodology and empirical findings in a reproducible R Markdown report

## Repository Files

- `ZIGP-project.Rmd` - Main report and analysis source
- `data.txt` - Data used by the report

The standalone `ZIGP-project.R` script was removed so the R Markdown report is the single source for the written analysis and code.

## How to Run

Open `ZIGP-project.Rmd` in RStudio and knit it to PDF. The report reads `data.txt` from the same repository folder:

```r
df = read.table("data.txt", header = TRUE)
```

Required R packages include:

```r
install.packages(c("mvtnorm", "HDInterval", "knitr", "rmarkdown"))
```

## Portfolio Summary

This project demonstrates custom Bayesian computation rather than use of a prebuilt ZIGP modeling function. The core technical work includes deriving the likelihood, implementing Metropolis updates, summarizing posterior estimates, comparing candidate models, and calculating influence diagnostics for individual observations.
