# Volatility Prediction of 2 ETFs: VT and EEM

## Project Overview

This project involves detailed financial time series analysis using various GARCH models to study the volatility and correlation of financial instruments' returns, specifically focusing on two stocks: VT and EEM. The analysis spans model fitting, forecasting, and correlation analysis, using a robust suite of statistical tools and libraries in R.

## Prerequisites

Before running the analysis, ensure the following libraries are installed in your R environment:

- `readxl`
- `fUnitRoots`
- `astsa`
- `tidyverse`
- `fImport`
- `fBasics`
- `rugarch`
- `highfrequency`
- `xts`
- `zoo`
- `stats`
- `chron`
- `forecast`
- `quantmod`
- `aTSA`
- `sandwich`
- `dplyr`
- `tidyr`
- `PortfolioAnalytics`
- `DEoptim`
- `ROI`
- `ROI.plugin.quadprog`
- `ROI.plugin.glpk`

These can be installed using R's `install.packages()` function.

## Data Acquisition

Historical daily adjusted prices for VT and EEM are fetched from 2008 and 2003 respectively up to March 2023 using the `quantmod` package:

```R
getSymbols('VT', from='2008-06-26', to='2023-03-31', auto.assign=TRUE, periodicity='daily')
getSymbols('EEM', from='2003-04-14', to='2023-03-31', auto.assign=TRUE, periodicity='daily')
```

## Analysis Tasks

### 1. **Volatility and Returns Analysis**
   - Graphical representation of daily prices, log prices, and log returns.
   - Calculation and graphical analysis of autocorrelation functions to identify any significant patterns or outliers.

### 2. **GARCH Model Fitting**
   - Fit GARCH(1,1) models to the yields of VT and EEM using a normal distribution. Assess the fit and check for the presence of heteroscedasticity, skewness, and the appropriateness of the distribution.

### 3. **Extended GARCH Modeling**
   - Extend the analysis by fitting EGARCH and APARCH models, employing both normal and skewed-t distributions. Evaluate improvements in fit using statistical tests and information criteria.

### 4. **Forecasting**
   - Conduct one-step-ahead forecasts for a year using all model specifications. Compare forecasts using the Diebold-Mariano test to rank the models.

### 5. **Correlation Analysis**
   - Analyze the correlations between the standardized residuals of the two financial instruments over different window sizes. Assess whether these correlations remain constant over time through graphical analysis.

## Findings Summary

### VT
- The APARCH model with a skewed-t distribution outperforms three other models, with the GJR model showing a statistically significant improvement when using the skewed-t distribution over the normal distribution.

### EEM
- The APARCH model with a skewed-t distribution is statistically superior to the same model with a normal distribution and outperforms other tested models. The GARCH model with a skewed-t distribution also shows significant improvements over its normal distribution counterpart.
