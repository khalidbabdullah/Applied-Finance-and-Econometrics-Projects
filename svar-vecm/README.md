# The Macroeconomic Impact of Fiscal Policy  
### Evidence from SVAR and VECM Analysis (U.S., 2000–2019)

## Overview
This project studies the dynamic effects of U.S. fiscal policy on real GDP using
Structural Vector Autoregression (SVAR) and Vector Error Correction Models (VECM).
The analysis focuses on identifying both short-run fiscal shocks and long-run equilibrium
relationships between government spending, tax revenue, and output.

## Research Question
How do government spending and tax shocks affect real GDP in the short run and long run,
and do these effects differ before and after the Global Financial Crisis?

## Data
Quarterly U.S. macroeconomic data (2000–2019) sourced from the Federal Reserve Economic
Data (FRED) database, including:
- Real GDP
- Government spending
- Tax revenue

Nominal variables are deflated using the GDP deflator and transformed into logarithms
for elasticity-based interpretation.

## Methodology
- Unit root testing (DF-GLS)
- Johansen cointegration analysis
- Vector Error Correction Model (VECM)
- Structural VAR (SVAR) with Blanchard–Perotti identification
- Impulse response functions (IRFs)
- Forecast error variance decomposition (FEVD)
- Structural break and subsample analysis (pre- vs post-GFC)
- Extensive residual diagnostics and robustness checks

## Key Findings
- Strong evidence of cointegration among GDP, spending, and taxes
- Government spending shocks produce modest and short-lived output responses
- Tax shocks generate larger and more persistent negative effects on GDP
- Fiscal multipliers appear stronger in the post-2008 period
- Results are robust across alternative specifications and diagnostics

## Files
- `analysis.R` — full R source code for data collection, SVAR and VECM estimation,
  diagnostics, and impulse response analysis
- `paper.pdf` — research paper presenting theory, methodology, results, and interpretation

## Tools
- R (vars, svars, urca, fredr, ggplot2)
