# Multiple monetary policy shocks from daily data: A heteroskedasticity IV approach

Replication files — Version 28 May 2026

Marc Burri and Daniel Kaufmann  
University of Neuchâtel, Institute of Economic Research  
*Economics Letters*, resubmission

---

## Overview

This replication package reproduces all figures and tables in Burri and Kaufmann (2026). The paper proposes a heteroskedasticity-based instrumental variables (HET-IV) estimator that identifies two-dimensional monetary policy shocks — a target shock and a path shock — from daily financial market data, without requiring intraday tick data or precise announcement timestamps.

All estimation routies are provided in the R package `hetiv`, which can be downloaded from: https://github.com/dankaufmann/hetiv. Full documentation of the package and an example application is available on: https://dankaufmann.github.io/hetiv/.

The replication code is split across two scripts that must be run in order:

1. `0_GetData.R` — *OPTIONAL* R script that creates the data set. Requires a FRED API key.
1. `1_MultiDimShocks.R` — R script that estimates local projection models, produces IRF figures, makes shock predictions, and conducts the weak instrument tests
2. `2_MultiDimWeakTests.m` *OPTIONAL* Octave weak-instrument test script to verify the results using the original Lewis and Mertens (2025) codes. The codes have been modified to run in Octave. The results have been verified in Matlab.

---

## File structure

| Path | Description |
|---|---|
| `0_GetData.R` | OPTIONAL Data replication file  |
| `1_MultiDimShocks.R` | Main R replication script |
| `2_MultiDimWeakTests.m` | OPTIONAL Octave weak-instrument test script |
| `Data/Data.RData` | Processed daily data (R format) |
| `Data/*.xlsx` | Various raw data files (Excel format) |
| `Functions/gweakivtest.m` | Lewis and Mertens (2025) weak-instrument test (Matlab/Octave) |
| `Functions/gweakivtest_critical_values.m` | Critical value computation for Lewis-Mertens test (Matlab/Octave) |
| `Results/` | All output files (figures, tables, weak instrument test results) |

---

## Data description

All data follow Burri and Kaufmann (2026, IRENE Working Paper 24-03). The sample covers 1988-01-04 to 2025-12-31.

### Financial market variables

| Variable | Description | Source | Transformation |
|---|---|---|---|
| `FFR` | Federal Funds Rate | Board of Governors (Table H.15, from FRED) | First difference |
| `IR3Mfed` | 3-month Treasury bill yield | Board of Governors (Table H.15, from FRED) | First difference |
| `IR6Mfed` | 6-month Treasury bill yield | Board of Governors (Table H.15, from FRED) | First difference |
| `IR2Yfed` | 2-year Treasury bond yield | Board of Governors (Table H.15, from FRED) | First difference |
| `IR3Yfed` | 3-year Treasury bond yield | Board of Governors (Table H.15, from FRED) | First difference |
| `IR5Yfed` | 5-year Treasury bond yield | Board of Governors (Table H.15, from FRED) | First difference |
| `IR7Yfed` | 7-year Treasury bond yield | Board of Governors (Table H.15, from FRED) | First difference |
| `IR10Yfed` | 10-year Treasury bond yield | Board of Governors (Table H.15, from FRED) | First difference |
| `IR30Yfed` | 30-year Treasury bond yield | Board of Governors (Table H.15, from FRED) | First difference |
| `IRSTfed` | Short-term rate (average of FFR, 3M, 6M) | Constructed | First difference |
| `IRMTfed` | Medium-term rate (average of 2Y, 3Y, 5Y) | Constructed | First difference |
| `IRLTfed` | Long-term rate (average of 7Y, 10Y, 30Y) | Constructed | First difference |
| `NEER` | Nominal effective exchange rate (USD per unit of foreign currency) | Board of Governors (Table H.10, from FRED); FRED series DTWEXM linked to DTWEXAFEGS | Log-difference × 100 |
| `Stocks` | Stock price index (weighted avg: 0.5 S&P 500 + 0.25 NASDAQ Composite + 0.25 NASDAQ 100) | S&P Dow Jones / Yahoo Finance, NASDAQ OMX (FRED series: NASDAQCOM, NASDAQ100) | Log-difference × 100 |
| `Spread` | Corporate bond spread (average of AAA and BAA, relative to 10Y Treasury) | Moody's (FRED series: BAA10Y, AAA10Y) | First difference |
| `VIX` | CBOE Volatility Index (average of VIX close and VXO close) | CBOE (FRED series : VIXCLS, VXOCLS) | First difference |
| `TRSkew` | Implied Treasury yield skewness | Bauer and Chernov (2024, *Journal of Finance*) | Level |

**Notes:** Holidays and weekends are excluded. Remaining missing values are linearly interpolated up to 4 working days. Transformations are applied after treatment of missing values. Differences are computed compared to the last available observation.

### High-frequency monetary policy surprises (instruments for HF-IV)

| Variable | Description | Source |
|---|---|---|
| `SFFR` | Federal funds rate surprise | Swanson (2021)  |
| `SFG` | Forward guidance surprise | Swanson (2021) |

**Notes:** Swanson (2021) surprises are extracted from an intraday event window around FOMC announcements. They include announcements from 5 July 1991 — 19 June 2019. For comparability, parameters of the local projections are estimated over the entire sample.

### FOMC announcement dates (event indicator)

| Period | Source |
|---|---|
| 1988–2019 | Bauer and Swanson (2022)|
| 2020–2025 | Acosta et al. (2025) |

**Notes:** All other weekdays, excluding holidays, are used as control days.

---

## Results produced

| File | Description | Paper output |
|---|---|---|
| `Results/IRF_HET_1.pdf` | HET-IV target shock IRFs | Figure 1(a) |
| `Results/IRF_HET_2.pdf` | HET-IV path shock IRFs | Figure 1(b) |
| `Results/IRF_HF_HET_1.pdf` | HET-IV vs. HF-IV, target shock | Figure 2(a) |
| `Results/IRF_HF_HET_2.pdf` | HET-IV vs. HF-IV, path shock | Figure 2(b) |
| `Results/IRF_HFRec_HET_1.pdf` | HET-IV vs. recursive HF-IV, target shock | Online Appendix |
| `Results/IRF_HFRec_HET_2.pdf` | HET-IV vs. recursive HF-IV, path shock | Online Appendix |
| `Results/PVals_HF_HET_Diff_1.pdf` | Bootstrap p-values, target shock | Online Appendix |
| `Results/PVals_HF_HET_Diff_2.pdf` | Bootstrap p-values, path shock | Online Appendix |
| `Results/CorrelationShocks2Dim.tex` | Shock correlation matrix | Online Appendix |
| `Results/WeakIVTest_HET1Dim.txt` | HET-IV target shock, 1 dimension weak-instrument test | Table 1(a) |
| `Results/WeakIVTest_HET2Dim.txt` | HET-IV path shock, 2 dimension weak-instrument test | Table 2(a) |
| `Results/WeakIVTest_HF1Dim.txt` | HF-IV target shock, 1 dimension weak-instrument test | Table 1(b) |
| `Results/WeakIVTest_HF12Dim.txt` | HF-IV path shock, 1 dimension weak-instrument test | Table 2(b) |
| `Results/WeakIVTest_HF2Dim.txt` | HF-IV target and path shock, 2 dimension weak-instrument test | Not reported |
|`WeakData2Dim_XXX.mat` | Matlab data files used by Octave/Matlab script to reproduce weak instrument tests using | Not reported |


---

## Software requirements

### R (scripts `1_MultiDimShocks.R`,  `0_GetData.R`)

Tested with R 4.x. Required packages:

```
hetiv, ggplot2, gridExtra, ivreg, sandwich, expm, 
xts, xtable, lmtest, tsbox, dplyr, matrixcalc, MASS, R.matlab, grid
```

The R package `hetiv` can be downloaded from: https://github.com/dankaufmann/hetiv.

### Matlab / Octave (script `2_MultiDimWeakTests.m`)

The script was adapted for GNU Octave (tested with Octave 8.x) with the `optim` and `statistics` packages. The original code by Lewis and Mertens (2025) is written for Matlab with the Optimization Toolbox.

The weak-instrument test functions (`gweakivtest.m`, `gweakivtest_critical_values.m`) are by Lewis and Mertens (2025) and are included in `Functions/`.

---

## How to replicate

1. Open R, install the `hetiv` package using 
```
install.packages("remotes")
remotes::install_github("dankaufmann/hetiv")
```
2. In R, set the working directory to `bk_2026_eclet_replication/`, and run `1_MultiDimShocks.R`. This produces all results. If desired, it also exports `.mat` files to `Results/`.

OPTIONAL: Run `0_GetData.R` to create the data set from scratch. Use the `exportMat` option to export Matlab data files for the weak instrument tests. Run `2_MultiDimWeakTests.m`. This reads the `.mat` files and writes the weak-instrument test summaries using the original Lewis and Mertens (2025) code to `Results/`.

**Notes:** The bootstrap in step 1 takes approximately 30–60 minutes for `B = 500` depending on hardware. Set to `B = 2000` for accurate results. Set `bootstrap = FALSE` to skip it.

---

## References

- Bauer, M.D. and Chernov, M. (2024). Interest rate skewness and biased beliefs. *Journal of Finance*, 79(1):173–217.
- Bauer, M.D. and Swanson, E.T. (2022). A reassessment of monetary policy surprises and high-frequency identification. In *NBER Macroeconomics Annual 2022*, volume 37.
- Burri, M. and Kaufmann, D. (2026). Measuring monetary policy shocks. IRENE Working Paper 24-03, University of Neuchâtel, https://ideas.repec.org/p/irn/wpaper/24-03.html.
- Lewis, D.J. and Mertens, K. (2025). A robust test for weak instruments with multiple endogenous regressors. *Review of Economic Studies*.
- Rigobon, R. and Sack, B. (2004). The impact of monetary policy on asset prices. *Journal of Monetary Economics*, 51(8):1553–1575.
- Swanson, E.T. (2021). Measuring the effects of Federal Reserve forward guidance and asset purchases on financial markets. *Journal of Monetary Economics*, 118:32–53.
