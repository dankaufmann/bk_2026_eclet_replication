#-------------------------------------------------------------------------------
# 1_MultiDimShocks.R
# Multiple monetary policy shocks from daily data: A heteroskedasticity IV approach
# by Marc Burri and Daniel Kaufmann
#-------------------------------------------------------------------------------
# Identifies two-dimensional monetary policy shocks (target and path) via
# heteroskedasticity-based IV (HET-IV) and high-frequency proxy IV (HF-IV),
# estimates impulse responses via local projections, and tests equality of
# the two approaches via a pairs bootstrap.
#
# Input:  ./DataRep/Data.RData         (daily financial market data and event indicators)
# Output: ./Results/WeakData2Dim_*.mat (data for Matlab weak-instrument tests; Section 2)
#         ./Results/IRF_HET_*.pdf      (baseline HET-IV IRFs; Section 3)
#         ./Results/IRF_HF_HET_*.pdf   (HET-IV vs. HF-IV comparison IRFs; Section 3)
#         ./Results/IRF_HFRec_HET_*.pdf(HET-IV vs. recursive HF-IV IRFs; Section 3)
#         ./Results/CorrelationShocks2Dim.tex (shock correlation table; Section 3)
#         ./Results/PVals_HF_HET_Diff_*.pdf (bootstrap p-values; Section 4)
#
# Run order: this script first, then 2_MultiDimWeakTests.m (Matlab/Octave)
#-------------------------------------------------------------------------------

# Load user-defined commands and packages
#install.packages("devtools")
#library(devtools)
#remove.packages("hetiv")
#install("C:\\Users\\daenu\\GitHub\\hetiv")
library("hetiv")
library("tsbox")
library("xts")
library("R.matlab")
library("grid")
library("gridExtra")
library("ggplot2")
library("dplyr")
library("xtable")

# source("./Functions/IdentificationFunctions.R")
# library("MultiRNG")
# library("ivreg")
# library("sandwich")
# library("lubridate")
# library("xts")
# library("xtable")
# library("lmtest")
# library("tsbox")
# library("dplyr")
# library("matrixcalc")
# library("MASS")


#-------------------------------------------------------------------------------
# 0) Load daily data
#-------------------------------------------------------------------------------
load("./DataRep/Data.RData")

#-------------------------------------------------------------------------------
# 1) Set baseline settings
#-------------------------------------------------------------------------------
# Sample settings (FOMC events from Bauer and Swanson, extended with our own event collection)
myStart <- "1988-02-01"
myEnd   <- "2022-12-31"

# Tolerance for SVD truncation of inverse in Kalman-filter prediction step
# Set to NA for default values of Moore-Penrose inverse
tolInv = NA

# Choose whether do bootstrap tests (takes a while)
bootstrap = TRUE
B = 20             # Number of bootstrap iterations (set to 500 for test purposes, use 2000 for final results)

# Baseline model specification
P        <- 1                   # 1 lags as controls 
E        <- 2                   # 2 dimensional shock
depVar   <- c("IRLTfed", "NEER", "Stocks")  # Variables used to extract shocks and compute the IRFs (Y)
conVar   <- c("IRSTfed", "IRMTfed", "NEER", "Stocks",  "Spread", "TRSkew")   # Control variables included with up to P lags (X)

# IRF graph settings
cum        <- TRUE  # Cumulative IRFs
NormFac    <- 0.25  # Normalization of initial response for IRFs
HStep      <- 1     # Only estimate every HStep IRF (more efficient with daily data; for test purposes)
HTick      <- 5     # Number of periods for each tick mark in graph (should be multiple of HStep)
H          <- 21    # Number of periods for IRF horizon (should be multiple of HStep)
figScaleW  <- 2.2   # Figure width scaling factor (per column)
figScaleH  <- 2.6   # Figure height scaling factor (per row)

shockLabs <- c()             # Empty vector, shocks are labeled in the plot functions
xLab      <- "Working days"  # X-axis label for IRF graphs
varLabs   <- data.frame(IRSTfed = "Short-term rate (in pp)",IRMTfed = "Medium-term rate (in pp)", IRLTfed = "Long-term rate (in pp)",
                        IR1Yfed = "1Y rate (in pp)", IR2Yfed = "2Y rate (in pp)", IR3Yfed = "3Y rate (in pp)", IR5Yfed = "5Y rate (in pp)", 
                        IR7Yfed = "7Y rate (in pp)",IR10Yfed = "10Y rate (in pp)", IR30Yfed = "30Y rate (in pp)",
                        NEER = "Exchange rate (in %)",
                        TRSkew  = "Treasury yield skew.",
                        Stocks = "Stock prices (in %)",
                        Spread = "Corp. spread (in pp)", 
                        VIX = "VIX (in pp)")

#-------------------------------------------------------------------------------
# 2) Check across term structure which instrument is strongest
# - Use a selection of short-term interest rates for first shock
# - Exports WeakData*.mat files used by 2_MultiDimWeakTests.m
#-------------------------------------------------------------------------------
depVar1  <- c("NEER")           # Only needed so that N > 1, also act as dummy control vars
N        <- 3                   # Number of dependent variables (only for weak instrument tests)

# First shock try short-term interest rates
for(firstShock in c("FFR", "IR3Mfed", "IR6Mfed", "IRSTfed")){
  
  # Second shock try medium- and long-term interest rates
  for(secondShock in c("IR2Yfed", "IR3Yfed", "IR5Yfed", "IRMTfed")){
    
    irfVars  <- c(firstShock, secondShock, depVar1)
    irfVars  <- unique(irfVars[irfVars != ""])
    
    infoVars <- c(conVar)
    infoVars <- unique(infoVars[infoVars != ""])
    
    y       <- Data[, irfVars]
    O       <- Data[, infoVars]
    
    ZTG   <- Data[, "SFFR"]
    ZFG   <- Data[, "SFG"]
    Z     <- ts_c(ZTG, ZFG)
    
    # Shorten to sample period 
    Dt   <- index(ts_span(y, myStart, myEnd))
    y    <- data.frame(ts_span(y, myStart, myEnd))
    O    <- data.frame(ts_span(O, myStart, myEnd))
    Ind  <- data.frame(ts_span(IndE, myStart, myEnd))
    Z    <- data.frame(ts_span(Z, myStart, myEnd))
    
    # Estimate the impact matrix (H = 1)
    resLP <- hetiv(y =         y[, 1:N], 
                   O =         O, 
                   Ind =       Ind, 
                   P =         P,
                   E =         E,
                   H =         1,
                   details =   TRUE)
    resHF <- proxyiv(y =       y[, 1:N], 
                   O =         O,
                   Z =         Z,
                   Ind =       Ind,
                   recursive = FALSE,
                   P =         P,
                   E =         E,
                   H =         1,
                   details =   TRUE)
    
    # Write data for weak instrument test for heteroscedasticity-based instruments
    writeMat(con= paste0("./Results/WeakData2Dim_", firstShock, "_", secondShock, ".mat"), myTable = resLP$WeakData)
    writeMat(con= paste0("./Results/WeakData2Dim_HF_", firstShock, "_", secondShock, ".mat"), myTable = resHF$WeakData)
    
  }
}

#-------------------------------------------------------------------------------
# 3) Estimate the baseline model for IRFs
# - Baseline: short-term (IRSTfed) and medium-term (IRMTfed) rate as normalization vars
# - Estimates HET-IV, HF-IV (non-recursive), and HF-IV (recursive) models
# - Extracts shock series, computes correlations
#-------------------------------------------------------------------------------
normVar <- c("IRSTfed", "IRMTfed")
depVar  <- c("Spread", "NEER", "Stocks", "VIX")  

irfVars  <- c(normVar, depVar)                # Variables used for IRF 
irfVars  <- unique(irfVars[irfVars != ""])
infoVars <- c(conVar)                         # Control variables included with up to P lags
infoVars <- unique(infoVars[infoVars != ""])

N        <- length(irfVars)     

# - y contains variables to compute IRF for (y1 response at h = 0 normalized to unity)
# - O contains the information set, included with 1..P lags
y     <- Data[, irfVars]
O     <- Data[, infoVars]

ZTG   <- Data[, "SFFR"]
ZFG   <- Data[, "SFG"]
Z     <- ts_c(ZTG, ZFG)

# Shorten to sample period 
Dt    <- index(ts_span(y, myStart, myEnd))
y     <- data.frame(ts_span(y, myStart, myEnd))
O     <- data.frame(ts_span(O, myStart, myEnd))
Ind   <- data.frame(ts_span(IndE, myStart, myEnd))
Z     <- data.frame(ts_span(Z, myStart, myEnd))

# Estimate the models
resHET <- hetiv(y =        y, 
               O =         O, 
               Ind =       Ind, 
               P =         P,
               E =         E,
               H =         H,
               cum =       cum,
               details =   TRUE)
resHF <- proxyiv(y =         y, 
                 O =         O,
                 Z =         Z,
                 Ind =       Ind,
                 recursive = FALSE,
                 P =         P,
                 E =         E,
                 H =         H,
                 cum =       cum,
                 details =   TRUE)
resHFrec <- proxyiv(y =      y, 
                 O =         O,
                 Z =         Z,
                 Ind =       Ind,
                 recursive = TRUE,
                 P =         P,
                 E =         E,
                 H =         H,
                 cum =       cum,
                 details =   TRUE)

# Plot baseline IRFs (HET-IV)
myGraphs <- plotirf(IRFest =  resHET$irf, 
                    IRFse =   resHET$se, 
                    HTick =   HTick, 
                    Labels =  varLabs[irfVars], 
                    ci     = c(0.9, 0.95))
for(e in 1:E){
  
  myXLab <- textGrob(xLab, gp=gpar(fontsize=9.5), vjust = -1)
  
  Ng <- dim(resHET$Psi)[1]
  if(is.null(Ng)){
    Ng = length(resHET$Psi)
  }
  
  p1       <- grid.arrange(grobs = myGraphs[((e-1)*Ng+1):(e*Ng)], nrow = 2, ncol = ceiling(Ng/2), bottom = myXLab)
  ggsave(p1, file = paste0("./Results/IRF_HET_", e, ".pdf"), width = ceiling(Ng/2)*figScaleW, height = 1.7*figScaleH)
  
}

# Plot comparison IRFs (HET-IV vs. HF-IV)
myGraphs <- plot2irf( IRF1 =        resHET$irf, 
                      IRF1se =      resHET$se, 
                      IRF2 =        resHF$irf, 
                      IRF2se =      resHF$se, 
                      HTick =       HTick, 
                      Labels =      varLabs[irfVars], 
                      ci =          c(0.9, 0.95))

for(e in 1:E){
  
  myXLab <- textGrob(xLab, gp=gpar(fontsize=9.5), vjust = -1)
  
  Ng <- dim(resHET$Psi)[1]
  if(is.null(Ng)){
    Ng = length(resHET$Psi)
  }
  
  p1       <- grid.arrange(grobs = myGraphs[((e-1)*Ng+1):(e*Ng)], nrow = 2, ncol = ceiling(Ng/2), bottom = myXLab)
  ggsave(p1, file = paste0("./Results/IRF_HF_HET_", e, ".pdf"), width = ceiling(Ng/2)*figScaleW, height = 1.7*figScaleH)
  
}

# Plot comparison IRFs (HET-IV vs. HF-IV rec.)
myGraphs <- plot2irf( IRF1 =        resHET$irf, 
                      IRF1se =      resHET$se, 
                      IRF2 =        resHFrec$irf, 
                      IRF2se =      resHFrec$se, 
                      HTick =       HTick, 
                      Labels =      varLabs[irfVars], 
                      ci =          c(0.9, 0.95))

for(e in 1:E){
  
  myXLab <- textGrob(xLab, gp=gpar(fontsize=9.5), vjust = -1)
  
  Ng <- dim(resHET$Psi)[1]
  if(is.null(Ng)){
    Ng = length(resHET$Psi)
  }
  
  p1       <- grid.arrange(grobs = myGraphs[((e-1)*Ng+1):(e*Ng)], nrow = 2, ncol = ceiling(Ng/2), bottom = myXLab)
  ggsave(p1, file = paste0("./Results/IRF_HFRec_HET_", e, ".pdf"), width = ceiling(Ng/2)*figScaleW, height = 1.7*figScaleH)
  
}

# Compute shocks and their correlations
# Extract shocks based on the Kalman filter for heteroscedasticity-based shocks
HETShocks <- kfpredict(resHET$Sig, resHET$SigR, as.matrix(resHET$Psi), as.matrix(resHET$et), tol = tolInv)
HETShocks <- xts(HETShocks, order.by = Dt)
colnames(HETShocks) <- c("Target (HET-IV)", "Path (HET-IV)")

HFShocks <- kfpredict(resHF$Sig, resHF$SigR, as.matrix(resHF$Psi), as.matrix(resHF$et), tol = tolInv)
HFShocks <- xts(HFShocks, order.by = Dt)
colnames(HFShocks) <- c("Target (HF-IV)", "Path (HF-IV)")

HFrecShocks <- kfpredict(resHFrec$Sig, resHFrec$SigR, as.matrix(resHFrec$Psi), as.matrix(resHFrec$et), tol = tolInv)
HFrecShocks <- xts(HFrecShocks, order.by = Dt)
colnames(HFrecShocks) <- c("Target (HF-IV rec.)", "Path (HF-IV rec.)")

# Correlation matrix all shocks
AllShocks <- ts_c(HETShocks, HFShocks, HFrecShocks)
myCor <- round(cor(AllShocks, use = "pairwise.complete.obs"), 2)
myCor[upper.tri(myCor)] <- NA
print(myCor, file = "./Results/CorrelationShocks2Dim.txt")
print(myCor %>% xtable(caption = "Correlation matrix of shocks", 
                       label = "tab:correlation_shocks") ,
      file = "./Results/CorrelationShocks2Dim.tex")


#-------------------------------------------------------------------------------
# 5) Runs bootstrap equality test
#-------------------------------------------------------------------------------
set.seed(42)


# Bootstrap equality of impulse responses
if(bootstrap == TRUE){
  
  # Reestimate baseline
  resHET_o <- hetiv(y =        y, 
                  O =         O, 
                  Ind =       Ind, 
                  P =         P,
                  E =         E,
                  H =         H,
                  cum =       cum,
                  details =   TRUE)
  resHF_o <- proxyiv(y =         y, 
                   O =         O,
                   Z =         Z,
                   Ind =       Ind,
                   recursive = FALSE,
                   P =         P,
                   E =         E,
                   H =         H,
                   cum =       cum,
                   details =   TRUE)
  
  # Dims: H x N X E
  boot_diffs <- array(NA, dim = c(H, N, E, B))
  for(b in 1:B){
    cat(paste0("Bootstrap iteration ", b, " of ", B, "\n"))
    
    # Do bootstrap resampling
    idx   <- sample(nrow(y), replace = TRUE)
    y_b   <- y[idx, ]
    O_b   <- O[idx, ]
    Ind_b <- Ind[idx, ]
    Z_b   <- Z[idx, ]
    
    # Estimate on bootstrap samples
    resHET_b <- hetiv(y =        y_b, 
                      O =         O_b, 
                      Ind =       Ind_b, 
                      P =         P,
                      E =         E,
                      H =         H,
                      cum =       cum,
                      details =   TRUE)
    resHF_b <- proxyiv(y =         y_b, 
                       O =         O_b,
                       Z =         Z_b,
                       Ind =       Ind_b,
                       recursive = FALSE,
                       P =         P,
                       E =         E,
                       H =         H,
                       cum =       cum,
                       details =   TRUE)
    
    boot_diffs[, , , b] <- resHET_b$irf - resHF_b$irf
 
  }
  
  # Dims: H x N X E
  obs_diff <- resHET_o$irf - resHF_o$irf
  
  p_boot <- array(NA, dim = c(H, N, E))
  for(e in 1:E){
    for(n in 1:N){
      for(h in 1:H){
        p_boot[h, n, e]   <- mean(abs(boot_diffs[h, n, e, ] - mean(boot_diffs[h, n, e, ])) >= abs(obs_diff[h, n, e] - mean(boot_diffs[h, n, e, ])))
      }
    }
  }
  
  dimnames(p_boot)[1] <- dimnames(resHET_o$irf )[1]
  
  myGraphs <- plotpval(p_boot, HTick, varLabs[irfVars], c(0.05, 0.10))
  for(e in 1:E){
    
    myXLab <- textGrob(xLab, gp=gpar(fontsize=9.5), vjust = -1)
    
    Ng <- dim(resHET$Psi)[1]
    if(is.null(Ng)){
      Ng = length(resHET$Psi)
    }
    
    p1       <- grid.arrange(grobs = myGraphs[((e-1)*Ng+1):(e*Ng)], nrow = 2, ncol = ceiling(Ng/2), bottom = myXLab)
    ggsave(p1, file = paste0("./Results/PVals_HF_HET_Diff_", e, ".pdf"), width = ceiling(Ng/2)*figScaleW, height = 1.7*figScaleH)
    
  }
}


#-------------------------------------------------------------------------------
# End of code
#-------------------------------------------------------------------------------



