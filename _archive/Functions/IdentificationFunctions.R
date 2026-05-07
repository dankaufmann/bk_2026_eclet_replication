#-------------------------------------------------------------------------------
# IdentificationFunctions.R
# Multiple monetary policy shocks from daily data: A heteroskedasticity IV approach
# by Marc Burri and Daniel Kaufmann
#-------------------------------------------------------------------------------
# Functions for HET-IV and HF-IV identification of monetary policy shocks
# via local projections, following Rigobon and Sack (2004) and Lewis (2022).
#
# Functions included:
#   estimLPHet       - LP-IV estimation using heteroscedasticity-based instruments (recursive)
#   estimLPProxy     - LP-IV estimation using high-frequency proxy instruments
#   estimLPHetFAST   - Fast version for bootstrap
#   estimLPProxyFAST - Fast version for bootstrap
#   extractShocksKF  - Kalman-filter shock extraction from LP-IV estimates
#   plotPvals        - Plot bootstrap p-values for IRF equality tests
#   plot2IRFs        - Plot two sets of IRFs for comparison (with 90% bands)
#   plotIRFs         - Plot a single set of IRFs (with 90% and 95% bands)
#
# Dependencies: tidyverse, dplyr, ggplot2, ivreg, sandwich, lmtest, MASS, matrixcalc
#-------------------------------------------------------------------------------

estimLPHet <- function(y, O, Ind, Interact, E, P, cumRep, H, HStep, NormFac){
  # Estimates impulse responses with recursive heteroscedasticity-IV
  # Ind:    is an event indicator, with values
  #         0: Control day (no event)
  #         1: Policy day (event)
  #         2: Contaminated control day (other event)
  #         Impulse responses are only estimated based on policy and control days
  # y:      Outcome variables, for which the impulse responses are computed.
  #         Note that for each dimension the first variable's response is
  #         normalized to unity. These are also the variables used for
  #         instruments and recursive ordering
  # O:      Information set which can, but does not have to be, identical to y.
  #         It will be included lagged 1...P
  # Interact: If TRUE, then O is interacted with regime dummies
  # P:      Maximum lag order for information set. If P = 0, no information set
  # H:      Maximum horizon up to which impulse response is computed
  # HStep:  If > 1, only every HStep impulse response is estimated
  #         (starting at H = 1)
  # cumRep: Vector of T/F of the same dimension as y indicating whether to compute
  #         the cumulative impulse response for every variable included in y
  # E:      Number of dimensions of shocks to be identified via Cholesky
  # NormFac:Normalizes the impulse response of the first variable to a specific value
  #         Set to 1 if no normalization is desired

  # Collect various properties of the data and observations to be used
  Nobs <- dim(y)[1]     # Number of observations in y
  beg  <- max(P+1)      # Start of the sample
  end  <- Nobs - H+1    # End of the sample
  N    <- dim(y)[2]     # Number of variables in y
  M    <- dim(O)[2]     # Number of variables in O

  if(sum(is.na(y))>0){
    warning("Missing values in y")
  }
  if(sum(is.na(O))>0){
    warning("Missing values in O")
  }
  if(sum(is.na(Ind))>0){
    warning("Missing values in Ind")
  }

  # Specify at which horizons IRF should be computed
  HSeries <- seq(1, H, HStep)
  HNum    <- length(HSeries)

  # Set up data set and and objects to save results
  DataM     <- data.frame(y, O, Ind)
  colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), "Ind")
  irfest    <- array(NA, dim = c(HNum, N, E))
  irfse     <- array(NA, dim = c(HNum, N, E))
  IVRes     <- list()

  # Compute lagged variables included in information set O(t-1)
  infoVars      <- c()
  if(P > 0){
    for(j in 1:M){
      for (p in 1:P){
        if(Interact == TRUE){
          for(e in 0:1){
            DataM[, paste0("Event", e, ".o", j, ".l", p)] <- dplyr::lag(DataM[, paste0("o", j)], p) * (Ind == e)
            infoVars <- c(infoVars, paste0("Event", e, ".o", j, ".l", p))
          }
        }else{
          DataM[, paste0("o", j, ".l", p)] <- dplyr::lag(DataM[, paste0("o", j)], p)
          infoVars <- c(infoVars, paste0("o", j, ".l", p))
        }

      }
    }
  }
  if(P == 0){
    # Otherwise, only regress on a constant
    if(Interact == TRUE){
      for(e in 0:1){
        DataM[, paste0("Event", e)] <- (Ind == e)
        infoVars <- c(infoVars, paste0("Event", e))
      }
    }else{
      infoVars <- "1"
    }

  }

  # Compute instrument separately for every dimension
  for(e in 1:E){

    # Set up the shock variable (instrumented variable), which is the e th variable
    # in y
    DataM$shockVar <- DataM[, paste0("y", e)]

    # Set up the instrument variable, which is also the e th variable in y
    DataM$Ze <- DataM[,paste0("y", e)]

    # Set up the control variables for recursive ordering
    # Note that we interact the dependent variable from previous
    # equation with the event dummy
    recVars <- c()
    recInst <- c()
    if(e>1){
      for(q in 1:(e-1)){
        recVars <- c(recVars, paste0("y", q))
        recInst <- c(recInst, paste0("Z", q))
      }
    }

    # Set up control variables for the local projections
    controls.info <- unique(c(infoVars))
    controls.info <- controls.info[controls.info !=""]

    controls.lp <- unique(c(infoVars, recVars))
    controls.lp <- controls.lp[controls.lp !=""]

    controls.iv <- unique(c(infoVars, recInst))
    controls.iv <- controls.iv[controls.iv !=""]

    # Account for missing data in instrument variable by setting it to an other event
    # which is not used in the estimation
    DataM$Ind[is.na(DataM$Ze)] = 2

    # Set up event days (Policy event, control, and other event)
    DataM$Event    <- (DataM$Ind == 1)
    DataM$NoEvent  <- (DataM$Ind == 0)
    DataM$OthEvent <- (DataM$Ind == 2)

    # Shorten data to subset of observations without missing values at beginning or end
    DataMSub <- DataM[beg:end, ]


    # Orthogonalize the instrument if we include lags of dependent variable
    # See Lewis (2022)
    if(P > 0){

      #print(controls.info)

      myFormula = paste0("Ze ~ ", paste(controls.info, collapse = "+"))
      orthModel <- lm(as.formula(myFormula), data = subset(DataMSub, Ind < 2), na.action="na.exclude")

      #print(summary(orthModel))
      DataMSub$Ze <- residuals(orthModel, na.action="na.exclude")

    }


    # Compute instrument and F-Statistic (see Lewis, 2022, ReStat, and
    # Rigobon, 2003, ReStat)
    # Compute number of events, control observations, and non-events
    Te <- sum(DataMSub$Event)
    Tn <- sum(DataMSub$NoEvent)
    To <- sum(DataMSub$OthEvent)
    Tt <- Te+Tn

    # Compute the instrument (see Lewis, 2022, ReStat and Rigobon and Sack, 2004, JME)
    DataMSub$Z <- (DataMSub$Event*Tt/Te - DataMSub$NoEvent*Tt/Tn)*DataMSub$Ze

    # Save instrument in original data for later use
    DataM[beg:end, paste0("Z", e)] <- DataMSub$Z
    DataM[beg:end, "Z"] <- DataMSub$Z

    # Estimate the impulse responses for every variable in y
    for (i in 1:N)
    {
      # Set up dependent variable
      DataM$depVar  <- DataM[, paste0("y", i)]
      cumRepi       <- cumRep[i]

      # Create object to save impulse responses and standard errors
      IRF <- array(data = NA, dim = c(HNum, 7))
      colnames(IRF) <- c("H", "irf", "se", "upper95", "lower95", "upper90", "lower90")
      rownames(IRF) <- HSeries

      # Compute the impulse response for every horizon in a paralleled for loop
      # The reason why this is possible is that we directly estimate cumulative
      # impulse responses
     for(h in HSeries){
       # Compute dependent variable (option for cumulative responses)
       # That is, we compute y(t+h) in order to compute the direct forecast later
       # (local projection)
       # NOTE: For parallel computing, need to do accumulation every time
       if(cumRepi == T){
         for(f in 1:h){
           if(f == 1){
             DataM$depVar.h <- dplyr::lead(DataM$depVar, f-1)
           }else{
             DataM$depVar.h <- DataM$depVar.h + dplyr::lead(DataM$depVar, f-1)
           }
         }
       }else{
         DataM$depVar.h <- dplyr::lead(DataM$depVar, h-1)
       }
  
       # Shorten data to subset which contains no missing values
       DataMSub <- DataM[beg:end, ]
  
       # LP (Jorda, 2005), dependent variable is just y(t+h)
       # Control variables for direct forecast y(t+h|t-1)
       # Compute direct forecast t+h|t-1, where shockVar is instrumented
       #print(controls.lp)
       #print(controls.iv)
  
       myFormula <- paste0("depVar.h ~ shockVar + ",
                           paste(controls.lp, collapse = "+"),
                           paste0("| Z + ",
                                  paste(controls.iv, collapse = "+")))
  
       IV.mod  <- ivreg(as.formula(myFormula), data = subset(DataMSub, Ind < 2))
       #IV.se   <- sqrt(diag(NeweyWest(IV.mod, prewhite = FALSE)))
       IV.se   <- sqrt(diag(vcovHC(IV.mod, type = "HC0")))
       #IV.sum  <- summary(IV.mod, vcov = NeweyWest(IV.mod, prewhite = FALSE), df = Inf, diagnostics = TRUE)
       IV.sum  <- summary(IV.mod, vcov = vcovHC(IV.mod, type = "HC0"), df = Inf, diagnostics = TRUE)
  
       # Compute residuals of dependent variables that we need later on for shock extraction
       # Need to do this only once for all shocks
       # Here, we calculate the residual for variable i
       # Only for event days
       if (e == 1 & h == 1){
         # Workaround to get the correct residuals
         # 1) set depVar.h to missing if Ind == 2
         DataMSub$depVar.h2 <- DataMSub$depVar.h
         DataMSub$depVar.h2[DataMSub$Ind == 2] <- NA
  
         # Save residuals on entire data set (missing depVar leads to missing residual)
         myFormula <- paste0("depVar.h2 ~", paste(controls.info, collapse = "+"))
         OLS.mod  <- lm(as.formula(myFormula), data = subset(DataMSub, Ind < 2), na.action="na.exclude")
  
         # Save residuals on event days for later use
         eti  <- residuals(OLS.mod)
         eti[DataMSub$Event != 1] <- NA
  
         # Save residuals on other days for later use
         vti  <- residuals(OLS.mod)
         vti[DataMSub$NoEvent != 1] <- NA
       }else{
         eti  <- NA
         vti  <- NA
  
         OLS.mod <- NA
       }
  
       # Account for missing values when saving back to original data set
       DataM$eti <- NA
       DataM$eti[beg:end] <- eti
       eti <- DataM$eti
  
       DataM$vti <- NA
       DataM$vti[beg:end] <- vti
       vti <- DataM$vti
  
       Result <- list(IV.mod, IV.se, IV.sum, eti, vti)

        # Save coefficient, standard error, as well as summary statistics
        IRF[h, "irf"]   <- IV.mod$coefficients["shockVar"]
        IRF[h, "se"]    <- IV.se["shockVar"]

        # Normalizes IRFs to a specific value. Because initial response
        # normalized to unity, just multiply by NormFac
        irfest[h, i, e] <- IRF[h, "irf"]*NormFac
        irfse[h, i, e]  <- IRF[h, "se"]*NormFac

        # Save residuals for later computation of variance-covariance matrix
        if(h == 1 & e == 1){
          if(i == 1){
            et <- data.frame(eti)
            vt <- data.frame(vti)

          }else{
            et <- data.frame(et, eti)
            vt <- data.frame(vt, vti)
          }
        }

        # Save IV results for every varialbe and every horizon
        IVRes[[paste0("IV.h", h, ".n", i, ".e", e)]] <- IV.mod

      }

      # Label the rows of the impulse responses to match those that have been
      # estimated and to start with 0 (immediate response)
      dimnames(irfest)[[1]] <- HSeries-1
      dimnames(irfse)[[1]]  <- HSeries-1

    }
  }

  # Compute variance-covariance matrix of residuals on event days, and impact matrix
  # needed for shock extraction
  # Compute VCOV on event days only, Note that this is correct, because residuals are set to missing before for non-event days
  Sig  <- var(et, use = "complete.obs")

  if(sum(!is.na(vt)>0)){
    SigR <- var(vt, use = "complete.obs")
  }else{
    SigR <- NA
  }

  Psi  <- irfest[1 , ,]


  # Save data for weak instruments test by Lewis Mertens (2024) later done in Matlab
  # Note that this does not work if there are observations with Ind = 2
  # Note that it assumes that hte first E variables are endogenous
  if(controls.info[1] != "1"){
    WeakData <- data.frame(DataM[, paste0("y", 1:E)], DataM[, paste0("Z", 1:E)], DataM[, controls.info])
  }else{
    WeakData <- data.frame(DataM[, paste0("y", 1:E)], DataM[, paste0("Z", 1:E)])
  }
  # Workaround for labelling
  if(E == 1){
    if(controls.info[1] != "1"){
      colnames(WeakData) <- c("y1", "Z1", controls.info)
    }else{
      colnames(WeakData) <- c("y1", "Z1")
    }
  }

  # Set data to missing if indicator is equal to 2
  WeakData[DataM$Ind == 2, ] <- NA

  Method = "Heteroscedasticity-IV"
  Obs <- data.frame(Tp = Te, Tc = Tn, To = To, Tt =Tt)

  return(list(irf = irfest, se = irfse,
              IVRes = IVRes,
              Obs = Obs, Method = Method,
              et = et, Sig = Sig, SigR = SigR, Psi = Psi, WeakData = WeakData))

}


estimLPProxy <- function(y, O, Z, recursive = F, Ind, E, P, cumRep, H, HStep, NormFac){
  # Estimates impulse responses with recursive or non-recursive proxy-LP
  # Ind:    is an event indicator, with values
  #         0: Control day (no event)
  #         1: Policy day (event)
  #         2: Contaminated control day (other event)
  #         Impulse responses are only estimated based on policy (and control) days
  # y:      Outcome variables, for which the impulse responses are computed.
  #         Note that for each dimension the first variable's response is
  #         normalized to unity. These are also the variables used for
  #         instruments and recursive ordering
  # Z:      Is a matrix of instruments for the first E endogenous variables in y
  # recursive T/F whether we additionally impose the recursive zero restrictions
  # O:      Information set which can, but does not have to be, identical to y.
  #         It will be included lagged 1...P
  # P:      Maximum lag order for information set. If P = 0, no information set
  # H:      Maximum horizon up to which impulse response is computed
  # HStep:  If > 1, only every HStep impulse response is estimated
  #         (starting at H = 1)
  # cumRep: Vector of T/F of the same dimension as y indicating whether to compute
  #         the cumulative impulse response for every variable included in y
  # E:      Number of dimensions of shocks to be identified via Cholesky
  # NormFac:Normalizes the impulse response of the first variable to a specific value
  #         Set to 1 if no normalization is desired

  # Collect various properties of the data and observations to be used
  Nobs <- dim(y)[1]     # Number of observations in y
  beg  <- max(P+1)      # Start of the sample
  end  <- Nobs - H+1    # End of the sample
  N    <- dim(y)[2]     # Number of variables in y
  M    <- dim(O)[2]     # Number of variables in O

  if(sum(is.na(y))>0){
    warning("Missing values in y")
  }
  if(sum(is.na(O))>0){
    warning("Missing values in O")
  }
  if(sum(is.na(Ind))>0){
    warning("Missing values in Ind")
  }

  # Specify at which horizons IRF should be computed
  HSeries <- seq(1, H, HStep)
  HNum    <- length(HSeries)

  # Set up data set and and objects to save results
  DataM     <- data.frame(y, O, Ind, Z)
  colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), "Ind", paste0("z", 1:E))
  irfest    <- array(NA, dim = c(HNum, N, E))
  irfse     <- array(NA, dim = c(HNum, N, E))
  IVRes     <- list()

  # Compute lagged variables included in information set O(t-1)
  infoVars      <- c()
  if(P > 0){
    for(j in 1:M){
      for (p in 1:P){
        DataM[, paste0("o", j, ".l", p)] <- dplyr::lag(DataM[, paste0("o", j)], p)
        infoVars <- c(infoVars, paste0("o", j, ".l", p))
      }
    }
  }
  if(P == 0){
    # Otherwise, only regress on a constant
    infoVars <- "1"
  }

  # Compute instrument separately for every dimension
  for(e in 1:E){

    # Set up the shock variable (instrumented variable), which is the e th variable
    # in y
    DataM$shockVar <- DataM[, paste0("y", e)]

    # Set up the instrument variable, which is the e th variable in z
    DataM$Ze <- DataM[,paste0("z", e)]

    # Set up the control variables for recursive ordering
    # Note that we interact the dependent variable from the previous
    # equation with the event dummy
    recVars <- c()
    recInst <- c()
    if(recursive == T){
      if(e>1){
        for(q in 1:(e-1)){
          recVars <- c(recVars, paste0("y", q))
          recInst <- c(recInst, paste0("Z", q))
        }
      }
    }

    # Set up control variables for the local projections
    controls.info <- unique(c(infoVars))
    controls.info <- controls.info[controls.info !=""]

    controls.lp <- unique(c(infoVars, recVars))
    controls.lp <- controls.lp[controls.lp !=""]

    controls.iv <- unique(c(infoVars, recInst))
    controls.iv <- controls.iv[controls.iv !=""]

    # Account for missing data in instrument variable by setting it to an other event
    # which is not used in the estimation
    # Set up event days (Policy event, control, and other event)
    DataM$Event    <- (DataM$Ind == 1)
    DataM$NoEvent  <- (DataM$Ind == 0)
    DataM$OthEvent <- (DataM$Ind == 2)

    # Shorten data to subset of observations without missing values
    DataMSub <- DataM[beg:end, ]

    # Collect some statistics about instrument and save instrument for later use
    Te <- sum(DataMSub$Event)
    Tn <- sum(DataMSub$NoEvent)
    To <- sum(DataMSub$OthEvent)
    Tt <- Te+Tn

    # Get the instrument for shock e
    DataMSub$Z <- DataMSub$Ze

    # Save instrument in original data for later use
    DataM[beg:end, paste0("Z", e)] <- DataMSub$Z
    DataM[beg:end, "Z"] <- DataMSub$Z

    # Estimate the impulse responses for every variable in y
    for (i in 1:N)
    {
      # Set up dependent variable
      DataM$depVar  <- DataM[, paste0("y", i)]
      cumRepi       <- cumRep[i]

      # Create object to save impulse responses and standard errors
      IRF <- array(data = NA, dim = c(HNum, 7))
      colnames(IRF) <- c("H", "irf", "se", "upper95", "lower95", "upper90", "lower90")
      rownames(IRF) <- HSeries

      # Compute the impulse response for every horizon in a paralleled for loop
      for(h in HSeries){
         # Compute dependent variable (option for cumulative responses)
         # That is, we compute y(t+h) in order to compute the direct forecast later
         # (local projection)
         # NOTE: For parallel computing, need to do accumulation every time
         if(cumRepi == T){
           for(f in 1:h){
             if(f == 1){
               DataM$depVar.h <- dplyr::lead(DataM$depVar, f-1)
             }else{
               DataM$depVar.h <- DataM$depVar.h + dplyr::lead(DataM$depVar, f-1)
             }
           }
         }else{
           DataM$depVar.h <- dplyr::lead(DataM$depVar, h-1)
         }
  
         # Shorten data to subset which contains no missing values
         DataMSub <- DataM[beg:end, ]
  
         # LP (Jorda, 2005), dependent variable is just y(t+h)
         # Control variables for direct forecast y(t+h|t-1)
         # Compute direct forecast t+h|t-1, where shockVar is instrumented
         myFormula <- paste0("depVar.h ~ shockVar + ",
                             paste(controls.lp, collapse = "+"),
                             paste0("| Z + ",
                                    paste(controls.iv, collapse = "+")))
  
         IV.mod  <- ivreg(as.formula(myFormula), data = subset(DataMSub, Ind < 2))
         #IV.se   <- sqrt(diag(NeweyWest(IV.mod, prewhite = FALSE)))
         IV.se   <- sqrt(diag(vcovHC(IV.mod, type = "HC0")))
         #IV.sum  <- summary(IV.mod, vcov = NeweyWest(IV.mod, prewhite = FALSE), df = Inf, diagnostics = TRUE)
         IV.sum  <- summary(IV.mod, vcov = vcovHC(IV.mod, type = "HC0"), df = Inf, diagnostics = TRUE)
  
         # Compute residuals of dependent variables that we need later on for shock extraction
         # Need to do this only once for all shocks
         # Here, we calculate the residual for variable i
         # Only for event days
         if (e == 1 & h == 1){
           # Workaround to get the correct residuals
           # 1) set depVar.h to missing if Ind == 2
           DataMSub$depVar.h2 <- DataMSub$depVar.h
           DataMSub$depVar.h2[DataMSub$Ind == 2] <- NA
  
           # Save residuals on entire data set (missing depVar leads to missing residual)
           myFormula <- paste0("depVar.h2 ~", paste(controls.info, collapse = "+"))
           OLS.mod  <- lm(as.formula(myFormula), data = DataMSub, na.action="na.exclude")
  
           # Save back for later use
           eti  <- residuals(OLS.mod)
           eti[DataMSub$Event != 1] <- NA
  
           # Save residuals on other days for later use
           vti  <- residuals(OLS.mod)
           vti[DataMSub$NoEvent != 1] <- NA
  
         }else{
           eti  <- NA
           vti  <- NA
         }
  
         # Account for missing values when saving back to original data set
         DataM$eti <- NA
         DataM$eti[beg:end] <- eti
         eti <- DataM$eti
  
         DataM$vti <- NA
         DataM$vti[beg:end] <- vti
         vti <- DataM$vti
  
        # Save coefficient, standard error, as well as summary statistics
        IRF[h, "irf"]   <- IV.mod$coefficients["shockVar"]
        IRF[h, "se"]    <- IV.se["shockVar"]

        # Normalizes IRFs to a specific value. Because initial response
        # normalized to unity, just multiply by NormFac
        irfest[h, i, e] <- IRF[h, "irf"]*NormFac
        irfse[h, i, e]  <- IRF[h, "se"]*NormFac

        # Save residuals for later computation of variance-covariance matrix
        if(h == 1 & e == 1){
          if(i == 1){
            et <- data.frame(eti)
            vt <- data.frame(vti)
          }else{
            et <- data.frame(et, eti)
            vt <- data.frame(vt, vti)
          }
        }

        # Save IV results for every varialbe and every horizon
        IVRes[[paste0("IV.h", h, ".n", i, ".e", e)]] <- IV.mod

      }

      # Label the rows of the impulse responses to match those that have been
      # estimated and to start with 0 (immediate response)
      dimnames(irfest)[[1]] <- HSeries-1
      dimnames(irfse)[[1]]  <- HSeries-1

    }
  }

  # Save data for weak instruments test by Lewis Mertens (2024) later done in Matlab
  # Note that this does not work if there are observations with Ind = 2
  # Note that it assumes that the first E variables are endogenous
  if(controls.info[1] != "1"){
    WeakData <- data.frame(DataM[, paste0("y", 1:E)], DataM[, paste0("Z", 1:E)], DataM[, controls.info])
  }else{
    WeakData <- data.frame(DataM[, paste0("y", 1:E)], DataM[, paste0("Z", 1:E)])
  }

  # Workaround for labelling
  if(E == 1){
    if(controls.info[1] != "1"){
      colnames(WeakData) <- c("y1", "Z1", controls.info)
    }else{
      colnames(WeakData) <- c("y1", "Z1")
    }
  }

  # Set data to missing if indicator is equal to 2
  WeakData[DataM$Ind == 2, ] <- NA

  Method = "Proxy-IV"
  Obs <- data.frame(Tp = Te, Tc = Tn, To = To, Tt =Tt)


  # Compute variance-covariance matrix of residuals on event days, and impact matrix
  # needed for shock extraction
  # Compute VCOV on event days only, Note that this is correct, because residuals are set to missing before for non-event days
  Sig <- var(et, use = "complete.obs")

  if(sum(!is.na(vt)>0)){
    SigR <- var(vt, use = "complete.obs")
  }else{
    SigR <- NA
  }

  Psi <- irfest[1 , ,]

  return(list(irf = irfest, se = irfse,
              IVRes = IVRes,
              Obs = Obs, Method = Method, WeakData = WeakData,
              et = et, Sig = Sig, SigR = SigR, Psi = Psi))
}


estimLPHetFast <- function(y, O, Ind, Interact, E, P, cumRep, H, HStep, NormFac){
  # FAST VERSION NOT COMPUTING EVERYTHING FOR BOOTSTRAP
  
  # Estimates impulse responses with recursive heteroscedasticity-IV
  # Ind:    is an event indicator, with values
  #         0: Control day (no event) 
  #         1: Policy day (event) 
  #         2: Contaminated control day (other event)
  #         Impulse responses are only estimated based on policy and control days
  # y:      Outcome variables, for which the impulse responses are computed.
  #         Note that for each dimension the first variable's response is 
  #         normalized to unity. These are also the variables used for 
  #         instruments and recursive ordering
  # O:      Information set which can, but does not have to be, identical to y. 
  #         It will be included lagged 1...P
  # Interact: If TRUE, then O is interacted with regime dummies
  # P:      Maximum lag order for information set. If P = 0, no information set
  # H:      Maximum horizon up to which impulse response is computed
  # HStep:  If > 1, only every HStep impulse response is estimated 
  #         (starting at H = 1)
  # cumRep: Vector of T/F of the same dimension as y indicating whether to compute
  #         the cumulative impulse response for every variable included in y
  # E:      Number of dimensions of shocks to be identified via Cholesky 
  # NormFac:Normalizes the impulse response of the first variable to a specific value
  #         Set to 1 if no normalization is desired
  
  # Collect various properties of the data and observations to be used
  Nobs <- dim(y)[1]     # Number of observations in y
  beg  <- max(P+1)      # Start of the sample
  end  <- Nobs - H+1    # End of the sample
  N    <- dim(y)[2]     # Number of variables in y
  M    <- dim(O)[2]     # Number of variables in O
  
  if(sum(is.na(y))>0){
    warning("Missing values in y")
  }
  if(sum(is.na(O))>0){
    warning("Missing values in O")
  }
  if(sum(is.na(Ind))>0){
    warning("Missing values in Ind")
  }
  
  # Specify at which horizons IRF should be computed
  HSeries <- seq(1, H, HStep)
  HNum    <- length(HSeries)
  
  # Set up data set and and objects to save results
  DataM     <- data.frame(y, O, Ind)
  colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), "Ind")
  irfest    <- array(NA, dim = c(HNum, N, E))
  
  # Compute lagged variables included in information set O(t-1) 
  infoVars      <- c()
  if(P > 0){
    for(j in 1:M){
      for (p in 1:P){
        if(Interact == TRUE){
          for(e in 0:1){
            DataM[, paste0("Event", e, ".o", j, ".l", p)] <- dplyr::lag(DataM[, paste0("o", j)], p) * (Ind == e)
            infoVars <- c(infoVars, paste0("Event", e, ".o", j, ".l", p))
          }
        }else{
          DataM[, paste0("o", j, ".l", p)] <- dplyr::lag(DataM[, paste0("o", j)], p)
          infoVars <- c(infoVars, paste0("o", j, ".l", p))  
        }
        
      }
    }
  }
  if(P == 0){
    # Otherwise, only regress on a constant
    if(Interact == TRUE){
      for(e in 0:1){
        DataM[, paste0("Event", e)] <- (Ind == e)
        infoVars <- c(infoVars, paste0("Event", e))
      }
    }else{
      infoVars <- "1"      
    }
    
  }
  
  # Compute instrument separately for every dimension
  for(e in 1:E){
    
    # Set up the shock variable (instrumented variable), which is the e th variable
    # in y
    DataM$shockVar <- DataM[, paste0("y", e)]
    
    # Set up the instrument variable, which is also the e th variable in y
    DataM$Ze <- DataM[,paste0("y", e)]
    
    # Set up the control variables for recursive ordering
    # Note that we interact the dependent variable from previous
    # equation with the event dummy
    recVars <- c()
    recInst <- c()
    if(e>1){
      for(q in 1:(e-1)){
        recVars <- c(recVars, paste0("y", q))  
        recInst <- c(recInst, paste0("Z", q))  
      }
    }
    
    # Set up control variables for the local projections
    controls.info <- unique(c(infoVars))
    controls.info <- controls.info[controls.info !=""]
    
    controls.lp <- unique(c(infoVars, recVars))
    controls.lp <- controls.lp[controls.lp !=""]
    
    controls.iv <- unique(c(infoVars, recInst))
    controls.iv <- controls.iv[controls.iv !=""]
    
    # Account for missing data in instrument variable by setting it to an other event
    # which is not used in the estimation
    DataM$Ind[is.na(DataM$Ze)] = 2
    
    # Set up event days (Policy event, control, and other event)
    DataM$Event    <- (DataM$Ind == 1)
    DataM$NoEvent  <- (DataM$Ind == 0)
    DataM$OthEvent <- (DataM$Ind == 2)
    
    # Shorten data to subset of observations without missing values at beginning or end
    DataMSub <- DataM[beg:end, ] 
    
    
    # Orthogonalize the instrument if we include lags of dependent variable
    # See Lewis (2022)
    if(P > 0){
      
      #print(controls.info)
      
      myFormula = paste0("Ze ~ ", paste(controls.info, collapse = "+"))
      orthModel <- lm(as.formula(myFormula), data = subset(DataMSub, Ind < 2), na.action="na.exclude")
      
      #print(summary(orthModel))
      DataMSub$Ze <- residuals(orthModel, na.action="na.exclude")
      
    }
    
    
    # Compute instrument and F-Statistic (see Lewis, 2022, ReStat, and 
    # Rigobon, 2003, ReStat)
    # Compute number of events, control observations, and non-events
    Te <- sum(DataMSub$Event)
    Tn <- sum(DataMSub$NoEvent)
    To <- sum(DataMSub$OthEvent)
    Tt <- Te+Tn
    
    # Compute the instrument (see Lewis, 2022, ReStat and Rigobon and Sack, 2004, JME)
    DataMSub$Z <- (DataMSub$Event*Tt/Te - DataMSub$NoEvent*Tt/Tn)*DataMSub$Ze
    
    # Save instrument in original data for later use
    DataM[beg:end, paste0("Z", e)] <- DataMSub$Z
    DataM[beg:end, "Z"] <- DataMSub$Z
    
    # Estimate the impulse responses for every variable in y
    for (i in 1:N)
    {
      # Set up dependent variable
      DataM$depVar  <- DataM[, paste0("y", i)]
      cumRepi       <- cumRep[i]
      
      
      # Compute the impulse response for every horizon in a paralleled for loop
      # The reason why this is possible is that we directly estimate cumulative
      # impulse responses
      for(h in HSeries){             
        # Compute dependent variable (option for cumulative responses)
        # That is, we compute y(t+h) in order to compute the direct forecast later
        # (local projection)
       if(cumRepi == T){
          for(f in 1:h){
            if(f == 1){
              DataM$depVar.h <- dplyr::lead(DataM$depVar, f-1)  
            }else{
              DataM$depVar.h <- DataM$depVar.h + dplyr::lead(DataM$depVar, f-1)    
            }
          }
        }else{
          DataM$depVar.h <- dplyr::lead(DataM$depVar, h-1)  
        }
        
        # Shorten data to subset which contains no missing values
        DataMSub <- DataM[beg:end, ]
        
        # LP (Jorda, 2005), dependent variable is just y(t+h)
        # Control variables for direct forecast y(t+h|t-1)
        # Compute direct forecast t+h|t-1, where shockVar is instrumented
        #print(controls.lp)
        #print(controls.iv)
        
        myFormula <- paste0("depVar.h ~ shockVar + ", 
                            paste(controls.lp, collapse = "+"),
                            paste0("| Z + ", 
                                   paste(controls.iv, collapse = "+")))
        
        IV.mod  <- ivreg(as.formula(myFormula), data = subset(DataMSub, Ind < 2))
        #IV.se   <- sqrt(diag(NeweyWest(IV.mod, prewhite = FALSE)))
        
        #  Result <- list(IV.mod)
        # }
        
        
        # Normalizes IRFs to a specific value. Because initial response
        # normalized to unity, just multiply by NormFac
        irfest[h, i, e] <- IV.mod$coefficients["shockVar"]*NormFac
        
      }
      
      # Label the rows of the impulse responses to match those that have been
      # estimated and to start with 0 (immediate response)
      dimnames(irfest)[[1]] <- HSeries-1
      
    }
  }
  
  return(list(irf = irfest))
  
}


estimLPProxyFast <- function(y, O, Z, recursive = F, Ind, E, P, cumRep, H, HStep, NormFac){
  # FAST VERSION NOT COMPUTING EVERYTHING FOR BOOTSTRAP
  
  # Estimates impulse responses with recursive or non-recursive proxy-LP
  # Ind:    is an event indicator, with values
  #         0: Control day (no event) 
  #         1: Policy day (event) 
  #         2: Contaminated control day (other event)
  #         Impulse responses are only estimated based on policy (and control) days
  # y:      Outcome variables, for which the impulse responses are computed.
  #         Note that for each dimension the first variable's response is 
  #         normalized to unity. These are also the variables used for 
  #         instruments and recursive ordering
  # Z:      Is a matrix of instruments for the first E endogenous variables in y
  # recursive T/F whether we additionally impose the recursive zero restrictions
  # O:      Information set which can, but does not have to be, identical to y. 
  #         It will be included lagged 1...P
  # P:      Maximum lag order for information set. If P = 0, no information set
  # H:      Maximum horizon up to which impulse response is computed
  # HStep:  If > 1, only every HStep impulse response is estimated 
  #         (starting at H = 1)
  # cumRep: Vector of T/F of the same dimension as y indicating whether to compute
  #         the cumulative impulse response for every variable included in y
  # E:      Number of dimensions of shocks to be identified via Cholesky 
  # NormFac:Normalizes the impulse response of the first variable to a specific value
  #         Set to 1 if no normalization is desired
  
  # Collect various properties of the data and observations to be used
  Nobs <- dim(y)[1]     # Number of observations in y
  beg  <- max(P+1)      # Start of the sample
  end  <- Nobs - H+1    # End of the sample
  N    <- dim(y)[2]     # Number of variables in y
  M    <- dim(O)[2]     # Number of variables in O
  
  if(sum(is.na(y))>0){
    warning("Missing values in y")
  }
  if(sum(is.na(O))>0){
    warning("Missing values in O")
  }
  if(sum(is.na(Ind))>0){
    warning("Missing values in Ind")
  }
  
  # Specify at which horizons IRF should be computed
  HSeries <- seq(1, H, HStep)
  HNum    <- length(HSeries)
  
  # Set up data set and and objects to save results
  DataM     <- data.frame(y, O, Ind, Z)
  colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), "Ind", paste0("z", 1:E))
  irfest    <- array(NA, dim = c(HNum, N, E))
  
  # Compute lagged variables included in information set O(t-1) 
  infoVars      <- c()
  if(P > 0){
    for(j in 1:M){
      for (p in 1:P){
        DataM[, paste0("o", j, ".l", p)] <- dplyr::lag(DataM[, paste0("o", j)], p)
        infoVars <- c(infoVars, paste0("o", j, ".l", p))
      }
    }
  }
  if(P == 0){
    # Otherwise, only regress on a constant
    infoVars <- "1"    
  }
  
  # Compute instrument separately for every dimension
  for(e in 1:E){
    
    # Set up the shock variable (instrumented variable), which is the e th variable
    # in y
    DataM$shockVar <- DataM[, paste0("y", e)]
    
    # Set up the instrument variable, which is the e th variable in z
    DataM$Ze <- DataM[,paste0("z", e)]
    
    # Set up the control variables for recursive ordering
    # Note that we interact the dependent variable from the previous
    # equation with the event dummy
    recVars <- c()
    recInst <- c()
    if(recursive == T){
      if(e>1){
        for(q in 1:(e-1)){
          recVars <- c(recVars, paste0("y", q))  
          recInst <- c(recInst, paste0("Z", q))  
        }
      }
    }
    
    # Set up control variables for the local projections
    controls.info <- unique(c(infoVars))
    controls.info <- controls.info[controls.info !=""]
    
    controls.lp <- unique(c(infoVars, recVars))
    controls.lp <- controls.lp[controls.lp !=""]
    
    controls.iv <- unique(c(infoVars, recInst))
    controls.iv <- controls.iv[controls.iv !=""]
    
    # Account for missing data in instrument variable by setting it to an other event
    # which is not used in the estimation
    # If instrument is missing, we can still predict the unobserved shocks for these periods!
    # So if there is an event, do not set to 2, only if there is no event
    # TODO: Check whether this really necessary. Why not include in estimation?
    #DataM$Ind[is.na(DataM$Ze & DataM$Ind == 0)] = 2
    
    # Set up event days (Policy event, control, and other event)
    DataM$Event    <- (DataM$Ind == 1)
    DataM$NoEvent  <- (DataM$Ind == 0)
    DataM$OthEvent <- (DataM$Ind == 2)
    
    # Shorten data to subset of observations without missing values
    DataMSub <- DataM[beg:end, ] 
    
    # Collect some statistics about instrument and save instrument for later use
    Te <- sum(DataMSub$Event)
    Tn <- sum(DataMSub$NoEvent)
    To <- sum(DataMSub$OthEvent)
    Tt <- Te+Tn
    
    # Get the instrument for shock e
    DataMSub$Z <- DataMSub$Ze
    
    # Save instrument in original data for later use
    DataM[beg:end, paste0("Z", e)] <- DataMSub$Z
    DataM[beg:end, "Z"] <- DataMSub$Z
    
    # Estimate the impulse responses for every variable in y
    for (i in 1:N)
    {
      # Set up dependent variable
      DataM$depVar  <- DataM[, paste0("y", i)]
      cumRepi       <- cumRep[i]
      
      # Compute the impulse response for every horizon in a paralleled for loop
       for(h in HSeries){
          if(cumRepi == T){
          for(f in 1:h){
            if(f == 1){
              DataM$depVar.h <- dplyr::lead(DataM$depVar, f-1)  
            }else{
              DataM$depVar.h <- DataM$depVar.h + dplyr::lead(DataM$depVar, f-1)    
            }
          }
        }else{
          DataM$depVar.h <- dplyr::lead(DataM$depVar, h-1)  
        }
        
        # Shorten data to subset which contains no missing values
        DataMSub <- DataM[beg:end, ]
        
        # LP (Jorda, 2005), dependent variable is just y(t+h)
        # Control variables for direct forecast y(t+h|t-1)
        # Compute direct forecast t+h|t-1, where shockVar is instrumented
        myFormula <- paste0("depVar.h ~ shockVar + ", 
                            paste(controls.lp, collapse = "+"),
                            paste0("| Z + ", 
                                   paste(controls.iv, collapse = "+")))
        
        IV.mod  <- ivreg(as.formula(myFormula), data = subset(DataMSub, Ind < 2))
        
        #             Result <- list(IV.mod)
        #           }
        
        # Save results from parallel computing for later use
        #for(h in 1:length(ResultsIRF)){
        #  IV.mod  <- ResultsIRF[[h]][[1]]
        
        # Normalizes IRFs to a specific value. Because initial response
        # normalized to unity, just multiply by NormFac
        irfest[h, i, e] <- IV.mod$coefficients["shockVar"]*NormFac
        
      }
      
      # Label the rows of the impulse responses to match those that have been
      # estimated and to start with 0 (immediate response)
      dimnames(irfest)[[1]] <- HSeries-1
      
    }
  }
  
  
  return(list(irf = irfest))
}





# Extract shocks using Kalman-filter
extractShocksKF <- function(Sig, SigR, Psi, et, tol, scale = TRUE){

  # Set tolerance for generalized inverse (maximum of specified value and the default value for ginv)
  tol <- max(sqrt(.Machine$double.eps), tol)

  if(is.na(tol)){
    tol <- sqrt(.Machine$double.eps)
  }

  eps <- as.matrix(et[, 1:dim(Psi)[2]])
  eps[, ] <- NA

  if(scale == TRUE){
    # Note that we have a unit impact normalization on Psi. SO we have to back
    # out the implied scale of the shocks, so that we can remormalize the shocks
    # to have a unit variance.
    # For this, we need R, the variance of shocks without a policy announcement
    q      <- vec(Sig - SigR)
    A      <- sapply(1:ncol(Psi), function(i) c(Psi[, i] %*% t(Psi[, i])))
    sig    <- c(ginv(A)%*%q)

    if(length(sig)>1){
      SigEps <- diag(sig)
      myScale <- solve(diag(sqrt(1/sig)))

      if(any(is.na(myScale))){
        myScale <- diag(abs(sig))
      }

    }else{
      SigEps  <- sig
      myScale <- 1/sqrt(1/sig)

      if(any(is.na(myScale))){
        myScale <- abs(sig)
      }
    }

  }else{
    myScale = 1
  }


  # Kalman filter formula:
  for(t in 1:dim(eps)[1]){
    eps[t, ] <- myScale%*%t(Psi)%*%ginv(Sig, tol = tol)%*%as.matrix(et[t, ])
  }

  return(eps)
}


plotPvals <- function(pvals, HTick, Labels, shockLabels){
  # Plot IRFS: First is the comparison response, second
  # point estimate, third standard error of point estimate

  myGraphs <- list()
  noDims   <- length(dim(pvals))
  HNum     <- dim(pvals)[1]
  HSeries  <- as.numeric(dimnames(pvals)[[1]])

  N <- dim(pvals)[2]
  n = 1
  if(noDims == 3){
    R <- dim(pvals)[3]
  }else{
    R <- 1
    dim(pvals) <- c(HNum, N, 1)
  }

  if(is.null(shockLabels)){
    shockLabels <- paste0("shock ", 1:R)
  }

  for(j in 1:R){
    for(i in 1:N){

      myIRF <- data.frame(HSeries, pvals[, i, j])
      colnames(myIRF) <- c("Horizon", "pval")

      g1 <- ggplot(myIRF, aes(Horizon))

      g1 <- g1 + theme_minimal() + xlab("") + ggtitle(paste0(Labels[i]))+theme(plot.title = element_text(size = 10))+
        ylab("")+scale_x_continuous(breaks=seq(HSeries[1], max(HSeries), HTick))
      #g1 <- g1 + geom_hline(yintercept = 0, size=0.2) + theme(legend.position = "none")

      g1 <- g1 + geom_line(aes(y = `pval`), colour = "steelblue", size = .7)

      g1 <- g1 + geom_hline(yintercept = 0.05, colour = "black", linetype = "solid", size=.4) + theme(legend.position = "none")
      g1 <- g1 + geom_hline(yintercept = 0.1, colour = "black", linetype = "dotted", size=.4) + theme(legend.position = "none")

      g1 <- g1 +  theme(panel.grid = element_line(color = "gray",
                                                  size = 0.2,
                                                  linetype = "dotted"),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(color = "black",
                                                    fill = NA,
                                                    size = 0.2,
                                                    linetype = "solid"))

      g1 <- g1 + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

      myGraphs[[n]] <- g1
      n = n+1
    }
  }
  return(myGraphs)
}


plot2IRFs <- function(IRF1, IRF1se, IRF2, IRF2se, HTick, Labels, shockLabels){
  # Plot IRFS: First is the comparison response, second
  # point estimate, third standard error of point estimate

  myGraphs <- list()
  noDims   <- length(dim(IRF1))
  HNum     <- dim(IRF1)[1]
  HSeries  <- as.numeric(dimnames(IRF1)[[1]])

  N <- dim(IRF1)[2]
  n = 1
  if(noDims == 3){
    R <- dim(IRF1)[3]
  }else{
    R <- 1
    dim(IRF1) <- c(HNum, N, 1)
    dim(IRF2) <- c(HNum, N, 1)
  }

  if(is.null(shockLabels)){
    shockLabels <- paste0("shock ", 1:R)
  }



  # Compute confidence intervals
  upper1_90 <- array(NA, dim = c(HNum, N, R))
  lower1_90 <- array(NA, dim = c(HNum, N, R))
  upper2_90 <- array(NA, dim = c(HNum, N, R))
  lower2_90 <- array(NA, dim = c(HNum, N, R))
  for(j in 1:R){
    for(i in 1:N){
      upper1_90[, i, j] <- IRF1[, i, j]+1.64*IRF1se[, i, j]
      lower1_90[, i, j] <- IRF1[, i, j]-1.64*IRF1se[, i, j]

      upper2_90[, i, j] <- IRF2[, i, j]+1.64*IRF2se[, i, j]
      lower2_90[, i, j] <- IRF2[, i, j]-1.64*IRF2se[, i, j]
    }
  }


  for(j in 1:R){
    for(i in 1:N){

      myIRF <- data.frame(HSeries, IRF1[, i, j], IRF2[, i, j], upper1_90[, i, j],
                            lower1_90[, i, j], upper2_90[, i, j], lower2_90[, i, j])
      colnames(myIRF) <- c("Horizon", "IRF1", "IRF2", "upper1_90", "lower1_90",
                             "upper2_90", "lower2_90")

      g1 <- ggplot(myIRF, aes(Horizon))

      g1 <- g1 + theme_minimal() + xlab("") + ggtitle(paste0(Labels[i]))+theme(plot.title = element_text(size = 10))+
        ylab("")+scale_x_continuous(breaks=seq(HSeries[1], max(HSeries), HTick))
      g1 <- g1 + geom_hline(yintercept = 0, size=0.2) + theme(legend.position = "none")

      g1 <- g1 + geom_line(aes(y = `IRF1`), colour = "steelblue", size = .7)
      g1 <- g1 + geom_line(aes(y = `IRF2`), size = 0.7, color = "darkred", linetype = "dotted")

      g1 <- g1 + geom_ribbon(aes_string(ymin = "lower1_90" , ymax = "upper1_90"), fill= "steelblue", alpha=1/10)
      g1 <- g1 + geom_ribbon(aes_string(ymin = "lower2_90" , ymax = "upper2_90"), fill= "darkred", alpha=1/10)

      g1 <- g1 +  theme(panel.grid = element_line(color = "gray",
                                                  size = 0.2,
                                                  linetype = "dotted"),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(color = "black",
                                                    fill = NA,
                                                    size = 0.2,
                                                    linetype = "solid"))

      g1 <- g1 + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

      myGraphs[[n]] <- g1
      n = n+1
    }
  }
  return(myGraphs)
}


plotIRFs <- function(IRFs, IRFest, IRFse, HTick, Labels, shockLabels){
  # Plot IRFS: First is the comparison response, second
  # point estimate, third standard error of point estimate

  # Workaround, if no "true" IRF for comparison, just set to missing
  if(is.array(IRFs) == FALSE){
    if(is.na(IRFs)){
      IRFs       <- IRFest
      IRFs[, , ] <- NaN
    }
  }
  myGraphs <- list()
  noDims   <- length(dim(IRFs))
  HNum     <- dim(IRFs)[1]
  HSeries  <- as.numeric(dimnames(IRFs)[[1]])

  N <- dim(IRFs)[2]
  n = 1
  if(noDims == 3){
    R <- dim(IRFs)[3]
  }else{
    R <- 1
    dim(IRFs) <- c(HNum, N, 1)
    dim(IRFest) <- c(HNum, N, 1)
  }

  if(is.null(shockLabels)){
    shockLabels <- paste0("shock ", 1:R)
  }


  if(is.array(IRFse) == TRUE){

    # Compute confidence intervals
    upper95 <- array(NA, dim = c(HNum, N, R))
    lower95 <- array(NA, dim = c(HNum, N, R))
    upper90 <- array(NA, dim = c(HNum, N, R))
    lower90 <- array(NA, dim = c(HNum, N, R))
    for(j in 1:R){
      for(i in 1:N){
        upper95[, i, j] <- IRFest[, i, j]+1.96*IRFse[, i, j]
        lower95[, i, j] <- IRFest[, i, j]-1.96*IRFse[, i, j]
        upper90[, i, j] <- IRFest[, i, j]+1.64*IRFse[, i, j]
        lower90[, i, j] <- IRFest[, i, j]-1.64*IRFse[, i, j]
      }
    }
  }

  for(j in 1:R){
    for(i in 1:N){

      if((is.array(IRFest) == TRUE & is.array(IRFse) == FALSE)){
        myIRF <- data.frame(HSeries, IRFs[, i, j], IRFest[, i, j])
        colnames(myIRF) <- c("Horizon", "Actual", "Estimate")
      }
      if((is.array(IRFest) == TRUE & is.array(IRFse) == TRUE)){
        myIRF <- data.frame(HSeries, IRFs[, i, j], IRFest[, i, j], upper95[, i, j],
                            lower95[, i, j], upper90[, i, j], lower90[, i, j])
        colnames(myIRF) <- c("Horizon", "Actual", "Estimate", "upper95", "lower95",
                             "upper90", "lower90")
      }
      if((is.array(IRFest) == FALSE & is.array(IRFse) == FALSE)){
        myIRF <- data.frame(HSeries, IRFs[, i, j])
        colnames(myIRF) <- c("Horizon", "Actual")
      }
      g1 <- ggplot(myIRF, aes(Horizon))

      g1 <- g1 + theme_minimal() + xlab("") + ggtitle(paste0(Labels[i]))+theme(plot.title = element_text(size = 10))+
                 ylab("")+scale_x_continuous(breaks=seq(HSeries[1], max(HSeries), HTick))
      g1 <- g1 + geom_hline(yintercept = 0, size=0.2) + theme(legend.position = "none")
      if(is.array(IRFest) == TRUE){
        g1 <- g1 + geom_line(aes(y = `Estimate`), colour = "steelblue", size = .7)
        g1 <- g1 + geom_line(aes(y = `Actual`), size = 0.7, color = "darkred", linetype = "dotted")
      }
      if(is.array(IRFse) == TRUE){
        g1 <- g1 + geom_ribbon(aes_string(ymin = "lower95" , ymax = "upper95"), fill= "steelblue", alpha=1/10)
        g1 <- g1 + geom_ribbon(aes_string(ymin = "lower90" , ymax = "upper90"), fill= "steelblue", alpha=2/10)
      }
      if((is.array(IRFest) == FALSE & is.array(IRFse) == FALSE)){
        g1 <- g1 + geom_line(aes(y = `Actual`), size = 0.7, colour = "darkred")
      }

      g1 <- g1 +  theme(panel.grid = element_line(color = "gray",
                                                  size = 0.2,
                                                  linetype = "dotted"),
                        panel.grid.minor = element_blank(),
                        panel.border = element_rect(color = "black",
                                                    fill = NA,
                                                    size = 0.2,
                                                    linetype = "solid"))

      g1 <- g1 + theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

      myGraphs[[n]] <- g1
      n = n+1
    }
  }
  return(myGraphs)
}
