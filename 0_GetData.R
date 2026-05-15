#-------------------------------------------------------------------------------
# 0_GetData.R
# Multiple monetary policy shocks from daily data: A heteroskedasticity IV approach
# by Marc Burri and Daniel Kaufmann
#-------------------------------------------------------------------------------
# OPTIONAL: Constructs data set
#
# Input:   FRED, Yahoo Finance, Excel files       
# Output: ./Data/Data.RData
#
# Run order: this script first, then 1_MultiDimShokcs.R
#-------------------------------------------------------------------------------

# In the first instance, install the package iv from github repository
# install.packages("devtools")
# library(devtools)
# remove.packages("hetiv")
# install("C:\\Users\\daenu\\GitHub\\hetiv")

# Load user-defined commands and packages
library("hetiv")
library("tsbox")
library("xts")
library("dplyr")
library("openxlsx")
library("fredr")
library("quantmod")

# Personal FRED API key (get your own key at https://fred.stlouisfed.org/docs/api/api_key.html)
fredr_set_key("YOUR_API_KEY_HERE")

# Sample settings
myStart <- "1988-01-01"
myEnd   <- "2025-12-31"
fillGap <- 4

#-------------------------------------------------------------------------------
# 1) Get event data from Bauer and Swanson (2022) and Acosta et al. (2025)
#-------------------------------------------------------------------------------
Date <- as.Date((Date = seq(as.Date(myStart), as.Date(myEnd), by = "days")))

BS <- read.xlsx("./Data/FOMC_Bauer_Swanson.xlsx", sheet = "High_Freq", detectDates = TRUE)
AC <- read.xlsx("./Data/USMP_database.xlsx", sheet = "Monetary Events", detectDates = TRUE)

BS <- data.frame(Date = as.Date(BS$Date))
AC <- data.frame(Date = as.Date(AC$Date))

# Create an Event Indicator, which is equal to 1 if there is an event in either
# of the two data sets and 0 otherwise
Events <- data.frame(Date = Date,
                     Ind  = Date %in% AC$Date | Date %in% BS$Date)
IndE    <- xts(as.numeric(Events$Ind), order.by = Events$Date)

#-------------------------------------------------------------------------------
# 2) Get multi-dimensional shocks from Swanson (2021)
#-------------------------------------------------------------------------------
SW <- read.xlsx("./Data/SwansonShocks.xlsx", sheet = "Data", detectDates = TRUE, startRow = 2)
SW <- data.frame(Date = as.Date(SW[, 1]), 
                 SFFR = SW[, 2],
                 SFG  = SW[, 3])
SW    <- xts(SW[, 2:3], order.by = SW[, 1])


#-------------------------------------------------------------------------------
# 3) Get financial market data from FRED
#-------------------------------------------------------------------------------
# Interest rate data
FFR      <- fredr("DFF")
IR3Mfed  <- fredr("DTB3")
IR6Mfed  <- fredr("DTB6")
IR1Yfed  <- fredr("DGS1")
IR2Yfed  <- fredr("DGS2")
IR3Yfed  <- fredr("DGS3")
IR5Yfed  <- fredr("DGS5")
IR7Yfed  <- fredr("DGS7")
IR10Yfed <- fredr("DGS10")
IR30Yfed <- fredr("DGS30")
             
# Create TS
FFR      <- xts(FFR$value, FFR$date)
IR3Mfed  <- xts(IR3Mfed$value, IR3Mfed$date)
IR6Mfed  <- xts(IR6Mfed$value, IR6Mfed$date)
IR1Yfed  <- xts(IR1Yfed$value, IR1Yfed$date)
IR2Yfed  <- xts(IR2Yfed$value, IR2Yfed$date)
IR3Yfed  <- xts(IR3Yfed$value, IR3Yfed$date)
IR5Yfed  <- xts(IR5Yfed$value, IR5Yfed$date)
IR7Yfed  <- xts(IR7Yfed$value, IR7Yfed$date)
IR10Yfed <- xts(IR10Yfed$value, IR10Yfed$date)
IR30Yfed <- xts(IR30Yfed$value, IR30Yfed$date)

# Create averages
IRSTfed <- (FFR + IR3Mfed + IR6Mfed)/3
IRMTfed <- (IR2Yfed + IR3Yfed + IR5Yfed)/3
IRLTfed <- (IR7Yfed + IR10Yfed + IR30Yfed)/3

# NEER data
NEER1 <- fredr("DTWEXM")
NEER2 <- fredr("DTWEXAFEGS")

# Create TS
NEER1 <- xts(NEER1$value, NEER1$date)
NEER2 <- xts(NEER2$value, NEER2$date)

# Link series and transform to units of US dollars per unit of foreign currency
linkDate <- ts_summary(NEER2)$start
NEER     <- 100*1/ts_bind(
                    ts_index(ts_span(NEER1, end = linkDate), linkDate),
                    ts_index(NEER2, linkDate)
                    )

# Stock price data
NASDAQ     <- fredr("NASDAQCOM")
NASDAQ100  <- fredr("NASDAQ100")
getSymbols("^GSPC", src = "yahoo", from = "1985-01-01")
SP500 <- GSPC[, "GSPC.Close"] 

# Create TS
NASDAQ    <- xts(NASDAQ$value, NASDAQ$date)
NASDAQ100 <- xts(NASDAQ100$value, NASDAQ100$date)
SP500     <- xts(SP500, index(SP500))

# Create weighted average
 Stocks <- (0.5*ts_index(SP500, "2000-12-01") + 
            0.25*ts_index(NASDAQ, "2000-12-01") + 
            0.25*ts_index(NASDAQ100, "2000-12-01"))

# Corporate bond spread data
SpreadBAA <- fredr("BAA10Y")
SpreadAAA <- fredr("AAA10Y")

# Create TS
SpreadBAA <- xts(SpreadBAA$value, SpreadBAA$date)
SpreadAAA <- xts(SpreadAAA$value, SpreadAAA$date)

# Create average
Spread <- (SpreadBAA + SpreadAAA)/2

# Get VIX data
VIX1 <- fredr("VXOCLS")
VIX2 <- fredr("VIXCLS")

# Create TS
VIX1 <- xts(VIX1$value, VIX1$date)
VIX2 <- xts(VIX2$value, VIX2$date)

# Link the series
linkDate <- ts_summary(VIX2)$start
VIX      <- ts_bind(ts_span(VIX1, end = linkDate),
                   VIX2)

#-------------------------------------------------------------------------------
# 4) Get treasury yield skewness from Bauer and Chernov (2024)
#-------------------------------------------------------------------------------
TRSkew <- read.xlsx("./Data/treasury-yield-skewness-data.xlsx", sheet = "ISK", detectDates = TRUE)
TRSkew <- xts(TRSkew$isk, order.by = TRSkew$date)

#-------------------------------------------------------------------------------
# 5) Create joint data set and do sampling decisions
#-------------------------------------------------------------------------------
yData <- ts_c(IndE, FFR, IR3Mfed, IR6Mfed, IR1Yfed, IR2Yfed, IR3Yfed, IR5Yfed, IR7Yfed, IR10Yfed, IR30Yfed,
              IRSTfed, IRMTfed, IRLTfed,
              NEER, Stocks, Spread, VIX, TRSkew, SW)
yData <- data.frame(date = index(yData), yData)

# Remove weekends and holidays
yData <- yData %>%
  filter(
    !weekdays(date) %in% c("Saturday", "Sunday") & !date %in% 
    as.Date(c(timeDate::holidayNYSE(1980:2026))) # Remove weekend and holidays)
) 

# Fill up to 4 missing values by linear interpolation
yData[, !colnames(yData) %in% c("date", "IndE")] <- filllinear(yData[,!colnames(yData) %in% c("date", "IndE")], fillGap)

# Transform data taking differences with respect to last working day
# Transform all variables so that they are stationary and measured in %/pp
for(xr in unique(c("NEER", "Stocks"))){
  yData[, xr] <- 100*logdiff(yData[, xr])
}

# Transform all variables so that they are stationary and measured in %/pp
for(ir in unique(c("FFR",  "IR3Mfed", "IR6Mfed", "IR1Yfed", 
                   "IR2Yfed", "IR3Yfed", "IR5Yfed", "IR7Yfed", "IR10Yfed", "IR30Yfed", 
                   "IRLTfed", "IRMTfed", "IRSTfed",
                   "Spread","VIX"))){
  yData[, ir] <- firstdiff(yData[, ir])
}

# Shorten to sample
Data <- subset(yData, date >= as.Date(myStart) & date <= as.Date(myEnd))

# Export as time series object without date vector
Data <- as.xts(Data[, -1])
save(Data, file="./Data/Data.RData")


#-------------------------------------------------------------------------------
# End of code
#-------------------------------------------------------------------------------



