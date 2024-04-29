##### Remote Sensing Final Project #####

library(lubridate)
library(dplyr)
library(hydroGOF)
library(wesanderson)

##### Lake Mead G-REALM remote sensing data (Birkett et al., 2017) #####

# Import ENVISAT data
ENVItbl <- read.delim("filename",header=FALSE,sep="",skip=49)
colnames(ENVItbl) <- c('Satellite mission name',
                             'Satellite repeat cycle',
                             'Calendar year/month/day of along track observations traversing target',
                             'Hour of day at mid point of along track pass traversing target',
                             'Minutes of hour at mid point of along track pass traversing target',
                             'Target height variation with respect to ENVISAT reference pass mean level (meters, default=999.99)',
                             'Estimated error of target height variation with respect to reference mean level (meters, default=99.999)',
                             'Mean along track Ku-band backscatter coefficient (decibels, default=999.99)',
                             'Wet tropospheric correction applied to range observation (RAD=ENVISAT radiometer, ECM=ECMWF Operational model, MIX=combination, U/A=unavailable, N/A=not applicable)',
                             'Ionosphere correction applied to range observation (GIM=GPS model, U/A=unavailable, N/A=not applicable)',
                             'Dry tropospheric correction applied to range observation (ECM=ECMWF Operational model, U/A=unavailable, N/A=not applicable)',
                             'Instrument operating mode 1 (default=9)',
                             'Instrument operating mode 2 (default=9)',
                             'Flag for potential frozen surface (ice-on=1, ice-off or unknown=0)',
                             'Target height variation in EGM2008 datum (meters above mean sea level, default=9999.99)',
                             'Flag for ENVISAT data source (GDR=0)'
                             )
ENVItbl <- subset(ENVItbl, `Target height variation in EGM2008 datum (meters above mean sea level, default=9999.99)` !=9999.99)
ENVItbl$Date <- ymd(ENVItbl$`Calendar year/month/day of along track observations traversing target`)
ENVItbl$Elevation_m <- ENVItbl$`Target height variation in EGM2008 datum (meters above mean sea level, default=9999.99)`
ENVItbl <- subset(ENVItbl, select=c(Date,Elevation_m))

# Import Sentinel-3A data
S3Atbl <- read.delim("filename",header=FALSE,sep="",skip=49)
colnames(S3Atbl) <- c('Satellite mission name',
                       'Satellite repeat cycle',
                       'Calendar year/month/day of along track observations traversing target',
                       'Hour of day at mid point of along track pass traversing target',
                       'Minutes of hour at mid point of along track pass traversing target',
                       'Target height variation with respect to ENVISAT reference pass mean level (meters, default=999.99)',
                       'Estimated error of target height variation with respect to reference mean level (meters, default=99.999)',
                       'Mean along track Ku-band backscatter coefficient (decibels, default=999.99)',
                       'Wet tropospheric correction applied to range observation (RAD=ENVISAT radiometer, ECM=ECMWF Operational model, MIX=combination, U/A=unavailable, N/A=not applicable)',
                       'Ionosphere correction applied to range observation (GIM=GPS model, U/A=unavailable, N/A=not applicable)',
                       'Dry tropospheric correction applied to range observation (ECM=ECMWF Operational model, U/A=unavailable, N/A=not applicable)',
                       'Instrument operating mode 1 (default=9)',
                       'Instrument operating mode 2 (default=9)',
                       'Flag for potential frozen surface (ice-on=1, ice-off or unknown=0)',
                       'Target height variation in EGM2008 datum (meters above mean sea level, default=9999.99)',
                       'Flag for Sen-3A data source (NTC-R=0, NTC=1, STC=2)'
                       )
S3Atbl <- subset(S3Atbl, `Target height variation in EGM2008 datum (meters above mean sea level, default=9999.99)` !=9999.99)
S3Atbl$Date <- ymd(S3Atbl$`Calendar year/month/day of along track observations traversing target`)
S3Atbl$Elevation_m <- S3Atbl$`Target height variation in EGM2008 datum (meters above mean sea level, default=9999.99)`
S3Atbl <- subset(S3Atbl, select=c(Date,Elevation_m))

# Convert masl to m depth
damheight <- 223 #(m)
damelev <- 376 #top of dam elevation (masl)
baseelev <- damelev-damheight
ENVItbl$Depth_m <- ENVItbl$Elevation_m-baseelev
S3Atbl$Depth_m <- S3Atbl$Elevation_m-baseelev

##### Lake Mead ResOpsUS data (Steyaert et al., 2022) #####

ResOpsUSdata <- read.csv("filename",header=TRUE)
ResOpsUSdata$Date <- ymd(ResOpsUSdata$date)
ResOpsUSdata$ROU_Elevation_m <- ResOpsUSdata$elevation
ResOpsUSdata$ROU_Storage_mcm <- ResOpsUSdata$storage
ResOpsUSdata <- subset(ResOpsUSdata, select=c(Date,ROU_Elevation_m,ROU_Storage_mcm))
ResOpsUSdata <- na.omit(ResOpsUSdata)

# Add ResOpsUS data to RS tables
ENVItbl <- merge(ENVItbl,ResOpsUSdata,by="Date",all.x=TRUE)
S3Atbl <- merge(S3Atbl,ResOpsUSdata,by="Date",all.x=TRUE)

##### Geometric Storage-Depth lookup table (Yigzaw et al., 2018) #####

Y_SDtbl <- read.csv('filename',skip=7)
colnames(Y_SDtbl) <- c('Depth_m','Area_skm','Storage_mcm')

# Round all depth values
Y_SDtbl$Depth_m <- round(Y_SDtbl$Depth_m, digits = 1)
ENVItbl$Depth_m <- round(ENVItbl$Depth_m, digits = 1)
S3Atbl$Depth_m <- round(S3Atbl$Depth_m, digits = 1)

# Add index column to time series tables
ENVItbl$index <- seq(1, length(ENVItbl$Depth_m))
S3Atbl$index <- seq(1, length(S3Atbl$Depth_m))

ENVItbl <- merge(ENVItbl,Y_SDtbl,by="Depth_m",all.x=TRUE)
ENVItbl <- ENVItbl[order(ENVItbl$index),]
ENVItbl$Y_EstStorage_mcm <- ENVItbl$Storage_mcm
ENVItbl <- subset(ENVItbl, select=-c(Storage_mcm,Area_skm))

S3Atbl <- merge(S3Atbl,Y_SDtbl,by="Depth_m",all.x=TRUE)
S3Atbl <- S3Atbl[order(S3Atbl$index),]
S3Atbl$Y_EstStorage_mcm <- S3Atbl$Storage_mcm
S3Atbl <- subset(S3Atbl, select=-c(Storage_mcm,Area_skm))

##### GRDL storage-depth lookup table (Hao et al., 2024) #####

GRDL_SDtbl <- read.csv('filename', skip=4,sep=';')
colnames(GRDL_SDtbl) <- c('Depth_m','Area_skm','Storage_mcm')

# GRDL lookup table does not include the range of values from G-REALM data

# GRDL_SDtbl$Depth_m <- round(GRDL_SDtbl$Depth_m)
# ENVItbl$Depth_m <- round(ENVItbl$Depth_m)
# S3Atbl$Depth_m <- round(S3Atbl$Depth_m)
# 
# ENVImergedGRDLSD <- merge(data.frame(Depth_m = ENVItbl$Depth_m,
#                                      index = ENVItbl$index),
#                           GRDL_SDtbl, by = "Depth_m", all.x = TRUE)
# ENVImergedGRDLSD <- ENVImergedGRDLSD[order(ENVImergedGRDLSD$index),]
# ENVItbl$GRDL_EstStorage_mcm <- ENVImergedGRDLSD$Storage_mcm
# 
# S3AmergedGRDLSD <- merge(data.frame(Depth_m = S3Atbl$Depth_m,
#                                     index = S3Atbl$index),
#                          GRDL_SDtbl, by = "Depth_m", all.x = TRUE)
# S3AmergedGRDLSD <- S3AmergedGRDLSD[order(S3AmergedGRDLSD$index),]
# S3Atbl$GRDL_EstStorage_mcm <- S3AmergedGRDLSD$Storage_mcm

##### My method of making S-E tables from ResOpsUS data #####

ENVIenddate <- max(ENVItbl$Date)
S3Astartdate <- min(S3Atbl$Date)

ENVIstartdate <- min(ENVItbl$Date)
S3Aenddate <- max(S3Atbl$Date)

my_SE <- data.frame(Date=ResOpsUSdata$Date,
                    Storage_mcm=ResOpsUSdata$ROU_Storage_mcm,
                    Elevation_m=ResOpsUSdata$ROU_Elevation_m)
my_SE <- na.omit(my_SE)

# Subset to only values between ENVI and S3A missions
my_SE <- my_SE[my_SE$Date>=ENVIenddate & my_SE$Date<=S3Astartdate,]
my_SE <- subset(my_SE,select=-c(Date))

my_SE <- my_SE[order(my_SE$Elevation_m),]

repeat {
  lengthbefore <- nrow(my_SE)
  my_SE <- my_SE[my_SE$Storage_mcm > lag(my_SE$Storage_mcm) &
                           my_SE$Storage_mcm < lead(my_SE$Storage_mcm),]
  lengthafter <- nrow(my_SE)
  if (lengthbefore==lengthafter){
    break
  }
}

my_SE <- na.omit(my_SE)

my_SE$Elevation_m <- round(my_SE$Elevation_m,digits=1)

my_SE <- aggregate(Storage_mcm ~ Elevation_m, data = my_SE, FUN = mean)
allelev <- seq(min(my_SE$Elevation_m), max(my_SE$Elevation_m), 0.1)
interpolated <- approx(x=my_SE$Elevation_m,y=my_SE$Storage_mcm,xout=allelev)
my_SE <- data.frame(Elevation_m=interpolated$x,
                        Storage_mcm=interpolated$y)

# Apply lookup table to G-REALM data

interp_SE_ENVI <- approx(my_SE$Elevation_m, my_SE$Storage_mcm, method = "linear", xout = ENVItbl$Elevation_m)
ENVItbl$MyEstStorage_mcm <- interp_SE_ENVI$y

interp_SE_S3A <- approx(my_SE$Elevation_m, my_SE$Storage_mcm, method = "linear", xout = S3Atbl$Elevation_m)
S3Atbl$MyEstStorage_mcm <- interp_SE_S3A$y

##### STATS #####

# NSE for RS elev vs ROU elev
NSE(sim=ENVItbl$Elevation_m,obs=ENVItbl$ROU_Elevation_m)
NSE(sim=S3Atbl$Elevation_m,obs=S3Atbl$ROU_Elevation_m, na.rm=TRUE)

# pbias for RS elev vs ROU elev
pbias(sim=ENVItbl$Elevation_m,obs=ENVItbl$ROU_Elevation_m)
pbias(sim=S3Atbl$Elevation_m,obs=S3Atbl$ROU_Elevation_m)

# NSE for Yigzaw vs ROU Storage
NSE(sim=ENVItbl$Y_EstStorage_mcm,obs=ENVItbl$ROU_Storage_mcm)
NSE(sim=S3Atbl$Y_EstStorage_mcm,obs=S3Atbl$ROU_Storage_mcm)

# pbias for Yigzaw vs ROU Storage
pbias(sim=ENVItbl$Y_EstStorage_mcm,obs=ENVItbl$ROU_Storage_mcm)
pbias(sim=S3Atbl$Y_EstStorage_mcm,obs=S3Atbl$ROU_Storage_mcm)

# NSE for My SD vs ROU Storage
NSE(sim=ENVItbl$MyEstStorage_mcm,obs=ENVItbl$ROU_Storage_mcm)
NSE(sim=S3Atbl$MyEstStorage_mcm,obs=S3Atbl$ROU_Storage_mcm)

# pbias for My SD vs ROU Storage
pbias(sim=ENVItbl$MyEstStorage_mcm,obs=ENVItbl$ROU_Storage_mcm)
pbias(sim=S3Atbl$MyEstStorage_mcm,obs=S3Atbl$ROU_Storage_mcm)

##### PLOTS #####

# Choose consistent colors
colors=wes_palette("Royal2")

# Plot G-REALM elevation vs ResOpsUS Elevation
plot(x=ymd(ResOpsUSdata$Date),
     y=ResOpsUSdata$ROU_Elevation_m,
     type='l',lwd=1.5,
     col=colors[2],
     xlim=c(ymd(ENVIstartdate),ymd(S3Aenddate)),
     xlab='Date',
     ylab='Elevation (m)')
lines(x=ymd(ENVItbl$Date),
      y=ENVItbl$Elevation_m,
      col=colors[1],lwd=1.5)
lines(x=ymd(S3Atbl$Date),
      y=S3Atbl$Elevation_m,
      col=colors[1],lwd=1.5)
legend('topright',legend=c('G-REALM','ResOpsUS'),col=c(colors[1],colors[2]),lty=c(1,1),lwd=c(1.5,1.5))

# Plot lookup table curves
plot(x=my_SE$Storage_mcm,
     y=my_SE$Elevation_m-baseelev,
     type='l',lwd=1.5,
     col=colors[5],
     ylim=c(0,220),xlim=c(0,56000),
     xlab='Storage (mcm)',ylab='Depth (m)')
lines(x=Y_SDtbl$Storage_mcm,
      y=Y_SDtbl$Depth_m,
      col=colors[3],lwd=1.5)
lines(x=GRDL_SDtbl$Storage_mcm,
      y=GRDL_SDtbl$Depth_m,
      col=colors[4],lwd=1.5)
abline(h=max(ENVItbl$Depth_m),lty=2,col='grey',lwd=1.5)
abline(h=min(S3Atbl$Depth_m),lty=2,col='grey',lwd=1.5)
legend('bottomright',legend=c('Geometric','GRDL','Mine','G-REALM range'),
       col=c(colors[3],colors[4],colors[5],'grey'),lty=c(1,1,1,2),lwd=c(1.5,1.5,1.5,1.5))

# Plot just my lookup table
plot(x=my_SE$Storage_mcm,
     y=my_SE$Elevation_m,
     type='l',lwd=1.5,
     col=colors[5],
     xlab='Storage (mcm)',ylab='Elevation (masl)')

# Plot ResOpsUS Storage; ENVI est storage and S3A est storage from Yigzaw and my method
plot(x=ymd(ResOpsUSdata$Date),
     y=ResOpsUSdata$ROU_Storage_mcm,
     type='l',lwd=1.5,
     col=colors[2],
     ylim=c(12000,30000),
     xlim=c(ymd(ENVIstartdate),ymd(S3Aenddate)),
     xlab='Date',
     ylab='Storage (mcm)')
lines(x=ymd(ENVItbl$Date),
      y=ENVItbl$Y_EstStorage_mcm,
      col=colors[3],lwd=1.5)
lines(x=ymd(S3Atbl$Date),
      y=S3Atbl$Y_EstStorage_mcm,
      col=colors[3],lwd=1.5)
lines(x=ymd(ENVItbl$Date),
      y=ENVItbl$MyEstStorage_mcm,
      col=colors[5],lwd=1.5)
lines(x=ymd(S3Atbl$Date),
      y=S3Atbl$MyEstStorage_mcm,
      col=colors[5],lwd=1.5)
legend('topright',legend=c('Geometric','Mine','ResOpsUS'),
       col=c(colors[3],colors[5],colors[2]),lty=c(1,1,1),lwd=c(1.5,1.5,1.5))

# Zoom in on just my storage with max/min lines
plot(x=ymd(ResOpsUSdata$Date),y=ResOpsUSdata$ROU_Storage_mcm,
     type='l',col=colors[2],lwd=1.5,
     ylim=c(min(my_SE$Storage_mcm),25000), # Work on this line
     xlim=c(ymd(ENVIstartdate),ymd(S3Aenddate)),
     xlab='Date',ylab='Storage (mcm)')
lines(x=ymd(ENVItbl$Date),y=ENVItbl$MyEstStorage_mcm,
      col=colors[5],lwd=1.5)
lines(x=ymd(S3Atbl$Date),y=S3Atbl$MyEstStorage_mcm,
      col=colors[5],lwd=1.5)
abline(h=max(my_SE$Storage_mcm),lty=2,col='grey',lwd=1.5)
abline(h=min(my_SE$Storage_mcm),lty=2,col='grey',lwd=1.5)
abline(v=ENVIstartdate,lty=2,col='lightblue',lwd=1.5)
abline(v=S3Aenddate,lty=2,col='lightblue',lwd=1.5)
legend('topright',legend=c('Mine','ResOpsUS','Max/Min Calibration Range','G-REALM start/end dates'),
       col=c(colors[5],colors[2],'grey','lightblue'),lty=c(1,1,2,2),lwd=c(1.5,1.5,1.5,1.5))

