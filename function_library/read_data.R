########################
#
#Functions to read data from MOV_UP drives
#
#####################

library(pacman)

p_load("tools","lubridate", "zoo")

#ensure columns are numeric if necessary
makenum <- function(x) {
  as.numeric(as.character(x))
}


#Name runs according to established conventions (car 1, car 2)
getrunname <- function (filenameval) {
  
  if(!is.null(grep("AE51", filenameval))) 
  {
    
    
    keepfilenames = filenameval
    
    finalnames = 
      
      lapply(keepfilenames, FUN = function(x){
        namevalues = NULL
        othervals = NULL
        othervalsfinal = NULL
        keepsingle = x
        othervals  = keepsingle[grep("AE51", substr(keepsingle, 1,4))]
        
        
        if(length(othervals)==0){
          x = sub("_car1", "car1", x)
          x = sub("_car2", "car2", x)
          x = sub("car1", "car1_", x)
          x = sub("car2", "car2_", x)
          namevalues = unlist(tstrsplit(x, "_", fixed=T, keep=1))
          #namevalues = unlist(tstrsplit(filenameval, " ", fixed=T, keep=1))
          namevalues = sub("car","_Car", namevalues)
          namevalues = grep("AE51", namevalues, invert=T, value = T)
        }
        
        
        if(length(othervals)>0){
          
          othervals = gsub("__","_", othervals )
          
          othervalsdate = unlist(tstrsplit(othervals, "_", fixed=T, keep=2))
          othervalsdate = unlist(tstrsplit(othervalsdate, "-", fixed=T, keep=1))
          othervalsdate = as.POSIXct(othervalsdate, format = "%Y%m%d", tz="America/Los_Angeles")
          othervalsdate = format(othervalsdate, format = "%Y%b%d")
          
          carvals = rep(NA, length.out = length(othervalsdate))
          carvals[grep("car1", othervals,ignore.case = T)] = "_Car1"
          carvals[grep("car2", othervals,ignore.case = T)] = "_Car2"
          
          
          othervalsfinal = paste0(othervalsdate, carvals)
          othervalsfinal = sub("NA","",othervalsfinal)
          
        }
        
        as.data.table(paste(othervalsfinal, namevalues))
        
        
      })
    
    return(gsub(" ", "", unlist(finalnames)))
    
    
  }
  
  filenameval = sub("_car1", "car1", filenameval)
  filenameval = sub("_car2", "car2", filenameval)
  filenameval = sub("car1", "car1_", filenameval)
  filenameval = sub("car2", "car2_", filenameval)
  namevalues = unlist(tstrsplit(filenameval, "_", fixed=T, keep=1))
  #namevalues = unlist(tstrsplit(filenameval, " ", fixed=T, keep=1))
  namevalues = sub("car","_Car", namevalues)
  namevalues
}

#Obtain correct fixed site name
getlocation <- function (filenameval) {
  
  namevalues = NA
  if(length(grep("Maywood", filenameval, value=T, ignore.case=T))>0)
    namevalues = "Maywood"
  if(length(grep("SandPt", filenameval, value=T, ignore.case=T))>0)
    namevalues = "SandPt"
  if(length(grep("SeaTacCC", filenameval, value=T, ignore.case=T))>0)
    namevalues = "SeaTacCC"
  if(length(grep("10th", filenameval, value=T, ignore.case=T))>0)
    namevalues = "10thWeller"
  if(length(grep("Weller", filenameval, value=T, ignore.case=T))>0)
    namevalues = "10thWeller"
  if(length(grep("W.dat", filenameval, value=T, ignore.case=T))>0)
    namevalues = "10thWeller"
  # if(length(grep("10th\+W", filenameval, value=T, ignore.case=T))>0)
  #   namevalues = "10thWeller"
  # if(length(grep("10thW", filenameval, value=T, ignore.case=T))>0)
  #   namevalues = "10thWeller"
  if(is.na(namevalues))
    namevalues = filenameval

  namevalues
}

#Value of average indicates the end of the sample interval
averageTime <- function (datetime, timeaverage) {
  #this establishes a common time scale for averaging
  as.POSIXct(ceiling(as.numeric(datetime)/(timeaverage*60))*(timeaverage*60),origin='1970-01-01')
}

fread2 <- function(filename) {
  tmp <- scan(file = filename, what = "character", sep="\t",quiet = TRUE,nlines=1)
  # remove empty lines caused by \r
  tmp <- tmp[tmp != ""]
  # paste lines back together together with \n character
  tmp <- paste(tmp, collapse = "\n")
  test = fread(tmp,fill=TRUE)
  test
}


read.gps <- function(datafile, runname, location, timeaverage, splineval = T) {
  
  dt <- fread(datafile)
  dt[, datetime.gps := 
       parse_date_time(paste(Date,`Time (local)`),
                       orders = c("%Y-%m-%d %H:%M:%S",
                                  "%m/%d/%Y %H:%M:%S"),
                       tz="America/Los_Angeles")]
  dt = dt[!is.na(datetime.gps),]
  dt[, runname := runname]
  dt[,Date:= NULL]
  dt[,`Time (local)` := NULL]
  #get Correct time average
  
  #get averating interval
  dt[,timeint := averageTime(datetime.gps, timeaverage)]
  
  #average data columns by timeint. Ensure that distinct runnames are kept seperated
  dt = 
    dt[,lapply(.SD, mean), by=c("timeint","runname"), 
       .SDcols=c("Latitude","Longitude","Altitude (m)","Speed (km/hr)")]
  
  if(splineval) {
    
    transcol = c("Latitude","Longitude","Altitude (m)","Speed (km/hr)")
    dt[, (transcol)  := 
         lapply(.SD, function(x)
           na.approx(x, na.rm=F, rule=1)), .SDcols = transcol ]
  }
  
  dt$location = location
  
  return(dt)
}

read.langan <- function(datafile, runname, location, timeaverage, splineval = T) {
  
  dt <- fread(datafile, skip=1,na.string="")
  dt[, `#` := NULL]
  #ask Time about the TimeZONE
  dt[, datetime.langan := as.POSIXct(`Date Time, GMT-07:00`, format="%m/%d/%y %I:%M:%S %p", tz="America/Los_Angeles")]
  dt[, runname := as.character(runname)]
  dt[,`Date Time, GMT-07:00`:= NULL]
  #get Correct time average
  setnames(dt, grep("CO", colnames(dt), value=T),"CO")
  setnames(dt, grep("Temp", colnames(dt), value=T), "TempC" )
  setnames(dt, grep("Stopped", colnames(dt), value=T), "Stopped")
  
  
  
  dt[,timeint := averageTime(datetime.langan, timeaverage)]
  dt = 
    dt[,lapply(.SD, mean), by=c("timeint","runname"), .SDcols=c("TempC","CO")]
  
  if(splineval) {
    
    transcol = c("CO","TempC")
    dt[, (transcol)  := 
         lapply(.SD, function(x)
           na.approx(x, na.rm=F, rule=1)), .SDcols = transcol ]
  }
  dt$location = location
  return(dt)
}

read.ptrak <- function(datafile, runname, location, timeaverage, screen=F,  splineval = T) {
  
  dt <- fread(datafile, skip="MM/dd/yyyy", fill=T)
  serial = fread(datafile, skip = "Serial Number:", nrows = 1, header = F)[,2]
  pncname = ifelse(screen, "pnc_screen","pnc_noscreen")
  serialname = ifelse(screen, "serial.screen", "serial.noscreen")
  colnames(dt) = c("Date", "Time", pncname)
  dt[, datetime.ptrak := as.POSIXct(paste0(Date,Time), format="%m/%d/%Y %T", tz="America/Los_Angeles")]
  dt = dt[!is.na(datetime.ptrak)]
  dt[, runname := runname]
  dt[, Date:= NULL]
  dt[, Time := NULL]
  #get Correct time average
  
  
  dt[,timeint := averageTime(datetime.ptrak, timeaverage)]
  
  dt = 
    dt[,lapply(.SD, FUN = function(x)
      mean(as.numeric(as.character(x)))), 
      by=c("timeint","runname"), 
      .SDcols=c(pncname)]
  
  if(splineval) {
    
    transcol = c(pncname)
    dt[, (transcol)  := 
         lapply(.SD, function(x)
           na.approx(x, na.rm=F, rule=1)), .SDcols = transcol ]
  }
  
  dt$location = location
  dt$serial = serial
  colnames(dt) = c("timeint", "runname", pncname, "location", serialname)
  return(dt)
}

read.co2 <- function(datafile, runname, location, timeaverage, splineval = T, serial = NA) {
  
  dt <- read.table(datafile, skip=1,na.string="?", sep="\t", fill=T, header=T)
  dt = data.table(dt)
  
  try(expr = dt[, "H?.O_.?C." := NULL], silent = T)
  try(expr = dt[, X := NULL])
  try(expr = dt[, `Input_Voltage_.V.` := NULL], silent = T)
  
  if(ncol(dt) ==9)
    colnames(dt)[1:9] = c("Date","Time","CO2 (umol/mol)","H2O (mmol/mol)", 
                          "Cell Temp (C)", "Cell Pressure (kPa)", 
                          "CO2 Absorp", "H2O Absorp", "Flow (L/min)")
  
  if(ncol(dt) ==8)
    colnames(dt)[1:8] = c("Date","Time","CO2 (umol/mol)","H2O (mmol/mol)", 
                          "Cell Temp (C)", "Cell Pressure (kPa)", 
                          "CO2 Absorp",  "Flow (L/min)")
  
  #ask Tim about the TimeZONE
  dt[, datetime.co2 := parse_date_time(paste(Date,Time), order= "Ymd HMS",
                                       tz="America/Los_Angeles")]
  
  dt = dt[!is.na(datetime.co2),]
  dt[, runname := as.character(runname)]
  dt[,`Date`:= NULL]
  dt[,`Time`:= NULL]
  
  dt[,timeint := averageTime(datetime.co2, timeaverage)]
  dt = 
    dt[,lapply(.SD, FUN = function(x){
      mean(makenum(x))}), by=c("timeint","runname"), .SDcols=c("CO2 (umol/mol)",
                                                               "CO2 Absorp",
                                                               "H2O (mmol/mol)")]
  
  if(splineval) {
    
    transcol = c("CO2 (umol/mol)",
                 "CO2 Absorp",
                 "H2O (mmol/mol)")
    dt[, (transcol)  := 
         lapply(.SD, function(x)
           na.approx(x, na.rm=F, rule=1)), .SDcols = transcol ]
  }
  
  dt$location = location
  dt$serial.co2 = serial
  return(dt)
  
}

read.cpc <- function(datafile, runname, location, timeaverage,  splineval = T) {
  
  #dt <- data.table(read.csv(datafile,sep=",", header=F))
  dt = data.table(fread(datafile, skip = 17, header = T))
  
  
  startdate =  fread(datafile, skip = "Start Date", nrows = 1, header = F)

  
  startdate = parse_date_time(startdate, c("%m/%d/%Y", "%m/%d/%y", "%Y-%m-%d"), 
                              tz="America/Los_Angeles")
  
  dt = lapply(1:floor(ncol(dt)/2), FUN = function(set){
    indexval = (1+2*(set-1)):(1+2*(set-1)+1)
    temp = dt[,..indexval]
    temp[,Time := as.POSIXct(paste(startdate[indexval[2]], Time), 
                             format="%Y-%m-%d %H:%M:%S")]
    temp
  })
  
  dt = rbindlist(dt)
  
  try(
    {serial = fread(datafile, skip = "Instrument ID", nrows = 1, header = F)[,2]})
  
  colnames(dt)[1] = "Time"
  colnames(dt)[2] = "Concentration (#/cc)"
  
  
  
  dt = dt[!is.na(Time),]
  
  dt = dt[,1:2,with=FALSE]
  
  dt[, runname := runname]
  #get Correct time average
  
  setnames(dt, "Time", "datetime.cpc")
  
  
  dt[,timeint := averageTime(datetime.cpc, timeaverage)]
  
  dt = 
    dt[,lapply(.SD, FUN = function(x) {
      mean(as.numeric(as.character(x)), na.rm=T)}), 
      by=c("timeint","runname"), 
      .SDcols=c("Concentration (#/cc)")]
  
  if(splineval) {
    
    transcol = c("Concentration (#/cc)")
    dt[, (transcol)  := 
         lapply(.SD, function(x)
           na.approx(x, na.rm=F,rule=1 )), .SDcols = transcol ]
  }
  
  try({
    dt$location = location
    dt$serial.cpc = serial
  })
  return(dt)
}


read.ae51 <- function(datafile, runname, location, timeaverage,  splineval = T) {
  
  if(file_ext(datafile) =="dat")
  {
    dt <- data.table(read.table(datafile,sep=";", skip=15, header=T))
    serial = fread(datafile,skip = "Device ID", nrows = 1, header = F)[,4]
  }
  if(file_ext(datafile) =="csv")
  {
    dt <- data.table(read.csv(datafile, sep=",", skip=15))
    serial = fread(datafile,skip = "Device ID", nrows = 1, header = F)[,4]
    
  }
  if(file_ext(datafile) =="txt")
  {
    dt <- data.table(read.delim(datafile, sep="\t", skip=15))
    serial = fread(datafile,skip = "Device ID", nrows = 1, header = F)[,1]
    serial = substr(serial, 13, nchar(serial))
  }
  
  #serial = fread(datafile,skip = "Device ID", nrows = 1, header = F)[,4]
  dt[, datetime.ae51 := parse_date_time(paste(Date,`Time`),c("Ymd HMS","mdY HMS","mdY HM"), tz="America/Los_Angeles")]
  dt[, runname := runname]
  dt[,Date:= NULL]
  dt[,`Time` := NULL]
  #get Correct time average
  
  
  dt[,timeint := averageTime(datetime.ae51, timeaverage)]
  
  dt = 
    dt[,lapply(.SD, mean), by=c("timeint","runname"), 
       .SDcols=c("Ref","Sen","ATN","Flow","PCB.temp","Status","Battery",
                 "BC")]
  
  #remove leading NA
  dt[, .SD[!!cumsum(!is.na(BC))], by = "runname"]
  
  if(splineval) {
    
    transcol = c("Ref","Sen","ATN","Flow","PCB.temp","Status","Battery",
                 "BC")
    dt[, (transcol)  := 
         lapply(.SD, function(x)
           
           na.approx(x, na.rm=F,rule = 1)), .SDcols = transcol ]
  }
  
  dt$location = location
  dt$serial.ae51 = serial
  
  #smooth data with ONA model, se Haggler paper
  dt[,ATN_ONA := round_any((ATN)-max(ATN)-0.005, 0.06, ceiling) + max(ATN)]
  dt[, BC_ONA := mean(BC, na.rm = T), by=c("runname","ATN_ONA")]
  dt[BC_ONA < 0, BC_ONA := NA]
  return(dt)
}


read.labview <- function(datafile, runname, location, timeaverage) {
  
  dt <- data.table(read.table(datafile,sep="\t", header=T))
  
  dt = dt[,which(unlist(lapply(dt, function(x)!all(is.na(x))))),with=F]
  
  cols.labview=c("Computer Time Stamp",
                 "Marker",
                 "GPS  Mode, 1=nofix, 2=2D, 3=3D",
                 "GPS  No. of Active Satellites",
                 "GPS  Hor. Precision (HDOP)",
                 "GPS  Time Stamp",
                 "GPS  Latitude (deg)",
                 "GPS  Longitude (deg)",
                 "GPS  Speed (km/h)",
                 "GPS  Direction (deg)",
                 "Time NO",
                 "Time NOx",
                 "AE52 Time",
                 "AE52Flow","AE52Status",
                 "Precon HS-2000 Temp (°C)",
                 "Precon HS-2000 RH (%)",
                 "SenseAir CO2 conc. (ppm)")
  
  keep.labview=c("Computer Time Stamp",
                 "Marker",
                 "GPS  Mode, 1=nofix, 2=2D, 3=3D",
                 "GPS  No. of Active Satellites",
                 "GPS  Hor. Precision (HDOP)",
                 "GPS  Time Stamp",
                 "GPS  Latitude (deg)",
                 "GPS  Longitude (deg)",
                 "GPS  Speed (km/h)",
                 "GPS  Direction (deg)",
                 "Precon HS-2000 Temp (°C)",
                 "Precon HS-2000 RH (%)",
                 "SenseAir CO2 conc. (ppm)")
  
  setnames(dt, colnames(dt),cols.labview)
  
  num.labview=c( "GPS  No. of Active Satellites",
                 "GPS  Time Stamp",
                 "GPS  Latitude (deg)",
                 "GPS  Longitude (deg)",
                 "GPS  Speed (km/h)",
                 "GPS  Direction (deg)",
                 "Precon HS-2000 Temp (°C)",
                 "Precon HS-2000 RH (%)",
                 "SenseAir CO2 conc. (ppm)")
  
  
  dt[, datetime.labview := as.POSIXct(`Computer Time Stamp`,
                                      format="%d/%m/%Y %I:%M:%S %p",
                                      tz="America/Los_Angeles")]
  
  dt[, `GPS  Time Stamp` := as.POSIXct(`GPS  Time Stamp`,
                                       format="%d/%m/%Y %I:%M:%S %p",
                                       tz="America/Los_Angeles")]
  dt[, runname := runname]
  dt[,`Computer Time Stamp`:= NULL]
  
  #get Correct time average
  
  dt[,timeint := averageTime(datetime.labview, timeaverage)]
  
  dt = 
    dt[,lapply(.SD, mean), by=c("timeint","runname"), 
       .SDcols=num.labview]
  
  if(splineval) {
    
    transcol = num.labview
    dt[, (transcol)  := 
         lapply(.SD, function(x)
           na.approx(x, na.rm=F, rule=1)), .SDcols = transcol ]
  }
  
  dt$location = location
  return(dt)
}

read.nano.scan <- function(datafile, runname, location, timeaverage, splineval=T) {
  
  dt <- fread(datafile,skip="File", header=T)
  dt[, datetime.nano.scan := parse_date_time(`Date Time`,c("Ymd HMS", "mdY HMS", "mdy HMS", "ymd HMS","%m/%d/%y %H:%M"), tz="America/Los_Angeles")]
  #correct so time is at the end of the interval
  dt[, datetime.nano.scan := datetime.nano.scan + 60]
  dt[, runname := runname]
  dt[,`Date Time`:= NULL]
  #get Correct time average
  
  
  dt[,timeint := averageTime(datetime.nano.scan, timeaverage)]
  
  dt = 
    dt[,lapply(.SD, mean), by=c("timeint","runname","Status")]
  
  if("154" %in% colnames(dt))
    setnames(dt, "154", "154.0")
  
  if(splineval) {
    
    transcol = c("11.5","15.4","20.5","27.4","36.5","48.7","64.9",
                 "86.6","115.5","154.0","205.4","273.8","365.2" ,
                 "Total Conc", "Median (nm)","Mean (nm)",
                 "Geo Mean (nm)", "Mode (nm)", "GSD",
                 "Particle Density (g/cc)")
    dt[, (transcol)  :=
         lapply(.SD, function(x)
           na.approx(x, na.rm=F, rule=1)),
       .SDcols = transcol ]
  }
  
  dt$location = location
  dt[,timeint := as.POSIXct(timeint)]
  return(dt)
}

read.nano.single <- function(datafile, runname, location, timeaverage, splineval=T) {
  
  #get datetime
  startval= read.csv(datafile, skip=7, nrow=1,colClasses="character")[2]
  timeval = start= read.csv(datafile, skip=8, nrow=1,colClasses="character")[2]
  start.time = parse_date_time(paste(startval, timeval), c("Ymd HMS", "mdY HMS"),tz="America/Los_Angeles")
  dt <- read.csv(datafile,skip=14, header=T)
  dt = data.table(dt)
  setnames(dt,names(dt), c("Time elapsed","single channel number/cc","status"))
  
  dt[,iteration := -999]
  dt$iteration[1] = 1
  dt[iteration==-999]$iteration = NA
  
  timecorr = data.table(Date = c(startval,t(dt[`Time elapsed`=="Date"][,2])),
                        Time = c(timeval,t(dt[`Time elapsed`=="Time"][,2])))
  
  timecorr[,datetime := parse_date_time(paste(as.character(timecorr$Date), 
                                              as.character(timecorr$Time)), c("Ymd HMS", "mdY HMS"),tz="America/Los_Angeles")]
  
  if(ncol(timecorr)>1){
    dt[`Time elapsed`=="Date", iteration := as.double((2:(nrow(timecorr))))]
    timecorr$iteration = 1:nrow(timecorr)
  }
  
  dt[,iteration := na.locf(iteration)]
  
  setkey(dt, iteration)
  setkey(timecorr, iteration)
  
  dt = dt[timecorr]
  
  dt[, runname := runname]
  
  dt[,`Time elapsed` := as.numeric(as.character(`Time elapsed`))]
  
  dt = dt[!is.na(`Time elapsed`),]
  
  dt[,datetime := datetime + `Time elapsed`]
  
  dt[,timeint := averageTime(datetime, timeaverage)]
  
  dt[,iteration := NULL]
  dt[, `Time elapsed` := NULL]
  dt[,Date := NULL]
  dt[,Time :=NULL]
  
  dt = 
    dt[,lapply(.SD, FUN = 
                 function(x) mean(as.numeric(as.character(x)))), 
       by=c("timeint", "status", "runname"),
       .SDcols = c("single channel number/cc")]
  
  if(splineval) {
    
    transcol = c("single channel number/cc")
    dt[, (transcol)  := 
         lapply(.SD, function(x)
           na.approx(x, na.rm=F, rule=1)), .SDcols = transcol ]
  }
  
  
  dt$location = location
  return(dt)
}

#round data.table
round_df <- function(df.table, digits) {
  
  nums <- names(df.table)[unlist(df.table[,lapply(.SD, is.numeric)])]
  
  df.table[, (nums = lapply(.SD, round, digits)), .SDcols=nums]
  
  df.table
}
