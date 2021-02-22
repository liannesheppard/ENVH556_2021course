# load dependencies
pacman::p_load(riem)

get_ASOS <- function(date_start, date_end, station) {

measures <- riem_measures(station = "SEA", 
                          date_start = date_start, 
                          date_end = date_end) 

measures = data.table(measures)

measurest = copy(measures)

transcol = c("tmpf","dwpf", "relh","drct","sknt","alti","mslp","vsby")
measurest[, (transcol)  := 
           lapply(.SD, function(x)
             na.spline(x, na.rm=F, maxgap=15 )), .SDcols = transcol ]

setnames(measurest, "valid","timeint")

measurest[,timeint := as.POSIXct(timeint)]
measurest
}

# station: three or four character site identifier
# valid: timestamp of the observation (UTC)
# tmpf: Air Temperature in Fahrenheit, typically @ 2 meters
# dwpf: Dew Point Temperature in Fahrenheit, typically @ 2 meters
# relh: Relative Humidity in \%
# drct: Wind Direction in degrees from north
# sknt: Wind Speed in knots
# p01i: One hour precipitation for the period from the observation time to the time of the previous hourly precipitation reset. This varies slightly by site. Values are in inches. This value may or may not contain frozen precipitation melted by some device on the sensor or estimated by some other means. Unfortunately, we do not know of an authoritative database denoting which station has which sensor.
# alti: Pressure altimeter in inches
# mslp: Sea Level Pressure in millibar
# vsby: Visibility in miles
# gust: Wind Gust in knots
# skyc1: Sky Level 1 Coverage
# skyc2: Sky Level 2 Coverage
# skyc3: Sky Level 3 Coverage
# skyc4: Sky Level 4 Coverage
# skyl1: Sky Level 1 Altitude in feet
# skyl2: Sky Level 2 Altitude in feet
# skyl3: Sky Level 3 Altitude in feet
# skyl4: Sky Level 4 Altitude in feet
# presentwx: Present Weather Codes (space seperated), see e.g. this manual for further explanations.
# metar: unprocessed reported observation in METAR format
