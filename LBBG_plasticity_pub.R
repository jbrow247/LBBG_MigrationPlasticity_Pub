##Load all packages ----

library(fossil)
library(adehabitatHR)
library(maps)
library(grid)
library(gridExtra)
library(lubridate)
library(tidyverse)
library(rgdal)
library(rworldmap)
library(rgeos)
library(RODBC)   ##connect to data base
library(sp)
select <- dplyr::select
##other packages used in this script: cowplot, data.table

##Load all user defined functions.----
norm_vec <- function(x) sqrt(sum(x^2))

new_point <- function(p0, p1, di) { # Finds point in distance di from point p0 in direction of point p1
  v = p1 - p0
  u = v / norm_vec(v)
  return (p0 + u * di)
}

## A function to calculate the distances between consecutive points ##
## Output in meters ##
pt2pt.distance <- function(latitude, longitude, lag = 1){
  require(fossil)
  distance <- NA
  for(i in 2:length(latitude)){
    distance[i] <- deg.dist(long1= longitude[i-lag],lat1 = latitude[i-lag], long2 = longitude[i], lat2 = latitude[i] )*1000 }
  return(distance)
}

## A function to calculate the time increment between consecutive points ##
## Default output in seconds, other options are "auto", "mins", "hours","days", or "weeks" ##
pt2pt.duration <- function(datetime, output.units='secs'){
  duration <- NA
  for(i in 2:length(datetime)){
    duration[i] <- difftime(datetime[i], datetime[i-1], units=output.units) }
  return(duration)
}

## A function to calculate the speed of movement between consecutive points ##
pt2pt.speed <- function(distance, duration){
  return(distance/duration)
}

###A function that finds the closest points to a regular (e.g. hourly) time series
trackSubSamp = function(df, int=1,unit='hours')
{
  id.i = unique(df$id)
  n.id = length(id.i)
  df.sub = list(n.id)
  timestep = paste(int,unit)
  
  # breakdown to datasets per bird
  for (i in 1:n.id)
  {
    df.i = df[df$id==id.i[i],]
    dat.seq <- seq(from = min(df.i$time, na.rm=T), to = max(df.i$time, na.rm=T), by = timestep)
    id.sub = sapply(dat.seq, function(x) which.min(abs(difftime(df.i$time, x, units='mins')))) #find gps points minimizing distance to each ts in dat.seq
    
    df.sub[[i]] = unique(df.i[id.sub,])     
    # the function unique makes sure that the rows in df.i[idx,] are unique - so no duplicate points
  }

  df.sub <- data.table::rbindlist(df.sub)
}
# 
# ##Same as trackSubSamp, but usable in a loop
# trackSubSamp.birdyear.loop = function(id.df,int=1,unit='hours'){
#   timestep = paste(int,unit)
#   dat.seq <- seq(from = min(id.df$ts, na.rm=T), to = max(id.df$ts, na.rm=T), by = timestep)
#   id.sub = sapply(dat.seq, function(x) which.min(abs(difftime(df.i$ts, x, units='secs')))) #find gps points minimizing distance to each ts in dat.seq
#   
#   df.sub[[i]] = unique(df.i[id.sub,])     
#   # the function unique makes sure that the rows in Dtemp[idx,] are unique - so no duplicate points
#   return(df.sub)
# }

#laea.proj <- "+proj=laea +lat_0=34.9 +lon_0=-4.16 +x_0=4321000 +y_0=3210000 +towgs84=0,0,0,0,0,0,0 +units=m"
laea.proj <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +towgs84=0,0,0,0,0,0,0 +units=m"

##defult project projection is lambert equal-area projection, which is centered on the spatialy center of the data (midpoint between most extreme values). This projection uses meters and is recommended for statistic analysis. 
mkSpatial = function(df, CRS.in = "+init=epsg:4326", 
                     CRS.out = laea.proj) {
  library(sp)
  library(rgdal)
  
  df$all = rep('all',nrow(df))
  df = SpatialPointsDataFrame(coords = cbind(df$lon,df$lat), data = df[,c('id','all','time', 'lon', 'lat', 'birdyear')], 
                              proj4string = CRS(CRS.in))        # the data is in WSG84 lat-lon
  
  df = spTransform(df, CRS(CRS.out)) 
  
}

mkgrid <- function(lat=NULL, lon = NULL, resolution=10000, buffer = 100000, 
                   projection = laea.proj  ){
  #lat is a vector of latitudes
  #lon is a vector of longitudes
  #resolution is the cell size of the raster. 
  #1 = ~ 133X165 km
  #.1 = ~ 13.9 x 17.1 km
  #.01 ~ 1.39 x 1.72 km     .... in WSG84
  #buffer adds space beyond the points in the grid. a buffer of 1 should equal ~ 111 km
  #projection: projection of the data.
  
  ##projection = laea, measured in m
  #resolution 1 000 = 1 km2...
  
  
  xmax <- max(lon)
  xmin <- min(lon)
  ymax <- max(lat)
  ymin <- min(lat)
  extent<- matrix(c(xmin - buffer, xmax + buffer, ymin - buffer, ymax + buffer), ncol=2)
  extent <- SpatialPoints(coords = extent, proj4string = CRS(projection))        # the data is in WSG84 lat-lon
  rast<- ascgen(extent, cellsize=resolution)
}


# ## Cleaning ----
gull <- readRDS("multitrack_gull_clean.RDS") ##running the above results in different daily.by??
gull <- gull %>% group_by(device_info_serial) %>%
  mutate(birdyear =  year(date_time-days(152)) - year(min(date_time)),  ## a bird year starts on year day 152 (usually jun1, except leap years), 
         id_birdyear = paste(device_info_serial, birdyear, sep = ".")) %>%
  filter(birdyear >=0 & date_time < ymd(paste0(year(max(date_time)),'-06-01'))) %>% #remove points in last year following jun 1 that don't contain migration/winter points
  ungroup()
gull <-  select(gull, id = device_info_serial, birdyear, id_birdyear, time = date_time, lat = latitude, lon = longitude, alt = altitude,
                start_date, start_lat = start_latitude, start_lon = start_longitude, mass, sex, dur, dist, speed.pts,
                col = key_name) %>% mutate(year = year(time))
original.n <- c(length(unique(gull$id_birdyear)), length(unique(gull$id)))
gull <- filter(gull, !(id == 757)) ## changed colonies!
gull <- filter(gull, !(id == 4023)) ## Doesn't migrate!
gull <- filter(gull, !(id_birdyear %in% c("325.1", "497.1","540.2"))) ## no arrival to colony (gps didn't finish downloading )

gull$col <- factor(gull$col)

gull.meta <- gull %>% select(id, start_date, start_lat, start_lon, mass, sex, col) %>% distinct()
gull <- gull %>% select(-mass, -sex, -col)

##Calculate distance to colony
gull <- gull %>% mutate(d2_col = deg.dist(lon, lat, start_lon, start_lat))

#Find first (=arrival time) and last (=departure time) visit to colony in a bird year
arrival.date <- gull %>%  filter(year != year(start_date) &   #no arrival in first year
                                   d2_col <= 10) %>%  ##select points within the 10k buffer of colony
  group_by(id, year) %>%
  summarise(col.arrival = min(time)) %>% ungroup() ##Find the first point within 10k of colony

departure.date <- gull %>% group_by(id) %>% mutate(end.year = max(year)) %>%
  ungroup() %>%
  filter(year != end.year &
           d2_col <= 10) %>%  ##select points within the 10k buffer of colony
  group_by(id, year) %>%
  summarise(col.depart = max(time)) %>% ungroup() ##Find the first point within 10k of colony

##join with gull to re-specify birdyear
arrival.date <- gull.meta %>% select(id, start_date) %>% full_join(arrival.date) %>%
  mutate(birdyear = year(col.arrival) - year(start_date) -1) %>% ## -1 because bird year of arrival corresponds to previous years departure
  select(id, col.arrival, birdyear)

departure.date <- gull.meta %>% select(id, start_date) %>% full_join(departure.date) %>%
  mutate(birdyear = year(col.depart) - year(start_date)) %>% select(id, col.depart, birdyear)

mig.date <- full_join(arrival.date, departure.date) %>% group_by(id) %>%
  arrange(birdyear) %>%
  mutate(year.start = (col.depart - lag(col.arrival))/2 + lag(col.arrival)) %>%
  left_join(gull.meta) %>%
  mutate(year.start = as.POSIXct(ifelse(birdyear == 0, start_date, year.start),
                                 origin = ymd("1970-01-01"))) %>%
  select(-col, -start_date, - start_lat, -start_lon) %>%
  mutate(year.end = lead(year.start)) %>%
  mutate(year.end = as.POSIXct(ifelse(birdyear == max(birdyear), col.arrival, year.end),
                               origin = ymd("1970-01-01"))) %>% ungroup() %>%
  filter(!is.na(year.start) & !is.na(year.end))  ##remove birds with missing departures/arrivals (e.g. year long gaps, no return data, etc. )

gull <- gull %>% select(-start_lon, -start_lat, -birdyear, -id_birdyear, -start_date)
gull <- gull %>% right_join(mig.date) %>%
   filter(time <= year.end & time >= year.start) %>%
  mutate(id_birdyear = paste(id, birdyear, sep = ".")) %>%
  select(-year.start, -year.end)

nopair <- names(table((gull %>% select(id, birdyear) %>% distinct())$id))[which(table((gull %>% select(id, birdyear) %>% distinct())$id)<2)]
gull <- filter(gull, ! id %in% nopair)

nonbreed <- gull %>% 
  filter((time >= floor_date(col.depart, unit = "days")) & 
           (time <= ceiling_date(col.arrival, unit = "days"))) %>% ##use rounded start and end dates to ensure arrival points are included in data. This is important for identifying gaps at the end of the year
  select(-col.depart, -col.arrival)

gull <- select(gull, -col.arrival, -col.depart)

rm(gull, mig.date, nopair)
## gull.glm, birdyear summary info----

tmp <- full_join(departure.date, arrival.date)
gull.glm <- left_join(gull.meta, tmp) %>% rename(start = start_date) 
rm(gull.meta)

## KDEs, daily ----
#subsample data
nonbreed.d <- as.data.frame(trackSubSamp(nonbreed,12)) ## 2 points per day

## point2point
nonbreed.d <- nonbreed.d %>% group_by(id_birdyear) %>%
  mutate(dur = pt2pt.duration(time), dist = pt2pt.distance(lat, lon), speed.pts = pt2pt.speed(dist, dur)) %>%
  ungroup()

##daily id
##Find birdyears with 1 pt per day
daily.by <- nonbreed %>% group_by(id, birdyear) %>% 
  summarise(max.gap = max(dur, na.rm=T)/3600) %>% 
  filter(max.gap <= 24)  ##131 birdyears with no gaps!!

id.v <-  unique(daily.by$id)

#Conver data to spatial points df
snonbreed.d <- mkSpatial(nonbreed.d)

#create grid
snonbreed.pt <- as.data.frame(snonbreed.d@coords)
ud.grid <- mkgrid(lat=snonbreed.pt$coords.x2, 
                  lon = snonbreed.pt$coords.x1, 
                  resolution = 10000, 
                  buffer = 450000, 
                  projection = laea.proj)

# get map
map <- getMap("coarse")
laea.map <- spTransform(map, CRS(laea.proj))
laea.map <- fortify(laea.map)    
map <- fortify(map)

##KDEs for ids with daily points
##calculate utilization distribution for each individual/year
multi.yr <- lapply(1:length(id.v), function(i) {  
  id.d <- filter(nonbreed.d,  id == id.v[i]) ## get all points for one id
  sid <- mkSpatial(id.d)  ##turn into spatial df
  ud <-  kernelUD(sid['birdyear'], h = 100000, grid = ud.grid)}) ##calculate UD per birdyear (sid)

## get 50% UD
multi.yr.50 <- lapply(1:length(id.v), function(i){
  name <- paste("ver50", id.v[i], sep = ".")
  assign(name, getverticeshr(multi.yr[[i]], 50))
})  

#Split 50% core area into seperate polygons
multi.yr.50p <- lapply(1:length(multi.yr.50), function(i) {
  p <- disaggregate(multi.yr.50[[i]])
  p@data$birdyear <- p@data$id #clarify variable name
  p@data$id <- rownames(p@data) #create id variable of rownames so data can be merged with fortified polygon
  p
}) 

## Gap: determin threshold and remove ----

##Filter points outside of any core area 
#or breeding colony- 
#this is so when a bird moves outside of the core area polygon 
#(e.g. long central placed foraging trip), 
#and returns to the same polygon, it does not count as a new visit to this core area. 

nonbreed.core <- data.frame()
for (i in 1:length(daily.by$id)) {
  id.n <- daily.by[i,]$id
  by <- daily.by[i,]$birdyear
  pts <- nonbreed.d[nonbreed.d$id_birdyear %in% paste(id.n, by, sep = "."),]
  pol <- multi.yr.50p[[which(id.n == id.v)]]
  s.pts <- mkSpatial(pts)
  pol.y <- pol[pol@data$birdyear == by,]
  overlap <- s.pts[pol.y,] ##id and year matched
  overlap@data$poly <- over(overlap, pol.y)$id #which polygon it overlapped with
  r.le <- rle(overlap@data$poly) #number of points in each polygon per each temporal visit
  overlap@data$p.visit <- rep(seq(1:length(r.le$lengths)), r.le$lengths) #each temporally distinct visit to a polygon has a unique id (e.g. if moved from poly 1, to poly 2, then back to poly 1, the corresponding p.visits will be 1, 2, and 3 )
  nonbreed.core <- rbind(nonbreed.core, overlap@data) ##overlap points from one id
  }  

## Select first and last point in each new polygon
core.ee <- nonbreed.core %>% 
  group_by(id, birdyear, poly, p.visit) %>% 
  mutate(ee = ifelse(time == min(time), "en", ifelse (time == max(time), "ex", NA))) %>%  ## a few p.visits only have 1 point in their last p visit (considered entrance). This helps filter the right points to calculate gap length
  filter(time == min(time) | time == max(time)) %>%
  ungroup() 

##Select birdyears with 1 pt per day
poly.dur <- core.ee %>% select(id, time, ee, poly, p.visit, birdyear) %>% spread(ee, time) %>%  
  mutate(dur = difftime(ex,en, units = "days")) %>% group_by(id, birdyear, poly) %>% 
  summarise(total_dur = as.numeric(sum(dur, na.rm = T))) 
poly.dur  %>% arrange(total_dur) ##shortest duration = 21 days (one of zero because migrated through fall poly in spring)

#ids with gap > 21 days
id.w.gap <-  filter(nonbreed.d, dur >= 3600*24*21) %>% select(id_birdyear) %>% 
  distinct() %>% pull(id_birdyear) 


## remove ids with gap, and ids that were then left without a pair
nonbreed.nogap.d <- filter(nonbreed.d, !(id_birdyear %in% id.w.gap))
pair.miss <- nonbreed.nogap.d %>% select(id, birdyear) %>% distinct() %>% group_by(id) %>% summarise(n=n()) %>% filter(n==1) %>% pull(id) #an additional 20 rows lost
nonbreed.nogap.d <- filter(nonbreed.nogap.d, ! (id %in% pair.miss)) ## 84 birds remain, with  237 years 
id.v <- unique(nonbreed.nogap.d$id)

nonbreed.nogap <- filter(nonbreed, !(id_birdyear %in% id.w.gap))
pair.miss <- nonbreed.nogap %>% select(id, birdyear) %>% distinct() %>% group_by(id) %>% summarise(n=n()) %>% filter(n==1) %>% pull(id) #an additional 20 rows lost
nonbreed.nogap <- filter(nonbreed.nogap, ! (id %in% pair.miss)) ## 84 birds remain, with  237 years 

saveRDS(nonbreed, "nonbreed.RDS")
rm(nonbreed, nonbreed.d)
## KDEs - no b.y. with gap ----
## create new KDE list
nogap_multi.yr<- lapply(1:length(id.v), function(i) {
  id.h <- filter(nonbreed.nogap.d,  id == id.v[i])
  sid <- mkSpatial(id.h)
  ud <-  kernelUD(sid['birdyear'], h = 100000, grid = ud.grid)})
names(nogap_multi.yr) <- id.v

nogap_multi.yr.50 <- lapply(1:length(id.v), function(i){
  name <- paste("ver50", id.v[i], sep = ".")
  assign(name, getverticeshr(nogap_multi.yr[[i]], 50))
}) ##core areas

nogap_multi.yr.50p <- lapply(1:length(nogap_multi.yr.50), function(i) {
  p <- disaggregate(nogap_multi.yr.50[[i]])
  p@data$birdyear <- p@data$id #clarify variable name
  p@data$id <- rownames(p@data) #create id variable of rownames so data can be merged with fortified polygon
  p
}) #Split 50% core area into seperate polygons

##nonbreeding Overlap----

##BA overlap of 95% KDE
UD_overlap <- lapply(1:length(nogap_multi.yr), function(i){
  kerneloverlaphr(nogap_multi.yr[[i]], meth = "BA", percent = 95, conditional = T)
})

UD_overlap_range  <- data.frame(matrix(ncol = 2, nrow = length(UD_overlap)))
for(i in 1:length(UD_overlap)){
  overlap <- UD_overlap[[i]]
  for(j in 1:length(overlap[1,])){ ###remove overlap within same years
    overlap[j,j] <- NA
  }
  range <-  range(overlap, na.rm = T)
  UD_overlap_range[i,] <- range
}
colnames(UD_overlap_range) <- c("min_overlap", "max_overlap")
UD_overlap_range$id <- names(nogap_multi.yr) ##order in loop
## max is 0.95 (overlap with)

UD_overlap_mn  <- data.frame(matrix(ncol = 1, nrow = length(UD_overlap)))
for(i in 1:length(UD_overlap)){
  overlap <- UD_overlap[[i]]
  for(j in 1:length(overlap[1,])){ ###remove overlap within same years
    overlap[j,j] <- NA
  }
  mean.overlap <-  mean(overlap, na.rm = T)
  UD_overlap_mn[i,] <- mean.overlap
}

colnames(UD_overlap_mn) <- c("mn_overlap")
UD_overlap_mn$id <- names(nogap_multi.yr) ##order in loop

## Gull.glm, add overlap ----
UD_overlap_range <- left_join(UD_overlap_range, UD_overlap_mn)
UD_overlap_range$id <- as.numeric(UD_overlap_range$id)
gull.glm <- left_join(gull.glm, UD_overlap_range)

## Seperate Core Areas ----

## convert polygons into a dataframe for plotting with ggplot
nogap_multi.yr.50p <- lapply(nogap_multi.yr.50p, 
                             function(x) spTransform(x, CRS("+init=epsg:4326")))

nogap_multi.yr.50df <- lapply(1:length(nogap_multi.yr.50p), function(i) {
  df <- fortify(nogap_multi.yr.50p[[i]]) ##convert to dataframe
  merge(df, nogap_multi.yr.50p[[i]]@data, by = "id") #add original data
}) ## convert polygons into a dataframe for plotting with ggplot, each id is one element of list. 
nogap_multi.yr.50df <- mapply(cbind, nogap_multi.yr.50df, "bird_id"=id.v, SIMPLIFY=F)

##Fragmented winter polygons 

ggplot(map, aes(long, lat, group = group)) + geom_polygon(fill = "white", col = "black") +
  geom_polygon(data = filter(nogap_multi.yr.50df[[which(id.v == 395)]], birdyear == 0), aes(col = birdyear, fill = id, group = id), alpha = .3) + ##birdyear = year, id = polygon number
  coord_fixed(xlim = c(-21, 17.8),
              ylim = c(10, 58.1))

ggplot(map, aes(long, lat, group = group)) + geom_polygon(fill = "white", col = "black") +
  geom_polygon(data = filter(nogap_multi.yr.50df[[which(id.v == 5296)]], birdyear == 0), aes(col = id, fill = id, group = id), alpha = .3) + ##birdyear = year, id = polygon number
  coord_fixed(xlim = c(-21, 17.8),
              ylim = c(10, 58.1))

ggplot(map, aes(long, lat, group = group)) + geom_polygon(fill = "white", col = "black") +
  geom_polygon(data = filter(nogap_multi.yr.50df[[which(id.v == 5296)]], birdyear == 1), aes(col = id, fill = id, group = id), alpha = .3) + ##birdyear = year, id = polygon number
  coord_fixed(xlim = c(-21, 17.8),
              ylim = c(10, 58.1))

ggplot(map, aes(long, lat, group = group)) + geom_polygon(fill = "white", col = "black") +
  geom_polygon(data = filter(nogap_multi.yr.50df[[which(id.v == 5335)]], birdyear == 1), aes(col = id, fill = id, group = id), alpha = .3) + ##birdyear = year, id = polygon number
  coord_fixed(xlim = c(-21, 17.8),
              ylim = c(10, 58.1))

ggplot(map, aes(long, lat, group = group)) + geom_polygon(fill = "white", col = "black") +
  geom_polygon(data = filter(nogap_multi.yr.50df[[which(id.v == 5535)]], birdyear == 1), aes(col = id, fill = id, group = id), alpha = .3) + ##birdyear = year, id = polygon number
  coord_fixed(xlim = c(-21, 17.8),
              ylim = c(10, 58.1))

# The following have fragmented polygons in 1 year that are joined in others. These years need to be merged:
#   
#   395.0: p1 + p2
# 5296.0: p1 + p2 + p3
# 5296.1: p5 + p6 
# 5335.1: p3 + p4
# 5535.1: p3 + p4

##fragmented WA poly
nogap_multi.yr.50df[[which(id.v == 395)]][nogap_multi.yr.50df[[which(id.v == 395)]]$id == 2,]$id <- 1
nogap_multi.yr.50df[[which(id.v == 5296)]][nogap_multi.yr.50df[[which(id.v == 5296)]]$id %in% c(2,3),]$id <- 1
nogap_multi.yr.50df[[which(id.v == 5296)]][nogap_multi.yr.50df[[which(id.v == 5296)]]$id == 6,]$id <- 5
nogap_multi.yr.50df[[which(id.v == 5335)]][nogap_multi.yr.50df[[which(id.v == 5335)]]$id == 4,]$id <- 3
nogap_multi.yr.50df[[which(id.v == 5535)]][nogap_multi.yr.50df[[which(id.v == 5535)]]$id == 4,]$id <- 3

###Create 1 dataframe for all polygons
poly.df <- lapply(1:length(id.v), function(i){
  df <- nogap_multi.yr.50df[[i]]
  df <- df %>% select( polygon = id,id = bird_id, birdyear, lat, lon = long, order) 
})
poly.df <- data.table::rbindlist(poly.df)
poly.df$uniq.p <- paste(poly.df$id, poly.df$birdyear, poly.df$polygon, sep = ".")

##create a list of polygons per core area (instead of a list per year)
poly.p <- lapply(split(poly.df[, c("lat","lon")], poly.df[, "uniq.p"]), Polygon)

##Find centroid of each core area
centroid <- lapply(1:length(poly.p), function(i){
  lat <- poly.p[[i]]@labpt[1]
  lon <- poly.p[[i]]@labpt[2]
  uniq.p <- split(poly.df[, 7], poly.df[, 7])[[i]][1]
  df <- data.frame(lat, lon, uniq.p)
  df
})
centroid <- data.table::rbindlist(centroid)
tmp <- poly.df %>% select(-lat, -lon, -order) %>% distinct()
centroid <-  left_join(centroid, tmp)
#centroid$birdyear <- as.numeric(centroid$birdyear) - 1  #birdyear was renumber starting at 1 instead of 0, so change back to join

centroid <- gull.glm %>% select(id, start_lat, start_lon) %>% distinct() %>% right_join(centroid)

## Define Winter area ----
## distance between colony and polygon centroid
centroid <- mutate(centroid, mig.dist = deg.dist(lon, lat, start_lon, start_lat)) 

#Calculate enter and exit times per polygon
##select only points in the polygons
gc()
remove(poly.pts)
poly.p <- poly.p[centroid$uniq.p] ##orders poly.p
poly.pts <- data.frame()
for (i in 1:length(centroid$uniq.p)) { ##find points within each core area
  ID <- centroid[i,]$id
  by <- centroid[i,]$birdyear
  p <- centroid[i,]$uniq.p
  pts <- nonbreed.nogap[nonbreed.nogap$id == ID & nonbreed.nogap$birdyear == by,] ## select points for the right id and bird year
  pts <- select(pts, id, time, lat, lon, birdyear, dur)
  pol <- poly.p[[p]] ## select the core area polygon
  overlap.log <- point.in.polygon(pts$lon, pts$lat, pol@coords[,2], pol@coords[,1]) ## Find points overlaping the polygon
  overlap.log <- overlap.log == 1 ##1 = overlaping
  overlap <- pts[overlap.log,]   ##overlapping points from the furthest polygon in i individual birdyear
  overlap$poly <- p
  poly.pts <- rbind(poly.pts, overlap)   ####list of all points in polygons. 
}

poly.pts <- ungroup(poly.pts)
poly.pts <- poly.pts %>% group_by(id) %>% arrange(time, .by_group = TRUE) %>% ungroup() #order points by time
poly.pts$p.index <- rep(1:length(rle(poly.pts$poly)$lengths), rle(poly.pts$poly)$lengths) 

#Find first and last point in a polygon
en.ex <- poly.pts %>% group_by(id, birdyear, poly, p.index) %>% summarise(enter = min(time), exit = max(time))
en.ex$time.in.poly <- difftime(en.ex$exit, en.ex$enter, units = "days") ## Find duration in polygon

# 
# Selection of wintering area:
# duration during the winter
winter.time.poly <- en.ex %>% mutate(enter = if_else(yday(enter) > 335 | yday(enter) < 90, yday(enter),
                                                     if_else(yday(enter) > 152, 335, 90)), 
                                     exit = if_else(yday(exit) > 335 | yday(exit) < 90, yday(exit),
                                                    if_else(yday(exit) > 152, 335, 90))) %>%
  mutate(enter = ifelse(enter < 152, enter+365, enter), exit = ifelse(exit < 152, exit +365, exit), wtime = exit-enter) %>% group_by(id, birdyear)%>% filter(wtime == max(wtime))

winter.poly.v <- winter.time.poly$poly

##89 days is the maximum possible duration duration (number of day from Dec 1 - Mar 1)


## furthest polygon. 

far_poly <- centroid %>% #in km
  group_by(id, birdyear) %>% filter(mig.dist == max(mig.dist)) %>% 
  ungroup() %>% select(id, birdyear, polygon, mig.dist, uniq.p)

far.poly.v <-  far_poly$uniq.p

#was the primary used between December - March in all but 2 cases (from same individual). 
#table(far.poly.v %in% winter.poly.v) 

## Based on maps, the resuls of time during winter make more sense - individual spent autumn on coast of portugal (furthest poly), moving to central spain for the winter months. 
## Select centroid from polygons in winter.poly.v, to use as migration distance
## Migration distance ----
winter_poly <- centroid %>% group_by(id, birdyear) %>% filter(uniq.p %in% winter.poly.v) %>% 
  ungroup() %>% select(id, birdyear, polygon, mig.dist, uniq.p)

##Winter arrival & Departure ----
## Select polygon of wintering area
winter.poly.p<- poly.p[winter_poly$uniq.p] ##In most cases,  two polygons in a similar region, so no change in migration distance (but see 4032.1).  

##select only points in wintering area polygon in each year
remove(winter.poly.pts)
winter.poly.pts <- data.frame()
nonbreed.nogap.d <- ungroup(nonbreed.nogap.d)
for (i in 1:length(winter_poly$id)) {
  ID <- winter_poly[i,]$id
  by <- winter_poly[i,]$birdyear
  pts <- nonbreed.nogap.d[nonbreed.nogap.d$id == ID & nonbreed.nogap.d$birdyear == by,]
  pol <- winter.poly.p[[i]]
  overlap.log <- point.in.polygon(pts$lon, pts$lat, pol@coords[,2], pol@coords[,1])
  overlap.log <- overlap.log == 1
  overlap <- pts[overlap.log,]   ##overlapping points from the furthest polygon in i individual birdyear
  winter.poly.pts <- rbind(winter.poly.pts, overlap)
}

## Select first and last point in each winter area (times at which a polygon was either entered or exited)
winter.poly.ee <- winter.poly.pts %>% 
  group_by(id, birdyear) %>% 
  mutate(ee = ifelse(time == min(time), "en", ifelse (time == max(time), "ex", NA))) %>%  
  filter(!is.na(ee)) %>%
  ungroup() 

##winter departure = last point in winter area
ex <- filter(winter.poly.ee, ee == "ex") %>% select(id, birdyear, winter.depart = time)

##Winter arrival = first point in winter area
en <- filter(winter.poly.ee, ee == "en") %>% select(id, birdyear, winter.arrive = time)

winter_poly$birdyear <- as.numeric(as.character(winter_poly$birdyear))
winter_poly<- left_join(winter_poly, ex) %>% left_join(en)
tmp <- winter_poly %>% select(id, birdyear, mig.dist, winter.arrive, winter.depart)

## Gull.glm, add mig.dist & winter time ----
gull.glm <- left_join(gull.glm, tmp)
###more obs in gull.glm than in winter_poly because data with gaps can still be used in analysis of arrival and departure from colony (if mig.dist not included as a fixed factor) ##

## remove arrivals and departures that occur during gap ----
gap <- nonbreed.nogap %>% group_by(id, birdyear) %>% mutate(start.gap = lag(time)) %>% ungroup() %>%
  filter(dur > 60*60*24*2) %>% select(id, birdyear, start.gap, end.gap = time)  

gap <- left_join(gap, select(gull.glm, id, birdyear, winter.arrive, winter.depart, col.depart, col.arrival))
filter(gap, winter.arrive >= start.gap & winter.arrive <= end.gap)
filter(gap, winter.depart >= start.gap & winter.depart <= end.gap)
filter(gap, col.depart >= start.gap & col.depart <= end.gap)
filter(gap, col.arrival >= start.gap & col.arrival <= end.gap)
##5215 by 0  col.depart, 5554 by 0 winter.arrive, 459 by 1 winter arrive
#removal done in .rmd

## Trajectory df ----
#For trajectories: we want all 'travel points' between colony and wintering area (autumn and spring).  
#Core areas enroute will be replaced with their polygon centroid so that the trajectory is smoothed through these areas.
id.v <- unique(nonbreed.nogap.d$id)

###want to use original polygons (not combined WA), because joins are otherwise wonky
nogap_multi.yr.50df.noagg <- lapply(1:length(nogap_multi.yr.50p), function(i) {
  df <- fortify(nogap_multi.yr.50p[[i]]) ##convert to dataframe
  merge(df, nogap_multi.yr.50p[[i]]@data, by = "id") #add original data
}) ## convert polygons into a dataframe for plotting with ggplot, each id is one element of list. 

nogap_multi.yr.50df.noagg <- mapply(cbind, nogap_multi.yr.50df.noagg, "bird_id"=id.v, SIMPLIFY=F)

###Create 1 dataframe for all polygons
poly.df.noagg <- lapply(1:length(id.v), function(i){
  df <- nogap_multi.yr.50df.noagg[[i]]
  df <- df %>% select( polygon = id,id = bird_id, birdyear, lat, lon = long, order) 
})

poly.df.noagg <- data.table::rbindlist(poly.df.noagg)

poly.df.noagg$uniq.p <- paste(poly.df.noagg$id, poly.df.noagg$birdyear, poly.df.noagg$polygon, sep = ".")

poly.p.noagg <- lapply(split(poly.df.noagg[, c("lat","lon")], poly.df.noagg[, "uniq.p"]), Polygon)
centroid.noagg <- lapply(1:length(poly.p.noagg), function(i){
  lat <- poly.p.noagg[[i]]@labpt[1]
  lon <- poly.p.noagg[[i]]@labpt[2]
  uniq.p <- split(poly.df.noagg[, 7], poly.df.noagg[, 7])[[i]][1]
  df <- data.frame(lat, lon, uniq.p)
  df
})
centroid.noagg <- data.table::rbindlist(centroid.noagg)
tmp <- poly.df.noagg %>% select(-lat, -lon, -order) %>% distinct()
centroid.noagg <-  left_join(centroid.noagg, tmp)

centroid.noagg$birdyear <- as.numeric(centroid.noagg$birdyear) #- 1  #birdyear was re-numbered starting at 1 instead of 0, so change back to join

##select only points in the polygons
nonbreed.nogap <- ungroup(nonbreed.nogap) ### use non subsampled for duration calculation. 
gc()
remove(poly.pts.noagg)
poly.p.noagg <- poly.p.noagg[centroid.noagg$uniq.p] ##orders poly.p.noagg 
poly.pts.noagg<- data.frame()
for (i in 1:length(centroid.noagg$uniq.p)) {
  ID <- centroid.noagg[i,]$id
  by <- centroid.noagg[i,]$birdyear
  p <- centroid.noagg[i,]$uniq.p
  pts <- nonbreed.nogap[nonbreed.nogap$id == ID & nonbreed.nogap$birdyear == by,]
  pts <- select(pts, id, time, lat, lon, birdyear, dur)
  pol <- poly.p.noagg[[i]]
  overlap.log <- point.in.polygon(pts$lon, pts$lat, pol@coords[,2], pol@coords[,1])
  overlap.log <- overlap.log == 1
  overlap <- pts[overlap.log,]   ##overlapping points from the furthest polygon in i individual birdyear
  overlap$poly <- p
  poly.pts.noagg <- rbind(poly.pts.noagg, overlap)   ####list of all points in polygons. 
}
poly.pts.noagg <- ungroup(poly.pts.noagg)
poly.pts.noagg <- poly.pts.noagg %>% group_by(id) %>% arrange(time, .by_group = TRUE) %>% ungroup()
poly.pts.noagg$p.index <- rep(1:length(rle(poly.pts.noagg$poly)$lengths), rle(poly.pts.noagg$poly)$lengths)

#find first and last point in core area
en.ex.noagg <- poly.pts.noagg %>% group_by(id, birdyear, poly, p.index) %>% summarise(enter = min(time), exit = max(time))
en.ex.noagg$time.in.poly <- difftime(en.ex.noagg$exit, en.ex.noagg$enter, units = "days")

##Replace points during time bird was in core area, replace with centroid.
#this stops distance from accumulating during central point foraging, as well as cleans up the actual migratory track. 
gull.glm$id_birdyear <- paste(gull.glm$id, gull.glm$birdyear, sep = ".")

smooth_pts <- filter(nonbreed.nogap, id_birdyear %in% gull.glm$id_birdyear) %>% ##hourly
  full_join(en.ex.noagg) %>% 
  filter(time > enter & time < exit) %>%
  select(id, birdyear, time, poly)

tmp <- centroid.noagg %>% select(c.lat = lat, c.lon = lon, uniq.p)
nonbreed.UDcentroid <- left_join(nonbreed.nogap, smooth_pts) %>% ##DF of points, core areas replaced by centroid
  left_join(tmp, by = c("poly" = "uniq.p")) %>% 
  mutate(lat = if_else(is.na(poly), lat, c.lat),
         lon = if_else(is.na(poly), lon, c.lon))

##Find point between colony exit and winter entrance (autumn), and winter exit and colony entrance (spring)
all_traj <- nonbreed.UDcentroid %>% left_join(gull.glm) %>% 
  mutate(direction = if_else(time >= col.depart & time <= winter.arrive, "Autumn", 
                             if_else(time >= winter.depart & time <= col.arrival, "Spring", "NA"))) %>%
  filter(!direction=="NA") %>% arrange(time) %>% 
  select(id, lon, lat, birdyear, id_birdyear, direction, time, poly) %>% 
  distinct()

##Only use trajectories with 1 point per day during travel periods (outside of core areas)
all_traj <- all_traj %>% group_by(id_birdyear, direction) %>% 
  mutate(dur = as.numeric(as.character(difftime(time, lag(time), units = "secs")))) %>% ungroup()
id.w.gap <- all_traj %>% filter(is.na(poly) & dur >= 3600*24) %>% select(id_birdyear, direction) %>% #points outside polygons have no polygon name (NA)
  distinct() %>% mutate(gap = T) #if there is a dur > 24 h, lable this traj with gap = T
all_traj <- all_traj %>% left_join(id.w.gap) %>% filter(is.na(gap)) %>% select(-gap) #remove points from traj where gap = T
pair.miss <- all_traj %>% select(id, birdyear,direction) %>% distinct() %>% ## Find traj left with no pair
  group_by(id, direction) %>% summarise(n=n()) %>% filter(n==1) %>% 
  select(id, direction) %>% mutate(gap = T) #an additional 18 rows lost
all_traj <- all_traj %>% left_join(pair.miss) %>% filter(is.na(gap)) %>% select(-gap) # Remove traj with no pair

saveRDS(all_traj, "all_traj.RDS")

## N core areas & core area overlap ----
##core areas with no overlap

### find winter areas with no overlap

##Fragmented winter area p
nogap_multi.yr.50p[[which(id.v == 395)]]@data[2,]$id <- "1"

nogap_multi.yr.50p[[which(id.v == 5296)]]@data[c(2,3),]$id <- "1"
nogap_multi.yr.50p[[which(id.v == 5296)]]@data[6,]$id <- "5"

nogap_multi.yr.50p[[which(id.v == 5335)]]@data[4,]$id <- "3"

nogap_multi.yr.50p[[which(id.v == 5535)]]@data[4,]$id <- "3"

no.wa.overlap <- vector()
for(i in 1:length(nogap_multi.yr.50p)){
  id.i <- id.v[i] #vector of bird ids
  sp <- nogap_multi.yr.50p[[i]] ## spatial polygon list of all core areas for a bird id
  p <- nogap_multi.yr.50p[[i]]@data ### all core areas across years
  p$poly <- paste(rep(id.i, times = length(p$id)) , p$birdyear, p$id, sep = ".")
  wa.id <- filter(winter.time.poly, id == id.i) %>% pull(poly) ## list of wa polys for that id
  
  for(j in 1:length(wa.id)){ ##for each wa polygon
  p1 <-  sp[which(p$poly == wa.id[j]),]
  other.wa <- wa.id[which(wa.id != wa.id[j])]
  
  no.overlap.v <- logical(length(other.wa))
  for(k in 1:length(other.wa)){ ## see if it overlaps with other wa polygons
    p2 <-sp[which(p$poly == other.wa[k]),]
    int <- gIntersection(p1,p2)
    no.overlap.v[k] <- is.null(int) 
    ## if there is 1 T, then wa in 1 yer doesn't overlap with wa in another year
  }
  if(!all(!no.overlap.v)){ ## if there is 1 or more True values (i.e. a wa that doesn't overlap with a wa in another year)
    no.wa.overlap <- rbind(no.wa.overlap, id.i)
 
  }
  }
}
no.wa.overlap <- unique(no.wa.overlap)

###534 spent year 3 in uk only

ggplot(map, aes(long, lat, group = group)) + geom_polygon(fill = "white", col = "black") +
  geom_polygon(data = filter(nogap_multi.yr.50df[[which(id.v == 606)]], birdyear == 1), aes(col = id, fill = id, group = id), alpha = .3) + ##birdyear = year, id = polygon number
  coord_fixed(xlim = c(-21, 17.8),
              ylim = c(10, 58.1))

no.so.overlap <- data.frame()
for (i in 1:length(nogap_multi.yr.50p)){
  id.i <- id.v[i] #vector of bird ids
  sp <- nogap_multi.yr.50p[[i]] ## spatial polygon list of all core areas for a bird id
  p <- nogap_multi.yr.50p[[i]]@data ### all core areas across years
  p$poly <- paste(rep(id.i, times = length(p$id)) , p$birdyear, p$id, sep = ".")
  wa.id <- filter(winter.time.poly, id == id.i) %>% pull(poly) ## list of wa polys for that id
  so.poly <- p[which(!(p$poly %in% wa.id)),]
  if(length(so.poly$id)==0) next ## if no stopovers, skip to next id
  for(j in 1:length(so.poly$id)){
    p.id <- so.poly[j,1]
    by <- so.poly[j,3]
    ## list all polygons in different years
    id.pair <- as.numeric(filter(p, birdyear != by) %>% pull(id))
    no.overlap.v <- logical(length(id.pair))
    p1<-sp[sp$id == p.id,]
    for(k in 1:length(id.pair)){
      p2<-sp[sp$id == id.pair[k],]
      int <- gIntersection(p1,p2)
      no.overlap.v[k] <- is.null(int)
    }
      if(all(no.overlap.v)){
        x<-data.frame(id.i, by, p.id)
        no.so.overlap <- rbind(no.so.overlap, x)  ##core area with no overlap in any year. 
      }
    }
  }
 
no.overlap <- select(no.so.overlap, id=id.i, birdyear=by) %>% distinct() %>% mutate(p.no.overlap = T)

## how many polygons per birdyear?
remove(n.winter.area)
n.winter.area <- data.frame()
for (i in 1:length(nogap_multi.yr.50p)){
  year <- unique(nogap_multi.yr.50p[[i]]@data$birdyear)  #vector of year ids
  id <- id.v[i] #vector of bird ids
  p <- nogap_multi.yr.50p[[i]]@data
  p <- select(p, id, birdyear) %>% distinct()
  n <- as.numeric(rle(as.numeric(p$birdyear))$lengths)
  df <- data.frame(id, year, n)
  n.winter.area <- rbind(n.winter.area, df)
}

# 
n.winter.area$year <- as.numeric(as.character(n.winter.area$year))
names(n.winter.area) <- c("id","year","n.winter.p")

## Gull.glm, add n core areas ----
gull.glm <- left_join(gull.glm, n.winter.area, by = c("id", "birdyear" = "year"))
no.overlap <- no.overlap %>% mutate(birdyear = as.numeric(as.character(birdyear))) 
gull.glm <- left_join(gull.glm, no.overlap)
gull.glm$change.wa <- ifelse(gull.glm$id %in% no.wa.overlap[,1], T, F)
## tracking summary ----

#migration duration
mig.dur <- nonbreed.nogap.d %>% group_by(id, birdyear) %>% 
  summarise(min = min(time), max = max(time), mig.dur = as.numeric(difftime(max, min, units = "days")))

#days with fix by birdyear
days.w.fix <- nonbreed.nogap.d %>%  mutate(yday = yday(time)) %>% select(id, birdyear, yday) %>% distinct() %>%
  group_by(id, birdyear)  %>% summarise(days.with.fix = n())
days.w.fix <- days.w.fix %>% left_join(mig.dur) %>% mutate(p.days.w.fix = days.with.fix/mig.dur) %>% ungroup() 
days.w.fix %>%  summarise(mean = mean(p.days.w.fix), min = min(p.days.w.fix), 
                          max = max(p.days.w.fix), median = median(p.days.w.fix)) 

days.w.fix <- nonbreed.nogap.d %>% group_by(id, birdyear) %>% 
  summarise(max.gap = max(as.numeric(difftime(time, lag(time), units = "days")), na.rm = T)) %>% left_join(days.w.fix) %>%
  select(id, birdyear, mig.dur, p.days.w.fix, max.gap) %>% ungroup()
days.w.fix %>%  summarise(mean = mean(max.gap), min = min(max.gap), 
                          max = max(max.gap), median = median(max.gap)) 

## WA site fidelity ----
 library(raster)
 nonbreed <- readRDS("nonbreed.RDS")
 nonbreed <- ungroup(nonbreed)
# 
 wa.pts <- data.frame()
 for (i in 1:length(gull.glm$id)) {
  s <- gull.glm[i,]$winter.arrive
  e <- gull.glm[i,]$winter.depart
  id.i <- gull.glm[i,]$id
  by.i <- gull.glm[i,]$birdyear
  wa <- filter(nonbreed,id == id.i& birdyear == by.i& time >= s & time <= e)
  wa.pts <- rbind(wa.pts, wa)
}
rm(nonbreed)

med.step <- quantile(wa.pts$dist, c(.5)) #50% = 43m, 90% = 1984m, 95% = 4602m

## wa points per id and by
id.v <- unique(wa.pts$id_birdyear)
wa.laea.l <- lapply(id.v, function(x){
  mkSpatial(wa.pts[wa.pts$id_birdyear == x,])
})
## all winter points of 1 individual
idu.v <- unique(wa.pts$id)
wau.laea.l <- lapply(idu.v, function(idu.v) mkSpatial(wa.pts[wa.pts$id == idu.v,]))

#make ltraj
wau.traj.l <- lapply(wau.laea.l, function(wa.laea){
  gc()
  adehabitatLT::as.ltraj(as.data.frame(wa.laea@coords), date = wa.laea@data$time,
                         id = wa.laea@data$birdyear,
                         slsp = "missing", proj4string = CRS(proj4string(wa.laea)))})
saveRDS(wau.laea.l, "wau.laea.l.RDS")
saveRDS(wau.traj.l, "wau.traj.l.RDS")
save.image("ln831.RData")
rm(list = ls())
wau.laea.l <- readRDS("wau.laea.l.RDS")
wau.traj.l <- readRDS("wau.traj.l.RDS")
laea.proj <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +towgs84=0,0,0,0,0,0,0 +units=m"

mkgrid <- function(lat=NULL, lon = NULL, resolution=10000, buffer = 100000, 
                   projection = laea.proj  ){
  #lat is a vector of latitudes
  #lon is a vector of longitudes
  #resolution is the cell size of the raster. 
  #1 = ~ 133X165 km
  #.1 = ~ 13.9 x 17.1 km
  #.01 ~ 1.39 x 1.72 km     .... in WSG84
  #buffer adds space beyond the points in the grid. a buffer of 1 should equal ~ 111 km
  #projection: projection of the data.
  
  ##projection = laea, measured in m
  #resolution 1 000 = 1 km2...
  
  
  xmax <- max(lon)
  xmin <- min(lon)
  ymax <- max(lat)
  ymin <- min(lat)
  extent<- matrix(c(xmin - buffer, xmax + buffer, ymin - buffer, ymax + buffer), ncol=2)
  extent <- SpatialPoints(coords = extent, proj4string = CRS(projection))        # the data is in WSG84 lat-lon
  rast<- ascgen(extent, cellsize=resolution)
}

##make grid
brbu.grid <- lapply(wau.laea.l[1:41], function(wa.grid) {
  gc()
  mkgrid(lat = wa.grid@coords[,2], lon = wa.grid@coords[,1],
         resolution = 500, buffer = 2500)
})

## diffusion coef
Du.l <- lapply(wau.traj.l[1:41], function(traj) BRB.D(traj, Tmax = 3*3600, Lmin= 20)) ## one bird has a point every 3 hours

wa_overlap.1 <- vector(mode = "list", length = 41)
for(i in 1:18){
  brb.i <- BRB(wau.traj.l[[i]],Du.l[[i]], filtershort = F, Tmax = 3*3600, Lmin= 20, hmin =150, type = "UD", grid = brbu.grid[[i]])
  gc()
  wa_overlap.1[[i]] <- kerneloverlaphr(brb.i, meth = "BA", percent = 95, conditional = T)
}

for(i in 20:41){
  brb.i <- BRB(wau.traj.l[[i]],Du.l[[i]], filtershort = F, Tmax = 3*3600, Lmin= 20, hmin =150, type = "UD", grid = brbu.grid[[i]])
  gc()
  wa_overlap.1[[i]] <- kerneloverlaphr(brb.i, meth = "BA", percent = 95, conditional = T)
}

saveRDS(wa_overlap.1, "wa_overlap.1.RDS")
##individual brbs
##make grid
brbu.grid <- lapply(wau.laea.l[42:82], function(wa.grid) {
  gc()
  mkgrid(lat = wa.grid@coords[,2], lon = wa.grid@coords[,1],
         resolution = 500, buffer = 2500)
})

## diffusion coef
Du.l <- lapply(wau.traj.l[42:82], function(traj) BRB.D(traj, Tmax = 3*3600, Lmin= 20)) ## one bird has a point every 3 hours

wa_overlap.2 <- vector(mode = "list", length = 41)
for(i in 1:41){
  brb.i <- BRB(wau.traj.l[[41+i]],Du.l[[i]], filtershort = F, Tmax = 3*3600, Lmin= 20, hmin =150, type = "UD", grid = brbu.grid[[i]])
  gc()
  wa_overlap.2[[i]] <- kerneloverlaphr(brb.i, meth = "BA", percent = 95, conditional = T)
}

##save environment  (ln861_2506.RData)
saveRDS(wa_overlap.2, "wa_overlap.2.RDS")
rm(list = ls())
wau.laea.l <- readRDS("wau.laea.l.RDS")
wau.traj.l <- readRDS("wau.traj.l.RDS")
laea.proj <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +towgs84=0,0,0,0,0,0,0 +units=m"

mkgrid <- function(lat=NULL, lon = NULL, resolution=10000, buffer = 100000, 
                   projection = laea.proj  ){
  #lat is a vector of latitudes
  #lon is a vector of longitudes
  #resolution is the cell size of the raster. 
  #1 = ~ 133X165 km
  #.1 = ~ 13.9 x 17.1 km
  #.01 ~ 1.39 x 1.72 km     .... in WSG84
  #buffer adds space beyond the points in the grid. a buffer of 1 should equal ~ 111 km
  #projection: projection of the data.
  
  ##projection = laea, measured in m
  #resolution 1 000 = 1 km2...
  
  
  xmax <- max(lon)
  xmin <- min(lon)
  ymax <- max(lat)
  ymin <- min(lat)
  extent<- matrix(c(xmin - buffer, xmax + buffer, ymin - buffer, ymax + buffer), ncol=2)
  extent <- SpatialPoints(coords = extent, proj4string = CRS(projection))        # the data is in WSG84 lat-lon
  rast<- ascgen(extent, cellsize=resolution)
}
##skipped 19 b/c takes too much memory!!!
brbu.grid <-  mkgrid(lat = wau.laea.l[[19]]@coords[,2], lon = wau.laea.l[[19]]@coords[,1],
                     resolution = 500, buffer = 500)
Du.l <-  BRB.D(wau.traj.l[[19]], Tmax = 3*3600, Lmin= 20) ## one bird has a point every 3 hours

gc()
x <- BRB(wau.traj.l[[19]],Du.l, filtershort = F, Tmax = 3*3600, Lmin= 20, hmin =150, type = "UD", grid = brbu.grid)
gc()

##null slot in overlap.1 b/c missing 19
##If memory limit, open kerneloverlaphr function and run through loop manually. 
#wa_overlap.19 <- kerneloverlaphr(x, meth = "BA", percent = 95, conditional = T)
vol <- getvolumeUD(x)
x <- lapply(x, function(y) {
  coo <- coordinates(y)
  y[order(coo[, 1], coo[, 2]), ]
})
vol <- lapply(vol, function(y) {
  coo <- coordinates(y)
  y[order(coo[, 1], coo[, 2]), ]
})
gp <- gridparameters(vol[[1]])
res <- matrix(0, ncol = length(x), nrow = length(x))
for (i in 1:length(x)) {
  for (j in 1:i) {
      vi <- x[[i]][[1]]
      vj <- x[[j]][[1]]
      ai <- vol[[i]][[1]]
      aj <- vol[[j]][[1]]
      ai[ai <= percent] <- 1
      ai[ai > percent] <- 0
      aj[aj <= percent] <- 1
      aj[aj > percent] <- 0
     
        vi <- vi * ai
        vj <- vj * aj
        res[j, i] <- res[i, j] <- sum(sqrt(vi) * sqrt(vj)) * 
          (gp[1, 2]^2)

  }
}
  
  rownames(res) <- names(x)
  colnames(res) <- names(x)

  wa_overlap.19 <- res

wa_overlap.1 <- readRDS("wa_overlap.1.RDS")
wa_overlap.1[[19]] <- wa_overlap.19
wa_overlap.2 <- readRDS("wa_overlap.2.RDS")
wa_overlap <- c(wa_overlap.1,wa_overlap.2)
saveRDS(wa_overlap, "wa_overlap_3035.RDS")

load("ln831.RData")
wa_overlap <- readRDS("wa_overlap_3035.RDS")

sf_range  <- data.frame(matrix(ncol = 4, nrow = length(idu.v)))
for(i in 1:length(wa_overlap)){
  wa_overlap.i <- wa_overlap[[i]]
  ##get range
  for(j in 1:length(wa_overlap.i[1,])){ ###remove overlap within same years
    wa_overlap.i[j,j] <- NA
  }
  range <-  range(wa_overlap.i, na.rm = T)
  sf_range[i,c(1,2)] <- range
  sf_range[i,3] <- mean(wa_overlap.i, na.rm = T)
  sf_range[i,4] <- idu.v[i]
}
colnames(sf_range) <- c("sf_min_overlap", "sf_max_overlap", "sf_mn_overlap", "id")
saveRDS(sf_range, "wa_overlap_range_3035.RDS")
gull.glm <- left_join(gull.glm, sf_range)

## Gull.glm export, add days with fix ----
gull.glm <- left_join(gull.glm, days.w.fix) %>% arrange(mig.dist)
saveRDS(gull.glm, "gull.glm.RDS")
save.image("ln964.RData")
## Overlap of random pairs ----
##### polygon overlap for randomization tests
#####Find pairings

start <- nonbreed.nogap.d %>% left_join(gull.glm) %>%
  select(id, id_birdyear, birdyear, slat = start_lat, slon = start_lon) %>% distinct()
end <- centroid %>% group_by(id, birdyear) %>% filter(uniq.p %in% winter.poly.v) %>% 
  ungroup() %>% select(id, birdyear, elat = lat, elon =lon) %>% 
  mutate(id_birdyear = paste(id, birdyear, sep = "."), birdyear = as.numeric(as.character(birdyear)))

##find distance between each start point, if id != id.  grouped within direction. 
focal <- full_join(start, end)
pair <- focal %>% rename(p.id = id, p.slon = slon, p.slat = slat, p.birdyear = birdyear,
                         p.id_birdyear = id_birdyear, p.elon = elon, p.elat = elat)
rand.pair.ol <- expand.grid(1:length(focal$id), 1:length(pair$p.id))

r.pair.ol <- data.frame()
for (i in seq_along(rand.pair.ol$Var1)){
  x <-  cbind(focal[rand.pair.ol[i,1],], pair[rand.pair.ol[i,2],])
  if(x$id == x$p.id) {
    next
  }
  x$s.dist <- deg.dist(x$slon, x$slat, x$p.slon, x$p.slat)
  x$e.dist <- deg.dist(x$elon, x$elat, x$p.elon, x$p.elat)
  if(x$s.dist > 250 | x$e.dist > 250){  ##start and end need to be within 250k
    next
  }
  x <- select(x, id, birdyear, id_birdyear, p.id, p.birdyear, p.id_birdyear)
  r.pair.ol <- rbind(r.pair.ol, x)
}

## Between individual UD pairs
r.pair.ol <- r.pair.ol[!duplicated(t(apply(r.pair.ol, 1, sort))), ]


r.pair.ol.norep <- r.pair.ol %>% select(-id_birdyear, -p.id_birdyear) %>%
  group_by(id, p.id) %>% ##for each unique id pair & direction
  sample_n(1) %>% ##randomly select 1 birdyear route per id
  ungroup() %>% mutate(id_birdyear = paste(id, birdyear, sep = "."), 
                       p.id_birdyear = paste(p.id, p.birdyear, sep = "."))

###calculate between id overlaps
#readRDS("nogap_multi.yr.RDS")
by.poly <- unlist(nogap_multi.yr)
by.poly.v <- names(by.poly)
rand.UD_overlap <- lapply(1:length(r.pair.ol.norep$id), function(i){
  f.poly <- by.poly[[which(by.poly.v == r.pair.ol.norep[i,]$id_birdyear)]]
  p.poly <- by.poly[[which(by.poly.v == r.pair.ol.norep[i,]$p.id_birdyear)]]
  ol <- list(f.poly, p.poly)
  names(ol) <- c(r.pair.ol.norep[i,]$id_birdyear,r.pair.ol.norep[i,]$p.id_birdyear)
  class(ol) <- "estUDm" 
  kerneloverlaphr(ol, meth = "BA", percent = 95, conditional = T)
})

saveRDS(rand.UD_overlap,"rand_overlap.RDS")
saveRDS(r.pair.ol.norep, "r.pair.ol.RDS")
save.image("ln1021.RData")
## Migration route variation ----

### trajectory averaging
##based on freeman et al. 

##df of single direction trajectories, multiples identified by birdyear
## df includes  id, birdyear, lon, lat, direction, date_time
## points within wintering area polygons (i.e. stopovers) were replaced with the polygon centroid to aid in  
#spacing points equally along the trajectory 

all_traj <- ungroup(all_traj) %>% arrange(id, time)

n_check <- all_traj %>% select(-time, -dur) %>% distinct() %>% arrange(id_birdyear)

##traj needs to be greater than 1 point
table(table((n_check %>% filter(direction == "Spring") %>% select(id_birdyear, direction, lon, lat) %>% distinct())$id_birdyear)<=1)
tmp <- n_check %>% filter(direction == "Spring") %>% select(id_birdyear, direction, lon, lat) %>% distinct()
nopath <- rle(tmp$id_birdyear)$values[which(rle(tmp$id_birdyear)$lengths <=1)]
all_traj <- filter(all_traj, !(id_birdyear %in% nopath & direction == "Spring"))
table(table((n_check %>% filter(direction == "Autumn") %>% select(id_birdyear, direction, lon, lat) %>% distinct())$id_birdyear)<=1)
tmp <- n_check %>% filter(direction == "Autumn") %>% select(id_birdyear, direction, lon, lat) %>% distinct()
nopath <- rle(tmp$id_birdyear)$values[which(rle(tmp$id_birdyear)$lengths <=1)]
all_traj <- filter(all_traj, !(id_birdyear %in% nopath & direction == "Autumn"))

## the three removed all still had multiple traj

pair.miss <- all_traj %>% select(id, id_birdyear, direction) %>% distinct()
pair.miss <- pair.miss %>% group_by(id, direction) %>% summarise(n.years = n()) %>% filter(n.years <=1)

all_traj <- all_traj %>% left_join(pair.miss) %>% filter(is.na(n.years)) %>% select(-n.years)

id.dir <- all_traj %>%  select(id, direction) %>% distinct() ##all traj
id.dir <- filter(id.dir, ! id %in% c(1402, 606, 5027, 5593, 534)) ## 534 - wa poly doesn't overlap, but explored that area

mn_traj <- data.frame(matrix(ncol = 5, nrow = 0))
names(mn_traj) <- c("mn.lon", "mn.lat", "within.var", "id", "direction")

for(a in 1:length(id.dir$id)){
  # for(a in a:length(id.dir$id)){  
  traj <- filter(all_traj, id == id.dir[a,]$id & direction == id.dir[a,]$direction)
  
  ###
  ###This creates a list of equally spaced points (currently n=11, specified in t1_eq), with each trajectory having it's own list element. 
  ntraj <- unique(traj$birdyear)  #trajectory id
  equi_n <- 500  ### number of points along the trajectory
  traj_l <- list()
  name <- c()
  for(i in 1:length(ntraj)){
    t1 <- filter(traj, birdyear == ntraj[i]) %>% 
      select(lon, lat) %>% distinct()
    if (length(t1$lon)<=10){ next }
    
    t1_l <-  SpatialLines(list(Lines(list(Line(cbind(t1$lon, t1$lat))), "id")))
    t1_eq <- spsample(t1_l, equi_n, type = "regular") ## 250 point equally placed along length of line
    traj_l[[i]] <- t1_eq@coords
    name <- c(name, ntraj[[i]])
    
  }
  if(length(name) <= 1 | is.null(name)) {next}
  if(length(name) != length(traj_l)) { traj_l <- traj_l[-which(sapply(traj_l, is.null))]}
  names(traj_l) <- name
  
  ###
  ###Create starting 'thread' for mean trajectory between start and end midpoints
  
  start <- cbind(x = mean(sapply(traj_l, '[[',1,1)), y = mean(sapply(traj_l, '[[',1,2))) ##mean of first elements in each list
  end <- cbind(x = mean(sapply(traj_l, '[[',equi_n,1)), y = mean(sapply(traj_l, '[[',equi_n,2)))
  
  ###
  ###Create starting 'thread' from mean trajectory points
  tmean <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(tmean) <- c("x","y")
  
  for(i in 1:equi_n){
    tmean[i,1] = mean(sapply(traj_l, '[[',i,1))
    tmean[i,2] = mean(sapply(traj_l, '[[',i,2))
  }
  
  
  ###
  ##Averaging
  
  iter_n <- 0
  repeat{
    iter_n <- iter_n + 1
    
    for(i in seq_along(tmean$x)){
      pmean <- as.numeric(tmean[i,]) ## select each  point
      
      nn <- data.frame()
      for(j in seq_along(traj_l)){
        x <- traj_l[[j]]
        x <- matrix(x[x[,2] < pmean[2]+1 & x[,2]  > pmean[2]-1,], ncol = 2)
        if (length(x) == 0) { x <- traj_l[[j]]}
        dist <- NA
        for(k in seq_along(x[,1])){
          point2 <- as.numeric(x[k,])
          dist[k] <- deg.dist(pmean[1], pmean[2], point2[1], point2[2])
        }
        nearest_neighbour1 <- x[dist == min(dist),]
        
        nn <- rbind(nn,nearest_neighbour1)
      }
      
      tmean[i,1] <- mean(nn[,1])
      tmean[i,2] <- mean(nn[,2])
    }
    
    tmean <- unique(tmean)
    tmean <- rbind(start, tmean, end)
    dist <- pt2pt.distance(tmean$y, tmean$x) ## dist in m
    dist_log <- dist > 25000 #25 km
    dist_log[1] <- FALSE
    
    if(iter_n == 100){      ###this will stop repeating after 100 iterations. option 2 is to stop repeating if distinces between tmean are all less than e.g. 25
      break
    }
    
    tmean1 <- data.frame()
    for(i in 1:length(dist_log)){
      if(dist_log[i]== T){
        newpt <- new_point(tmean[i-1,], tmean[i,], di = norm_vec(tmean[i-1,]- tmean[i,])/2) #new point at midpoint 
        tmean1 <- rbind(tmean1, newpt)  ## order matters
      }
      tmean1 <- rbind(tmean1, tmean[i,])
    }
    tmean <- tmean1
    
  }
  
  tmean <-  SpatialLines(list(Lines(list(Line(cbind(tmean$x, tmean$y))), "id")))
  tmean <- spsample(tmean, equi_n, type = "regular") ## 500 point equally placed along length of line
  tmean <- tmean@coords
  

  traj_l <- lapply(traj_l, FUN = data.frame)
  for(i in 1:length(traj_l)){
    names(traj_l[[i]]) <- c("lon", "lat")
  }

    nn_dist <- NA
  for(l in seq_along(traj_l)){
    x <- traj_l[[l]]
    for(i in seq_along(tmean[,1])){
      point <- tmean[i,]
      dist <- NA
      for(j in seq_along(x[,1])){
        point2 <- x[j,]
        dist[j] <- deg.dist(point[1], point[2], point2[1,1], point2[1,2])
      }
      nn_dist[i] <- min(dist)
    }
    tmean <- cbind(tmean, nn_dist)
  }
  
  tmean <- data.frame(tmean)
  tmean$within.var <- rowSums(tmean[,3:ncol(tmean)])/(ncol(tmean)-2-1)  ##var = sum of (x-x.mn)/n-1.  here n = num traj = (ncol - 2(=tmean coords))
  tmean <- tmean[,c(1:2,ncol(tmean))]
  names(tmean) <- c("mn.lon", "mn.lat", "within.var")
  tmean$id <- id.dir[[a, "id"]]
  tmean$direction <- id.dir[[a, "direction"]]
  mn_traj <- rbind(mn_traj, tmean)

  
}  
saveRDS(mn_traj, "mn_traj_limitlat.RDS")


## Route variation of random pairs ----
#####Pairing of individuals with close start and end points

start <- all_traj %>% group_by(id, id_birdyear, direction) %>% filter(row_number()==1) %>% 
  ungroup() %>% arrange(direction, id, id_birdyear) %>% 
  rename(slon = lon, slat = lat) %>% select(-time, -poly, -dur)
end <- all_traj %>% group_by(id, id_birdyear, direction) %>% filter(row_number()== n()) %>% 
  ungroup() %>% arrange(direction, id, id_birdyear)%>% rename(elon = lon, elat = lat) %>% 
  select(-time, -poly, -dur)
x <- all_traj %>% group_by(id, id_birdyear, direction) %>% summarise(n= n()) %>% filter(n <= 10)

all_traj <- left_join(all_traj, x) %>% filter(is.na(n)) %>% select(-n)
start <-left_join(start, x) %>% filter(is.na(n)) %>% select(-n)
end <-left_join(end, x) %>% filter(is.na(n)) %>% select(-n)
##find distance between each start point, if id != id.  grouped within direction. 
focal <- full_join(start, end) #start and end points for each bird year
pair <- focal %>% rename(p.id = id, p.slon = slon, p.slat = slat, p.birdyear = birdyear, #copy of focal, renamed with pair
                         p.id_birdyear = id_birdyear, p.direction = direction,
                         p.elon = elon, p.elat = elat)
rand.pair <- expand.grid(1:length(focal$id), 1:length(pair$p.id)) ## create df size of every possible id_birdyear pair

r.pair <- data.frame()
for (i in seq_along(rand.pair$Var1)){
  x <-  cbind(focal[rand.pair[i,1],], pair[rand.pair[i,2],])
  if(x$id == x$p.id | x$direction != x$p.direction) {
    next
  }
  x$s.dist <- deg.dist(x$slon, x$slat, x$p.slon, x$p.slat)
  x$e.dist <- deg.dist(x$elon, x$elat, x$p.elon, x$p.elat)
  if(x$s.dist > 250 | x$e.dist > 250){ next }
  x <- select(x, id, birdyear, p.id, p.birdyear, direction)
  r.pair<- rbind(r.pair,x)
}

r.pair <- r.pair[!duplicated(t(apply(r.pair, 1, sort))),]
###for each id pair, only one bird year should be used

set.seed(15)
r.pair.norep <- r.pair %>% group_by(id, p.id, direction) %>% ##for each unique id pair & direction
  sample_n(1) %>% ##randomly select 1 birdyear route per id
  ungroup()
rm(r.pair)
## calculate average trajectory and variance from paired roots. 


rand_traj <- data.frame(matrix(ncol = 8, nrow = 0))
names(rand_traj) <- c("mn.lon", "mn.lat", "within.var", "id","birdyear", "p.id", "p.birdyear", "direction")
for(a in seq_along(r.pair.norep$id)){
  # for(a in a:length(id.dir$id)){
  traj <- filter(all_traj, id == r.pair.norep[a,]$id & birdyear == r.pair.norep[a,]$birdyear & 
                   direction == r.pair.norep[a,]$direction)
  
  p.traj <- filter(all_traj, id == r.pair.norep[a,]$p.id & birdyear == r.pair.norep[a,]$p.birdyear & direction == r.pair.norep[a,]$direction)
  
  if(length(traj$id) == 0 | length(p.traj$id) == 0) {next}
  
  ###
  ###This creates a list of equally spaced points (currently n=11, specified in t1_eq), with each trajectory having it's own list element.

  equi_n <- 500  ### number of points along the trajectory
  traj_l <- list()
  
  t1 <- traj %>% select(lon, lat) %>% distinct()
  t1_l <-  SpatialLines(list(Lines(list(Line(cbind(t1$lon, t1$lat))), "id")))
  t1_eq <- spsample(t1_l, equi_n, type = "regular") ## 250 point equally placed along length of line
  traj_l[[1]] <- t1_eq@coords
  
  t1 <- p.traj %>% select(lon, lat) %>% distinct()
  t1_l <-  SpatialLines(list(Lines(list(Line(cbind(t1$lon, t1$lat))), "id")))
  t1_eq <- spsample(t1_l, equi_n, type = "regular") ## 250 point equally placed along length of line
  traj_l[[2]] <- t1_eq@coords
  
  names(traj_l) <- c(1,2)
  
  ###
  ###Create starting 'thread' for mean trajectory between start and end midpoints
  
  start <- cbind(x = mean(sapply(traj_l, '[[',1,1)), y = mean(sapply(traj_l, '[[',1,2))) ##mean of first elements in each list
  end <- cbind(x = mean(sapply(traj_l, '[[',equi_n,1)), y = mean(sapply(traj_l, '[[',equi_n,2)))
  
  ###
  ###Create starting 'thread' from mean trajectory points
  tmean <- data.frame(matrix(ncol = 2, nrow = length(1:equi_n)))
  colnames(tmean) <- c("x","y")
  
  for(i in 1:equi_n){
    tmean[i,1] = mean(sapply(traj_l, '[[',i,1))
    tmean[i,2] = mean(sapply(traj_l, '[[',i,2))
  }
  
  ###
  ##Averaging
    iter_n <- 0
  repeat{
    iter_n <- iter_n + 1
    
    for(i in seq_along(tmean$x)){
      pmean <- as.numeric(tmean[i,]) ## select each  point
      
      nn <- data.frame()
      for(j in seq_along(traj_l)){
        x <- traj_l[[j]]
        x <- matrix(x[x[,2] < pmean[2]+1 & x[,2]  > pmean[2]-1,], ncol = 2)
        if (length(x) == 0 ) { x <- traj_l[[j]]}
        dist <- NA
        for(k in seq_along(x[,1])){
          point2 <- as.numeric(x[k,])
          dist[k] <- deg.dist(pmean[1], pmean[2], point2[1], point2[2])
        }
        nearest_neighbour1 <- x[dist == min(dist),]
        
        nn <- rbind(nn, nearest_neighbour1)
      }
      
      tmean[i,1] <- mean(nn[,1])
      tmean[i,2] <- mean(nn[,2])
    }
    
    tmean <- unique(tmean)
    tmean <- rbind(start, tmean, end)
    dist <- pt2pt.distance(tmean$y, tmean$x) ## dist in m
    dist_log <- dist > 25000
    dist_log[1] <- FALSE
    
    if(iter_n == 100){      ###this will stop repeating after 100 iterations. option 2 is to stop repeating if distinces between tmean are all less than e.g. 25
      break
    }
    
    tmean1 <- data.frame()
    for(i in 1:length(dist_log)){
      if(dist_log[i]== T){
        newpt <- new_point(tmean[i-1,], tmean[i,], di = norm_vec(tmean[i-1,]- tmean[i,])/2) #new point at midpoint
        tmean1 <- rbind(tmean1, newpt)  ## order matters
      }
      tmean1 <- rbind(tmean1, tmean[i,])
    }
    tmean <- tmean1
    
  }
  
  tmean <-  SpatialLines(list(Lines(list(Line(cbind(tmean$x, tmean$y))), "id")))
  tmean <- spsample(tmean, equi_n, type = "regular") ## 100 point equally placed along length of line
  tmean <- tmean@coords
  

  traj_l <- lapply(traj_l, FUN = data.frame)
  for(i in 1:length(traj_l)){
    names(traj_l[[i]]) <- c("lon", "lat")
  }

  nn_dist <- NA
  for(l in seq_along(traj_l)){
    x <- traj_l[[l]]
    for(i in seq_along(tmean[,1])){
      point <- tmean[i,]
      dist <- NA
      for(j in seq_along(x[,1])){
        point2 <- x[j,]
        dist[j] <- deg.dist(point[1], point[2], point2[1,1], point2[1,2])
      }
      nn_dist[i] <- min(dist)
    }
    tmean <- cbind(tmean, nn_dist)
  }
  
  tmean <- data.frame(tmean)
  tmean$within.var <- rowSums(tmean[,3:ncol(tmean)])/(ncol(tmean)-2-1)  ##var = sum of (x-x.mn)/n-1.  here n = num traj = (ncol - 2(=tmean coords))
  tmean <- tmean[,c(1:2,ncol(tmean))]
  names(tmean) <- c("mn.lon", "mn.lat", "within.var")
  
  tmean$id <- r.pair.norep[[a, "id"]]
  tmean$birdyear <- r.pair.norep[[a, "birdyear"]]
  tmean$p.id <- r.pair.norep[[a, "p.id"]]
  tmean$p.birdyear <- r.pair.norep[[a, "p.birdyear"]]
  tmean$direction <- r.pair.norep[[a, "direction"]]
  rand_traj <- rbind(rand_traj,tmean)

}
saveRDS(rand_traj, "rand_traj_limlat.RDS")


#find variance
rand_traj_sum <- rand_traj %>% group_by(id, birdyear, p.id, p.birdyear, direction) %>% 
  summarise(mn_var = mean(within.var)) %>% ungroup()
rm(rand_traj)

rand_traj_sum$pair <- "Between" 

saveRDS(rand_traj_sum, "rand_traj_sum.RDS") ##Randomization test carried out in results_d4
saveRDS(r.pair.norep, "r_pair_traj.RDS")


##Route var examples ----
library(cowplot)
map.world <- rworldmap::getMap("low")
map.world <- fortify(map.world)
example.id <- c(5554,833,608, 782, 537, 4024,5337,540, 5060,5134,344)
mn.ex <- filter(mn_traj, id %in% example.id)

ID <- example.id[1]
##UDs
# UD <- nogap_multi.yr[[as.character(ID)]]
# p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear

pts <- filter(all_traj, id == ID) %>% mutate(birdyear = factor(birdyear))
mn <- filter(mn.ex, id == ID)

pv5554<- ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) +
  #geom_polygon(data = p50, aes(fill = id,col = id), size = 1, alpha = .1) +
  geom_path(data=pts, aes(lon, lat, group = birdyear), size = 1) +
  geom_point(data=pts, aes(lon, lat, group = NA), size = 1) +
  geom_point(data=mn, aes(mn.lon, mn.lat, group = NA, col = within.var), size = 1) +
  theme(legend.position = c(.85,.3), legend.background = element_rect(fill="white"),axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  scale_colour_viridis_c() +
  labs(col = "Variation \n (km)") +
  ggtitle("ID = 5554, Mean Variation = 13 km")+
  coord_fixed(xlim = c(-10, 10), 
              ylim = c(38.5, 53),ratio = 1.2) ## add core areas

ID <- example.id[10]
##UDs
# UD <- nogap_multi.yr[[as.character(ID)]]
# p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear

pts <- filter(all_traj, id == ID & direction == "Autumn") %>% mutate(birdyear = factor(birdyear))
mn <- filter(mn.ex, id == ID & direction == "Autumn")

pv5134<- ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) +
  # geom_polygon(data = p50, aes(fill = id,col = id), size = 1, alpha = .1) +
  geom_path(data=pts, aes(lon, lat, group = birdyear), size = 1) +
  geom_point(data=pts, aes(lon, lat, group = NA), size = 1) +
  geom_point(data=mn, aes(mn.lon, mn.lat, group = NA, col = within.var), size = 1) +
  theme(legend.position = c(.85,.3),legend.background = element_rect(fill="white"),axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  scale_colour_viridis_c() +
  ggtitle("ID = 5134, Mean Variation = 20 km") +
  labs(col = "Variation \n (km)") +
  coord_fixed(xlim = c(-7.5, 7.5), 
              ylim = c(40.5, 51.5),ratio = 1.2) ## add core areas



ID <- example.id[2]
##UDs
# UD <- nogap_multi.yr[[as.character(ID)]]
# p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear

pts <- filter(all_traj, id == ID & direction == "Autumn") %>% mutate(birdyear = factor(birdyear))
mn <- filter(mn.ex, id == ID & direction == "Autumn")

pv833<- ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) +
  # geom_polygon(data = p50, aes(fill = id,col = id), size = 1, alpha = .1) +
  geom_path(data=pts, aes(lon, lat, group = birdyear), size = 1) +
  geom_point(data=pts, aes(lon, lat, group = NA), size = 1) +
  geom_point(data=mn, aes(mn.lon, mn.lat, group = NA, col = within.var), size = 1) +
  theme(legend.position = c(.85,.3),legend.background = element_rect(fill="white"),axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  scale_colour_viridis_c() +
  ggtitle("ID = 833, Mean Variation = 17 km") +
  labs(col = "Variation \n (km)") +
  coord_fixed(xlim = c(-30, 22), 
              ylim = c(14, 52),ratio = 1.2) ## add core areas


ID <- example.id[3]
##UDs
# UD <- nogap_multi.yr[[as.character(ID)]]
# p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear

pts <- filter(all_traj, id == ID & direction == "Spring") %>% mutate(birdyear = factor(birdyear))
mn <- filter(mn.ex, id == ID & direction == "Spring")

pv608<- ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) +
  # geom_polygon(data = p50, aes(fill = id,col = id), size = 1, alpha = .1) +
  geom_path(data=pts, aes(lon, lat, group = birdyear), size = 1) +
  geom_point(data=pts, aes(lon, lat, group = NA), size = 1) +
  geom_point(data=mn, aes(mn.lon, mn.lat, group = NA, col = within.var), size = 1) +
  theme(legend.position = c(.85,.3),legend.background = element_rect(fill="white"),axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  scale_colour_viridis_c() +
  ggtitle("ID = 608, Mean Variation = 13 km") +
  labs(col = "Variation \n (km)") +
  coord_fixed(xlim = c(-2, 5), 
              ylim = c(49, 54),ratio = 1.2) ## add core areas



ID <- example.id[5]
##UDs
# UD <- nogap_multi.yr[[as.character(ID)]]
# p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear

pts <- filter(all_traj, id == ID & direction == "Autumn") %>% mutate(birdyear = factor(birdyear))
mn <- filter(mn.ex, id == ID & direction == "Autumn")

pv537<- ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) +
  geom_path(data=pts, aes(lon, lat, group = birdyear), size = 1) +
  geom_point(data=pts, aes(lon, lat, group = NA), size = 1) +
  geom_point(data=mn, aes(mn.lon, mn.lat, group = NA, col = within.var), size = 1) +
  theme(legend.position = c(.13,.3), legend.background = element_rect(fill="white"),axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  scale_colour_viridis_c() + 
  ggtitle("ID = 537, Mean Variation = 33 km") +
  labs(col = "Variation \n (km)") +
  coord_fixed(xlim = c(-2, 5), 
              ylim = c(49, 54),ratio = 1.2) ## add core areas 

ID <- example.id[7]
##UDs
# UD <- nogap_multi.yr[[as.character(ID)]]
# p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear

pts <- filter(all_traj, id == ID & direction == "Autumn") %>% mutate(birdyear = factor(birdyear))
mn <- filter(mn.ex, id == ID & direction == "Autumn")

pv5337<- ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) +
  # geom_polygon(data = p50, aes(fill = id,col = id), size = 1, alpha = .1) +
  geom_path(data=pts, aes(lon, lat, group = birdyear), size = 1) +
  geom_point(data=pts, aes(lon, lat, group = NA), size = 1) +
  geom_point(data=mn, aes(mn.lon, mn.lat, group = NA, col = within.var), size = 1) +
  theme(legend.position = c(.85,.3), legend.background = element_rect(fill="white"),axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  scale_colour_viridis_c() + scale_fill_viridis_d(begin = .1, end = .8) +
  ggtitle("ID = 5337, Mean Variation = 106 km") +
  labs(col = "Variation \n (km)") +
  coord_fixed(xlim = c(-10, 5), 
              ylim = c(40, 51),ratio = 1.2) ## add core areas

ID <- example.id[8]
##UDs
# UD <- nogap_multi.yr[[as.character(ID)]]
# p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear

pts <- filter(all_traj, id == ID & direction == "Spring") %>% mutate(birdyear = factor(birdyear))
mn <- filter(mn.ex, id == ID & direction == "Spring")

pv540<- ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) +
  geom_path(data=pts, aes(lon, lat, group = birdyear), size = 1) +
  geom_point(data=pts, aes(lon, lat, group = NA), size = 1) +
  geom_point(data=mn, aes(mn.lon, mn.lat, group = NA, col = within.var), size = 1) +
  theme(legend.position = c(.85,.3), legend.background = element_rect(fill="white"),axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  scale_colour_viridis_c() +
  ggtitle("ID = 540, Mean Variation = 111 km") +
  labs(col = "Variation \n (km)") +
  coord_fixed(xlim = c(-26, 22), 
              ylim = c(18, 53),ratio = 1.2) ## add core areas

#plot_grid(pv537, pv608, pv5337, pv5554, pv540, pv833, nrow = 3) 

png("tvar_ex.png",  width = 20, height = 30, units = "cm", res = 400)
plot_grid(pv537, pv608, pv5337, pv5134, pv540, pv833, nrow = 3) 
dev.off()



## Maps - overlap ----


library(cowplot)
theme_set(theme_bw())
example.id <- c(478,608,1400,483,5027,5296, 871,860,4047,5068)

## 478
ID <- example.id[1]
##UDs
UD <- nogap_multi.yr[[as.character(ID)]]
p25 <- fortify(spTransform(getverticeshr(UD, 25), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p75 <- fortify(spTransform(getverticeshr(UD, 75), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p95 <- fortify(spTransform(getverticeshr(UD, 95), CRS("+init=epsg:4326"))) ##@data$id is birdyear
##points
pts.d <- filter(nonbreed.nogap.d, id == ID) %>% select(lat, lon, time, birdyear) %>% mutate(birdyear = factor(birdyear))

OL <- gull.glm %>% filter(id == ID) %>% select(mn_overlap) %>% distinct() %>% pull(mn_overlap)

p478 <-
  ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) + ## map
  ##UDS per year
  geom_polygon(data = p95, aes(col = id, fill = NA), size = .2,  fill = NA) +
  geom_polygon(data = p75, aes(col = id, fill = NA), size = .5) +
  geom_polygon(data = p50, aes(fill = id, col = id), size = 1.2, alpha = .2) +
  geom_polygon(data = p25, aes(col = id, fill = NA), size = .5) +
  #points
  geom_point(data = pts.d, aes(lon, lat, group = birdyear, col = birdyear), size = .8)+
  ##appearance
  scale_colour_viridis_d(end = .9, option = "C") + scale_fill_viridis_d(begin = .1, end = .8) +

  theme(legend.position = "none", axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  coord_fixed(xlim = c(-4, 4.5), 
              ylim = c(49,55) ,ratio = 1.2) + 
  ggtitle(paste( " ID = ", ID, ", Overlap = ", round(OL, digits = 2), ", Seasons = ", 
                 gull.glm %>% filter(id == ID & !is.na(mig.dist)) %>% summarise(n = n()) %>% pull(n), sep = ""))

## 608
ID <- example.id[2]

UD <- nogap_multi.yr[[as.character(ID)]]
p25 <- fortify(spTransform(getverticeshr(UD, 25), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p75 <- fortify(spTransform(getverticeshr(UD, 75), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p95 <- fortify(spTransform(getverticeshr(UD, 95), CRS("+init=epsg:4326"))) ##@data$id is birdyear
pts.d <- filter(nonbreed.nogap.d, id == ID) %>% select(lat, lon, time, birdyear) %>% mutate(birdyear = factor(birdyear))
OL <- gull.glm %>% filter(id == ID) %>% select(mn_overlap) %>% distinct() %>% pull(mn_overlap)

p608 <-
  ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) + ## map
  geom_polygon(data = p95, aes(col = id), size = .2, fill = NA) +
  geom_polygon(data = p75, aes(fill = NA, col = id), size = .5) +
  geom_polygon(data = p50, aes(fill = id, col = id), size = 1.2, alpha = .2) +
  geom_polygon(data = p25, aes(fill = NA, col = id), size = .5) +
  geom_point(data = pts.d, aes(lon, lat, group = birdyear, col = birdyear), size = .8)+
  scale_colour_viridis_d(end = .9, option = "C") + scale_fill_viridis_d(begin = .1, end = .7) +
  theme(legend.position = "none", axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  coord_fixed(xlim = c(-5, 3.5), 
              ylim = c(49,55) ,ratio = 1.2) + 
  ggtitle(paste( " ID = ", ID, ", Overlap = ", round(OL, digits = 2), ", Seasons = ", 
                 gull.glm %>% filter(id == ID & !is.na(mig.dist)) %>% summarise(n = n()) %>% pull(n), sep = ""))

#1400
ID <- example.id[3]
UD <- nogap_multi.yr[[as.character(ID)]]
p25 <- fortify(spTransform(getverticeshr(UD, 25), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p75 <- fortify(spTransform(getverticeshr(UD, 75), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p95 <- fortify(spTransform(getverticeshr(UD, 95), CRS("+init=epsg:4326"))) ##@data$id is birdyear
pts.d <- filter(nonbreed.nogap.d, id == ID) %>% select(lat, lon, time, birdyear) %>% mutate(birdyear = factor(birdyear))
OL <- gull.glm %>% filter(id == ID) %>% select(mn_overlap) %>% distinct() %>% pull(mn_overlap)
p1400 <-
  ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) + ## map
  geom_polygon(data = p95, aes(col = id), size = .2, fill = NA) +
  geom_polygon(data = p75, aes(fill = NA, col = id), size = .5) +
  geom_polygon(data = p50, aes(fill = id, col = id), size = 1.2, alpha = .2) +
  geom_polygon(data = p25, aes(fill = NA, col = id), size = .5) +
  geom_point(data = pts.d, aes(lon, lat, group = birdyear, col = birdyear), size = .8)+
  scale_colour_viridis_d(end = .9, option = "C") + scale_fill_viridis_d(begin = .1, end = .7) +
  theme(legend.position = "none", axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  coord_fixed(xlim = c(-14, 7), 
              ylim = c(36,52) ,ratio = 1.2) + 
  ggtitle(paste( " ID = ", ID, ", Overlap = ", round(OL, digits = 2), ", Seasons = ", 
                 gull.glm %>% filter(id == ID & !is.na(mig.dist)) %>% summarise(n = n()) %>% pull(n), sep = ""))

#483
ID <- example.id[4]
UD <- nogap_multi.yr[[as.character(ID)]]
p25 <- fortify(spTransform(getverticeshr(UD, 25), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p75 <- fortify(spTransform(getverticeshr(UD, 75), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p95 <- fortify(spTransform(getverticeshr(UD, 95), CRS("+init=epsg:4326"))) ##@data$id is birdyear
pts.d <- filter(nonbreed.nogap.d, id == ID) %>% select(lat, lon, time, birdyear) %>% mutate(birdyear = factor(birdyear))
OL <- gull.glm %>% filter(id == ID) %>% select(mn_overlap) %>% distinct() %>% pull(mn_overlap)
p483 <-
  ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) + ## map
  geom_polygon(data = p95, aes(col = id), size = .2, fill = NA) +
  geom_polygon(data = p75, aes(fill = NA, col = id), size = .5) +
  geom_polygon(data = p50, aes(fill = id, col = id), size = 1.2, alpha = .2) +
  geom_polygon(data = p25, aes(fill = NA, col = id), size = .5) +
  geom_point(data = pts.d, aes(lon, lat, group = birdyear, col = birdyear), size = .8)+
  scale_colour_viridis_d(end = .9, option = "C") + scale_fill_viridis_d(begin = .1, end = .7) +
  theme(legend.position = "none", axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  coord_fixed(xlim = c(-15, 7), 
              ylim = c(37,54) ,ratio = 1.2) + 
  ggtitle(paste( " ID = ", ID, ", Overlap = ", round(OL, digits = 2), ", Seasons = ", 
                 gull.glm %>% filter(id == ID & !is.na(mig.dist)) %>% summarise(n = n()) %>% pull(n), sep = ""))

#5027
ID <- example.id[5]
UD <- nogap_multi.yr[[as.character(ID)]]
p25 <- fortify(spTransform(getverticeshr(UD, 25), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p75 <- fortify(spTransform(getverticeshr(UD, 75), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p95 <- fortify(spTransform(getverticeshr(UD, 95), CRS("+init=epsg:4326"))) ##@data$id is birdyear
pts.d <- filter(nonbreed.nogap.d, id == ID) %>% select(lat, lon, time, birdyear) %>% mutate(birdyear = factor(birdyear))
OL <- gull.glm %>% filter(id == ID) %>% select(mn_overlap) %>% distinct() %>% pull(mn_overlap)
p5027 <-
  ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) + ## map
  geom_polygon(data = p95, aes(col = id), size = .2, fill = NA) +
  geom_polygon(data = p75, aes(fill = NA, col = id), size = .5) +
  geom_polygon(data = p50, aes(fill = id, col = id), size = 1.2, alpha = .2) +
  geom_polygon(data = p25, aes(fill = NA, col = id), size = .5) +
  geom_point(data = pts.d, aes(lon, lat, group = birdyear, col = birdyear), size = .8)+
  scale_colour_viridis_d(end = .9, option = "C") + scale_fill_viridis_d(begin = .1, end = .7) +
  theme(legend.position = "none", axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  coord_fixed(xlim = c(-20, 11), 
              ylim = c(32,55.5) ,ratio = 1.2) + 
  ggtitle(paste( " ID = ", ID, ", Overlap = ", round(OL, digits = 2), ", Seasons = ", 
                 gull.glm %>% filter(id == ID & !is.na(mig.dist)) %>% summarise(n = n()) %>% pull(n), sep = ""))

## 5296
ID <- example.id[6]
UD <- nogap_multi.yr[[as.character(ID)]]
p25 <- fortify(spTransform(getverticeshr(UD, 25), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p75 <- fortify(spTransform(getverticeshr(UD, 75), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p95 <- fortify(spTransform(getverticeshr(UD, 95), CRS("+init=epsg:4326"))) ##@data$id is birdyear
pts.d <- filter(nonbreed.nogap.d, id == ID) %>% select(lat, lon, time, birdyear) %>% mutate(birdyear = factor(birdyear))
OL <- gull.glm %>% filter(id == ID) %>% select(mn_overlap) %>% distinct() %>% pull(mn_overlap)
p5296 <-
  ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) + ## map
  geom_polygon(data = p95, aes(col = id), size = .2, fill = NA) +
  geom_polygon(data = p75, aes(fill = NA, col = id), size = .5) +
  geom_polygon(data = p50, aes(fill = id, col = id), size = 1.2, alpha = .2) +
  geom_polygon(data = p25, aes(fill = NA, col = id), size = .5) +
  geom_point(data = pts.d, aes(lon, lat, group = birdyear, col = birdyear), size = .8)+
  scale_colour_viridis_d(end = .9, option = "C") + scale_fill_viridis_d(begin = .1, end = .7) +
  theme(legend.position = "none", axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  coord_fixed(xlim = c(-32, 17), 
              ylim = c(17,55) ,ratio = 1.2) + 
  ggtitle(paste( " ID = ", ID, ", Overlap = ", round(OL, digits = 2), ", Seasons = ", 
                 gull.glm %>% filter(id == ID & !is.na(mig.dist)) %>% summarise(n = n()) %>% pull(n), sep = ""))


## 871
ID <- example.id[7]
UD <- nogap_multi.yr[[as.character(ID)]]
p25 <- fortify(spTransform(getverticeshr(UD, 25), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p75 <- fortify(spTransform(getverticeshr(UD, 75), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p95 <- fortify(spTransform(getverticeshr(UD, 95), CRS("+init=epsg:4326"))) ##@data$id is birdyear
pts.d <- filter(nonbreed.nogap.d, id == ID) %>% select(lat, lon, time, birdyear) %>% mutate(birdyear = factor(birdyear))
OL <- gull.glm %>% filter(id == ID) %>% select(mn_overlap) %>% distinct() %>% pull(mn_overlap)
p871 <-
  ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) + ## map
  geom_polygon(data = p95, aes(col = id), size = .2, fill = NA) +
  geom_polygon(data = p75, aes(fill = NA, col = id), size = .5) +
  geom_polygon(data = p50, aes(fill = id, col = id), size = 1.2, alpha = .2) +
  geom_polygon(data = p25, aes(fill = NA, col = id), size = .5) +
  geom_point(data = pts.d, aes(lon, lat, group = birdyear, col = birdyear), size = .8)+
    scale_colour_viridis_d(end = .9, option = "C") + scale_fill_viridis_d(begin = .1, end = .7) +
    theme(legend.position = "none", axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
          panel.grid.major = element_line(colour= "white")) +
    coord_fixed(xlim = c(-5, 3.5), 
                ylim = c(49,55) ,ratio = 1.2) + 
  ggtitle(paste( " ID = ", ID, ", Overlap = ", round(OL, digits = 2), ", Seasons = ", 
                 gull.glm %>% filter(id == ID & !is.na(mig.dist)) %>% summarise(n = n()) %>% pull(n), sep = ""))

#5068
ID <- example.id[10]
UD <- nogap_multi.yr[[as.character(ID)]]
p25 <- fortify(spTransform(getverticeshr(UD, 25), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p75 <- fortify(spTransform(getverticeshr(UD, 75), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p95 <- fortify(spTransform(getverticeshr(UD, 95), CRS("+init=epsg:4326"))) ##@data$id is birdyear
pts.d <- filter(nonbreed.nogap.d, id == ID) %>% select(lat, lon, time, birdyear) %>% mutate(birdyear = factor(birdyear))
OL <- gull.glm %>% filter(id == ID) %>% select(mn_overlap) %>% distinct() %>% pull(mn_overlap)
p5068 <-
  ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) + ## map
  geom_polygon(data = p95, aes(col = id), size = .2, fill = NA) +
  geom_polygon(data = p75, aes(fill = NA, col = id), size = .5) +
  geom_polygon(data = p50, aes(fill = id, col = id), size = 1.2, alpha = .2) +
  geom_polygon(data = p25, aes(fill = NA, col = id), size = .5) +
  geom_point(data = pts.d, aes(lon, lat, group = birdyear, col = birdyear), size = .8)+
  scale_colour_viridis_d(end = .9, option = "C") + scale_fill_viridis_d(begin = .1, end = .7) +
  theme(legend.position = "none", axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  coord_fixed(xlim = c(-20, 11), 
              ylim = c(30,53.5) ,ratio = 1.2) + 
  ggtitle(paste( " ID = ", ID, ", Overlap = ", round(OL, digits = 2), ", Seasons = ", 
                 gull.glm %>% filter(id == ID & !is.na(mig.dist)) %>% summarise(n = n()) %>% pull(n), sep = ""))


#860
ID <- example.id[8]
UD <- nogap_multi.yr[[as.character(ID)]]
p25 <- fortify(spTransform(getverticeshr(UD, 25), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p50 <- fortify(spTransform(getverticeshr(UD, 50), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p75 <- fortify(spTransform(getverticeshr(UD, 75), CRS("+init=epsg:4326"))) ##@data$id is birdyear
p95 <- fortify(spTransform(getverticeshr(UD, 95), CRS("+init=epsg:4326"))) ##@data$id is birdyear
pts.d <- filter(nonbreed.nogap.d, id == ID) %>% select(lat, lon, time, birdyear) %>% mutate(birdyear = factor(birdyear))
OL <- gull.glm %>% filter(id == ID) %>% select(mn_overlap) %>% distinct() %>% pull(mn_overlap)
p860 <-
  ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) + ## map
  geom_polygon(data = p95, aes(col = id), size = .2, fill = NA) +
  geom_polygon(data = p75, aes(fill = NA, col = id), size = .5) +
  geom_polygon(data = p50, aes(fill = id, col = id), size = 1.2, alpha = .2) +
  geom_polygon(data = p25, aes(fill = NA, col = id), size = .5) +
  geom_point(data = pts.d, aes(lon, lat, group = birdyear, col = birdyear), size = .8)+
  scale_colour_viridis_d(end = .9, option = "C") + scale_fill_viridis_d(begin = .1, end = .7) +
  theme(legend.position = "none", axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  coord_fixed(xlim = c(-15, 7), 
              ylim = c(37,54) ,ratio = 1.2) + 
  ggtitle(paste( " ID = ", ID, ", Overlap = ", round(OL, digits = 2), ", Seasons = ", 
                 gull.glm %>% filter(id == ID & !is.na(mig.dist)) %>% summarise(n = n()) %>% pull(n), sep = ""))



ol.page <- plot_grid(p478,p871,p1400,p860,p5027,p5068, nrow = 3)

png("ol_ex.png",  width = 20, height = 30, units = "cm", res = 400)
ol.page
dev.off()

## Mean traj computation example ----
all_traj <- ungroup(all_traj) %>% arrange(id, time)
traj5524 <- filter(all_traj, id == 5524 & direction == "Autumn")
p5524 <- nogap_multi.yr.50p[[74]]
p5524 <-  spTransform(p5524, CRS("+init=epsg:4326"))
p5524 <- fortify(p5524)
p5524.0 <- filter(p5524, id %in% c("1","2"))
p5524.1 <- filter(p5524, id %in% c("3","4"))

### first need plot with UD + real points
t.0 <- filter(traj5524, birdyear == 0) %>% select(lon, lat)
t.1<- filter(traj5524, birdyear == 1) %>% select(lon, lat)

SFa <- ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) +
  geom_polygon(data = p5524.0, fill = "#0D0887FF",col = "#0D0887FF", size = 1, alpha = .1) +
  geom_polygon(data = p5524.1, fill = "#CC4678FF",col = "#CC4678FF", size = 1, alpha = .1) +
  geom_path(data=t.0, aes(lon, lat, group = NA), col = "#0D0887FF", size = 1) +
  geom_path(data = t.1, aes(lon, lat,group = NA),col = "#CC4678FF", size = 1) +
  geom_point(data = t.0, aes(lon, lat,group = NA), col = "#0D0887FF", size = 2) +
  geom_point(data = t.1, aes(lon, lat,group = NA),col = "#CC4678FF", size = 2) +
  annotate(xmin = -11, xmax = -1, ymin = 37, ymax = 44, 
           geom = "rect", alpha = 0, col = "black") +
  theme(legend.position = "none", axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
        panel.grid.major = element_line(colour= "white")) +
  coord_fixed(xlim = c(-15, 14), 
              ylim = c(32.5, 53),ratio = 1.2) ## add core areas

##df of single direction trajectories, multiples identified by birdyear
## df includes  id, birdyear, lon, lat, direction, date_time
## points within wintering area polygons (i.e. stopovers) were replaced with the polygon centroid to aid in  
#spacing points equally along the trajectory 
id.dir <- traj5524 %>%  select(id, direction) %>% distinct() 
mn_traj <- data.frame(matrix(ncol = 5, nrow = 0))
names(mn_traj) <- c("mn.lon", "mn.lat", "within.var", "id", "direction")
a <- 1

traj <- filter(traj5524, id == id.dir[a,]$id & direction == id.dir[a,]$direction)

###This creates a list of equally spaced points (currently n=11, specified in t1_eq), with each trajectory having it's own list element. 
ntraj <- unique(traj$birdyear)  #trajectory id
equi_n <- 100  ### number of points along the trajectory
traj_l <- list()
for(i in 1:length(ntraj)){
  t1 <- filter(traj, birdyear == ntraj[i]) %>% select(lon, lat) %>% distinct()
  if (length(t1$lon)<=10){ ntraj <- ntraj[-i] 
  } else {
    t1_l <-  SpatialLines(list(Lines(list(Line(cbind(t1$lon, t1$lat))), "id")))
    t1_eq <- spsample(t1_l, equi_n, type = "regular") ## 250 point equally placed along length of line
    traj_l[[i]] <- t1_eq@coords
    
  }
}

names(traj_l) <- ntraj

traj_l_sample <- lapply(traj_l, function(x) {
  x <- as.data.frame(x)
  x$seq <- 1:equi_n
  filter(x, coords.x1 >= -12 & coords.x1 <=0 & coords.x2 >= 36 & coords.x2 <= 45)
})    

t0.0 <- data.frame(traj_l_sample[[1]])  
t0.1 <- data.frame(traj_l_sample[[2]]) 

###Create starting 'thread' for mean trajectory between start and end midpoints

start <- cbind(x = mean(sapply(traj_l, '[[',1,1)), y = mean(sapply(traj_l, '[[',1,2))) ##mean of first elements in each list
end <- cbind(x = mean(sapply(traj_l, '[[',equi_n,1)), y = mean(sapply(traj_l, '[[',equi_n,2)))

###Create starting 'thread' from mean trajectory points
tmean <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(tmean) <- c("x","y")

for(i in 1:equi_n){
  tmean[i,1] = mean(sapply(traj_l, '[[',i,1))
  tmean[i,2] = mean(sapply(traj_l, '[[',i,2))
}

t1.0 <- tmean %>% mutate(seq = 1:equi_n) %>% filter(x > -12 & x < 0 & y > 36 & y < 45)

SFb <- ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) +
  geom_path(data=t0.0, aes(coords.x1, coords.x2, group = NA), size = 1, alpha = .5, linetype = '21') +
  geom_path(data = t0.1, aes(coords.x1, coords.x2,group = NA), size = 1, alpha = .5, linetype = '21') +
geom_point(data = t0.0, aes(coords.x1, coords.x2,group = NA, fill = seq), size = 2, shape = 21) +
  geom_point(data = t0.1, aes(coords.x1, coords.x2,group = NA, fill = seq), size = 2, shape = 21) +
 geom_path(data = t1.0, aes(x, y,group = NA), size = 1) +
  geom_point(data = t1.0, aes(x, y,group = NA, fill = seq), size = 2, shape = 21) +
  geom_point(data = t0.0[t0.0$seq %in% c(50,59,64,81),], aes(coords.x1, coords.x2,group = NA, fill = seq), size = 4, shape = 21) +
  geom_point(data = t0.1[t0.1$seq %in% c(50,59,64,81),], aes(coords.x1, coords.x2,group = NA, fill = seq), size = 4, shape = 21) +
geom_point(data = t1.0[t1.0$seq %in% c(50,59,64,81),], aes(x, y,group = NA, fill = seq), size = 4, shape = 21,) +
  scale_fill_viridis_c(option = "C") + theme(legend.position = "none", axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
                                             panel.grid.major = element_line(colour= "white")) +
  coord_fixed(xlim = c(-11, -1), 
              ylim = c(37, 44),ratio = 1.2)  

### nearest neighbour points for 50,59,64,81:
nn.ex <- data.frame()
for(i in c(50,59,64,81)){
  pmean <- as.numeric(tmean[i,]) ## select each  point
  
  nn <- data.frame()
  for(j in seq_along(traj_l)){
    x <- traj_l[[j]]
    x<- data.frame(x)
    x$seq <- 1:equi_n
    x <- x[x[,2] < pmean[2]+1 & x[,2]  > pmean[2]-1,]
    if (length(x) == 0) { x <- traj_l[[j]]}
    dist <- NA
    for(k in seq_along(x[,1])){
      point2 <- as.numeric(x[k,])
      dist[k] <- deg.dist(pmean[1], pmean[2], point2[1], point2[2])
    }
    nearest_neighbour1 <- x[dist == min(dist),]
    
    nn <- rbind(nn,nearest_neighbour1)
    
  }
  
  names(nn) <- c("x", "y", "seq")
  x <- mean(nn[,1])
  y <- mean(nn[,2])
  seq <- i
  nn <- rbind(nn, data.frame(x,y, seq))
  nn.ex <- rbind(nn, nn.ex)
  
}   

##Averaging

iter_n <- 0

iter_n <- iter_n + 1

for(i in seq_along(tmean$x)){
  pmean <- as.numeric(tmean[i,]) ## select each  point
  
  nn <- data.frame()
  for(j in seq_along(traj_l)){
    x <- traj_l[[j]]
    x <- matrix(x[x[,2] < pmean[2]+1 & x[,2]  > pmean[2]-1,], ncol = 2)
    if (length(x) == 0) { x <- traj_l[[j]]}
    dist <- NA
    for(k in seq_along(x[,1])){
      point2 <- as.numeric(x[k,])
      dist[k] <- deg.dist(pmean[1], pmean[2], point2[1], point2[2])
    }
    nearest_neighbour1 <- x[dist == min(dist),]
    
    nn <- rbind(nn,nearest_neighbour1)
  }
  
  tmean[i,1] <- mean(nn[,1])
  tmean[i,2] <- mean(nn[,2])
}   

tmean <- unique(tmean)
tmean <- rbind(start, tmean, end)

t1.1 <- tmean %>% mutate(seq = 1:length(tmean$x)) %>%  filter(x > -12 & x < 0 & y > 36 & y < 45)

dist <- pt2pt.distance(tmean$y, tmean$x) ## dist in m
dist_log <- dist > 25
dist_log[1] <- FALSE

tmean1 <- data.frame()
for(i in 1:length(dist_log)){
  if(dist_log[i]== T){
    newpt <- new_point(tmean[i-1,], tmean[i,], di = norm_vec(tmean[i-1,]- tmean[i,])/2) #new point at midpoint 
    tmean1 <- rbind(tmean1, newpt)  ## order matters
  }
  tmean1 <- rbind(tmean1, tmean[i,])
}

new.pts <- anti_join(tmean1, tmean)
tmean <- tmean1

SFc <- ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) +
  geom_path(data=t0.0, aes(coords.x1, coords.x2, group = NA), size = 1, alpha = .5, linetype = '21') +
  geom_path(data = t0.1, aes(coords.x1, coords.x2,group = NA), size = 1, alpha = .5, linetype = '21') +
 geom_point(data = t0.0, aes(coords.x1, coords.x2,group = NA, fill = seq), size = 2, shape = 21) +
  geom_point(data = t0.1, aes(coords.x1, coords.x2,group = NA, fill = seq), size = 2, shape = 21) +
 geom_path(data = t1.0, aes(x, y,group = NA), size = 1, alpha = .25) +
  geom_point(data = t1.0, aes(x, y,group = NA, fill = seq), size = 2, shape = 21, alpha = .25) +
  geom_path(data = t1.1, aes(x, y,group = NA), size = 1) +
  geom_point(data = new.pts, aes(x, y,group = NA), fill = '#20A387FF', size = 2, shape = 21) +
  geom_point(data = t1.1, aes(x, y,group = NA, fill = seq), size = 2, shape = 21) +
  geom_point(data = nn.ex,aes(x, y,group = NA, fill = seq), size = 4, shape = 21) +
  
  scale_fill_viridis_c(option = "C") + theme(legend.position = "none", axis.title = element_blank(),
                                             panel.background = element_rect(fill = "gray93"),
                                             panel.grid.major = element_line(colour= "white")) + 
  coord_fixed(xlim = c(-11, -1), 
              ylim = c(37, 44),ratio = 1.2)  

###
###Create starting 'thread' for mean trajectory between start and end midpoints

start <- cbind(x = mean(sapply(traj_l, '[[',1,1)), y = mean(sapply(traj_l, '[[',1,2))) ##mean of first elements in each list
end <- cbind(x = mean(sapply(traj_l, '[[',equi_n,1)), y = mean(sapply(traj_l, '[[',equi_n,2)))
#tmean <- data.frame(x = seq(start[,1], end[,1], length.out = 250), y = seq(start[,2], end[,2], length.out = 250))

###
###Create starting 'thread' from mean trajectory points
tmean <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(tmean) <- c("x","y")

for(i in 1:equi_n){
  tmean[i,1] = mean(sapply(traj_l, '[[',i,1))
  tmean[i,2] = mean(sapply(traj_l, '[[',i,2))
}


###
##Averaging

iter_n <- 0
repeat{
  iter_n <- iter_n + 1
  
  for(i in seq_along(tmean$x)){
    pmean <- as.numeric(tmean[i,]) ## select each  point
    
    nn <- data.frame()
    for(j in seq_along(traj_l)){
      x <- traj_l[[j]]
      x <- matrix(x[x[,2] < pmean[2]+1 & x[,2]  > pmean[2]-1,], ncol = 2)
      if (length(x) == 0) { x <- traj_l[[j]]}
      dist <- NA
      for(k in seq_along(x[,1])){
        point2 <- as.numeric(x[k,])
        dist[k] <- deg.dist(pmean[1], pmean[2], point2[1], point2[2])
      }
      nearest_neighbour1 <- x[dist == min(dist),]
      
      nn <- rbind(nn,nearest_neighbour1)
    }
    
    tmean[i,1] <- mean(nn[,1])
    tmean[i,2] <- mean(nn[,2])
  }
  
  tmean <- unique(tmean)
  tmean <- rbind(start, tmean, end)
  dist <- pt2pt.distance(tmean$y, tmean$x) ## dist in m
  dist_log <- dist > 25
  dist_log[1] <- FALSE
  
  if(iter_n == 100){      ###this will stop repeating after 100 iterations. option 2 is to stop repeating if distinces between tmean are all less than e.g. 25
    break
  }
  
  tmean1 <- data.frame()
  for(i in 1:length(dist_log)){
    if(dist_log[i]== T){
      newpt <- new_point(tmean[i-1,], tmean[i,], di = norm_vec(tmean[i-1,]- tmean[i,])/2) #new point at midpoint 
      tmean1 <- rbind(tmean1, newpt)  ## order matters
    }
    tmean1 <- rbind(tmean1, tmean[i,])
  }
  tmean <- tmean1

}

tmean <-  SpatialLines(list(Lines(list(Line(cbind(tmean$x, tmean$y))), "id")))
tmean <- spsample(tmean, equi_n, type = "regular") ## 100 point equally placed along length of line
tmean <- tmean@coords

traj_l <- lapply(traj_l, FUN = data.frame)
for(i in 1:length(traj_l)){
  names(traj_l[[i]]) <- c("lon", "lat")
}

ex2 <- c(52,59,65,79)

nn.ex2 <- data.frame()
for(i in ex2){
  pmean <- as.numeric(tmean[i,]) ## select each  point
  
  nn <- data.frame()
  for(j in seq_along(traj_l)){
    x <- traj_l[[j]]
    x <-x[x[,2] < pmean[2]+1 & x[,2]  > pmean[2]-1,]
    if (length(x$lon) == 0) { x <- traj_l[[j]]}
    dist <- NA
    for(k in seq_along(x[,1])){
      point2 <- as.numeric(x[k,])
      dist[k] <- deg.dist(pmean[1], pmean[2], point2[1], point2[2])
    }
    nearest_neighbour1 <- x[dist == min(dist),]
    
    nn <- rbind(nn,nearest_neighbour1)
    
  }
  names(nn) <- c("x", "y")
  nn$mn_seq <- i
  nn.ex2 <- rbind(nn, nn.ex2)
  
} 

nn_dist <- NA
for(l in seq_along(traj_l)){
  x <- traj_l[[l]]
  for(i in seq_along(tmean[,1])){
    point <- tmean[i,]
    dist <- NA
    for(j in seq_along(x[,1])){
      point2 <- x[j,]
      dist[j] <- deg.dist(point[1], point[2], point2[1,1], point2[1,2])
    }
    nn_dist[i] <- min(dist)
  }
  tmean <- cbind(tmean, nn_dist)
}

tmean <- data.frame(tmean)
tmean$within.var <- rowSums(tmean[,3:ncol(tmean)])/(ncol(tmean)-2-1)  ##var = sum of (x-x.mn)/n-1.  here n = num traj = (ncol - 2(=tmean coords))
tmean <- tmean[,c(1:2,ncol(tmean))]
names(tmean) <- c("mn.lon", "mn.lat", "within.var")
nn.ex2$wvar <- rep(tmean[rev(ex2),]$within.var, each = length(traj_l))

SFd<-  ggplot(map.world, aes(long, lat, group = group)) + 
  geom_polygon(fill = "white", col = "black", size = .25) +
  geom_path(data=t0.0, aes(coords.x1, coords.x2, group = NA), size = 1, alpha = .5, linetype = '21') +
  geom_path(data = t0.1, aes(coords.x1, coords.x2,group = NA), size = 1, alpha = .5, linetype = '21') +
  geom_point(data = t0.0, aes(coords.x1, coords.x2,group = NA), size = 2, shape = 21, fill = "black") +
  geom_point(data = t0.1, aes(coords.x1, coords.x2,group = NA), size = 2, shape = 21, fill = "black") +

  geom_point(data = nn.ex2, aes(x,y,group = NA, fill = wvar), size = 4, shape = 21) +
  geom_path(data = tmean, aes(mn.lon, mn.lat, group = NA), size = 1) +
  geom_point(data = tmean, aes(mn.lon, mn.lat, group = NA, fill = within.var), size = 2, shape = 21) +
  geom_point(data = tmean[ex2,], aes(mn.lon, mn.lat, group = NA, fill = within.var), size = 4, shape = 21) +
  scale_fill_viridis_c(name = "Variance \n (km)") + theme(axis.title = element_blank(), panel.background = element_rect(fill = "gray93"),
                                                          panel.grid.major = element_line(colour= "white"), legend.position = c(.68,.33)) +
  coord_fixed(xlim = c(-11, -1), 
              ylim = c(37, 44),ratio = 1.2)   

png("FigS3.png",  width = 20, height = 16, units = "cm", res = 400)
cowplot::plot_grid(SFa,SFb,SFc,SFd, nrow = 2)
dev.off()
