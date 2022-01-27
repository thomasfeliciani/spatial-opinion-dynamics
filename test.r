#


rm(list = ls())

# ______________________________________________________________________________
# Looking for ways of creating arbitrary spatial distributions.
#
# fake district:
library(ggplot2)
library(ggpubr)

vk <- expand.grid(
  east = 1:50,
  north = 1:50
)
# and removing some of its cells:
vk <- vk[sample(1:nrow(vk), size = 0.8 * nrow(vk), replace = FALSE),]


vk$x <- rbeta(n = nrow(vk), shape1 = 3, shape2 = 3)

vk$xx <- vk$x ^ 2



makeGradient <- function (hotspot, coords, steepness = 0) {
  x <- c()
  
  # Calculating the Euclidean distance of each square from the hotspot:
  for (c in 1:nrow(coords)) {
    distance = sqrt(sum((hotspot - c(coords$east[c], coords$north[c])) ^ 2))
    x[c] <- distance
  }
  
  # Transofrming the distance values so that (1) they are inverted (squares
  # closer to the hotspot get higher values); (2) they range between 0 and 1;
  # (3) we introduce a "steepness" parameter that determines how close to the
  # hotspot the values drop to zero.
  x <- 1 - ((x - min(x)) / ((max(x) - min(x)) * (1 - steepness)))
  #x <- 1 - ((x - min(x)) / (max(x) - min(x)) * steepness)#(1 - steepness))
  
  # Ensuring range is correct:
  x[x < 0] <- 0
  
  return(x)
}

vk$xx <- makeGradient(
  hotspot = c(east = 10, north = 10),
  coords = vk,
  steepness = 0.3
)


logAmplify <- function(x, steepness = 10) {
  L = 1 # max value
  #steepness = 10 # steepness
  x0 = 0.5 # sigmoid's midpoint
  e = exp(1) # Euler's constant e
  return(
    L / (1 + (e ^ (- steepness * (x - x0))))
  )
}
#plot(logAmplify(0:100/100, steepness = 10))
#cor(vk$x, vk$xx, method = "spearman")
vk$xx <- logAmplify(vk$x, steepness = 2)


xplot <- ggplot(data = vk, aes(x = east, y = north, color = x)) +
  geom_point(shape = 15, size = 2) +
  scale_color_viridis_c(option = "C", limits = c(0,1))
xxplot <- ggplot(data = vk, aes(x = east, y = north, color = xx)) +
  geom_point(shape = 15, size = 2) +
  scale_color_viridis_c(option = "C", limits = c(0,1))
xhist <- ggplot(data = vk, aes(x = x)) + geom_histogram() + xlim(c(0,1))
xxhist <- ggplot(data = vk, aes(x = xx)) + geom_histogram() + xlim(c(0,1))

print(ggpubr::ggarrange(
  xplot, xxplot, xhist, xxhist, ncol = 2, nrow = 2,
  hjust = -1, common.legend = FALSE
))




plot(logAmplify(0:100/100, steepness = 10), ylim = c(0,1))


#








logAmplify <- function(x) {
  L = 1 # max value
  k = 10 # steepness
  x0 = 0.5 # sigmoid's midpoint
  e = exp(1) # Euler's constant e
  return(
    L / (1 + (e ^ (-k * (x - x0))))
  )
}

x <- 1:100 / 100
xx <- sapply(test, FUN = logamplify)
plot(x, xx)





















library("geosphere")










# Calculating the matrix of distances between all pairs of squares. We
# save them in a matrix "distmat":
coords <- as.matrix(cbind(x$east, x$north))
distmat <- matrix(NA, nrow = nrow(x), ncol = nrow(x))
for (i in 1:nrow(x)) {
  #distmat[i,] <- geosphere::distVincentyEllipsoid( # more precise
  distmat[i,] <- geosphere::distMeeus( # faster
    p1 = coords[i,],
    p2 = coords
  )
}

# We calculate the resulting proximity matrix.
proxmat <- matrix(NA, nrow = nrow(x), ncol = nrow(x))
proxmat <- exp(-distmat / input$distanceDecaySlope)

# We also need to decide how we treat a square's distance from itself.
# Here we set it to about 52 (meters), i.e. the average distance between
# random points in a 100 by 100 meters square.
diag(proxmat) <- exp(-52.140543316 / input$distanceDecaySlope)

# We can finally calculate the LISA:
moranI <- moranI(x = x$pnw, proxmat = proxmat, dens = x$density)




################################################################################

# experimenting with segregating by line

x$test <- apply(
  as.matrix(cbind(x$east, x$north)),
  MARGIN = 1,
  FUN = function (p) {
    geosphere::dist2Line(
      p = p,
      line = as.matrix(rbind(
        c(min(x$east), min(x$north)),
        c(max(x$east), max(x$north))
      )),
      distfun = distMeeus
    )
  }
) [1,] # We only need the first raw, which is the one containing the distances.

ggplot(x, aes(x = east, y = north, color = test)) + geom_point()


####################3

bearing = 0

linep <- geosphere::destPoint(
  p = c(median(x$east), median(x$north)), # approximate center
  b = c(bearing -5, bearing + 5), # bearing; tangent to the target line
  d = 30000 # Distance (meters) from center to outskirts. Just a big number.
)

x$distToLine <- apply(
  as.matrix(cbind(x$east, x$north)),
  MARGIN = 1,
  FUN = function (p) {
    geosphere::dist2Line(
      p = p,
      line = linep,
      distfun = distMeeus
    )
  }
) [1,]

x$distToLine <- normalize(x$distToLine)

ggplot(x, aes(x = east, y = north, color = distToLine)) +
  geom_point() +
  geom_point(data = as.data.frame(linep), aes(x = lon, y = lat), color = "red")












################################################################################
################################################################################
################################################################################
################################################################################
#
# Exploring noise
rm(list = ls())
load("./geoData/geoData.RData")
#vks <- rgdal::readOGR(
#  dsn = "./geoData/2021-cbs_vk100_2020_v1",
#  layer = "CBS_vk100_2020_v1"
#)



wijk = 1
sample = 1#0.01


vkx <- vk[vk$district == districts$name[wijk],]
npop <- sum(vkx$density)

# Creating blank agentset...
w <- data.frame(
  sq = rep(NA, times = npop),
  g = rep(0, times = npop),
  o = rep(NA, times = npop),
  east = rep(NA, times = npop),
  north = rep(NA, times = npop),
  rdCoords = rep(NA, times = npop)
)

# ... and adding attributes (coordinates, group) according to the density and
# group composition of each square.
counter <- 1
for (sq in 1:nrow(vkx)) {
  #print(sq)
  residents <- counter:(counter + (vkx$density[sq] - 1))
  w$sq[residents] <- vkx$id[sq]
  w$east[residents] <- vkx$east[sq]
  w$north[residents] <- vkx$north[sq]
  w$rdCoords[residents] <- vkx$rdCoords[sq]
  
  # Reassigning the group w$g of non-west-background residents, if there are any
  if (vkx$pnw[sq] != 0) {
    numberNonWest <- round(vkx$density[sq] * vkx$pnw[sq])
    w$g[counter:(counter + numberNonWest - 1)] <- 1
  }
  
  counter <- counter + vkx$density[sq]
}

# Sampling agentset
w <- w[runif(n = round(nrow(w) * sample), min = 1, max = nrow(w)),]




noiseCoord <- function (w, rast, dist, proportion = 0.5, s = 100) {
  sample <- sample(1:nrow(w), size = round(nrow(w) * proportion))
  
  for (i in sample) { # For each agent i that we are going to move
    
    # In which row of the raster dataframe is the agent?
    startCellKey <- which(rast$rdCoords == w$rdCoords[i])
    
    # Sampling from its neighboring cells where it should go.
    # Note that the probability that a cell is chosen depends on its proximity
    # to the starting cell:
    targetCell <- sample(
      x = rast$rdCoords,
      size = 1,
      prob = exp(-dist[startCellKey,] / s)
    )
    targetCellKey <- which(rast$rdCoords == targetCell)
    
    # Now we can update location and coordinates of the agent i:
    w$rdCoords[i] <- targetCell
    w$east[i] <- rast$east[targetCellKey]
    w$north[i] <- rast$north[targetCellKey]
  }
  
  return(w)
}

noiseLocation <- function (w, rast, dist, proportion = 0.5, s = 100) {
  
  # We choose a cell in the center of the district:
  distanceToMedianCoord <- geosphere::distMeeus(
    p1 = c(median(w$east), median(w$north)),
    p2 = as.matrix(cbind(rast$east, rast$north))
  )
  startCellKey <- which(distanceToMedianCoord == min(distanceToMedianCoord))
  
  # And calculate its proximity to all cells:
  proximity <- exp(-dist[startCellKey,] / s)
  
  # We now move (a sample of) agents to random locations around the center cell
  sample <- sample(1:nrow(w), size = round(nrow(w) * proportion))
  w$rdCoords[sample] <- sample(
    x = rast$rdCoords,
    size = length(sample),
    prob = proximity,
    replace = TRUE
  )
  
  # ... and we update their coordinates accordingly.
  for (i in sample) {
    targetCellKey <- which(rast$rdCoords == w$rdCoords[i])
    w$east[i] <- rast$east[targetCellKey]
    w$north[i] <- rast$north[targetCellKey]
  }
  
  return(w)
}

noiseSwap <- function (w, proportion = 0.5) {
  gt <- w$g
  ot <- w$o
  sample <- sample(1:nrow(w), size = round(nrow(w) * proportion))
  sampleA <- sample[1:floor(length(sample) / 2)]
  sampleB <- sample[(floor(length(sample) / 2) + 1):(length(sampleA) * 2)]
  for (i in 1:length(sampleA)) {
    w$g[sampleA[i]] <- gt[sampleB[i]]
    w$o[sampleA[i]] <- ot[sampleB[i]]
    w$g[sampleB[i]] <- gt[sampleA[i]]
    w$o[sampleB[i]] <- ot[sampleA[i]]
  }
  return(w)
  #return(data.frame(g = g, o = o))
}



ww <- noiseCoord(
  w,
  rast = raster[[wijk]],
  dist = distances[[wijk]],
  proportion = 0.5,
  s = 100
)
ww <- noiseLocation(
  w,
  rast = raster[[wijk]],
  dist = distances[[wijk]],
  proportion = 0.5,
  s = 200
)
ww <- noiseSwap(
  w, 
  proportion = 0.5
)
plot(w$east, w$north)
plot(ww$east, ww$north)


# Restructuring the raster: this entails calculating, for each cell, its density
# and proportion of non-western, and eliminating the empty cells from the
# raster and distance matrix.
rast <- raster[[wijk]]
dist <- distances[[wijk]]
rast$density <- rast$pnw <- NA

for(cell in 1:nrow(rast)){
  residents <- subset(w, w$rdCoords == rast$rdCoords[cell])
  rast$density[cell] <- nrow(residents)
  rast$pnw[cell] <- sum(residents$g) / nrow(residents)
}
emptyCells <- which(rast$density == 0)

rast <- rast[-emptyCells,]
dist <- dist[-emptyCells, -emptyCells]
#

library("ggplot2")
library("gridExtra")

densityplot <- ggplot(rast, aes(x = east, y = north, color = density)) +
  geom_point(shape = 15, size = 2) +
  scale_color_viridis_c(option = "A") +
  #scale_x_continuous(
  #  limits = c(min(raster[[wijk]]$east), max(raster[[wijk]]$east))) +
  #scale_y_continuous(
  #  limits = c(min(raster[[wijk]]$north), max(raster[[wijk]]$north))
  #) + 
  theme(legend.position = "top")

pnwplot <- ggplot(rast, aes(x = east, y = north, color = pnw)) +
  geom_point(shape = 15, size = 2) +
  scale_color_viridis_c(option = "A") +
  #scale_x_continuous(
  #  limits = c(min(raster[[wijk]]$east), max(raster[[wijk]]$east))) +
  #scale_y_continuous(
  #  limits = c(min(raster[[wijk]]$north), max(raster[[wijk]]$north))
  #) + 
  theme(legend.position = "top")

grid.arrange(densityplot, pnwplot, ncol = 2)


#







################################################################################

sp::coordinates(w) = c("east","north")

library("rgdal")

buurtShape <- rgdal::readOGR(
  dsn = "./geoData/WijkBuurtkaart_2019_v3",
  layer = "buurt_2019_v3"
)
buurtShape <- sp::spTransform(buurtShape, CRS("+proj=longlat +ellps=WGS84"))

rast <- rgdal::readOGR(
  dsn = "./geoData/NL_vierkant_100meter_bij_100meter",
  layer = "NL_vierkant100m"
)
# NL_vierkant_100meter_bij_100meter
#temp <- rast

rast$E <- as.numeric(substr(
  as.character(rast$C28992R100),
  start = 2,
  stop = 5
))
rast$N <- as.numeric(substr(
  as.character(rast$C28992R100),
  start = 7,
  stop = length(as.character(rast$C28992R100))
))

rast$east <- (rast$E * 100) + 50 # (meters)
rast$north <- (rast$N * 100) + 50

# Transforming the cbs coordinate system to WGS-84
wgs84coords <- rd2wgs84(rast$east, rast$north)
rast$E <- wgs84coords[,2]
rast$N <- wgs84coords[,1]

rast <- rast@data
sp::coordinates(rast) = c("E","N")
llCRS <- sp::CRS("+proj=longlat +ellps=WGS84")
##rast <- sp::SpatialPoints(rast, proj4string = llCRS)
proj4string(rast) <- proj4string(buurtShape)

rastbu <- sp::over(rast, buurtShape)

rast$BU_CODE <- rastbu$BU_CODE
rast$BU_NAAM <- rastbu$BU_NAAM
rast$WK_CODE <- rastbu$WK_CODE
rast$GM_CODE <- rastbu$GM_CODE
rast$GM_NAAM <- rastbu$GM_NAAM
rast$POSTCODE <- rastbu$POSTCODE





test <- w$east[order(w$east)][c(1,2)]
test[2] - test[1]
test <- w$east[order(w$east, decreasing = TRUE)][c(5,6)]
test[1] - test[2]
mean(c(0.001138123, 0.001400141))

test <- w$north[order(w$north)][c(7,8)]
test[2] - test[1]
test <- w$north[order(w$north, decreasing = TRUE)][c(3,4)]
test[1] - test[2]
mean(c(0.00091, test[1] - test[2]))
51.90581 - 51.90672
51.90672 - 51.90762
51.92567 - 51.92643
mean(c(51.90581 - 51.90672, 51.90672 - 51.90762, 51.92567 - 51.92643))

c(0.001269132, 0.0008566667)

rast <- raster::raster(
  resolution = c(0.001269132, 0.0008566667),
  xmn = min(w$east), xmx = max(w$east),
  ymn = min(w$north), ymx = max(w$north),
  crs = "+proj=longlat +ellps=WGS84 +no_defs"
)
rast <- raster::setValues(rast, sample(1:10, size = length(rast), replace = T))
plot(rast)
raster::plot(rast)
points(w$east, w$north)

ggplot(rast) +
  geom_tile(aes(fill = value))


####
rast <- raster::raster(as.matrix(cbind(w$east, w$north)))

rast <- raster::raster(
  #x = as.matrix(cbind(w$east, w$north)),
  xmn = min(w$east), xmx = max(w$east),
  ymn = min(w$north), ymx = max(w$north),
  crs = "+proj=longlat +ellps=WGS84"#,
)

rast <- raster::raster(
  resolution = c(0.001269132, 0.0008566667),
  xmn = min(w$east), xmx = max(w$east),
  ymn = min(w$north), ymx = max(w$north),
  crs = "+proj=longlat +ellps=WGS84 +no_defs"
)

test <- raster::raster(nrows=108, ncols=21, xmn=0, xmx=10)



sp::coordinates(w) = c("east","north")
grd <- sp::points2grid(
  points = w, tolerance = 0.77
)

grd <- sp::SpatialPixels(w)

# makegrid





## spatial points for agents
sp::coordinates(w) = c("east","north")
wsp <- sp::SpatialPoints(
  w,
  proj4string = sp::CRS("+proj=longlat +ellps=WGS84")
)
#proj4string(wsp) <- sp::CRS("+proj=longlat +ellps=WGS84")


#vks <- subset(vks, vks$)
vks$E <- as.numeric(substr(
  as.character(vks$c28992r100),
  start = 2,
  stop = 5
))
vks$N <- as.numeric(substr(
  as.character(vks$c28992r100),
  start = 7,
  stop = length(as.character(vks$c28992r100))
))


# The coordinates of each square (vkd$E and vkd$N) point at their south-west
# corners. The format of these coordinates makes it easy to find the centroid:
vkd$east <- (vkd$E * 100) + 50 # (meters)
vkd$north <- (vkd$N * 100) + 50

# Transforming the cbs coordinate system to WGS-84
wgs84coords <- rd2wgs84(vkd$east, vkd$north)
vkd$E <- wgs84coords[,2]
vkd$N <- wgs84coords[,1]

# Make sure same coordinate system
buurtShape <- sp::spTransform(buurtShape, CRS("+proj=longlat +ellps=WGS84")) 
sp::coordinates(vkd) = c("E","N")









