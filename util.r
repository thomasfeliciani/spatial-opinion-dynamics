#  Utility functions

# This script contains general purpose, auxiliary functions
# for the main script, "simulation.r".
library("ggplot2")
library("ggpubr")
library("ggmap")
library("ggsn")


# Opinion formation models______________________________________________________
#
# This calculates the interaction weight (or 'similarity') between two
# interacting agents. "agents" refers to the dataframe containing the
# simulation agentset; "ego" and "alter" are two row indices that identify the
# pair of agents of which we wish to know the similarity; "H" is the model
# parameter.
computeWeight <- function(agents, ego, alter, H){
  return(
    1 - (abs(agents$opinion[ego] - agents$opinion[alter]) * H +
           abs(agents$group[ego] - agents$group[alter]) * (1 - H)) / 1
  )
}


# Function that returns the new opinion of the interacting agent(s), resulting
# from the interaction between ego and alter.
#
NIcomputeOpinion <- function(
  agents,
  ego,
  alter,
  H,
  typeInteraction = "two-way",
  rateOpinionChange = 1
) {
  
  # We store the interaction weight and the opinion difference
  w <- computeWeight(agents, ego, alter, H)
  opinionDiff <- abs (agents$opinion[alter] - agents$opinion[ego])
  
  # We run a convergence test
  hasConverged <- TRUE
  if (w != 0 & opinionDiff != 0 & opinionDiff != 2 ) hasConverged <- FALSE
  
  # We update the opinion of ego
  oEgo <- agents$opinion[ego] +
    rateOpinionChange * (agents$opinion[alter] - agents$opinion[ego]) * w / 2
  if (oEgo < -1) oEgo <- -1
  else if (oEgo > 1) oEgo <- 1
  
  # If interactions are two-way (i.e. if alter influences ego and at the same
  # time ego influences alter), then we also determine the new opinion of alter.
  if (typeInteraction == "two-way"){
    oAlter <- agents$opinion[alter] +
      rateOpinionChange * (agents$opinion[ego] - agents$opinion[alter]) * w / 2
    if (oAlter < -1) {oAlter <- -1}
    else if (oAlter > 1) {oAlter <- 1}
    return(list(value = c(oEgo, oAlter), hasConverged = hasConverged))
  } else {
    return(list(value = oEgo, hasConverged = hasConverged))
  }
}


# Functions that run the persuasive argument (PA) model
#
PAcomputeOpinion <- function(ego, alter){
  opinion <- ego$opinion + argument(ego$opinion, alter$opinion)
  if (opinion < -1) opinion <- -1
  else if (opinion > 1) opinion <- 1
  return(opinion)
}

# For the PA model, this function returns the effect of a "pseudo-argument":
# that is, the amount of influence on the opinion of ego that an interaction
# with alter would have produced, if alter had communicated an argument to ego.
argument <-function(opinion, j_opinion){
  if (rbinom(1, 1, (j_opinion + 1) / 2) == 1) {  
    # If j picks a pro argument...
    if (rbinom(1, 1, (opinion + 1) / 2) == 1) { 
      # ...and i drops a pro argument, then a=0 (ineffective argument exchange)
      return(0)
    } else {                                                
      # ...and i drops a con argument, then i's opinion gets a positive push
      return(2 / S)
    }
  } else {                                                     
    # If j picks a con argument
    if (rbinom(1, 1, (opinion + 1) / 2) == 1){
      # ...and i drops a pro argument, then i's opinion gets a negative push
      return(-2 / S)
    } else {                                                   
      # ...and i drops a con argument, then a=0 (ineffective argument exchange)
      return(0)
    }
  }
}


# Spatial tools ________________________________________________________________
#
# Moran's I
moranI <- function(x, y = NULL, proxmat, dens = NULL, N = length(x)) {
  # Adapted from Anselin (1995, eq. 7, 10, 11)
  # https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1538-4632.1995.tb00338.x
  #
  #dens: the proportion of individuals in each cell over the district population
  #if individual level data dens is.null and N is simply length of input
  #if we have aggregate data then N should be total population size (or actually
  #just a large number)
  if(is.null(y)){y <- x}
  if(is.null(dens)){dens <- rep(1/N, times = N)}
  
  #correct scaling of opinions for densities
  v1dens_ind <- rep(x, times = (dens * N))
  v1dens <- (x - mean(v1dens_ind))/sd(v1dens_ind)
  v2dens_ind <- rep(y, times = (dens * N))
  v2dens <- (y - mean(v2dens_ind)) / sd(v2dens_ind)
  
  # (density) weighted proximity matrix
  w <- proxmat
  wdens <- t(dens * t(w))
  wdens <- wdens / rowSums(wdens)
  
  # density and proximity weighted locals
  localI <- (v1dens * wdens %*% v2dens) #formula 7
  
  # correct the normalization constants
  m2 <- sum(v1dens^2 * dens)
  S0 <- N #we know the weight matrix for the individual level should add up to N
  ydens <- S0 * m2
  globalI <- sum(localI * dens * N) / ydens # formula 10/11
  
  return(list(
    globalI = globalI,
    localI = as.numeric(localI)
  ))
}


# Function to convert between coordinates systems
rd2wgs84 <- function(X, Y)
{
  #http://www.dekoepel.nl/pdf/Transformatieformules.pdf
  #basispunten definieren
  X0 <- 155000.00 
  Y0 <- 463000.00 
  j0 <- 52.15517440
  l0 <- 5.38720621
  #coefficienten definieren
  K01<- 3235.65389 
  K20<- -32.58297
  K02<- -0.24750 
  K21<- -0.84978 
  K03<- -0.06550 
  K22<- -0.01709 
  K10<- -0.00738 
  K40<- 0.00530 
  K23<- -0.00039 
  K41<- 0.00033 
  K11<- -0.00012 
  
  L10<- 5260.52916
  L11<- 105.94684
  L12<- 2.45656
  L30<- -0.81885
  L13<- 0.05594
  L31<- -0.05607
  L01<- 0.01199
  L32<- -0.00256
  L14<- 0.00128
  L02<- 0.00022
  L20<- -0.00022
  L50<- 0.00026
  
  dX <- (X - X0)*10^-5
  dY <- (Y - Y0)*10^-5 
  {
    j <- j0 + 
      (
        K01*dX^0*dY^1 +
          K02*dX^0*dY^2 +
          K03*dX^0*dY^3 +
          K10*dX^1*dY^0 +
          K20*dX^2*dY^0 +
          K21*dX^2*dY^1 +
          K22*dX^2*dY^2 +
          K23*dX^1*dY^3 +
          K40*dX^2*dY^0 +
          K41*dX^2*dY^1 
      )/3600
  }
  
  {
    l <- l0 + 
      (
        L10*dX^1*dY^0 +
          L11*dX^1*dY^1 +
          L12*dX^1*dY^2 +
          L30*dX^3*dY^0 +
          L13*dX^1*dY^3 +
          L31*dX^3*dY^1 +
          L01*dX^0*dY^1 +
          L32*dX^3*dY^2 +
          L14*dX^1*dY^4 +
          L02*dX^0*dY^2 +
          L20*dX^2*dY^0 +
          L50*dX^5*dY^0 
      )/3600
  }
  wgs84<-cbind(j,l)
  return(wgs84)
}


# Signal processing ____________________________________________________________
# These functions are useful for the artificial manipulation of the distribution
# across space of e.g. an ethnic group, or any other variable.

logAmplify <- function(x, steepness = 10, L = 1, x0 = 0.5) {
  # L is the max output value;
  # x0 is the sigmoid's midpoint.
  
  e = exp(1) # Euler's constant e
  return(
    L / (1 + (e ^ (- steepness * (x - x0))))
  )
}

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



# Graphics______________________________________________________________________
# This function plots the agents on a map and shows how groups and opinions
# are distributed across space.
plotMap <- function(agents, zoom = 15) {
  myMap <- ggmap::get_stamenmap( # downloads map tiles by Stamen Design 2021
    bbox = c(
      left = min(agents$y_coor),
      bottom = min(agents$x_coor),
      right = max(agents$y_coor),
      top = max(agents$x_coor)),
    maptype = "toner-background",
    crop = FALSE, zoom = zoom
  )
  
  mapTheme <- ggplot2::theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
  
  groupPlot <- ggmap(myMap, darken = c(0.5, "white")) + 
    geom_point( # group density (tile / square unit)
      aes(x = y_coor, y = x_coor, color = as.factor(group)),
      data = agents, size = 0.6, alpha = 0.4,
      position = position_jitter(width = 0.00065, height = 0.00045, seed = 1)
    ) +
    scale_color_manual(
      values = c("-1" = "darkorange", "1" = "darkmagenta"),
      labels = c("natives and western", "non-western")
    ) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    ggtitle("group")
  
  opinionPlot <- ggmap(myMap, darken = c(0.5, "white")) + 
    geom_point( # group density (tile / square unit)
      aes(x = y_coor, y = x_coor, color = opinion),
      data = agents, size = 0.6, alpha = 0.3,
      position = position_jitter(width = 0.00065, height = 0.00045, seed = 1)
    ) +
    scale_colour_gradientn(colors = c("red", "gray", "blue")) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    ggtitle("opinion")
  
  return(ggpubr::ggarrange(
    groupPlot + mapTheme, opinionPlot + mapTheme,
    ncol = 2, hjust = -1, common.legend = FALSE
  ))
}




# Miscellanea___________________________________________________________________
#
writeToCSV <- function(variablesToExport) {
  write.table(
    variablesToExport,
    file = "Raw output.csv",
    row.names = FALSE,
    na = "",
    col.names = FALSE,
    sep = ","
  )
}


normalize <- function(x) {
  min <- min(x, na.rm = TRUE)
  max <- max(x, na.rm = TRUE)
  ifelse(
    min - max == 0,
    x[!is.na(x)] <- 1, # if no variability in x, all values are assumed to be 1.
    x <- (x - min) / (max - min)
  )
  return(x)
}