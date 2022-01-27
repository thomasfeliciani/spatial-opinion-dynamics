# Shiny app to visualize how district features can be manipulated.

#rm(list = ls())

#source("util.r")
#source("../util.r")
#library("reshape2")
library("scales")
#library("sp")
library("geosphere")
library("ggplot2")
library("ggmap")
library("ggpubr")
library("shiny")


#load("./geoData/geoData.RData")
load("./geoData.RData")



# Utilities ____________________________________________________________________
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


# Functions for signal processing:
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


makeGradient_linear <- function (bearing, coords, steepness = 0) {
  linep <- geosphere::destPoint(
    p = c(median(coords$east), median(coords$north)), # approximate center
    b = c(bearing -5, bearing + 5), # bearing; tangent to the target line
    d = 30000 # Distance (meters) from center to outskirts. Just a big number.
  )
  
  x <- apply(
    as.matrix(cbind(coords$east, coords$north)),
    MARGIN = 1,
    FUN = function (p) {
      geosphere::dist2Line(
        p = p,
        line = linep,
        distfun = distMeeus
      )
    }
  ) [1,]
  
  #x <- (1 - normalize(x))
  x <- 1 - normalize(x) - steepness
  x[x < 0] <- 0
  x <- normalize (x)
  
  return(normalize(x))
  #return(normalize(x) * (steepness - 1) + 1)
}

makeGradient_radial <- function (hotspot, coords, steepness = 0) {
  x <- c()
  
  # Calculating the Euclidean distance of each square from the hotspot:
  x <- geosphere::distMeeus( # faster
    p1 = hotspot,
    p2 = cbind(coords$east, coords$north)
  )
  #for (c in 1:nrow(coords)) {
  #  distance = sqrt(sum((hotspot - c(coords$east[c], coords$north[c])) ^ 2))
  #  x[c] <- distance
  #}
  
  # Transforming the distance values so that (1) they are inverted (squares
  # closer to the hotspot get higher values); (2) they range between 0 and 1;
  # (3) we introduce a "steepness" parameter that determines how close to the
  # hotspot the values drop to zero.
  x <- 1 - ((x - min(x)) / ((max(x) - min(x)) * (1 - steepness)))
  #x <- 1 - ((x - min(x)) / (max(x) - min(x)) * steepness)#(1 - steepness))
  
  # Ensuring range is correct:
  x[x < 0] <- 0
  
  return(x)
}

normalize <- function(x) {
  min <- min(x, na.rm = TRUE)
  max <- max(x, na.rm = TRUE)
  return((x - min) / (max - min))
}





# UI ___________________________________________________________________________

ui <- fluidPage(
  
  fluidRow(column(
    12,
    mainPanel(
      width = 12,
      plotOutput(outputId = "densityMap", height = "450px")#,
      #textOutput(outputId = "globalI")
    )
  )),
  
  fluidRow(
    column(
      12,
      column(
        4, h2("Location"), style = "background-color:#e0e0e0;",
        selectInput(
          inputId = "place",
          label = "choose a location:",
          choices = c(
            "'s-Gravenhage",
            "'s-Hertogenbosch",
            "Almere",
            "Amsterdam",
            "Apeldoorn",
            "Breda",
            "Eindhoven",
            "Enschede",
            "Groningen",
            "Leeuwarden",
            "Leiden",
            "Lelystad",
            "Maastricht",
            "Nijmegen",
            "Rotterdam",
            "Terschelling",
            "Tilburg",
            "Utrecht",
            "Zwolle"
          ),
          selected = "Leiden"
        ),
        
        sliderInput(
          inputId = "resize", label = "sample size", min = 1, max = 100,
          post = "%", value = 100
        ),
        
        helpText(
          "The distance decay function is used to calculate the LISA scores
          (Local Indicator of Spatial Association), which depend on the
          definition of proximity betweneen points on the map.
          Use this parameter to define how steeply proximity decays over 
          the distance between points. As an indication, 10 is very steep, and
          1000 is very mild. The default value is 100."
        ),
        
        column(
          5, style = "background-color:#e0e0e0;",
          numericInput(
          inputId = "distanceDecaySlope", label = "distance decay parameter",
          value = 100, min = 1, max = 10000
          )
        ),
        
        column(
          7, style = "background-color:#e0e0e0;",
          helpText(
          "Note: due to memory requirements (and the long computing time)
          LISA scores will not be calculated for the larger locations."
          )
        ),
        
      ),
      
      column(
        4, h2("Rearrange groups"), style = "background-color:#cfcfcf;",
        
        helpText(
          "By increasing the overlay above 0% we can artificially reposition
          the two groups following a linear gradient. The angle allows to
          rotate the gradient to the desired position; and the steepness
            parameter defines how far the gradient radiates."
        ),
        
        column(
          12, style = "background-color:#cfcfcf;",
          
          sliderInput(
            inputId = "overlayGradient", label = "overlay gradient",
            min = 0, max = 100, post = "%", value = 0
          ),
          
          sliderInput(
            inputId = "bearing", label = "gradient angle",
            min = 0, max = 360, value = 45, post = "°"
          ),
          
          sliderInput(
            inputId = "hotspotSteepness", label = "gradient steepness",
            min = 0, max = 0.99, value = 0.5
          )
          
        )#,
        
        #column(
        #  6, "hotspot coordinates:", style = "background-color:#c6c6c6;",
        #  
        #  sliderInput(
        #    inputId = "hotspotLong", label = "longitude",
        #    min = 0, max = 100, post = "%", value = 80
        #  ),
        #  
        #  sliderInput(
        #    inputId = "hotspotLat", label = "latitude",
        #    min = 0, max = 100, post = "%", value = 80
        #  )
        #)
        
      ),
      
      column(
        4, h2("Amplification"), style = "background-color:#e0e0e0;",
        
        helpText(
          "Amplification applies a logisitic transformation to increase the
            density of the ethnic groups in each square. A group's density is
            increased where it's above 50% and decreased where it's below 50%."
        ),
        helpText(
          "Increasing amplification strength creates more homogeneous
            squares and sharper boundaries between local ethnic clusters."
        ),
        
        column(
          3, style = "background-color:#e0e0e0;",
          radioButtons(
            inputId = "amplification_", label = "amplify?",
            choices = c("on", "off"), selected = "off"
          )
        ),
        column(
          9, style = "background-color:#e0e0e0;",
          sliderInput(
            inputId = "amplSteepness", label = "amplification strength",
            min = 0, max = 20, value = 10
          )
        )
      )
    )
  )
)


# Server _______________________________________________________________________

server <- function(input, output) {
  
  if (FALSE) { ### for debugging
    input <- data.frame(
      place="Zwolle", resize=100, distanceDecaySlope=100,
      overlayHotspot=80,# hotspotLong=70, hotspotLat=70,
      bearing = 45, hotspotSteepness=0.5,
      amplification_="on", amplSteepness=10
    )
  } ### end debug tool
  
  world <- reactive({
    x <- subset(vk, vk$municipality == input$place)
    x$density <- x$density / 100 * input$resize
    
    # Transforming signal (x$pnw) as needed, i.e. by adding a hotspot or by
    # amplifying it.
    gradient <- makeGradient_linear(
      bearing = input$bearing,
      coords = x,
      steepness = input$hotspotSteepness
    )
    #hotspot <- makeGradient_radial(
    #  
    #  # (here we find the requested hotspot coordinates)
    #  hotspot = c(
    #    east = min(x$east) + (
    #      (max(x$east) - min(x$east)) / 100 * input$hotspotLong
    #    ), 
    #    north = min(x$north) + (
    #      (max(x$north) - min(x$north)) / 100 * input$hotspotLat
    #    )
    #  ),
    #  coords = x,
    #  steepness = input$hotspotSteepness
    #)
    
    # Then we overlay the original data and the hotspot:
    x$pnw <- sapply(
      1:length(x$pnw),
      FUN = function(i) weighted.mean(
        c(x$pnw[i], gradient[i]),
        w = c(1 - (input$overlayGradient / 100), input$overlayGradient / 100)
      )
    )
    
    
    # Last we amplify the ethnic composition signal, if needed. 
    # We do this transformation by using a logistic equation that increases
    # a group's density in squares where it's above 0.5, and decreasing it 
    # in squares where it is below 0.5.
    if (input$amplification_ == "on") x$pnw <- logAmplify(
      x$pnw, steepness = input$amplSteepness
    )
    
    
    # Based on the density (overall and of each group), we can now get the 
    # absolute number of residents for each group:
    x$pw <- 1 - x$pnw
    x$countNW <- x$density * x$pnw
    x$countW <- x$density * x$pw
    
    
    
    # Last, we calculate the LISA (local indicator of spatial association).
    # This requires the calculation of two very large matrices storing the 
    # distance (distmat) and resulting proximity (proxmat) between all pairs
    # of squares. Because of RAM and time constraint, we only do this for
    # sufficiently small locations:
    if (nrow(x) < 3000) { # Threshold is set by the number of squares.
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
      N <- sum(x$density)
      moranI <- moranI(
        x = x$pnw,
        proxmat = proxmat,
        dens = x$density / N,
        N = N
      )
      x$localI <- moranI$localI
      Istatistics <- paste(
        "Global I =", round(moranI$globalI, 3), "\nLISA scores:")
      
    } else {
      # Else, i.e. if there are too many squares, then we'd rather not 
      # calculate these spatial measures:
      x$localI <- NA
      Istatistics <- paste0(
        "Could not calculate LISA and Moran's I\n",
        "because the chosen location is too large."
      )
    }
    
    
    
    
    # Downloading map tiles
    mapTiles <- ggmap::get_stamenmap( # map tiles by Stamen Design 2021
      bbox = c(
        left = min(x$east),
        bottom = min(x$north),
        right = max(x$east),
        top = max(x$north)),
      maptype = "toner-background",
      crop = FALSE, zoom = 12#input$zoom
    )
    
    mapTheme <- ggplot2::theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(angle = 35, vjust = 1, hjust = 0.5),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
    
    pointerSize <- 0.13 / max( # only used to approximate the right tile size
      c(max(x$east) - min(x$east)),
      c(max(x$north) - min(x$north))
    )
    
    list(
      x = x,
      Istatistics = Istatistics,
      mapTiles = mapTiles,
      mapTheme = mapTheme,
      pointerSize = pointerSize
    )
  })
  
  
  
  
  output$densityMap <- renderPlot({
    out <- world()
    densityPlot <- ggmap(out$mapTiles, darken = c(0.4, "white")) + 
      geom_point( # group density (tile / square unit)
        aes(x = east, y = north, color = density),
        data = out$x, size = out$pointerSize, alpha = 0.7, shape = 15
      ) +
      scale_colour_viridis_c(
        option = "C", trans = scales::log10_trans()) +
      ggtitle("\nnumber of residents") +
      out$mapTheme
    
    nwPlot <- ggmap(out$mapTiles, darken = c(0.4, "white")) + 
      geom_point( # proportion non-western
        aes(x = east, y = north, color = pnw),#countNW),
        data = out$x, size = out$pointerSize, alpha = 0.9, shape = 15
      ) +
      scale_color_viridis_c(option = "C", limits = c(0,1)) +
      ggtitle("non-western migr. backgr.\n(proportion of residents)") +
      out$mapTheme
    
    nwCountPlot <- ggmap(out$mapTiles, darken = c(0.4, "white")) + 
      geom_point( # number of non-western
        aes(x = east, y = north, color = countNW + 0.01),
        data = out$x, size = out$pointerSize, alpha = 0.9, shape = 15
      ) +
      scale_color_viridis_c(
        option = "C",
        #limits = c(0, NA)#,
        breaks = c(0, 1, 10,100, 1000),
        trans = scales::log10_trans()
      ) +
      ggtitle("non-western migr. backgr.\n(number of residents)") +
      out$mapTheme
    
    if (! all(is.na(out$x$localI))) {
      LISAplot <- ggmap(out$mapTiles, darken = c(0.4, "white")) + 
        geom_point( # number of non-western
          aes(x = east, y = north, color = localI),
          data = out$x, size = out$pointerSize, alpha = 0.9, shape = 15
        ) +
        scale_color_viridis_c(option = "C") +
        ggtitle(out$Istatistics) + out$mapTheme
    } else {
      LISAplot <- ggmap(out$mapTiles, darken = c(0.9, "white")) + 
        ggtitle(out$Istatistics) + out$mapTheme
    }
    
    
    #wPlot <- ggmap(mapTiles, darken = c(0.4, "white")) + 
    #  geom_point( # group density (tile / square unit)
    #    aes(x = east, y = north, color = pw),
    #    data = x, size = pointerSize, alpha = 0.9, shape = 15
    #  ) +
    #  scale_color_viridis_c(option = "C", limits = c(0,1)) +
    #  ggtitle("native or western background\n(proportion of residents)")
    
    print(ggpubr::ggarrange(
      densityPlot,# + out$mapTheme,
      nwPlot,# + out$mapTheme,
      nwCountPlot,# + out$mapTheme,
      LISAplot,# + out$mapTheme,
      ncol = 4, hjust = -1, common.legend = FALSE
    ))
    
  })
  
}


shinyApp(ui = ui, server = server)
