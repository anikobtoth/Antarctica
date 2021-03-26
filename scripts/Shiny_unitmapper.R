library(shiny)
library(leaflet)
library(rgdal)

ui <- fluidPage(
  fileInput(inputId = "unit", label = "Choose file", multiple = TRUE, accept = c(".geoTiff")), 
  
  leafletOutput(outputId = "mymap", width = "700px", 
                height = "700px")
)

resolutions <- c(8192, 4096, 2048, 1024, 512, 256)
zoom <- 0
maxZoom <- 5

crsAntartica <- leafletCRS(
  crsClass = "L.Proj.CRS",
  code = "EPSG:3031",
  proj4def = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs",
  resolutions = resolutions,
  origin = c(-4194304, 4194304),
  bounds =  list( c(-4194304, -4194304), c(4194304, 4194304) ))

antarticaTilesURL <- "//map1{s}.vis.earthdata.nasa.gov/wmts-antarctic/MODIS_Aqua_CorrectedReflectance_TrueColor/default/2014-12-01/EPSG3031_250m/{z}/{y}/{x}.jpg"
border <- readOGR("../Data/Base", "Antarctic_landpoly") %>% spTransform(CRS("+init=epsg:4326"))


server <- function(input, output){
  
  
  output$mymap <- renderLeaflet({
    
   
    # initialise the leaflet object
    basemap <- leaflet(options = leafletOptions(
      crs = crsAntartica, minZoom = zoom, maxZoom = maxZoom, worldCopyJump = FALSE)) %>%
      setView(0, -90, 0) %>%
      addPolygons(data = border, color = "#ff0000", weight = 2, fill = FALSE) %>%
      #addCircleMarkers(data = points, label = ~Name) %>%
      addTiles(urlTemplate = antarticaTilesURL,
               layerId = "antartica_tiles",
               attribution = "<a href='https://earthdata.nasa.gov/gibs'> NASA EOSDIS GIBS</a>&nbsp;&nbsp;&nbsp; <a href='https://github.com/nasa-gibs/web-examples/blob/release/leaflet/js/antarctic-epsg3031.js'> View Source </a>",
               options = tileOptions(
                 tileSize = 512,
                 subdomains = "abc",
                 noWrap = TRUE,
                 continuousWorld = TRUE,
                 format = "image%2Fjpeg"
               )) %>%
      addGraticule() %>%
      htmlwidgets::onRender(
        "function(el, t){
       var myMap = this;
       debugger;
       var tileLayer = myMap.layerManager._byLayerId['tile\\nantartica_tiles'];

       // HACK: BEGIN
       // Leaflet does not yet handle these kind of projections nicely. Monkey
       // patch the getTileUrl function to ensure requests are within
       // tile matrix set boundaries.
       var superGetTileUrl = tileLayer.getTileUrl;

       tileLayer.getTileUrl = function(coords) {
         debugger;
         var max = Math.pow(2, tileLayer._getZoomForUrl() + 1);
         if ( coords.x < 0 ) { return ''; }
         if ( coords.y < 0 ) { return ''; }
         if ( coords.x >= max ) { return ''; }
         if ( coords.y >= max ) { return ''; }
           return superGetTileUrl.call(tileLayer, coords);
       };
       // HACK: END
    }")
    
    
    
  })
  
}

shinyApp(ui, server)
