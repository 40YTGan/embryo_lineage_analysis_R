#' drawEmb
#' This file contains several functions relevant to drawing a 3D (interactive) embryo, 
#' able to represent expression levels and highlight lineages
#' if deliver multiple "lineages" must also deliver a list "lineageColors" with elements corresponding to "lineages"
#' 
#' @import plotly
source('LineageProcessing.R')
if(!require('plotly')){stop("Package 'plotly' could not be loaded")}
if (!requireNamespace("webshot", quietly = TRUE)) {
  stop("Package 'webshot' is required. Install with install.packages('webshot')")
}
if (!requireNamespace("magick", quietly = TRUE)) {
  stop("Package 'magick' is required. Install with install.packages('magick')")
}
if (!requireNamespace("animation", quietly = TRUE)) {
  stop("Package 'animation' is required. Install with install.packages('animation')")
}

#' drawEmbExp 
#'
#' @param embryoCD the embryo dataframe with cell, time, and blot(expression) column
#' @param time which time to plot the embryo
#' @param lineages a list of lineages to highlight (plot with the given symbol in "shapes" parameter)
#' @param shapes a list of shapes to plot each lineage given in "lineages"
#' @param ReporterForAll color the nucleus points by expression or not (default TRUE)
#' @param colorScheme a list of at least 2 vectors made of a number between 0 and 1 to specify expression fraction, and a color value to specify color fro that fraction
#' @param maxBlot upper limit for expression coloring 
#' @param minBlot lower limit for expresion coloring
#' @param xSize size of each x unit in embryoCD dataframe
#' @param ySize size of each y unit in embryoCD dataframe
#' @param zSize size of each z unit in embryoCD dataframe
#' @param whiteSpace relative extra space on the 3 dimensions to leave beyond the actual range of cell loactions
#' @param center which point in the 3d space does the camera aim at
#' @param viewPoint where to put the imaginative camera initially
#'
#' @return a list made of: the plotly figure, a dataframe of the highlighted cells, and a dataframe of cells not highlighted
#' @export
#'
#' @examples
drawEmbExp<-function(embryoCD, time, lineages=c("ABa", "ABp","MS", "E", "C", "D", "P4"), 
                     shapes = NULL, 
                     ReporterForAll=F, colorScheme = NULL, 
                     maxBlot = NULL, minBlot = NULL, 
                     xSize=0.08651, ySize=0.08651, zSize=0.5, whiteSpace = c(0.05,2), 
                     center = list(x=0,y=0,z=0), viewPoint = list(x=0,y=0,z=1.8),
                     showAxis = F){
  embDat <- embryoCD|>grepCells(lineages = "ALL", times = time)
  row.names(embDat)<- NULL
  embDat$x <- embDat$x*xSize
  embDat$y <- embDat$y*ySize
  embDat$z <- embDat$z*zSize
  expressions <- aggregate(embryoCD[,"blot"], by = list(embryoCD[,"cell"]), FUN = "mean")
  names(expressions)<-c("cell", "blot")
  embDat[,"blot"]<-NULL
  embDat<- left_join(embDat, expressions, by = join_by(cell))
  if(is.null(maxBlot)){maxBlot <- max(embDat$blot)}
  if(is.null(minBlot)){minBlot <- min(embDat$blot)}
  if(identical(lineages, NULL)){lineages <- list("ABa", "ABp","MS", "E", "C", "D", "P4")} #default lineages
  if(length(shapes) != length(lineages)){
    print("\'shapes\' parameters not properly specified, using circles for all")
    # shapesSet <- schema(F)$traces$scatter$attributes$marker$symbol$values
    # shapesSet <- grep("[a-z]", shapesSet, value = TRUE)
    # shapesSet <- shapesSet[-grep("dot", shapesSet)]
    # shapesSet <- shapesSet[-grep("open", shapesSet)]
    # shapes <- seq(from = 1, by = 1, length.out = length(lineages))
    # shapes = shapesSet[shapes]
    shapes <- rep("circle", times = length(lineages))
  }
  if(is.null(colorScheme)){colorScheme =  "Viridis"}
  else if(colorScheme=="blackBlueOrange"){
    colorScheme<-list(c(0,rgb(0,0,0)), c(0.3,rgb(0,0.2,1)), c(0.7,rgb(1,0,0.2)), c(1,rgb(1,0.7,0.3)))}
  else if(colorScheme=="blackOrange"){colorScheme<-list(c(0,rgb(0,0,0)), c(0.6,rgb(1,0.1,0.1)), c(1,rgb(1,0.7,0.3)))}
  
  # rangeDim <- function(dimVals){
  #   span <- max(dimVals)-min(dimVals)
  #   return(span*c(-whiteSpace[1],whiteSpace[1]) + c(-whiteSpace[2],whiteSpace[2]) + range(dimVals))
  # }
  # xrange <- rangeDim(embryoCD[,"x"])
  # yrange <- rangeDim(embryoCD[,"y"])
  # zrange <- rangeDim(embryoCD[,"z"])
  if(showAxis){
    xtitle <- "x"
    ytitle <- "y"
    ztitle <- "z"
  }
  else{
    xtitle <- ""
    ytitle <- ""
    ztitle <- ""
  }
  fig<-plot_ly() |> #initiate a plotly figure object
    layout(
      scene = list(
        aspectmode = "data", 
        xaxis = list(title = xtitle, showgrid = showAxis, showticklabels=showAxis),
        yaxis = list(title = ytitle, showgrid = showAxis, showticklabels=showAxis),
        zaxis = list(title = ztitle, showgrid = showAxis, showticklabels=showAxis), 
        camera = list(eye = viewPoint, center=center, up = list(x = 0, y = 1, z = 0)
                      #,projection=list(type="orthographic")
        )), 
      paper_bgcolor = rgb(1,1,1)
    )
  
  selectCells <- NULL
  for (i in seq_along(lineages)) { #plot highlighted cells
    lineage <- lineages[[i]]
    thisSymbol <- shapes[[i]]
    thisCells <- grepCells(CDData = embDat, lineages = lineage, dataReturn = F)
    fig <- fig%>%AddGroupExp(groupName = lineage, data = embDat, selectedCells = thisCells, symbol = thisSymbol, 
                          colorScheme = colorScheme, colorMin = minBlot, colorMax = maxBlot)
    selectCells <- union(selectCells, thisCells)
  }
  
  otherCells <- embDat|>rownames()|>as.integer()
  otherCells <- otherCells[!otherCells%in%selectCells]
  if(is.null(lineages)){alpha <- 1}
  else{alpha <- 0.4}
  fig<-fig%>%add_trace(name="other cells",
                       x=embDat[otherCells,"x"], 
                       y=embDat[otherCells,"y"], 
                       z=embDat[otherCells,"z"], #cells that are not selected will be plotted transparent
                       type = "scatter3d", mode="markers",
                       marker = list(color=embDat[otherCells,"blot"], size = 6, opacity = alpha, 
                                     symbol = "circle-open",
                                     #line=list(color = "grey", width=3), 
                                     colorscale=colorScheme, cmin=minBlot, cmax=maxBlot)
  )
  
  return(list(fig, embDat[selectCells,], embDat[otherCells,]))
}

AddGroupExp <-function(plotlyFig, groupName, data, selectedCells, symbol, 
                   colorScheme, colorMin, colorMax, opacity=1){
  plotlyFig<-plotlyFig|>add_trace(name=groupName,
             x=data[selectedCells,"x"], 
             y=data[selectedCells,"y"], 
             z=data[selectedCells,"z"], #cells that are not selected will be plotted transparent
             type = "scatter3d", mode="markers",
             marker = list(color=data[selectedCells,"blot"], size = 6, opacity = 1, 
                           symbol = symbol, 
                           colorscale=colorScheme, cmin=colorMin, cmax=colorMax, 
                           colorbar = list(title="expression", x = 0)
             )
  )
  return(plotlyFig)
}

drawEmbLine <- function(embryoCD, time, lineages=NULL, colors = NULL, opacitys = NULL,
                        xSize=0.08651, ySize=0.08651, zSize=0.5, 
                        otherOpacity = 0.2, cellSize = 10,
                        center = list(x=0,y=0,z=0), viewPoint = list(x=0,y=0,z=1.8), 
                        showAxis = F){
  embDat <- embryoCD|>grepCells(lineages = "ALL", times = time)
  row.names(embDat)<- NULL
  embDat$x <- embDat$x*xSize
  embDat$y <- embDat$y*ySize
  embDat$z <- embDat$z*zSize
  if(identical(lineages, NULL)){lineages <- list("ABa", "ABp","MS", "E", "C", "D", "P4")} #default lineages
  traceCount <- length(lineages)
  if(length(colors) != traceCount){
    print("\'colors\' argument not properly specified")
    library("viridis")
    colors <- viridis_pal(option = "H")(traceCount)
  }
  if(length(opacitys) != traceCount){
    print("\'opacitys\' argument not properly specified")
    opacitys <- rep(1, traceCount)
  }
  if(showAxis){
    xtitle <- "x"
    ytitle <- "y"
    ztitle <- "z"
  }
  else{
    xtitle <- ""
    ytitle <- ""
    ztitle <- ""
  }
  fig<-plot_ly() |> #initiate a plotly figure object
    layout(
      scene = list(
        aspectmode = "data", 
        xaxis = list(title = xtitle, showgrid = showAxis, showticklabels=showAxis),
        yaxis = list(title = ytitle, showgrid = showAxis, showticklabels=showAxis),
        zaxis = list(title = ztitle, showgrid = showAxis, showticklabels=showAxis), 
        camera = list(eye = viewPoint, center=center, up = list(x = 0, y = 1, z = 0)
                      #,projection=list(type="orthographic")
                      )), 
      paper_bgcolor = rgb(1,1,1)
    )
  selectCells <- NULL
  for (i in seq_along(lineages)) { #add each lineage as a trace
    lineage <- lineages[[i]]
    color <- colors[[i]]
    opacity <- opacitys[[i]]
    thisCells <- grepCells(CDData = embDat, lineages = lineage, dataReturn = F)
    fig <- fig|>AddGroupLine(groupName = lineage, 
                             data = embDat, selectedCells = thisCells,
                             color=color, opacity = opacity, cellSize=cellSize)
    selectCells <- union(selectCells, thisCells)
  }
  otherCells <- embDat|>rownames()|>as.integer()
  otherCells <- otherCells[!otherCells%in%selectCells]
  fig<-fig%>%add_trace(
    name="other cells",
    x=embDat[otherCells,"x"], 
    y=embDat[otherCells,"y"], 
    z=embDat[otherCells,"z"], #cells that are not selected will be plotted transparent
    type = "scatter3d", mode="markers",
    marker = list(color="black", size = cellSize, opacity = otherOpacity, symbol = "circle")
  )
  return(list(fig, embDat[selectCells,], embDat[-selectCells,]))
}

AddGroupLine <- function(plotlyFig, groupName, data, selectedCells, color, opacity, cellSize){
  out <- plotlyFig|>add_trace(
    name = groupName, 
    x=data[selectedCells,"x"], 
    y=data[selectedCells,"y"], 
    z=data[selectedCells,"z"], 
    type = "scatter3d", mode="markers",
    marker = list(color=color, size = cellSize, opacity = opacity, symbol = 'circle')
  )
 out
}

#' printEmbImg Generates a .png image from drawEmbExp or drawEmbLine output figure
#' Require Reticulate that connects to a python distribution with "kaleido" and "plotly" installed
#'  recommend setting up a python venv for this and specify the python kernel by:
#'  reticulate::use_virtualenv("the/directory/holding/python/executable")
saveEmbImg <- function(
    fig, output_file = "embryo_DV_view.png",center = list(0,0,0),
    legend = F, colorscale = F,
    width = 600, height = 450, keepLabels = F
    ){
  fig <- fig|>
    layout(
      margin = list(l = 0, r = 0, t = 0, b = 0),
      scene = list(domain = list(x = c(0, 1), y = c(0, 1)))
    )
  if(!keepLabels){fig<-fig|> hide_colorbar() |> hide_legend()}
  plotly::save_image(fig, file = output_file, width = width, height = height, scale = 2)
}

#' createRotationGIF
#' NOT IMPLEMENTED!!
#' Generates a GIF showing rotation around the x-axis.
#' need Reticulate that connects to a python distribution with "kaleido" and "plotly" installed
#' @param fig The plotly figure object generated by drawEmbExp or drawEmbLine.
#' @param output_file Path to save the output GIF. Default: "animation.gif".
#' @param radius Camera distance from the center. Default: 2.5.
#' @param num_frames Number of frames in the animation. Default: 36.
#' @param duration Total duration of the GIF in seconds. Default: 5.
#' @param center the coordinate that camera rotates around, default (0,0,0)
#' @param width Width of the output GIF in pixels. Default: 800.
#' @param height Height of the output GIF in pixels. Default: 800.
#' @param dpi DPI for the output images. Default: 150.
#' 
#' @note Requires 'webshot', 'magick', and 'animation' packages. Ensure ImageMagick and PhantomJS are installed and accessible in the system PATH.
#'
# RotatingGIF <- function(fig, output_file = "animation.gif", 
#                               radius = 2.5, num_frames = 36, 
#                               duration = 5, center = list(0,0,0),
#                               width = 800, height = 800, dpi = 150) {
#   # Create temporary directory for frames
#   temp_dir <- tempdir()
#   dir.create(temp_dir, showWarnings = FALSE)
#   # Generate angles for rotation
#   angles <- seq(0, 2 * pi, length.out = num_frames)
#   frame_files <- vector(mode='list', length=length(angles))
#   
#   # Generate each frame
#   for (i in seq_along(angles)) {
#     theta <- angles[i]
#     eye_x <- center["x"]
#     eye_y <- radius * sin(theta)
#     eye_z <- radius * cos(theta)
#     fig_upd <- fig %>% layout(
#       scene = list
#       (camera = list(
#         center = center,
#         eye = list(x = eye_x, y = eye_y, z = eye_z),
#         up = list(x = 0, y = 0, z = 1)
#         )
#       )
#     )
#     # Save frame as PNG
#     if(i<10){fileNO<-paste0("00",i)}
#     else if(i<100){fileNO<-paste0("0",i)}
#     else{fileNO<-paste0("0",i)}
#     frame_files[i]<-paste0(temp_dir,"/img",fileNO,".png")
#     plotly::save_image(
#       fig_upd, 
#       file = frame_files[i], 
#       width = width, 
#       height = height, 
#       scale = dpi / 72
#     )
#   }
#   # Use magick to create GIF
    #saw memory overflow here
#   imgs <- list.files(temp_dir, full.names = TRUE)
#   img_list <- lapply(imgs, magick::image_read)
#   img_joined <- magick::image_join(img_list)
#   img_animated <- magick::image_animate(img_joined, fps = num_frames/duration)
#   magick::image_write(image = img_animated, path = output_file)
#   message("GIF saved to ", output_file)
#   file.remove(frame_files)
#   return(img_animated)
# }