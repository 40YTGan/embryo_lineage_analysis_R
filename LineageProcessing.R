#' LineageProcessing
#' This file contains several utility functions that helps processing data in CD-like dataFrames
#' @import dplyr
#' @import 
if(!require('dplyr')){stop("Package 'dplyr' could not be loaded")}

#' timeAvg averages the selected attribute for all cells of the same time, can add corresponding standard deviation attribute
#' @param datList list of CD-like dataFrames
#' @param attribute the attribute to be averaged
#' @return
#' @export
#'
#' @examples
timeAvg <- function(datList, attribute){
  datList |> group_by(time) |> 
      summarise(stdev:=sd(!!sym(attribute)), !!attribute:=mean(!!sym(attribute)))
}

#' alignTime: 'align' the cell time of a list of CD-like dataFrame by a specific cell stage 
#' so that there is a "common" reference point (like common start time of lineage)
#' @param datList list of CD-like dataFrames
#' @param alignCell the cell used to find common time point to align
#' @param alignPoint a 0 to 1 ratio number showing what ratio of the alignCell life to align together (1 means the end, 0 means the start)
#' @param alignTime the time of the `alignCell` at `alignPoint` after alignment, if the input is `"mean"` the align time will be the mean value of `alignPoint` across all embryos
#' @param alignBlot whether and method to align the blot values together default "no" for not align, "mean" for aligning the mean values together (by addition/subtraction)
#' @return modified list of CD-like dataFrames
#' @export
#' @examples `alignTime(datList, alignCell="MSa", alignPoint = 100, alignTime = "mean", alignBlot = "no")` 
#' to have every CD-like dataFrame have the end time of MSa cell equals the mean end time of MSa cell
#' all other cell time points will be moved in parallel
alignTime <- function(datList, alignCell, alignPoint = 1, alignTime = 0, alignBlot = "No"){
  rawTimes <- datList|>lapply(function(df){df[with(df, cell==alignCell),'time']|>quantile(alignPoint)})
  if(alignTime=='mean'){alignTime<-mean(unlist(rawTimes))}
  else if(alignTime=='max'){alignTime<-max(unlist(rawTimes))}
  embryos <- names(datList)
  times<-embryos |> lapply(function(embName){
    displace <- alignTime-rawTimes[[embName]] #move all times by the same amount
    datList[[embName]][,'time'] = datList[[embName]][,'time']+displace
  })
  if(alignBlot=="mean"){
    rawBlots <- datList|>lapply(function(df){df[with(df, expr = cell==alignCell),'blot']|>mean()})
    blotMean <- mean(unlist(rawBlots))
    blotDiff <- rawBlots|>lapply(function(val){blotMean-val})
    names(blotDiff)<-embryos
    datList <- embryos|>lapply(function(embName){
      df <- datList[[embName]]
      df$blot <- df$blot+blotDiff[[embName]]
      df
    })
    names(datList) <- embryos
  }
  for(i in 1:length(datList)){datList[[i]][,'time']<-times[[i]]}
  return(datList)
}

#' grepCells extract the data of given cells from the given CD table
grepCells <- function(CDData, cells=NULL, lineages=NULL, 
                      times = "ALL", timesDiffThreshold = 1.5, dataReturn = TRUE){
  selectCells <- integer(0)
  CellREs <- NULL
  if((is.null(cells) && is.null(lineages)) || is.null(times)){
    if (dataReturn) {return(CDData[0, ])}
    return(integer(0))
  }
  else if("ALL"%in%cells || "ALL"%in%lineages){
    selectCells<-CDData|>rownames()|>as.integer()}
  else{
    CellREs <- cells|> lapply(function(line){paste0("^", gsub("x", "[a-z]", line),"$")})
    CellREs <- c(CellREs, lineageRE(lineages))
  }
  for (RE in CellREs){selectCells <- union(selectCells, grep(RE, CDData$cell))}
  
  if(!"ALL"%in%times){#by default, no filter by times
    timeMatches <- NULL
    for(t in times){
      if (t %in% CDData$time) { # Exact match found
        tMatches <- which(CDData$time == t)
      } else {# Find nearest time
        timeDiffs <- abs(CDData$time - t)
        minDiff <- min(timeDiffs)
        if(minDiff>timesDiffThreshold){tMatches=NULL}
        else{tMatches <- which(timeDiffs-minDiff < 0.001)}
      }
      timeMatches <- union(timeMatches, tMatches)
    }
    selectCells <- intersect(selectCells,timeMatches)
  }
  if(dataReturn){return(CDData[selectCells, ])}
  else{return(selectCells)}
}

#' RePosition orientate one embryo at a given time point so that 
#'  the embry is centered at (0,0,0) coordinate
#'  the AB axis is aligned to the x-axis, DV to z-axis and RL to y-axis
#'
#' @param CDFrame a CD-like dataFrame
#' @param time 
#' @param indicatorP lineage on the relative posterior of embryo center 
#' @param indicatorD lineage on the relative dorsal of embryo center, use NULL not NaN if not considering
#' @param indicatorV lineage on the relative ventral of embryo center, use NULL not NaN if not considering
#' @param indicatorL lineage on the relative left of embryo center, use NULL not NaN if not considering
#' @param indicatorR lineage on the relative right of embryo center, use NULL not NaN if not considering
#' @param xSize 
#' @param ySize 
#' @param zSize 
#'
#' @return
#' @export
#'
#' @examples
RePosition <- function(CDFrame, time, indicatorP = "C", indicatorD = "Cxa",
                       indicatorV = "MSxxp", indicatorL = NULL, indicatorR = NULL,
                       xSize=0.08651, ySize=0.08651, zSize=0.5){
  CDFrame <- CDFrame |> grepCells(lineages = "ALL", times = time)
  CDFrame$x <- CDFrame$x*xSize
  CDFrame$y <- CDFrame$y*ySize
  CDFrame$z <- CDFrame$z*zSize
  #move the embryo to 0,0,0 by averaging the top 95% and bottom 5% values
  CDFrame$x <- CDFrame$x-mean(quantile(CDFrame$x, probs = c(0.05, 0.95)))
  CDFrame$y <- CDFrame$y-mean(quantile(CDFrame$y, probs = c(0.05, 0.95)))
  CDFrame$z <- CDFrame$z-mean(quantile(CDFrame$z, probs = c(0.05, 0.95)))
  
  transformMatrix <- rotationVec(CDFrame, indicatorP, indicatorD, indicatorV, indicatorL, indicatorR)
  
  coords <- as.matrix(CDFrame[, c("x", "y", "z")])
  transformed_coords <- coords %*% transformMatrix
  CDFrame$x <- transformed_coords[, 1]
  CDFrame$y <- transformed_coords[, 2]
  CDFrame$z <- transformed_coords[, 3]
  row.names(CDFrame)<- NULL
  return(CDFrame)
}

totalRePosition <- function(CDFrame, time, indicatorP = "C", indicatorD = "Cxa",
                            indicatorV = "MSxxp", indicatorL = NULL, indicatorR = NULL,
                            xSize=0.08651, ySize=0.08651, zSize=0.5){
  CDFrame$x <- CDFrame$x*xSize
  CDFrame$y <- CDFrame$y*ySize
  CDFrame$z <- CDFrame$z*zSize
  
  #align the embryo center to (0,0,0)
  positionAgg_FUN <- function(X){mean(quantile(X, probs = c(0.05, 0.95)))}
  centerCord <- aggregate(cbind(x,y,z)~time, data = CDFrame, FUN = positionAgg_FUN)
  cellCount <- count(CDFrame, time)
  t_start <- min(cellCount[cellCount[,"n"]>=50,"time"]) #embryos with >50 cells have more stable center coordination calculated by this method
  coord_start <- centerCord[centerCord[,"time"]==t_start,c("x","y","z")]
  centerCord[centerCord[,"time"]<t_start,c("x","y","z")] <- coord_start
  CDFrame <- left_join(CDFrame, centerCord, by = "time", suffix = c("", "_c"))
  CDFrame[,c("x","y","z")] <- CDFrame[,c("x","y","z")] - CDFrame[,c("x_c","y_c","z_c")]
  CDFrame <- CDFrame[,!names(CDFrame)%in%c("x_c","y_c","z_c")]
  CD_timeRef <- CDFrame |> grepCells(lineages = "ALL", times = time)
  
  transformMatrix <- rotationVec(CD_timeRef, indicatorP, indicatorD, indicatorV, indicatorL, indicatorR)
  coords <- as.matrix(CDFrame[, c("x", "y", "z")])
  transformed_coords <- coords %*% transformMatrix
  CDFrame$x <- transformed_coords[, 1]
  CDFrame$y <- transformed_coords[, 2]
  CDFrame$z <- transformed_coords[, 3]
  return(CDFrame)
}


rotationVec <- function(CDFrame, indicatorP, indicatorD, indicatorV, indicatorL, indicatorR){
  # 1) find the AB axis for this embryo
  pca1 <- prcomp(CDFrame[,c('x','y','z')]) #prcomp() always returns unit vectors
  APVect <- pca1$rotation[c(1,2,3)] #get the PC1 transformation UNIT vector
  PCells <- grepCells(CDFrame, lineages = indicatorP) #make sure that positive value for principle component 1 transformation vector is toward the A side
  PCellsVect <- c(mean(PCells$x), mean(PCells$y), mean(PCells$z))
  if(sum(PCellsVect*APVect)>0){APVect <- -APVect}
  # 2) find the unit transformation vector toward D side
  DVVect <- c(0,0,0)
  if(!is.null(indicatorD)){ #directional contribution from each indicator lineage
    DCells <- grepCells(CDFrame, lineages = indicatorD)
    DVVect_D <- c(mean(DCells$x), mean(DCells$y), mean(DCells$z)) |> OrthgonalProj(APVect)
    DVVect <- DVVect + DVVect_D
  }
  if(!is.null(indicatorV)){
    VCells <- grepCells(CDFrame, lineages = indicatorV) 
    VDVect_V <- c(mean(VCells$x), mean(VCells$y), mean(VCells$z)) |> OrthgonalProj(APVect)
    DVVect <- DVVect - VDVect_V
  }
  if(!is.null(indicatorR)){
    RCells <- grepCells(CDFrame, lineages = indicatorR)
    DVVect_R <- c(mean(RCells$x), mean(RCells$y), mean(RCells$z)) |> cross_product(APVect)
    DVVect <- DVVect + DVVect_R
  }
  if(!is.null(indicatorL)){
    LCells <- grepCells(CDFrame, lineages = indicatorL)
    VDVect_L <- c(mean(LCells$x), mean(LCells$y), mean(LCells$z)) |> cross_product(APVect)
    DVVect <- DVVect - VDVect_L
  }
  #do not scale before adding together so that lineages further from AP make more contributions
  DVVect <- DVVect/sqrt(sum(DVVect^2))#scale into a unit vector
  # 3) use vector cross product to get RL transformation vector (+ side is the right side)
  RLVect <- cross_product(APVect, DVVect)
  transformMatrix <- cbind(APVect,RLVect,DVVect)
  return(transformMatrix)
}
#'Gram-Schmidt Process gives the orthogonal component to AP axis
OrthgonalProj <- function(v,u){
  vOrth <- v - (sum(v*u)*u)
}
#' vector cross-product for orthogonal vector
cross_product <- function(v, u) {
  c(v[2] * u[3] - v[3] * u[2],
    v[3] * u[1] - v[1] * u[3],
    v[1] * u[2] - v[2] * u[1])
}

#' blotRange return the range (min, max) of blot values in a given CD-like dataframe
#' @param datList 
#' @return a named list with a value 'max' and a value 'min' 
#' @export
blotRange <- function(datList){
  maxBlot <- datalist$blot|>max()
  minBlot <- datalist$blot|>min()
  list(max = maxBlot, min = minBlot)
}

#' aggEmbryos merge a list of CD-like dataFrame into one dataFrame, with options to keep embryo identifier
#' @param datList list of CD-like dataFrames
#' @param labelEmbryo Boolean, whether to add a column to identify the embryo each cell come from
#' @return a dataframe combining all CD-like dataframes in the first dimension
#' @export
aggEmbryos <- function(datList, labelEmbryo = TRUE){
  if(!is.data.frame(datList[[1]])){ #will do aggregation if a LargeList of DataFrame(s) is delivered
    return(datList)
  }
  if(labelEmbryo){
    out <- lapply(names(datList), 
                  function(name) {
      datList[[name]] |> mutate(embryo = name)
    })  
  }
  else{out<-datList}
  out|>bind_rows()
}

#' blotAvgFun get the average blot value for each cell, can choose a cutoff to select for top blot values only
#' @param datList 
#' @param cutoff 
blotAvgFun <- function(datList, cutoff = 0.5){
  CD_merged <- aggEmbryos(datList, labelEmbryo = TRUE)
  output <- CD_merged |>
    group_by(cell, embryo) |>
    summarise(
      avgBlot = mean(blot[blot >= quantile(blot, cutoff)], na.rm = TRUE)
    )
  output
}


#' lineageNames returns all possible starts of daughter lineages (including the mother cell itself) under a given lineage/mother cell
#'
#' @param lineage a string or list of lineage name / mother cell name
#'
#' @return a list of regular expressions for cell
#' @export
#'
#' @examples lineageNames(c("EMS", "C"))
lineageNames <- function(lineage){
  P4 <- c("P4", "Z2", "Z3")
  P3 <- c("P3", "D", P4)
  P2 <- c("P2", "C", P3)
  EMS <- c("EMS", "MS", "E")
  E <- c("E$", "Ex")
  P1 <- c("P1", EMS, P2)
  P0 <- c("AB", P1)
  lineageDict <- list("P4" = P4, "P3" = P3, "P2" = P2, "P1" = P1, "P0" = P0, "EMS" = EMS, "E" = E)
  out<- lineage
  if(lineage %in% names(lineageDict)){out<-unlist(lineageDict[lineage])}
  else{out <- lineage}
  out
}

#' lineageRE returns all regular expressions needed to grep cells under a given lineage/mother cell
#'
#' @param lineage a string or list of lineage name / mother cell name, x as single-letter wildcard
#'
#' @return a list of regular expressions for cell
#' @export
#' @examples lineageRE(c("EMS", "C"))
lineageRE <- function(lineage){
  lineage<- lineage |> lapply(lineageNames)|>unlist()
  REs <- lineage |> lapply(function(line){paste0("^", gsub("x", "[a-z]", line))})
  return(REs)
}