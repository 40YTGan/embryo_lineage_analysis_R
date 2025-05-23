library(tidyverse)
library(ggtree)
library("ggplot2")

plot_trees <- function(inDir, root = "P0", 
                       plotExpression = T, lowerBound = 0.4, upperBound = 0.97, 
                       min_gain=100, max_gain = 10000, transform = "identity",
                       start_time = FALSE, end_time=FALSE,
                       tip_lab=TRUE,node_lab=FALSE, outDir = NaN, 
                       exp_col="blot", cell_col="cell", time_col="time",
                       width = 7, height = 7, dpi = 300){
  if(is.na(outDir)){outDir = inDir}
  files <- list.files(inDir, pattern=".+\\.csv$")
  embNames <- gsub("(.*)\\.csv$", "\\1", files)
  for (i in seq_along(files)) {
    directory <- paste0(inDir,files[i])
    tree <- CD_tree_plot(directory, root = root, transform = transform,
                         start_time = start_time, end_time = end_time, 
                         exp_col=exp_col, cell_col=cell_col, time_col=time_col,
                         plotExpression = plotExpression, min_gain=min_gain,max_gain=max_gain,  
                         tip_lab=tip_lab,node_lab=node_lab)
    ggsave(filename = paste0(outDir, embNames[i], root,".png"), plot = tree, 
           width = width, height = height, dpi = dpi, units = "in", device='png')
  }
}

GetSister <- function(x){
    if(is.na(x)) return(NA)		
    if(substr(x,nchar(x),nchar(x))=="a"){
        return(paste(substr(x,1,nchar(x)-1),"p",sep=""))
        }
    if(substr(x,nchar(x),nchar(x))=="d"){
        return(paste(substr(x,1,nchar(x)-1),"v",sep=""))
        }
    if(substr(x,nchar(x),nchar(x))=="l"){
        return(paste(substr(x,1,nchar(x)-1),"r",sep=""))
        }	
    if(substr(x,nchar(x),nchar(x))=="p"){
        return(paste(substr(x,1,nchar(x)-1),"a",sep=""))
        }
    if(substr(x,nchar(x),nchar(x))=="v"){
        return(paste(substr(x,1,nchar(x)-1),"d",sep=""))
        }
    if(substr(x,nchar(x),nchar(x))=="r"){
        return(paste(substr(x,1,nchar(x)-1),"l",sep=""))
        }	



    if(x=="AB") return("P1")
    if(x=="EMS") return("P2")
    if(x=="C") return("P3")
    if(x=="D") return("P4")
    if(x=="Z2") return("Z3")
    if(x=="E") return("MS")

    if(x=="P1") return("AB")
    if(x=="P2") return("EMS")
    if(x=="P3") return("C")
    if(x=="P4") return("D")
    if(x=="Z3") return("Z2")
    if(x=="MS") return("E")

    return(NA)
}

GetParent <- function(x){

    if(is.na(x) || x=="P") return(NA)
    if(is.element(substr(x,nchar(x),nchar(x)),c("a","p","d","v","l","r"))) return (substr(x,1,nchar(x)-1))
    if(x=="AB" || x=="P1") return("P0")
    if(x=="EMS" || x=="P2") return("P1")
    if(x=="E" || x=="MS") return("EMS")
    if(x=="C" || x=="P3") return("P2")
    if(x=="D" || x=="P4") return("P3")
    if(x=="Z2" || x=="Z3") return("P4")
    return("P")
}

##recursively create Newick format string from complete data frame (columns: 
MakeNewick <- function(data, root="P0", outer=TRUE,verbose=FALSE, tiplab = TRUE, nodelab=FALSE){
    if(verbose){
        message(root)
    }
    if(data[root,"hasDaughter"]){
        daughters <- sort(rownames(data[data[,"parents"]==root,]))
        if(outer){
            return(paste("((",MakeNewick(data,daughters[1],outer=FALSE),",",MakeNewick(data,daughters[2],outer=FALSE),")",root,":",data[root,"length"],");", sep=""))
        }else{
            return(paste("(",MakeNewick(data,daughters[1],outer=FALSE),",",MakeNewick(data,daughters[2],outer=FALSE),")",root,":",data[root,"length"], sep=""))
        }
        
        ##        message(paste("(",MakeNewick(data,daughters[1]),",",MakeNewick(data,daughters[2]),")",root,data[root,"length"], sep=""))
    }else{
        return(paste(root,data[root,"length"],sep=":"))
        ##        message(paste(root,data[root,"length"],sep=":"))
    }    
}

##arguments to add:
## color scale choice (or customize after the fact)
## error bars??? maybe this is an option for a separate related function

##example: Waterston/Murray ACD
##CD_tree_plot("~/Dropbox/Epic/ACD/ACD20080418_pha-4_3E3C5_1.csv", min_gain=0, max_gain=8000, end_time=300, node_lab=FALSE)
CD_tree_plot <- function(file, root="P0", 
                         min_gain=-1000, max_gain=FALSE, transform = "identity",
                         start_time = FALSE, end_time=FALSE, 
                         plotExpression=TRUE, lowerBound = 0.4, upperBound = 0.97,
                         exp_col="blot", tip_lab=TRUE,node_lab=FALSE, cell_col="cell", time_col="time", 
                         extra_branch_l=15, exp_legend="Expression"){

    CD <- read.csv(file, header=T)
    if(is.numeric(start_time)){
      CD <- CD[CD$time >= start_time,]
    }
    if(is.numeric(end_time)){
      CD <- CD[CD$time <= end_time,]
      X_Max <- (end_time-start_time)+15
    }
    else{X_Max <- NA}

    if(cell_col != "cell"){
        colnames(CD)[colnames(CD) == cell_col] <- "cell"
    }
    if(time_col != "time"){
        colnames(CD)[colnames(CD) == time_col] <- "time"
    }
    CD <- CD|>drop_na(time)
    
    CD <- CD[order(CD$cell),]

    cells<-unique(CD$cell)
    birthTime <- aggregate(CD$time,by=list(CD$cell), FUN=min)[,2]
    names(birthTime) <- cells
    endTime <- aggregate(CD$time,by=list(CD$cell), FUN=max)[,2]
    names(endTime) <- cells
    length <- (endTime-birthTime) + 1
    parents <- sapply(cells,GetParent)
    hasParent <- cells %in% cells[parents %in% cells]
    names(hasParent) <- cells

    ##debugging code
    ##message(paste("Eal", birthTime["Eal"], endTime["Eal"], length["Eal"], "Ear", birthTime["Ear"], endTime["Ear"],length["Ear"],sep=" "))

    ##add dummy branches to connect trees to root P0 or P. Each gets the same minute/tp length
    while(sum(!hasParent & !(cells %in% c("P0", "P"))) > 0){
        message("ADDING PARENTS FOR:", sum(!hasParent), "CELLS", sep=" ")
        message(paste(cells[!hasParent], collapse=" "))
        newCells <- unique(parents[cells[!hasParent]])
        cells <- c(cells, newCells)
        length <- c(length, rep(extra_branch_l,length(newCells)))
        parents <- c(parents, sapply(newCells,GetParent))    
        names(length) <- cells
        names(parents) <- cells
        hasParent <- cells %in% cells[parents %in% cells]
        names(hasParent) <- cells
    }

    hasDaughter <- cells %in% parents
    names(hasDaughter) <- cells

    ##QC
    sisters <- sapply(cells,GetSister)
    hasSister <- cells %in% sisters
    message(paste(sum(!hasSister),"cells don't have sisters - should equal one for P0"))

    tree_table <- data.frame(parents,sisters,length,hasDaughter)
    newick_tree <- MakeNewick(tree_table, root=root)
 ##   message(newick_tree)
    ## load into ggtree
    tree <- read.tree(text=newick_tree)
    
##    message(exp_col)
##    message(colnames(CD))
    if(plotExpression){
        funAvg <- function(X, lower, upper){
          lowerV <- quantile(X, lower)
          upperV <- quantile(X, upper)
          Xfiltered <- X[X>=lowerV[[1]]]
          Xfiltered <- Xfiltered[Xfiltered<=upperV[[1]]]
          return(mean(Xfiltered))
        }
        blot <- aggregate(CD[,exp_col], by=list(CD$cell), FUN=funAvg, lower=lowerBound, upper = upperBound)
        rownames(blot) <- blot[,1]
        blot[is.na(blot[,2]), 2] <- min_gain
        ##blot <- blot[-1]
        ##fixme consider NA
        blot[setdiff(union(tree$tip.label, tree$node.label),rownames(blot)),2] <- 0
        blot[,1] <- rownames(blot)
        colnames(blot) <- c("label", "blot")

        if(is.numeric(min_gain)){blot[blot[,2] < min_gain,2] <- min_gain}
        if(is.numeric(max_gain)){blot[blot[,2] > max_gain,2] <- max_gain
          message("adjusting max gain")}
        
        g <- ggtree(tree, ladderize=F) %<+% blot
        g <- g + xlim(0, X_Max) + geom_tree(aes(color=blot))
        
        #define color scaling and labels according to transform argument
        if(transform == "identity"){
          newMin <- min_gain
          newMax <- max_gain
          midpoint <- (newMin+newMax)/2
        }
        else if(transform == "log2"){
          newMin <- log2(min_gain)
          newMax <- log2(max_gain)
          midpoint <- 2**((newMin+newMax)/2)
        }
        else if(transform == "log10"){
          newMin <- log10(min_gain)
          newMax <- log10(max_gain)
          midpoint <- 10**((newMin+newMax)/2)
        }
        else if(transform == "log"){
          ewMin <- log(min_gain)
          newMax <- log(max_gain)
          midpoint <- exp((newMin+newMax)/2)
        }

        g <- g + scale_color_gradient2(name=exp_legend, 
                                       transform = transform,
                                       low=rgb(0,0,0.9),mid=rgb(0.9,0.85,0.8),high="red", 
                                       midpoint = midpoint, lim = c(min_gain,max_gain))
    }else{
        g <- ggtree(tree, ladderize=F)+ xlim(0, X_Max)
    }
    
    if(tip_lab){
        g <- g + xlim(0, X_Max*1.1) + geom_tiplab()     
    }
    if(node_lab){
        g <- g + geom_nodelab()     
    }
    g <- g + theme_tree2()
    return(g)
}

