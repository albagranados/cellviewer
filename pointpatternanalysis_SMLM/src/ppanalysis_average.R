library(spatstat)
library(extrafont)
loadfonts()
setwd("/home/alba/Dropbox (CRG)/postdoc_CRG/coding/cellviewer/point_pattern_analysis/src")

# define functions
build_pp <- function(pptype='marked', fn1, fn2=NULL, pn1, pn2=NULL, units='pixel', storm_file=FALSE)
{
  if (pptype %in% c("marked"))
  {
    if (storm_file){  # if original STORM file
      data1 <- read.table(fn1, skip=1, header=FALSE,
                          colClasses = c(rep("NULL", 3), rep("numeric", 2), rep("NULL", 13)))
      data2 <- read.table(fn2, skip=1, header=FALSE,
                          colClasses = c(rep("NULL", 3), rep("numeric", 2), rep("NULL", 13)))
    }
    else{
      data1 <- read.table(fn1, skip=0)
      data2 <- read.table(fn2, skip=0)
    }
    data <- rbind(data1, data2)
    m1 <- factor(rep(pn1, dim(data1)[1]))
    m2 <- factor(rep(pn2, dim(data2)[1]))
    m <- factor(c(rep(pn1, dim(data1)[1]),rep(pn2, dim(data2)[1])))
    rnahistone <- cbind(data,m)
    
    points1 <- ppp(data1[,1], data1[,2], c(min(data1[,1]),max(data1[,1])),
                   c(min(data1[,2]),max(data1[,2])),marks=m1)
    points2 <- ppp(data2[,1], data2[,2], c(min(data2[,1]),max(data2[,1])),
                   c(min(data2[,2]),max(data2[,2])),marks=m2)
    points <- ppp(data[,1], data[,2], c(min(data[,1]),max(data[,1])),
                  c(min(data[,2]),max(data[,2])),marks=m)
    points_unmarked <- ppp(data[,1], data[,2], c(min(data[,1]),max(data[,1])),
                           c(min(data[,2]),max(data[,2])))
    unitname(points) <- units # 1 pixel = 160 nm (dSTORM)
    unitname(points_unmarked) <- units
    unitname(points1) <- units
    unitname(points2) <- units
    
    return(list(first = points1, second = points2, all = points, all_unmarked = points_unmarked))
  }
  else
  {
    if (storm_file){  # if original STORM file
      data1 <- read.table(fn1, skip=1, header=FALSE,
                          colClasses = c(rep("NULL", 3), rep("numeric", 2), rep("NULL", 13)))
    }
    else{
      data1 <- read.table(fn1, skip=0)
    }
    m1 <- factor(rep(pn1, dim(data1)[1]))
    points1 <- ppp(data1[,1], data1[,2], c(min(data1[,1]),max(data1[,1])),
                   c(min(data1[,2]),max(data1[,2])),marks=m1)
    unitname(points1) <- units

    return(list(first = points1))
  }
}

# select experiment
levels1 = 'SMC1'; levels2 = 'CTCF'
units = 'pixels'; pptype='marked'
path <- "../../../../../Shared files Victoria-Alba/2017-07-17_HeLa_DualColor_SMC1_CTCF"
dirs <- list.dirs(path = path, full.names = TRUE, recursive = FALSE)

# Run twice here. Select control/treated
dir = dirs[1]; col = "black"; lty = 3  
dir = dirs[2]; col = "red"; lty = 1  
dir = dirs[3]; col = "black"; lty = 3

fileNames <- list.files(path = dir, pattern = "\\.txt$", all.files = FALSE,
                full.names = FALSE, recursive = FALSE,
                ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

r_eval = seq(0, 4, length.out=100)  # select radius range
colfunc <- colorRampPalette(c("black", "#D3D3D3")); col_palette <- colfunc(length(fileNames))
legend_plot <- list(); cell_g12 <- list(); cell_p12 <- list(); cell_r0 <- c()
pdf("plot.pdf", family="CM Roman")
plot(0, 0, ylim=c(0,6), xlim=c(0, r_eval[length(r_eval)]),type = "l", 
     xlab=paste(c('r ','[', units,']'), collapse = ''), ylab="pair correlation function",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
i <- 1; cell_no <- 1; image_no <- 1
while (i <= 2){#length(fileNames)){  # loop in number of images/files (i). Each image may have multiple cells (j)
  cat("i = ", i, '\n')
  if (is.na(pmatch("Ch", unlist(strsplit(fileNames[i], "\\_"))))){i <- i + 1; print('skip file.')} 
  else{
    split = unlist(strsplit(fileNames[i], "\\_"))
    imageName = paste(split[1:pmatch("Ch", split)-1], collapse='_')
    image_no <- image_no + 1
    
    cellNamesCh1 = fileNames[grep(paste(c(imageName, "Ch1"), collapse='_'), fileNames)] 
    cellNamesCh2 = fileNames[grep(paste(c(imageName, "Ch2"), collapse='_'), fileNames)]
    cat('cellNamesCh1 = ', cellNamesCh1, '\ncellNamesCh2 = ', cellNamesCh2, '\n');
    
    for (j in 1:max(length(cellNamesCh1),length(cellNamesCh2))){  # loop in cell number
      cat("j = ", j, '\n')
      path_to_file1 = paste(c(dir,cellNamesCh1[j]), collapse='/'); path_to_file2 = paste(c(dir,cellNamesCh2[j]), collapse='/')

      if (pptype %in% c("marked")){
        res <- build_pp(pptype=pptype, fn1=path_to_file1, fn2=path_to_file2, pn1=levels1, pn2=levels2, units=units)
        points1 = res$first; points2 = res$second; points = res$all; points_unmarked = res$all_unmarked
        area_window <- (points$window$xrange[2]-points$window$xrange[1])*(points$window$yrange[2]-points$window$yrange[1])
        cat('av. density 1 = ', points1$n/area_window, '\t av. density 2 = ', points2$n/area_window, '\n')
        # lambda1 <- intensity(points1); lambda2 <- intensity(points2)  # equivalent to
      }
      else{
        res <- build_pp(pptype=pptype, fn1=path_to_file1, pn1=levels1, units=units)
        points = res$first; 
        area_window <- (points$window$xrange[2]-points$window$xrange[1])*(points$window$yrange[2]-points$window$yrange[1])
        cat('av. density = ', points$n/area_window, '\n')
      }
      g_12 <- pcfcross(points, i=levels1, j=levels2, r=r_eval, breaks=NULL, correction="isotropic")
      # g_22 <- pcfcross(points, i=levels2, j=levels2, r=r_eval, breaks=NULL, correction="isotropic")
      # g_11 <- pcfcross(points, i=levels1, j=levels1, r=NULL, breaks=NULL, correction="isotropic")
      # g_cross_eval <- approx(g_12$r, g_12$iso, xout=r_eval, method='linear')
      # g_pair_eval <- approx(g_11$r, g_11$iso, xout=r_eval, method='linear')
      
      # M_12 <- markconnect(points, levels1, levels2)
      # M_11 <- markconnect(points, levels1, levels1)
      # M_22 <- markconnect(points, levels2, levels2)
      # p_cross_eval <- approx(M_12$r, M_12$iso, xout=r_eval, method='linear')

      # delta_12 <- lambda2*g_cross_eval$y - lambda1*g_pair_eval$y
      
      # find first zero
      # g_env <- envelope(points, pcfcross, i=levels1, j=levels2, nsim=30, global=F)
      cluster_size = r_eval[which(c(0,diff(sign(g_12$iso-1)))!=0)[1]]  # brut force
      cat('estimated cluster size r0 = ', cluster_size, '\n')
      
      # save in .RData; skip positions 1 (or more) in g_12 due to 'inf' value
      cell_g12[[cell_no]] <- data.frame(r=r_eval[10:length(r_eval)], g_cross=g_12$iso[10:length(r_eval)],
                                        name=cellNamesCh1[j])
      cell_r0[cell_no] = cluster_size
      # cell_g12[[cell_no]] <- data.frame(cbind(r=r_eval, g_cross=g_cross_eval$y))
      # cell_p12[[cell_no]] <- data.frame(cbind(r=r_eval, p_cross=p_cross_eval$y))
      
      legend_plot[cell_no] <- cellNamesCh1[j]; col <- col_palette[cell_no]  # or comment
      lines(g_12$r, g_12$iso, col=col, lty=lty)
      cell_no <- cell_no + 1
      # Sys.sleep(0.1) 
    }
    if (pptype %in% c('marked')){ i <- i + 2*max(length(cellNamesCh1),length(cellNamesCh2))}  # 2 channels
    else{i <- i + 1*max(length(cellNamesCh1),length(cellNamesCh2))}  # 1 channel (Ch1 files only)
  }
}
abline(1,0, col="black", lty=2)
legend("topright", legend=legend_plot, col=col_palette, lty=lty, cex=1.2)  # or col=col
dev.off(); embed_fonts("plot.pdf", outfile="plot.pdf")

# g_12_av <- rowMeans(sapply(cell_g12, `[[`, 'g_cross')); g_12_std <- apply(sapply(cell_g12, `[[`, 'g_cross'), 1, sd, na.rm = TRUE)

#### PLOTS
# loop 1
plot(cell_g12[[1]]$r, g_12_av, ylim=c(0.5,6), xlim=c(0, 3), type="l", col = col, xlab=paste(c('r ','[', units,']'), collapse=''), 
     ylab="av. cross-pair correlation function")
# loop 2
lines(cell_g12[[1]]$r, g_12_av, col=col, lty=1)

arrows(cell_g12[[1]]$r, g_12_av-g_12_std, r_eval, g_12_av+g_12_std, col = col, length=0.05, angle=90, code=3)
abline(h=1, col="blue")

# point pattern
pdf('pointpattern.pdf', family="CM Roman")
plot(points, pch = 1, cex=0.1, cols="black", use.marks=FALSE, main='')
dev.off(); embed_fonts("pointpattern.pdf", outfile="pointpattern.pdf")
##########

#### SAVING:
# would you like to save it as e.g., K4me2SMC1_DMSO.RData to run longitudinal_data_analysis.R?
save(cell_g12, cell_r0, r_eval, file=paste(c(paste(c(dir,paste(c(levels1,levels2), collapse = '')),
                                                   collapse='/'), '.RData'), collapse = ''))
