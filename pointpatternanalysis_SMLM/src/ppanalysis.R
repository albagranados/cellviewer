# check out (Ba10_Rspatialcourse_CMIS_PDF-Standard.pdf), (spatstatManual.pdf)
# point pattern analysis for single files. 
# create /output directory at same level of /src
# input files must be .txt. If bin->txt run in matlab bin2txt.m

### import libraries
library(spatstat)  # statistical point pattern analysis
library(lme4)  # class. statistics
library(nlme) # class. statistics
library(MASS) # class. statistics
library(car) # class. statistics
library(extrafont)  # fonts for latex pdf
loadfonts()
library(knitr) # class. statistics
# library(ggplot2)
graphics.off()

# set working directory
setwd("/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/pointpatternanalysis_SMLM/src") 

source("utilities.R")

### select data

# # dual
# pptype='marked'; units = 'pixels'; units_out = 'nm'; unit_size=160 #nm/px
# exp_name1 <- "DMSO"; exp_name2 <- "ActD"  # NA
# levels1 = 'SMC1'; levels2 = 'PolII'
# path_to_experiment1 = '/home/alba/ISIS/nfs/users/jsolon/agranados/data/vicky/2017-12-14_HeLa_DualColor_RNApolII_SMC1/RNApolII_SMC1 in HeLa DMSO Controls'
# path_to_experiment2 = "/home/alba/ISIS/nfs/users/jsolon/agranados/data/vicky/2017-12-14_HeLa_DualColor_RNApolII_SMC1/RNApolII_SMC1 in HeLa ActD treated"
# 
# storm_file=0
# exp_names <- c(exp_name1, exp_name2)
# path_to_experiments <- c(path_to_experiment1, path_to_experiment2)

# # mono
pptype='unmarked'; units = 'nm'; units_out = 'nm'; unit_size=1
levels1 = '1'; levels2 = ''
exp_name1 <- "r50"; exp_name2 <- "r100"
path_to_experiment1 = "/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/pointpatternanalysis_SMLM/output/exp1"
path_to_experiment2 = "/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/pointpatternanalysis_SMLM/output/exp2"

storm_file = 0
exp_names <- c(exp_name1, exp_name2)
path_to_experiments <- c(path_to_experiment1, path_to_experiment2)

### compute
compute_ppsummaries <- data.frame(run=1, 
                                  nearestneighbour=0,
                                  Kfunction=1,
                                    Lfunction=1,
                                  crosscorrelation=1, 
                                  markconnection=0,
                                  envelopes=0,
                                  plot_functions=0)
longitudinal_dataset <- data.frame(run=1,
                                    generate = 1,
                                      stat_rrange_probs = 0.1,
                                      units = units_out,
                                      plot_2experiments = 1,
                                      save=1,
                                    statanalysis = 1,
                                      fitting = 'exp'
                                   )

### parameters
# r_eval = seq(0, 2, length.out=200) # 300
r_eval = seq(0, 300, length.out=200) 

# -------------------- Dual color: MARKED POINT PATTERN ------
# ------------------------------------------------------------
if (pptype == "marked"){

  if (compute_ppsummaries$run){
    for (ii in 1:length(path_to_experiments)){
    
      exp_name <- exp_names[ii]
      path_to_experiment <- path_to_experiments[ii]
      dirs_channels <- list.dirs(path = path_to_experiment, full.names = TRUE, recursive = FALSE)  # Ch1&Ch2 of experiment
      
      cat('\n********* Point pattern analysis of experiment ', exp_name, '************\n')
    
      # read files in each channel
      fileNamesCh1 <- list.files(path = dirs_channels[1], pattern = "\\.txt$", all.files = FALSE,
                                 full.names = FALSE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
      fileNamesCh2 <- list.files(path = dirs_channels[2], pattern = "\\.txt$", all.files = FALSE,
                                 full.names = FALSE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
      num_files = 1:length(fileNamesCh1)
      
      # initialize
      g12_r0_all <- c(); g22_r0_all <- c(); g11_r0_all <- c(); g_r0_all <- c()
      
      G12_all <- list(); K12_all <- list(); L12_all <- list();
      g12_all <- list(); p12_all <- list(); M12_all <- list();
      G22_all <- list(); K22_all <- list(); L22_all <- list();
      g22_all <- list(); p22_all <- list(); M22_all <- list();
      g11_all <- list();
      
      intensities <- matrix(nrow = length(fileNamesCh1), ncol = 2)
      
      cat('level 1 = ', levels1, "\nlevel 2 = ", levels2, '\n')
      
      # if multiple crops of a cell, this variable should contain the reference to the cell.
      cellName_previous <- paste(unlist(strsplit(fileNamesCh1[1], '_'))[1:(length(unlist(strsplit(fileNamesCh1[1], '_')))-2)], 
                                collapse='_')
      cell_nos <- c(); no <- 1
      for (j in num_files){
        
        cellName_current <- paste(unlist(strsplit(fileNamesCh1[j], '_'))[1:(length(unlist(strsplit(fileNamesCh1[1], '_')))-2)], 
                                   collapse='_')
        ### read & build point pattern
        cat("Analyzing file ", j, ' [name =', fileNamesCh1[j], ']... \n')
        path_to_file1 = paste(c(dirs_channels[1],fileNamesCh1[j]), collapse='/')
        path_to_file2 = paste(c(dirs_channels[2],fileNamesCh2[j]), collapse='/')
          
        res <- build_pp(pptype=pptype, fn1=path_to_file1, fn2=path_to_file2, pn1=levels1, pn2=levels2, storm_file=storm_file)
        points1 = res$first; points2 = res$second; points = res$all; points_unmarked = res$all_unmarked
      
        summary(points)
      
        ### plot pp
        if(j==num_files[1]){  # plot only first point pattern
          openpdf(paste("pointpattern_", exp_name, ".pdf", sep = ''))
          plot(points, cex=0.2, main='', pch=16, cols=c(rgb(173,216,230,max = 255), rgb(255,165,0,max = 255,alpha=125)), 
               use.marks=TRUE, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
          closepdf(paste("pointpattern_", exp_name, ".pdf", sep = ''))
        }
        
        if(compute_ppsummaries$nearestneighbour){
          cat('Computing nearest-neighbour distance d.f. ... ')
          G_11 <- Gcross(points, i=levels1, j=levels1, r=NULL, breaks=NULL, correction="rs")
          G_22 <- Gcross(points, i=levels2, j=levels2, r=NULL, breaks=NULL, correction="rs")
          G_12 <- Gcross(points, i=levels1, j=levels2, r=NULL, breaks=NULL, correction="rs")
          if(compute_ppsummaries$envelopes){ 
            G_env <- envelope(points, Gcross, i=levels1, j=levels2, r=NULL, correction="rs", nsim=99, global=F, savepatterns=TRUE)
          }
          G12_all[[j]] <- data.frame(r=G_12$r, G12=G_12$iso)
          G22_all[[j]] <- data.frame(r=G_22$r, G22=G_22$iso)
          cat('Done.\n')
      
          if (compute_ppsummaries$plot_functions){
            openpdf("nearestneighbour.pdf")
            plot(G_11$r, G_11$rs, col="black", xlab=paste('r ','[', units,']', sep = ''), 
                 ylab=expression(G(r)), type='l', lty=1)
            lines(G_22$r, G_22$rs, col="black", lty=6)
            lines(G_12$r, G_12$rs, col="black", lty=3)
            lines(G_12$r, G_12$theo, col="grey50", lty=2)
            if(compute_ppsummaries$envelopes){
              polygon(c(0, G_env$r, rev(G_env$r), 0), c(c(0, G_env$lo), rev(c(0,G_env$hi))), col=rgb(0.5, 0.5, 0.5,0.5), lty=0)
            }
            legend("bottomright", legend=c(expression(G[11]), expression(G[22]), expression(G[12]),
                                           expression(G[pois])), 
                   lty=c(1,6,3,2), col=c("black", "black","black","grey50"))
            closepdf("nearestneighbour.pdf")
          }
        }
      
        if(compute_ppsummaries$Kfunction){
          cat('Computing K-function... ')
          K_12 <- Kcross(points, i=levels1, j=levels2, r=r_eval, breaks=NULL, correction="border", ratio=FALSE)
          # K_11 <- Kcross(points, i=levels1, j=levels1, r=r_eval, breaks=NULL, correction="border", ratio=FALSE) 
          K_22 <- Kcross(points, i=levels2, j=levels2, r=r_eval, breaks=NULL, correction="border", ratio=FALSE)
          if(compute_ppsummaries$envelopes){ 
              K_env <- envelope(points, Kcross, i=levels1, j=levels2, r=r_eval, nsim=20, global=F)
          }
          K12_all[[j]] <- data.frame(r=r_eval, K12=K_12$border)
          K22_all[[j]] <- data.frame(r=r_eval, K22=K_22$border)
          cat('Done.\n')
          
          if (compute_ppsummaries$plot_functions){
            openpdf("Kfunction.pdf")
            plot(r_eval, K_11$border, col="black", xlab=paste('r ','[', units,']', sep = ''), 
                 ylab=expression(K(r)), type='l', lty=1)
            lines(r_eval, K_22$border, col="black", lty=6)
            lines(r_eval, K_12$border, col="black", lty=3)
            lines(r_eval, K_12$theo, col="grey50", lty=2)
            if(compute_ppsummaries$envelopes){ 
                polygon(c(0, K_env$r, rev(K_env$r), 0), c(c(0, K_env$lo), rev(c(0,K_env$hi))), col=rgb(0.5, 0.5, 0.5,0.5), lty=0)
            }
            legend("bottomright", legend=c(expression(K[11]), expression(K[22]), expression(K[12]),expression(K[pois])), 
                   lty=c(1,6,3,2), col=c("black", "black","black","grey50"))
            closepdf('Kfunction.pdf')
          }
        }
      
        if(compute_ppsummaries$Kfunction && compute_ppsummaries$Lfunction){
          cat('Computing L-function... ')
          if(compute_ppsummaries$envelopes){ 
            L_env <- envelope(points, Lcross, i=levels1, j=levels2, correction="isotropic", nsim=20, global=F)
          }
          L12_all[[j]] <- data.frame(r=r_eval, L12=sqrt(K_12$border/pi)-r_eval)
          L22_all[[j]] <- data.frame(r=r_eval, L22=sqrt(K_22$border/pi)-r_eval)
          cat('Done.\n')
          
          if (compute_ppsummaries$plot_functions){
            openpdf("Lfunction.pdf")
            plot(r_eval, sqrt(K_11$border/pi)-r_eval, col="black", xlab=paste('r ','[', units,']', sep = ''), 
                 ylab=expression(L(r)-r), type='l', lty=1) 
            lines(r_eval, sqrt(K_22$border/pi)-r_eval, col="black", lty=6)
            lines(r_eval, sqrt(K_12$border/pi)-r_eval, col="black", lty=3)
            lines(r_eval, sqrt(K_12$theo/pi)-r_eval, col="grey50", lty=2)
            if(compute_ppsummaries$envelopes){ 
                polygon(c(0, K_env$r, rev(K_env$r), 0), c(c(0, sqrt(K_env$lo/pi)- K_env$r), rev(c(0,sqrt(K_env$hi/pi)-K_env$r))),
                    col=rgb(0.5, 0.5, 0.5,0.5), lty=0)
            }
            legend("topright", legend=c(expression(L[11]-r), expression(L[22]-r), expression(L[12]-r),
                                        expression(L[pois]-r)), lty=c(1,6,3,2), 
                   col=c("black", "black","black","grey50"))
            closepdf('Lfunction.pdf')
          }
        }
        
        if(compute_ppsummaries$crosscorrelation){
          cat('Computing pair-correlation function... ')
          g_11 <- pcfcross(points, i=levels1, j=levels1, r=r_eval, breaks=NULL, correction="isotropic")
          g_22 <- pcfcross(points, i=levels2, j=levels2, r=r_eval, breaks=NULL, correction="isotropic")
          g_12 <- pcfcross(points, i=levels1, j=levels2, r=r_eval, breaks=NULL, correction="isotropic")
          if(compute_ppsummaries$envelopes){ 
            g_env <- envelope(points, pcfcross, i=levels1, j=levels2, nsim=5, global=F)
          }
          cluster_size_12 = r_eval[which(c(0,diff(sign(g_12$iso-1)))!=0)[1]]  # brut force
          cluster_size_22 = r_eval[which(c(0,diff(sign(g_22$iso-1)))!=0)[1]]  # brut force
          cluster_size_11 = r_eval[which(c(0,diff(sign(g_11$iso-1)))!=0)[1]]  # brut force
          g12_all[[j]] <- data.frame(r=r_eval, g12=g_12$iso)
          g22_all[[j]] <- data.frame(r=r_eval, g22=g_22$iso)
          g11_all[[j]] <- data.frame(r=r_eval, g11=g_11$iso)
          g12_r0_all[j] = cluster_size_12
          g22_r0_all[j] = cluster_size_22
          g11_r0_all[j] = cluster_size_11
          cat('Done.\n')
          
          if (compute_ppsummaries$plot_functions){
            openpdf("crosscorrelation.pdf")
            plot(r_eval, g_12$iso, col="black", xlab=paste('r ','[', units,']', sep = ''), 
                 ylab=expression(g(r)), type='l', lty=3) 
            # lines(r_eval, g_22$iso, col="black", lty=6)
            # lines(r_eval, g_11$iso, col="black", lty=1)
            lines(r_eval, g_12$theo, col="grey50", lty=2)
            if(compute_ppsummaries$envelopes){ 
              polygon(c(0, g_env$r, rev(g_env$r), 0), c(c(0, g_env$lo), rev(c(0,g_env$hi))), col=rgb(0.5, 0.5, 0.5,0.5), lty=0)
            }
            legend("topright", legend=c(expression(g[12]), expression(g[11]), expression(g[22]),
                                        expression(g[pois])), lty=c(3,1,6,2), 
                   col=c("black", "black","black","grey50"))
            closepdf("crosscorrelation.pdf")
          }
        }
        
        if(compute_ppsummaries$markconnection){
          cat('Computing mark connection function... ')
          M_11 <- markconnect(points, i=levels1, j=levels1, r=r_eval, breaks=NULL, correction="isotropic")
          M_22 <- markconnect(points, i=levels2, j=levels2, r=r_eval, breaks=NULL, correction="isotropic")
          M_12 <- markconnect(points, i=levels1, j=levels2, r=r_eval, breaks=NULL, correction="isotropic")
          M12_all[[j]] <- data.frame(r=r_eval, M12=M_12$iso)
          M22_all[[j]] <- data.frame(r=r_eval, M22=M_12$iso)
          cat('Done.\n')
          
          if (compute_ppsummaries$plot_functions){
            openpdf("markconnection.pdf")
            plot(r_eval, M_11$iso, ylim=c(0,0.5), col="black", xlab=paste('r ','[', units,']', sep = ''), 
                 ylab=expression(p(r)), type='l', lty=1) 
            lines(r_eval, M_22$iso, col="black", lty=6)
            lines(r_eval, M_12$iso, col="black", lty=3)
            lines(r_eval, M_12$theo, col="grey50", lty=2)
            if(compute_ppsummaries$envelopes){ 
              polygon(c(0, M_env$r, rev(M_env$r), 0), c(c(0, M_env$lo), rev(c(0,M_env$hi))), col=rgb(0.5, 0.5, 0.5,0.5), lty=0)
            }
            legend("topright", legend=c(expression(p[11]), expression(p[22]), expression(p[12]),
                                        expression(p[pois])), lty=c(1,6,3,2), 
                   col=c("black", "black","black","grey50"))
            closepdf("markconnection.pdf")
          }
        }
        
        intensities[j,] <- intensity(points)
        
        if (cellName_current != cellName_previous){ no = no + 1 }
        cell_nos <- append(cell_nos, no); cellName_previous <- cellName_current
      }
      # save txt file with name of analyzed file
      textfile <- file(paste(c(paste(c(path_to_experiment,paste(c(levels1,levels2,'_',exp_name, '_Ch1_fileNames'), collapse = '')),
                                     collapse='/'), '.txt'), collapse = ''), "w")
      cat(fileNamesCh1, file=textfile, sep='\n')
      close(textfile)
      textfile <- file(paste(c(paste(c(path_to_experiment,paste(c(levels1,levels2,'_',exp_name, '_Ch2_fileNames'), collapse = '')),
                                     collapse='/'), '.txt'), collapse = ''), "w")
      cat(fileNamesCh2, file=textfile, sep='\n')
      close(textfile)
      
      # save:
      path_to_RData <- paste(c(paste(c(path_to_experiment,paste(c(levels1,levels2,'_',exp_name), collapse = '')),
                                     collapse='/'), '.RData'), collapse = '')
      save(path_to_experiments, exp_name, levels1, levels2, 
           G12_all, G22_all, K12_all, K22_all, L12_all, L22_all, g12_all, g22_all, g11_all, 
           g12_r0_all, g22_r0_all, g11_r0_all,  
           r_eval, cell_nos, intensities,
           file=path_to_RData)
    }
  }

  if(longitudinal_dataset$run){
    if (longitudinal_dataset$generate){
      generate_longitudinal_data(path_to_experiments, pptype, exp_names, levels1, levels2, 
                                 units=longitudinal_dataset$units, unit_size = unit_size,
                                 plot_2experiments=longitudinal_dataset$plot_2experiments,
                                 save=longitudinal_dataset$save, stat_rrange_probs=longitudinal_dataset$stat_rrange_probs)
    }
    if (longitudinal_dataset$statanalysis){
      statistical_analysis(path_to_experiments, pptype, fitting='exp')
    }
  }
}


# ----------------------- UNMARKED POINT PATTERN ------
# -----------------------------------------------------
if (pptype %in% c("unmarked")){
  if (compute_ppsummaries$run){
    for (ii in 1:length(path_to_experiments)){
      
      exp_name <- exp_names[ii]
      path_to_experiment <- path_to_experiments[ii]
      dirs_channels <- list.dirs(path = path_to_experiment, full.names = TRUE, recursive = FALSE)  # Ch1&Ch2 of experiment
      
      cat('\n********* Point pattern analysis of experiment ', exp_name, '************\n')

      # read files in each channel
      fileNamesCh1 <- list.files(path = dirs_channels[1], pattern = "\\.txt$", all.files = FALSE,
                                 full.names = FALSE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
      num_files = 1:length(fileNamesCh1)
      
      # initialize
      g11_r0_all <- c()
      G11_all <- list(); K11_all <- list(); L11_all <- list();
      g11_all <- list(); p11_all <- list(); M11_all <- list();
      
      intensities <- matrix(nrow = length(fileNamesCh1), ncol = 1)
      
      cat('level 1 = ', levels1, "\nlevel 2 = ", levels2, '\n')
      
      # if multiple crops of a cell, this variable should contain the reference to the cell.
      cellName_previous <- paste(unlist(strsplit(fileNamesCh1[1], '_'))[1:(length(unlist(strsplit(fileNamesCh1[1], '_')))-2)], 
                                 collapse='_')
      cell_nos <- c(); no <- 1
      for (j in num_files){
        
        cellName_current <- paste(unlist(strsplit(fileNamesCh1[j], '_'))[1:(length(unlist(strsplit(fileNamesCh1[1], '_')))-2)], 
                                  collapse='_')
        
        ### read & build point pattern
        cat("Analyzing file ", j, ' [name =', fileNamesCh1[j], ']... \n')
        path_to_file1 = paste(c(dirs_channels[1],fileNamesCh1[j]), collapse='/')
        
        res <- build_pp(pptype=pptype, fn1=path_to_file1, pn1=levels1, units=units, storm_file=storm_file)
        points = res$first

        ### plot pp
        if(j==num_files[1]){  # plot only first point pattern
          openpdf(paste("pointpattern_", exp_name, ".pdf", sep = ''))
          plot(points, cex=0.2, main='', pch=16, cols="black", 
               use.marks=TRUE, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
          closepdf(paste("pointpattern_", exp_name, ".pdf", sep = ''))
        }
        
        if(compute_ppsummaries$nearestneighbour){
          cat('Computing nearest-neighbour distance d.f. ... ')
          G_11 <- Gcross(points, i=levels1, j=levels1, r=NULL, breaks=NULL, correction="rs")
          if(compute_ppsummaries$envelopes){ 
            G_env <- envelope(points, Gcross, i=levels1, j=levels1, r=NULL, correction="rs", nsim=99, global=F, savepatterns=TRUE)
          }
          G11_all[[j]] <- data.frame(r=G_11$r, G11=G_11$rs)
          cat('Done.\n')
          
          if (compute_ppsummaries$plot_functions){
            openpdf("nearestneighbour.pdf")
            plot(G_11$r, G_11$rs, col="black", xlab=paste('r ','[', units,']', sep = ''), 
                 ylab=expression(G(r)), type='l', lty=1)
            lines(G_11$r, G_11$theo, col="grey50", lty=2)
            if(compute_ppsummaries$envelopes){ 
              polygon(c(0, G_env$r, rev(G_env$r), 0), c(c(0, G_env$lo), rev(c(0,G_env$hi))), col=rgb(0.5, 0.5, 0.5,0.5), lty=0)
            }
            legend("bottomright", legend=c(expression(G[11]), expression(G[pois])), 
                   lty=c(1,2), col=c("black", "grey50"))
            closepdf("nearestneighbour.pdf")
          }
        }
        
        if(compute_ppsummaries$Kfunction){
          cat('Computing K-function... ')
          K_11 <- Kcross(points, i=levels1, j=levels1, r=r_eval, breaks=NULL, correction="border", ratio=FALSE)
          if(compute_ppsummaries$envelopes){ 
            K_env <- envelope(points, Kcross, i=levels1, j=levels1, r=r_eval, nsim=20, global=F)
          }
          K11_all[[j]] <- data.frame(r=r_eval, K11=K_11$border)
          cat('Done.\n')
          
          if (compute_ppsummaries$plot_functions){
            openpdf("Kfunction.pdf")
            plot(r_eval, K_11$border, col="black", xlab=paste('r ','[', units,']', sep = ''), 
                 ylab=expression(K(r)), type='l', lty=1) 
            lines(r_eval, K_11$theo, col="grey50", lty=2)
            if(compute_ppsummaries$envelopes){ 
              polygon(c(0, K_env$r, rev(K_env$r), 0), c(c(0, K_env$lo), rev(c(0,K_env$hi))), col=rgb(0.5, 0.5, 0.5,0.5), lty=0)
            }
            legend("bottomright", legend=c(expression(K[11]),expression(K[pois])), 
                   lty=c(1,2), col=c("black", "grey50"))
            closepdf('Kfunction.pdf')
          }
        }
        
        if(compute_ppsummaries$Kfunction && compute_ppsummaries$Lfunction){
          cat('Computing L-function... ')
          if(compute_ppsummaries$envelopes){ 
            L_env <- envelope(points, Lcross, i=levels1, j=levels2, correction="isotropic", nsim=20, global=F)
          }
          cat('Done.\n')
          
          if (compute_ppsummaries$plot_functions){
            openpdf("Lfunction.pdf")
            plot(r_eval, sqrt(K_11$border/pi)-r_eval, col="black", xlab=paste('r ','[', units,']', sep = ''), 
                 ylab=expression(L(r)-r), type='l', lty=1) 
            lines(r_eval, sqrt(K_11$theo/pi)-r_eval, col="grey50", lty=2)
            if(compute_ppsummaries$envelopes){ 
              polygon(c(0, K_env$r, rev(K_env$r), 0), c(c(0, sqrt(K_env$lo/pi)- K_env$r), rev(c(0,sqrt(K_env$hi/pi)-K_env$r))),
                      col=rgb(0.5, 0.5, 0.5,0.5), lty=0)
            }
            legend("topright", legend=c(expression(L[11]-r), expression(L[pois]-r)), lty=c(1,2), 
                   col=c("black", "grey50"))
            closepdf('Lfunction.pdf')
          }
        }
        
        if(compute_ppsummaries$crosscorrelation){
          cat('Computing pair-correlation function... ')
          # g_11 <- pcfcross(points, i=levels1, j=levels1, r=r_eval, breaks=NULL, correction="isotropic")
          g_11 <- pcf(points, r=r_eval, correction="isotropic")
          if(compute_ppsummaries$envelopes){ 
            g_env <- envelope(points, pcfcross, i=levels1, j=levels1, nsim=5, global=F)
          }
          cluster_size = r_eval[which(c(0,diff(sign(g_11$iso-1)))!=0)[1]]  # brut force
          g11_all[[j]] <- data.frame(r=r_eval, g11=g_11$iso)
          g11_r0_all[j] = cluster_size
          cat('Done.\n')
          
          if (compute_ppsummaries$plot_functions){
            openpdf("crosscorrelation.pdf")
            plot(r_eval, g_11$iso, col="black", xlab=paste('r ','[', units,']', sep = ''), 
                 ylab=expression(g(r)), type='l', lty=1) 
            lines(r_eval, g_11$theo, col="grey50", lty=2)
            if(compute_ppsummaries$envelopes){ 
              polygon(c(0, g_env$r, rev(g_env$r), 0), c(c(0, g_env$lo), rev(c(0,g_env$hi))), col=rgb(0.5, 0.5, 0.5,0.5), lty=0)
            }
            legend("topright", legend=c(expression(g[11]), expression(g[pois])), lty=c(1,2), 
                   col=c("black", "grey50"))
            closepdf("crosscorrelation.pdf")
          }
        }
        
        if(compute_ppsummaries$markconnection){
          cat('Computing mark connection function... ')
          M_11 <- markconnect(points, i=levels1, j=levels1, r=r_eval, breaks=NULL, correction="isotropic")
          M11_all[[j]] <- data.frame(r=r_eval, M11=M_11$iso)
          cat('Done.\n')
          
          if (compute_ppsummaries$plot_functions){
            openpdf("markconnection.pdf")
            plot(r_eval, M_11$iso, ylim=c(0,0.5), col="black", xlab=paste('r ','[', units,']', sep = ''), 
                 ylab=expression(p(r)), type='l', lty=1) 
            lines(r_eval, M_11$theo, col="grey50", lty=2)
            if(compute_ppsummaries$envelopes){ 
              polygon(c(0, M_env$r, rev(M_env$r), 0), c(c(0, M_env$lo), rev(c(0,M_env$hi))), col=rgb(0.5, 0.5, 0.5,0.5), lty=0)
            }
            legend("topright", legend=c(expression(p[11]), expression(p[pois])), lty=c(1,2), 
                   col=c("black", "grey50"))
            closepdf("markconnection.pdf")
          }
        }
        intensities[j,] <- intensity(points)
        
        if (cellName_current != cellName_previous){ no = no + 1 }
        cell_nos <- append(cell_nos, no); cellName_previous <- cellName_current
      }
      # save txt file with name of analyzed file
      textfile <- file(paste(c(paste(c(path_to_experiment,paste(c(levels1,levels2,'_',exp_name, '_Ch1_fileNames'), collapse = '')),
                                     collapse='/'), '.txt'), collapse = ''), "w")
      cat(fileNamesCh1, file=textfile, sep='\n')
      close(textfile)

      # save:
      path_to_RData <- paste(c(paste(c(path_to_experiment,paste(c(levels1,levels2,'_',exp_name), collapse = '')),
                                     collapse='/'), '.RData'), collapse = '')
      save(path_to_experiments, exp_name, levels1, levels2, 
           G11_all, K11_all, L11_all, g11_all, p11_all, g11_r0_all, r_eval, intensities,
           file=path_to_RData)
    }
  }
  
  if(longitudinal_dataset$run){
    if (longitudinal_dataset$generate){
      generate_longitudinal_data(path_to_experiments, pptype, exp_names, levels1, levels2, 
                                 units=longitudinal_dataset$units, unit_size = unit_size,
                                 plot_2experiments=longitudinal_dataset$plot_2experiments,
                                 save=longitudinal_dataset$save, stat_rrange_probs=longitudinal_dataset$stat_rrange_probs)
    }
    if (longitudinal_dataset$statanalysis){
      statistical_analysis(path_to_experiments, pptype, fitting=longitudinal_dataset$fitting)
    }
  }
}
