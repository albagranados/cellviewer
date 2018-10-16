### define functions
build_pp <- function(pptype='marked', fn1, fn2, pn1, pn2, units='pixel', storm_file=TRUE)
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
                   c(min(data1[,2]),max(data1[,2])), marks=m1)
    unitname(points1) <- units
    
    return(list(first = points1))
  }
}

openpdf <- function(file, dir="../output/", width=7, height=7, pointsize=18)
{
  dir.create(file.path(dir), showWarnings = FALSE)
  fname <- paste(dir,file, sep="")
  pdf(fname, family="CM Roman", width = width, height = height, pointsize = pointsize)
}
closepdf <- function(file, dir="../output/")
{
  fname <- paste(dir,file, sep="")
  dev.off(); embed_fonts(fname, outfile=fname)
  plot_crop(fname, quiet = TRUE)
}

generate_longitudinal_data <- function(path_to_experiments, pptype, exp_names, levels1, levels2, 
                                       units = 'pixels', unit_size = 1, reso_r = 20, 
                                       plot_2experiments=TRUE, save=TRUE, stat_rrange_probs=0.5)
{
  cat('Generating longitudinal data...\n')
  exp_name1 <- exp_names[1]; exp_name2 <- exp_names[2]
  path_to_experiment1 <- path_to_experiments[1];   path_to_experiment2 <- path_to_experiments[2]
  filename_experiments <- "_2experiments"
  if (length(path_to_experiments) < 2){  # onlly one experiment
    exp_name1 <- exp_names[1]; exp_name2 <- NULL
    path_to_experiment1 <- path_to_experiments[1]; path_to_experiment2 <- NULL
    filename_experiments <- "_1experiment"
  }
  path_to_RData1 <- paste(c(path_to_experiment1,paste(c(levels1,levels2,'_',exp_name1, '.RData'), collapse = '')), collapse='/')
  path_to_RData2 <- paste(c(path_to_experiment2,paste(c(levels1,levels2,'_',exp_name2, '.RData'), collapse = '')), collapse='/')
  
  if (pptype == "marked"){
    load(path_to_RData1)
    
    r_eval <- r_eval*unit_size
    data.g12.exp1 <- t(sapply(g12_all, `[[`, 'g12'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]
    data.g22.exp1 <- t(sapply(g22_all, `[[`, 'g22'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]
    data.g11.exp1 <- t(sapply(g11_all, `[[`, 'g11'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]
    data.K12.exp1 <- t(sapply(K12_all, `[[`, 'K12'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]*(unit_size^2)
    data.K22.exp1 <- t(sapply(K22_all, `[[`, 'K22'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]*(unit_size^2)
    g12_r0_exp1 <- unit_size*g12_r0_all
    g12_rend_exp1 <- which(r_eval[which(r_eval>=reso_r)]>=quantile(g12_r0_exp1, probs=stat_rrange_probs, na.rm=TRUE))[1]
    g22_r0_exp1 <- unit_size*g22_r0_all
    g22_rend_exp1 <- which(r_eval[which(r_eval>=reso_r)]>=quantile(g22_r0_exp1, probs=stat_rrange_probs, na.rm=TRUE))[1]
    g11_r0_exp1 <- unit_size*g11_r0_all
    g11_rend_exp1 <- which(r_eval[which(r_eval>=reso_r)]>=quantile(g11_r0_exp1, probs=stat_rrange_probs, na.rm=TRUE))[1]
    intensities.exp1 <- intensities/(unit_size^2)
    npoints.exp1 <- npoints
    
    cat('for statistical analysis (exp1), fill g12 with NA function values above r = ', r_eval[which(r_eval>=reso_r)][g12_rend_exp1], '\n')
    cat('for statistical analysis (exp1), fill g22 with NA function values above r = ', r_eval[which(r_eval>=reso_r)][g22_rend_exp1], '\n')
    cat('for statistical analysis (exp1), fill g11 with NA function values above r = ', r_eval[which(r_eval>=reso_r)][g11_rend_exp1], '\n')
    
    load(path_to_RData2)
    r_eval <- r_eval*unit_size
    data.g12.exp2 <- t(sapply(g12_all, `[[`, 'g12'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]
    data.g22.exp2 <- t(sapply(g22_all, `[[`, 'g22'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]
    data.g11.exp2 <- t(sapply(g11_all, `[[`, 'g11'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]
    data.K12.exp2 <- t(sapply(K12_all, `[[`, 'K12'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]*(unit_size^2)
    data.K22.exp2 <- t(sapply(K22_all, `[[`, 'K22'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]*(unit_size^2)
    g12_r0_exp2 <- unit_size*g12_r0_all
    g12_rend_exp2 <- which(r_eval[which(r_eval>=reso_r)]>=quantile(g12_r0_exp2, probs=stat_rrange_probs, na.rm=TRUE))[1]
    g22_r0_exp2 <- unit_size*g22_r0_all
    g22_rend_exp2 <- which(r_eval[which(r_eval>=reso_r)]>=quantile(g22_r0_exp2, probs=stat_rrange_probs, na.rm=TRUE))[1]
    g11_r0_exp2 <- unit_size*g11_r0_all
    g11_rend_exp2 <- which(r_eval[which(r_eval>=reso_r)]>=quantile(g11_r0_exp2, probs=stat_rrange_probs, na.rm=TRUE))[1]
    intensities.exp2 <- intensities/(unit_size^2)
    npoints.exp2 <- npoints
    
    cat('for statistical analysis (exp2), fill g12 with NA function values above r = ', r_eval[which(r_eval>=reso_r)][g12_rend_exp2], '\n')
    cat('for statistical analysis (exp2), fill g22 with NA function values above r = ', r_eval[which(r_eval>=reso_r)][g22_rend_exp2], '\n')
    cat('for statistical analysis (exp2), fill g11 with NA function values above r = ', r_eval[which(r_eval>=reso_r)][g11_rend_exp2], '\n')
    
    data_temp <- data.frame( g12 = (rbind(data.g12.exp1, data.g12.exp2)), g22 = (rbind(data.g22.exp1, data.g22.exp2)),
                             g11 = (rbind(data.g11.exp1, data.g11.exp2)),
                             K12 = (rbind(data.K12.exp1, data.K12.exp2)), K22 = (rbind(data.K22.exp1, data.K22.exp2)))
    data_wide <- data.frame(id = factor(1:nrow(data_temp)),
                            experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(data.g12.exp1), nrow(data.g12.exp2)))),
                            data_temp)
    data_wide$experiment <- relevel(data_wide$experiment, ref=exp_name1)
    
    parameters_g12 <- data.frame(id = factor(1:nrow(data_temp)),
                                 experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(data.g12.exp1),
                                                                                    nrow(data.g12.exp2)))), 
                                 r0 = c(g12_r0_exp1, g12_r0_exp2))
    
    parameters_g22 <- data.frame(id = factor(1:nrow(data_temp)),
                                 experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(data.g22.exp1),
                                                                                    nrow(data.g22.exp2)))), 
                                 r0 = c(g22_r0_exp1, g22_r0_exp2))
    
    parameters_g11 <- data.frame(id = factor(1:nrow(data_temp)),
                                 experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(data.g11.exp1),
                                                                                    nrow(data.g11.exp2)))), 
                                 r0 = c(g11_r0_exp1, g11_r0_exp2))
    
    density <- data.frame(id = factor(1:nrow(data_temp)),
                          experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(intensities.exp1),
                                                                             nrow(intensities.exp2)))), 
                          level1 = c(intensities.exp1[,1], intensities.exp2[,1]), 
                          level2 = c(intensities.exp1[,2], intensities.exp2[,2]))
    npoints <- data.frame(id = factor(1:nrow(data_temp)),
                          experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(npoints.exp1),
                                                                             nrow(npoints.exp2)))),
                          level1 = c(npoints.exp1[,1], npoints.exp2[,1]),
                          level2 = c(npoints.exp1[,2], npoints.exp2[,2]))
    # npoints <- 'temp'
    
    # convert wide-formatted data into long
    num_measurements <- dim(data.g12.exp1)[2]
    data_allrange <- reshape(data_wide, 
                             varying=list(names(data_wide)[3:(3+num_measurements-1)],
                                          names(data_wide)[(3+num_measurements):(3+2*num_measurements-1)], 
                                          names(data_wide)[(3+2*num_measurements):(3+3*num_measurements-1)],
                                          names(data_wide)[(3+3*num_measurements):(3+4*num_measurements-1)],
                                          names(data_wide)[(3+4*num_measurements):(3+5*num_measurements-1)]), 
                             idvar=c("id", "experiment"),
                             direction="long", timevar='r', 
                             times=r_eval[(length(r_eval)-(dim(data.g12.exp1)[2])+1):length(r_eval)], 
                             v.names=c('g12', 'g22', 'g11', 'K12', 'K22'))
    data_allrange$experiment <- relevel(data_allrange$experiment,ref=exp_name1)
    
    # data_wide[which(data_wide$experiment==exp_name1), (g12_rend_exp1+2):length(data_wide[1,])] <- NA
    # data_wide[which(data_wide$experiment==exp_name2), (g12_rend_exp2+2):length(data_wide[1,])] <- NA
    data_wide[which(data_wide$experiment==exp_name1), (2+g12_rend_exp1):(2+num_measurements)] <- NA
    data_wide[which(data_wide$experiment==exp_name1), (2+num_measurements+g22_rend_exp1):(2+2*num_measurements)] <- NA
    data_wide[which(data_wide$experiment==exp_name1), (2+2*num_measurements+g11_rend_exp1):(2+3*num_measurements)] <- NA
    data_wide[which(data_wide$experiment==exp_name2), (2+g12_rend_exp2):(2+num_measurements)] <- NA
    data_wide[which(data_wide$experiment==exp_name2), (2+num_measurements+g22_rend_exp2):(2+2*num_measurements)] <- NA
    data_wide[which(data_wide$experiment==exp_name2), (2+2*num_measurements+g11_rend_exp2):(2+3*num_measurements)] <- NA
    data <- reshape(data_wide, 
                    varying=list(names(data_wide)[3:(3+num_measurements-1)],
                                 names(data_wide)[(3+num_measurements):(3+2*num_measurements-1)], 
                                 names(data_wide)[(3+2*num_measurements):(3+3*num_measurements-1)],
                                 names(data_wide)[(3+3*num_measurements):(3+4*num_measurements-1)],
                                 names(data_wide)[(3+4*num_measurements):(3+5*num_measurements-1)]), 
                    idvar=c("id", "experiment"),
                    direction="long", timevar='r', 
                    times=r_eval[(length(r_eval)-(dim(data.g12.exp1)[2])+1):length(r_eval)], 
                    v.names=c('g12', 'g22', 'g11', 'K12', 'K22'))
    data$experiment <- relevel(data_allrange$experiment,ref=exp_name1)
    
    col_palette <-  c(rgb(173,216,230,max = 255), rgb(255,165,0,max = 255,alpha=125)); lty=c(1,2)
    r_eval_plot <- which(r_eval>=reso_r)[1]:(0.8*length(r_eval))
    ylim=c(0.5,10); xlim <- c(r_eval[r_eval_plot][1],r_eval[r_eval_plot][length(r_eval[r_eval_plot])])
    if (plot_2experiments & (length(g12_all)>0)){
      openpdf("g12_2experiments.pdf")
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''), 
           ylab=eval(bquote(expression(g[.(levels1)][','][.(levels2)](r)))), type='n', ylim = ylim)
      for (ii in 1:length(exp_names)){
        subset <- subset(data_allrange, data$experiment %in% exp_names[ii])
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, g12, lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
      }
      abline(h=1, lty=2, col='black'); legend("topright", legend=exp_names, lty=lty, col=col_palette)
      closepdf("g12_2experiments.pdf");  cat('\t\'g12_2experiments.pdf\' created.\n')
    }
    if (plot_2experiments & (length(g22_all)>0)){
      openpdf("g22_2experiments.pdf")
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''), 
           ylab=eval(bquote(expression(g[.(levels2)][','][.(levels2)](r)))), type='n', ylim = ylim)
      for (ii in 1:length(exp_names)){
        subset <- subset(data_allrange, data$experiment %in% exp_names[ii])
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, g22, lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
      }
      abline(h=1, lty=2, col='black'); legend("topright", legend=exp_names, lty=lty, col=col_palette)
      closepdf("g22_2experiments.pdf");  cat('\t\'g22_2experiments.pdf\' created.\n')
    }
    if (plot_2experiments & (length(g11_all)>0)){
      openpdf("g11_2experiments.pdf")
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''), 
           ylab=eval(bquote(expression(g[.(levels1)][','][.(levels1)](r)))), type='n', ylim = ylim)
      for (ii in 1:length(exp_names)){
        subset <- subset(data_allrange, data$experiment %in% exp_names[ii])
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, g11, lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
      }
      abline(h=1, lty=2, col='black'); legend("topright", legend=exp_names, lty=lty, col=col_palette)
      closepdf("g11_2experiments.pdf"); cat('\t\'g11_2experiments.pdf\' created.\n')
    }
    ylim=c(0,0.5*max(data$K12, na.rm=TRUE))
    if (plot_2experiments & (length(K12_all)>0)){
      openpdf("K12_2experiments.pdf")
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''), 
           ylab=eval(bquote(expression(K[.(levels1)][','][.(levels2)](r)))), type='n', ylim = ylim)
      for (ii in 1:length(exp_names)){
        subset <- subset(data_allrange, data$experiment %in% exp_names[ii])
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, K12, lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
      }
      lines(r_eval, pi*r_eval^2, col="black", lty=2); legend("topleft", legend=exp_names, lty=lty, col=col_palette)
      closepdf("K12_2experiments.pdf"); cat('\t\'K12_2experiments.pdf\' created.\n')
      
      openpdf("cdf12_empirical_npoints_2experiments.pdf")
      plot(1, xlim=c(r_eval[r_eval_plot][1],2*150),#min(parameters_g11$r0, na.rm=TRUE)), 
           ylim=c(0,100), xlab=paste('2r ','[', units,']', sep = ''),
           ylab=paste('Av. number of points (', levels1, ',', levels2, ') [%]', sep=''), type='n')
      for (ii in 1:length(exp_names)){
        subset <- data$experiment %in% exp_names[ii]
        K12_median <- apply(matrix(data_allrange$K12[subset], nrow=length(unique(data$id[subset]))), 2, median)
        lines(2*r_eval[r_eval_plot], 100*K12_median[r_eval_plot]/K12_median[which(r_eval>100)[1]], 
              lty=lty[ii], main='', lwd=1, col=col_palette[ii])
        
        # subset <- subset(data_allrange, data$experiment %in% exp_names[ii])
        # for (jj in unique(subset$id)){
        #   with(subset[which(subset$id == jj), ],
        #        lines(2*r, 100*K12/K12[which(r>200)[1]], lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        # }
      }
      legend("topleft", legend=exp_names, lty=lty, col=col_palette);
      closepdf("cdf12_empirical_npoints_2experiments.pdf"); cat('\t\'cdf12_empirical_npoints_2experiments.pdf\' created.\n')
    }
    if (save){
      save(data, data_allrange, units, parameters_g12, parameters_g22, parameters_g11, density, npoints, 
           r_eval, levels1, levels2, 
           file=paste(c(paste(c(path_to_experiments[1], paste(c(levels1,levels2), collapse = '')), collapse='/'), 
                        '_2experiments.RData'), collapse = ''))
    }
  }
  if (pptype == "unmarked"){
    load(path_to_RData1)
    r_eval <- r_eval*unit_size
    data.g11.exp1 <- t(sapply(g11_all, `[[`, 'g11'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]
    data.K11.exp1 <- t(sapply(K11_all, `[[`, 'K11'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]*(unit_size^2)
    g11_r0_exp1 <- unit_size*g11_r0_all
    g11_rend_exp1 <- which(r_eval[which(r_eval>=reso_r)]>=quantile(g11_r0_exp1, probs=stat_rrange_probs, na.rm=TRUE))[1]
    intensities.exp1 <- intensities/(unit_size^2)
    npoints.exp1 <- npoints
    cat('for statistical analysis (exp1), fill with NA function values above r = ', r_eval[which(r_eval>=reso_r)][g11_rend_exp1], '\n')
    
    if (length(path_to_experiments) > 1){
      load(path_to_RData2)
      r_eval <- r_eval*unit_size
      data.g11.exp2 <- t(sapply(g11_all, `[[`, 'g11'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]
      data.K11.exp2 <- t(sapply(K11_all, `[[`, 'K11'))[,which(r_eval>=reso_r)[1]:(length(r_eval))]*(unit_size^2)
      g11_r0_exp2 <- unit_size*g11_r0_all
      g11_rend_exp2 <- which(r_eval[which(r_eval>=reso_r)]>=quantile(g11_r0_exp2, probs=stat_rrange_probs, na.rm=TRUE))[1]
      intensities.exp2 <- intensities/(unit_size^2)
      npoints.exp2 <- npoints
      cat('for statistical analysis (exp2), fill with NA function values above r = ', r_eval[which(r_eval>=reso_r)][g11_rend_exp2], '\n')
    } else {
      data.g11.exp2 <- NULL
      data.K11.exp2 <- NULL
      g11_r0_exp2 <- NULL
      g11_rend_exp2 <- NULL
      intensities.exp2 <- NULL
      npoints.exp2 <- NULL
    }
    data_temp <- data.frame( g11 = (rbind(data.g11.exp1, data.g11.exp2)), 
                             K11 = (rbind(data.K11.exp1, data.K11.exp2)))
    data_wide <- data.frame(id = factor(1:nrow(data_temp)),
                            experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(data.g11.exp1), nrow(data.g11.exp2)))),
                            data_temp)
    data_wide$experiment <- relevel(data_wide$experiment, ref=exp_name1)
    
    parameters_g11 <- data.frame(id = factor(1:nrow(data_temp)),
                                 experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(data.g11.exp1),
                                                                                    nrow(data.g11.exp2)))), 
                                 r0 = c(g11_r0_exp1, g11_r0_exp2))
    density <- data.frame(id = factor(1:nrow(data_temp)),
                          experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(intensities.exp1),
                                                                             nrow(intensities.exp2)))), 
                          level1 = c(intensities.exp1[,1], intensities.exp2[,1]))
    npoints <- data.frame(id = factor(1:nrow(data_temp)),
                          experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(npoints.exp1),
                                                                             nrow(npoints.exp2)))), 
                          level1 = c(npoints.exp1[,1], npoints.exp2[,1]))
    
    # convert wide-formatted data into long
    num_measurements <- dim(data.g11.exp1)[2]
    data_allrange <- reshape(data_wide, varying=list(names(data_wide)[3:(3+num_measurements-1)],
                                                     names(data_wide)[(3+num_measurements):(3+2*num_measurements-1)]), 
                             idvar=c("id", "experiment"),
                             direction="long", timevar='r', times=r_eval[(length(r_eval)-(dim(data.g11.exp1)[2])+1):length(r_eval)], 
                             v.names=c('g11', 'K11'))
    data_allrange$experiment <- relevel(data_allrange$experiment,ref=exp_name1)
    
    # data_wide[which(data_wide$experiment==exp_name1), (g11_rend_exp1+2):length(data_wide[1,])] <- NA
    # data_wide[which(data_wide$experiment==exp_name2), (g11_rend_exp2+2):length(data_wide[1,])] <- NA
    data_wide[which(data_wide$experiment==exp_name1), (2+g11_rend_exp1):(2+num_measurements)] <- NA
    if (length(path_to_experiments) > 1){
      data_wide[which(data_wide$experiment==exp_name2), (2+g11_rend_exp2):(2+num_measurements)] <- NA
    }
    data <- reshape(data_wide, varying=list(names(data_wide)[3:(3+num_measurements-1)],
                                            names(data_wide)[(3+num_measurements):(3+2*num_measurements-1)]), 
                    idvar=c("id", "experiment"),
                    direction="long", timevar='r', times=r_eval[(length(r_eval)-(dim(data.g11.exp1)[2])+1):length(r_eval)], 
                    v.names=c('g11', 'K11'))
    data$experiment <- relevel(data_allrange$experiment,ref=exp_name1)
    
    r_eval_plot <- which(r_eval>=reso_r)[1]:(0.8*length(r_eval))
    ylim=c(0.5,10)
    # r_eval_plot <- which(r_eval>=reso_r)[1]:(0.5*length(r_eval)); ylim=c(0.5,2.8)
    if (plot_2experiments & (length(g11_all)>0)){
      openpdf(paste(c("g11", filename_experiments, '.pdf'),collapse=''))
      # col_palette <- colorRampPalette(c("black", "grey60"))(length(exp_names)); lty=c(1,2)
      col_palette <-  c(rgb(173,216,230,max = 255), rgb(255,165,0,max = 255,alpha=125)); lty=c(1,2)
      xlim <- c(r_eval[r_eval_plot][1],r_eval[r_eval_plot][length(r_eval[r_eval_plot])])
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''), 
           ylab=eval(bquote(expression(g[.(levels1)][','][.(levels1)](r)))), type='n', ylim = ylim)
      for (ii in 1:length(exp_names)){
        exp_name <- exp_names[ii]
        subset <- subset(data_allrange, data$experiment %in% exp_name)
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, g11, lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
      }
      abline(h=1, lty=2, col='black')
      legend("topright", legend=exp_names, lty=lty, col=col_palette)
      closepdf(paste(c("g11", filename_experiments, '.pdf'),collapse=''))
      cat('\t\'g11_experiment.pdf\' created.\n')
    }
    ylim=c(0,60)
    if (plot_2experiments & (length(K11_all)>0)){
      openpdf(paste(c("K11", filename_experiments, '.pdf'),collapse=''))
      # col_palette <- colorRampPalette(c("black", "grey60"))(length(exp_names)); lty=c(1,2)
      col_palette <-  c(rgb(173,216,230,max = 255), rgb(255,165,0,max = 255,alpha=125)); lty=c(1,2)
      xlim <- c(r_eval[r_eval_plot][1],r_eval[r_eval_plot][length(r_eval[r_eval_plot])])
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''), 
           ylab=eval(bquote(expression(K[.(levels1)][','][.(levels1)](r)))), type='n', ylim = ylim)
      for (ii in 1:length(exp_names)){
        exp_name <- exp_names[ii]
        subset <- subset(data_allrange, data$experiment %in% exp_name)
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, K11, lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
      }
      lines(r_eval[r_eval_plot], pi*r_eval[r_eval_plot]^2, col="black", lty=2)
      legend("topleft", legend=exp_names, lty=lty, col=col_palette)
      closepdf(paste(c("K11", filename_experiments, '.pdf'),collapse=''))
      cat('\t\'K11_experiment.pdf\' created.\n')
      
      openpdf(paste(c("cdf_empirical_npoints", filename_experiments, '.pdf'), collapse=''))
      col_palette <-  c(rgb(173,216,230,max = 255), rgb(255,165,0,max = 255,alpha=125)); lty=c(1,2)
      xlim <- c(r_eval[r_eval_plot][1],2*min(parameters_g11$r0, na.rm=TRUE))
      # r_eval[r_eval_plot][length(r_eval[r_eval_plot])])
      plot(1, xlim=xlim, ylim=c(0,100), xlab=paste('r ','[', units,']', sep = ''), 
           ylab=paste('Av. number of points (', levels1, ',', levels1, ') [%]', sep=''), type='n')
      for (ii in 1:length(exp_names)){
        exp_name <- exp_names[ii]
        subset <- subset(data_allrange, data$experiment %in% exp_name)
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], 
               lines(2*r, 100*K11/K11[which(r>100)[1]], lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
      }
      legend("bottomright", legend=exp_names, lty=lty, col=col_palette)
      closepdf(paste(c("cdf_empirical_npoints", filename_experiments, '.pdf'), collapse=''))
      cat('\t\'cdf_empirical_npoints_2experiments.pdf\' created.\n')
    }
    if (save){
      save(data, data_allrange, units, parameters_g11, density, npoints, r_eval, levels1, levels2, 
           file=paste(c(paste(c(path_to_experiments[1], paste(c(levels1,levels2), collapse = '')), collapse='/'), 
                        filename_experiments, '.RData'), collapse = ''))
    }
  }
}

statistical_analysis <- function(path_to_experiments, pptype, fitting='exp', save=TRUE){
  
  library(lme4)
  library(nlme)
  library(MASS)
  library(car)
  library(lattice)
  library(pracma) # error function
  
  cat('Running statistical analysis...\n')
  if (length(path_to_experiments) > 1){
    filename_experiments <- "_2experiments"
    RData_2experiments <- list.files(path = path_to_experiments[1], pattern = "\\_2experiments.RData$", all.files = FALSE,
                                     full.names = FALSE, recursive = FALSE,
                                     ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  } else{
    filename_experiments <- "_1experiment"
    RData_2experiments <- list.files(path = path_to_experiments[1], pattern = "\\_1experiment.RData$", all.files = FALSE,
                                     full.names = FALSE, recursive = FALSE,
                                     ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  }
  load(paste(c(path_to_experiments[1], RData_2experiments), collapse='/'))
  attach(data)
  
  reso_r <- data$r[1]
  
  if (length(path_to_experiments) > 1){
    exp1_name <- levels(data$experiment)[1]; exp2_name <- levels(data$experiment)[2];
  } else{
    exp1_name <- levels(data$experiment)[1]; exp2_name <- NULL;
  }
  
  data$r <- data$r
  if(fitting %in% c('lin')){
    attach(data)
    cat('linear transformation... Done\n')
  }
  if (pptype == "marked"){
    
    ## exponential fitting
    if(fitting %in% c('exp', 'expsq')){
      data$g12 <- log(data$g12-1)
      data$g22 <- log(data$g22-1)
      data$g11 <- log(data$g11-1)
      attach(data)
      cat('Exponential transformation... Done\n')
    }
    
    ## response feature analysis (Faraway, 2006, p.206)
    # fit ordinary regressions separately to each group
    parameters_g12 <- data.frame(parameters_g12, 
                                 beta0=matrix(nrow=nlevels(data$id),ncol=1), 
                                 beta1=matrix(nrow=nlevels(data$id),ncol=1), 
                                 R2=matrix(nrow=nlevels(data$id),ncol=1))
    K_r0_g12 <- numeric(nlevels(data$id))
    
    parameters_g22 <- data.frame(parameters_g22, 
                                 beta0=matrix(nrow=nlevels(data$id),ncol=1), 
                                 beta1=matrix(nrow=nlevels(data$id),ncol=1), 
                                 R2=matrix(nrow=nlevels(data$id),ncol=1))
    K_r0_g22 <- numeric(nlevels(data$id))
    
    parameters_g11 <- data.frame(parameters_g11, 
                                 beta0=matrix(nrow=nlevels(data$id),ncol=1), 
                                 beta1=matrix(nrow=nlevels(data$id),ncol=1), 
                                 R2=matrix(nrow=nlevels(data$id),ncol=1))
    for(i in 1:nlevels(data$id)){
      if (fitting %in% c("lin", "exp")){lmod <- lm(g12 ~ r, subset=(id==i), data=data)}
      if (fitting == "expsq"){lmod <- lm(g12 ~ 1 + I(r^2), subset=(id==i), data=data)}
      parameters_g12$R2[i] <- summary(lmod)$adj.r.squared
      parameters_g12$beta0[i] <- coef(lmod)[1]
      parameters_g12$beta1[i] <- coef(lmod)[2]
      a <- tail(which(r_eval < 0.5*parameters_g12$r0[which(parameters_g12$id == i)]),n=1)
      if (length(a)==0){ K_r0_g12[i] <- NA}
      else{
        K_r0_g12[i] <- K12[id == i][a]
      }
      if (fitting %in% c("lin", "exp")){lmod <- lm(g22 ~ r, subset=(id==i), data=data)}
      if (fitting == "expsq"){lmod <- lm(g22 ~ 1 + I(r^2), subset=(id==i), data=data)}
      parameters_g22$R2[i] <- summary(lmod)$adj.r.squared
      parameters_g22$beta0[i] <- coef(lmod)[1]
      parameters_g22$beta1[i] <- coef(lmod)[2]
      a <- tail(which(r_eval < 0.5*parameters_g22$r0[which(parameters_g22$id == i)]),n=1)
      if (length(a)==0){ K_r0_g22[i] <- NA}
      else{
        K_r0_g22[i] <- K22[id == i][a]
      }
      if (fitting %in% c("lin", "exp")){lmod <- lm(g11 ~ r, subset=(id==i), data=data)}
      if (fitting == "expsq"){lmod <- lm(g11 ~ 1 + I(r^2), subset=(id==i), data=data)}
      parameters_g11$R2[i] <- summary(lmod)$adj.r.squared
      parameters_g11$beta0[i] <- coef(lmod)[1]
      parameters_g11$beta1[i] <- coef(lmod)[2]
      a <- tail(which(r_eval < 0.5*parameters_g11$r0[which(parameters_g11$id == i)]),n=1)
    }
    cat('Checking quality of fitting... ')
    fileName <- paste(c("../output/", unlist(strsplit(RData_2experiments, "\\_"))[1], "_summaryfitting_", fitting, '.txt'), collapse = "")
    if (file.exists(fileName)) file.remove(fileName)    
    if (fitting %in% c("exp","lin")){
      lmod <- lm(g11 ~ r, subset=(experiment==levels(data$experiment)[1]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g11_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g11_plots.pdf'), collapse=""))
      lmod <- lm(g11 ~ r, subset=(experiment==levels(data$experiment)[2]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[2], '_g11_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[2], '_g11_plots.pdf'), collapse=""))
      lmod <- lm(g22 ~ r, subset=(experiment==levels(data$experiment)[1]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g22_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g22_plots.pdf'), collapse=""))
      lmod <- lm(g22 ~ r, subset=(experiment==levels(data$experiment)[2]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[2], '_g22_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[2], '_g22_plots.pdf'), collapse=""))
      lmod <- lm(g12 ~ r, subset=(experiment==levels(data$experiment)[1]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g12_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g12_plots.pdf'), collapse=""))
      lmod <- lm(g12 ~ r, subset=(experiment==levels(data$experiment)[2]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[2], '_g12_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[2], '_g12_plots.pdf'), collapse=""))
      cat('Done.\n')
    }
    if (fitting=="expsq"){
      lmod <- lm(g11 ~ 1+I(r^2), subset=(experiment==levels(data$experiment)[1]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g11_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g11_plots.pdf'), collapse=""))
      lmod <- lm(g11 ~ 1+I(r^2), subset=(experiment==levels(data$experiment)[2]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[2], '_g11_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[2], '_g11_plots.pdf'), collapse=""))
      lmod <- lm(g22 ~ 1+I(r^2), subset=(experiment==levels(data$experiment)[1]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g22_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g22_plots.pdf'), collapse=""))
      lmod <- lm(g22 ~ 1+I(r^2), subset=(experiment==levels(data$experiment)[2]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[2], '_g22_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[2], '_g22_plots.pdf'), collapse=""))
      lmod <- lm(g12 ~ 1+I(r^2), subset=(experiment==levels(data$experiment)[1]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g12_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g12_plots.pdf'), collapse=""))
      lmod <- lm(g12 ~ 1+I(r^2), subset=(experiment==levels(data$experiment)[2]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[2], '_g12_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[2], '_g12_plots.pdf'), collapse=""))
      cat('Done.\n')
    }
    if(fitting %in% c('exp', 'expsq')){
      amplitude_g12 <- exp(parameters_g12$beta0)
      amplitudeAtr0_g12 <- exp(parameters_g12$beta0*parameters_g12$beta1*reso_r)+1
      length_scale_g12 <- -parameters_g12$beta1^(-1)
      area_g12 <- -length_scale_g12*amplitude_g12*(exp(-1/(length_scale_g12)*parameters_g12$r0)-1)
      
      amplitude_g22 <- exp(parameters_g22$beta0)
      amplitudeAtr0_g22 <- exp(parameters_g22$beta0*parameters_g22$beta1*reso_r)+1
      length_scale_g22 <- -parameters_g22$beta1^(-1)
      area_g22 <- -length_scale_g22*amplitude_g22*(exp(-1/(length_scale_g22)*parameters_g22$r0)-1)
      
      amplitude_g11 <- exp(parameters_g11$beta0)
      amplitudeAtr0_g11 <- exp(parameters_g11$beta0*parameters_g11$beta1*reso_r)+1
      length_scale_g11 <- -parameters_g11$beta1^(-1)
      area_g11 <- -length_scale_g11*amplitude_g11*(exp(-1/(length_scale_g11)*parameters_g11$r0)-1)
      
      if (fitting == 'expsq'){
        clusterradiusR_g12 <- sqrt(length_scale_g12)
        clusterradiusR_g22 <- sqrt(length_scale_g22)
        clusterradiusR_g11 <- sqrt(length_scale_g11)
        kappa_g12 <- 1/(amplitude_g12*pi*length_scale_g12)  # number of clusters per area
        kappa_g22 <- 1/(amplitude_g22*pi*length_scale_g22)  # number of clusters per area
        kappa_g11 <- 1/(amplitude_g11*pi*length_scale_g11)  # number of clusters per area
        phi_g12 <- amplitude_g12/4   # rho_cluster/rho_average
        phi_g22 <- amplitude_g22/4   # rho_cluster/rho_average
        phi_g11 <- amplitude_g11/4   # rho_cluster/rho_average
        rho_g12 <- phi_g12*(density$level1 + density$level2)
        rho_g22 <- phi_g22*density$level2
        rho_g11 <- phi_g11*density$level1
        Nclusters_g22 <- density$level2/kappa_g22  #  average number of points per cluster
        Nclusters_g11 <- density$level1/kappa_g11  #  average number of points per cluster
      }
      if (fitting == 'exp'){
        clusterradiusR_g12 <- length_scale_g12
        clusterradiusR_g22 <- length_scale_g22
        clusterradiusR_g11 <- length_scale_g11
        phi_g12 <- 2*amplitude_g12   # rho_cluster/rho_average
        phi_g22 <- 2*amplitude_g22   # rho_cluster/rho_average
        phi_g11 <- 2*amplitude_g11   # rho_cluster/rho_average
        rho_g12 <- phi_g12*(density$level1 + density$level2)
        rho_g22 <- phi_g22*density$level2
        rho_g11 <- phi_g11*density$level1
        Nclusters_g22 <- 2*amplitude_g22*pi*length_scale_g22^2*density$level2  #  average number of points per cluster
        Nclusters_g11 <- 2*amplitude_g11*pi*length_scale_g11^2*density$level1  #  average number of points per cluster
      }  
      
      # tests (experiments)
      cat('Running tests... \n')
      pexp <- data$experiment[match(1:nlevels(data$id),data$id)]
      
      col_palette <-  c(rgb(173,216,230,max = 255,alpha=125), rgb(255,165,0,max = 255,alpha=125)); lty=c(1,2)
      r_eval_plot <- which(r_eval>=reso_r)[1]:(0.8*length(r_eval))
      xlim <- c(r_eval[r_eval_plot][1],r_eval[r_eval_plot][length(r_eval[r_eval_plot])])
      ylim=c(0.5,10)
      openpdf(paste(c("g12_2experiments_", fitting, '.pdf'), collapse=""))
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''),
           ylab=eval(bquote(expression(g[.(levels1)][','][.(levels2)](r)))), type='n', ylim = ylim)
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        subset <- subset(data_allrange, data$experiment %in% exp_name)
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, g12, lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
        abline(h=1, lty=2, col='black')
        legend("topright", legend=exp_names, lty=lty, col=col_palette)
      }
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        medianr0 <- median(parameters_g11$r0[which(parameters_g12$experiment == exp_name)], na.rm=TRUE)
        medianlengthscale <- median(length_scale_g12[pexp==exp_name], na.rm=TRUE)
        medianamplitude <- median(amplitude_g12[pexp==exp_name], na.rm=TRUE)
        if (fitting=="exp"){fitted_line <- 1+medianamplitude*exp(-seq(0, medianr0, by=0.5)/medianlengthscale)}
        if (fitting=="expsq"){fitted_line <- 1+medianamplitude*exp(-seq(0, medianr0, by=0.5)^2/medianlengthscale)}
        if (fitting=="lin"){fitted_line <- medianamplitude - 1/medianlengthscale*seq(0, medianr0, by=0.5)}
        lines(seq(0, medianr0, by=0.5),fitted_line , col='dimgray', lty=1)
      }
      closepdf(paste(c("g12_2experiments_", fitting, '.pdf'), collapse=""))
      cat('\t\'g12_2experiments.pdf\' created.\n')
      openpdf(paste(c("g11_2experiments_", fitting, '.pdf'), collapse=""))
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''),
           ylab=eval(bquote(expression(g[.(levels1)][','][.(levels1)](r)))), type='n', ylim = ylim)
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        subset <- subset(data_allrange, data$experiment %in% exp_name)
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, g11, lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
        abline(h=1, lty=2, col='black')
        legend("topright", legend=exp_names, lty=lty, col=col_palette)
      }
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        medianr0 <- median(parameters_g11$r0[which(parameters_g11$experiment == exp_name)], na.rm=TRUE)
        medianlengthscale <- median(length_scale_g11[pexp==exp_name], na.rm=TRUE)
        medianamplitude <- median(amplitude_g11[pexp==exp_name], na.rm=TRUE)
        if (fitting=="exp"){fitted_line <- 1+medianamplitude*exp(-seq(0, medianr0, by=0.5)/medianlengthscale)}
        if (fitting=="expsq"){fitted_line <- 1+medianamplitude*exp(-seq(0, medianr0, by=0.5)^2/medianlengthscale)}
        if (fitting=="lin"){fitted_line <- medianamplitude - 1/medianlengthscale*seq(0, medianr0, by=0.5)}
        lines(seq(0, medianr0, by=0.5),fitted_line , col='dimgray', lty=1)
      }
      closepdf(paste(c("g11_2experiments_", fitting, '.pdf'), collapse=""))
      cat('\t\'g11_2experiments.pdf\' created.\n')
      openpdf(paste(c("g22_2experiments_", fitting, '.pdf'), collapse=""))
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''),
           ylab=eval(bquote(expression(g[.(levels2)][','][.(levels2)](r)))), type='n', ylim = ylim)
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        subset <- subset(data_allrange, data$experiment %in% exp_name)
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, g22, lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
        abline(h=1, lty=2, col='black')
        legend("topright", legend=exp_names, lty=lty, col=col_palette)
      }
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        medianr0 <- median(parameters_g11$r0[which(parameters_g22$experiment == exp_name)], na.rm=TRUE)
        medianlengthscale <- median(length_scale_g22[pexp==exp_name], na.rm=TRUE)
        medianamplitude <- median(amplitude_g22[pexp==exp_name], na.rm=TRUE)
        if (fitting=="exp"){fitted_line <- 1+medianamplitude*exp(-seq(0, medianr0, by=0.5)/medianlengthscale)}
        if (fitting=="expsq"){fitted_line <- 1+medianamplitude*exp(-seq(0, medianr0, by=0.5)^2/medianlengthscale)}
        if (fitting=="lin"){fitted_line <- medianamplitude - 1/medianlengthscale*seq(0, medianr0, by=0.5)}
        lines(seq(0, medianr0, by=0.5),fitted_line , col='dimgray', lty=1)
      }
      closepdf(paste(c("g22_2experiments_", fitting, '.pdf'), collapse=""))
      cat('\t\'g22_2experiments.pdf\' created.\n')
      openpdf(paste(c("g11_log_2experiments_", fitting, '.pdf'), collapse=""))
      plot(1, xlim=c(0,median(parameters_g12$r0, na.rm=TRUE)), 
           ylim=range(data$g11, na.rm=TRUE), xlab=paste('r ','[', units,']', sep = ''),
           ylab=eval(bquote(expression(ln(g[.(levels1)][','][.(levels1)](r)-1)))), type='n')
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        subset <- subset(data, data$experiment %in% exp_name)
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, g11, lty=lty[ii], main='', pch= 20, lwd=1, col=col_palette[ii]))
        }
        legend("topright", legend=exp_names, lty=lty, col=col_palette)
      }
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        medianr0 <- median(parameters_g11$r0[which(parameters_g11$experiment == exp_name)], na.rm=TRUE)
        medianlengthscale <- median(length_scale_g11[pexp==exp_name], na.rm=TRUE)
        medianamplitude <- median(amplitude_g11[pexp==exp_name], na.rm=TRUE)
        if (fitting=="exp"){fitted_line <- log(medianamplitude)-seq(0, medianr0, by=0.5)/medianlengthscale}
        if (fitting=="expsq"){fitted_line <- log(medianamplitude)-seq(0, medianr0, by=0.5)^2/medianlengthscale}
        if (fitting=="lin"){fitted_line <- medianamplitude - 1/medianlengthscale*seq(0, medianr0, by=0.5)}
        lines(seq(0, medianr0, by=0.5),fitted_line , col='dimgray', lty=1)
      }
      closepdf(paste(c("g11_log_2experiments_", fitting, '.pdf'), collapse=""))
      openpdf(paste(c("g22_log_2experiments_", fitting, '.pdf'), collapse=""))
      plot(1, xlim=c(0,median(parameters_g22$r0, na.rm=TRUE)), 
           ylim=range(data$g22, na.rm=TRUE), xlab=paste('r ','[', units,']', sep = ''),
           ylab=eval(bquote(expression(ln(g[.(levels2)][','][.(levels2)](r)-1)))), type='n')
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        subset <- subset(data, data$experiment %in% exp_name)
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, g22, lty=lty[ii], main='', pch= 20, lwd=1, col=col_palette[ii]))
        }
        legend("topright", legend=exp_names, lty=lty, col=col_palette)
      }
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        medianr0 <- median(parameters_g22$r0[which(parameters_g22$experiment == exp_name)], na.rm=TRUE)
        medianlengthscale <- median(length_scale_g22[pexp==exp_name], na.rm=TRUE)
        medianamplitude <- median(amplitude_g22[pexp==exp_name], na.rm=TRUE)
        if (fitting=="exp"){fitted_line <- log(medianamplitude)-seq(0, medianr0, by=0.5)/medianlengthscale}
        if (fitting=="expsq"){fitted_line <- log(medianamplitude)-seq(0, medianr0, by=0.5)^2/medianlengthscale}
        if (fitting=="lin"){fitted_line <- medianamplitude - 1/medianlengthscale*seq(0, medianr0, by=0.5)}
        lines(seq(0, medianr0, by=0.5),fitted_line , col='dimgray', lty=1)
      }
      closepdf(paste(c("g22_log_2experiments_", fitting, '.pdf'), collapse=""))
      openpdf(paste(c("g12_log_2experiments_", fitting, '.pdf'), collapse=""))
      plot(1, xlim=c(0,median(parameters_g12$r0, na.rm=TRUE)), 
           ylim=range(data$g12, na.rm=TRUE), xlab=paste('r ','[', units,']', sep = ''),
           ylab=eval(bquote(expression(ln(g[.(levels1)][','][.(levels2)](r)-1)))), type='n')
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        subset <- subset(data, data$experiment %in% exp_name)
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, g12, lty=lty[ii], main='', pch= 20, lwd=1, col=col_palette[ii]))
        }
        legend("topright", legend=exp_names, lty=lty, col=col_palette)
      }
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        medianr0 <- median(parameters_g12$r0[which(parameters_g12$experiment == exp_name)], na.rm=TRUE)
        medianlengthscale <- median(length_scale_g12[pexp==exp_name], na.rm=TRUE)
        medianamplitude <- median(amplitude_g12[pexp==exp_name], na.rm=TRUE)
        if (fitting=="exp"){fitted_line <- log(medianamplitude)-seq(0, medianr0, by=0.5)/medianlengthscale}
        if (fitting=="expsq"){fitted_line <- log(medianamplitude)-seq(0, medianr0, by=0.5)^2/medianlengthscale}
        if (fitting=="lin"){fitted_line <- medianamplitude - 1/medianlengthscale*seq(0, medianr0, by=0.5)}
        lines(seq(0, medianr0, by=0.5),fitted_line , col='dimgray', lty=1)
      }
      closepdf(paste(c("g12_log_2experiments_", fitting, '.pdf'), collapse=""))
      cat('\t\'g12/g11/g22_2experiments(FITTING).pdf\' created.\n')
      openpdf('cdf11_npoints_2experiments.pdf')
      plot(1, xlim=c(0,2*min(parameters_g11$r0, na.rm=TRUE)), ylim=c(0,100), xlab=paste('2r ','[', units,']', sep = ''),
           ylab=paste('Av. number of points (', levels1, ',', levels1, ') [%]', sep=''), type='n')
      for (ii in  1:length(levels(data$experiment))){
        if (fitting=="expsq"){
          exp_name <- levels(data$experiment)[ii]
          medianclusterradiusR <- median(clusterradiusR_g11[pexp==exp_name], na.rm=TRUE)
          r <- seq(0,2*min(parameters_g11$r0[pexp==exp_name], na.rm=TRUE), by=5)
          cdf_line <- 100*(pnorm(r, mean=0, sd=medianclusterradiusR/4)-pnorm(-r, mean=0, sd=medianclusterradiusR/4))
          lines(2*r, cdf_line , col=col_palette[ii], lty=1)
          legend("bottomright", legend=exp_names, lty=1, col=col_palette)
        }
      }
      closepdf("cdf11_npoints_2experiments.pdf")
      cat('\t\'cdf11_npoints_2experiments(FITTING).pdf\' created.\n')
      openpdf('cdf22_npoints_2experiments.pdf')
      plot(1, xlim=c(0,2*min(parameters_g11$r0, na.rm=TRUE)), ylim=c(0,100), xlab=paste('2r ','[', units,']', sep = ''),
           ylab=paste('Av. number of points (', levels2, ',', levels2, ') [%]', sep=''), type='n')
      for (ii in  1:length(levels(data$experiment))){
        if (fitting=="expsq"){
          exp_name <- levels(data$experiment)[ii]
          medianclusterradiusR <- median(clusterradiusR_g22[pexp==exp_name], na.rm=TRUE)
          r <- seq(0,2*min(parameters_g22$r0[pexp==exp_name], na.rm=TRUE), by=5)
          cdf_line <- 100*(pnorm(r, mean=0, sd=medianclusterradiusR/4)-pnorm(-r, mean=0, sd=medianclusterradiusR/4))
          lines(2*r, cdf_line , col=col_palette[ii], lty=1)
          legend("bottomright", legend=exp_names, lty=1, col=col_palette)
        }
      }
      closepdf("cdf22_npoints_2experiments.pdf")
      cat('\t\'cdf22_npoints_2experiments(FITTING).pdf\' created.\n')
      openpdf('cdf12_npoints_2experiments.pdf')
      plot(1, xlim=c(0,2*min(parameters_g11$r0, na.rm=TRUE)), ylim=c(0,100), xlab=paste('2r ','[', units,']', sep = ''),
           ylab=paste('Av. number of points (', levels1, ',', levels2, ') [%]', sep=''), type='n')
      for (ii in  1:length(levels(data$experiment))){
        if (fitting=="expsq"){
          exp_name <- levels(data$experiment)[ii]
          medianclusterradiusR <- median(clusterradiusR_g12[pexp==exp_name], na.rm=TRUE)
          r <- seq(0,2*min(parameters_g12$r0[pexp==exp_name], na.rm=TRUE), by=1)
          cdf_line <- 100*(pnorm(r, mean=0, sd=medianclusterradiusR/4)-pnorm(-r, mean=0, sd=medianclusterradiusR/4))
          lines(2*r, cdf_line , col=col_palette[ii], lty=1)
          legend("bottomright", legend=exp_names, lty=1, col=col_palette)
        }
      }
      closepdf("cdf12_npoints_2experiments.pdf")
      cat('\t\'cdf12_npoints_2experiments(FITTING).pdf\' created.\n')
      
      parameters_g12_r0.test <- t.test(parameters_g12$r0[which(parameters_g12$experiment == exp1_name)], parameters_g12$r0[which(parameters_g12$experiment == exp2_name)])
      density_level1.test <- t.test(density$level1[pexp==exp1_name],density$level1[pexp==exp2_name])
      density_level2.test <- t.test(density$level2[pexp==exp1_name],density$level2[pexp==exp2_name])
      amplitude_g12.test <- t.test(amplitude_g12[pexp==exp1_name],amplitude_g12[pexp==exp2_name])
      amplitudeAtr0_g12.test <- t.test(amplitudeAtr0_g12[pexp==exp1_name],amplitudeAtr0_g12[pexp==exp2_name])
      length_scale_g12.test <- t.test(length_scale_g12[pexp==exp1_name],length_scale_g12[pexp==exp2_name])
      clusterradiusR_g12.test <- t.test(clusterradiusR_g12[pexp==exp1_name],clusterradiusR_g12[pexp==exp2_name])
      parameters_g12_R2.test <- t.test(parameters_g12$R2[pexp==exp1_name],parameters_g12$R2[pexp==exp2_name])
      area_g12.test <- t.test(area_g12[pexp==exp1_name],area_g12[pexp==exp2_name])
      phi_g12.test <- t.test(phi_g12[pexp==exp1_name],phi_g12[pexp==exp2_name])
      rho_g12.test <- t.test(rho_g12[pexp==exp1_name],rho_g12[pexp==exp2_name])
      
      parameters_g22_r0.test <- t.test(parameters_g22$r0[which(parameters_g22$experiment == exp1_name)], parameters_g22$r0[which(parameters_g22$experiment == exp2_name)])
      amplitude_g22.test <- t.test(amplitude_g22[pexp==exp1_name],amplitude_g22[pexp==exp2_name])
      amplitudeAtr0_g22.test <- t.test(amplitudeAtr0_g22[pexp==exp1_name],amplitudeAtr0_g22[pexp==exp2_name])
      length_scale_g22.test <- t.test(length_scale_g22[pexp==exp1_name],length_scale_g22[pexp==exp2_name])
      clusterradiusR_g22.test <- t.test(clusterradiusR_g22[pexp==exp1_name],clusterradiusR_g22[pexp==exp2_name])
      parameters_g22_R2.test <- t.test(parameters_g22$R2[pexp==exp1_name],parameters_g22$R2[pexp==exp2_name])
      area_g22.test <- t.test(area_g22[pexp==exp1_name],area_g22[pexp==exp2_name])
      phi_g22.test <- t.test(phi_g22[pexp==exp1_name],phi_g22[pexp==exp2_name])
      rho_g22.test <- t.test(rho_g22[pexp==exp1_name],rho_g22[pexp==exp2_name])
      Nclusters_g22.test <- t.test(Nclusters_g22[pexp==exp1_name],Nclusters_g22[pexp==exp2_name])
      
      parameters_g11_r0.test <- t.test(parameters_g11$r0[which(parameters_g11$experiment == exp1_name)], parameters_g11$r0[which(parameters_g11$experiment == exp2_name)])
      amplitude_g11.test <- t.test(amplitude_g11[pexp==exp1_name],amplitude_g11[pexp==exp2_name])
      amplitudeAtr0_g11.test <- t.test(amplitudeAtr0_g11[pexp==exp1_name],amplitudeAtr0_g11[pexp==exp2_name])
      length_scale_g11.test <- t.test(length_scale_g11[pexp==exp1_name],length_scale_g11[pexp==exp2_name])
      clusterradiusR_g11.test <- t.test(clusterradiusR_g11[pexp==exp1_name],clusterradiusR_g11[pexp==exp2_name])
      parameters_g11_R2.test <- t.test(parameters_g11$R2[pexp==exp1_name],parameters_g11$R2[pexp==exp2_name])
      area_g11.test <- t.test(area_g11[pexp==exp1_name],area_g11[pexp==exp2_name])
      phi_g11.test <- t.test(phi_g11[pexp==exp1_name],phi_g11[pexp==exp2_name])
      rho_g11.test <- t.test(rho_g11[pexp==exp1_name],rho_g11[pexp==exp2_name])
      Nclusters_g11.test <- t.test(Nclusters_g11[pexp==exp1_name],Nclusters_g11[pexp==exp2_name])
      
      if (fitting=='expsq'){
        kappa_g12.test <- t.test(kappa_g12[pexp==exp1_name],kappa_g12[pexp==exp2_name])
        kappa_g22.test <- t.test(kappa_g22[pexp==exp1_name],kappa_g22[pexp==exp2_name])
        kappa_g11.test <- t.test(kappa_g11[pexp==exp1_name],kappa_g11[pexp==exp2_name])
      }
      cat('Done.\n')
      
      # # plot results
      cat('Plotting results of statistical tests...\n')
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests.pdf", collapse = "_"), pointsize=11)
      # old <- par(mfrow=c(2, 2))
      
      ## G12
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests", levels1, "density.pdf"), collapse = "_"))
      ylim <- c(0,1.2*max(density$level1, na.rm=TRUE))
      bp <- boxplot(split(density$level1, pexp), names=c(exp1_name, exp2_name), ylim=ylim, ylab=bquote('Average density ' ~ .(levels1) ~ ' [points/' ~ nm^2 ~']'))
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(density$level1, na.rm=TRUE), y1 = 1.04*max(density$level1, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(density$level1, na.rm=TRUE),  get_asterisk(density_level1.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests", levels1, "density.pdf"), collapse = "_"))
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests", levels2, "density.pdf"), collapse = "_"))
      ylim <- c(0,1.2*max(density$level2, na.rm=TRUE))
      bp <- boxplot(split(density$level2, pexp), 
                    names=c(exp1_name, exp2_name), ylim=ylim, ylab=bquote('Average density ' ~ .(levels2) ~ ' [points/' ~ nm^2 ~']'))
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(density$level2, na.rm=TRUE), y1 = 1.04*max(density$level2, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(density$level2, na.rm=TRUE),  get_asterisk(density_level2.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests", levels2, "density.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_density.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_r0_g12.pdf"), collapse = "_"))
      bp <- boxplot(split(parameters_g12$r0,pexp), ylab=expression(paste('Correlation range ', r[0], ' [nm]')), ylim=c(50,1.2*max(parameters_g12$r0, na.rm=TRUE)))
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(parameters_g12$r0, na.rm=TRUE), y1 = 1.04*max(parameters_g12$r0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(parameters_g12$r0, na.rm=TRUE), get_asterisk(parameters_g12_r0.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_r0_g12.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_r0_g12.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitude_g12.pdf"), collapse = "_"))
      ylim <- c(0,1.2*max(amplitude_g12, na.rm=TRUE))
      bp <- boxplot(split(amplitude_g12,pexp), ylab=expression("Amplitude A"), names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(amplitude_g12, na.rm=TRUE), y1 = 1.04*max(amplitude_g12, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(amplitude_g12, na.rm=TRUE),  get_asterisk(amplitude_g12.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitude_g12.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_amplitude_g12.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitudeAtr0_g12.pdf"), collapse = "_"))
      ylim <- c(0,1.2*max(amplitudeAtr0_g12, na.rm=TRUE))
      bp <- boxplot(split(amplitudeAtr0_g12,pexp), ylab=expression("Correlation at 20 nm [a.u.]"), names=c(exp1_name, exp2_name), ylim=ylim)
      # text(bp$group, bp$out, parameters_g11$id[which(amplitudeAtr0_g12 %in% bp$out)], cex=0.5, pos = 4)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(amplitudeAtr0_g12, na.rm=TRUE), y1 = 1.04*max(amplitudeAtr0_g12, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(amplitudeAtr0_g12, na.rm=TRUE),  get_asterisk(amplitudeAtr0_g12.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitudeAtr0_g12.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_amplitudeAtr0_g12.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_lengthscale_g12.pdf"), collapse = "_"))
      ylim <- c(10,1.2*max(length_scale_g12, na.rm=TRUE))
      bp <- boxplot(split(length_scale_g12,pexp), ylab=expression(paste("Length scale ", lambda, " [nm]")), 
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(length_scale_g12, na.rm=TRUE), y1 = 1.04*max(length_scale_g12, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(length_scale_g12, na.rm=TRUE) ,  get_asterisk(length_scale_g12.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_lengthscale_g12.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_lengthscale_g12.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_clusterradiusR_g12.pdf"), collapse = "_"))
      ylim <- c(0.9*min(clusterradiusR_g12, na.rm=TRUE),1.2*max(clusterradiusR_g12, na.rm=TRUE))
      bp <- boxplot(split(clusterradiusR_g12,pexp), ylab="Cluster radius R [nm]", 
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(clusterradiusR_g12, na.rm=TRUE), y1 = 1.04*max(clusterradiusR_g12, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(clusterradiusR_g12, na.rm=TRUE) ,  get_asterisk(clusterradiusR_g12.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_clusterradiusR_g12.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_clusterradiusR_g12.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_phi_g12.pdf"), collapse = "_"))
      ylim <- c(0.9*min(phi_g12, na.rm=TRUE),1.2*max(phi_g12, na.rm=TRUE))
      bp <- boxplot(split(phi_g12,pexp), ylab=expression(phi1^cluster),
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(phi_g12, na.rm=TRUE), y1 = 1.04*max(phi_g12, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(phi_g12, na.rm=TRUE) ,  get_asterisk(phi_g12.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_phi_g12.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_phi_g12.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_rho_g12.pdf"), collapse = "_"))
      ylim <- c(0.9*min(rho_g12, na.rm=TRUE),1.2*max(rho_g12, na.rm=TRUE))
      bp <- boxplot(split(rho_g12,pexp), ylab=expression(rho^cluster),
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(rho_g12, na.rm=TRUE), y1 = 1.04*max(rho_g12, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(rho_g12, na.rm=TRUE) ,  get_asterisk(rho_g12.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_rho_g12.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_rho_g12.pdf\' created.\n')
      
      if (fitting == 'expsq'){
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_kappa_g12.pdf"), collapse = "_"))
        ylim <- c(0.9*min(kappa_g12, na.rm=TRUE),1.2*max(kappa_g12, na.rm=TRUE))
        bp <- boxplot(split(kappa_g12,pexp), ylab=expression(paste(kappa, " [clusters/", nm^2, "]")), 
                      names=c(exp1_name, exp2_name), ylim=ylim)
        segments(x0 = 1, x1 = 2, y0 = 1.04*max(kappa_g12, na.rm=TRUE), y1 = 1.04*max(kappa_g12, na.rm=TRUE), col = "black")
        text( 1.5 , 1.09*max(kappa_g12, na.rm=TRUE) ,  get_asterisk(kappa_g12.test) , cex=1)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_kappa_g12.pdf"), collapse = "_"))
        cat('\t\'statistical_tests_kappa_g12.pdf\' created.\n')
      }
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_R2_g12.pdf"), collapse = "_"))
      # ylim <- c(10,1.2*max(parameters_g12$R2, na.rm=TRUE))
      bp <- boxplot(split(parameters_g12$R2,pexp), ylab=expression(paste(R^2,"(adj)")), names=c(exp1_name, exp2_name))
      text(bp$group, bp$out, parameters_g12$id[which(parameters_g12$R2 %in% bp$out)], cex=0.5, pos = 4)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(parameters_g12$R2, na.rm=TRUE), y1 = 1.04*max(parameters_g12$R2, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(parameters_g12$R2, na.rm=TRUE) ,  get_asterisk(parameters_g12_R2.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_R2_g12.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_R2_g12.pdf\' created.\n')
      
      
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_area_g12.pdf", collapse = "_"))
      # ylim <- c(0,1.2*max(area_g12, na.rm=TRUE))
      # bp <- boxplot(split(area_g12,pexp), main=expression("area"), names=c(exp1_name, exp2_name),  ylim=ylim)
      # segments(x0 = 1, x1 = 2, y0 = 1.04*max(area_g12, na.rm=TRUE), y1 = 1.04*max(area_g12, na.rm=TRUE), col = "black")
      # text( 1.5 , 1.09*max(area_g12, na.rm=TRUE) , get_asterisk(area_g12.test), cex=1.5)
      # closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_area.pdf", collapse = "_"))
      # cat('\t\'statistical_tests_area_g12.pdf\' created.\n')
      # # par(old); closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_g12.pdf", collapse = "_"))
      
      ## G22
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_r0_g22.pdf"), collapse = "_"))
      bp <- boxplot(split(parameters_g22$r0,pexp), ylab=expression(paste('Correlation range ', r[0], ' [nm]')), ylim=c(50,1.2*max(parameters_g22$r0, na.rm=TRUE)))
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(parameters_g22$r0, na.rm=TRUE), y1 = 1.04*max(parameters_g22$r0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(parameters_g22$r0, na.rm=TRUE), get_asterisk(parameters_g22_r0.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_r0_g22.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_r0_g22.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitude_g22.pdf"), collapse = "_"))
      ylim <- c(0,1.2*max(amplitude_g22, na.rm=TRUE))
      bp <- boxplot(split(amplitude_g22,pexp), ylab=expression("Amplitude A"), names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(amplitude_g22, na.rm=TRUE), y1 = 1.04*max(amplitude_g22, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(amplitude_g22, na.rm=TRUE),  get_asterisk(amplitude_g22.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitude_g22.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_amplitude_g22.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitudeAtr0_g22.pdf"), collapse = "_"))
      ylim <- c(0,1.2*max(amplitudeAtr0_g22, na.rm=TRUE))
      bp <- boxplot(split(amplitudeAtr0_g22,pexp), ylab=expression("Correlation at 20 nm [a.u.]"), names=c(exp1_name, exp2_name), ylim=ylim)
      # text(bp$group, bp$out, parameters_g22$id[which(amplitudeAtr0_g22 %in% bp$out)], cex=0.5, pos = 4)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(amplitudeAtr0_g22, na.rm=TRUE), y1 = 1.04*max(amplitudeAtr0_g22, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(amplitudeAtr0_g22, na.rm=TRUE),  get_asterisk(amplitudeAtr0_g22.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitudeAtr0_g22.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_amplitudeAtr0_g22.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_lengthscale_g22.pdf"), collapse = "_"))
      ylim <- c(10,1.2*max(length_scale_g22, na.rm=TRUE))
      bp <- boxplot(split(length_scale_g22,pexp), ylab=expression(paste("Length scale ", lambda, " [nm]")), 
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(length_scale_g22, na.rm=TRUE), y1 = 1.04*max(length_scale_g22, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(length_scale_g22, na.rm=TRUE) ,  get_asterisk(length_scale_g22.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_lengthscale_g22.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_lengthscale_g22.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_clusterradiusR_g22.pdf"), collapse = "_"))
      ylim <- c(0.9*min(clusterradiusR_g22, na.rm=TRUE) ,1.2*max(clusterradiusR_g22, na.rm=TRUE))
      bp <- boxplot(split(clusterradiusR_g22,pexp), ylab="Cluster radius R [nm]", 
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(clusterradiusR_g22, na.rm=TRUE), y1 = 1.04*max(clusterradiusR_g22, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(clusterradiusR_g22, na.rm=TRUE) ,  get_asterisk(clusterradiusR_g22.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_clusterradiusR_g22.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_clusterradiusR_g22.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_phi_g22.pdf"), collapse = "_"))
      ylim <- c(0.9*min(phi_g22, na.rm=TRUE),1.2*max(phi_g22, na.rm=TRUE))
      bp <- boxplot(split(phi_g22,pexp), ylab=expression(phi1^cluster),
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(phi_g22, na.rm=TRUE), y1 = 1.04*max(phi_g22, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(phi_g22, na.rm=TRUE) ,  get_asterisk(phi_g22.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_phi_g22.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_phi_g22.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_rho_g22.pdf"), collapse = "_"))
      ylim <- c(0.9*min(rho_g22, na.rm=TRUE),1.2*max(rho_g22, na.rm=TRUE))
      bp <- boxplot(split(rho_g22,pexp), ylab=expression(rho^cluster),
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(rho_g22, na.rm=TRUE), y1 = 1.04*max(rho_g22, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(rho_g22, na.rm=TRUE) ,  get_asterisk(rho_g22.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_rho_g22.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_rho_g22.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_Nclusters_g22.pdf"), collapse = "_"))
      ylim <- c(0.9*min(Nclusters_g22, na.rm=TRUE),1.2*max(Nclusters_g22, na.rm=TRUE))
      bp <- boxplot(split(Nclusters_g22,pexp), ylab=expression(paste(N^cluster, ' [points/cluster]')),
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(Nclusters_g22, na.rm=TRUE), y1 = 1.04*max(Nclusters_g22, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(Nclusters_g22, na.rm=TRUE) ,  get_asterisk(Nclusters_g22.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_Nclusters_g22.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_Nclusters_g22.pdf\' created.\n')
      
      if (fitting == 'expsq'){
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_kappa_g22.pdf"), collapse = "_"))
        ylim <- c(0.9*min(kappa_g22, na.rm=TRUE),1.2*max(kappa_g22, na.rm=TRUE))
        bp <- boxplot(split(kappa_g22,pexp), ylab=expression(paste(kappa, " [clusters/", nm^2, "]")), 
                      names=c(exp1_name, exp2_name), ylim=ylim)
        segments(x0 = 1, x1 = 2, y0 = 1.04*max(kappa_g22, na.rm=TRUE), y1 = 1.04*max(kappa_g22, na.rm=TRUE), col = "black")
        text( 1.5 , 1.09*max(kappa_g22, na.rm=TRUE) ,  get_asterisk(kappa_g22.test) , cex=1)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_kappa_g22.pdf"), collapse = "_"))
        cat('\t\'statistical_tests_kappa_g22.pdf\' created.\n')
      }
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_R2_g22.pdf"), collapse = "_"))
      # ylim <- c(10,1.2*max(parameters_g22$R2, na.rm=TRUE))
      bp <- boxplot(split(parameters_g22$R2,pexp), ylab=expression(paste(R^2,"(adj)")), names=c(exp1_name, exp2_name))
      # text(bp$group, bp$out, parameters_g22$id[which(parameters_g22$R2 %in% bp$out)], cex=0.5, pos = 4)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(parameters_g22$R2, na.rm=TRUE), y1 = 1.04*max(parameters_g22$R2, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(parameters_g22$R2, na.rm=TRUE) ,  get_asterisk(parameters_g22_R2.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_R2_g22.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_R2_g22.pdf\' created.\n')
      
      ## G11
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_r0_g11.pdf"), collapse = "_"))
      bp <- boxplot(split(parameters_g11$r0,pexp), ylab=expression(paste('Correlation range ', r[0], ' [nm]')), ylim=c(50,1.2*max(parameters_g11$r0, na.rm=TRUE)))
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(parameters_g11$r0, na.rm=TRUE), y1 = 1.04*max(parameters_g11$r0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(parameters_g11$r0, na.rm=TRUE), get_asterisk(parameters_g11_r0.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_r0_g11.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_r0_g11.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitude_g11.pdf"), collapse = "_"))
      ylim <- c(0,1.2*max(amplitude_g11, na.rm=TRUE))
      bp <- boxplot(split(amplitude_g11,pexp), ylab=expression("Amplitude A"), names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(amplitude_g11, na.rm=TRUE), y1 = 1.04*max(amplitude_g11, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(amplitude_g11, na.rm=TRUE),  get_asterisk(amplitude_g11.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitude_g11.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_amplitude_g11.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitudeAtr0_g11.pdf"), collapse = "_"))
      ylim <- c(0,1.2*max(amplitudeAtr0_g11, na.rm=TRUE))
      bp <- boxplot(split(amplitudeAtr0_g11,pexp), ylab=expression("Correlation at 20 nm [a.u.]"), names=c(exp1_name, exp2_name), ylim=ylim)
      # text(bp$group, bp$out, parameters_g11$id[which(amplitudeAtr0_g11 %in% bp$out)], cex=0.5, pos = 4)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(amplitudeAtr0_g11, na.rm=TRUE), y1 = 1.04*max(amplitudeAtr0_g11, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(amplitudeAtr0_g11, na.rm=TRUE),  get_asterisk(amplitudeAtr0_g11.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitudeAtr0_g11.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_amplitudeAtr0_g11.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_lengthscale_g11.pdf"), collapse = "_"))
      ylim <- c(10,1.2*max(length_scale_g11, na.rm=TRUE))
      bp <- boxplot(split(length_scale_g11,pexp), ylab=expression(paste("Length scale ", lambda, " [nm]")), 
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(length_scale_g11, na.rm=TRUE), y1 = 1.04*max(length_scale_g11, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(length_scale_g11, na.rm=TRUE) ,  get_asterisk(length_scale_g11.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_lengthscale_g11.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_lengthscale_g11.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_clusterradiusR_g11.pdf"), collapse = "_"))
      ylim <- c(0.9*min(clusterradiusR_g11, na.rm=TRUE),1.2*max(clusterradiusR_g11, na.rm=TRUE))
      bp <- boxplot(split(clusterradiusR_g11,pexp), ylab="Cluster radius R [nm]", 
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(clusterradiusR_g11, na.rm=TRUE), y1 = 1.04*max(clusterradiusR_g11, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(clusterradiusR_g11, na.rm=TRUE) ,  get_asterisk(clusterradiusR_g11.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_clusterradiusR_g11.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_clusterradiusR_g11.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_phi_g11.pdf"), collapse = "_"))
      ylim <- c(0.9*min(phi_g11, na.rm=TRUE),1.2*max(phi_g11, na.rm=TRUE))
      bp <- boxplot(split(phi_g11,pexp), ylab=expression(phi1^cluster),
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(phi_g11, na.rm=TRUE), y1 = 1.04*max(phi_g11, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(phi_g11, na.rm=TRUE) ,  get_asterisk(phi_g11.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_phi_g11.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_phi_g11.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_rho_g11.pdf"), collapse = "_"))
      ylim <- c(0.9*min(rho_g11, na.rm=TRUE),1.2*max(rho_g11, na.rm=TRUE))
      bp <- boxplot(split(rho_g11,pexp), ylab=expression(rho^cluster),
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(rho_g11, na.rm=TRUE), y1 = 1.04*max(rho_g11, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(rho_g11, na.rm=TRUE) ,  get_asterisk(rho_g11.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_rho_g11.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_rho_g11.pdf\' created.\n')
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_Nclusters_g11.pdf"), collapse = "_"))
      ylim <- c(0.9*min(Nclusters_g11, na.rm=TRUE),1.2*max(Nclusters_g11, na.rm=TRUE))
      bp <- boxplot(split(Nclusters_g11,pexp), ylab=expression(paste(N^cluster, ' [points/cluster]')),
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(Nclusters_g11, na.rm=TRUE), y1 = 1.04*max(Nclusters_g11, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(Nclusters_g11, na.rm=TRUE) ,  get_asterisk(Nclusters_g11.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_Nclusters_g11.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_Nclusters_g11.pdf\' created.\n')
      
      if (fitting == 'expsq'){
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_kappa_g11.pdf"), collapse = "_"))
        ylim <- c(0.9*min(kappa_g11, na.rm=TRUE),1.2*max(kappa_g11, na.rm=TRUE))
        bp <- boxplot(split(kappa_g11,pexp), ylab=expression(paste(kappa, " [clusters/", nm^2, "]")), 
                      names=c(exp1_name, exp2_name), ylim=ylim)
        segments(x0 = 1, x1 = 2, y0 = 1.04*max(kappa_g11, na.rm=TRUE), y1 = 1.04*max(kappa_g11, na.rm=TRUE), col = "black")
        text( 1.5 , 1.09*max(kappa_g11, na.rm=TRUE) ,  get_asterisk(kappa_g11.test) , cex=1)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_kappa_g11.pdf"), collapse = "_"))
        cat('\t\'statistical_tests_kappa_g11.pdf\' created.\n')
      }
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_R2_g11.pdf"), collapse = "_"))
      # ylim <- c(10,1.2*max(parameters_g11$R2, na.rm=TRUE))
      bp <- boxplot(split(parameters_g11$R2,pexp), ylab=expression(paste(R^2,"(adj)")), names=c(exp1_name, exp2_name))
      # text(bp$group, bp$out, parameters_g11$id[which(parameters_g11$R2 %in% bp$out)], cex=0.5, pos = 4)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(parameters_g11$R2, na.rm=TRUE), y1 = 1.04*max(parameters_g11$R2, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(parameters_g11$R2, na.rm=TRUE) ,  get_asterisk(parameters_g11_R2.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_R2_g11.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_R2_g11.pdf\' created.\n')
    }
    
    if(fitting %in% c('lin')){
      area <- 0.5*beta0*parameters_g12$r0
      
      # tests (experiments)
      cat('Running tests... ')
      pexp <- data$experiment[match(1:nlevels(data$id),data$id)]
      parameters_g12_r0.test <- t.test(parameters_g12$r0[which(parameters_g12$experiment == exp1_name)], parameters_g12$r0[which(parameters_g12$experiment == exp2_name)])
      beta0.test <- t.test(beta0[pexp==exp1_name],beta0[pexp==exp2_name])
      beta1.test <- t.test(beta1[pexp==exp1_name],beta1[pexp==exp2_name])
      area.test <- t.test(area[pexp==exp1_name],area[pexp==exp2_name])
      K_r0_g12.test <- t.test(K_r0_g12[which(parameters_g12$experiment == exp1_name)], K_r0_g12[which(parameters_g12$experiment == exp2_name)])
      cat('Done.\n')
      
      # # plot results
      cat('Plotting results of statistical tests...\n')
      openpdf("statistical_tests.pdf", pointsize=11)
      old <- par(mfrow=c(2, 2))
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_r0.pdf", collapse = "_"))
      ylim <- NULL
      bp <- boxplot(split(parameters_g12$r0,pexp), main=expression(r[0]), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.01*max(parameters_g12$r0, na.rm=TRUE), y1 = 1.01*max(parameters_g12$r0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.02*max(parameters_g12$r0, na.rm=TRUE), get_asterisk(parameters_g12_r0.test) , cex=1.2)
      # closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_r0.pdf", collapse = "_"))
      # cat('\t\'statistical_tests_r0.pdf\' created.\n')
      
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_beta0.pdf", collapse = "_"))
      ylim <- NULL
      bp <- boxplot(split(beta0,pexp), main="beta0", names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.01*max(beta0, na.rm=TRUE), y1 = 1.01*max(beta0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.02*max(beta0, na.rm=TRUE),  get_asterisk(beta0.test) , cex=1.2)
      # closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_beta0.pdf", collapse = "_"))
      # cat('\t\'statistical_tests_beta0.pdf\' created.\n')
      
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_beta1.pdf", collapse = "_"))
      ylim <- NULL
      bp <- boxplot(split(beta1,pexp), main="beta1", 
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.01*max(beta1, na.rm=TRUE), y1 = 1.01*max(beta1, na.rm=TRUE), col = "black")
      text( 1.5 , 1.02*max(beta1, na.rm=TRUE) ,  get_asterisk(beta1.test) , cex=1.2)
      # closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_beta1.pdf", collapse = "_"))
      # cat('\t\'statistical_tests_beta1.pdf\' created.\n')
      
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_area.pdf", collapse = "_"))
      ylim <- NULL
      bp <- boxplot(split(area,pexp), main=expression("area"), names=c(exp1_name, exp2_name),  ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.01*max(area, na.rm=TRUE), y1 = 1.01*max(area, na.rm=TRUE), col = "black")
      text( 1.5 , 1.02*max(area, na.rm=TRUE) , get_asterisk(area.test), cex=1.2)
      # closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_area.pdf", collapse = "_"))
      # cat('\t\'statistical_tests_area.pdf\' created.\n')
      par(old); closepdf("statistical_tests.pdf")
      cat('\t\'statistical_tests.pdf\' created.\n')
      
      openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_Kr0.pdf", collapse = "_"))
      ylim <- NULL
      bp <- boxplot(split(K_r0,pexp), main=expression(K(r[0])), names=c(exp1_name, exp2_name),  ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.01*max(K_r0, na.rm=TRUE), y1 = 1.01*max(K_r0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.02*max(K_r0, na.rm=TRUE) , get_asterisk(K_r0.test), cex=1.2)
      closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_Kr0.pdf", collapse = "_"))
    }
    if (save){
      parameters_g11_aux <- parameters_g11; parameters_g22_aux <- parameters_g22; parameters_g12_aux <- parameters_g12
      load(paste(c(path_to_experiments[1], RData_2experiments), collapse='/'))
      parameters_g11 <- parameters_g11_aux; parameters_g22 <- parameters_g22_aux; parameters_g12 <- parameters_g12_aux
      save(data, data_allrange, units, parameters_g11, parameters_g22, parameters_g12, density, npoints,
           r_eval, levels1, levels2, 
           file=paste(c(paste(c(path_to_experiments[1], paste(c(levels1,levels2), collapse = '')), collapse='/'), 
                        '_', fitting, '_2experiments.RData'), collapse = ''))
    }
  }
  if (pptype=="unmarked"){
    ## exponential fitting
    if(fitting %in% c('exp','expsq')){
      data$g11 <- log(data$g11-1)
      attach(data)
      cat('Exponential transformation... Done\n')
    }
    ## response feature analysis (Faraway, 2006, p.206)
    # fit ordinary regressions separately to each group
    parameters_g11 <- data.frame(parameters_g11, 
                                 beta0=matrix(nrow=nlevels(data$id),ncol=1), 
                                 beta1=matrix(nrow=nlevels(data$id),ncol=1), 
                                 R2=matrix(nrow=nlevels(data$id),ncol=1))  # initialize
    K_r0_g11 <- numeric(nlevels(data$id))
    
    for(i in 1:nlevels(data$id)){
      if (fitting %in% c("lin", "exp")){lmod <- lm(g11 ~ r, subset=(id==i), data=data)}
      if (fitting == "expsq"){lmod <- lm(g11 ~ 1 + I(r^2), subset=(id==i), data=data)}
      parameters_g11$R2[i] <- summary(lmod)$adj.r.squared
      parameters_g11$beta0[i] <- coef(lmod)[1]
      parameters_g11$beta1[i] <- coef(lmod)[2]
      a <- tail(which(r_eval < 0.5*parameters_g11$r0[which(parameters_g11$id == i)]),n=1)
      if (length(a)==0){ K_r0_g11[i] <- NA}
      else{
        K_r0_g11[i] <- K11[id == i][a]
      }
    }
    cat('Checking quality of fitting... ')
    fileName <- paste(c("../output/", unlist(strsplit(RData_2experiments, "\\_"))[1], "_summaryfitting_", fitting, '.txt'), collapse = "")
    if (file.exists(fileName)) file.remove(fileName)    
    if (fitting %in% c("exp","lin")){
      lmod <- lm(g11 ~ r, subset=(experiment==levels(data$experiment)[1]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g11_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g11_plots.pdf'), collapse=""))
      if (!is.empty(exp2_name)){
        lmod <- lm(g11 ~ r, subset=(experiment==levels(data$experiment)[2]), data = data)
        capture.output(summary(lmod), file=fileName, append=TRUE)
      }
      cat('Done.\n')
    }
    if (fitting=="expsq"){
      lmod <- lm(g11 ~ 1+I(r^2), subset=(experiment==levels(data$experiment)[1]), data = data)
      capture.output(summary(lmod), file=fileName, append=TRUE)
      openpdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g11_plots.pdf'), collapse=""))
      plot(lmod)
      closepdf(paste(c(strsplit(fileName, "\\.txt"), '_', levels(data$experiment)[1], '_g11_plots.pdf'), collapse=""))
      if (!is.empty(exp2_name)){
        lmod <- lm(g11 ~ 1+I(r^2), subset=(experiment==levels(data$experiment)[2]), data = data)
        capture.output(summary(lmod), file=fileName, append=TRUE)
      }
      cat('Done.\n')
    }
    if(fitting %in% c('exp', 'expsq')){
      amplitude_g11 <- exp(parameters_g11$beta0)
      length_scale_g11 <- -parameters_g11$beta1^(-1)
      amplitudeAtr0_g11 <- exp(parameters_g11$beta0*parameters_g11$beta1*reso_r)+1
      area_g11 <- -length_scale_g11*amplitude_g11*(exp(-1/(length_scale_g11)*parameters_g11$r0)-1)
      
      if (fitting == 'expsq'){
        clusterradiusR_g11 <- sqrt(length_scale_g11)
        kappa_g11 <- 1/(amplitude_g11*pi*length_scale_g11)  # number of clusters per area
        phi_g11 <- amplitude_g11/4   # rho_cluster/rho_average
        # Nclusters_g11 <- <- amplitude_g11*pi*length_scale_g11*density$level1
        Nclusters_g11 <- density$level1/kappa_g11  #  average number of points per cluster
      }      
      if (fitting == 'exp'){
        clusterradiusR_g11 <- length_scale_g11
        phi_g11 <- 2*amplitude_g11   # rho_cluster/rho_average
        Nclusters_g11 <- 2*amplitude_g11*pi*length_scale_g11^2*density$level1  # average number of points per cluster
      }      
      
      # tests (experiments)
      cat('Running tests... \n')
      pexp <- data$experiment[match(1:nlevels(data$id),data$id)]
      
      col_palette <-  c(rgb(173,216,230,max = 255,alpha=125), rgb(255,165,0,max = 255,alpha=125)); lty=c(1,2)
      r_eval_plot <- which(r_eval>=reso_r)[1]:(0.8*length(r_eval))
      xlim <- c(r_eval[r_eval_plot][1],r_eval[r_eval_plot][length(r_eval[r_eval_plot])])
      ylim=c(0.5,10)
      openpdf(paste(c("g11", filename_experiments, "_", fitting, '.pdf'), collapse=""))
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''),
           ylab=eval(bquote(expression(g[.(levels1)][','][.(levels1)](r)))), type='n', ylim = ylim)
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        subset <- subset(data_allrange, data$experiment %in% exp_name)
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, g11, lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
        abline(h=1, lty=2, col='black')
        legend("topright", legend=exp_names, lty=lty, col=col_palette)
      }
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        medianr0 <- median(parameters_g11$r0[which(parameters_g11$experiment == exp_name)], na.rm=TRUE)
        medianlengthscale <- median(length_scale_g11[pexp==exp_name], na.rm=TRUE)
        medianamplitude <- median(amplitude_g11[pexp==exp_name], na.rm=TRUE)
        if (fitting=="exp"){fitted_line <- 1+medianamplitude*exp(-seq(0, medianr0, by=0.5)/medianlengthscale)}
        if (fitting=="expsq"){fitted_line <- 1+medianamplitude*exp(-seq(0, medianr0, by=0.5)^2/medianlengthscale)}
        if (fitting=="lin"){fitted_line <- medianamplitude - 1/medianlengthscale*seq(0, medianr0, by=0.5)}
        lines(seq(0, medianr0, by=0.5),fitted_line , col='dimgray', lty=1)
      }
      closepdf(paste(c("g11", filename_experiments, "_", fitting, '.pdf'), collapse=""))
      cat('\t\'g11_2experiments.pdf\' created.\n')
      openpdf(paste(c("g11_log", filename_experiments, "_", fitting, '.pdf'), collapse=""))
      plot(1, xlim=c(0,median(parameters_g11$r0, na.rm=TRUE)), 
           ylim=range(data$g11, na.rm=TRUE), xlab=paste('r ','[', units,']', sep = ''),
           ylab=eval(bquote(expression(ln(g[.(levels1)][','][.(levels1)](r)-1)))), type='n')
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        subset <- subset(data, data$experiment %in% exp_name)
        for (jj in unique(subset$id)){
          with(subset[which(subset$id == jj), ], lines(r, g11, lty=lty[ii], main='', pch= 20, lwd=1, col=col_palette[ii]))
        }
        legend("topright", legend=exp_names, lty=lty, col=col_palette)
      }
      for (ii in  1:length(levels(data$experiment))){
        exp_name <- levels(data$experiment)[ii]
        medianr0 <- median(parameters_g11$r0[which(parameters_g11$experiment == exp_name)], na.rm=TRUE)
        medianlengthscale <- median(length_scale_g11[pexp==exp_name], na.rm=TRUE)
        medianamplitude <- median(amplitude_g11[pexp==exp_name], na.rm=TRUE)
        if (fitting=="exp"){fitted_line <- log(medianamplitude)-seq(0, medianr0, by=0.5)/medianlengthscale}
        if (fitting=="expsq"){fitted_line <- log(medianamplitude)-seq(0, medianr0, by=0.5)^2/medianlengthscale}
        if (fitting=="lin"){fitted_line <- medianamplitude - 1/medianlengthscale*seq(0, medianr0, by=0.5)}
        lines(seq(0, medianr0, by=0.5),fitted_line , col='dimgray', lty=1)
      }
      closepdf(paste(c("g11_log", filename_experiments, "_", fitting, '.pdf'), collapse=""))
      cat('\t\'g11_log_2experiments(FITTING).pdf\' created.\n')
      # openpdf(paste(c('cdf11_npoints', filename_experiments, '.pdf'), collapse=''))
      # col_palette <-  c(rgb(173,216,230,max = 255), rgb(255,165,0,max = 255))
      # plot(1, xlim=c(0,2*min(parameters_g11$r0, na.rm=TRUE)), ylim=c(0,100), xlab=paste('2r ','[', units,']', sep = ''),
      #      ylab=paste('Av. number of points (', levels1, ',', levels1, ') [%]', sep=''), type='n')
      # for (ii in  1:length(levels(data$experiment))){
      #   exp_name <- levels(data$experiment)[ii]
      #   medianclusterradiusR <- median(clusterradiusR_g11[pexp==exp_name], na.rm=TRUE)
      #   if (fitting=="expsq"){
      #     r <- seq(0,2*min(parameters_g11$r0[pexp==exp_name], na.rm=TRUE), by=5)
      #     cdf_line <- 100*(pnorm(r, mean=0, sd=medianclusterradiusR/4)-pnorm(-r, mean=0, sd=medianclusterradiusR/4))
      #   }
      #   lines(2*r, cdf_line , col=col_palette[ii], lty=1)
      #   legend("bottomright", legend=exp_names, lty=1, col=col_palette)
      # }
      # closepdf(paste(c('cdf11_npoints', filename_experiments, '.pdf'), collapse=''))
      # cat('\t\'cdf11_npoints_2experiments(FITTING).pdf\' created.\n')
      
      if (!is.empty(exp2_name)){
        # tests (experiments)
        cat('Running tests... ')
        parameters_g11_r0.test <- t.test(parameters_g11$r0[which(parameters_g11$experiment == exp1_name)], parameters_g11$r0[which(parameters_g11$experiment == exp2_name)])
        amplitude_g11.test <- t.test(amplitude_g11[pexp==exp1_name],amplitude_g11[pexp==exp2_name])
        amplitudeAtr0_g11.test <- t.test(amplitudeAtr0_g11[pexp==exp1_name],amplitudeAtr0_g11[pexp==exp2_name])
        length_scale_g11.test <- t.test(length_scale_g11[pexp==exp1_name],length_scale_g11[pexp==exp2_name])
        clusterradiusR_g11.test <- t.test(clusterradiusR_g11[pexp==exp1_name],clusterradiusR_g11[pexp==exp2_name])
        parameters_g11_R2.test <- t.test(parameters_g11$R2[pexp==exp1_name],parameters_g11$R2[pexp==exp2_name])
        area_g11.test <- t.test(area_g11[pexp==exp1_name],area_g11[pexp==exp2_name])
        density_level1.test <- t.test(density$level1[pexp==exp1_name],density$level1[pexp==exp2_name])
        length_scale_g11.test <- t.test(length_scale_g11[pexp==exp1_name],length_scale_g11[pexp==exp2_name])
        phi_g11.test <- t.test(phi_g11[pexp==exp1_name],phi_g11[pexp==exp2_name])
        Nclusters_g11.test <- t.test(Nclusters_g11[pexp==exp1_name],Nclusters_g11[pexp==exp2_name])
        if (fitting == 'expsq'){
          kappa_g11.test <- t.test(kappa_g11[pexp==exp1_name],kappa_g11[pexp==exp2_name])
        }
        cat('Done.\n')
        
        # # plot results
        cat('Plotting results of statistical tests...\n')
        # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests.pdf", collapse = "_"), pointsize=11)
        # old <- par(mfrow=c(2, 2))
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_density.pdf"), collapse = "_"))
        ylim <- c(0,1.2*max(density$level1, na.rm=TRUE))
        bp <- boxplot(split(density$level1, pexp),
                      names=c(exp1_name, exp2_name), ylim=ylim, ylab=bquote('Average density ' ~ .(levels1) ~ ' [points/' ~ nm^2 ~']'))
        # text(bp$group, bp$out, density$id[which(density$level1 == bp$out, arr.ind=TRUE)], cex=0.5, pos = 4)
        segments(x0 = 1, x1 = 2, y0 = 1.04*max(density$level1, na.rm=TRUE), y1 = 1.04*max(density$level1, na.rm=TRUE), col = "black")
        text( 1.5 , 1.09*max(density$level1, na.rm=TRUE),  get_asterisk(density_level1.test) , cex=1)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_density.pdf"), collapse = "_"))
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_r0_g11.pdf"), collapse = "_"))
        bp <- boxplot(split(parameters_g11$r0,pexp), ylab=expression(paste('Correlation range ', r[0], ' [nm]')), ylim=c(50,1.2*max(parameters_g11$r0, na.rm=TRUE)))
        # text(bp$group, bp$out, parameters_g11$id[which(parameters_g11$r0 %in% bp$out)], cex=0.5, pos = 4)
        segments(x0 = 1, x1 = 2, y0 = 1.04*max(parameters_g11$r0, na.rm=TRUE), y1 = 1.04*max(parameters_g11$r0, na.rm=TRUE), col = "black")
        text( 1.5 , 1.09*max(parameters_g11$r0, na.rm=TRUE), get_asterisk(parameters_g11_r0.test) , cex=1)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_r0_g11.pdf"), collapse = "_"))
        cat('\t\'statistical_tests_r0_g11.pdf\' created.\n')
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitudeAtr0_g11.pdf"), collapse = "_"))
        ylim <- c(0,1.2*max(amplitudeAtr0_g11, na.rm=TRUE))
        bp <- boxplot(split(amplitudeAtr0_g11,pexp), ylab=expression("Correlation at 20 nm [a.u.]"), names=c(exp1_name, exp2_name), ylim=ylim)
        # text(bp$group, bp$out, parameters_g11$id[which(amplitudeAtr0_g11 %in% bp$out)], cex=0.5, pos = 4)
        segments(x0 = 1, x1 = 2, y0 = 1.04*max(amplitudeAtr0_g11, na.rm=TRUE), y1 = 1.04*max(amplitudeAtr0_g11, na.rm=TRUE), col = "black")
        text( 1.5 , 1.09*max(amplitudeAtr0_g11, na.rm=TRUE),  get_asterisk(amplitudeAtr0_g11.test) , cex=1)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitudeAtr0_g11.pdf"), collapse = "_"))
        cat('\t\'statistical_tests_amplitudeAtr0_g11.pdf\' created.\n')
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitude_g11.pdf"), collapse = "_"))
        ylim <- c(0,1.2*max(amplitude_g11, na.rm=TRUE))
        bp <- boxplot(split(amplitude_g11,pexp), ylab="Amplitude A", names=c(exp1_name, exp2_name), ylim=ylim)
        # text(bp$group, bp$out, parameters_g11$id[which(amplitude_g11 %in% bp$out)], cex=0.5, pos = 4)
        segments(x0 = 1, x1 = 2, y0 = 1.04*max(amplitude_g11, na.rm=TRUE), y1 = 1.04*max(amplitude_g11, na.rm=TRUE), col = "black")
        text( 1.5 , 1.09*max(amplitude_g11, na.rm=TRUE),  get_asterisk(amplitude_g11.test) , cex=1)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_amplitude_g11.pdf"), collapse = "_"))
        cat('\t\'statistical_tests_amplitude_g11.pdf\' created.\n')
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_lengthscale_g11.pdf"), collapse = "_"))
        ylim <- c(10,1.2*max(length_scale_g11, na.rm=TRUE))
        bp <- boxplot(split(length_scale_g11,pexp), ylab=expression(paste("Length scale ", lambda, " [nm]")), 
                      names=c(exp1_name, exp2_name), ylim=ylim)
        segments(x0 = 1, x1 = 2, y0 = 1.04*max(length_scale_g11, na.rm=TRUE), y1 = 1.04*max(length_scale_g11, na.rm=TRUE), col = "black")
        text( 1.5 , 1.09*max(length_scale_g11, na.rm=TRUE) ,  get_asterisk(length_scale_g11.test) , cex=1)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_lengthscale_g11.pdf"), collapse = "_"))
        cat('\t\'statistical_tests_lengthscale_g11.pdf\' created.\n')
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_clusterradiusR_g11.pdf"), collapse = "_"))
        ylim <- c(0.9*min(clusterradiusR_g11, na.rm=TRUE),1.2*max(clusterradiusR_g11, na.rm=TRUE))
        bp <- boxplot(split(clusterradiusR_g11,pexp), ylab="cluster radius R [nm]", 
                      names=c(exp1_name, exp2_name), ylim=ylim)
        segments(x0 = 1, x1 = 2, y0 = 1.04*max(clusterradiusR_g11, na.rm=TRUE), y1 = 1.04*max(clusterradiusR_g11, na.rm=TRUE), col = "black")
        text( 1.5 , 1.09*max(clusterradiusR_g11, na.rm=TRUE) ,  get_asterisk(clusterradiusR_g11.test) , cex=1)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_clusterradiusR_g11.pdf"), collapse = "_"))
        cat('\t\'statistical_tests_clusterradiusR_g11.pdf\' created.\n')
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_phi_g11.pdf"), collapse = "_"))
        ylim <- c(0.9*min(phi_g11, na.rm=TRUE),1.2*max(phi_g11, na.rm=TRUE))
        bp <- boxplot(split(phi_g11,pexp), ylab=expression(phi1^cluster),
                      names=c(exp1_name, exp2_name), ylim=ylim)
        segments(x0 = 1, x1 = 2, y0 = 1.04*max(phi_g11, na.rm=TRUE), y1 = 1.04*max(phi_g11, na.rm=TRUE), col = "black")
        text( 1.5 , 1.09*max(phi_g11, na.rm=TRUE) ,  get_asterisk(phi_g11.test) , cex=1)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_phi_g11.pdf"), collapse = "_"))
        cat('\t\'statistical_tests_phi_g11.pdf\' created.\n')
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_Nclusters_g11.pdf"), collapse = "_"))
        ylim <- c(0.9*min(Nclusters_g11, na.rm=TRUE),1.2*max(Nclusters_g11, na.rm=TRUE))
        bp <- boxplot(split(Nclusters_g11,pexp), ylab=expression(paste(N^cluster, ' [points/cluster]')),
                      names=c(exp1_name, exp2_name), ylim=ylim)
        segments(x0 = 1, x1 = 2, y0 = 1.04*max(Nclusters_g11, na.rm=TRUE), y1 = 1.04*max(Nclusters_g11, na.rm=TRUE), col = "black")
        text( 1.5 , 1.09*max(Nclusters_g11, na.rm=TRUE) ,  get_asterisk(Nclusters_g11.test) , cex=1)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_Nclusters_g11.pdf"), collapse = "_"))
        cat('\t\'statistical_tests_Nclusters_g11.pdf\' created.\n')
        
        if (fitting == 'expsq'){
          openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_kappa_g11.pdf"), collapse = "_"))
          ylim <- c(0.9*min(kappa_g11, na.rm=TRUE),1.2*max(kappa_g11, na.rm=TRUE))
          bp <- boxplot(split(kappa_g11,pexp), ylab=expression(paste(kappa, " [clusters/", nm^2, "]")), 
                        names=c(exp1_name, exp2_name), ylim=ylim)
          segments(x0 = 1, x1 = 2, y0 = 1.04*max(kappa_g11, na.rm=TRUE), y1 = 1.04*max(kappa_g11, na.rm=TRUE), col = "black")
          text( 1.5 , 1.09*max(kappa_g11, na.rm=TRUE) ,  get_asterisk(kappa_g11.test) , cex=1)
          closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_kappa_g11.pdf"), collapse = "_"))
          cat('\t\'statistical_tests_kappa_g11.pdf\' created.\n')
        }
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_R2_g11.pdf"), collapse = "_"))
        # ylim <- c(10,1.2*max(parameters_g11$R2, na.rm=TRUE))
        bp <- boxplot(split(parameters_g11$R2,pexp), ylab=expression(paste(R^2,"(adj)")), names=c(exp1_name, exp2_name))
        # text(bp$group, bp$out, parameters_g11$id[which(parameters_g11$R2 %in% bp$out)], cex=0.5, pos = 4)
        segments(x0 = 1, x1 = 2, y0 = 1.04*max(parameters_g11$R2, na.rm=TRUE), y1 = 1.04*max(parameters_g11$R2, na.rm=TRUE), col = "black")
        text( 1.5 , 1.09*max(parameters_g11$R2, na.rm=TRUE) ,  get_asterisk(parameters_g11_R2.test) , cex=1)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests_R2_g11.pdf"), collapse = "_"))
        cat('\t\'statistical_tests_R2_g11.pdf\' created.\n')
        
      } else{
        parameters_g11$r0[which(parameters_g11$experiment == exp1_name)]
        amplitude_g11[pexp==exp1_name]
        amplitudeAtr0_g11[pexp==exp1_name]
        length_scale_g11[pexp==exp1_name]
        clusterradiusR_g11[pexp==exp1_name]
        parameters_g11$R2[pexp==exp1_name]
        area_g11[pexp==exp1_name]
        density$level1[pexp==exp1_name]
        length_scale_g11[pexp==exp1_name]
        phi_g11[pexp==exp1_name]
        Nclusters_g11[pexp==exp1_name]
        if (fitting == 'expsq'){
          kappa_g11[pexp==exp1_name]
        }
        
        # # plot results
        cat('Plotting results...\n')
        # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests.pdf", collapse = "_"), pointsize=11)
        # old <- par(mfrow=c(2, 2))
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "density.pdf"), collapse = "_"))
        ylim <- c(0,1.2*max(density$level1, na.rm=TRUE))
        bp <- boxplot(split(density$level1, pexp),
                      names=c(exp1_name, exp2_name), ylim=ylim, ylab=bquote('Average density ' ~ .(levels1) ~ ' [points/' ~ nm^2 ~']'))
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "density.pdf"), collapse = "_"))
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "r0_g11.pdf"), collapse = "_"))
        bp <- boxplot(split(parameters_g11$r0,pexp), ylab=expression(paste('Correlation range ', r[0], ' [nm]')), ylim=c(50,1.2*max(parameters_g11$r0, na.rm=TRUE)))
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "r0_g11.pdf"), collapse = "_"))
        cat('\t\'r0_g11.pdf\' created.\n')
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "amplitudeAtr0_g11.pdf"), collapse = "_"))
        ylim <- c(0,1.2*max(amplitudeAtr0_g11, na.rm=TRUE))
        bp <- boxplot(split(amplitudeAtr0_g11,pexp), ylab=expression("Correlation at 20 nm [a.u.]"), names=c(exp1_name, exp2_name), ylim=ylim)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "amplitudeAtr0_g11.pdf"), collapse = "_"))
        cat('\t\'amplitudeAtr0_g11.pdf\' created.\n')
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "amplitude_g11.pdf"), collapse = "_"))
        ylim <- c(0,1.2*max(amplitude_g11, na.rm=TRUE))
        bp <- boxplot(split(amplitude_g11,pexp), ylab="Amplitude A", names=c(exp1_name, exp2_name), ylim=ylim)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "amplitude_g11.pdf"), collapse = "_"))
        cat('\t\'statistical_tests_amplitude_g11.pdf\' created.\n')
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "lengthscale_g11.pdf"), collapse = "_"))
        ylim <- c(10,1.2*max(length_scale_g11, na.rm=TRUE))
        bp <- boxplot(split(length_scale_g11,pexp), ylab=expression(paste("Length scale ", lambda, " [nm]")), 
                      names=c(exp1_name, exp2_name), ylim=ylim)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "lengthscale_g11.pdf"), collapse = "_"))
        cat('\t\'lengthscale_g11.pdf\' created.\n')
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "clusterradiusR_g11.pdf"), collapse = "_"))
        ylim <- c(0.9*min(clusterradiusR_g11, na.rm=TRUE),1.2*max(clusterradiusR_g11, na.rm=TRUE))
        bp <- boxplot(split(clusterradiusR_g11,pexp), ylab="cluster radius R [nm]", 
                      names=c(exp1_name, exp2_name), ylim=ylim)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "clusterradiusR_g11.pdf"), collapse = "_"))
        cat('\t\'clusterradiusR_g11.pdf\' created.\n')
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "phi_g11.pdf"), collapse = "_"))
        ylim <- c(0.9*min(phi_g11, na.rm=TRUE),1.2*max(phi_g11, na.rm=TRUE))
        bp <- boxplot(split(phi_g11,pexp), ylab=expression(phi1^cluster),
                      names=c(exp1_name, exp2_name), ylim=ylim)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "phi_g11.pdf"), collapse = "_"))
        cat('\t\'phi_g11.pdf\' created.\n')
        
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "Nclusters_g11.pdf"), collapse = "_"))
        ylim <- c(0.9*min(Nclusters_g11, na.rm=TRUE),1.2*max(Nclusters_g11, na.rm=TRUE))
        bp <- boxplot(split(Nclusters_g11,pexp), ylab=expression(paste(N^cluster, ' [points/cluster]')),
                      names=c(exp1_name, exp2_name), ylim=ylim)
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "Nclusters_g11.pdf"), collapse = "_"))
        cat('\t\'Nclusters_g11.pdf\' created.\n')
        
        if (fitting == 'expsq'){
          openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "kappa_g11.pdf"), collapse = "_"))
          ylim <- c(0.9*min(kappa_g11, na.rm=TRUE),1.2*max(kappa_g11, na.rm=TRUE))
          bp <- boxplot(split(kappa_g11,pexp), ylab=expression(paste(kappa, " [clusters/", nm^2, "]")), 
                        names=c(exp1_name, exp2_name), ylim=ylim)
          closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "kappa_g11.pdf"), collapse = "_"))
          cat('\t\'kappa_g11.pdf\' created.\n')
        }
        openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "R2_g11.pdf"), collapse = "_"))
        # ylim <- c(10,1.2*max(parameters_g11$R2, na.rm=TRUE))
        bp <- boxplot(split(parameters_g11$R2,pexp), ylab=expression(paste(R^2,"(adj)")), names=c(exp1_name, exp2_name))
        closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "R2_g11.pdf"), collapse = "_"))
        cat('\t\'R2_g11.pdf\' created.\n')
      }
    }
    
    if(fitting %in% c('lin')){
      area <- 0.5*beta0*parameters_g11$r0
      
      # tests (experiments)
      cat('Running tests... ')
      pexp <- data$experiment[match(1:nlevels(data$id),data$id)]
      parameters_g11:r0.test <- t.test(parameters_g11$r0[which(parameters_g11$experiment == exp1_name)], parameters_g11$r0[which(parameters_g11$experiment == exp2_name)])
      beta0.test <- t.test(beta0[pexp==exp1_name],beta0[pexp==exp2_name])
      beta1.test <- t.test(beta1[pexp==exp1_name],beta1[pexp==exp2_name])
      area.test <- t.test(area[pexp==exp1_name],area[pexp==exp2_name])
      K_r0.test <- t.test(K_r0[which(parameters_g11$experiment == exp1_name)], K_r0[which(parameters_g11$experiment == exp2_name)])
      cat('Done.\n')
      
      # # plot results
      cat('Plotting results of statistical tests...\n')
      openpdf("statistical_tests.pdf", pointsize=11)
      old <- par(mfrow=c(2, 2))
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_r0.pdf", collapse = "_"))
      ylim <- NULL
      bp <- boxplot(split(parameters_g11$r0,pexp), main=expression(r[0]), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.01*max(parameters_g11$r0, na.rm=TRUE), y1 = 1.01*max(parameters_g11$r0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.02*max(parameters_g11$r0, na.rm=TRUE), get_asterisk(parameters_g11_r0.test) , cex=1.2)
      # closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_r0.pdf", collapse = "_"))
      # cat('\t\'statistical_tests_r0.pdf\' created.\n')
      
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_beta0.pdf", collapse = "_"))
      ylim <- NULL
      bp <- boxplot(split(beta0,pexp), main="beta0", names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.01*max(beta0, na.rm=TRUE), y1 = 1.01*max(beta0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.02*max(beta0, na.rm=TRUE),  get_asterisk(beta0.test) , cex=1.2)
      # closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_beta0.pdf", collapse = "_"))
      # cat('\t\'statistical_tests_beta0.pdf\' created.\n')
      
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_beta1.pdf", collapse = "_"))
      ylim <- NULL
      bp <- boxplot(split(beta1,pexp), main="beta1", 
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.01*max(beta1, na.rm=TRUE), y1 = 1.01*max(beta1, na.rm=TRUE), col = "black")
      text( 1.5 , 1.02*max(beta1, na.rm=TRUE) ,  get_asterisk(beta1.test) , cex=1.2)
      # closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_beta1.pdf", collapse = "_"))
      # cat('\t\'statistical_tests_beta1.pdf\' created.\n')
      
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_area.pdf", collapse = "_"))
      ylim <- NULL
      bp <- boxplot(split(area,pexp), main=expression("area"), names=c(exp1_name, exp2_name),  ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.01*max(area, na.rm=TRUE), y1 = 1.01*max(area, na.rm=TRUE), col = "black")
      text( 1.5 , 1.02*max(area, na.rm=TRUE) , get_asterisk(area.test), cex=1.2)
      # closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_area.pdf", collapse = "_"))
      # cat('\t\'statistical_tests_area.pdf\' created.\n')
      par(old); closepdf("statistical_tests.pdf")
      cat('\t\'statistical_tests.pdf\' created.\n')
      
      openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_Kr0.pdf", collapse = "_"))
      ylim <- NULL
      bp <- boxplot(split(K_r0,pexp), main=expression(K(r[0])), names=c(exp1_name, exp2_name),  ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.01*max(K_r0, na.rm=TRUE), y1 = 1.01*max(K_r0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.02*max(K_r0, na.rm=TRUE) , get_asterisk(K_r0.test), cex=1.2)
      closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_Kr0.pdf", collapse = "_"))
    }
    
    if (save){
      parameters_g11_aux <- parameters_g11
      load(paste(c(path_to_experiments[1], RData_2experiments), collapse='/'))
      parameters_g11 <- parameters_g11_aux
      save(data, data_allrange, units, parameters_g11, density, npoints, r_eval, levels1, levels2, 
           file=paste(c(paste(c(path_to_experiments[1], paste(c(levels1,levels2), collapse = '')), collapse='/'), 
                        '_', fitting, filename_experiments, '.RData'), collapse = ''))
    }
  }
}

get_asterisk <- function(test, low_sig=0.05, mid_sig=0.01, high_sig=0.001){
  if (test$p.value > low_sig){
    return("")
  }
  else {
    if (test$p.value <= high_sig){
      return("***")
    }
    if (test$p.value <= mid_sig){
      return("**")
    }
    else{
      return("*")
    }
  }
}

rotate <- function(alpha){
  alpha <- alpha*pi/180
  R <- matrix(c(cos(alpha), -sin(alpha), sin(alpha), cos(alpha)), # the data elements 
              nrow=2,              # number of rows 
              ncol=2,              # number of columns 
              byrow = TRUE)  
  return(R)
}

savetext <- function(pp, file, dir="../output/", ncolumns=2)
{
  dir.create(file.path(dir), showWarnings = FALSE)
  fname <- paste(dir,file, sep="")
  write(t(pp),fname,ncolumns)
}

triangle <- function(cenx, ceny, rot, height, gamma)
{
  Ax <- - height*tan(gamma*pi/180); Ay <- - height/3
  Bx <- height*tan(gamma*pi/180); By <- - height/3
  Cx <- 0; Cy <- height - height/3
  vertex <- list(x=c(Ax,Bx,Cx),y=c(Ay, By, Cy))
  vertex_rot <- rotate(alpha=rot)%*%t(matrix(c(vertex$x,vertex$y), nrow=3, ncol=2))
  vertex_rot_trans <- vertex_rot + c(cenx, ceny) 
  shape <- owin(poly=list(x=vertex_rot_trans[1,], y=vertex_rot_trans[2,]))
  return(shape)  
}

quadrilateral <- function(cenx, ceny, rot, height, width)
{
  Ax <- - width/2; Ay <- - height/2
  Bx <- width/2; By <- - height/2
  Cx <- width/2; Cy <- height/2
  Dx <- - width/2; Dy <- height/2
  vertex <- list(x=c(Ax,Bx,Cx,Dx),y=c(Ay, By, Cy, Dy))
  vertex_rot <- rotate(alpha=rot)%*%t(matrix(c(vertex$x,vertex$y), nrow=4, ncol=2))
  vertex_rot_trans <- vertex_rot + c(cenx,ceny)
  shape <- owin(poly=list(x=vertex_rot_trans[1,], y=vertex_rot_trans[2,]))
  return(shape)  
}

circle <- function(distribution, n=Nclusters, radius=radius, centre=centre, sig=sig){
  
  if (distribution=='uniform'){
    points <- runifdisc(n=n, radius=radius, centre=centre) 
    return(points)
  }
  if (distribution=='matern'){
    # npoints <- rpoispp(n)
    npoints <- rpois(1, lambda=Nclusters)
    points <- runifdisc(n=npoints, radius=radius, centre=centre) 
    return(points)
  }
  if (distribution=='thomas'){
    npoints <- rpois(1, lambda=Nclusters)
    p <- mvrnorm(n=npoints, mu=c(0,0), Sigma=matrix(c(sig^2,0,0,sig^2), 2, 2))
    if (npoints==1){
      p[1] <- p[1] + centre[1]
      p[2] <- p[2] + centre[2]
    }
    else{
      p[,1] <- p[,1] + centre[1]
      p[,2] <- p[,2] + centre[2]
    }
    points <- as.ppp(p, c(centre[1]-radius, centre[1]+radius, centre[2]-radius, centre[2]+radius))
    return(points)
  }
  if (distribution=='gaussian_truncated'){
    p <- mvrnorm(n=npoints, mu=c(0,0), Sigma=matrix(c(sig^2,0,0,sig^2), 2, 2))
    out <- which((abs(p[,1])>radius) | (abs(p[,2])>radius))
    if (length(out)>0){
      p_constrained <- p[-out,] 
    }
    else{ p_constrained <- p }
    if (n < dim(p_constrained)[1]){ p_reduced <- p_constrained[1:n,]}
    else{ p_reduced <- p_constrained}
    p_reduced[,1] <- p_reduced[,1] + centre[1]
    p_reduced[,2] <- p_reduced[,2] + centre[2]
    points <- as.ppp(p_reduced, c(centre[1]-radius, centre[1]+radius, centre[2]-radius, centre[2]+radius))
    return(points)
  }
  if (distribution=='gaussian'){
    p <- mvrnorm(n=n, mu=c(0,0), Sigma=matrix(c(sig^2,0,0,sig^2), 2, 2))
    p[,1] <- p[,1] + centre[1]
    p[,2] <- p[,2] + centre[2]
    points <- as.ppp(p, c(centre[1]-radius, centre[1]+radius, centre[2]-radius, centre[2]+radius))
    return(points)
  }
}