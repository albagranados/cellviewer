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
                   c(min(data1[,2]),max(data1[,2])),marks=m1)
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
                                       units = 'pixels', unit_size = 1,
                                       plot_2experiments=TRUE, save=TRUE, stat_rrange_probs=0.5)
{
  cat('Generating longitudinal data...\n')
  exp_name1 <- exp_names[1]; exp_name2 <- exp_names[2]
  path_to_experiment1 <- path_to_experiments[1];   path_to_experiment2 <- path_to_experiments[2]
  
  path_to_RData1 <- paste(c(path_to_experiment1,paste(c(levels1,levels2,'_',exp_name1, '.RData'), collapse = '')), collapse='/')
  path_to_RData2 <- paste(c(path_to_experiment2,paste(c(levels1,levels2,'_',exp_name2, '.RData'), collapse = '')), collapse='/')
  
  if (pptype == "marked"){
    load(path_to_RData1)
    r_eval <- r_eval*unit_size
    data.g12.exp1 <- t(sapply(g12_all, `[[`, 'g12'))[,which(r_eval>=20)[1]:(length(r_eval))]
    data.g22.exp1 <- t(sapply(g22_all, `[[`, 'g22'))[,which(r_eval>=20)[1]:(length(r_eval))]
    data.K12.exp1 <- t(sapply(K12_all, `[[`, 'K12'))[,which(r_eval>=20)[1]:(length(r_eval))]
    data.K22.exp1 <- t(sapply(K22_all, `[[`, 'K22'))[,which(r_eval>=20)[1]:(length(r_eval))]
    g12_r0_exp1 <- unit_size*g12_r0_all
    g12_rend_exp1 <- which(r_eval>=quantile(g12_r0_exp1, probs=stat_rrange_probs, na.rm=TRUE))[1]
    g22_r0_exp1 <- unit_size*g22_r0_all
    g22_rend_exp1 <- which(r_eval>=quantile(g22_r0_exp1, probs=stat_rrange_probs, na.rm=TRUE))[1]
    intensities.exp1 <- intensities/(unit_size^2)
    
    cat('for statistical analysis, fill with NA function values above r = ', r_eval[g12_rend_exp1], '\n')
    load(path_to_RData2)
    r_eval <- r_eval*unit_size
    data.g12.exp2 <- t(sapply(g12_all, `[[`, 'g12'))[,which(r_eval>=20)[1]:(length(r_eval))]
    data.g22.exp2 <- t(sapply(g22_all, `[[`, 'g22'))[,which(r_eval>=20)[1]:(length(r_eval))]
    data.K12.exp2 <- t(sapply(K12_all, `[[`, 'K12'))[,which(r_eval>=20)[1]:(length(r_eval))]
    data.K22.exp2 <- t(sapply(K22_all, `[[`, 'K22'))[,which(r_eval>=20)[1]:(length(r_eval))]
    g12_r0_exp2 <- unit_size*g12_r0_all
    g12_rend_exp2 <- which(r_eval>=quantile(g12_r0_exp2, probs=stat_rrange_probs, na.rm=TRUE))[1]
    g22_r0_exp2 <- unit_size*g22_r0_all
    g22_rend_exp2 <- which(r_eval>=quantile(g22_r0_exp2, probs=stat_rrange_probs, na.rm=TRUE))[1]
    intensities.exp2 <- intensities/(unit_size^2)
  
    cat('for statistical analysis, fill with NA function values above r = ', r_eval[g12_rend_exp2], '\n')
    
    data_temp <- data.frame( g12 = (rbind(data.g12.exp1, data.g12.exp2)), g22 = (rbind(data.g22.exp1, data.g22.exp2)), 
                             K12 = (rbind(data.K12.exp1, data.K12.exp2)), K22 = (rbind(data.K22.exp1, data.K22.exp2)))
    data_wide <- data.frame(id_cell = factor(1:nrow(data_temp)),
                            experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(data.g12.exp1), nrow(data.g12.exp2)))),
                            data_temp)
    data_wide$experiment <- relevel(data_wide$experiment, ref=exp_name1)
    
    g12_r0 <- data.frame(id_cell = factor(1:nrow(data_temp)),
                          experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(data.g12.exp1),
                                                                             nrow(data.g12.exp2)))), 
                          r0 = c(g12_r0_exp1, g12_r0_exp2))
    
    g22_r0 <- data.frame(id_cell = factor(1:nrow(data_temp)),
                         experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(data.g22.exp1),
                                                                            nrow(data.g22.exp2)))), 
                         r0 = c(g22_r0_exp1, g22_r0_exp2))
    
    density <- data.frame(id_cell = factor(1:nrow(data_temp)),
                      experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(intensities.exp1),
                                                                         nrow(intensities.exp2)))), 
                      label1 = c(intensities.exp1[,1], intensities.exp2[,1]), 
                      label2 = c(intensities.exp1[,2], intensities.exp2[,2]))
    
    # convert wide-formatted data into long
    num_measurements <- dim(data.g12.exp1)[2]
    data_allrange <- reshape(data_wide, varying=list(names(data_wide)[3:(3+num_measurements-1)],
                                                     names(data_wide)[(3+num_measurements):(3+2*num_measurements-1)], 
                                                     names(data_wide)[(3+2*num_measurements):(3+3*num_measurements-1)],
                                                     names(data_wide)[(3+3*num_measurements):(3+4*num_measurements-1)]), 
                             idvar=c("id_cell", "experiment"),
                             direction="long", timevar='r', times=r_eval[(length(r_eval)-(dim(data.g12.exp1)[2])+1):length(r_eval)], 
                             v.names=c('g12', 'g22', 'K12', 'K22'))
    data_allrange$experiment <- relevel(data_allrange$experiment,ref=exp_name1)
    
    # data_wide[which(data_wide$experiment==exp_name1), (g12_rend_exp1+2):length(data_wide[1,])] <- NA
    # data_wide[which(data_wide$experiment==exp_name2), (g12_rend_exp2+2):length(data_wide[1,])] <- NA
    data_wide[which(data_wide$experiment==exp_name1), (2+g12_rend_exp1):(2+num_measurements)] <- NA
    data_wide[which(data_wide$experiment==exp_name1), (2+num_measurements+g22_rend_exp1):(2+2*num_measurements)] <- NA
    data_wide[which(data_wide$experiment==exp_name2), (2+g12_rend_exp2):(2+num_measurements)] <- NA
    data_wide[which(data_wide$experiment==exp_name2), (2+num_measurements+g22_rend_exp2):(2+2*num_measurements)] <- NA
    data <- reshape(data_wide, varying=list(names(data_wide)[3:(3+num_measurements-1)],
                                            names(data_wide)[(3+num_measurements):(3+2*num_measurements-1)], 
                                            names(data_wide)[(3+2*num_measurements):(3+3*num_measurements-1)],
                                            names(data_wide)[(3+3*num_measurements):(3+4*num_measurements-1)]), 
                             idvar=c("id_cell", "experiment"),
                             direction="long", timevar='r', times=r_eval[(length(r_eval)-(dim(data.g12.exp1)[2])+1):length(r_eval)], 
                             v.names=c('g12', 'g22', 'K12', 'K22'))
    data$experiment <- relevel(data_allrange$experiment,ref=exp_name1)
    
    r_eval_plot <- which(r_eval>=20)[1]:(0.8*length(r_eval))
    ylim=c(0.5,5)
    # r_eval_plot <- which(r_eval>=20)[1]:(0.5*length(r_eval)); ylim=c(0.5,2.8)
    if (plot_2experiments & (length(g12_all)>0)){
      openpdf("g12_2experiments.pdf")
      # col_palette <- colorRampPalette(c("black", "grey60"))(length(exp_names)); lty=c(1,2)
      col_palette <-  c(rgb(173,216,230,max = 255), rgb(255,165,0,max = 255,alpha=125)); lty=c(1,2)
      xlim <- c(r_eval[r_eval_plot][1],r_eval[r_eval_plot][length(r_eval[r_eval_plot])])
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''), 
           ylab=eval(bquote(expression(g[.(levels1)][','][.(levels2)](r)))), type='n', ylim = ylim)
      for (ii in 1:length(exp_names)){
        exp_name <- exp_names[ii]
        subset <- subset(data_allrange, data$experiment %in% exp_name)
        for (jj in unique(subset$id_cell)){
          with(subset[which(subset$id_cell == jj), ], lines(r[r_eval_plot], g12[r_eval_plot], lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
      }
      abline(h=1, lty=2, col='black')
      legend("topright", legend=exp_names, lty=lty, col=col_palette)
      closepdf("g12_2experiments.pdf")
      cat('\t\'g12_2experiments.pdf\' created.\n')
    }
    if (plot_2experiments & (length(g22_all)>0)){
      openpdf("g22_2experiments.pdf")
      # col_palette <- colorRampPalette(c("black", "grey60"))(length(exp_names)); lty=c(1,2)
      col_palette <-  c(rgb(173,216,230,max = 255), rgb(255,165,0,max = 255,alpha=125)); lty=c(1,2)
      xlim <- c(r_eval[r_eval_plot][1],r_eval[r_eval_plot][length(r_eval[r_eval_plot])])
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''), 
           ylab=eval(bquote(expression(g[.(levels1)][','][.(levels2)](r)))), type='n', ylim = ylim)
      for (ii in 1:length(exp_names)){
        exp_name <- exp_names[ii]
        subset <- subset(data_allrange, data$experiment %in% exp_name)
        for (jj in unique(subset$id_cell)){
          with(subset[which(subset$id_cell == jj), ], lines(r[r_eval_plot], g22[r_eval_plot], lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
      }
      abline(h=1, lty=2, col='black')
      legend("topright", legend=exp_names, lty=lty, col=col_palette)
      closepdf("g22_2experiments.pdf")
      cat('\t\'g22_2experiments.pdf\' created.\n')
    }
    ylim=c(0,60)
    if (plot_2experiments & (length(K12_all)>0)){
      openpdf("K12_2experiments.pdf")
      # col_palette <- colorRampPalette(c("black", "grey60"))(length(exp_names)); lty=c(1,2)
      col_palette <-  c(rgb(173,216,230,max = 255), rgb(255,165,0,max = 255,alpha=125)); lty=c(1,2)
      xlim <- c(r_eval[r_eval_plot][1],r_eval[r_eval_plot][length(r_eval[r_eval_plot])])
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''), 
           ylab=eval(bquote(expression(K[.(levels1)][','][.(levels2)](r)))), type='n', ylim = ylim)
      for (ii in 1:length(exp_names)){
        exp_name <- exp_names[ii]
        subset <- subset(data_allrange, data$experiment %in% exp_name)
        for (jj in unique(subset$id_cell)){
          with(subset[which(subset$id_cell == jj), ], lines(r[r_eval_plot], K12[r_eval_plot], lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
      }
      lines(r_eval[r_eval_plot], pi*r_eval[r_eval_plot]^2, col="black", lty=2)
      legend("topleft", legend=exp_names, lty=lty, col=col_palette)
      closepdf("K12_2experiments.pdf")
      cat('\t\'K12_2experiments.pdf\' created.\n')
    }
    if (save){
      save(data, data_allrange, g12_r0, g22_r0, density, r_eval, levels1, levels2, 
           file=paste(c(paste(c(path_to_experiments[1], paste(c(levels1,levels2), collapse = '')), collapse='/'), 
                                                        '_2experiments.RData'), collapse = ''))
    }
  }
  if (pptype == "unmarked"){
    load(path_to_RData1)
    r_eval <- r_eval*unit_size
    print(g11_all)
    print(which(r_eval>=20)[1]:(length(r_eval)))
    data.g11.exp1 <- t(sapply(g11_all, `[[`, 'g11'))[,which(r_eval>=20)[1]:(length(r_eval))]
    data.K11.exp1 <- t(sapply(K11_all, `[[`, 'K11'))[,which(r_eval>=20)[1]:(length(r_eval))]
    g11_r0_exp1 <- unit_size*g11_r0_all
    g11_rend_exp1 <- which(r_eval>=quantile(g11_r0_exp1, probs=stat_rrange_probs, na.rm=TRUE))[1]
    intensities.exp1 <- intensities/(unit_size^2)
    cat('for statistical analysis, fill with NA function values above r = ', r_eval[g11_rend_exp1], '\n')

    load(path_to_RData2)
    r_eval <- r_eval*unit_size
    data.g11.exp2 <- t(sapply(g11_all, `[[`, 'g11'))[,which(r_eval>=20)[1]:(length(r_eval))]
    data.K11.exp2 <- t(sapply(K11_all, `[[`, 'K11'))[,which(r_eval>=20)[1]:(length(r_eval))]
    g11_r0_exp2 <- unit_size*g11_r0_all
    g11_rend_exp2 <- which(r_eval>=quantile(g11_r0_exp2, probs=stat_rrange_probs, na.rm=TRUE))[1]
    intensities.exp2 <- intensities/(unit_size^2)
    cat('for statistical analysis, fill with NA function values above r = ', r_eval[g11_rend_exp2], '\n')
    
    data_temp <- data.frame( g11 = (rbind(data.g11.exp1, data.g11.exp2)), 
                             K11 = (rbind(data.K11.exp1, data.K11.exp2)))
    data_wide <- data.frame(id_cell = factor(1:nrow(data_temp)),
                            experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(data.g11.exp1), nrow(data.g11.exp2)))),
                            data_temp)
    data_wide$experiment <- relevel(data_wide$experiment, ref=exp_name1)
    
    g11_r0 <- data.frame(id_cell = factor(1:nrow(data_temp)),
                         experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(data.g11.exp1),
                                                                            nrow(data.g11.exp2)))), 
                         r0 = c(g11_r0_exp1, g11_r0_exp2))
    density <- data.frame(id_cell = factor(1:nrow(data_temp)),
                          experiment = factor(rep(c(exp_name1, exp_name2), c(nrow(intensities.exp1),
                                                                             nrow(intensities.exp2)))), 
                          label1 = c(intensities.exp1[,1], intensities.exp2[,1]))
    
    # convert wide-formatted data into long
    num_measurements <- dim(data.g11.exp1)[2]
    data_allrange <- reshape(data_wide, varying=list(names(data_wide)[3:(3+num_measurements-1)],
                                                     names(data_wide)[(3+num_measurements):(3+2*num_measurements-1)]), 
                             idvar=c("id_cell", "experiment"),
                             direction="long", timevar='r', times=r_eval[(length(r_eval)-(dim(data.g11.exp1)[2])+1):length(r_eval)], 
                             v.names=c('g11', 'K11'))
    data_allrange$experiment <- relevel(data_allrange$experiment,ref=exp_name1)
    
    # data_wide[which(data_wide$experiment==exp_name1), (g11_rend_exp1+2):length(data_wide[1,])] <- NA
    # data_wide[which(data_wide$experiment==exp_name2), (g11_rend_exp2+2):length(data_wide[1,])] <- NA
    data_wide[which(data_wide$experiment==exp_name1), (2+g11_rend_exp1):(2+num_measurements)] <- NA
    data_wide[which(data_wide$experiment==exp_name2), (2+g11_rend_exp2):(2+num_measurements)] <- NA
    data <- reshape(data_wide, varying=list(names(data_wide)[3:(3+num_measurements-1)],
                                            names(data_wide)[(3+num_measurements):(3+2*num_measurements-1)]), 
                    idvar=c("id_cell", "experiment"),
                    direction="long", timevar='r', times=r_eval[(length(r_eval)-(dim(data.g11.exp1)[2])+1):length(r_eval)], 
                    v.names=c('g11', 'K11'))
    data$experiment <- relevel(data_allrange$experiment,ref=exp_name1)
    
    r_eval_plot <- which(r_eval>=20)[1]:(0.8*length(r_eval))
    ylim=c(0.5,5)
    # r_eval_plot <- which(r_eval>=20)[1]:(0.5*length(r_eval)); ylim=c(0.5,2.8)
    if (plot_2experiments & (length(g11_all)>0)){
      openpdf("g11_2experiments.pdf")
      # col_palette <- colorRampPalette(c("black", "grey60"))(length(exp_names)); lty=c(1,2)
      col_palette <-  c(rgb(173,216,230,max = 255), rgb(255,165,0,max = 255,alpha=125)); lty=c(1,2)
      xlim <- c(r_eval[r_eval_plot][1],r_eval[r_eval_plot][length(r_eval[r_eval_plot])])
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''), 
           ylab=eval(bquote(expression(g[.(levels1)][','][.(levels2)](r)))), type='n', ylim = ylim)
      for (ii in 1:length(exp_names)){
        exp_name <- exp_names[ii]
        subset <- subset(data_allrange, data$experiment %in% exp_name)
        for (jj in unique(subset$id_cell)){
          with(subset[which(subset$id_cell == jj), ], lines(r[r_eval_plot], g11[r_eval_plot], lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
      }
      abline(h=1, lty=2, col='black')
      legend("topright", legend=exp_names, lty=lty, col=col_palette)
      closepdf("g11_2experiments.pdf")
      cat('\t\'g11_2experiments.pdf\' created.\n')
    }
    ylim=c(0,60)
    if (plot_2experiments & (length(K11_all)>0)){
      openpdf("K11_2experiments.pdf")
      # col_palette <- colorRampPalette(c("black", "grey60"))(length(exp_names)); lty=c(1,2)
      col_palette <-  c(rgb(173,216,230,max = 255), rgb(255,165,0,max = 255,alpha=125)); lty=c(1,2)
      xlim <- c(r_eval[r_eval_plot][1],r_eval[r_eval_plot][length(r_eval[r_eval_plot])])
      plot(1, xlim=xlim, xlab=paste('r ','[', units,']', sep = ''), 
           ylab=eval(bquote(expression(K[.(levels1)][','][.(levels2)](r)))), type='n', ylim = ylim)
      for (ii in 1:length(exp_names)){
        exp_name <- exp_names[ii]
        subset <- subset(data_allrange, data$experiment %in% exp_name)
        for (jj in unique(subset$id_cell)){
          with(subset[which(subset$id_cell == jj), ], lines(r[r_eval_plot], K11[r_eval_plot], lty=lty[ii], main='', lwd=1, col=col_palette[ii]))
        }
      }
      lines(r_eval[r_eval_plot], pi*r_eval[r_eval_plot]^2, col="black", lty=2)
      legend("topleft", legend=exp_names, lty=lty, col=col_palette)
      closepdf("K11_2experiments.pdf")
      cat('\t\'K11_2experiments.pdf\' created.\n')
    }
    
    if (save){
      save(data, data_allrange, g11_r0, density, r_eval, levels1, levels2, 
           file=paste(c(paste(c(path_to_experiments[1], paste(c(levels1,levels2), collapse = '')), collapse='/'), 
                        '_2experiments.RData'), collapse = ''))
    }
  }
}

statistical_analysis <- function(path_to_experiments, pptype, fitting='exp'){

  library(lme4)
  library(nlme)
  library(MASS)
  library(car)

  cat('Running statistical analysis...\n')  
  RData_2experiments <- list.files(path = path_to_experiments[1], pattern = "\\_2experiments.RData$", all.files = FALSE,
                                   full.names = FALSE, recursive = FALSE,
                                   ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  load(paste(c(path_to_experiments[1], RData_2experiments), collapse='/'))
  attach(data)
  
  exp1_name <- levels(data$experiment)[1]; exp2_name <- levels(data$experiment)[2];
  
  data$r <- data$r
  if(fitting %in% c('lin')){
    attach(data)
    cat('linear transformation... Done\n')
  }
  if (pptype == "marked"){
    
    ## exponential fitting
    if(fitting %in% c('exp')){
      data$g12 <- log(data$g12-1)
      data$g22 <- log(data$g22-1)
      attach(data)
      cat('Exponential transformation... Done\n')
    }
    
    ## response feature analysis (Faraway, 2006, p.206)
    # fit ordinary regressions separately to each group
    beta0_g12 <- numeric(nlevels(data$id_cell)); beta1_g12 <- numeric(nlevels(data$id_cell))
    K_r0_g12 <- numeric(nlevels(data$id_cell))
    beta0_g22 <- numeric(nlevels(data$id_cell)); beta1_g22 <- numeric(nlevels(data$id_cell))
    K_r0_g22 <- numeric(nlevels(data$id_cell))
    for(i in 1:nlevels(data$id_cell)){
      lmod <- lm(g12 ~ r, subset=(id_cell==i), data=data)
      beta0_g12[i] <- coef(lmod)[1]
      beta1_g12[i] <- coef(lmod)[2]
      a <- tail(which(r_eval < 0.5*g12_r0$r0[which(g12_r0$id_cell == i)]),n=1)
      if (length(a)==0){ K_r0_g12[i] <- NA}
      else{
        K_r0_g12[i] <- K12[id_cell == i][a]
      }
      lmod <- lm(g22 ~ r, subset=(id_cell==i), data=data)
      beta0_g22[i] <- coef(lmod)[1]
      beta1_g22[i] <- coef(lmod)[2]
      a <- tail(which(r_eval < 0.5*g12_r0$r0[which(g22_r0$id_cell == i)]),n=1)
      if (length(a)==0){ K_r0_g22[i] <- NA}
      else{
        K_r0_g22[i] <- K22[id_cell == i][a]
      }
    }
    if(fitting %in% c('exp')){
      amplitude_g12 <- exp(beta0_g12)
      length_scale_g12 <- -beta1_g12^(-1)
      area_g12 <- -length_scale_g12*amplitude_g12*(exp(-1/(length_scale_g12)*g12_r0$r0)-1)
    
      amplitude_g22 <- exp(beta0_g22)
      length_scale_g22 <- -beta1_g22^(-1)
      area_g22 <- -length_scale_g22*amplitude_g22*(exp(-1/(length_scale_g12)*g22_r0$r0)-1)
      
      # tests (experiments)
      cat('Running tests... ')
      pexp <- data$experiment[match(1:nlevels(data$id_cell),data$id_cell)]
      g12_r0.test <- t.test(g12_r0$r0[which(g12_r0$experiment == exp1_name)], g12_r0$r0[which(g12_r0$experiment == exp2_name)])
      density_label1.test <- t.test(density$label1[pexp==exp1_name],density$label1[pexp==exp2_name])
      density_label2.test <- t.test(density$label2[pexp==exp1_name],density$label2[pexp==exp2_name])
      amplitude_g12.test <- t.test(amplitude_g12[pexp==exp1_name],amplitude_g12[pexp==exp2_name])
      length_scale_g12.test <- t.test(length_scale_g12[pexp==exp1_name],length_scale_g12[pexp==exp2_name])
      area_g12.test <- t.test(area_g12[pexp==exp1_name],area_g12[pexp==exp2_name])
      
      g22_r0.test <- t.test(g22_r0$r0[which(g22_r0$experiment == exp1_name)], g22_r0$r0[which(g22_r0$experiment == exp2_name)])
      density_label1.test <- t.test(density$label1[pexp==exp1_name],density$label1[pexp==exp2_name])
      density_label2.test <- t.test(density$label2[pexp==exp1_name],density$label2[pexp==exp2_name])
      amplitude_g22.test <- t.test(amplitude_g22[pexp==exp1_name],amplitude_g22[pexp==exp2_name])
      length_scale_g22.test <- t.test(length_scale_g22[pexp==exp1_name],length_scale_g22[pexp==exp2_name])
      area_g22.test <- t.test(area_g22[pexp==exp1_name],area_g22[pexp==exp2_name])
      cat('Done.\n')
      
      # # plot results
      cat('Plotting results of statistical tests...\n')
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests.pdf", collapse = "_"), pointsize=11)
      # old <- par(mfrow=c(2, 2))
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests", levels1, "density.pdf"), collapse = "_"))
      ylim <- c(0,1.2*max(density$label1, na.rm=TRUE))
      bp <- boxplot(split(density$label1, pexp), main=levels1, 
                    names=c(exp1_name, exp2_name), ylim=ylim, ylab=expression(paste('Average density [points/',nm^2,']')))
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(density$label1, na.rm=TRUE), y1 = 1.04*max(density$label1, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(density$label1, na.rm=TRUE),  get_asterisk(density_label1.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests", levels1, "density.pdf"), collapse = "_"))
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests", levels2, "density.pdf"), collapse = "_"))
      ylim <- c(0,1.2*max(density$label2, na.rm=TRUE))
      bp <- boxplot(split(density$label2, pexp), main=levels2, 
                    names=c(exp1_name, exp2_name), ylim=ylim, ylab=expression(paste('Average density [points/',nm^2,']')))
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(density$label2, na.rm=TRUE), y1 = 1.04*max(density$label2, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(density$label2, na.rm=TRUE),  get_asterisk(density_label2.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests", levels2, "density.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_density.pdf\' created.\n')
      
      openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_r0_g12.pdf", collapse = "_"))
      bp <- boxplot(split(g12_r0$r0,pexp), ylab=expression(paste('Correlation range ', r[0], ' [nm]')), ylim=c(50,1.2*max(g12_r0$r0, na.rm=TRUE)))
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(g12_r0$r0, na.rm=TRUE), y1 = 1.04*max(g12_r0$r0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(g12_r0$r0, na.rm=TRUE), get_asterisk(g12_r0.test) , cex=1)
      closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_r0_g12.pdf", collapse = "_"))
      cat('\t\'statistical_tests_r0_g12.pdf\' created.\n')
      
      openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_amplitude_g12.pdf", collapse = "_"))
      ylim <- c(0,1.2*max(amplitude_g12, na.rm=TRUE))
      bp <- boxplot(split(amplitude_g12,pexp), ylab=expression("Amplitude A"), names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(amplitude_g12, na.rm=TRUE), y1 = 1.04*max(amplitude_g12, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(amplitude_g12, na.rm=TRUE),  get_asterisk(amplitude_g12.test) , cex=1)
      closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_amplitude_g12.pdf", collapse = "_"))
      cat('\t\'statistical_tests_amplitude_g12.pdf\' created.\n')
      
      openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_lengthscale_g12.pdf", collapse = "_"))
      ylim <- c(10,1.2*max(length_scale_g12, na.rm=TRUE))
      bp <- boxplot(split(length_scale_g12,pexp), ylab=expression(paste("Length scale ", lambda, " [nm]")), 
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(length_scale_g12, na.rm=TRUE), y1 = 1.04*max(length_scale_g12, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(length_scale_g12, na.rm=TRUE) ,  get_asterisk(length_scale_g12.test) , cex=1)
      closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_lengthscale_g12.pdf", collapse = "_"))
      cat('\t\'statistical_tests_lengthscale_g12.pdf\' created.\n')
      
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_area_g12.pdf", collapse = "_"))
      # ylim <- c(0,1.2*max(area_g12, na.rm=TRUE))
      # bp <- boxplot(split(area_g12,pexp), main=expression("area"), names=c(exp1_name, exp2_name),  ylim=ylim)
      # segments(x0 = 1, x1 = 2, y0 = 1.04*max(area_g12, na.rm=TRUE), y1 = 1.04*max(area_g12, na.rm=TRUE), col = "black")
      # text( 1.5 , 1.09*max(area_g12, na.rm=TRUE) , get_asterisk(area_g12.test), cex=1.5)
      # closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_area.pdf", collapse = "_"))
      # cat('\t\'statistical_tests_area_g12.pdf\' created.\n')
      # # par(old); closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_g12.pdf", collapse = "_"))
      
      openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_r0_g22.pdf", collapse = "_"))
      bp <- boxplot(split(g22_r0$r0,pexp), ylab=expression(paste('Correlation range ', r[0], ' [nm]')), ylim=c(50,1.2*max(g12_r0$r0, na.rm=TRUE)))
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(g22_r0$r0, na.rm=TRUE), y1 = 1.04*max(g22_r0$r0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(g22_r0$r0, na.rm=TRUE), get_asterisk(g22_r0.test) , cex=1)
      closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_r0_g22.pdf", collapse = "_"))
      cat('\t\'statistical_tests_r0_g22.pdf\' created.\n')
      
      openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_amplitude_g22.pdf", collapse = "_"))
      ylim <- c(0,1.2*max(amplitude_g22, na.rm=TRUE))
      bp <- boxplot(split(amplitude_g22,pexp), ylab=expression("Amplitude A"), names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(amplitude_g22, na.rm=TRUE), y1 = 1.04*max(amplitude_g22, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(amplitude_g22, na.rm=TRUE),  get_asterisk(amplitude_g22.test) , cex=1)
      closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_amplitude_g22.pdf", collapse = "_"))
      cat('\t\'statistical_tests_amplitude_g22.pdf\' created.\n')
      
      openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_lengthscale_g22.pdf", collapse = "_"))
      ylim <- c(10,1.2*max(length_scale_g22, na.rm=TRUE))
      bp <- boxplot(split(length_scale_g22,pexp), ylab=expression(paste("Length scale ", lambda, " [nm]")), 
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(length_scale_g22, na.rm=TRUE), y1 = 1.04*max(length_scale_g22, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(length_scale_g22, na.rm=TRUE) ,  get_asterisk(length_scale_g22.test) , cex=1)
      closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_lengthscale_g22.pdf", collapse = "_"))
      cat('\t\'statistical_tests_lengthscale_g22.pdf\' created.\n')
    }
  
    if(fitting %in% c('lin')){
      area <- 0.5*beta0*g12_r0$r0
      
      # tests (experiments)
      cat('Running tests... ')
      pexp <- data$experiment[match(1:nlevels(data$id_cell),data$id_cell)]
      g12_r0.test <- t.test(g12_r0$r0[which(g12_r0$experiment == exp1_name)], g12_r0$r0[which(g12_r0$experiment == exp2_name)])
      beta0.test <- t.test(beta0[pexp==exp1_name],beta0[pexp==exp2_name])
      beta1.test <- t.test(beta1[pexp==exp1_name],beta1[pexp==exp2_name])
      area.test <- t.test(area[pexp==exp1_name],area[pexp==exp2_name])
      K_r0.test <- t.test(K_r0[which(g12_r0$experiment == exp1_name)], K_r0[which(g12_r0$experiment == exp2_name)])
      cat('Done.\n')
      
      # # plot results
      cat('Plotting results of statistical tests...\n')
      openpdf("statistical_tests.pdf", pointsize=11)
      old <- par(mfrow=c(2, 2))
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_r0.pdf", collapse = "_"))
      ylim <- NULL
      bp <- boxplot(split(g12_r0$r0,pexp), main=expression(r[0]), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.01*max(g12_r0$r0, na.rm=TRUE), y1 = 1.01*max(g12_r0$r0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.02*max(g12_r0$r0, na.rm=TRUE), get_asterisk(g12_r0.test) , cex=1.2)
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
  }
  if (pptype=="unmarked"){
    ## exponential fitting
    if(fitting %in% c('exp')){
      data$g11 <- log(data$g11-1)
      attach(data)
      cat('Exponential transformation... Done\n')
    }
    
    ## response feature analysis (Faraway, 2006, p.206)
    # fit ordinary regressions separately to each group
    beta0_g11 <- numeric(nlevels(data$id_cell)); beta1_g11 <- numeric(nlevels(data$id_cell))
    K_r0_g11 <- numeric(nlevels(data$id_cell))
    for(i in 1:nlevels(data$id_cell)){
      lmod <- lm(g11 ~ r, subset=(id_cell==i), data=data)
      beta0_g11[i] <- coef(lmod)[1]
      beta1_g11[i] <- coef(lmod)[2]
      a <- tail(which(r_eval < 0.5*g11_r0$r0[which(g11_r0$id_cell == i)]),n=1)
      if (length(a)==0){ K_r0_g11[i] <- NA}
      else{
        K_r0_g11[i] <- K11[id_cell == i][a]
      }
    }
    if(fitting %in% c('exp')){
      amplitude_g11 <- exp(beta0_g11)
      length_scale_g11 <- -beta1_g11^(-1)
      area_g11 <- -length_scale_g11*amplitude_g11*(exp(-1/(length_scale_g11)*g11_r0$r0)-1)
      
      # tests (experiments)
      cat('Running tests... ')
      pexp <- data$experiment[match(1:nlevels(data$id_cell),data$id_cell)]
      g11_r0.test <- t.test(g11_r0$r0[which(g11_r0$experiment == exp1_name)], g11_r0$r0[which(g11_r0$experiment == exp2_name)])
      density_label1.test <- t.test(density$label1[pexp==exp1_name],density$label1[pexp==exp2_name])
      density_label2.test <- t.test(density$label2[pexp==exp1_name],density$label2[pexp==exp2_name])
      amplitude_g11.test <- t.test(amplitude_g11[pexp==exp1_name],amplitude_g11[pexp==exp2_name])
      length_scale_g11.test <- t.test(length_scale_g11[pexp==exp1_name],length_scale_g11[pexp==exp2_name])
      area_g11.test <- t.test(area_g11[pexp==exp1_name],area_g11[pexp==exp2_name])
      cat('Done.\n')
      
      # # plot results
      cat('Plotting results of statistical tests...\n')
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests.pdf", collapse = "_"), pointsize=11)
      # old <- par(mfrow=c(2, 2))
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests", levels1, "density.pdf"), collapse = "_"))
      ylim <- c(0,1.2*max(density$label1, na.rm=TRUE))
      bp <- boxplot(split(density$label1, pexp), main=levels1, 
                    names=c(exp1_name, exp2_name), ylim=ylim, ylab=expression(paste('Average density [points/',nm^2,']')))
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(density$label1, na.rm=TRUE), y1 = 1.04*max(density$label1, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(density$label1, na.rm=TRUE),  get_asterisk(density_label1.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests", levels1, "density.pdf"), collapse = "_"))
      
      openpdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests", levels2, "density.pdf"), collapse = "_"))
      ylim <- c(0,1.2*max(density$label2, na.rm=TRUE))
      bp <- boxplot(split(density$label2, pexp), main=levels2, 
                    names=c(exp1_name, exp2_name), ylim=ylim, ylab=expression(paste('Average density [points/',nm^2,']')))
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(density$label2, na.rm=TRUE), y1 = 1.04*max(density$label2, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(density$label2, na.rm=TRUE),  get_asterisk(density_label2.test) , cex=1)
      closepdf(paste(c(unlist(strsplit(RData_2experiments, "\\_"))[1], "statistical_tests", levels2, "density.pdf"), collapse = "_"))
      cat('\t\'statistical_tests_density.pdf\' created.\n')
      
      openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_r0_g11.pdf", collapse = "_"))
      bp <- boxplot(split(g11_r0$r0,pexp), ylab=expression(paste('Correlation range ', r[0], ' [nm]')), ylim=c(50,1.2*max(g11_r0$r0, na.rm=TRUE)))
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(g11_r0$r0, na.rm=TRUE), y1 = 1.04*max(g11_r0$r0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(g11_r0$r0, na.rm=TRUE), get_asterisk(g11_r0.test) , cex=1)
      closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_r0_g11.pdf", collapse = "_"))
      cat('\t\'statistical_tests_r0_g11.pdf\' created.\n')
      
      openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_amplitude_g11.pdf", collapse = "_"))
      ylim <- c(0,1.2*max(amplitude_g11, na.rm=TRUE))
      bp <- boxplot(split(amplitude_g11,pexp), ylab=expression("Amplitude A"), names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(amplitude_g11, na.rm=TRUE), y1 = 1.04*max(amplitude_g11, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(amplitude_g11, na.rm=TRUE),  get_asterisk(amplitude_g11.test) , cex=1)
      closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_amplitude_g11.pdf", collapse = "_"))
      cat('\t\'statistical_tests_amplitude_g11.pdf\' created.\n')
      
      openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_lengthscale_g11.pdf", collapse = "_"))
      ylim <- c(10,1.2*max(length_scale_g11, na.rm=TRUE))
      bp <- boxplot(split(length_scale_g11,pexp), ylab=expression(paste("Length scale ", lambda, " [nm]")), 
                    names=c(exp1_name, exp2_name), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.04*max(length_scale_g11, na.rm=TRUE), y1 = 1.04*max(length_scale_g11, na.rm=TRUE), col = "black")
      text( 1.5 , 1.09*max(length_scale_g11, na.rm=TRUE) ,  get_asterisk(length_scale_g11.test) , cex=1)
      closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_lengthscale_g11.pdf", collapse = "_"))
      cat('\t\'statistical_tests_lengthscale_g11.pdf\' created.\n')
      
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_area_g11.pdf", collapse = "_"))
      # ylim <- c(0,1.2*max(area_g11, na.rm=TRUE))
      # bp <- boxplot(split(area_g11,pexp), main=expression("area"), names=c(exp1_name, exp2_name),  ylim=ylim)
      # segments(x0 = 1, x1 = 2, y0 = 1.04*max(area_g11, na.rm=TRUE), y1 = 1.04*max(area_g11, na.rm=TRUE), col = "black")
      # text( 1.5 , 1.09*max(area_g11, na.rm=TRUE) , get_asterisk(area_g11.test), cex=1.5)
      # closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_area.pdf", collapse = "_"))
      # cat('\t\'statistical_tests_area_g11.pdf\' created.\n')
      # # par(old); closepdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_g11.pdf", collapse = "_"))
    }
    
    if(fitting %in% c('lin')){
      area <- 0.5*beta0*g11_r0$r0
      
      # tests (experiments)
      cat('Running tests... ')
      pexp <- data$experiment[match(1:nlevels(data$id_cell),data$id_cell)]
      g11_r0.test <- t.test(g11_r0$r0[which(g11_r0$experiment == exp1_name)], g11_r0$r0[which(g11_r0$experiment == exp2_name)])
      beta0.test <- t.test(beta0[pexp==exp1_name],beta0[pexp==exp2_name])
      beta1.test <- t.test(beta1[pexp==exp1_name],beta1[pexp==exp2_name])
      area.test <- t.test(area[pexp==exp1_name],area[pexp==exp2_name])
      K_r0.test <- t.test(K_r0[which(g11_r0$experiment == exp1_name)], K_r0[which(g11_r0$experiment == exp2_name)])
      cat('Done.\n')
      
      # # plot results
      cat('Plotting results of statistical tests...\n')
      openpdf("statistical_tests.pdf", pointsize=11)
      old <- par(mfrow=c(2, 2))
      # openpdf(paste(unlist(strsplit(RData_2experiments, "\\_"))[1], "_statistical_tests_r0.pdf", collapse = "_"))
      ylim <- NULL
      bp <- boxplot(split(g11_r0$r0,pexp), main=expression(r[0]), ylim=ylim)
      segments(x0 = 1, x1 = 2, y0 = 1.01*max(g11_r0$r0, na.rm=TRUE), y1 = 1.01*max(g11_r0$r0, na.rm=TRUE), col = "black")
      text( 1.5 , 1.02*max(g11_r0$r0, na.rm=TRUE), get_asterisk(g11_r0.test) , cex=1.2)
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

