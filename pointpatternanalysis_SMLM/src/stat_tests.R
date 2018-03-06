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

### Compare set1 and set2 (monoVSdual, noisyVSnoisefree, ...)
### select data

# SET1
set1_name <- "obs"
pptype='unmarked'; units = 'pixels'; units_out = 'nm'; unit_size=160
levels1 = 'CTCF'; levels2 = ''
exp_name1 <- "DMSO"; exp_name2 <- "ActD"  # NA
path_to_experiment1 = '/home/alba/ISIS/nfs/users/jsolon/agranados/data/vicky/2017-06-15_HeLa_antiCTCF_DMSO_ActD/noisy/2017-06-15_HeLa_antiCTCF_DMSO'
path_to_experiment2 = '/home/alba/ISIS/nfs/users/jsolon/agranados/data/vicky/2017-06-15_HeLa_antiCTCF_DMSO_ActD/noisy/2017-06-15_HeLa_antiCTCF_ActD'

storm_file = 0
exp_names <- c(exp_name1, exp_name2)
path_to_experiments <- c(path_to_experiment1, path_to_experiment2)

## compute
RData_2experiments <- list.files(path = path_to_experiments[1], pattern = "\\_2experiments.RData$", 
                                 all.files = FALSE, full.names = FALSE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
load(paste(c(path_to_experiments[1], RData_2experiments), collapse='/'))

parameters1 <- parameters_g11
density1 <- density$level1

# SET2
set2_name <- "sim"
pptype='marked'; units = 'pixels'; units_out = 'nm'; unit_size=160 #nm/px
exp_name1 <- "DMSO"; exp_name2 <- "ActD"  # NA
levels1 = 'CTCF'; levels2 = ''
path_to_experiment1 <- '/home/alba/ISIS/nfs/users/jsolon/agranados/data/vicky/2017-06-15_HeLa_antiCTCF_DMSO_ActD/noisy/simulations/expsq'
path_to_experiment2 <- '/home/alba/ISIS/nfs/users/jsolon/agranados/data/vicky/2017-06-15_HeLa_antiCTCF_DMSO_ActD/noisy/simulations/expsq'

storm_file=0
exp_names <- c(exp_name1, exp_name2)
path_to_experiments <- c(path_to_experiment1, path_to_experiment2)

RData_2experiments <- list.files(path = path_to_experiments[1], pattern = "\\_2experiments.RData$", 
                                 all.files = FALSE, full.names = FALSE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
load(paste(c(path_to_experiments[1], RData_2experiments), collapse='/'))

parameters2 <- parameters_g11
density2 <- density$level1

# experiment 1 comparison:
amplitude_g1 <- exp(parameters1$beta0)  # set1
length_scale_g1 <- -parameters1$beta1^(-1)
clusterradiusR_g1 <- 2*sqrt(length_scale_g1)
phi_g1 <- amplitude_g1/4   # rho_cluster/rho_average
Nclusters_g1 <- amplitude_g1*pi*length_scale_g1*density1  # average number of points per cluster
kappa_g1 <- 1/(amplitude_g1*pi*length_scale_g1)  

amplitude_g2 <- exp(parameters2$beta0)  # set2
length_scale_g2 <- -parameters2$beta1^(-1)
clusterradiusR_g2 <- 2*sqrt(length_scale_g2)
phi_g2 <- amplitude_g2/4   # rho_cluster/rho_average
Nclusters_g2 <- amplitude_g2*pi*length_scale_g2*density2  #  average number of points per cluster
kappa_g2 <- 1/(amplitude_g2*pi*length_scale_g2)  

pexp_1 <- parameters1$experiment  # set1
pexp_2 <- parameters2$experiment  # set2

parameters_r0.test <- t.test(parameters1$r0[which(parameters1$experiment == exp_name1)],
                              parameters2$r0[which(parameters2$experiment == exp_name1)])
density.test <- t.test(density1[pexp_1==exp_name1],density2[pexp_2==exp_name1])
amplitude.test <- t.test(amplitude_g1[pexp_1==exp_name1],amplitude_g2[pexp_2==exp_name1])
length_scale.test <- t.test(length_scale_g1[pexp_1==exp_name1],length_scale_g2[pexp_2==exp_name1])
clusterradiusR.test <- t.test(clusterradiusR_g1[pexp_1==exp_name1],clusterradiusR_g2[pexp_2==exp_name1])
phi.test <- t.test(phi_g1[pexp_1==exp_name1],phi_g2[pexp_2==exp_name1])
Nclusters.test <- t.test(Nclusters_g1[pexp_1==exp_name1],Nclusters_g2[pexp_2==exp_name1])
kappa.test <- t.test(kappa_g1[pexp_1==exp_name1],kappa_g2[pexp_2==exp_name1])

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_density.pdf"), collapse = ""))
set1 <- density1[pexp_1==exp_name1]
set2 <- density2[pexp_2==exp_name1]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim,
              names=c(set1_name, set2_name), ylab=bquote('Average density ' ~ .(levels1) ~ ' [points/' ~ nm^2 ~']'))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), 
         y1 =  1.04*max(density1, density2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(density1, density2, na.rm=TRUE), get_asterisk(density.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_density.pdf"), collapse = ""))

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_r0.pdf"), collapse = ""))
set1 <- parameters1$r0[which(pexp_1 == exp_name1)]
set2 <- parameters2$r0[which(pexp_2 == exp_name1)]
ylim <- c(50,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, names=c(set1_name, set2_name), ylab=expression(paste('Correlation range ', r[0], ' [nm]')), 
              ylim=ylim)
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE), get_asterisk(parameters_r0.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_r0.pdf"), collapse = ""))
cat('\t\'statistical_tests_r0.pdf\' created.\n')

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_amplitude.pdf"), collapse = ""))
set1 <- amplitude_g1[which(pexp_1 == exp_name1)]
set2 <- amplitude_g2[which(pexp_2 == exp_name1)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, names=c(set1_name, set2_name), ylab=expression("Amplitude A"))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_amplitude.pdf"), collapse = ""))
cat('\t\'statistical_tests_amplitude.pdf\' created.\n')

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_length_scale.pdf"), collapse = ""))
set1 <- length_scale_g1[which(pexp_1 == exp_name1)]
set2 <- length_scale_g2[which(pexp_2 == exp_name1)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim, 
              names=c(set1_name, set2_name), ylab=expression("Length scale [nm]"))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_length_scale.pdf"), collapse = ""))
cat('\t\'statistical_tests_length_scale.pdf\' created.\n')

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_clusterradiusR.pdf"), collapse = ""))
set1 <- clusterradiusR_g1[which(pexp_1 == exp_name1)]
set2 <- clusterradiusR_g2[which(pexp_2 == exp_name1)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim, 
              names=c(set1_name, set2_name), ylab=expression("cluster radius R [nm]"))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_clusterradiusR.pdf"), collapse = ""))
cat('\t\'statistical_clusterradiusR_scale.pdf\' created.\n')

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_phi.pdf"), collapse = ""))
set1 <- phi_g1[which(pexp_1 == exp_name1)]
set2 <- phi_g2[which(pexp_2 == exp_name1)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim, 
              names=c(set1_name, set2_name), ylab=expression(phi^(cluster)))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_phi.pdf"), collapse = ""))
cat('\t\'statistical_phi_scale.pdf\' created.\n')

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_Nclusters.pdf"), collapse = ""))
set1 <- Nclusters_g1[which(pexp_1 == exp_name1)]
set2 <- Nclusters_g2[which(pexp_2 == exp_name1)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim, 
              names=c(set1_name, set2_name), ylab=expression(paste(N^cluster, ' [points/cluster]')))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_Nclusters.pdf"), collapse = ""))
cat('\t\'statistical_Nclusters_scale.pdf\' created.\n')

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_kappa.pdf"), collapse = ""))
set1 <- kappa_g1[which(pexp_1 == exp_name1)]
set2 <- kappa_g2[which(pexp_2 == exp_name1)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim, 
              names=c(set1_name, set2_name), ylab=expression(paste(kappa, " [clusters/", nm^2, "]")))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name1, "_kappa.pdf"), collapse = ""))
cat('\t\'statistical_kappa_scale.pdf\' created.\n')

parameters_r0.test <- t.test(parameters1$r0[which(parameters1$experiment == exp_name2)],
                             parameters2$r0[which(parameters2$experiment == exp_name2)])
density.test <- t.test(density1[pexp_1==exp_name2],density2[pexp_2==exp_name2])
amplitude.test <- t.test(amplitude_g1[pexp_1==exp_name2],amplitude_g2[pexp_2==exp_name2])
length_scale.test <- t.test(length_scale_g1[pexp_1==exp_name2],length_scale_g2[pexp_2==exp_name2])
clusterradiusR.test <- t.test(clusterradiusR_g1[pexp_1==exp_name2],clusterradiusR_g2[pexp_2==exp_name2])
phi.test <- t.test(phi_g1[pexp_1==exp_name2],phi_g2[pexp_2==exp_name2])
Nclusters.test <- t.test(Nclusters_g1[pexp_1==exp_name2],Nclusters_g2[pexp_2==exp_name2])
kappa.test <- t.test(kappa_g1[pexp_1==exp_name2],kappa_g2[pexp_2==exp_name2])

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_density.pdf"), collapse = ""))
set1 <- density1[pexp_1==exp_name2]
set2 <- density2[pexp_2==exp_name2]
ylim <- c(0, 1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim,
              names=c(set1_name, set2_name), ylab=bquote('Average density ' ~ .(levels1) ~ ' [points/' ~ nm^2 ~']'))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), 
         y1 =  1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE), get_asterisk(density.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_density.pdf"), collapse = ""))

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_r0.pdf"), collapse = ""))
set1 <- parameters1$r0[which(pexp_1 == exp_name2)]
set2 <- parameters2$r0[which(pexp_2 == exp_name2)]
ylim <- c(50,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, names=c(set1_name, set2_name), ylab=expression(paste('Correlation range ', r[0], ' [nm]')), 
              ylim=ylim)
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE), get_asterisk(parameters_r0.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_r0.pdf"), collapse = ""))
cat('\t\'statistical_tests_r0.pdf\' created.\n')

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_amplitude.pdf"), collapse = ""))
set1 <- amplitude_g1[which(pexp_1 == exp_name2)]
set2 <- amplitude_g2[which(pexp_2 == exp_name2)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim, 
              names=c(set1_name, set2_name), ylab=expression("Amplitude A"))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_amplitude.pdf"), collapse = ""))
cat('\t\'statistical_tests_amplitude.pdf\' created.\n')

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_length_scale.pdf"), collapse = ""))
set1 <- length_scale_g1[which(pexp_1 == exp_name2)]
set2 <- length_scale_g2[which(pexp_2 == exp_name2)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim,
              names=c(set1_name, set2_name), ylab=expression("Length scale [nm]"))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_length_scale.pdf"), collapse = ""))
cat('\t\'statistical_tests_length_scale.pdf\' created.\n')

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_clusterradiusR.pdf"), collapse = ""))
set1 <- clusterradiusR_g1[which(pexp_1 == exp_name2)]
set2 <- clusterradiusR_g2[which(pexp_2 == exp_name2)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim,
              names=c(set1_name, set2_name), ylab=expression("Cluster radius [nm]"))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_clusterradiusR.pdf"), collapse = ""))
cat('\t\'statistical_tests_clusterradiusR.pdf\' created.\n')

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_phi.pdf"), collapse = ""))
set1 <- phi_g1[which(pexp_1 == exp_name2)]
set2 <- phi_g2[which(pexp_2 == exp_name2)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim, 
              names=c(set1_name, set2_name), ylab=expression(phi^(cluster)))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_phi.pdf"), collapse = ""))
cat('\t\'statistical_phi_scale.pdf\' created.\n')

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_Nclusters.pdf"), collapse = ""))
set1 <- Nclusters_g1[which(pexp_1 == exp_name2)]
set2 <- Nclusters_g2[which(pexp_2 == exp_name2)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim, 
              names=c(set1_name, set2_name), ylab=expression(paste(N^cluster, ' [points/cluster]')))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_Nclusters.pdf"), collapse = ""))
cat('\t\'statistical_Nclusters_scale.pdf\' created.\n')

openpdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_kappa.pdf"), collapse = ""))
set1 <- kappa_g1[which(pexp_1 == exp_name2)]
set2 <- kappa_g2[which(pexp_2 == exp_name2)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim, 
              names=c(set1_name, set2_name), ylab=expression(paste(kappa, " [clusters/", nm^2, "]")))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c(paste(c(set1_name, set2_name, '_'), collapse=''), "statistical_tests_", levels1, '_', exp_name2, "_kappa.pdf"), collapse = ""))
cat('\t\'statistical_kappa_scale.pdf\' created.\n')


