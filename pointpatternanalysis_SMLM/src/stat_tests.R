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

# mono
pptype='unmarked'; units = 'pixels'; units_out = 'nm'; unit_size=160
levels1 = 'SMC3'; levels2 = '_'
exp_name1 <- "DMSO"; exp_name2 <- "ActD"  # NA
path_to_experiment1 = "/home/alba/ISIS/nfs/users/jsolon/agranados/data/vicky/2017-06-15_HeLa_antiCTCF_DMSO_ActD/2017-06-15_HeLa_antiCTCF_DMSO"
path_to_experiment2 = "/home/alba/ISIS/nfs/users/jsolon/agranados/data/vicky/2017-06-15_HeLa_antiCTCF_DMSO_ActD/2017-06-15_HeLa_antiCTCF_ActD"

storm_file = 0
exp_names <- c(exp_name1, exp_name2)
path_to_experiments <- c(path_to_experiment1, path_to_experiment2)

## compute
RData_2experiments <- list.files(path = path_to_experiments[1], pattern = "\\_2experiments.RData$", 
                                 all.files = FALSE, full.names = FALSE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
load(paste(c(path_to_experiments[1], RData_2experiments), collapse='/'))

parameters1 <- parameters_g11
density1 <- density

# dual
pptype='marked'; units = 'pixels'; units_out = 'nm'; unit_size=160 #nm/px
exp_name1 <- "DMSO"; exp_name2 <- "ActD"  # NA
levels1 = 'SMC3'; levels2 = 'SMC3'
path_to_experiment1 <- '/home/alba/ISIS/nfs/users/jsolon/agranados/data/vicky/2017-07-17_HeLa_DualColor_SMC1_CTCF/SMC1_CTCF in DMSO Controls'
path_to_experiment2 <- '/home/alba/ISIS/nfs/users/jsolon/agranados/data/vicky/2017-07-17_HeLa_DualColor_SMC1_CTCF/SMC1_CTCF in ActD Treated'

storm_file=0
exp_names <- c(exp_name1, exp_name2)
path_to_experiments <- c(path_to_experiment1, path_to_experiment2)

RData_2experiments <- list.files(path = path_to_experiments[1], pattern = "\\_2experiments.RData$", 
                                 all.files = FALSE, full.names = FALSE, recursive = FALSE,
                                 ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
load(paste(c(path_to_experiments[1], RData_2experiments), collapse='/'))

parameters2 <- parameters_g22
density2 <- density

# experiment 1 comparison:
amplitude_g1 <- exp(parameters1$beta0)  # mono
amplitudeAtr0_g1 <- exp(parameters1$beta0*parameters1$beta1*20)+1
length_scale_g1 <- -parameters1$beta1^(-1)

amplitude_g2 <- exp(parameters2$beta0)  # dual
amplitudeAtr0_g2 <- exp(parameters2$beta0*parameters2$beta1*20)+1
length_scale_g2 <- -parameters2$beta1^(-1)

pexp_1 <- parameters1$experiment  # mono
pexp_2 <- parameters2$experiment  # dual

parameters_r0.test <- t.test(parameters1$r0[which(parameters1$experiment == exp_name1)],
                              parameters2$r0[which(parameters2$experiment == exp_name1)])
density.test <- t.test(density1$level1[pexp_1==exp_name1],density2$level2[pexp_2==exp_name1])
amplitude.test <- t.test(amplitude_g1[pexp_1==exp_name1],amplitude_g2[pexp_2==exp_name1])
amplitudeAtr0.test <- t.test(amplitudeAtr0_g1[pexp_1==exp_name1],amplitudeAtr0_g2[pexp_2==exp_name1])
length_scale.test <- t.test(length_scale_g1[pexp_1==exp_name1],length_scale_g2[pexp_2==exp_name1])

openpdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name1, "density.pdf"), collapse = "_"))
set1 <- density1$level1[pexp_1==exp_name1]
set2 <- density2$level2[pexp_2==exp_name1]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim,
              names=c('mono', 'dual'), ylab=bquote('Average density ' ~ .(levels2) ~ ' [points/' ~ nm^2 ~']'))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), 
         y1 =  1.04*max(density1$level1, density2$level2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(density1$level1, density2$level2, na.rm=TRUE), get_asterisk(density.test) , cex=1)
closepdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name1, "density.pdf"), collapse = "_"))

openpdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name1, "r0.pdf"), collapse = "_"))
set1 <- parameters1$r0[which(pexp_1 == exp_name1)]
set2 <- parameters2$r0[which(pexp_2 == exp_name1)]
ylim <- c(50,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, names=c('mono', 'dual'), ylab=expression(paste('Correlation range ', r[0], ' [nm]')), 
              ylim=ylim)
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE), get_asterisk(parameters_r0.test) , cex=1)
closepdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name1, "r0.pdf"), collapse = "_"))
cat('\t\'statistical_tests_r0.pdf\' created.\n')

openpdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name1, "amplitude.pdf"), collapse = "_"))
set1 <- amplitude_g1[which(pexp_1 == exp_name1)]
set2 <- amplitude_g2[which(pexp_2 == exp_name1)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, names=c('mono', 'dual'), ylab=expression("Amplitude A"))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name1, "amplitude.pdf"), collapse = "_"))
cat('\t\'statistical_tests_amplitude.pdf\' created.\n')

openpdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name1, "amplitudeAtr0.pdf"), collapse = "_"))
set1 <- amplitudeAtr0_g1[which(pexp_1 == exp_name1)]
set2 <- amplitudeAtr0_g2[which(pexp_2 == exp_name1)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, names=c('mono', 'dual'), ylab=expression("Correlation at 20 nm [a.u.]"))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitudeAtr0.test) , cex=1)
closepdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name1, "amplitudeAtr0.pdf"), collapse = "_"))
cat('\t\'statistical_tests_amplitudeAtr0.pdf\' created.\n')

openpdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name1, "length_scale.pdf"), collapse = "_"))
set1 <- length_scale_g1[which(pexp_1 == exp_name1)]
set2 <- length_scale_g2[which(pexp_2 == exp_name1)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim, 
              names=c('mono', 'dual'), ylab=expression("Length scale [nm]"))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name1, "length_scale.pdf"), collapse = "_"))
cat('\t\'statistical_tests_length_scale.pdf\' created.\n')

parameters_r0.test <- t.test(parameters1$r0[which(parameters1$experiment == exp_name2)],
                             parameters2$r0[which(parameters2$experiment == exp_name2)])
density.test <- t.test(density1$level1[pexp_1==exp_name2],density2$level2[pexp_2==exp_name2])
amplitude.test <- t.test(amplitude_g1[pexp_1==exp_name2],amplitude_g2[pexp_2==exp_name2])
amplitudeAtr0.test <- t.test(amplitudeAtr0_g1[pexp_1==exp_name2],amplitudeAtr0_g2[pexp_2==exp_name2])
length_scale.test <- t.test(length_scale_g1[pexp_1==exp_name2],length_scale_g2[pexp_2==exp_name2])

print(exp_name2)
print(density.test$p.value)
print(get_asterisk(density.test))

openpdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name2, "density.pdf"), collapse = "_"))
set1 <- density1$level1[pexp_1==exp_name2]
set2 <- density2$level2[pexp_2==exp_name2]
ylim <- c(0, 1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim,
              names=c('mono', 'dual'), ylab=bquote('Average density ' ~ .(levels2) ~ ' [points/' ~ nm^2 ~']'))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), 
         y1 =  1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE), get_asterisk(density.test) , cex=1)
closepdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name2, "density.pdf"), collapse = "_"))

openpdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name2, "r0.pdf"), collapse = "_"))
set1 <- parameters1$r0[which(pexp_1 == exp_name2)]
set2 <- parameters2$r0[which(pexp_2 == exp_name2)]
ylim <- c(50,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, names=c('mono', 'dual'), ylab=expression(paste('Correlation range ', r[0], ' [nm]')), 
              ylim=ylim)
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE), get_asterisk(parameters_r0.test) , cex=1)
closepdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name2, "r0.pdf"), collapse = "_"))
cat('\t\'statistical_tests_r0.pdf\' created.\n')

openpdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name2, "amplitude.pdf"), collapse = "_"))
set1 <- amplitude_g1[which(pexp_1 == exp_name2)]
set2 <- amplitude_g2[which(pexp_2 == exp_name2)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim, 
              names=c('mono', 'dual'), ylab=expression("Amplitude A"))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name2, "amplitude.pdf"), collapse = "_"))
cat('\t\'statistical_tests_amplitude.pdf\' created.\n')

openpdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name2, "amplitudeAtr0.pdf"), collapse = "_"))
set1 <- amplitudeAtr0_g1[which(pexp_1 == exp_name2)]
set2 <- amplitudeAtr0_g2[which(pexp_2 == exp_name2)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim, 
              names=c('mono', 'dual'), ylab=expression("Correlation at 20 nm [a.u.]"))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitudeAtr0.test) , cex=1)
closepdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name2, "amplitudeAtr0.pdf"), collapse = "_"))
cat('\t\'statistical_tests_amplitudeAtr0.pdf\' created.\n')

openpdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name2, "length_scale.pdf"), collapse = "_"))
set1 <- length_scale_g1[which(pexp_1 == exp_name2)]
set2 <- length_scale_g2[which(pexp_2 == exp_name2)]
ylim <- c(0,1.2*max(set1, set2, na.rm=TRUE))
bp <- boxplot(set1, set2, ylim=ylim,
              names=c('mono', 'dual'), ylab=expression("Length scale [nm]"))
segments(x0 = 1, x1 = 2, y0 = 1.04*max(set1, set2, na.rm=TRUE), y1 = 1.04*max(set1, set2, na.rm=TRUE), col = "black")
text( 1.5 , 1.09*max(set1, set2, na.rm=TRUE),  get_asterisk(amplitude.test) , cex=1)
closepdf(paste(c("MonoDual", "statistical_tests", levels2, exp_name2, "length_scale.pdf"), collapse = "_"))
cat('\t\'statistical_tests_length_scale.pdf\' created.\n')




