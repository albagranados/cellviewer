# read: http://rpsychologist.com/r-guide-longitudinal-lme-lmer; 
# http://tutorials.iq.harvard.edu/R/Rstatistics/Rstatistics.html#orgheadline36
# Statistical analysis of the pair cross-correlation functions. 
# Run after generate_longitudinal_dataset.R

# rm(list=ls()) # clean the workspace
library(lattice) # for Lattice graphics
trellis.device(color=FALSE) # to get black-and-white figures
library(lme4)
library(nlme)
library(MASS)
library(car)
library(texreg)  # latex table results model
library(extrafont)
loadfonts()
setwd("/home/alba/Dropbox (CRG)/postdoc_CRG/coding/cellviewer/point_pattern_analysis/src")

# first: run generate_longitudinal_dataset.R. Then, load:
load("../../data/victoria/Cohesin_CTCF/K4me2SMC1.RData")
load("../../data/victoria/Cohesin_CTCF/SMC3SMC1.RData")
load("../../data/victoria/Cohesin_CTCF/CTCFSMC1.RData")

attach(data)

exp1_name <- "DMSO"; exp2_name <- "ActD"; unit_size <- 160
exp1_name <- "exp1"; unit_size <- 1


# --------------- primary analysis ----
# ---------------------------------------

example.lm <- lm(log(g12-1) ~ r, data=data[which(data$id_cell == 3), ])
with(data[which(data$id_cell == 3), ], plot(r, g12)); abline(1,0)
with(data[which(data$id_cell == 3), ], plot(r, log(g12-1))); abline(coef(example.lm))

xyplot(g12 ~ r|experiment, data, type=c("p", "r"), ylab = 'log(y-1)', trip=FALSE) 
xyplot(g12 ~ r |id_cell, data, type=c("p", "r"), trip=FALSE) 
exp1.lm <- lm(g12 ~ r, data=data[which(data$experiment == exp1_name), ]) 
with(data[which(data$experiment == exp1_name), ], plot(r, g12, ylab='pair correlation', pch=20)); abline(coef(exp1.lm))
exp2.lm <- lm(g12 ~ r, data=data[which(data$experiment == exp2_name), ]) 
with(data[which(data$experiment == exp2_name), ], plot(r, g12, ylab='pair correlation')); abline(coef(exp2.lm))

## regression model
# linear fitting
g12 <- g12
# exponential fitting
data$g12 <- log(data$g12-1)
# power law fitting
data$g12 <- log(data$g12-1)
data$r <- log(data$r)

attach(data)

## response feature analysis (Faraway, 2006, p.206)

# fit ordinary regressions separately to each group
slopes <- numeric(nlevels(data$id_cell)); intercepts <- numeric(nlevels(data$id_cell))
ratio <- numeric(nlevels(data$id_cell)); area <- numeric(nlevels(data$id_cell))
# pdf("plot.pdf", family="CM Roman")
plot(0, 0, ylim=c(log(1.1-1),log(5-1)), xlim=c(0, r_eval[length(r_eval)]),
     xlab='r', ylab="correlation function", cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
for(i in 1:nlevels(data$id_cell)){
  lmod <- lm(g12 ~ r, subset=(id_cell==i), data=data)
  intercepts[i] <- exp(coef(lmod)[1])
  slopes[i] <- -unit_size*coef(lmod)[2]^(-1)
  ratio[i] <- intercepts[i]/slopes[i]
  # area[i] <- -coef(lmod)[2]*data_r0$r0[which(data_r0$id_cell==i)]
  area[i] <- intercepts[i]*(-slopes[i])*(exp(coef(lmod)[2]*data_r0$r0[which(data_r0$id_cell==i)])-
                                          intercepts[i])+r_eval[1]
  # intercepts[i] <- coef(lmod)[1]
  # slopes[i] <- coef(lmod)[2]
  # area[i] <- intercepts[i]*(1-intercepts[i])/slopes[i]
  cat('r0 = ', data_r0$r0[which(data_r0$id_cell==i)], '\t a = ', intercepts[i], '\t slopes', slopes[i],
      '\t area = ', area[i], '\n')
  
  with(data[which(data$id_cell == i), ], lines(r, g12, type="p", pch=20, col="#A8A8A8")); 
  #abline(coef(lmod))
  Sys.sleep(0.05) 
}
# dev.off(); embed_fonts("plot.pdf", outfile="plot.pdf")

pexp <- data$experiment[match(1:nlevels(data$id_cell),data$id_cell)]
plot(intercepts, slopes, xlab="Intercept", ylab="Slope")
pdf("plot.pdf", family="CM Roman")
old <- par(mfrow=c(2, 2))
boxplot(split(data_r0$r0,pexp), main='r0') #, ylim=c(70,320))
boxplot(split(intercepts,pexp), main="amplitude", names=c(exp1_name, exp2_name), ylim=c(0,15))
boxplot(split(slopes,pexp), main="length scale [nm]", names=c(exp1_name, exp2_name), ylim=c(10,100))
# boxplot(split(ratio,pexp), main="A/length scale", names=c(exp1_name, exp2_name))
boxplot(split(area,pexp), main="area", names=c(exp1_name, exp2_name))
par(old)
dev.off(); embed_fonts("plot.pdf", outfile="plot.pdf")

t.test(data_r0$r0[which(data_r0$experiment == exp1_name)],
       data_r0$r0[which(data_r0$experiment == exp2_name)]) 
t.test(intercepts[pexp==exp1_name],intercepts[pexp==exp2_name])
t.test(slopes[pexp==exp1_name],slopes[pexp==exp2_name])  
t.test(ratio[pexp==exp1_name],ratio[pexp==exp2_name])
t.test(area[pexp==exp1_name],area[pexp==exp2_name])
#  the t-statistic for testing whether the corresponding regression coefficient is different from 0.
# OBS: It requires choosing an important characteristic. We have chosen two here: the slope and the 
# intercept. For many datasets, this is not an easy choice and at least some information is lost by doing this.


# --------------- secondary analysis ----
# ---------------------------------------

# histograms of intercepts/slope distributions
old <- par(mfrow=c(1, 2))
low = min(intercepts[pexp==exp1_name], intercepts[pexp==exp2_name])
high = max(intercepts[pexp==exp1_name], intercepts[pexp==exp2_name])
breaks = seq(low, high, l=10)
hist(intercepts[pexp==exp1_name], col=rgb(0.1,0.1,0.1,0.5), xlab='intercepts', main='', breaks=breaks)
hist(intercepts[pexp==exp2_name], col=rgb(0.8,0.8,0.8,0.5), breaks=breaks, add=T); box()
low = min(slopes[pexp==exp1_name], slopes[pexp==exp2_name])
high = max(slopes[pexp==exp1_name], slopes[pexp==exp2_name])
breaks = seq(low, high, l=10)
hist(slopes[pexp==exp1_name], col=rgb(0.1,0.1,0.1,0.5), xlab='slopes', main='', breaks=breaks)
hist(slopes[pexp==exp2_name], col=rgb(0.8,0.8,0.8,0.5), breaks=breaks, add=T); box()
par(old)

# breaks=seq(min(slopes[pexp==exp2_name]), max(slopes[pexp==exp2_name]), 

# idem:
exp1.list <- lmList(g12 ~ r | id_cell, data=data[pexp==exp1_name,], na.action=na.omit)
boxplot(coef(exp1.list)[,1],coef(exp1.list)[,2], names=c("intercept", "slope"), main=exp1_name)
exp2.list <- lmList(g12 ~ r | id_cell, data=data[pexp==exp2_name,], na.action=na.omit)
boxplot(coef(exp2.list)[,1],coef(exp2.list)[,2], names=c("intercept", "slope"), main=exp2_name)
plot(intervals(exp1.list), main=exp1_name)
plot(intervals(exp2.list), main=exp2_name)
# fitting a linear model to the observations in each group, returning a list of linear-model objects

## violations linear model
# basic residual plots
par(mfrow = c(2, 2))
plot(exp1.lm, pch = 23)
par(mfrow = c(2, 2))
plot(exp2.lm, pch = 23)
plot(exp1.lm)
# bwplot(pexp[pexp==exp1_name] ~ residuals(exp1.lm))
# We now plot the standardized residual against the observed values of the variable r.
# Leverage (or influential) points and outliers I (Leverage points are those which have great influence on the fitted
# model, that is, those whose x-value is distant from the other x-values.)
example.lm <- lm(log(g12-1) ~ r, subset=(id_cell==1), data=data)
plot(r_eval, rstandard(example.lm), ylab="Standardized Residuals", xlab="r") 
abline(0, 0) 
plot ( resid ( example. lm ) ~ fitted ( example. lm ) )# Residuals vs. fitted values
plot(rstandard(example.lm)) # Good leverage points have their standardized residuals within the interval [−2, 2]
plot(example.lm) # Diagnostic plots
plot(g12 ~ r, subset=(id_cell==1), data=data)
# Remove invalid data points. Good leverage points have their standardized residuals within
# the interval [−2, 2]

## models (see Stephen_Mbunai_Sonal_Nagda):
# First    a  generalised  linear  model  (fixed  effects  model)  was   fitted  to  check  the  
# significance  of  each  of  the  fixed   effects
# OBS: The  r  was  highly  significant.    With  every  increase, there would be a decrease 
# of 0.38 (±0.008) g12's. ActD  had  a  significantly  higher  g12  by  0.06(±0.01) g12's  than  
# exp1. 
#  In contrast, if the observations in a \group" represent longitudinal data on a single individual,
# then the structure of the s may be speci ed to capture autocorrelation among the errors, as
# is common in observations collected over time.  The
# lme function in the nlme package can handle autocorrelated and heteroscedastic error. 
# Intercept varies by id_cell and r.
# random effects: how much variability there is due to id_cell and r (the two random effects); 
# 'Residual' which stands for the variability that’s ndata=data, ot due to either r or id_cell (= \epsilon). 
# fixed effects: t value is simply the estimate divided by the standard err. For the effect of one
# of the independent variables.  It’s  the verage of our data for one of the values of this indep. variable

## model
data.model0 <- lmer(g12 ~ 1 + 1|experiment, data=data, REML = FALSE) 
summary(data.model0)
# OBS: a null model (i.e., a model with a random effects grouping structure, but no fixed-effects predictors

## model
data.model1 <- lmer(g12 ~ r + 1|experiment, data=data) 
summary(data.model1)
data.model1 <- lme(g12 ~ 1 + r, random =  ~ 1|experiment, data=data, na.action=na.omit)
summary(data.model1)
data.model1ML <- lme(g12 ~ 1 + r, random =  ~ 1|experiment, data=data, na.action=na.omit, method='ML')
summary(data.model1)

fixef(data.model1)
ranef(data.model1)

data.lm <- lm(g12 ~ 1 + r, data = data, na.action=na.omit)
anova(data.lm, data.model1)

## model
data.model2 <- lme(g12 ~ 1 + r, random = ~ r | experiment, data=data, na.action=na.omit)
summary(data.model2)

## model
data.model3 <- lmer(g12 ~ r + (r|experiment), data=data, REML = FALSE) 
data.model3 <- lme(g12 ~ 1 + r, random = ~ 1 + r | experiment, data=data, na.action=na.omit)
summary(data.model3)
# OBS: in addition to estimating the distribution of intercepts across experiments, we also estimate 
# the  distribution of the slope of g12 on r. 

anova(data.model1, data.model3)
# p- values for mixed models aren’t as straightforward as they are for the linear model. The  logic  
# of  the  likelihood  ratio  test  is  to  compare  the  likelihood  of  two  models  with  each  other.  
# First,  the  model  without the factor  that  you’re  interested  in  (the null  model),  then  the
# model with the  factor  that  you’re  interested
# REML=FALSE: It is necessary to do this when you compare models using the likelihood ratio test
# experiment affected g12 (x2=, p=) lowering it by ...\pm std error

## model
data.lme <- lme(g12 ~ 1 + r, random = ~ 1 | experiment,
                data = data, method = "ML")
summary(data.lme)data.model5ML <- lmer(g12 ~ 1 + r + experiment + 1|id_cell, data=data, na.action=na.omit, REML=FALSE)
data.lm <- lm(g12 ~ 1 + r, data = data)
summary(data.lm)
# OBS: The parameter estimates for the two models seem quite close, however:
#   • the coefficient of IQ differs by about 3 standard errors between the two models
#   • the standard error of the intercept is twice as large in lang.lme –
# the ordinary regression analysis produces an over-optimistic impression of the precision of this estimate.

## model
data.model4 <- lmer(g12 ~ r + 1|id_cell, data=data, na.action=na.omit)
data.model4MLcor <- lme(g12 ~ r, random = ~ 1|id_cell, data=data, na.action=na.omit,method='ML', 
                           correlation=corAR1(value=0))
summary(data.model4MLcor)
summary(data.model4)

## model
data.model5 <- lmer(g12 ~ r + experiment + 1|id_cell, data=data, na.action=na.omit)
data.model5MLcor <- lme(g12 ~ r + experiment, random = ~ 1|id_cell, data=data, na.action=na.omit,
                        method='ML', correlation=corAR1(value=0))
summary(data.model5MLcor)
summary(data.model5)

anova(data.model4, data.model5)
anova(data.model4MLcor, data.model5MLcor)

## model
data.model6 <- lmer(g12 ~ r + (1+r|id_cell), data=data, na.action=na.omit)
data.model7 <- lmer(g12 ~ r + experiment + (1+r|id_cell), data=data, na.action=na.omit)
summary(data.model7)
data.model7 <- lme(g12 ~ 1 + r + experiment, random = ~ 1+r|id_cell, data=data, na.action=na.omit)
summary(data.model7)
fixef(data.model7)
ranef(data.model7)

summary(lme(g12 ~ 1 + r + experiment, random = ~ 1|id_cell, data=data, 
                          na.action=na.omit))
summary(lme(g12 ~ 1 + r + experiment, random = ~ 1|id_cell, data=data, 
                          na.action=na.omit, correlation=corAR1(value=0.8)))

anova(data.model7, data.model5)
anova(data.model4, data.model5, data.model7)
anova(data.model6, data.model7)

## model
data.model8ML <- lmer(g12 ~ r*experiment + (1 + r | id_cell), data=data, na.action=na.omit, REML=FALSE)
summary(data.model8ML)
anova(data.model7ML, data.model8ML)
plot(data.model8ML)

## examine residual
plot(data.model3, experiment ~ resid(.), abline = 0, ylab = "experiment")
qqnorm(data.model3, ~ resid(., type = "p"))

res_lme=residuals(model1)
plot(res_lme)
qqnorm(res_lme)
qqline(res_lme)
plot(model1)

# each r and each id_cell is assigned a different intercept (of course because model |)
# we used a random intercept model: whatever the effect of experiment is, it’s going to be the same for 
# all id_cell and r. Otherwise we need a random slop model. They  are  also  allowed  to  have  
# different  slopes  for  the  effect  of  experiment. 

## Testing the immediate post-intervention difference
with(data, t.test(g_12.2~experiment))
# Obs: The t.test() function is a standard one for one-sample, paired, and two-sample t-tests. 
#     By assigning a formula-alike measure.2~tx argument into the function, R will
#     assume it’s a two-sample t-test problem and calculate the result for us. Ooops, no significant 
#     difference were detected at the immediate post-intervention measurement.
#     This data needs a more in-depth analysis, and we would provide a walkthrough step-by-step.

## output model to latex
# > texreg(data.model,
#                caption="The Importance of Clustering Standard Errors",
#                dcolumn=FALSE,
#                model.names=c("M1"),
#                override.se=list(summary(data.model)$coef[,2]),
#                override.pval=list(summary(data.model)$coef[,4]))
