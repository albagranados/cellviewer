# # see http://ase.tufts.edu/gsc/gradresources/guidetomixedmodelsinr/mixed%20model%20guide.html
# run after ppanalysis_average (cell_g12: g_12 for each cell in the experiment)

# ------------------------- TOW EXPERIMENTS------
# ----------------------------------------------
exp1_name <- "DMSO"
load("../../data/victoria/Cohesin_CTCF/K4me2SMC1_DMSO.RData")
data.exp1 <- t(sapply(cell_g12, `[[`, 'g_cross'))
r0_exp1 <- cell_r0
rend_exp1 <- which(r_eval>=quantile(r0_exp1, probs=0.5, na.rm=TRUE))[1]

exp2_name <- "ActD"
load("../../data/victoria/Cohesin_CTCF/K4me2SMC1_actd.RData")
data.exp2 <- t(sapply(cell_g12, `[[`, 'g_cross'))
r0_exp2 <- cell_r0
rend_exp2 <- which(r_eval>=quantile(r0_exp2, probs=0.5, na.rm=TRUE))[1]

r_eval <- cell_g12[[1]]$r

data_temp <- data.frame( r = (rbind(data.exp1, data.exp2)))
data_wide <- data.frame(
  id_cell = factor(1:nrow(data_temp)),
  # col.names=c('cell','experiment',’Age8′,’Age10′,’Age1),
  # colClasses=c('factor','factor', rep('numeric', 4))),
  experiment = factor(rep(c(exp1_name, exp2_name), c(nrow(data.exp1),nrow(data.exp2)))),
  data_temp
  )
data_wide$experiment <- relevel(data_wide$experiment, ref=exp1_name)

## first zero-crossing of the cross pair-correlation functions
data_r0 <- data.frame(id_cell = factor(1:nrow(data_temp)),
                           experiment = factor(rep(c(exp1_name, exp2_name), c(nrow(data.exp1),
                                               nrow(data.exp2)))), r0 = c(r0_exp1, r0_exp2))

data_wide[which(data_wide$experiment==exp1_name), rend_exp1:length(data_wide[1,])] <- NA
data_wide[which(data_wide$experiment==exp2_name), rend_exp2:length(data_wide[1,])] <- NA

# for(i in 1:nlevels(data_wide$id_cell)){
#   # print(i); print(data_r0$r0[i]); 
#   if (!is.na(data_r0$r0[i]))
#   {
#     print((2+which(r_eval == data_r0$r0[i])):length(data_wide[i,]))
#     data_wide[i, (2+which(r_eval == data_r0$r0[i])):length(data_wide[i,])] <- NA
#   }
# }

# convert wide-formatted data into long
data <- reshape(data_wide, varying=names(data_wide)[3:ncol(data_wide)] , idvar=c("id_cell", "experiment"), 
                     direction="long", timevar='r', times=r_eval, v.names='g12')
data$experiment <- relevel(data$experiment,ref=exp1_name)

# write.table(data, "../../data/victoria/Cohesin_CTCF/K4me2SMC1_DMSOActD_long.txt", sep="\t")
save(data_wide, data, data_r0, r_eval, file="../../data/victoria/Cohesin_CTCF/K4me2SMC1.RData")


# ------------------------- ONE EXPERIMENT------
# ----------------------------------------------
load("../")
exp1_name <- "exp1"
data.exp1 <- t(sapply(cell_g12, `[[`, 'g_cross'))
r0_exp1 <- cell_r0
rend_exp1 <- which(r_eval>=quantile(r0_exp1, probs=0.5, na.rm=TRUE))[1]  # probs = 0.5

r_eval <- cell_g12[[1]]$r

data_temp <- data.frame( r = data.exp1)
data_wide <- data.frame(
  id_cell = factor(1:nrow(data_temp)),
  experiment = factor(rep(exp1_name, nrow(data.exp1))),
  data_temp
)
data_wide$experiment <- relevel(data_wide$experiment, ref=exp1_name)

## first zero-crossing of the cross pair-correlation functions
data_r0 <- data.frame(id_cell = factor(1:nrow(data_temp)),
                      experiment = factor(rep(exp1_name, nrow(data.exp1))), r0 = r0_exp1)

# data_wide[which(data_wide$experiment==exp1_name), rend_exp1:length(data_wide[1,])] <- NA
for (i in which(data_wide$experiment==exp1_name)){
  data_wide[which(data_wide$experiment==exp1_name)[i], which(r_eval==cell_r0[i]):length(data_wide[1,])] <- NA
}

# convert wide-formatted data into long
data <- reshape(data_wide, varying=names(data_wide)[3:ncol(data_wide)] , idvar=c("id_cell", "experiment"), 
                direction="long", timevar='r', times=r_eval, v.names='g12')
data$experiment <- relevel(data$experiment,ref=exp1_name)

save(data_wide, data, data_r0, r_eval, file="../../data/test/synthetic_pp/test.RData")
