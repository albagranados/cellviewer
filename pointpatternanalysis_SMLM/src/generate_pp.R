# create point patterns
# http://mapas.mma.gov.br/i3geo/pacotes/rlib/win/spatstat/html/owin.html

### import libraries
library(spatstat)  # statistical point pattern analysis
library(car) # class. statistics
library(MASS)
library(extrafont)  # fonts for latex pdf
library(RandomFields)
loadfonts()
library(knitr) # class. statistics
# library(ggplot2)
graphics.off()

# set working directory
setwd("/home/alba/ownCloud/postdoc_CRG/coding/github/cellviewer/pointpatternanalysis_SMLM/src") 

source("utilities.R")

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
    p[,1] <- p[,1] + centre[1]
    p[,2] <- p[,2] + centre[2]
    points <- as.ppp(p, c(centre[1]-radius, centre[1]+radius, centre[2]-radius, centre[2]+radius))
    return(points)
  }
  if (distribution=='gaussian_truncated'){
    p <- mvrnorm(n=3*n, mu=c(0,0), Sigma=matrix(c(sig^2,0,0,sig^2), 2, 2))
    out <- which((abs(p[,1])>radius) | (abs(p[,2])>radius))
    print(out)
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

## input parameters
num_realizations <- 10
cluster_type <- 'thomas'
# exp fitting, estimated parameters
A <- 10
lambda <- 100
# background point pattern
ntotal <- 1 # 1000  # total number of points
area <- 20*20*160^2 # 20*20*160^2  # nm2
rho_av <- ntotal/area  # average density
# number of clusters
num_clusters <- 10  # number of clusters

# cluster point pattern
rcluster <- lambda  # average radius clusters [nm]
Nclusters <- 50# round(2*A*pi*lambda^2*rho_av)  # average number of proteins per cluster
rho_cluster <-  Nclusters/(pi*rcluster^2)
phi <- rho_cluster/rho_av  # increased density of points in clusters

# simulate
for (j in 1:num_realizations){
  # # Mosaic pp
  # grid <- dirichlet(runifpoint(3000, win = square(20)))
  # logLambda <- rMosaicField(grid, rnorm, dimyx=512, rgenargs=list(mean=1, sd=1))  # mean points per pixel
  # Lambda <- exp(logLambda)
  # points <- rpoispp(Lambda)

  radii <- rnorm(num_clusters, mean = rcluster, sd = 0)
  centreclusters <- cbind(runif(num_clusters, min = 3*max(radii), max = sqrt(area)-3*max(radii)),
                          runif(num_clusters, min = 3*max(radii), max = sqrt(area)-3*max(radii)))
  pp_background <- rpoint(ntotal, win=owin(c(0,sqrt(area)), c(0,sqrt(area))))

  points <- pp_background
  for (i in 1:num_clusters){
    pp_cluster <- circle(cluster_type, n=Nclusters, radius=rcluster, centre=centreclusters[i,], sig=rcluster/3)
    points <- superimpose(pp_cluster,points)
  }
  plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  
  fileName <- paste(c('synthetic_', cluster_type, '_A', A, '_lambda', lambda, '_nt', ntotal, '_numclusters', 
                      num_clusters, '_', j), collapse='')
  openpdf(paste(fileName, ".pdf", sep = ''))
  plot(points, cex=0.13, main='', pch=16, cols="black", cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
  closepdf(paste(fileName, ".pdf", sep = ''))
  savetext(matrix(c(points$x,points$y), nrow=points$n, ncol=2), paste(fileName, ".txt", sep = ''))
}


# # Tools in spatstat can be used to simulate such models by hand. Figure 14.29 shows examples
# # of linked Cox and balanced Cox processes, generated by
# P <- dirichlet(runifpoint(100, win = square(2)))
# logLambda <- rMosaicField(P, rnorm, dimyx=512, rgenargs=list(mean=4, sd=1))
# Lambda <- exp(logLambda)
# X <- rpoispp(Lambda)
# Xlinked <- X %mark% factor(sample(c("a","b"), npoints(X), replace=TRUE, prob=c(2,3)/5))
# Y <- rpoispp(100)
# Xbalanced <- Y %mark% factor(ifelse(Lambda[Y] > exp(4), "a", "b"))
# # and between them a multitype log-Gaussian Cox process, generated by
# X1 <- rLGCP("exp", mu=4, var=0.2, scale=.1)
# Z1 <- log(attr(X1, "Lambda"))
# Z2 <- log(attr(rLGCP("exp", 4, var=0.1, scale=.1), "Lambda"))
# Lam1 <- exp(Z1)
# Lam2 <- exp((Z1 + Z2)/2)
# Xlgcp <- superimpose(a=rpoispp(Lam1), b=rpoispp(Lam2))

# # # Thomas pp
# points <- rThomas(kappa=num_clusters/area, scale=round(rcluster/3), mu=Nclusters, 
#                   win = owin(c(0,sqrt(area)), c(0,sqrt(area))),
#                   drop=TRUE, saveLambda=FALSE, saveparents=TRUE)

# # triangle
# X <- rpoint(50, win=triangle(cenx=25, ceny=25, gamma=20, rot=0,height=15))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)

# # a square
# X <- rpoint(40, win=quadrilateral(cenx=20, ceny=10, rot=130, height=1, width=1))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)

# # a rectange
# X <- rpoint(80, win=quadrilateral(cenx=10, ceny=30, rot=10, height=3, width=1))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)

# 
# # pixelate
# im <- pixellate(points, step=0.1)
# im <- as.im(points, dimyx=100)
# openpdf("pppixelated_triangle_square_low.pdf")
# plot(im, main='')
# closepdf("pppixelated_triangle_square_low.pdf")
# savetext(matrix(im$v, nrow=im$dim[1], ncol=im$dim[2]), 'pixelated_triangle_square.txt', ncolumns=im$dim[2])
