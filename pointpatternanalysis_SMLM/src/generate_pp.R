# create point patterns
# http://mapas.mma.gov.br/i3geo/pacotes/rlib/win/spatstat/html/owin.html

### import libraries
library(spatstat)  # statistical point pattern analysis
library(car) # class. statistics
library(extrafont)  # fonts for latex pdf
loadfonts()
library(knitr) # class. statistics
# library(ggplot2)
graphics.off()

# set working directory
setwd("/home/alba/DropboxCRG/postdoc_CRG/coding/cellviewer/point_pattern_analysis/src") 

source("../src_clean/utilities.R")

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

# exp fitting, estimated parameters
A <- 30
lambda <- 200

# observation window
ntotal <- 1000  # total number of points
area <- 16*10^6  # nm2
rho_av <- ntotal/area  # average density

# clusters
rclusters <- lambda  # average radius clusters [nm]
nclusters <- 2*A*pi*lambda^2*rho_av  # average number of proteins per cluster
rho_cluster <-  nclusters/(pi*rclusters^2)

phi <- rho_cluster/rho_av  # increased density of points in clusters

# simulate
clusters <- 5  # number of clusters
rclusters <- rnorm(clusters, mean = xi_clusters, sd = 1)

background <- rpoint(ntotal, win=owin(c(0,sqrt(area)), c(0,sqrt(area))))
plot(background, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)

# a circle
runifdisc(30)

# # triangle
# X <- rpoint(50, win=triangle(cenx=25, ceny=25, gamma=20, rot=0,height=15))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# 
# # triangle
# X <- rpoint(60, win=triangle(cenx=13, ceny=15, gamma=10, rot=-20,height=3))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# 
# X <- rpoint(60, win=triangle(cenx=20, ceny=20, gamma=20, rot=10,height=3))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# 
# X <- rpoint(60, win=triangle(cenx=30, ceny=30, gamma=10, rot=150,height=6))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# 
# X <- rpoint(60, win=triangle(cenx=37, ceny=37, gamma=40, rot=150,height=2))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# 
# # a square
# X <- rpoint(40, win=quadrilateral(cenx=20, ceny=10, rot=130, height=1, width=1))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# 
# X <- rpoint(60, win=quadrilateral(cenx=35, ceny=15, rot=-10, height=2, width=2))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# 
# X <- rpoint(80, win=quadrilateral(cenx=15, ceny=35, rot=120, height=3, width=3))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# 
# X <- rpoint(60, win=quadrilateral(cenx=15, ceny=35, rot=120, height=3, width=3))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# 
# # a rectange
# X <- rpoint(80, win=quadrilateral(cenx=10, ceny=30, rot=10, height=3, width=1))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# 
# X <- rpoint(80, win=quadrilateral(cenx=20, ceny=30, rot=-50, height=2.5, width=0.8))
# points <- superimpose(points,X)
# plot(points, cex=0.2, main='', pch=16, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# 
# fileName <- 'pp_sample'
# openpdf(paste(fileName, ".pdf", sep = ''))
# plot(points, cex=0.2, main='', pch=16, cols="black", 
#      use.marks=TRUE, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5)
# closepdf(paste(fileName, ".pdf", sep = ''))
# 
# savetext(matrix(c(points$x,points$y), nrow=points$n, ncol=2), paste(fileName, ".txt", sep = ''))
# 
# 
# 
# # pixelate
# im <- pixellate(points, step=0.1)
# im <- as.im(points, dimyx=100)
# openpdf("pppixelated_triangle_square_low.pdf")
# plot(im, main='')
# closepdf("pppixelated_triangle_square_low.pdf")
# savetext(matrix(im$v, nrow=im$dim[1], ncol=im$dim[2]), 'pixelated_triangle_square.txt', ncolumns=im$dim[2])
