#Landscapemodel_2.R
#Script for generating fitness landscape model - 

set.seed(1812)
LEN=150
x <- seq(-5, 5, len = LEN)                          #define resolution of landscape
y <- seq(-5, 5, len = LEN)
Xh <- expand.grid(x = x, y = y)
Xm <- expand.grid(x = x, y = y)

Xh <- transform(Xh, 
                z= 1*dnorm(x,-2,sd = 1.5)* dnorm(y,-2,sd=1.5)+
                  .5*dnorm(x,-1,sd = 1)* dnorm(y,-3,sd=1)+
                  .4*dnorm(x,.2,sd = 1.1)* dnorm(y,0.2,sd=1.2)
                
)

Xm <- transform(Xm, 
                z=1*dnorm(x, 2,sd = 1.5)*dnorm(y,2,sd=1.5)+
                  .5*dnorm(x,3,sd = 1)* dnorm(y,1.5,sd=1)+
                  .3*dnorm(x,0,sd = 1.1)* dnorm(y,0,sd=1.2)
)

Xall <- expand.grid(x = x, y = y)
Xall <- transform(Xall, 
                  z= 1*dnorm(x,-2,sd = 1.5)* dnorm(y,-2,sd=1.5)+
                    .5*dnorm(x,-1,sd = 1)* dnorm(y,-3,sd=1)+
                    .3*dnorm(x,.2,sd = 1.1)* dnorm(y,0.2,sd=1.2)+(
                      1*dnorm(x, 2,sd = 1.5)*dnorm(y,2,sd=1.5)+
                        .5*dnorm(x,3,sd = 1)* dnorm(y,1.5,sd=1)+
                        .3*dnorm(x,0,sd = 1.1)* dnorm(y,0,sd=1.2)
                    )
)

XallC <- expand.grid(x = x, y = y)
XallC <- transform(Xall, 
                   z= 1*dnorm(x,-2,sd = 1.5)* dnorm(y,-2,sd=1.5)+
                     .5*dnorm(x,-1,sd = 1)* dnorm(y,-3,sd=1)+
                     .3*dnorm(x,.2,sd = 1.1)* dnorm(y,0.2,sd=1.2)-
                     1*dnorm(x, 2,sd = 1.5)*dnorm(y,2,sd=1.5)-
                     .5*dnorm(x,3,sd = 1)* dnorm(y,1.5,sd=1)-
                     .3*dnorm(x,0,sd = 1.1)* dnorm(y,0,sd=1.2)
)


mpal.colors <- colorRampPalette( c('lightgrey','lightgrey',"#9E5471","#993366") ) #define gradient colors
hpal.colors <- colorRampPalette( c('lightgrey','lightgrey',"#25858E","#195D60") ) #define gradient colors
SurfaceColors <- SurfaceColors <- colorRampPalette( c(RColorBrewer::brewer.pal(n = 4,"RdYlBu")[4],'lightgrey',RColorBrewer::brewer.pal(n = 4,"RdYlBu")[1]) ) #define gradient colors

nbcol <- 500
mcolor <- mpal.colors(nbcol)
hcolor <- hpal.colors(nbcol)
Ccolor <- SurfaceColors(nbcol)

##HUMAN
z <- matrix(Xh$z, nrow = LEN)
nrz <- nrow(z)
ncz <- ncol(z)

nbcol <- 500
zfaceta <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfaceta, nbcol)

pdf("~/Desktop/landscape_output_OVERLAP.pdf")
persp(x, y, z,xlab = "Sequence Dimension 1",ylab = "Sequence Dimension 2",zlab = "Fitness",
      theta = 50, phi = 30,     #adjust viewing angle
      ltheta = 90,              #adjust lighting angle
      shade = 0.6,d=10,border=NA,lwd = 0.1,col=hcolor[facetcol],r=3,box = F,expand = 0.8)


#####All
z <- matrix(Xall$z, nrow = LEN)
nrz <- nrow(z)
ncz <- ncol(z)

nbcol <- 500
zfacetALL <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
facetcol <- cut(zfaceta, nbcol)

#colors for all
zC <- matrix(XallC$z, nrow = LEN)
nrzC <- nrow(zC)
nczC <- ncol(zC)

nbcolC <- 500
zfacetaC <- zC[-1, -1] + zC[-1, -nczC] + zC[-nrzC, -1] + zC[-nrzC, -nczC]
facetcolC <- cut(zfacetaC, nbcol)


pdf("~/Desktop/landscape_output_OVERLAP.pdf")
persp(x, y, z,xlab = "Sequence Dimension 1",ylab = "Sequence Dimension 2",zlab = "Fitness",
      theta = 50, phi = 30,     #adjust viewing angle
      ltheta = 90,              #adjust lighting angle
      shade = 0.6,d=10,border=NA,lwd = 0.1,col=Ccolor[facetcolC],r=3,box = F,expand = 0.8)

contour(zfaceta,lwd=2,col = "#195D60",drawlabels = F,labcex = 0.01)
contour(zfacetb,lwd = 2,col = "#993366",drawlabels = F,add = T,nlevels = 5)
#contour(zfacetc,col = viridis::viridis(10),drawlabels = F,nlevels = 10)

#image(z,col=color[facetcol])
dev.off()

