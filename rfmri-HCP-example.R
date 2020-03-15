

rm(list=ls())

## PCA fators are saved in the following file
load("rsfmri-example.Rdata");
## Note that the data dimension is dim*length (5*1200) 

source("dpd_VAR_library.R");

# Best VAR order
VAR.best(z, 5)$opt

## Critical value
#n.est = 5*5 + 5*6/2;
#critical = q.Brown.bridge2(n.est, n=2000, r=5000, level=.95, n.cl=10)
## Critical value is given by 16.08714
critical = 16.08714;

## Now find change-points using MDPDE
c0 = score.chol.binary(z, p=1, alpha=0, critical=critical);
c0$khat

c1 = score.chol.binary(z, p=1, alpha=0.1, critical=critical);
c1$khat

c2 = score.chol.binary(z, p=1, alpha=0.2, critical=critical);
c2$khat

c3 = score.chol.binary(z, p=1, alpha=0.3, critical=critical);
c3$khat

c4 = score.chol.binary(z, p=1, alpha=0.4, critical=critical);
c4$khat


## Draw change-points in the paper
## Draw line plot for changes
e1=e2=e3=e4=e0 =list()
e0$khat = c0$khat[c0$khat >1 & c0$khat < 1200];
e1$khat = c1$khat[c1$khat >1 & c1$khat < 1200];
e2$khat = c2$khat[c2$khat >1 & c2$khat < 1200];
e3$khat = c3$khat[c3$khat >1 & c3$khat < 1200];
e4$khat = c4$khat[c4$khat >1 & c4$khat < 1200];

par(mfrow=c(1,1))
par(mar=c(2,0,2,0.1), oma=c(0,3.2,0,0))
par(cex=.5)
title=c("Change-points");
plot(e0$khat, rep(1, length(e0$khat)), type="p", cex=.7, yaxt='n', xlim=c(1, 1200), ylim=c(0.5, 5.5), xlab="", ylab="",main="Change-points");	
abline(h=1, col="grey", lwd=.5)
abline(h=2, col="grey", lwd=.5)
points(e1$khat, rep(2, length(e1$khat)), type="p", cex=.7, xlim=c(1, 1200)); 
abline(h=3, col="grey", lwd=.5)
points(e2$khat, rep(3, length(e2$khat)), type="p", cex=.7, xlim=c(1, 1200)); 
abline(h=4, col="grey", lwd=.5)
points(e3$khat, rep(4, length(e3$khat)), type="p", cex=.7, xlim=c(1, 1200)); 
abline(h=5, col="grey", lwd=.5)
points(e4$khat, rep(5, length(e4$khat)), type="p", cex=.7, xlim=c(1, 1200)); 
abline(v=c(292, 435, 882), lwd=1, lty=2);
axis(side = 2, tck = -.015, at=c(1, 2, 3,4,5), labels = c(expression(paste(alpha, "=0")), 
                                                          expression(paste(alpha, "=.1")),
                                                          expression(paste(alpha, "=.2")),
                                                          expression(paste(alpha, "=.3")),
                                                          expression(paste(alpha, "=.4"))))

