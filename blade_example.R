### reproducing Section 5.3: Thermal Stress Analysis of Jet Engine Turbine Blade ###
library(matlabr)
library(randtoolbox)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(lhs)
source("GP.R")
source("matern.R")
source("stacking_design.R")

# a function for FEM simulations through matlab
f <- function(X, l, MM, tt, return.time=FALSE){
  if(is.null(dim(X))) n <- length(X) else n <- nrow(X)
  dfr <- data.frame(X*0.5+0.25, rep(MM/tt^l, n)) # scale X to [-1,1]
  write.csv(dfr, "Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("Rmatlab_files/SolveJetBlade.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  dfm <- read.table("Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  if(sum(abs(dfm[,1:3]-dfr))>1e-5){print("wrong data!"); break}
  
  if(return.time) return(list(value=dfm[,4],time=dfm[,5])) else return(dfm[,4])
}

cost <- c(0.75,1.07,2.13,11.51) # this should be NULL, but for the sake of reproducibilty, I used the saved costs
d <- 2              # d: dimension of X (scalar)
n.init <- 5*d
alpha <- NULL
tt <- 2             # mesh = 0.4/(tt^l), lowest level = 0.2
MM <- 0.1
Lmax <- 5
k <- NULL           # k: the parameter for the Matern kernel function (NULL means it'll be estimated by LOOCV) 
xnew <- expand.grid(seq(0,1,0.02),seq(0,1,0.02)) # xnew.ori: test data (original scale)
xnew.ori <- xnew*0.5+0.25
n.max <- 300        # the maximum number of sample size
log.fg <- TRUE

set.seed(1)
X.test <- maximinLHS(20,2)
y.test <- f(X.test,l=5, MM, tt, return.time=TRUE)
# load("2d_jet_blade_true.Rdata")
# y.test <- true.save

ML.out <- stacking_design(cost, d, n.init, epsilon=5, alpha, 
                          tt, MM, Lmax, f, k, xnew, n.max, 
                          log.fg=TRUE, norm="L2", model.save=TRUE, save.n = TRUE,
                          init=rep(2.5,d), lower=0.5, upper=5)

P <- matrix(0, nrow=nrow(X.test), ncol=ML.out$L)
for(l in 1:ML.out$L){
  out <- pred.GP(ML.out$model[[l]], X.test)
  P[,l] <- out$mu
}
y.pred <- rowSums(P)

# empirical L2 error
sqrt(mean((y.pred - y.test$value)^2))

# summarize sample sizes
sample.size.df <- data.frame(matrix(0,ncol=4,nrow=16))
colnames(sample.size.df) <- c("name", "stage", "value", "total")
sample.size.df[1,4] <- ML.out$n.save[[1]]
sample.size.df[1,3] <- ML.out$n.save[[1]]
sample.size.df[5:6,4] <- ML.out$n.save[[2]]
sample.size.df[5:6,3] <- ML.out$n.save[[2]] - c(ML.out$n.save[[1]],0)
sample.size.df[9:11,4] <- ML.out$n.save[[3]]
sample.size.df[9:11,3] <- ML.out$n.save[[3]] - c(ML.out$n.save[[2]],0)
sample.size.df[13:16,4] <- ML.out$n.save[[4]]
sample.size.df[13:16,3] <- ML.out$n.save[[4]] - c(ML.out$n.save[[3]],0)

sample.size.df[,1] <- rep(paste0("n",1:4),4)
sample.size.df[,1] <- factor(sample.size.df[,1], levels = paste0("n",1:4))
sample.size.df[,2] <- rep(1:4,each=4)
sample.size.df[,2] <- factor(sample.size.df[,2], levels = c(4, 3, 2, 1))

# Figrue 1
pdf("blade_design.pdf", width=8.5, height=2.5)
layout(matrix(c(rep(1:4,each=2),5),nrow=1))
par(mar=c(2, 2, 4, 0.2))
plot(sobol(sample.size.df[1,4], 2, scrambling = 2),xlab="",ylab="",col=1, yaxt = "n", xaxt = "n",cex=2, pch=0, main=expression(italic(L)==1))
title(xlab="x1 (pressure side)", ylab="x2 (suction side)", line=0.5, cex.lab=1.2)
plot(sobol(sample.size.df[5,4], 2, scrambling = 2),xlab="",ylab="",col=1, yaxt = "n", xaxt = "n",cex=2, pch=0, main=expression(italic(L)==2))
points(sobol(sample.size.df[6,4], 2, scrambling = 2),col=2,pch=1,cex=3)
title(xlab="x1 (pressure side)", ylab="x2 (suction side)", line=0.5, cex.lab=1.2)
plot(sobol(sample.size.df[9,4], 2, scrambling = 2),xlab="",ylab="",col=1, yaxt = "n", xaxt = "n",cex=2, pch=0, main=expression(italic(L)==3))
points(sobol(sample.size.df[10,4], 2, scrambling = 2),col=2,pch=1,cex=3)
points(sobol(sample.size.df[11,4], 2, scrambling = 2),col=4,pch=2,cex=2)
title(xlab="x1 (pressure side)", ylab="x2 (suction side)", line=0.5, cex.lab=1.2)
plot(sobol(sample.size.df[13,4], 2, scrambling = 2),xlab="",ylab="",col=1, yaxt = "n", xaxt = "n",cex=2, pch=0, main=expression(italic(L)==4))
points(sobol(sample.size.df[14,4], 2, scrambling = 2),col=2,pch=1,cex=3)
points(sobol(sample.size.df[15,4], 2, scrambling = 2),col=4,pch=2,cex=2)
points(sobol(sample.size.df[16,4], 2, scrambling = 2),col=5,pch=4,cex=2)
title(xlab="x1 (pressure side)", ylab="x2 (suction side)", line=0.5, cex.lab=1.2)

par(mar=c(0, 0, 0, 0))
plot(1,1,xaxt="n",type="n",yaxt="n", bty="n")
legend("center", 
       legend = c(expression(italic(l)==1), expression(italic(l)==2), expression(italic(l)==3), expression(italic(l)==4)), 
       pch=c(0,1,2,4), col=c(1,2,4,5), cex=1.5, pt.cex=c(2,3,2,2), bty = "n", y.intersp=2)
dev.off()

pdf("blade_design_size.pdf", width=10, height=2)
g1 <- ggplot(sample.size.df[1:4,], aes(x = name,y = value ,fill = stage)) +
  geom_bar(stat = "identity") +
  labs(y = "sample size", x = "") + 
  ylim(c(0,43)) + 
  scale_x_discrete(labels=c(expression(n[1]), expression(n[2]), expression(n[3]), expression(n[4]))) +
  scale_fill_manual(breaks=c(1:4), values = c("lightblue1","lightblue2", "lightblue3", "lightblue4")) +
  coord_flip() + theme( panel.background = element_blank(), legend.position = "none")

g1 <- g1 + geom_text(aes(name, total + 3, label = total, fill = NULL), size = 4, data = sample.size.df[sample.size.df$stage==1 & sample.size.df$total > 0,])
g1 <- g1 + theme(axis.text.y = element_text(size = 12))

g2 <- ggplot(sample.size.df[1:8,], aes(x = name,y = value ,fill = stage)) +
  geom_bar(stat = "identity") +
  labs(y = "sample size", x = "") + 
  ylim(c(0,43)) + 
  scale_x_discrete(labels=c(expression(n[1]), expression(n[2]), expression(n[3]), expression(n[4]))) +
  scale_fill_manual(breaks=c(1:4), values = c("lightblue1","lightblue2", "lightblue3", "lightblue4")) +
  coord_flip() + theme( panel.background = element_blank(), legend.position = "none")

g2 <- g2 + geom_text(aes(name, total + 3, label = total, fill = NULL), size = 4,
                     data = sample.size.df[sample.size.df$stage==2 & sample.size.df$total > 0,])
g2 <- g2 + theme(axis.text.y = element_text(size = 12))

g3 <- ggplot(sample.size.df[1:12,], aes(x = name,y = value ,fill = stage)) +
  geom_bar(stat = "identity") +
  labs(y = "sample size", x = "") + 
  ylim(c(0,43)) + 
  scale_x_discrete(labels=c(expression(n[1]), expression(n[2]), expression(n[3]), expression(n[4]))) +
  scale_fill_manual(breaks=c(1:4), values = c("lightblue1","lightblue2", "lightblue3", "lightblue4")) +
  coord_flip() + theme( panel.background = element_blank(), legend.position = "none")

g3 <- g3 + geom_text(aes(name, total + 3, label = total, fill = NULL), size = 4,
                     data = sample.size.df[sample.size.df$stage==3 & sample.size.df$total > 0,])
g3 <- g3 + theme(axis.text.y = element_text(size = 12))

g4 <- ggplot(sample.size.df, aes(x = name,y = value ,fill = stage)) +
  geom_bar(stat = "identity") +
  labs(y = "sample size", x = "") + 
  ylim(c(0,43)) + 
  scale_x_discrete(labels=c(expression(n[1]), expression(n[2]), expression(n[3]), expression(n[4]))) +
  scale_fill_manual(breaks=c(1:4), values = c("lightblue1","lightblue2", "lightblue3", "lightblue4")) +
  coord_flip() + theme( panel.background = element_blank(), legend.position = "none")

g4 <- g4 + geom_text(aes(name, total + c(3,3,3,3), label = total, fill = NULL), size = 4,
                     data = sample.size.df[sample.size.df$stage==4,])
g4 <- g4 + theme(axis.text.y = element_text(size = 12))

grid.arrange(g1, g2, g3, g4, ncol = 4)
dev.off()


# Figure 12 (no longer needed)
# pdf("alpha_estimation_blade.pdf", width=6, height=2.5)
# z.val <- c(ML.out$z[[2]][1:ML.out$n.save[[3]][2]], ML.out$z[[3]][1:ML.out$n.save[[3]][3]])
# df_alpha <- data.frame(mesh_l=rep(1/tt^(2:3), ML.out$n.save[[3]][2:3]),
#                        z=abs(z.val))
# lm.fit <- lm(log(z) ~ log(mesh_l), data=df_alpha)
# 
# g1 <- ggplot(df_alpha, aes(x=log(mesh_l), y=log(z), group=log(mesh_l))) + 
#   geom_boxplot(fill='#A4A4A4', color="black", width=0.25)+
#   xlab(expression(log(italic(T^-l))))+ylab(bquote("log|" ~ f[italic(l)] ~ "(x) -" ~ f[italic(l)-1] ~ "(x)|"))+#ylab(bquote("log|" ~ z[italic(l)] ~  "|"))+
#   theme_bw()+ggtitle(expression(italic(L)==3)) + theme(plot.title = element_text(hjust = 0.5))
# 
# g1 <- g1 + geom_abline(intercept = lm.fit$coefficients[1], 
#                        slope = lm.fit$coefficients[2], color = "red", linetype = "dashed") 
# 
# z.val <- unlist(ML.out$z[2:4])
# df_alpha <- data.frame(mesh_l=rep(1/tt^(2:4), ML.out$n[2:4]),
#                        z=abs(z.val))
# lm.fit <- lm(log(z) ~ log(mesh_l), data=df_alpha)
# 
# g2 <- ggplot(df_alpha, aes(x=log(mesh_l), y=log(z), group=log(mesh_l))) + 
#   geom_boxplot(fill='#A4A4A4', color="black", width=0.25)+
#   xlab(expression(log(italic(T^-l))))+ylab(bquote("log|" ~ f[italic(l)] ~ "(x) -" ~ f[italic(l)-1] ~ "(x)|"))+#ylab(bquote("log|" ~ z[italic(l)] ~  "|"))+
#   theme_bw() + ggtitle(expression(italic(L)==4)) + theme(plot.title = element_text(hjust = 0.5))
# 
# g2 <- g2 + geom_abline(intercept = lm.fit$coefficients[1], 
#                        slope = lm.fit$coefficients[2], color = "red", linetype = "dashed") 
# 
# grid.arrange(g1, g2, ncol = 2)
# dev.off()


# Figure 13
pdf("blade_prediction.pdf", width=7.5, height = 3.3)
par(mfrow=c(1,2))
par(mar=c(2,1,0,0))
pmat <- persp(seq(0,1,0.02)*0.5+0.25, seq(0,1,0.02)*0.5+0.25, matrix(ML.out$mean, 51, 51), col="lightblue", 
              theta = 300, phi = 20, ticktype = "detailed", 
              xlab="", ylab="", zlab="f(x1,x2)",zlim=c(0,50))

label.pos <- trans3d(0.5, 1.1, -50, pmat)
text(label.pos$x, label.pos$y, labels="suction (MPa)", adj=c(0, NA), srt=-25, cex=0.8)

label.pos <- trans3d(4, 1, -400, pmat)
text(label.pos$x, label.pos$y, labels="pressure (MPa)", adj=c(0, NA), srt=62, cex=0.8)

mypoints <- trans3d(X.test[,1]*0.5+0.25,X.test[,2]*0.5+0.25,y.test$value, pmat = pmat)
points(mypoints, pch = 16,col = 2)
par(mar=c(4,4,1,1.5))
fields::image.plot(seq(0,1,0.02)*0.5+0.25, seq(0,1,0.02)*0.5+0.25, 
                   matrix(rowSums(sqrt(t(t(ML.out$V[,1:ML.out$L])*ML.out$n[1:ML.out$L])))+ML.out$B[,ML.out$L], 51, 51),
                   xlab="pressure (MPa)", ylab="suction (MPa)", 
                   legend.mar=6.5,
                   col = hcl.colors(12, "YlOrRd", rev = TRUE))
dev.off()


