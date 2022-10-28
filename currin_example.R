### reproducing Section 5.1: Multi-fidelity Currin function ###
library(randtoolbox)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
source("GP.R")
source("matern.R")
source("stacking_design.R")

# input setting 
d <- 2             # d: dimension of X (scalar)
n.init <- 5*d
alpha <- NULL
tt <- 2             # mesh = 0.4/(tt^l), lowest level = 0.2
MM <- 16
Lmax <- 6
beta <- 2
cost <- tt^(beta*(1:Lmax))
k <- NULL           # k: the parameter for the Matern kernel function (NULL means it'll be estimated by LOOCV) 
xnew <- expand.grid(seq(0,1,0.02),seq(0,1,0.02)) # xnew: test data (a data frame)
n.max <- 100*d        # the maximum number of sample size
log.fg <- TRUE

# synthetic function
f <- function(X, l, MM, tt){
  fact1 <- 1 - exp(-1/(2*X[,2]))
  fact2 <- 2300*X[,1]^3 + 1900*X[,1]^2 + 2092*X[,1] + 60
  fact3 <- 100*X[,1]^3 + 500*X[,1]^2 + 4*X[,1] + 20
  
  y <- fact1 * fact2/fact3
  y <- y + (MM*tt^(-l))*exp(-1.4*X[,1])*cos(3.5*pi*X[,2])
  return(y)
}

# currin function
f_true <- function(X) {
  fact1 <- 1 - exp(-1/(2*X[,2]))
  fact2 <- 2300*X[,1]^3 + 1900*X[,1]^2 + 2092*X[,1] + 60
  fact3 <- 100*X[,1]^3 + 500*X[,1]^2 + 4*X[,1] + 20
  
  y <- fact1 * fact2/fact3
  return(y)
}

# true values
y.true <- f_true(xnew)

# Figure 3
pdf("synthetic_illustration.pdf", width=10, height=3)
par(mfrow=c(1,4))
par(mar=c(0,1,2,0))
z <- outer(seq(0,1,0.02), seq(0,1,0.02), function(x,y) f(cbind(x,y), 1, MM, tt))
persp(z=z, col="lightblue", theta = 300, phi = 20, ticktype = "detailed", 
      xlab="x1", ylab="x2", zlab="f(x1,x2)", main=expression(italic(l)==1), zlim=c(-5,20))

z <- outer(seq(0,1,0.02), seq(0,1,0.02), function(x,y) f(cbind(x,y), 2, MM, tt))
persp(z=z, col="lightblue", theta = 300, phi = 20, ticktype = "detailed", 
      xlab="x1", ylab="x2", zlab="f(x1,x2)", main=expression(italic(l)==2), zlim=c(-5,20))

z <- outer(seq(0,1,0.02), seq(0,1,0.02), function(x,y) f(cbind(x,y), 3, MM, tt))
persp(z=z, col="lightblue", theta = 300, phi = 20, ticktype = "detailed", 
      xlab="x1", ylab="x2", zlab="f(x1,x2)", main=expression(italic(l)==3), zlim=c(-5,20))

z <- outer(seq(0,1,0.02), seq(0,1,0.02), function(x,y) f_true(cbind(x,y)))
persp(z=z, col="lightblue", theta = 300, phi = 20, ticktype = "detailed", 
      xlab="x1", ylab="x2", zlab="f(x1,x2)", main=expression(italic(l)%->%infinity), zlim=c(-5,20))
dev.off()

# MLGP with stacking design
ML.out <- stacking_design(cost, d, n.init, epsilon=1, alpha, 
                             tt, MM, Lmax, f, k, xnew, n.max, 
                             log.fg=TRUE, norm="L2", model.save=TRUE, save.n = TRUE,
                             init=rep(2.5,d), lower=0.5, upper=5)
# L2 error
print(sqrt(mean((y.true - ML.out$mean)^2)))

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

pdf("synthetic_design.pdf", width=8.5, height=2.5)
layout(matrix(c(rep(1:4,each=2),5),nrow=1))
# par(mfrow=c(1,4))
par(mar=c(2, 2, 4, 0.2))
plot(sobol(sample.size.df[1,4], 2, scrambling = 2),xlab="",ylab="",col=1, yaxt = "n", xaxt = "n",cex=2, pch=0, main=expression(italic(L)==1))
title(xlab="x1", ylab="x2", line=0.5, cex.lab=1.2)
plot(sobol(sample.size.df[5,4], 2, scrambling = 2),xlab="",ylab="",col=1, yaxt = "n", xaxt = "n",cex=2, pch=0, main=expression(italic(L)==2))
points(sobol(sample.size.df[6,4], 2, scrambling = 2),col=2,pch=1,cex=3)
title(xlab="x1", ylab="x2", line=0.5, cex.lab=1.2)
plot(sobol(sample.size.df[9,4], 2, scrambling = 2),xlab="",ylab="",col=1, yaxt = "n", xaxt = "n",cex=2, pch=0, main=expression(italic(L)==3))
points(sobol(sample.size.df[10,4], 2, scrambling = 2),col=2,pch=1,cex=3)
points(sobol(sample.size.df[11,4], 2, scrambling = 2),col=4,pch=2,cex=2)
title(xlab="x1", ylab="x2", line=0.5, cex.lab=1.2)
plot(sobol(sample.size.df[13,4], 2, scrambling = 2),xlab="",ylab="",col=1, yaxt = "n", xaxt = "n",cex=2, pch=0, main=expression(italic(L)==4))
points(sobol(sample.size.df[14,4], 2, scrambling = 2),col=2,pch=1,cex=3)
points(sobol(sample.size.df[15,4], 2, scrambling = 2),col=4,pch=2,cex=2)
points(sobol(sample.size.df[16,4], 2, scrambling = 2),col=5,pch=4,cex=2)
title(xlab="x1", ylab="x2", line=0.5, cex.lab=1.2)

par(mar=c(0, 0, 0, 0))
plot(1,1,xaxt="n",type="n",yaxt="n", bty="n")
legend("center", 
       legend = c(expression(italic(l)==1), expression(italic(l)==2), expression(italic(l)==3), expression(italic(l)==4)), 
       pch=c(0,1,2,4), col=c(1,2,4,5), cex=1.5, pt.cex=c(2,3,2,2), bty = "n", y.intersp=2)
dev.off()

# Figure 4
pdf("synthetic_design_size.pdf", width=10, height=2)
g1 <- ggplot(sample.size.df[1:4,], aes(x = name,y = value ,fill = stage)) +
  geom_bar(stat = "identity") +
  labs(y = "sample size", x = "") + 
  ylim(c(0,116)) + 
  scale_x_discrete(labels=c(expression(n[1]), expression(n[2]), expression(n[3]), expression(n[4]))) +
  scale_fill_manual(breaks=c(1:4), values = c("lightblue1","lightblue2", "lightblue3", "lightblue4")) +
  coord_flip() + theme( panel.background = element_blank(), legend.position = "none")

g1 <- g1 + geom_text(aes(name, total + 3, label = total, fill = NULL), size = 4, data = sample.size.df[sample.size.df$stage==1 & sample.size.df$value > 0,])
g1 <- g1 + theme(axis.text.y = element_text(size = 12))

g2 <- ggplot(sample.size.df[1:8,], aes(x = name,y = value ,fill = stage)) +
  geom_bar(stat = "identity") +
  labs(y = "sample size", x = "") + 
  ylim(c(0,116)) + 
  scale_x_discrete(labels=c(expression(n[1]), expression(n[2]), expression(n[3]), expression(n[4]))) +
  scale_fill_manual(breaks=c(1:4), values = c("lightblue1","lightblue2", "lightblue3", "lightblue4")) +
  coord_flip() + theme( panel.background = element_blank(), legend.position = "none")

g2 <- g2 + geom_text(aes(name, total + 3, label = total, fill = NULL), size = 4,
                     data = sample.size.df[sample.size.df$stage==2 & sample.size.df$value > 0,])
g2 <- g2 + theme(axis.text.y = element_text(size = 12))

g3 <- ggplot(sample.size.df[1:12,], aes(x = name,y = value ,fill = stage)) +
  geom_bar(stat = "identity") +
  labs(y = "sample size", x = "") + 
  ylim(c(0,116)) + 
  scale_x_discrete(labels=c(expression(n[1]), expression(n[2]), expression(n[3]), expression(n[4]))) +
  scale_fill_manual(breaks=c(1:4), values = c("lightblue1","lightblue2", "lightblue3", "lightblue4")) +
  coord_flip() + theme( panel.background = element_blank(), legend.position = "none")

g3 <- g3 + geom_text(aes(name, total + 3, label = total, fill = NULL), size = 4,
                     data = sample.size.df[sample.size.df$stage==3 & sample.size.df$value > 0,])
g3 <- g3 + theme(axis.text.y = element_text(size = 12))

g4 <- ggplot(sample.size.df, aes(x = name,y = value ,fill = stage)) +
  geom_bar(stat = "identity") +
  labs(y = "sample size", x = "") + 
  ylim(c(0,116)) + 
  scale_x_discrete(labels=c(expression(n[1]), expression(n[2]), expression(n[3]), expression(n[4]))) +
  scale_fill_manual(breaks=c(1:4), values = c("lightblue1","lightblue2", "lightblue3", "lightblue4")) +
  coord_flip() + theme( panel.background = element_blank(), legend.position = "none")

g4 <- g4 + geom_text(aes(name, total + c(1,5,5,5), label = total, fill = NULL), size = 4,
                     data = sample.size.df[sample.size.df$stage==4 & sample.size.df$value > 0,])
g4 <- g4 + theme(axis.text.y = element_text(size = 12))

grid.arrange(g1, g2, g3, g4, ncol = 4)
dev.off()


# Figure 5
pdf("alpha_estimation.pdf", width=6, height=2.5)
z.val <- c(ML.out$z[[2]][1:ML.out$n.save[[3]][2]], ML.out$z[[3]][1:ML.out$n.save[[3]][3]])
df_alpha <- data.frame(mesh_l=rep(1/tt^(2:3), ML.out$n.save[[3]][2:3]),
                       z=abs(z.val))
lm.fit <- lm(log(z) ~ log(mesh_l), data=df_alpha)

g1 <- ggplot(df_alpha, aes(x=log(mesh_l), y=log(z), group=log(mesh_l))) + 
  geom_boxplot(fill='#A4A4A4', color="black", width=0.25)+
  xlab(expression(log(italic(T^-l))))+ylab(bquote("log|" ~ f[italic(l)] ~ "(x) -" ~ f[italic(l)-1] ~ "(x)|"))+#ylab(bquote("log|" ~ z[italic(l)] ~  "|"))+
  theme_bw()+ggtitle(expression(italic(L)==3)) + theme(plot.title = element_text(hjust = 0.5))

g1 <- g1 + geom_abline(intercept = lm.fit$coefficients[1], 
                       slope = lm.fit$coefficients[2], color = "red", linetype = "dashed") 

z.val <- unlist(ML.out$z[2:4])
df_alpha <- data.frame(mesh_l=rep(1/tt^(2:4), ML.out$n[2:4]),
                       z=abs(z.val))
lm.fit <- lm(log(z) ~ log(mesh_l), data=df_alpha)

g2 <- ggplot(df_alpha, aes(x=log(mesh_l), y=log(z), group=log(mesh_l))) + 
  geom_boxplot(fill='#A4A4A4', color="black", width=0.25)+
  xlab(expression(log(italic(T^-l))))+ylab(bquote("log|" ~ f[italic(l)] ~ "(x) -" ~ f[italic(l)-1] ~ "(x)|"))+#ylab(bquote("log|" ~ z[italic(l)] ~  "|"))+
  theme_bw() + ggtitle(expression(italic(L)==4)) + theme(plot.title = element_text(hjust = 0.5))

g2 <- g2 + geom_abline(intercept = lm.fit$coefficients[1], 
                       slope = lm.fit$coefficients[2], color = "red", linetype = "dashed") 

grid.arrange(g1, g2, ncol = 2)
dev.off()



# Run stacking designs with various epsilon choices
# L2 norm
epsilon <- c(0.5, 1, 2)

true.error <- rep(0,length(epsilon))
n.all <- matrix(0,nrow=length(epsilon),ncol=Lmax)
cost.all <- rep(0,length(epsilon))

for(i in 1:length(epsilon)){
  ML.out <- stacking_design(cost, d, n.init, epsilon=epsilon[i], alpha, 
                            tt, MM, Lmax, f, k, xnew, n.max, 
                            log.fg=TRUE, norm="L2", model.save=TRUE, save.n = TRUE,
                            init=rep(2.5,d), lower=0.5, upper=5)

  # compute the true error
  true.error[i] <- sqrt(mean((y.true - ML.out$mean)^2))
  # sample size
  n.all[i,1:length(ML.out$n)] <- ML.out$n
  cost.all[i] <- sum(ML.out$n *  ML.out$cost)
}

df1 <- data.frame("epsilon"=rep(epsilon,each=5),
                  "error"=rep(true.error,each=5),
                  "size"=c(n.all[1,1:5], n.all[2,1:5], n.all[3,1:5]),
                  "level"=factor(rep(1:5,length(epsilon))))

# uniform norm
epsilon <- c(1, 2.5, 5)

true.error <- rep(0,length(epsilon))
n.all <- matrix(0,nrow=length(epsilon),ncol=Lmax)
cost.all <- rep(0,length(epsilon))

for(i in 1:length(epsilon)){
  # perform forward ML emulator
  ML.out <- stacking_design(cost, d, n.init, epsilon=epsilon[i], alpha, 
                            tt, MM, Lmax, f, k, xnew, n.max, 
                            log.fg=TRUE, norm="uniform", model.save=TRUE, save.n = TRUE,
                            init=rep(2.5,d), lower=0.5, upper=5)

  # compute the true error
  true.error[i] <- sqrt(mean((y.true - ML.out$mean)^2))
  # sample size
  n.all[i,1:length(ML.out$n)] <- ML.out$n
  cost.all[i] <- sum(ML.out$n *  ML.out$cost)
}

df2 <- data.frame("epsilon"=rep(epsilon,each=5),
                  "error"=rep(true.error,each=5),
                  "size"=c(n.all[1,1:5], n.all[2,1:5], n.all[3,1:5]),
                  "level"=factor(rep(1:5,length(epsilon))))

colnames(df1)[4] <- colnames(df2)[4] <- "meshsize"
df1$meshsize <- factor(MM/tt^as.numeric(df1$meshsize))
df2$meshsize <- factor(MM/tt^as.numeric(df2$meshsize))
df1$meshsize <- factor(df1$meshsize, levels = c(8, 4, 2, 1, 0.5))
df2$meshsize <- factor(df2$meshsize, levels = c(8, 4, 2, 1, 0.5))

# Figure 6
my.cols <- brewer.pal(6, "Blues")
pdf(paste0("illustration epsilon.pdf"), width=11, height=2.8)
g1 <- ggplot(data=df1, aes(x=factor(epsilon))) +
  geom_bar(aes(y=size, fill=meshsize), stat="identity", position=position_dodge())+
  geom_point(aes(y=epsilon*100+6,group=1,shape="epsilon"),size=3)+
  geom_point(aes(y=error*100+6,group=1,shape="test error"),size=3)+
  scale_shape_manual(bquote(L[2] ~ "norm"), values = c(3,5)) +
  scale_y_continuous( 
    # Features of the first axis
    name = "sample size",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~(.-6)/100, name="error"))+ 
  scale_fill_manual(values = my.cols[2:6])+theme_bw() + 
  guides(shape = guide_legend(order = 1),col = guide_legend(order = 2))+
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "(a)")

g2 <- ggplot(data=df2, aes(x=factor(epsilon))) +
  geom_bar(aes(y=size, fill=meshsize), stat="identity", position=position_dodge())+
  geom_point(aes(y=epsilon*40+6,group=1,shape="epsilon"),size=3)+
  geom_point(aes(y=error*40+6,group=1,shape="test error"),size=3)+
  scale_shape_manual(bquote(L[infinity] ~ "norm"), values = c(3,5)) +
  scale_y_continuous( 
    # Features of the first axis
    name = "sample size",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~(.-6)/40, name="error"))+ 
  scale_fill_manual(values = my.cols[2:6])+theme_bw() + 
  guides(shape = guide_legend(order = 1),col = guide_legend(order = 2))+
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "(b)")

grid.arrange(g1, g2, ncol = 2)

dev.off()
