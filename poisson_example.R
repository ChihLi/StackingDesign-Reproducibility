### reproducing Section 5.2: Poisson's equation ###
library(matlabr)
library(randtoolbox)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
source("GP.R")
source("matern.R")
source("stacking_design.R")

# input setting 
tt <- 2             # mesh = 0.4/(tt^l), lowest level = 0.2
MM <- 0.4

# true function
f_true0 <- function(z1,z2,x) {exp(x*z1)*sin(pi*z1)*sin(pi*z2)}
f_true <- function(x) {
  x <- 2*x-1
  integrate(function(z2) { 
    sapply(z2, function(z2) {
      integrate(function(z1) f_true0(z1,z2,x), 0, 1)$value
    })
  }, 0, 1)$value
}

# a function for FEM simulations through matlab
f <- function(X, l, MM, tt, return.time=FALSE){
  if(is.null(dim(X))) n <- length(X) else n <- nrow(X)
  dfr <- data.frame(2*X-1, rep(MM/tt^l, n)) # scale X to [-1,1]
  write.csv(dfr, "Rmatlab_files/generate_text/temp_to_matlab.txt", row.names=F)
  run_matlab_script("Rmatlab_files/SolvePoissonEqn.m", verbose = FALSE, desktop = FALSE, 
                    splash = FALSE, display = FALSE, wait = TRUE, single_thread = FALSE,
                    intern = TRUE)
  dfm <- read.table("Rmatlab_files/generate_text/temp_to_r.txt", sep = ",")
  if(sum(abs(dfm[,1:2]-dfr))>1e-5){print("wrong data!"); break}
  
  if(return.time) return(list(value=dfm[,3],time=dfm[,4])) else return(dfm[,3])
}

# settings 
cost <- c(0.18,0.19,0.23,0.27,0.55) # this should be NULL, but for the sake of reproducibilty, I used the saved costs
d <- 1              # d: dimension of X (scalar)
n.init <- 5*d
alpha <- NULL
tt <- 2             # mesh = 0.4/(tt^l), lowest level = 0.2
MM <- 0.4
Lmax <- 6
k <- NULL           # k: the parameter for the Matern kernel function (NULL means it'll be estimated by LOOCV) 
xnew.ori <- matrix(seq(-1,1,0.001),ncol=1) # xnew.ori: test data (original scale)
xnew <- (xnew.ori+1)/2
n.max <- 100        # the maximum number of sample size
log.fg <- TRUE

# true values
y.true <- apply(xnew,1,f_true)

# perform forward ML emulator
ML.out <- stacking_design(cost, d, n.init, epsilon=0.05, alpha, 
                          tt, MM, Lmax, f, k, xnew, n.max, 
                          log.fg=TRUE, norm="uniform", model.save=TRUE, save.n = TRUE,
                          init=rep(1,d), lower=0.1, upper=5)

# empirical uniform-norm error
print(max(abs((y.true - ML.out$mean))))

# total computational cost
print(sum(ML.out$n * ML.out$cost))

# Figure 7
pdf("poisson_design.pdf", width=11.2, height=2.4)
layout(matrix(c(rep(1:5,each=2),6),nrow=1))
par(mar=c(4, 4, 2, 1))
for(l in 1:length(ML.out$n)){
  title.l <- bquote(italic(L) ~ "=" ~ .(l))
  plot(xnew.ori, y.true,col=1,type="n",lty=2,lwd=1, xlab="x", ylab="y", main=title.l, ylim=c(0.2,0.7))
  polygon(c(rev(xnew.ori), xnew.ori), 
          c(rev(rowSums(ML.out$P[,1:l,drop=FALSE])-rowSums(sqrt(t(t(ML.out$V[,1:l,drop=FALSE])*ML.out$n[1:l])))- ML.out$B[,l]), 
            rowSums(ML.out$P[,1:l,drop=FALSE])+rowSums(sqrt(t(t(ML.out$V[,1:l,drop=FALSE])*ML.out$n[1:l])))+ML.out$B[,l]), col = 'grey80', border = NA)
  
  
  lines(xnew.ori, rowSums(ML.out$P[,1:l,drop=FALSE]),col=2,lwd=3)
  lines(xnew.ori, y.true,col=1,type="l",lty=2,lwd=2)
  for(i in 1:l) points(ML.out$X[[i]][1:ML.out$n.save[[l]][i]]*2-1, ML.out$Y[[i]][1:ML.out$n.save[[l]][i]],pch=c(0,1,2,4,15)[i],col=4,cex=1.5)
}
par(mar=c(0, 0, 0, 0))
plot(1,1,xaxt="n",type="n",yaxt="n", bty="n")
legend("center", legend = c(expression(italic(l)==1), expression(italic(l)==2), expression(italic(l)==3), expression(italic(l)==4), expression(italic(l)==5), expression(f[infinity]), "ML emulator"), 
       pch=c(0,1,2,4,15,NA,NA), col=c(rep(4,5),1,2), lty=c(rep(NA,5),2,1),lwd=c(rep(NA,5),1,2),cex=1.5, bty = "n", 
       pt.cex=c(2,3,2,2,2,NA,NA), y.intersp=2)
dev.off()


sample.size.df <- data.frame(matrix(0,ncol=4,nrow=25))
colnames(sample.size.df) <- c("name", "stage", "value", "total")
sample.size.df[1,4] <- ML.out$n.save[[1]]
sample.size.df[1,3] <- ML.out$n.save[[1]]
sample.size.df[6:7,4] <- ML.out$n.save[[2]]
sample.size.df[6:7,3] <- ML.out$n.save[[2]] - c(ML.out$n.save[[1]],0)
sample.size.df[11:13,4] <- ML.out$n.save[[3]]
sample.size.df[11:13,3] <- ML.out$n.save[[3]] - c(ML.out$n.save[[2]],0)
sample.size.df[16:19,4] <- ML.out$n.save[[4]]
sample.size.df[16:19,3] <- ML.out$n.save[[4]] - c(ML.out$n.save[[3]],0)
sample.size.df[21:25,4] <- ML.out$n.save[[5]]
sample.size.df[21:25,3] <- ML.out$n.save[[5]] - c(ML.out$n.save[[4]],0)

sample.size.df[,1] <- rep(paste0("n",1:5),5)
sample.size.df[,1] <- factor(sample.size.df[,1], levels = paste0("n",1:5))
sample.size.df[,2] <- rep(1:5,each=5)
sample.size.df[,2] <- factor(sample.size.df[,2], levels = c(5,4, 3, 2, 1))

pdf("poisson_design_size.pdf", width=10, height=2)
g1 <- ggplot(sample.size.df[1:5,], aes(x = name,y = value ,fill = stage)) +
  geom_bar(stat = "identity") +
  labs(y = "sample size", x = "") + 
  ylim(c(0,8)) + 
  scale_x_discrete(labels=c(expression(n[1]), expression(n[2]), expression(n[3]), expression(n[4]), expression(n[5]))) +
  scale_fill_manual(breaks=c(1:5), values = c("lightblue", "lightblue1","lightblue2", "lightblue3", "lightblue4")) +
  coord_flip() + theme( panel.background = element_blank(), legend.position = "none")

g1 <- g1 + geom_text(aes(name, total + 1, label = total, fill = NULL), size = 4, data = sample.size.df[sample.size.df$stage==1 & sample.size.df$value > 0,])
g1 <- g1 + theme(axis.text.y = element_text(size = 12))

g2 <- ggplot(sample.size.df[1:10,], aes(x = name,y = value ,fill = stage)) +
  geom_bar(stat = "identity") +
  labs(y = "sample size", x = "") + 
  ylim(c(0,8)) + 
  scale_x_discrete(labels=c(expression(n[1]), expression(n[2]), expression(n[3]), expression(n[4]), expression(n[5]))) +
  scale_fill_manual(breaks=c(1:5), values = c("lightblue", "lightblue1","lightblue2", "lightblue3", "lightblue4")) +
  coord_flip() + theme( panel.background = element_blank(), legend.position = "none")

g2 <- g2 + geom_text(aes(name, total + 1, label = total, fill = NULL), size = 4,
                     data = sample.size.df[sample.size.df$stage==2,][1:2,])
g2 <- g2 + theme(axis.text.y = element_text(size = 12))

g3 <- ggplot(sample.size.df[1:15,], aes(x = name,y = value ,fill = stage)) +
  geom_bar(stat = "identity") +
  labs(y = "sample size", x = "") + 
  ylim(c(0,8)) + 
  scale_x_discrete(labels=c(expression(n[1]), expression(n[2]), expression(n[3]), expression(n[4]), expression(n[5]))) +
  scale_fill_manual(breaks=c(1:5), values = c("lightblue", "lightblue1","lightblue2", "lightblue3", "lightblue4")) +
  coord_flip() + theme( panel.background = element_blank(), legend.position = "none")

g3 <- g3 + geom_text(aes(name, total + 1, label = total, fill = NULL), size = 4,
                     data = sample.size.df[sample.size.df$stage==3,][1:3,])
g3 <- g3 + theme(axis.text.y = element_text(size = 12))

g4 <- ggplot(sample.size.df[1:20,], aes(x = name,y = value ,fill = stage)) +
  geom_bar(stat = "identity") +
  labs(y = "sample size", x = "") + 
  ylim(c(0,8)) + 
  scale_x_discrete(labels=c(expression(n[1]), expression(n[2]), expression(n[3]), expression(n[4]), expression(n[5]))) +
  scale_fill_manual(breaks=c(1:5), values = c("lightblue", "lightblue1","lightblue2", "lightblue3", "lightblue4")) +
  coord_flip() + theme( panel.background = element_blank(), legend.position = "none")

g4 <- g4 + geom_text(aes(name, total + 1, label = total, fill = NULL), size = 4,
                     data = sample.size.df[sample.size.df$stage==4,][1:4,])
g4 <- g4 + theme(axis.text.y = element_text(size = 12))

g5 <- ggplot(sample.size.df, aes(x = name,y = value ,fill = stage)) +
  geom_bar(stat = "identity") +
  labs(y = "sample size", x = "") + 
  ylim(c(0,8)) + 
  scale_x_discrete(labels=c(expression(n[1]), expression(n[2]), expression(n[3]), expression(n[4]), expression(n[5]))) +
  scale_fill_manual(breaks=c(1:5), values = c("lightblue", "lightblue1","lightblue2", "lightblue3", "lightblue4")) +
  coord_flip() + theme( panel.background = element_blank(), legend.position = "none")

g5 <- g5 + geom_text(aes(name, total + 1, label = total, fill = NULL), size = 4,
                     data = sample.size.df[sample.size.df$stage==5,])
g5 <- g5 + theme(axis.text.y = element_text(size = 12))

grid.arrange(g1, g2, g3, g4, g5, ncol = 5)
dev.off()


# Run stacking designs with various epsilon choices
# L2 norm 
epsilon <- c(0.075, 0.05, 0.015) 

true.error <- rep(0,length(epsilon))
n.all <- matrix(0,nrow=length(epsilon),ncol=Lmax)

for(i in 1:length(epsilon)){
  ML.out <- stacking_design(cost, d, n.init, epsilon=epsilon[i], alpha, 
                            tt, MM, Lmax, f, k, xnew, n.max, 
                            log.fg=TRUE, norm="L2", model.save=TRUE, save.n = TRUE,
                            init=rep(1,d), lower=0.1, upper=5)
  
  # compute the true error
  true.error[i] <- sqrt(mean((y.true - ML.out$mean)^2))
  # sample size
  n.all[i,1:length(ML.out$n)] <- ML.out$n
}

df1 <- data.frame("epsilon"=rep(epsilon,each=5),
                  "error"=rep(true.error,each=5),
                  "size"=c(n.all[1,1:5], n.all[2,1:5], n.all[3,1:5]),
                  "level"=factor(rep(1:5,length(epsilon))))


# uniform norm
epsilon <- c(0.1, 0.075, 0.05) 

true.error <- rep(0,length(epsilon))
n.all <- matrix(0,nrow=length(epsilon),ncol=Lmax)

for(i in 1:length(epsilon)){
  ML.out <- stacking_design(cost, d, n.init, epsilon=epsilon[i], alpha, 
                            tt, MM, Lmax, f, k, xnew, n.max, 
                            log.fg=TRUE, norm="uniform", model.save=TRUE, save.n = TRUE,
                            init=rep(1,d), lower=0.1, upper=5)
  # compute the true error
  true.error[i] <- max(abs((y.true - ML.out$mean)))
  # sample size
  n.all[i,1:length(ML.out$n)] <- ML.out$n
}

df2 <- data.frame("epsilon"=rep(epsilon,each=5),
                  "error"=rep(true.error,each=5),
                  "size"=c(n.all[1,1:5], n.all[2,1:5], n.all[3,1:5]),
                  "level"=factor(rep(1:5,length(epsilon))))


df1$level <- factor(MM/tt^as.numeric(df1$level))
df2$level <- factor(MM/tt^as.numeric(df2$level))
colnames(df1)[4] <- colnames(df2)[4] <- "meshsize"
df1$meshsize <- factor(df1$meshsize, levels = c(0.2, 0.1, 0.05, 0.025, 0.0125))
df2$meshsize <- factor(df2$meshsize, levels = c(0.2, 0.1, 0.05, 0.025, 0.0125))

# Figure 8
my.cols <- brewer.pal(6, "Blues")
pdf(paste0("poisson epsilon.pdf"), width=11, height=2.8)
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
  geom_point(aes(y=epsilon*80+6,group=1,shape="epsilon"),size=3)+
  geom_point(aes(y=error*80+6,group=1,shape="test error"),size=3)+
  scale_shape_manual(bquote(L[infinity] ~ "norm"), values = c(3,5)) +
  scale_y_continuous( 
    # Features of the first axis
    name = "sample size",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~(.-6)/80, name="error"))+ 
  scale_fill_manual(values = my.cols[2:6])+theme_bw() + 
  guides(shape = guide_legend(order = 1),col = guide_legend(order = 2))+
  scale_x_discrete(labels = NULL, breaks = NULL) + labs(x = "(b)")

grid.arrange(g1, g2, ncol = 2)

dev.off()
