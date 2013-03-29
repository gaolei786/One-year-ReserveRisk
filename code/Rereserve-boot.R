library(ChainLadder)

#############################׼����һ���ڷ��յļ��㡢�봫ͳ���յıȽ�##########################################
#���ȣ������ݷ���׼����
data <- read.table("data3.txt", sep = "")
rownames(data) <- 0:9
colnames(data) <- 0:9
data <- cbind(expand.grid(0:9, 0:9), as.vector(as.matrix(data)))
data[, 3] <- ifelse(data[, 3] == 0, NA, data[, 3])
data <- data[!is.na(data[[3]]), ]
names(data) <- c("accyear", "devyear", "value")
data <- as.triangle(data, origin = "accyear", dev = "devyear", value = "value")
data <- incr2cum(data)
Mack.data <- MackChainLadder(data, est.sigma="Mack")
f <- Mack.data$f[1:9]
R0 <- summary(Mack.data)$Totals[4, 1]
sum(data[1:9, 2]) / sum(data[1:9, 1])


Mack.summary <- summary(Mack.data)
write.table(Mack.summary$ByOrigin[, c(1, 3, 4)], "clipboard", sep = "\t")

#�ù�һ����ģ�ͣ���ɢ�����ɷֲ�ģ�ͣ�����׼����
mle <- function() {

  y <- NULL
  X0 <- array(0, dim = c(55, 20))
  for (i in 1:10) {
    for (j in 1:(11 - i)) {
      y <- c(y, X[i, j])
      X0[length(y), c(i, 10 + j)] <- c(1, 1)
    }
  }
  coef <- exp(coefficients(glm(y ~ X0 - 1, family = poisson(link = "log"))))
  u <- coef[1:10] * sum(X[1, ]) / coef[1]
  names(u) <- NULL
  r <- c((coef[11:20] * coef[1] / sum(X[1, ]))[1:9], 1 - sum(coef[11:19] * coef[1] / sum(X[1, ])))
  names(r) <- NULL
  fit <- fitted(glm(y ~ X0 - 1, family = poisson(link = "log")))
  phi <- sum(((y - fit) / fit^0.5)^2)/(55-19)
  return(list(u = u, r = r, phi = phi))
}
X <- cum2incr(data)
result <- mle()
u <- mle()$u
r <- mle()$r
phi <- mle()$phi

###############################����bootstrap���δ���⸶��ģ��ֲ�
dev <- f[9:1]
dev <- cumprod(dev)
Mack.boot.result1 <- list() #�����µ�ģ�������⸶����������
Mack.boot.result2 <- list() #����ÿ�εõ���׼����
Mack.boot.result3 <- NULL  #����õ���Ԥ�������⸶֮�ͣ���Ԥ���⸶������׼����ͬ��
set.seed(23)               #׼����������������Ԥ���⸶��������Ϊ���������㳬ɢ���ֲ����������
for (k in 1:5000) {
  data.fitted <- data
  for (i in 1:9) {
      dev <- f[(10 - i):1]
      dev <- cumprod(dev)
      data.fitted[i, 1:(10-i)] <- (data[i, 11-i]/dev[1:(10-i)])[(10-i):1]
  }
  data.fitted.incr <- cum2incr(data.fitted)
  data.incr <- cum2incr(data)
  data.unscaled.resi1 <- (data.incr - data.fitted.incr) / data.fitted.incr^0.5
  data.unscaled.resi <- as.vector(data.unscaled.resi1)
  data.unscaled.resi <- data.unscaled.resi[!is.na(data.unscaled.resi)]
 
  data.unscaled.resi2 <- data.unscaled.resi1
  for (i in 1:10) {
    for (j in 1:(11 - i)) {
      data.unscaled.resi2[i, j] <- sample(data.unscaled.resi, 1)
    }
  }
  data.boot <- data.unscaled.resi2 * (data.fitted.incr)^0.5 +data.fitted.incr
  for (i in 1:10) {
    for (j in 1:(11-i)) {
      if(data.boot[i, j] < 0) {
      data.boot[i, j] <- data.fitted.incr[i, j]
      }
    }
  }
  Mack.boot <- MackChainLadder(incr2cum(data.boot), est.sigma="Mack")
  
  Mack.boot.result1[[k]] <- cum2incr(data)
  RR <- 0
  for (i in 2:10) {
    for (j in (12-i):10) {
      print(c(i,j,cum2incr(Mack.boot$FullTriangle)[i, j]))
       sim.c <- rgamma(1, shape=cum2incr(Mack.boot$FullTriangle)[i, j]/phi,scale=phi) #��٤���ֲ����������ɢ�����ɷֲ�

       Mack.boot.result1[[k]][i, j] <- sim.c
       RR <- RR + sim.c
    }
  }

  Mack.boot.result2[[k]] <- summary(Mack.boot)$ByOrigin$IBNR
  Mack.boot.result3 <- c(Mack.boot.result3, RR)
}

#Mack.boot.result3[1] <- sapply(Mack.boot.result2, sum)[1]


################����Re-reserving������û���һ���ڷ��յ�δ���⸶ģ��ֲ�#############################

T <- NULL  #�����ظ�����׼����
C.result <- NULL  #�¶Խ����ϵ���������Щ����ģ��Ϊbootstrap�Ľ��
All.C <- list()  #��������������
set.seed(123)
for (k in 1:5000) { 
  C <- NULL
  for (i in 1:9) {
    
   # u <- data[i + 1, 10 - i] * f[10-i]
    #sigma2 <-  (Mack.data$sigma)^2 * data[i + 1, 10 - i]
    #data[i + 1, 11 - i] <- rnorm(1, u, (sigma2)^0.5)
    #data[i+1, 11 - i] <- data[i + 1, 10 - i] + phi * rpois(1, u[i+1]*r[11-i]/phi)
    data[i+1, 11 - i] <- data[i + 1, 10 - i] + Mack.boot.result1[[k]][i+1, 11-i]  
    C <- c(C, data[i + 1, 11 - i] - data[i + 1, 10 - i])
  }
  C.result <- rbind(C.result, C)
  Mack.1 <- MackChainLadder(data, est.sigma="Mack")
  All.C[[k]] <-  Mack.1$FullTriangle
  R1 <- summary(Mack.1)$Totals[4, 1]
  T <- c(T, sum(C) + R1)
}
windows()

###########��7�¹���ģ��############
plot(0:3, cum2incr(All.C[[1]])[8, 1:4], type = "b", xlim = c(0, 11), ylim = c(0, 3.2 * 10^6),
  xlab = "��չ��",  ylab = "�����⸶")
X83 <- NULL
fun1 <- function(x) {
 lines(0:3, cum2incr(x)[8, 1:4], col = rgb(0, 0,0, 0.05))
  points(rnorm(1, 4.5, 0.3), rnorm(1, cum2incr(x)[8, 4], 1000), cex = 0.2, col = rgb(0, 0,0, 0.05)) 
  X83 <- c(X83, cum2incr(x)[8, 4])
}
X83 <- sapply(All.C, fun1)
x <- X83
lines(4.5 + 4096528 * density(x)$y, density(x)$x, col = rgb(0, 0, 0, 0.2))
#lines(4.5 + 2996528 *dpois( round(seq(500000, 2500000, length = 400)/phi, 0), u[8]*r[4]/(phi))/phi, 
 # lty = 2, round(seq(500000, 2500000, length = 400), 0),col = rgb(0, 0, 0, 0.2))
abline(v =4.5, col = rgb(0, 0, 0, 0.2), lwd = 2)

####################################################################################

#########################׼�����չ�������#############################################

#########ɢ��ͼ�����ܶȹ���##################
Devres <- summary(Mack.data)$Totals[4, 1] - T
windows()
plot(Devres, xlab = "", main = "Scatterplot of CDR", col = rgb(0, 0, 0, 0.2), ylab = "", xlim = c(0, 8000))
points(c(1:5000)[Devres < quantile(Devres, probs = c(0.015))],Devres[Devres < quantile(Devres, probs = c(0.015))], 
pch = 7, col = 2, cex = 1)
abline(h = quantile(Devres, probs = c(0.015)), col = rgb(0, 0, 0, 0.2), lwd = 3)
abline(h = -quantile(Devres, probs = c(0.015)), col = rgb(0, 0, 0, 0.2), lwd = 3)
abline(h = 0, col = rgb(0, 0, 0, 0.5), lwd = 3)
lines(5100 +400*26884720 * density(Devres)$y, density(Devres)$x, col = rgb(0, 0, 0, 0.5))
abline(v =5100, col = rgb(0, 0, 0, 0.2), lwd = 2)
polygon(c(5100 +400*26884720 *density(Devres)$y[density(Devres)$x<quantile(Devres, probs = c(0.015))], 5100),
 c(density(Devres)$x[density(Devres)$x<quantile(Devres, probs = c(0.015))],quantile(Devres, probs = c(0.015))),
  col = "orange")


########��Ӧ�����⸶�仯###################

par(mfrow=c(2, 2))
for (i in c(6, 7, 8, 9)) {
#hist(C.result[, i])
u1 <- u[i+1]*r[11-i]
plot(C.result[, i],xlab = "", ylab = "", main = paste("I(", i, ")"), col = rgb(0, 0, 0, 0.05))
points(c(1:5000)[Devres < quantile(Devres, probs = c(0.015))],C.result[Devres < quantile(Devres, probs = c(0.015)),i],
 pch = 7, col = 2, cex = 1)
abline(h = u1)
#abline(v = c(1:500)[Devres < quantile(Devres, probs = c(0.025))])
}

#######################################################################################

#####################Reserving ��bootstrap������###########################################

############����ֲ�����##########
plot(0, 0, xlim = c(10000000, 30000000), ylim = c(0, 2.9*10^(-7)), type = "n", xlab = "", ylab = "" )
lines(density(T, na.rm = T),  lwd = 3, col = rgb(0, 0, 0, 0.4))
lines(density(Mack.boot.result3, na.rm = T), lwd = 3, col = rgb(0, 0, 0, 0.1))
polygon( c(density(T)$x[density(T)$x<quantile(T, probs = c(0.05))],
quantile(T, probs = c(0.05))),c(density(T)$y[density(T)$x<quantile(T, probs = c(0.05))], 0),
, col = rgb(0, 0, 0, 0.4), border = NA)
polygon( c(density(T)$x[density(T)$x>quantile(T, probs = c(0.95))],
quantile(T, probs = c(0.95))),c(density(T)$y[density(T)$x>quantile(T, probs = c(0.95))], 0),
, col = rgb(0, 0, 0, 0.4), border = NA)
polygon( c(density(Mack.boot.result3)$x[density(Mack.boot.result3)$x<quantile(Mack.boot.result3, probs = c(0.05))],
quantile(Mack.boot.result3, probs = c(0.05))),c(density(Mack.boot.result3)$y[density(Mack.boot.result3)$x<quantile(Mack.boot.result3, probs = c(0.05))], 0),
, col = rgb(0, 0, 0, 0.1), border = NA)
polygon( c(density(Mack.boot.result3)$x[density(Mack.boot.result3)$x>quantile(Mack.boot.result3, probs = c(0.95))],
quantile(Mack.boot.result3, probs = c(0.95))),c(density(Mack.boot.result3)$y[density(Mack.boot.result3)$x>quantile(Mack.boot.result3, probs = c(0.95))], 0),
, col = rgb(0, 0, 0, 0.1), border = NA)
legend("topright", c("R", "RR"), lty = 1, col = c(rgb(0, 0, 0, 0.1),  rgb(0, 0, 0, 0.4)), lwd = 3)
arrows(2.3e07, 0.1e-7, 2.5e07, 1e-7,col = rgb(0, 0, 0, 0.4), length = 0.15)
text(2.5e07, 1.2e-7, expression(alpha==0.05))
###############################

##########�ܵı䶯################
windows()
par(mfrow = c(1, 2))
set.seed(123)
plot(0:9, All.C[[1]][8, ], type = "n", ylim = c(3*10^6, 16 * 10^6), xlab = "��չ��",ylab=expression(C[paste(7,",",j)]))
abline(v = 3, lwd = 2, col = rgb(0, 0, 0, 0.1))
fun1 <- function(x) {
 lines(0:9, x[8, ], col = rgb(0, 0, 0, 0.01), lwd = 0.01)
}
sapply(All.C, fun1)
lines(2:9, (All.C[[100]])[8, 3:10], col = "white", lwd = 2)
 
set.seed(123)
plot(0:9, c(data[8, 1:3],incr2cum(Mack.boot.result1[[1]])[8, -(1:3)]), xlab = "��չ��", type = "n",
 ylim = c(3*10^6, 16* 10^6),ylab=expression(C[paste(7,",",j)]))
fun1 <- function(x) {
 lines(0:9, c(data[8, 1:3],incr2cum(x)[8, -(1:3)]),col = rgb(0, 0, 0, 0.01), lwd = 0.01)
}
abline(v = 3, lwd = 2, col = rgb(0, 0, 0, 0.1))
sapply(Mack.boot.result1, fun1)
lines(2:9, c((data)[8, 3],incr2cum(Mack.boot.result1[[100]])[8, -(1:3)]), col = "white", lwd = 2)
##############################

##############�����䶯############
par(mfrow= c(1, 2))
plot(0:9, cum2incr(All.C[[1]])[8, ], type = "n", ylim = c(0, 2.5 * 10^6), xlab = "��չ��", ylab =expression(X[paste(7,",",j)]))
fun1 <- function(x) {
  lines(0:9, cum2incr(x)[8, ], col = rgb(0, 0, 0, 0.01), lwd = 0.01)
}
sapply(All.C, fun1)
lines(2:9, cum2incr(All.C[[100]])[8, 3:10], col = "white", lwd = 2)
plot(0:9, c(cum2incr(data)[8, 1:3],(Mack.boot.result1[[1]])[8, -(1:3)]), type = "n", ylim = c(0, 2.5* 10^6), xlab = "��չ��", ylab =expression(X[paste(7,",",j)]))
fun1 <- function(x) {
 lines(0:9, c(cum2incr(data)[8, 1:3],x[8, -(1:3)]), col = rgb(0, 0, 0, 0.01), lwd = 0.01)
}
sapply(Mack.boot.result1, fun1)
lines(2:9, c(cum2incr(data)[8, 3],(Mack.boot.result1[[100]])[8, -(1:3)]), col = "white", lwd = 2)
###############################

########################################################################################



