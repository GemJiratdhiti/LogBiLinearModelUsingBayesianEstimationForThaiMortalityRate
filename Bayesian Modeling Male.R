library(MASS)
library(forecast)
library(TruncatedNormal)
library(ggplot2)

# import and transform data
setwd("...") # Set your own working directory
ABmle <- read.csv("Male data/output_male.csv")
Ktmle <- read.csv("Male data/output_kt_male.csv")
Dxt <- as.matrix(read.csv("Male data/Mortality Male.csv"))[,-1:-2]
Ext <- as.matrix(read.csv("Male data/Population Male.csv"))[,-1:-2]
Amle <- as.numeric(ABmle[,2])
Bmle <- as.numeric(ABmle[,3])
Kmle <- data.frame(cbind(1:20, t(Ktmle)))

# set loops
set.seed(2500)
N = 20000
M = 101
t = 20

### initial values ###

# fit kt in SLR to find g0 and S0
SLR <- lm(Kmle[,2]~Kmle[,1], data = Kmle)
g0 = t(t(SLR$coefficients))
S0 = vcov(SLR)

# fit transformed kt in AR(1) to find roh and s2k
X <- matrix(1, nrow = t, ncol = 2)
X[, 2] = 1:t
Eta <- X %*% g0
plot.ts(Kmle[2] - Eta)
AR <- arima(Kmle[2] - Eta, order = c(1,0,0))
roh.ar = AR$coef[1]
s2k.ar = AR$sigma2

# other initial values
s2b.mle = var(Bmle)
alpha.mle = Amle
beta.mle = Bmle
kappa.mle = Kmle[,2]
s2p = 1
Ab = 2.1
Bb = (Ab - 1) * s2b.mle
Ak = 2.1
Bk = (Ak - 1) * s2k.ar
Ba = 0.001
Aa = Ba * exp(alpha.mle)

P <- matrix(0, nrow = t, ncol = t)
for(i in 1:t - 1){
  P[i + 1, i] = roh.ar
}

Q <- matrix(0, nrow = t, ncol = t)
for(i in 1:t - 1){
  Q[i, i + 1] = -roh.ar
  Q[i, i] = 1 + roh.ar^2
  Q[i + 1, i] = -roh.ar
}
Q[t, t] = 1

gamma <- matrix(0, nrow = 2, ncol = N)
roh <- rep(0,N)
s2k <- rep(0,N)
s2b <- rep(0,N)
alpha <- matrix(0, nrow = M, ncol = N)
beta <- matrix(0, nrow = M, ncol = N)
kappa <- matrix(0, nrow = t, ncol = N)

gamma[,1] <- g0
roh[1] <- roh.ar
s2k[1] <- s2k.ar
s2b[1] <- s2b.mle
alpha[,1] <- alpha.mle
beta[,1] <- beta.mle
kappa[,1] <- kappa.mle

Pgamma <- c()
Proh <- c()
Ps2k <- c()
Ps2b <- c()
Palpha <- c()
Pbeta <- c()
Pkappa <- c()

### MCMC ###

countk = 0
countb = 0
sigmax = sqrt(s2b.mle / 32)
sigmat = 1/2
for(i in 2:N){
  # gamma
  S. <- solve(t(X) %*% Q %*% X + s2k[i - 1] * solve(S0))
  g. <- S. %*% (t(X) %*% Q %*% kappa[,i - 1] + s2k[i - 1] * solve(S0) %*% g0)
  gamma[,i] <- t(mvrnorm(1, g., s2k[i - 1] * S.))
  
  # s2b
  s2b[i] <- 1 / rgamma(1, Ab + M / 2, Bb + (t(beta[,i - 1]) %*% beta[,i - 1] / 2))
  
  # s2k
  eta <- X %*% gamma[,i]
  S1 = (kappa[1,i - 1] - eta[1]) ^ 2
  for(j in 2:t){
    s1 <- (kappa[j,i - 1] - eta[j] - roh[i - 1] * (kappa[j - 1,i - 1] - eta[j - 1])) ^ 2
    S1 = S1 + s1
  }
  s2k[i] <- 1 / rgamma(1, Ak + t / 2, Bk + S1 / 2)
  
  # roh
  Ap = 0
  Bp = 0
  for(j in 2:t){
    ap <- (kappa[j - 1,i - 1] - eta[j - 1]) ^ 2
    Ap = Ap + ap
    bp <- (kappa[j,i - 1] - eta[j]) * (kappa[j - 1,i - 1] - eta[j - 1])
    Bp = Bp + bp
  }
  mup. <- Bp / (Ap + s2k[i] / s2p)
  s2p. <- s2k[i - 1] / (Ap + s2k[i] / s2p)
  roh[i] <- rtnorm(1, mup., s2p., -1, 1)
  
  # kappa
  # x = kappa[j,i - 1]
  # y = 1
  # k = kappa[,i - 1]
  # a = alpha[,i - 1]
  # b = beta[,i - 1]
  llh.kt <- function(x, y, k, a, b){
    f.Dt = 0
    fkt <- -(x - eta[1]) ^ 2 / (2 * s2k[i])
    for(m in 1:M){
      fDt <- as.numeric(-Ext[m,y] * exp(a[m] + b[m] * x) + b[m] * x * Dxt[m,y])
      f.Dt = f.Dt + fDt
    }
    if(y < t){
      fkt1.kt <- -(k[y + 1] - eta[y + 1] - roh[i] * (x - eta[y])) ^ 2 / (2 * s2k[i])
    }
    if(y > 1){
      fkt.kt1 <- -(x - eta[y] - roh[i] * (k[y - 1] - eta[y - 1])) ^ 2 / (2 * s2k[i])
    }
    if(y == 1){
      llh <- f.Dt + fkt + fkt1.kt
    }
    else if(y == t){
      llh <- f.Dt + fkt.kt1
    }
    else{
      llh <- f.Dt + fkt.kt1 + fkt1.kt
    }
    return(llh)
  }
  kap <- kappa[,i - 1]
  alp <- alpha[,i - 1]
  bet <- beta[,i - 1]
  j = 1
  for(j in 1:t){
    pp.kt <- rnorm(1, kappa[j,i - 1], sigmat)
    psi1 <- min(1, exp(llh.kt(pp.kt, j, kap, alp, bet) - llh.kt(kappa[j,i - 1], j, kap, alp, bet)))
    u <- runif(1)
    if(u > psi1 | is.nan(psi1) == T){
      kap[j] <- kappa[j,i - 1]
    }
    else{
      kap[j] <- pp.kt
      countk = countk + 1
    }
    kbar <- mean(kap)
    kap <- kap - kbar
    alp <- alp + bet * kbar
  }
  kappa[,i] <- kap 
  alpha[,i] <- alp
  
  # beta
  # x = pp.bx
  # y = 68
  # k = kap
  # a = alp
  # b = bet
  llh.bx <- function(x, y, k, a, b){
    f.Dt = 0
    fbx <- -x ^ 2 / (2 * s2b[i])
    for(n in 1:t){
      fDt <- as.numeric(-Ext[y,n] * exp(a[y] + x * k[n]) + x * k[n] * Dxt[y,n])
      f.Dt = f.Dt + fDt
    }
    llh <- f.Dt + fbx
    return(llh)
  }
  kap <- kappa[,i]
  alp <- alpha[,i]
  bet <- beta[,i - 1]
  for(j in 1:M){
    pp.bx <- rnorm(1, beta[j,i - 1], sigmax)
    psi2 <- min(1, exp(llh.bx(pp.bx, j, kap, alp, bet) - llh.bx(beta[j,i - 1], j, kap, alp, bet)))
    v <- runif(1)
    if(v > psi2 | is.nan(psi2) == T){
      bet[j] <- beta[j,i - 1]
    }
    else{
      bet[j] <- pp.bx
      countb = countb + 1
    }
    bsum <- sum(bet)
    bet <- bet / bsum
    kap <- kap * bsum
  }
  kappa[,i] <- kap
  beta[,i] <- bet
  
  # alpha
  c <- rep(0,M)
  D. <- rep(0,M)
  for(j in 1:M){
    c[j] <- Ext[j,] %*% exp(beta[j,i] * kappa[,i])
    D.[j] <- rowSums(Dxt)[j]
    alpha[j,i] <- log(rgamma(1, Aa[j] + D.[j], Ba + c[j]))
  }
  
  # posterior distribution (burn-in)
  if(i > N / 2 & i %% 10 == 1){
    Pgamma <- cbind(Pgamma, gamma[,i])
    Proh <- cbind(Proh, roh[i])
    Ps2k <- cbind(Ps2k, s2k[i])
    Ps2b <- cbind(Ps2b, s2b[i])
    Palpha <- cbind(Palpha, alpha[,i])
    Pbeta <- cbind(Pbeta, beta[,i])
    Pkappa <- cbind(Pkappa, kappa[,i])
  }
}

### Checking ###

rateK <- countk / (t * N)
rateB <- countb / (M * N)
rateK
rateB
hist(Pgamma)
plot(Pgamma[1,], type = "b")
plot(Pgamma[2,], type = "b")
hist(Proh)
plot(Proh[1,], type = "b")
hist(Ps2k)
plot(Ps2k[1,], type = "b")
hist(Ps2b)
plot(Ps2b[1,], type = "b")

# Convergence of Parameters
par(mfrow = c(2,2))
for(i in 1:4){
  a <- 30 * (i - 1)
  ma <- paste("Alpha ", a)
  plot(1:200, alpha[a + 1,1:200], main = ma, type = "l")
  abline(h = alpha.mle[i], col = "blue")
}
for(i in 1:4){
  a <- 30 * (i - 1)
  mb <- paste("Beta ", a)
  plot(1:200, beta[a + 1,1:200], main = mb, type = "l")
  abline(h = beta.mle[i], col = "blue")
}
for(i in 1:t){
  mk <- paste("Kappa ", i + 1997)
  plot(1:200, kappa[i,1:200], main = mk, type = "l")
  abline(h = kappa.mle[i], col = "blue")
}
plot(1:100, kappa[1,1:100], type = "b")

# Posterior Distribution
par(mfrow = c(1,1))
par(mfrow = c(2,2))
for(i in 1:4){
  a <- 30 * (i - 1)
  m1 <- paste("Alpha ", a)
  hist(Palpha[a + 1,], main = m1)
  abline(v = mean(Palpha[a + 1,]), col = "red")
  abline(v = alpha.mle[a + 1], col = "blue")
}
for(i in 1:4){
  a <- 30 * (i - 1)
  m1 <- paste("Beta ", a)
  hist(Pbeta[a + 1,], main = m1)
  abline(v = mean(Pbeta[a + 1,]), col = "red")
  abline(v = beta.mle[a + 1], col = "blue")
}
for(i in 1:t){
  m2 <- paste("Kappa ", (i + 1997))
  hist(Pkappa[i,], main = m2)
  abline(v = mean(Pkappa[i,]), col = "red")
  abline(v = kappa.mle[i], col = "blue")
}
for(i in 1:4){
  a <- 30 * (i - 1)
  m1 <- paste("Alpha ", a)
  plot(Palpha[a + 1,], main = m1, type = "b")
  abline(h = mean(Palpha[a + 1,]), col = "red")
  abline(h = alpha.mle[a + 1], col = "blue")
}
for(i in 1:4){
  a <- 30 * (i - 1)
  m1 <- paste("Beta ", a)
  plot(Pbeta[a + 1,], main = m1, type = "b")
  abline(h = mean(Pbeta[a + 1,]), col = "red")
  abline(h = beta.mle[a + 1], col = "blue")
}
for(i in 1:t){
  m2 <- paste("Kappa ", (i + 1997))
  plot(Pkappa[i,], main = m2, type = "b")
  abline(h = mean(Pkappa[i,]), col = "red")
  abline(h = kappa.mle[i], col = "blue")
}

### Bayesian Estimation ###

PostmeanA <- c()
PostmeanB <- c()
PostmeanK <- c()
for(i in 1:M){
  PostmeanA <- cbind(PostmeanA, mean(Palpha[i,]))
  PostmeanB <- cbind(PostmeanB, mean(Pbeta[i,]))
  if(i <= t){
    PostmeanK <- cbind(PostmeanK, mean(Pkappa[i,]))
  }
}

par(mfrow = c(1,1))
plot(x = 0:(M - 1), y = PostmeanA, main = "Estimated Alpha (M)",
     xlab = "Age", ylab = "Alpha", type = "b", col = "red", pch = 4, lty = 4)
lines(x = 0:(M - 1), y = alpha.mle, type = "b", col = "blue", pch = 3, lty = 3)
legend(x = 0, y = -2, legend = c("MLE","Bayes"), col = c("blue","red"),
       lty = c(3,4), pch = c(3,4), cex = 0.7)
plot(x = 0:(M - 1), y = PostmeanB, main = "Estimated Beta (M)",
     xlab = "Age", ylab = "Beta", type = "b", col = "red", pch = 4)
lines(x = 0:(M - 1), y = beta.mle, type = "b", col = "blue", pch = 3)
legend(x = 72, y = 0.075, legend = c("MLE","Bayes"), col = c("blue","red"),
       lty = c(3,4), pch = c(3,4), cex = 0.7)
plot(x = 2540:2559, y = PostmeanK, main = "Estimated Kappa (M)",
     xlab = "Year", ylab = "Kappa", type = "b", col = "red", pch = 4)
lines(x = 2540:2559, y =  kappa.mle, type = "b", col = "blue", pch = 3)
legend(x = 2553.5, y = 11, legend = c("MLE","Bayes"), col = c("blue","red"),
       lty = c(3,4), pch = c(3,4), cex = 0.7)
PostmeanA - alpha.mle
PostmeanB - beta.mle
PostmeanK - kappa.mle

### Log mortailty rate ###

logm.mle <- function(i){
  return(alpha.mle + beta.mle * kappa.mle[i])
}
plot.logm.mle <- function(a, b){
  plot(x = 0:(M - 1), y = logm.mle(b - 2539), main = "Male log mortal rate (MLE)",
       xlab = "Age", ylab = "Log mortality rate", type = "l", col = "red")
  for(i in (b - a):1){
    lines(x = 0:(M - 1), y = logm.mle(b - 2539 - i), type = "l",
          col = rainbow(1, start = (i + 0.226 * (b - a - i) - 4.3) / (b - a)))
  }
}
plot.logm.mle(2540,2559)

logm.bayes <- function(i){
  return(PostmeanA + PostmeanB * PostmeanK[i])
}
plot.logm.bayes <- function(a, b){
  plot(x = 0:(M - 1), y = logm.bayes(b - 2539), main = "Male log mortal rate (Bayes)",
       xlab = "Age", ylab = "Log mortality rate", type = "l", col = "red")
  for(i in (b - a):1){
    lines(x = 0:(M - 1), y = logm.bayes(b - 2539 - i), type = "l",
          col = rainbow(1, start = (i + 0.226 * (b - a - i) - 4.3) / (b - a)))
  }
}
plot.logm.bayes(2540,2559)

### Print value ###

write.csv(PostmeanA, "output_axbay_male.csv", row.names = F)
write.csv(PostmeanB, "output_bxbay_male.csv", row.names = F)
write.csv(PostmeanK, "output_ktbay_male.csv", row.names = F)