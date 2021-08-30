# This script answer the following question:
# QUESTION: Is gglasso fitted w/o intercept on centered Y
#           the same as gglasso on Y with intercept?
#---------------------------------------------------------
set.seed(123)

# Genrate X and expand 
#(not scaled pre expansion as it doesn't make a difference)
#----------------------------------------------------------
X0 <- matrix(runif(n = 20*6,min = -2.5,max = 2.5),nrow = 20,ncol = 6) %*% diag(exp(1:6))
grp <- rep(1:ncol(X0),each = 3)
Xp <- matrix(nrow = nrow(X0), ncol= ncol(X0)*3)
for (j in 1:ncol(X0)) {
  
  tmp <- ns(x = X0[, j],df = 3)
  Xp[,1:3 + (j-1)*3 ] <- tmp
}

# Generate beta
#--------------
sparsity <- .5 #rate of non-zero betas
active <- rep(sample(0:1,size = ncol(X0),replace = T,prob = c(1-sparsity,sparsity)),each = 3)
betas <- matrix(rnorm(n = ncol(X0)*3 + 1, mean = 5,sd = 1) * c(1,active),ncol=1)

# Generate observed data ==> y
#-----------------------------
sigma.obs <- 1
y <- cbind(1,Xp)%*% betas + sigma.obs *rnorm(nrow(X0))
mod1 <- gglasso(x = Xp,y = y,group = grp,loss = "ls",intercept = T,dfmax = )
mod2 <- gglasso(x = Xp,y = y-mean(y),group = grp,loss = "ls",intercept = F)

cf1 <- coef(mod1)
cf2 <- coef(mod2)
cf2[1,] <- mean(y)

plot(log(mod1$lambda),1:length(mod1$lambda),col=0,ylim = c(0,1),main = "Cor of soln to true beta\n Seem better off fitting on y not on y-mean(y)"
     , xlab = "log(lambda)", ylab = "Cor")
lines(log(mod1$lambda),c(cor(betas,cf1)))
lines(log(mod2$lambda),c(cor(betas,cf2)),col=2)

legend("bottomleft",legend = c("Fit with intercept", "Add y.bar as intercept"),lty = c(1,1),col=1:2)

