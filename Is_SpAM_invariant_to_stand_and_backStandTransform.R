# This script answer the following question:
# QUESTION: Does SpAM soln matches the soln
#           standardized and transformed back
#         
# Answer: The solutions follow but don't match
#         The intercept is constant for untransformed SpAM
#         Maybe that's clue
#---------------------------------------------------------

#=========================================================
#                 Utility Functions
#=========================================================


# center
#--------
centerX <- function(X){
  col.mean <- apply(X,2,mean)
  X <- sweep(X,2,col.mean)
  return(list(X=X,col.mean=col.mean))
}

# block standardize the regression
#----------------------------------
blk.standardize <- function(X,group.id,spline.deg){
  require(Matrix)
  
  p <- ncol(X)/spline.deg
  if(p!=floor(p))
    stop("number of columns of x should be a integer multiple of spline.deg")
  
  
  # initializging
  #--------------
  Qs <- array(dim = dim(X),dimnames = dimnames(X))
  R.inv <- matrix(0, spline.deg*p , spline.deg*p)
  
  
  for (gp in 1:length(unique(group.id))){
    rows <- cols <- group.id == gp
    decomp <- qr(X[, cols])
    if (decomp$rank < length(sum(cols))) 
      stop("Block belonging to columns ", paste(cols, collapse = ", "), 
           " has not full rank! \n")
    R.inv[rows,cols] <- solve(qr.R(decomp) * 1/sqrt(nrow(X)))
    Qs[, cols] <- qr.Q(decomp) * sqrt(nrow(X))
  }
  
  return(list(Qs = Qs, R.inv = R.inv))
}

# SpAM function
#dyn.load("C:/Users/Afshin/Documents/Ali/Thesis/ThLasso/SAM/libs/x64/SAM.dll")
SpAM <- function(X, y, lambda, smooth.deg){
  lmbs.len <- length(lambda)
  #centX <- centerX(X)
  #colmean.X <- centX$col.mean
  #X <- centX$X
  mn.y <- mean(y)
  y <- y - mn.y
  out = .C("grplasso", y = as.double(y), X = as.double(X), 
           lambda = as.double(lambda) , nnlambda = as.integer(lmbs.len), 
           nn = as.integer(nrow(X)), dd = as.integer(ncol(X0)), pp = as.integer(smooth.deg), 
           ww = as.double(matrix(0,ncol(X), lmbs.len)), mmax_ite = as.integer(3e+08), 
           tthol = as.double(1e-08), iinput = as.integer(1), 
           df = as.integer(rep(0, lmbs.len)), sse = as.double(rep(0, 
                                                                  lmbs.len)), func_norm = as.double(matrix(0, ncol(X)/smooth.deg, lmbs.len)), 
           package = "SAM")
  cf <- matrix(out$w, ncol = lmbs.len)
  #cf0 <- mn.y - t(cf)%*% colmean.X
  cf0 <- rep(mn.y,ncol(cf))
  cf <-rbind(t(cf0),cf)
  return(cf)
}

#=========================================================
#                 Evaluation of Qeustion
#=========================================================-

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

sparsity <- .4 #rate of non-zero betas
active <- rep(sample(0:1,size = ncol(X0),replace = T,prob = c(1-sparsity,sparsity)),each = 3)
betas <- matrix(rnorm(n = ncol(X0)*3 + 1, mean = 2,sd = 1) * c(1,active),ncol=1)
sigma.obs <- 0.5
y.obs.prenoise <- cbind(1, Xp) %*% betas


#~~~~~~~~~~~~~~~~~~~~~Standardizaiton Related~~~~~~~~~~~~~~~~~~~~~~~~~~
# Center
#--------
cent <- centerX(X = Xp)
X <- cent$X
col.mean <- cent$col.mean

# Standardize
#------------
stand <- blk.standardize(X = X,group.id = rep(1:4,each=3),spline.deg = 3)
X <- stand$Qs
R.inv <- stand$R.inv

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sim param
#----------
sim.num <- 1
set.seed(126)
seeds <- sample(1:(10*sim.num), size = sim.num)

# Fit param
#----------
lmbs.len <- 20
lmb.range <- c(0.2, 2)
lmbs <- exp(seq(from=log(lmb.range[2]),to=log(lmb.range[1]),length.out=lmbs.len))


# Result storage init
#--------------------
#DOF.est.gg.i <- DOF.est.sm.i <- array(NA,dim = c(sim.num,lmbs.len))
obs.noise <- array(NA,dim = c(sim.num,nrow(X0)))
y.hat.gg <- y.hat.sam <- array(NA,dim = c(sim.num,lmbs.len,nrow(X0)))


sn<- 1

# Data generation
#----------------
set.seed(seeds[sn])
obs.noise[sn,] <-  rnorm(n=nrow(X0)) * sigma.obs

y.obs <- y.obs.prenoise + obs.noise[sn,]


# SpAM on Orthogonalized Data
#----------------------------
cf <- SpAM(X = X, y = y.obs,lambda = sqrt(3)* lmbs,smooth.deg = 3)
rownames(cf)<- c("(Intercept)",paste("V",1:(nrow(cf)-1),sep=""))

#gglasso on Orthogonalized Data
#------------------------------
mod.gg <- gglasso(x = X,y = y.obs-mean(y.obs), group = grp, loss="ls", lambda = lmbs/sqrt(nrow(X0)),intercept = T,eps=1e-8)

cf.gg <- coef(mod.gg)
cf.gg[1,] <- mean(y.obs) - t(coef(mod.gg)[-1,])%*% apply(X,MARGIN = 2,mean)

# Transfrom back the Solns
#-------------------------
cf[-1,] <- R.inv %*% cf[-1,]
cf.gg[-1,] <- R.inv %*% cf.gg[-1,]

cf[1, ] <- cf[1,] - apply(cf[-1, , drop = FALSE] * col.mean, 2, sum)
cf.gg[1, ] <- cf.gg[1,] - apply(cf.gg[-1, , drop = FALSE] * col.mean, 2, sum)

# SpAM on original Data
#----------------------
#mod <- samQL(X = X0,y = y.obs,lambda = sqrt(3)*lmbs,thol = 1e-08)
#cf.sam <- (rbind(mod$intercept,mod$w))
cf.sam <- SpAM(X = Xp, y = y.obs,lambda = sqrt(3)* lmbs,smooth.deg = 3)
rownames(cf.sam) <- c("(Intercept)",paste("V",1:(nrow(cf)-1),sep=""))

# Comparison plot
#----------------
plot(c(cf),c(cf.sam),main="SpAM stz coefs vs unstz SpAM \n (Intercept might be a cluprit)",xlab = "Spam",ylab="gglasso")

abline(coef(lm(c(cf)~c(cf.sam))))
coefs <- as.character(signif(coef(lm(c(cf)~c(cf.sam))),3))
text(x = 2,y=3,labels = paste(coefs[1],"+",coefs[2]))

cors <- rep(NA,ncol(cf))
for(i in 1:ncol(cf)){
  cors[i] <- cor(cf[,i],cf.sam[,i]) 
}

plot(log(lmbs),cors,main="Soln cor vs. lambda\n Overall follows but deteriorates for small lambdas")
