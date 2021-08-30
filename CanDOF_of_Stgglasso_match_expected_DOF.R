# This script answer the following question:
# QUESTION: Can we get a SpAM like agreement with DOF using
#           gglasso with standardization and correct
#           (i.e. 1/sqrt(n)*lambda) lambda scaling
#         
# Answer: solve SpAM as usual with pen = sqrt(3)*lambda
#         orthogonalize, solve with gglasso with pen = lambda/sqrt(n)
#         then transform back and compute SGL DOF with expanded X
#         and lambda for both ==> We get really good agreement!!
#---------------------------------------------------------

#=========================================================
#                 Utility Functions
#=========================================================

get.GL.DOF <- function(X, beta.hat, lambda, group.id, SGL = F){
  require(Matrix)
  if(length(beta.hat) != length(group.id))
    stop("Length of beta.hat needs to be equal to length of group.id")
  tb.group <- table(group.id)
  if(!all(tb.group ==tb.group[1]))
    stop("all groups need to be the same size")
  if(any(c(group.id,group.id[length(group.id)])-c(1,group.id) < 0))
    stop("all groups need to be in blocks (not inter-dispersed) and in increasing order")
  
  active.set <- which(beta.hat != 0)
  if(length(active.set) == 0) # DOF is zero if all betas are zero
    return(0)
  groups.active <- group.id[active.set]
  bHat.active <- beta.hat[active.set]
  sz.active <- length(bHat.active)
  sz.group <- tb.group[1]
  num.group <- sz.active/sz.group
  
  # If the solver was regular group lasso
  #--------------------------------------
  if(!SGL){
    block.diagonizer <- bdiag( rep(list(matrix(1,sz.group,sz.group)),num.group))
    
    # Normalizing operator
    norm.blocks <- tapply(beta.hat^2,INDEX = group.id,FUN = sum,simplify = TRUE)[groups.active]
    normalizer.mat <- matrix(1/sqrt(norm.blocks),nrow = sz.group*num.group,ncol=sz.group * num.group,byrow = TRUE)
    
    # Projection block diagonal operator
    P_bet <- bHat.active%*%t(bHat.active) * block.diagonizer * normalizer.mat^2
    Id <- diag(1,nrow = sz.group*num.group, ncol = sz.group *num.group)
    D <- normalizer.mat * (Id - P_bet)
  }
  
  # If the solver was standardized group lasso
  #-------------------------------------------
  else{
    D <- matrix(0, nrow = sz.active, ncol = sz.active)
    for(gr in 1:num.group){
      blk.ndx <- groups.active %in% unique(groups.active)[gr]
      group.b <- which(group.id %in% groups.active[blk.ndx])
      X.gr <- X[,group.b]
      beta.gr <- beta.hat[group.b]
      y.hat.gr <- X.gr %*% beta.gr
      
      nrm2.y.hat.gr <- as.numeric(t(y.hat.gr) %*% y.hat.gr)
      D[blk.ndx, blk.ndx] <- t(X.gr) %*% X.gr / sqrt(nrm2.y.hat.gr)  - 
        t(X.gr) %*% y.hat.gr %*% t(y.hat.gr) %*% X.gr / (nrm2.y.hat.gr)^1.5
      
    }
  }
  
  
  # cov(y.hat, y)
  #--------------
  inv.term <- tryCatch(expr = solve(t(X[,active.set,drop=F]) %*% X[,active.set,drop=F] + lambda * D )
                       ,error = function(exception){
                         return(exception)})
  
  if(inherits(inv.term,"simpleError"))
    stop(paste(inv.term$message,"\n  ***Singularity usually occurs when regression was done with too smal of a lambda.If system is singular try larger lambda***"))
  
  cov.yHat.y <- X[,active.set,drop=F] %*%  inv.term %*% t(X[,active.set,drop=F])
  GL.DOF <-sum(diag(cov.yHat.y))
  return(GL.DOF)
}

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
  require(SAM)
  lmbs.len <- length(lambda)
  centX <- centerX(X)
  colmean.X <- centX$col.mean
  X <- centX$X
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
  cf0 <- mn.y - t(cf)%*% colmean.X
  #cf0 <- rep(mn.y,ncol(cf))
  cf <-rbind(t(cf0),cf)
  return(cf)
}

#=========================================================
#                 Evaluation of Qeustion
#=========================================================-
require(gglasso)
require(splines)
set.seed(126)

# Genrate X and expand 
#(not scaled pre expansion as it doesn't make a difference)
#----------------------------------------------------------
X0 <- matrix(runif(n = 30*20,min = -2.5,max = 2.5),nrow = 30,ncol = 20) %*% diag(1:20)
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
stand <- blk.standardize(X = X,group.id = grp,spline.deg = 3)
X <- stand$Qs
R.inv <- stand$R.inv

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sim param
#----------
sim.num <- 100
set.seed(126)
seeds <- sample(1:(10*sim.num), size = sim.num)

# Fit param
#----------
lmbs.len <- 20
lmb.range <- c(0.2, 5)
lmbs <- exp(seq(from=log(lmb.range[2]),to=log(lmb.range[1]),length.out=lmbs.len))


# Result storage init
#--------------------
DOF.est.gg.i <- DOF.est.sm.i <- array(NA,dim = c(sim.num,lmbs.len))
obs.noise <- array(NA,dim = c(sim.num,nrow(X0)))
y.hat.gg <- y.hat.sam <- array(NA,dim = c(sim.num,lmbs.len,nrow(X0)))


for(sn in 1:sim.num){
  # Data generation
  #----------------
  set.seed(seeds[sn])
  obs.noise[sn,] <-  rnorm(n=nrow(X0)) * sigma.obs
  
  y.obs <- y.obs.prenoise + obs.noise[sn,]
  
  
  # SpAM on Orthogonalized Data
  #----------------------------
  #cf <- SpAM(X = X, y = y.obs,lambda = sqrt(3)* lmbs,smooth.deg = 3)
  #rownames(cf)<- c("(Intercept)",paste("V",1:(nrow(cf)-1),sep=""))
  
  #gglasso on Orthogonalized Data
  #------------------------------
  mod.gg <- gglasso(x = X,y = y.obs-mean(y.obs), group = grp, loss="ls", lambda = lmbs/sqrt(nrow(X0)),intercept = T,eps=1e-8)
  
  cf.gg <- coef(mod.gg)
  cf.gg[1,] <- mean(y.obs) - t(coef(mod.gg)[-1,])%*% apply(X,MARGIN = 2,mean)
  
  # Transfrom back the Solns
  #-------------------------
  cf.gg[-1,] <- R.inv %*% cf.gg[-1,]
  cf.gg[1, ] <- cf.gg[1,] - apply(cf.gg[-1, , drop = FALSE] * col.mean, 2, sum)
  
  
  # SpAM on original Data
  #----------------------
  #mod <- samQL(X = X0,y = y.obs,lambda = sqrt(3)*lmbs,thol = 1e-08)
  #cf.sam <- (rbind(mod$intercept,mod$w))
  cf.sam <- SpAM(X = Xp, y = y.obs,lambda = sqrt(3)* lmbs,smooth.deg = 3)
  rownames(cf.sam) <- c("(Intercept)",paste("V",1:(nrow(cf.sam)-1),sep=""))
  
  y.hat.gg[sn,,] <- t(cbind(1,Xp) %*% cf.gg)
  #y.hat.sam[sn,,] <- t(predict(mod,newdata = X0)$values)
  y.hat.sam[sn,,] <- t(cbind(1,Xp) %*% cf.sam)
  
  ## Evaluating the concordance of cf.gg and cf.sam
  ##-----------------------------------------------
  # plot(c(cf.gg),c(cf.sam))
  
  
  # Fitting and DOF
  #----------------
  for(i in 1:lmbs.len){
    
    #DOF estimate
    #-------------
    DOF.est.gg.i[sn,i] <- 1 + get.GL.DOF(X = Xp, beta.hat = cf.gg[-1,i], lambda = lmbs[i],group.id = grp,SGL = T)
    DOF.est.sm.i[sn,i] <- 1 + get.GL.DOF(X = Xp, beta.hat = cf.sam[-1,i], lambda = lmbs[i],group.id = grp,SGL = T)
  }# over lambda
  
}

# Compute y.hat.mean and set dim to c(sim.num,lmbs.len,n.obs)
#-------------------------------------------------------------
y.hat.mean <- aperm(array(y.obs.prenoise,dim = c(sim.num,nrow(X0),lmbs.len)),perm = c(1,3,2))


# Set dim to c(sim.num,lmbs.len,n.obs)
#------------------------------------
obs.noise <- aperm(array(obs.noise,dim = c(sim.num,nrow(X0),lmbs.len)),perm = c(1,3,2))


# Compute DOFs (ture)
#--------------------
mult.gg <- (y.hat.gg - y.hat.mean)*(obs.noise)
cov.gg.i <- apply(mult.gg,MARGIN = c(2,3),mean)               # Mean over simulations
DOFs.gg <- apply(cov.gg.i,MARGIN = 1,FUN = sum)/(sigma.obs^2) # Sum over observation number n


mult.sm <- (y.hat.sam - y.hat.mean)*(obs.noise)
cov.sm.i <- apply(mult.sm,MARGIN = c(2,3),mean)               # Mean over simulations
DOFs.sm <- apply(cov.sm.i,MARGIN = 1,FUN = sum)/(sigma.obs^2) # Sum over observation number n


# Compute Est. DOF
#-----------------
DOF.ests.gg <- apply(DOF.est.gg.i,MARGIN =2,FUN = mean )
DOF.ests.sm <- apply(DOF.est.sm.i,MARGIN =2,FUN = mean )


# Plot
#-----
plot(DOFs.gg,DOFs.gg,ylim=range(c(DOFs.gg,DOF.ests.gg,DOFs.sm,DOF.ests.sm)),col=0,
     main = "Est. DOF with standardized gglasso", xlab = "True DOF", ylab = "Estimated DOF")
lines(DOFs.gg,DOFs.gg,lty=2)
points(DOFs.gg,DOF.ests.gg,col=2, pch = 4,lwd =2)
points(DOFs.sm,DOF.ests.sm,col=1, pch = 1,cex=2)
legend("bottomright",legend = c("Standardized gglasso", "SpAM","45 degree line")
       ,pt.cex = c(1,2,1),pch=c(4,1,NA),lwd=c(2,1,1),col=c(2,1,1),lty=c(NA,NA,2))


#plot(apply(cf.gg,2,function(col){sum(col[betas!=0]!=0)/length(betas)}))
#lines(apply(cf.gg,2,function(col){sum(col[betas==0]!=0)/length(betas)}),col=2)