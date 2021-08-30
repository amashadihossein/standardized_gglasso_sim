# This script answer the following question:
# QUESTION: Does scaling pre expansion matter?
#----------------------------------------------

X0 <- matrix(runif(n = 5*6,min = -2.5,max = 2.5),nrow = 5,ncol = 6) %*% diag(exp(1:6))

# z scale
#--------
X.zscale <- scale(X0)
  
# unit scale
#-----------
rng <- apply(apply(X0,2,range),2,diff)
mn <- apply(X0,2,min)
X.uscale <- sweep(sweep(X0,MARGIN = 2,STATS = mn),MARGIN = 2,STATS = rng,FUN = "/")

Xp0 <- Xpz <- Xpu <- matrix(nrow = nrow(X0), ncol= ncol(X0)*3)
for (j in 1:ncol(X0)) {
  
  tmp0 <- ns(x = X0[, j],df = 3)
  tmpz <- ns(x = X.zscale[, j], df =  3)
  tmpu <- ns(x = X.uscale[, j],df = 3)
  
  Xp0[,1:3 + (j-1)*3 ] <- tmp0
  Xpz[,1:3 + (j-1)*3 ] <- tmpz
  Xpu[,1:3 + (j-1)*3 ] <- tmpu
}


if(all.equal(Xp0, Xpz,check.attributes = F) & all.equal(Xp0, Xpu,check.attributes = F)){
  cat("scale pre expansion doesn't make a difference")
}else{
  cat("scale pre expansion makes a difference")
}

