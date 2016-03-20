library(lattice)
library(glmnet)

###########################
####Deep's LP Functions####
###########################

#This is an edited version of the LP.Score.fun() that fills in NAs for 
#discrete data when m is larger than the number of unique data points
LP.Score.fun2 <- function(x,m){ 
  u <- (rank(x,ties.method = c("average")) - .5)/length(x) 
  n <- min(length(unique(u ))-1, m )
  S.mat <- as.matrix(poly(u ,df=n)) 
  if (ncol(S.mat) < m) {
    cols <- matrix(rep(0, (length(x)*(m-n))), nrow= length(x))
    S.mat <- cbind(S.mat, cols)
  }
  return(as.matrix(scale(S.mat)))
}


LP.smooth <- function(CR,n,method){ ###--"AIC" or "BIC"
  CR.s <- sort(CR^2,decreasing=TRUE,index=TRUE)$x
  aa <- rep(0,length(CR.s))
  if(method=="AIC"){ penalty <- 2}
  if(method=="BIC"){ penalty <- log(n)}
  aa[1] <- CR.s[1] - penalty/n
  if(aa[1]< 0){ return(rep(0,length(CR))) }
  for(i in 2: length(CR.s)){
    aa[i] <- aa[(i-1)] + (CR.s[i] - penalty/n)
  }
  #plot(aa,type="b",ylab=method,cex.axis=1.2,cex.lab=1.2)
  CR[CR^2<CR.s[which(aa==max(aa))]] <- 0
  return(CR)
}

#---------------------------------------------------------------------------


############################
####Brendan's Functions#####
############################



#This function creates the T matrix, but does not use pairwise interactions
LPT <- function(x, m) {
  xMat <- x
  colnames(xMat) <- seq(1:ncol(xMat))
  tvals <- lapply(xMat, LP.Score.fun2, m=m)
  for (i in 1:ncol(xMat)) {
    colnames(tvals[[i]])<- paste0(colnames(xMat)[i], "_", 1:m)
  }
  tmat <- do.call(cbind, tvals)
  return(tmat)
}




#This function will create a levelplot() of the LP Feature correlations 
#for significant features found after LP.Smooth. The function takes for inputs
#your data, your target vector, and m to be used in LP.Score.fun2

#This plot may be illegible if you have many significant features after smoothing

LPplot <- function(x, y, m) {
  t <- LPT(x, m)
  LPcorr <- cor(y, t)
  Tsmooth <- LP.smooth(LPcorr, length(y), "AIC")
  reducedT <- t[, which(Tsmooth != 0)]
  reducedsplit <- strsplit(colnames(reducedT), "_")
  tmatsplit <- strsplit(colnames(t), "_")
  tmatnames <- sapply(tmatsplit, "[", 1)
  reducednames <- sapply(reducedsplit, "[", 1)
  sigfeat <- t[,tmatnames %in% reducednames]
  sigcor <- cor(y, sigfeat)
  matcor <- matrix(sigcor^2, nrow=m, ncol=length(unique(reducednames)))
  matcor[is.na(matcor)] <- 0
  cp <- colorRampPalette(c("red", "orange", "white"))
  levelplot(matcor, ylab = "Significant Features", xlab = "LP Moments", 
            col.regions=cp, main="LP Feature Correlations", aspect="fill", scales=list(
              y=list(
                at=seq(1:length(unique(reducednames))),
                labels=colnames(x[,as.numeric(unique(reducednames))]),
                cex=.5
              ),
              x=list(
                at=seq(1:m),
                labels= c(1:m)
              )
            )
  )
}


#This function runs cross validated lasso on your data. Returns a list with 
#the training data Tmatrix, a vector of predictions, and Tmatrix after smoothing

LP.logistic <- function(train, target, test, m) {
  t_train <- LPT(train, m)
  LPcorr <- cor(target, t_train)
  Tsmooth <- LP.smooth(LPcorr, length(target), "AIC")
  reducedT <- t_train[, which(Tsmooth != 0)]
  t_test <- LPT(test, m)
  testsmooth <- t_test[, which(Tsmooth != 0)]
  logistic <- cv.glmnet(reducedT, target, family="multinomial")
  prediction <- as.factor(predict(logistic, testsmooth,
                                  type="class", s="lambda.min"))
  output <- list("Tmatrix" = t_train, "Prediction" = prediction, 
                 "Smooth" = reducedT)
  return(output)
}