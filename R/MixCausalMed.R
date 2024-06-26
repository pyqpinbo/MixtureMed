#' Integration of decision trees with mixture mediation model
#'
#' @param data The data frame contains treatment variable, confounding variables, mediator variable and outcome variable.
#' @param K The number of mixing components
#' @keywords Model estimation for mixture of mediation model
#' @import e1071
#' @import rpart
#' @import stats
#'
#' @export
#' 
#' @examples
#' R code examples of how to use your function
#'

MixCausalMed <- function(data, K) {
  X = data$X
  M = data$M
  Z = data$Z
  Y = data$Y
  n = length(Y)
  U<-rep(1, n)
  V<-matrix(c(U, data$X,data$Z), ncol=3)
  V1<-matrix(c(U, data$X, data$M,data$Z), ncol=4)
  beta1<-runif(K)
  beta2<-runif(K)
  a<-runif(K)
  b<-runif(K)
  c<-runif(K)
  ceta1 <-  runif(K)
  ceta2 <-  runif(K)
  sigma1<-runif(K, 0.5, 1)
  sigma2<-runif(K, 0.5, 1)
  piz_est <-  stats::runif(n,0,1)
  
  iter = 0
  epsilon <- 1e-4
  library(rpart)
  library(e1071)
  
   repeat {
    iter = iter + 1
    ## E-step
    #
    mean1 <- beta1[1]+a[1]*X+ceta1[1]*Z
    f1 = stats::dnorm(M, mean = mean1, sd = sigma1[1])
    mean11 <-  beta2[1]+c[1]*X+b[1]*M+ceta2[1]*Z
    f11 = stats::dnorm(Y, mean = mean11, sd = sigma2[1])
    F1 = f1*f11
    #
    mean2 <- beta1[2]+a[2]*X+ceta1[2]*Z
    f2 =stats::dnorm(M, mean = mean2, sd = sigma1[2])
    mean22 <-  beta2[2]+c[2]*X+b[2]*M+ceta2[2]*Z
    f22 = stats::dnorm(Y, mean = mean22, sd = sigma2[2])
    F2 = f2*f22
    
    w <- piz_est*F1 / (piz_est*F1 + (1-piz_est)*F2)
    H = cbind(F1,F2)
    
    oldbeta1<-beta1
    oldbeta2<-beta2
    olda<-a
    oldb<-b
    oldc<-c
    oldpiz_est<-piz_est
    oldsigma1<-sigma1
    oldsigma2<-sigma2
    oldceta1 <- ceta1
    oldceta2 <- ceta2
    
    ## M-step based on DT method
    
    multiplepiz = matrix(0,nrow = n,ncol = 1)
    for (l in 1:n){
      multiplepiz[l,] <- rbinom(1,size = 1,prob = w[l])
      }
    
    pizprob1 <-  as.factor(multiplepiz[,1])
    update_piz <- c(1,1) #average step
    
    data1 <- as.data.frame(cbind(X,Z,pizprob1))
    
    #tune parameter 
    obj3 <- e1071::tune.rpart(pizprob1 ~ X+Z, data = data1, minsplit = c(11, 20, 25), cp = c(0.001, 0.005, 0.1))
    bc <- obj3$best.parameters[1]
    bg <- obj3$best.parameters[2]
    
    #tree1
    mod1 <- rpart::rpart(pizprob1 ~ X + Z , data = data1, method = 'class', control = rpart.control(minsplit = bc[[1]],
                  minbucket = round(bc[[1]]/3), cp = bg[[1]]), xval = 10, parms = list(split = 'gini'))
    
    cp.min1 <-  mod1$cptable[which.min(mod1$cptable[,"xerror"]),"CP"]
    tree1 <- rpart::prune(mod1, cp = cp.min1)
    proba1 <- predict(tree1, newdata = data1, type = 'prob')
    
    
    update_piz1<- c(1,1)
    #update_pred1  <- c(1,1)
    
    for (i in 1:n) {
      update_piz1[i] <- proba1[i,colnames(proba1)== 1] 
    }
    
    pizprob1 <- as.numeric(as.character(pizprob1))
    
    #average
    for (i in 1:n) {
      update_piz <- update_piz1
    }
    
   
    
    piz_est <- update_piz
    #
    model_m <- stats::lm(M ~ X + Z, weights = w, data = data)
    beta1[1] <- stats::coef(model_m)[1]
    a[1]<- stats::coef(model_m)[2]
    ceta1[1]<- stats::coef(model_m)[3]
    
    model_y <- stats::lm(Y ~ X + M + Z, weights = w, data = data)
    beta2[1] <- stats::coef(model_y)[1]
    c[1] <- stats::coef(model_y)[2]
    b[1] <- stats::coef(model_y)[3]
    ceta2[1] <- stats::coef(model_y)[4]
    sigma1[1]<-{(M-c(beta1[1], a[1],ceta1[1])%*%t(V))%*%diag(w)%*%t(M-c(beta1[1], a[1],ceta1[1])%*%t(V))/sum(w)}^0.5
    sigma2[1]<-{(Y-c(beta2[1], c[1], b[1],ceta2[1])%*%t(V1))%*%diag(w)%*%t((Y-c(beta2[1], c[1], b[1],ceta2[1])%*%t(V1)))/sum(w)}^0.5
    
    #*
    model_m2 <- stats::lm(M ~ X + Z, weights = 1-w, data = data)
    beta1[2] <- stats::coef(model_m2)[1]
    a[2] <- stats::coef(model_m2)[2]
    ceta1[2] <- stats::coef(model_m2)[3]
    
    model_y2 <- stats::lm(Y ~ X + M + Z, weights = 1-w, data = data)
    beta2[2] <- stats::coef(model_y2)[1]
    c[2] <- stats::coef(model_y2)[2]
    b[2] <- stats::coef(model_y2)[3]
    ceta2[2] <- stats::coef(model_y2)[4]
    sigma1[2]<-{(M-c(beta1[2], a[2],ceta1[2])%*%t(V))%*%diag(1-w)%*%t(M-c(beta1[2], a[2],ceta1[2])%*%t(V))/sum(1-w)}^0.5
    sigma2[2]<-{(Y-c(beta2[2], c[2], b[2],ceta2[2])%*%t(V1))%*%diag(1-w)%*%t((Y-c(beta2[2], c[2], b[2],ceta2[2])%*%t(V1)))/sum(1-w)}^0.5
    
    
    
    
    em <- sum((abs(piz_est- oldpiz_est) <  epsilon) )
    
    em <- as.numeric(em)
    
    if ( sum(abs(b-oldb)< epsilon) & sum(abs(a-olda)< epsilon)  || iter > 10000) { 
      break 
    }
    
    cat('iter', iter, 'a', a, 'b', b, 'c', c, 'sigma1',sigma1,'sigma2',sigma2,'ceta1',ceta1,'ceta2',ceta2,'\n')
  } 
  out = list(NIE = cbind(a[1]*b[1],a[2]*b[2]), piz = piz_est)
} 
