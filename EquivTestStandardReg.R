#########################################################################  
### function DEL()  from: Jones and Waller (2013) 
### 	"Computing Confidence Intervals for Standardized Regression Coefficients"
###		Psychological Methods
###		2013, Vol. 18, No. 4, 435– 453

DEL <- function(X = NULL, y = NULL,
cov.x = NULL, cov.xy = NULL,
                    var.y = NULL, Nobs = NULL,
                    alpha = .05, digits = 3) {
     # vech function
       vech <- function(x) t(x[!upper.tri(x)])
     # Transition or Duplicator Matrix
       Dn <- function(x){
		mat <- diag(x)
		index <- seq(x*(x+1) / 2) 
		mat[lower.tri(mat, TRUE)] <- index 
		mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
        outer(c(mat), index, function(x, y) ifelse(x == y, 1, 0))
        }
        
        
    DIAG <- function(x = 1, nrow, ncol) {
        if (length(x) == 1) 
            x <- as.matrix(x)
        if (is.matrix(x)) 
            return(diag(x))
        else return(diag(x, nrow, ncol))
    }        
        
########
#       Error Checking                      #
 if(is.null(X) & !is.null(y))
    stop("\n y is not defined\n Need to specify both X and y\n")
 if(!is.null(X) & is.null(y))
    stop("\n X is not defined\n Need to specify both X and y\n")
 if(is.null(X) & is.null(y)) {
    if(is.null(cov.x) | is.null(cov.xy) | is.null(var.y) | is.null(Nobs))
      stop("\nYou need to specify covariances and sample size\n")
    scov <- rbind(cbind(cov.x, cov.xy), c(cov.xy, var.y))
    N <- Nobs
    p <- nrow(cov.x)
  } else {
	X <- as.matrix(X)
    y <- as.matrix(y)
    scov <- cov(cbind(X, y))
    N <- length(y)
    p <- ncol(X)
  }        
        
        
# Create covariance matrix of covariances under normality
# See Browne (1984) Eqn 4.6
Kp.lft <- solve(t(Dn(p + 1)) %*% Dn(p + 1)) %*% t(Dn(p + 1))
cov.cov <- 2  * Kp.lft %*% (scov %x% scov) %*% t(Kp.lft)
param <- c(vech(scov))			# sigma_v
ncovs <- length(param)


# Find the vector element numbers for the variances of X
  v.x.pl <-c(1, rep(0, p - 1))
  for(i in 2:p) v.x.pl[i] <-v.x.pl[i - 1] + p - (i - 2)

# Store covariances and variances to use in the derivatives
  cx  <- scov[1:p, 1:p]
  cxy <- scov[1:p, p+1]
  vy  <- scov[p+1, p+1]

	sx  <- sqrt(diag(cx))
	sy  <- sqrt(vy)
	bu  <- solve(cx) %*% cxy 
	ncx <- length(vech(cx)) 

# Derivatives of standardized regression coefficients wrt the covariances.
# These are based on Yuan & Chan's (2011) equation 13
  db <- matrix(0, p, ncovs)
  V <-  matrix(0, p, ncx)
  V[as.matrix(cbind(1:p, v.x.pl))] <- 1

db[, 1:ncx] <- (											# eq (13) part 1: h_2(sigma)
		DIAG(c(solve(DIAG(2  * sx *  sy)) %*% bu)) %*% V - 
		DIAG(sx/sy) %*% (t(bu) %x% solve(cx)) %*% (Dn(p))
		)		
		
db[, (ncx+1):(ncx+p)] <- diag(sx / sy) %*% solve(cx) 	# eq (13) part 2: h_2(sigma) 
db[,ncovs] <- -diag(sx / (2  * sy^3)) %*% bu				# eq (13) part 3: h_3(sigma)


db					# hdot(s)

# Re-order the derivatives
  cx.nms <- matrix(0, p, p)
  cxy.nms <- c(rep(0, p), "var_y")

for(i in 1:p){ for(j in 1:p){ cx.nms[i, j] <- paste("cov_x", i, "x", j, sep='') }}
for(i in 1:p){ cxy.nms[i] <- paste("cov_x", i, "y", sep='') }

cx.nms
cxy.nms

old.ord <- c(vech(cx.nms), cxy.nms)
new.ord <- vech(rbind(cbind(cx.nms, cxy.nms[1:p]), c(cxy.nms)))


db <- db[, match(new.ord, old.ord)]
  
  
# Compute covariance matrix of standardized
# regression coefficients using the Delta Method
DEL.cmat <- db %*% cov.cov %*% t(db) / (N - 3) 			# hdot(s) %*% S_u %*% hdot'(s)  # eq (23)
DEL.cmat <- db %*% cov.cov %*% t(db) / (N) 
b.nms <- NULL

  for(i in 1:p) b.nms[i] <- paste("beta_", i, sep='')
  rownames(DEL.cmat) <- colnames(DEL.cmat) <- b.nms

# compute standard errors and confidence intervals
  DELse <- sqrt(diag(DEL.cmat))
  CIs <- as.data.frame(matrix(0, p, 3))
  colnames(CIs) <- c("lbound", "estimate", "ubound")
  for(i in 1:p) rownames(CIs)[i] <- paste("beta_", i, sep='')
  
  
tc <- qt(alpha / 2, N - p - 1, lower = F)
beta <- diag(sx) %*% bu * sy^-1
for(i in 1:p) CIs[i,] <- c(beta[i] - tc  * DELse[i], beta[i], beta[i] + tc *  DELse[i])
#cat("\n", 100 * (1 - alpha),"% CIs for Standardized Regression Coefficients:\n\n", sep='')
  return(list(cov.mat=DEL.cmat,SEs=DELse,alpha=alpha,CIs=CIs))
}




#################################################################################
### function equivBeta()


## Y is the response variable (vector of length N)
## Xmatrix a matrix of dimension N by K, where each column is a covariate
## DELTA is either 	- a single number for equiv. margin = [-DELTA, DELTA] for k in 1,..K
##					- a vector of length K for equiv. margin = [-DELTA[k], DELTA[k]] for k in 1,..K
##					- a matrix of dim K by 2 for equiv. margin = [DELTA[k,1], DELTA[k,2]] for k in 1,..K


equivBeta <- function(Y= rnorm(100), Xmatrix= cbind(rnorm(100),rnorm(100)), DELTA= 0.1){

X <- cbind(1,Xmatrix)
N <- dim(X[,-1])[1]
K <- dim(X[,-1])[2]

if(dim(cbind(DELTA))[1]==1){DELTA = rep(DELTA, dim(Xmatrix)[2])}
if(dim(cbind(DELTA))[2]==1){DELTA = cbind(-DELTA,DELTA)}


colnames(DELTA) <- c("Delta_1", "Delta_2")
rownames(DELTA) <- paste("k", c(0:(dim(X)[2]-1)), sep="_")

print(DELTA)

lmmod <- summary(lm(Y~X[,-1]))
beta_hat <- lmmod$coef[,1]
SE_beta_hat <- lmmod$coef[,2]
CI <- confint(lm(Y~X[,-1]), level = 0.90)

pval <- p1 <- p2 <- rep(0,K)

for(k in 1:(K+1)){ 

p1[k] <- pt((beta_hat[k] - DELTA[k,1])/SE_beta_hat[k], N-K-1, 0, lower.tail=FALSE)
p2[k] <- pt((-beta_hat[k] + DELTA[k,2])/SE_beta_hat[k], N-K-1, 0, lower.tail=FALSE)
pval[k] <- max(c(p1[k],p2[k]))

}

names(beta_hat)<-paste("beta", c(1:dim(X)[2])-1, sep="_")
return(list(beta= beta_hat, pval= pval, DELTA= DELTA, CI=CI))
}

# equivBeta(DELTA=rbind(c(-0.2,0.2),c(-0.3,0.15)))


#################################################################################
### function equivstandardBeta()


## Y is the response variable (vector of length N)
## Xmatrix a matrix of dimension N by K, where each column is a covariate
## DELTA is either 	- a single number for equiv. margin = [-DELTA, DELTA] for k in 1,..K
##					- a vector of length K for equiv. margin = [-DELTA[k], DELTA[k]] for k in 1,..K
##					- a matrix of dim K by 2 for equiv. margin = [DELTA[k,1], DELTA[k,2]] for k in 1,..K
## random is TRUE or FALSE to indicate if regressors are assumed fixed or random.

equivstandardBeta <- function(Y= rnorm(100), Xmatrix= cbind(rnorm(100),rnorm(100)), DELTA= rep(0.1, dim(Xmatrix)[2]), random= FALSE){

X <- cbind(1,Xmatrix)
N <- dim(X[,-1])[1]
K <- dim(X[,-1])[2]

if(dim(cbind(DELTA))[1]==1){DELTA = rep(DELTA, dim(Xmatrix)[2])}

if(dim(cbind(DELTA))[2]==1){DELTA = cbind(-DELTA,DELTA)}

colnames(DELTA) <- c("Delta_1", "Delta_2")
rownames(DELTA) <- paste("k", c(1:dim(Xmatrix)[2]), sep="_")

unstandard_beta <- lm(Y~X[,-1])$coefficients[-1]
sigma2_Y <- var(Y)
sigma2 <- apply(X[,-1],2, var)
standard_beta<-unstandard_beta*(c(sqrt(sigma2))/c(sqrt(sigma2_Y)))


p1 <- p2 <- pval <- R2YdotX <- R2XkdotXminK <- lambda_U2 <- lambda_L2 <- upperCI2 <- lowerCI2 <-upperCI <- lowerCI <- SE_beta_FIX <- lambda_U <- lambda_L <- rep(0,K)

b_vec <- standard_beta

if(random){ SE_std_beta_RDM <- DEL(X=Xmatrix, y=Y)$SEs }

for(k in 1:K){ 
	R2XkdotXminK[k] <- (summary(lm(X[,-1][,k]~X[,-1][,-k])))$r.squared
	R2YdotX[k]      <- (summary(lm(Y~X[,-1])))$r.squared
	SE_beta_FIX[k]  <- sqrt( (1-R2YdotX[k])/( (1-R2XkdotXminK[k])*(N-K-1)  ) )  # see Kelley2007 eq80.
#	t_vec[k] 		<- b_vec[k]/s_b_[k]	# see Kelley2007 eq79.

}

for(k in 1:K){ 
if(!random){

	p1[k] <- pt(b_vec[k]/SE_beta_FIX[k], N-K-1, DELTA[k,1]*sqrt(N*(1-R2XkdotXminK[k]))/sqrt(1-R2YdotX[k]), lower.tail=FALSE)
	p2[k] <- pt(-b_vec[k]/SE_beta_FIX[k], N-K-1, -DELTA[k,2]*sqrt(N*(1-R2XkdotXminK[k]))/sqrt(1-R2YdotX[k]), lower.tail=FALSE)
	pval[k] <- max(c(p1[k], p2[k]))

	}

if(random){
	
	p1[k] <- pt((b_vec[k] - DELTA[k,1])/SE_std_beta_RDM[k], df=N-K-1, 0, lower.tail=FALSE)
	p2[k] <- pt((DELTA[k,2] - b_vec[k])/SE_std_beta_RDM[k], df=N-K-1, 0, lower.tail=FALSE)	
	pval[k] <- max(c(p1[k], p2[k])) 	
	
	}
}

CI <- confint(lm( scale(Y)~-1+scale(X[,-1])), level=0.90)

return(list(standard_beta=standard_beta, pval=pval, DELTA= DELTA, CI=CI))
}

#equivstandardBeta(DELTA=c(0.2,0.3), random=FALSE)
#equivstandardBeta(DELTA=rbind(c(-0.2,0.3),c(-0.5,0.1)), random=TRUE)



#################################################################################
### function equivdiffP2()


## Y is the response variable (vector of length N)
## Xmatrix a matrix of dimension N by K, where each column is a covariate
## DELTA is either 	- a single number for non-inf. margin = [-Inf, DELTA] for k in 1,..K
##					- a vector of length K for non-inf. margin = [-Inf, DELTA[k]] for k in 1,..K
## random is TRUE or FALSE to indicate if regressors are assumed fixed or random.


equivdiffP2 <- function(Y= rnorm(100), Xmatrix= cbind(rnorm(100),rnorm(100)), DELTA= 0.1, random=FALSE){
	
X <- cbind(1,Xmatrix)	
N <- dim(Xmatrix)[1]
K <- dim(Xmatrix)[2]

if(length(DELTA)!=K){DELTA <- rep(DELTA[1], K) }


lmmod <- summary(lm(Y~X[,-1]))
R2 <- lmmod$r.squared
diffR2k <- unlist(lapply(c(2:(K+1)), function(k) {R2-summary(lm(Y~X[,-k]))$r.squared}))

if(random==FALSE){
pval <- rep(0, K)
for(k in 1:K){
	pval[k] <- pt(sqrt((N-K-1)*diffR2k[k])/sqrt(1-R2), N-K-1, sqrt(N*DELTA[k])/sqrt(1-R2), lower.tail=TRUE)
	}
}

if(random==TRUE){
SE_std_beta_RDM <- DEL(X = Xmatrix, y = Y)$SEs
R2XkdotXminK <- pval <- rep(0, K)
for(k in 1:K){
	R2XkdotXminK[k] <- (summary(lm(X[,-1][,k]~X[,-1][,-k])))$r.squared
	pval[k] <- pt((sqrt(diffR2k[k]) - sqrt(DELTA[k]))/ (SE_std_beta_RDM[k]*sqrt(1-R2XkdotXminK[k])), N-K-1, lower.tail=TRUE)
	}
}

return(list(diffR2k= diffR2k, pval= pval, DELTA= DELTA))
}
	      
#################################################################################	      
library("BayesFactor")
library("mvtnorm")


BFstandardBeta<-function(Y= yvec, Xmatrix= Xmat, BFthres=3, random=FALSE){

K<-dim(Xmatrix)[2]
mydata<-data.frame(Y, Xmatrix)
colnames(mydata)<- c(c("yvector"),paste("X",1:K,sep=""))
head(mydata)
BFmod <- regressionBF(yvector~. , data= mydata)

BF<-result<-rep(0,K)
for(k in 1:K){

whichk<-paste("X",k,sep="")
BF_without_k <-BFmod[!grepl(whichk,names(BFmod)$numerator)][
which.max(nchar(names(BFmod)$numerator[!grepl(whichk,names(BFmod)$numerator)]))]
BF_full <- BFmod[which.max(nchar(names(BFmod)$numerator))]
#BF[k] <- as.numeric(slot(BF_without_k/BF_full,"bayesFactor")[1])

BF[k] <- exp(as.numeric(slot(BF_without_k,"bayesFactor")[1]))/exp(as.numeric(slot(BF_full,"bayesFactor")[1]))

if(BF[k]<= 1/BFthres){result[k]<-"positive" }
if(BF[k]>BFthres){result[k]<-"negative"}	
if(BF[k]> 1/BFthres & BF[k]<BFthres){result[k]<-"inconclusive"}
}
return(list(BF=c(BF), BFthres=c(BFthres),conclusion= result))
}

	      
	      
#################################################################################
std_beta_CET<-function(Y = rnorm(100), Xmatrix = cbind(rnorm(100),rnorm(100)), 
		       k = 1, alpha1 = 0.05, random = FALSE, DELTA = rep(0.1, dim(Xmatrix)[2])){

N <- length(Y); K <- dim(Xmatrix)[2]
lmmod <- summary(lm(Y ~ Xmatrix)); R2 <- lmmod$r.squared
beta_hat <- lmmod$coef[-1,1];
R2XkXk <- unlist(lapply(c(1:K), function(k) {summary(lm(Xmatrix[,k]~ Xmatrix[,-k]))$r.squared}))
std_beta_hat <- beta_hat*(apply(Xmatrix,2,sd)/sd(Y))
SE_std_beta_FIX <- sqrt((1-R2)/((1-R2XkXk)*(N-K-1)))

pval1 <- 2*pt(abs(std_beta_hat/SE_std_beta_FIX), N-K-1, lower.tail=FALSE)[k]
pval2 <- equivstandardBeta(Y = Y , Xmatrix = Xmatrix , DELTA = DELTA , random = random)$pval[k]
if(pval1 < alpha1){return(list(pval=pval1,result=1))}
if(pval1 >= alpha1){return(list(pval=pval2,result=-1)) }
	}
	
