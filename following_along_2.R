
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/harlanhappydog/EquivTestStandardReg/master/EquivTestStandardReg.R", ssl.verifypeer = FALSE)

eval(parse(text = script))


## Example data: ##
set.seed(123)
y <- rnorm(100)
X <- cbind(1, rnorm(100), rpois(100,4))


## PART 1 ##
lmmod <- summary(lm(y~X[,-1]))
N <- length(y); K <- dim(X[,-1])[2]
beta_hat <- lmmod$coef[,1]; SE_beta_hat <- lmmod$coef[,2]
pval <- 2*pt(abs(beta_hat/SE_beta_hat), N-K-1, 0, lower.tail=FALSE)
pval

# alternatively:
summary(lm(y~X[,-1]))$coef[,4]

## PART 2 ##
R2 <- lmmod$r.squared
diffR2k <- unlist(lapply(c(2:(K+1)), function(k) {R2-summary(lm(y~X[,-k]))$r.squared}))
pval <- pf((N-K-1)*(diffR2k/(1-R2)), 1, N-K-1, 0, lower.tail=FALSE)
pval

# alternatively:
summary(lm(y~X[,-1]))$coef[,4][-1]

## PART 3 ##
DELTA <- rbind(c(-0.3, 0.25), c(-0.2, 0.2), c(-0.3, 0.1))
pval <- p1 <- p2 <- rep(0, K+1)
for(k in 1:(K+1)){
	p1[k] <- pt((beta_hat[k] - DELTA[k,1])/SE_beta_hat[k], N-K-1, 0, lower.tail=FALSE)
	p2[k] <- pt((-beta_hat[k] + DELTA[k,2])/SE_beta_hat[k], N-K-1, 0, lower.tail=FALSE)
	pval[k] <- max(c(p1[k],p2[k]))
}
pval

# alternatively:
equivBeta(Y = y, Xmatrix = X[,-1], DELTA = DELTA)$pval

## PART 4 ##
b_vec <- (beta_hat*(apply(X,2,sd)/sd(y)))[-1]
b_vec

# alternatively:
lm(scale(y) ~ scale(X[,-1]))$coef[-1]

## PART 5 ##
DELTA <- rbind(c(-0.2, 0.2), c(-0.3, 0.1))
SE_beta_FIX <- R2YdotX <- R2XkdotXminK <- pval <- p1 <- p2 <- rep(0, K)
for(k in 1:K){
	R2XkdotXminK[k] <- (summary(lm(X[,-1][,k]~X[,-1][,-k])))$r.squared
	R2YdotX[k]      <- (summary(lm(y~X[,-1])))$r.squared
	SE_beta_FIX[k]  <- sqrt( (1-R2YdotX[k])/( (1-R2XkdotXminK[k])*(N-K-1)  ) )  

	p1[k] <- pt(b_vec[k]/SE_beta_FIX[k], N-K-1,
	             DELTA[k,1]*sqrt(N*(1-R2XkdotXminK[k]))/sqrt(1-R2YdotX[k]), lower.tail=FALSE)
	p2[k] <- pt(-b_vec[k]/SE_beta_FIX[k], N-K-1,
              	-DELTA[k,2]*sqrt(N*(1-R2XkdotXminK[k]))/sqrt(1-R2YdotX[k]), lower.tail=FALSE)
	pval[k] <- max(c(p1[k], p2[k]))
}	
pval
	
# alternatively:
equivstandardBeta(Y = y, Xmatrix = X[,-1], DELTA = DELTA, random = FALSE)$pval

## PART 6 ##
SE_std_beta_RDM <- DEL(X=X[,-1], y=y)$SEs
pval <- p1 <- p2 <- rep(0, K)
for(k in 1:K){
	p1[k] <- pt((b_vec[k] - DELTA[k,1])/SE_std_beta_RDM[k], N-K-1, 0, lower.tail=FALSE)
	p2[k] <- pt((DELTA[k,2] - b_vec[k])/SE_std_beta_RDM[k], N-K-1, 0, lower.tail=FALSE)	
	pval[k] <- max(c(p1[k], p2[k]))
}
pval

# alternatively:
equivstandardBeta(Y = y, Xmatrix = X[,-1], DELTA = DELTA, random = TRUE)$pval



## PART 7 ##
DELTA <- rep(0.05,2)
pval <- rep(0, K)
for(k in 1:K){
	pval[k] <- pt(sqrt((N-K-1)*diffR2k[k])/sqrt(1-R2), N-K-1, sqrt(N*DELTA[k])/sqrt(1-R2), lower.tail=TRUE)
}
pval

# alternatively:
equivdiffP2(Y = y, Xmatrix = X[,-1], DELTA = DELTA, random = FALSE)$pval

## PART 8 ##
DELTA <- rep(0.05,2)
pval <- rep(0, K)
for(k in 1:K){
pval[k] <- pt((sqrt(diffR2k[k]) - sqrt(DELTA[k]))/ (SE_std_beta_RDM[k]*sqrt(1-R2XkdotXminK[k])), N-K-1, lower.tail=TRUE)
}
pval

# alternatively:
equivdiffP2(Y = y, Xmatrix = X[,-1], DELTA = DELTA, random = TRUE)$pval

#########################################################################  
