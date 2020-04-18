




#################################################################################
### function equivBeta()


## Y is the response variable (vector of length N)
## Xmatrix a matrix of dimension N by K, where each column is a covariate
## DELTA is either 	- a single number for equiv. margin = [-DELTA, DELTA] for k in 1,..K
##					- a vector of length K for equiv. margin = [-DELTA[k], DELTA[k]] for k in 1,..K



equivBeta <- function(Y= rnorm(100), Xmatrix= cbind(rnorm(100),rnorm(100)), DELTA= 0.1){
  
  Xmatrix <-cbind(Xmatrix)	
  DELTA<-cbind(DELTA)
  X <- cbind(1,Xmatrix)
  N <- dim(cbind(X[,-1]))[1]
  K <- dim(cbind(X[,-1]))[2]
  
  if(sum(dim(matrix(DELTA)))==2){ DELTA = cbind(-DELTA,DELTA) }
  DELTA <- t(matrix(rep(DELTA,(K+1)),2,))
  colnames(DELTA) <- c("Delta_1", "Delta_2")
  rownames(DELTA) <- paste("k", c(0:(dim(X)[2]-1)), sep="_")
  
  
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


# equivBeta(DELTA=rbind(c(-0.2,0.2),c(-0.2,0.2),c(-0.3,0.15)))




#################################################################################
### function equivstandardBeta()


## Y is the response variable (vector of length N)
## Xmatrix a matrix of dimension N by K, where each column is a covariate
## DELTA is either 	- a single number for equiv. margin = [-DELTA, DELTA] for k in 1,..K
##					- a vector of length K for equiv. margin = [-DELTA[k], DELTA[k]] for k in 1,..K

equivstandardBeta <- function(Y= rnorm(100), Xmatrix= cbind(rnorm(100),rnorm(100)), DELTA= rep(0.1, dim(Xmatrix)[2]), kvec=1:K){
Xmatrix <-cbind(Xmatrix)	
X <- cbind(1,Xmatrix)
N <- dim(cbind(X[,-1]))[1]
K <- dim(cbind(X[,-1]))[2]

if(sum(dim(matrix(DELTA)))==2){ DELTA = cbind(-DELTA,DELTA) }
DELTA <- t(matrix(rep(DELTA,K),2,))

colnames(DELTA) <- c("Delta_1", "Delta_2")
rownames(DELTA) <- paste("k", c(1:dim(Xmatrix)[2]), sep="_")

unstandard_beta <- lm(Y~ cbind(X[,-1]))$coefficients[-1]
sigma2_Y <- var(Y)
sigma2_X <- apply(cbind(X[,-1]),2, var)
standard_beta<-unstandard_beta*(c(sqrt(sigma2_X))/c(sqrt(sigma2_Y)))
#  lm(scale(Y)~ scale(cbind(X[,-1])))$coefficients[-1]

pval_fix<- pval_rdm <- Kupper<- Klower <- p1 <- p2 <- pval <- R2YdotXmink <- R2YdotX <- R2XkdotXminK <- lambda_U2 <- lambda_L2 <- upperCI2 <- lowerCI2 <-upperCI <- lowerCI <- SE_beta_FIX <- lambda_U <- lambda_L <- rep(0,K)

b_vec <- standard_beta

for(k in kvec){ 
	
if(K>1){	Xmink <- cbind(cbind(X[,-1])[,-k])}
if(K==1){	Xmink <- rep(1,N)}

	R2YdotXmink[k]      <- summary(lm(Y~ Xmink))$r.squared	
	R2XkdotXminK[k] <- (summary(lm(cbind(X[,-1])[,k]~ Xmink)))$r.squared
	R2YdotX[k]      <- (summary(lm(Y~ cbind(X[,-1]))))$r.squared
	SE_beta_FIX[k]  <- sqrt( (1-R2YdotX[k])/( (1-R2XkdotXminK[k])*(N-K-1)  ) )  # see Kelley2007 eq80.
#	t_vec[k] 		<- b_vec[k]/s_b_[k]	# see Kelley2007 eq79.


}

for(k in kvec){ 


P2_1 = (1-R2XkdotXminK[k])*DELTA[k,1]^{2} + R2YdotXmink[k]
ncp_1 = sqrt(N*(1-R2XkdotXminK[k])) * (DELTA[k,1]/sqrt( 1 - P2_1 ))

P2_2 = (1-R2XkdotXminK[k])*DELTA[k,2]^{2} + R2YdotXmink[k]
ncp_2 = sqrt(N*(1-R2XkdotXminK[k])) * (-DELTA[k,2]/sqrt( 1 - P2_1 ))

	p1[k] <- pt(b_vec[k]/SE_beta_FIX[k], N-K-1, ncp=ncp_1, lower.tail=FALSE)
	p2[k] <- pt(-b_vec[k]/SE_beta_FIX[k], N-K-1, ncp=ncp_2, lower.tail=FALSE)
	pval[k] <- max(c(p1[k], p2[k]))

}

CI <- confint(lm( scale(Y)~-1+scale(cbind(X[,-1]))), level=0.90)
return(list(standard_beta=standard_beta, pval = pval, DELTA= DELTA, CI=CI))
}

#################################################################################
### function equivdiffP2()


## Y is the response variable (vector of length N)
## Xmatrix a matrix of dimension N by K, where each column is a covariate
## DELTA is either 	- a single number for non-inf. margin = [-Inf, DELTA] for k in 1,..K
##					- a vector of length K for non-inf. margin = [-Inf, DELTA[k]] for k in 1,..K


equivdiffP2 <- function(Y= rnorm(100), Xmatrix= cbind(rnorm(100),rnorm(100)), DELTA= 0.1, kvec=1:K){
Xmatrix <-cbind(Xmatrix)		
X <- cbind(1,Xmatrix)	
N <- dim(Xmatrix)[1]
K <- dim(Xmatrix)[2]

if(length(DELTA)!=K){DELTA <- rep(DELTA[1], K) }


lmmod <- summary(lm(Y~X[,-1]))
R2 <- lmmod$r.squared
diffR2k <- unlist(lapply(c(2:(K+1)), function(k) {R2-summary(lm(Y~X[,-k]))$r.squared}))

R2XkdotXminK <- pval <- rep(0, K)
for(k in kvec){
	if(K>1){		Xmink <- cbind(cbind(X[,-1])[,-k])}
	if(K==1){	Xmink <- rep(1,N)}
	R2XkdotXminK[k] <- (summary(lm(cbind(X[,-1])[,k]~ Xmink)))$r.squared
	ncp_1<-sqrt(N*DELTA[k])/sqrt(1-DELTA[k]+  R2XkdotXminK[k])
	pval[k] <- pt(sqrt((N-K-1)*diffR2k[k])/sqrt(1-R2), N-K-1, ncp=ncp_1, lower.tail=TRUE)
	}

return(list(diffR2k= diffR2k, pval= pval, DELTA= DELTA))
}



#################################################################################

noninfP2 <- function(Y= rnorm(100), Xmatrix= cbind(rnorm(100),rnorm(100)), DELTA= 0.1, random=FALSE, tol=1.0e-12){
Xmatrix <-cbind(Xmatrix)	
X <- cbind(1,Xmatrix)	
n<-N <- dim(Xmatrix)[1]
k<-K <- dim(Xmatrix)[2]
delta<-DELTA

lmmod <- summary(lm(Y~X[,-1]))
Rsq <- R2 <- lmmod$r.squared

if(random==FALSE){
	Fstat 	<- (Rsq/k)/((1-Rsq)/(n-k-1))
	pval 	<- pf(Fstat, df1 = k, df2 = n-k-1, ncp=(n*delta)/(1-delta), lower.tail = TRUE) 
	
}

if(random==TRUE){
Psq <- Rsq; Psq_last <- 1; # initial value

    F_num 	 <- (n-k-1)*Rsq*(delta-1)
	F_den 	 <- ((Rsq-1) * (delta*(n-k-1) + k))
	Fstat 	 <- F_num/F_den

while(abs(Psq_last - Psq) > tol){
    Psq_last <- Psq
    v 		 <- (((n-k-1)*Psq + k)^2)/(n-1-(n-k-1)*(1-Psq)^2)
    Psq_num  <- (n-k-1)*Rsq - (1-Rsq)*k*Fstat
	Psq_den  <- (n-k-1)*(Rsq + (1-Rsq)*Fstat)
	Psq 		 <- Psq_num/Psq_den
}
pval <- pf(Fstat, v, n-k-1, lower.tail=TRUE)}




return(list(R2 = R2, pval= pval, DELTA= DELTA))
}



	
	
	      
#################################################################################	      
library("BayesFactor")
library("mvtnorm")


BFstandardBeta<-function(Y= yvec, Xmatrix= Xmat, BFthres=3, kvec=c(1:K)){
  Xmatrix<-cbind(Xmatrix)
  K<-dim(Xmatrix)[2]; print(K)
  mydata<-data.frame(Y, Xmatrix)
  colnames(mydata)<- c(c("yvector"),paste("X",1:K,sep=""))
  head(mydata)
  BFmod <- regressionBF(yvector~. , data= mydata)
  BF<-result <- 0*kvec
  if(K>1){
  for(k in kvec){
    
    whichk<-paste("X",k,sep="")
    BF_without_k <-BFmod[!grepl(whichk,names(BFmod)$numerator)][
      which.max(nchar(names(BFmod)$numerator[!grepl(whichk,names(BFmod)$numerator)]))]
    BF_full <- BFmod[which.max(nchar(names(BFmod)$numerator))]
    #BF[k] <- as.numeric(slot(BF_without_k/BF_full,"bayesFactor")[1])
    
    BF[k] <- exp(as.numeric(slot(BF_without_k,"bayesFactor")[1]))/exp(as.numeric(slot(BF_full,"bayesFactor")[1]))
  } 
 }
  
  if(K==1){
    BF<-exp(slot(regressionBF(yvector~. , data= mydata), "bayesFactor")[1])
  }
  for(k in kvec){
    if(BF[k]<= 1/BFthres){result[k]<-"positive" }
    if(BF[k]>BFthres){result[k]<-"negative"}	
    if(BF[k]> 1/BFthres & BF[k]<BFthres){result[k]<-"inconclusive"}
  }
    return(list(BF=c(BF), BFthres=c(BFthres),conclusion= result))
  }


	      
	      
#################################################################################
std_beta_CET<-function(Y = rnorm(100), Xmatrix = cbind(rnorm(100),rnorm(100)), 
		       k = 1, alpha1 = 0.05,  DELTA = rep(0.1, dim(Xmatrix)[2])){

N <- length(Y); K <- dim(Xmatrix)[2]
lmmod <- summary(lm(Y ~ Xmatrix)); R2 <- lmmod$r.squared
beta_hat <- lmmod$coef[-1,1];
R2XkXk <- unlist(lapply(c(1:K), function(k) {summary(lm(Xmatrix[,k]~ Xmatrix[,-k]))$r.squared}))
std_beta_hat <- beta_hat*(apply(Xmatrix,2,sd)/sd(Y))
SE_std_beta_FIX <- sqrt((1-R2)/((1-R2XkXk)*(N-K-1)))

pval1 <- 2*pt(abs(std_beta_hat/SE_std_beta_FIX), N-K-1, lower.tail=FALSE)[k]
pval2 <- equivstandardBeta(Y = Y , Xmatrix = Xmatrix , DELTA = DELTA )$pval[k]
if(pval1 < alpha1){return(list(pval=pval1,result=1))}
if(pval1 >= alpha1){return(list(pval=pval2,result=-1)) }
	}




###################################################################################
linear_reg<- function(y, X, DELTA=0.1, BF_thres=6, alpha = 0.05){

lmsummary<-summary((lm(salary ~ rank + discipline + yrs.since.phd + yrs.service + sex, data=Salaries)))

	CETfunction<-function(x) { if(x[1]<0.05){return("postive")}
				if(x[1]>=0.05 & x[2]<0.05){return("negative")}
					if(x[1]>=0.05 & x[2]>=0.05){return("inconclusive")}
				}
	
equivBeta <- suppressWarnings(equivstandardBeta(Y = y, Xmatrix = X[,-1], DELTA = DELTA))
CET_summary_table <- data.frame(beta=lmsummary$coef[,1], std_beta= c(NA,equivBeta$standard_beta), NHST_p=lmsummary$coef[,4], equiv_p=c(NA,equivBeta$pval) ,
CET_con=c(NA,apply(cbind(lmsummary$coef[,4], c(NA,equivBeta$pval))[-1,],1, CETfunction)))									
				
BFsumm <- BFstandardBeta(Y= y, Xmatrix=X[,-1], BFthres= BF_thres)
BF_summary_table <- data.frame(c(NA,BFsumm$BF), BF_con=c(NA,BFsumm$conclusion))

return(print(data.frame(CET_summary_table ,BF_summary_table), digits=3))}
	
