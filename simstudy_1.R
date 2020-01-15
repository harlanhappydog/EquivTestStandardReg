###  Code for Simulation Study 1 ###



######################################## 
simstudy1 <- function(randomreg=FALSE, nSim=100){

resultsmat <- data.frame(
JJJindex=NA,
true_std_beta= NA, 
power=NA, 
as.matrix(expand.grid(Nvar=c(2,4), 
	sigma2=c(0.05, 0.15, 0.5, 0.50001), 
	Nsample=c(180, 540, 1000, 3500))))

resultsmat$true_std_beta <- NA
resultsmat$alpha_sig <- 0.05

resultsmat<-resultsmat[order(resultsmat$Nvar, resultsmat$sigma2),]
pval_list<-list()
problist<-list()


for(jjj in 1:dim(resultsmat)[1]){

print(jjj)

resultsmat[jjj, "JJJindex"] <- paste("jjj",jjj,sep="")
sigma2 	<- resultsmat[jjj, "sigma2"] 
N 		<- resultsmat[jjj, "Nsample"]
nVar 	<- resultsmat[jjj, "Nvar"]

basematrix <- data.frame(expand.grid(
				X1=c(0,1),
				X2=c(0,1),
				X3=c(0,1),
				X4=c(0,1)))

X1 <- rep(basematrix$X1,900000)
X2 <- rep(basematrix$X2,900000)
X3 <- rep(basematrix$X3,900000)
X4 <- rep(basematrix$X4,900000)


if(nVar==2){
	epsilon <- rnorm(length(X1), 0, sqrt(sigma2))
	X  <- as.matrix(cbind(1, X1, X2))
	betavec <- c(-0.2, 0.1, 0.2)
	if(sigma2==0.50001){betavec <- c(-0.2, 0.0, 0.2)}

	Y <- c(X%*%betavec) + epsilon
	std_beta_hat <- betavec[-1]*(apply(X[,-1],2,sd)/sd(Y))
	true_std_beta <- std_beta_hat[1]
	resultsmat[jjj, "true_std_beta"] <- round(true_std_beta,3)
	Xmatrix_ <- X[1:N,]
}

if(nVar==4){
	epsilon <- rnorm(length(X1), 0, sqrt(sigma2))
	X  <- as.matrix(cbind(1, X1, X2, X3, X4))
	betavec <- c(0.2, 0.1, 0.14, -0.1, -0.1)
	if(sigma2== 0.50001){betavec <- c(0.2, 0.00, 0.14, -0.1, -0.1)}

	Y <- c(X%*%betavec) + epsilon
	std_beta_hat <- betavec[-1]*(apply(X[,-1],2,sd)/sd(Y))
	true_std_beta <- std_beta_hat[1]
	resultsmat[jjj, "true_std_beta"] <- round(true_std_beta, 3)
	Xmatrix_ <- X[1:N,]
}


pval_list[[jjj]] <- list()
Deltavec <- seq(0.01, 0.25, 0.005)

for(iii in 1:nSim){
	if(randomreg){Xmatrix_ <- X[sample(1:dim(X)[1], N, replace=TRUE),]}
	epsilon <- rnorm(N, 0, sqrt(sigma2))
	Y <- Xmatrix_%*%betavec + epsilon

	lmmod <- summary(lm(Y ~ Xmatrix_[,-1]))
	R2 <- lmmod$r.squared
	N <- length(Y); K <- dim(Xmatrix_[,-1])[2]
	beta_hat <- lmmod$coef[-1,][1,1]; SE_beta_hat <- lmmod$coef[-1,][1,2];
	
	std_beta_hat 	<- beta_hat*(apply(Xmatrix_[,-1],2,sd)/sd(Y))[1]
	R2XkXk 	<- unlist(lapply(c(1:K), function(k) {
				summary(lm(Xmatrix_[,-1][,k]~Xmatrix_[,-1][,-k]))$r.squared}))[1]
	SE_std_beta_FIX <- sqrt((1-R2)/((1-R2XkXk)*(N-K-1)))
	SE_std_beta_RDM <- DEL(X=Xmatrix_[,-1], y=Y)$SEs[1]



	pval_list[[jjj]][[iii]]<-vector() 

	for(kk in 1:length(Deltavec)){
		
		Delta <- Deltavec[kk]
		
		# Fixed regressors:		
		if(!randomreg){		
		pval <- equivstandardBeta(Y=Y, Xmatrix= Xmatrix_[,-1], DELTA= Delta, random=FALSE)$pval[1]
			}

		# Random regressors:	
		if(randomreg){
		pval <- equivstandardBeta(Y=Y, Xmatrix= Xmatrix_[,-1], DELTA= Delta, random=TRUE)$pval[1]
			}
		pval_list[[jjj]][[iii]][kk] <- pval
		
		}
	
	}

pvalmat <- as.data.frame(pval_list[[jjj]], col.names=NA, row.names=Deltavec)	

problist[[jjj]] <- rowMeans(pvalmat<resultsmat[jjj, "alpha_sig"])

print(resultsmat[jjj,])
}



problist_df <- as.data.frame(problist, col.names=c(1:dim(resultsmat)[1]))
colnames(problist_df)<-paste("jjj",c(1:dim(resultsmat)[1]), sep="")
problist_df$Delta<-Deltavec
library(tidyr)
problist_long<-gather(problist_df,JJJindex,pr_less_alpha,jjj1:jjj32, factor_key=TRUE)


resultsmatall<-merge(resultsmat,problist_long, by="JJJindex",all=TRUE)

resultsmatall$true_std_beta <- as.factor(resultsmatall$true_std_beta)
resultsmatall$N <- as.factor(resultsmatall$Nsample)

resultsmatall$g <- as.factor(paste(as.character(resultsmatall$true_std_beta),as.character(resultsmatall$N), sep="_"))

resultsmatall <- transform(resultsmatall,
  Nvar = factor(Nvar, levels = sort(unique(resultsmatall$Nvar)), c( (("K = 2")),(("K = 4")))  ))

resultsmatall$KK <- as.numeric(resultsmatall$Nvar)*2
resultsmatall$NN <- as.numeric(as.character(resultsmatall$N))

return(resultsmatall)}

########################################


simtest <- simstudy1(FALSE, 50)


resultsmatall<-simtest
####  plot with full axis:
library(ggplot2)
qplot(x=Delta, y= pr_less_alpha, group= g, pch=N, lty=N, col= true_std_beta, data= resultsmatall)+geom_line()+  geom_hline(yintercept = unique(resultsmat$alpha_sig))+ facet_grid(Nvar ~. ) + scale_x_continuous(breaks = seq(0, 0.2, by = 0.02))  + scale_y_continuous(breaks = seq(0, 1, by = 0.1))+ labs(x = expression(Delta),y = expression("probability of p" < alpha))+ coord_cartesian(ylim = c(0, 1)) 


####  plot with truncated axis:
library(ggplot2)
qplot(x=Delta, y= pr_less_alpha, group= g, pch=N, lty=N, col= true_std_beta, data= resultsmatall)+geom_line()+  geom_hline(yintercept = unique(resultsmat$alpha_sig))+ facet_grid(Nvar ~ . ) + scale_x_continuous(breaks = seq(0, 0.2, by = 0.02))  + scale_y_continuous(breaks = seq(-0.01, 1, by = 0.05))+  coord_cartesian(ylim = c(0, 0.20))+ labs(x = expression(Delta), y = expression("probability of p" < alpha))






