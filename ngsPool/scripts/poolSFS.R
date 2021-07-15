
calcSFS <- function(n, daf, norm=TRUE) {
        # from relative frequencies [0,1], both population and sample
        # including non-snps
        x <- seq(0,1,1/(2*n))
        sfs <- rep(0, length(x))
        # daf is [0,1]
        for (i in daf) {
                ind <- which.min(abs(i-x))
                sfs[ind] <- sfs[ind] + 1
        }
        if (norm) sfs <- sfs/sum(sfs)
        sfs
}

getExpSFS <- function(K,S=0) {
        ee <- (1/(1:((S*2)-1))^(1/K))
        ee/sum(ee)
}

getFit <- function(vec,Q=NA,th=1e-8) {
        P <- vec + th; #P <- P/sum(P) + th
        Q <- Q + th
        KL <- sum(P*log(P/Q), na.rm=T)
        KL2 <-sum(Q*log(Q/P), na.rm=T)
        MSE <- sum((Q-P)^2,na.rm=T)/length(P)
        RMSE <- sum(sqrt((Q-P)^2),na.rm=T)/length(P)
        res <- c(KL, KL2, MSE, RMSE)
        res
}

calcFit <- function(par, data, metric="KL", func=getExpSFS, th=1e-8) {

        if(is.null(dim(data))) data <- t(data)
        n <- (ncol(data)+1)/2

        Q <- func(S=n, K=par[1])

        metrics <- t(apply(FUN=getFit, X=data, Q=Q, MAR=1, th=th))

        res <- apply(FUN=sum, X=metrics, MAR=2, na.rm=T)

        if (metric=="KL") return(res[1]); if (metric=="KL2") return(res[2])
        if (metric=="MSE") return(res[3]); if (metric=="RMSE") return(res[4])
}

norm <- function(a,th=0.05) {
        a[which(a<th)] <- 0
        a/sum(a)
}


args <- commandArgs(T)

fin_saf <- args[1]
norm_th <- as.numeric(args[2]) # 0.05
offset <- as.numeric(args[3]) # 1e-6

if (length(args)<3) {
	if (length(args)==1) {
		norm_th <- 0.5
		offset <- 1e-6
	}
	if (length(args)==2) offset <- 1e-6
}

# read saf file
saf <- exp(read.table(fin_saf))

# sample size for diploids
S <- (ncol(saf)-1)/2

# only SNPs (it should have been done a SNP calling)
daf_mle <- (apply(X=saf, MAR=1, FUN=which.max)-1) / (S*2)
ind <- which(daf_mle>0 & daf_mle<1)
nsnps_sample <- length(ind)


# normalise and set likelihoods to 0 if below norm_th
saf_th <- t(apply(X=saf[ind,2:(S*2)], FUN=norm, MAR=1, th=norm_th))

# 1) SFS from MLE(freq) (max saf) 
sfs_mle <- calcSFS(S, daf_mle[ind], F)
cat("count",sfs_mle,"\n")

# 2) fit SFS from MLE(freq) (max saf)
sfs_mle_snps <- sfs_mle[2:(S*2)]
sfs_mle_snps <- sfs_mle_snps / sum(sfs_mle_snps)
opt_mle <- optim(par=c(1), fn=calcFit, data=sfs_mle_snps, lower=c(0), upper=c(1e2), method="L-BFGS-B", metric="KL", th=offset)
sfs_opt_mle <- c(NA, getExpSFS(S=S, K=opt_mle$par[1])*nsnps_sample, NA)
cat("fit_count",sfs_opt_mle,"\n")

# 3) fit SFS from SAF
opt_th <- optim(par=opt_mle$par, fn=calcFit, data=saf_th, lower=c(0), upper=c(1e2), method="L-BFGS-B", metric="KL", th=offset)
sfs_opt_th <- c(NA, getExpSFS(S=S, K=opt_th$par[1])*nsnps_sample, NA)
cat("fit_saf",sfs_opt_th,"\n")






