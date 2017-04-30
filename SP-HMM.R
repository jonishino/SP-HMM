#---------- calc. heritability from beta:effect sizes(vector), p:derived allele frequencies(vector), K:prevalence-----------
get_herit <- function(beta,p,K,f00,itr){
 OR1<-exp(beta) 
 OR2<-exp(2*beta)
 f0 <- f00
 for(i in 1:itr){
  RR1 <- OR1/(1+f0*(OR1-1))
  RR2 <- OR2/(1+f0*(OR2-1))
  f0 <- K/((1-p)^2+2*p*(1-p)*RR1+p^2*RR2)
 }
 f1 <- f0*RR1
 f2 <- f0*RR2
 T <- qnorm(1-f0) 
 myu0 <- 0 
 myu1 <- T-qnorm(1-f1) 
 myu2 <- T-qnorm(1-f2)
 myuall <- myu1*2*p*(1-p) + myu2*p^2
 v_exp <- (1-p)^2*(myu0-myuall)^2 + 2*p*(1-p)*(myu1-myuall)^2 + p^2*(myu2-myuall)^2
 he <- v_exp/(1+v_exp)
 return(he)
}

get_herit2 <- function(RR1,p,K){ # approximation version
 T <- qnorm(1-K)
 nyu <- dnorm(T)/K
 he <- 2*p*(1-p)*(RR1-1)^2/(nyu^2)
 return(he)
}






hmm_em <- function(y, v_y, init_pi, init_p_b, a, emitr,prev=NA,p=NA,brk=1,fixedpi=0){
  pi <- init_pi;pi0 <- 1-pi
  B <- length(a)
  if(length(init_p_b)==1) p_b <- rep(1,B)/B
  if(length(init_p_b)>1) p_b <- init_p_b
  dnorm <- matrix(0,nrow=B,ncol=length(y))
  for(b in 1:B) dnorm[b,] <- dnorm(y, mean=a[b], sd=sqrt(v_y))
  mlike0 <- dnorm(y, mean=0, sd=sqrt(v_y))
  pre_pi=0
  for(itr in 1:emitr){
    mlike1 <- t(as.matrix(p_b)) %*% dnorm
    tau1 <- 1- pi0 * mlike0 / (pi0 * mlike0  + pi * mlike1)
    if(fixedpi==1 & itr <= 200){
    } else{ pi <- mean(tau1) }
    pi0 <- 1-pi
    p_b <- (dnorm%*%t(tau1/mlike1))*(p_b/sum(tau1))
#    plot(p_b~a,type="h",ylim=c(0,0.6),main=paste(itr,"th iteration"))
    plot(p_b~a,type="h",main=paste(itr,"th iteration"))
    cat(itr, pi, "\n")
    if(brk==1 & itr > 200 & abs((pre_pi-pi)/pre_pi) < 0.00005) break
    pre_pi <- pi
  }
  pm <- t(t(as.matrix(p_b*a)) %*% dnorm /mlike1*tau1) #calc. posterior mean

  # SNP heritability from estimated gentic structure
  if(is.na(prev)!=1 & (is.na(p))[1]!=1){
    h<-0
    for( b in 1:B ) h <- h+(p_b[b]*pi)*get_herit(beta=a[b],p=p,K=prev,f00=prev,itr=10)
    return(list(pi=pi, p_b=p_b,pm=pm, a=a,herit=sum(h),loglike=sum(log(pi0*mlike0  + pi*mlike1))))
  }
  else {
    return(list(pi=pi, p_b=p_b,pm=pm, a=a,loglike=sum(log(pi0*mlike0  + pi*mlike1))))
  }

}

# automatically determined initial values (pi & g)
hmm_em_aut_pre<- function(y, v_y,init_pi=seq(0.1,0.9,by=0.1),emitr00=200,emitr01=200,emitr=2000,
                        a0=c(seq(-0.3,-0.005,by=0.005),seq(0.005,0.3,by=0.005)),
                        a=c(seq(-0.3,-0.005,by=0.005),seq(0.005,0.3,by=0.005)),prev=NA,p=NA,brk=1,nbs=0){
#-----1. Determine initial pi by selecting pi giving maximum likelihood
  loglikes <- pi <- c() 
  for (i in 1:length(init_pi)){
    hmm_out1<-hmm_em(y=y,v_y=v_y,init_pi=init_pi[i],init_p_b=0,a=a0,emitr=emitr00,prev=prev,p=p,brk=0)
    loglikes[i]<-hmm_out1$loglike
    pi[i]<-hmm_out1$pi
    print(paste("loglike:",loglikes[i]));print(paste("pi:",pi[i]))
  }
#-----2. Determine initial g using initial pi (pi is fixed in the process)
  hmm_out2<-hmm_em(y=y,v_y=v_y,init_pi=pi[which.max(loglikes)],init_p_b=0,a=a,emitr=emitr01,prev=prev,p=p,brk=0,fixedpi=1)
#-----3. Run EM algorithm using initial pi & g determined above
  hmm_out3<-hmm_em(y=y,v_y=v_y,init_pi=pi[which.max(loglikes)],init_p_b=hmm_out2$p_b,a=a,emitr=emitr,prev=prev,p=p,brk=brk)
  return(hmm_out3)
}






hmm_em_aut<- function(y, v_y,init_pi=seq(0.1,0.9,by=0.1),emitr00=200,emitr01=200,emitr=2000,
                        a0=c(seq(-0.3,-0.005,by=0.005),seq(0.005,0.3,by=0.005)),
                        a=c(seq(-0.3,-0.005,by=0.005),seq(0.005,0.3,by=0.005)),prev=NA,p=NA,brk=1,nbs=0){
   hmm_out <- hmm_em_aut_pre(y, v_y,init_pi=init_pi,emitr00=emitr00,emitr01=emitr01,emitr=emitr,
                        a0=a0,a=a,prev=prev,p=p,brk=brk,nbs=nbs)
  #----- Bootstrap
  if(nbs>0){
    library(doSNOW)
    y_bs<-matrix(0,nrow=length(y),ncol=nbs)
    for(i in 1:nbs) y_bs[,i] <- rnorm(length(y), mean=sample(c(a,0),length(y),
     replace=TRUE,prob=c(hmm_out$p_b*hmm_out$pi,(1-hmm_out$pi))), sd=sqrt(v_y))
     registerDoSNOW(makeCluster(4, type = "SOCK"))
     aaa<- foreach(i = 1:nbs,.export=c("hmm_em_aut_pre","hmm_em","get_herit")) %dopar% hmm_em_aut_pre(y=y_bs[,i], v_y=v_y,init_pi=init_pi,emitr00=emitr00,emitr01=emitr01,emitr=emitr,
                        a0=a0,a=a,prev=prev,p=p,brk=brk,nbs=0)

     bs_pi <- bs_herit <- c()
    if(is.na(prev)==1){      
      for(i in 1:nbs) bs_pi <- append(bs_pi,aaa[[i]]$pi)
      hmm_out <- append(hmm_out,list(bs_pi))    
      names(hmm_out)[length(hmm_out)] <- "bs_pi"
    }
    if(is.na(prev)!=1){      
      for(i in 1:nbs) bs_pi <- append(bs_pi,aaa[[i]]$pi)
      for(i in 1:nbs) bs_herit <- append(bs_herit,aaa[[i]]$herit)
      hmm_out <- append(hmm_out,list(bs_pi))
      hmm_out <- append(hmm_out,list(bs_herit))
      names(hmm_out)[(length(hmm_out)-1)] <- "bs_pi"
      names(hmm_out)[length(hmm_out)] <- "bs_herit"
    }
  }
  return(hmm_out)
}






