## R functions of continuous-time markov chain(CTMC)

## Enter a transition rate matrix

R=matrix(c(0,3,0,0,
           2,0,3,0,
           0,4,0,3,
           0,0,4,0),nrow=4,byrow=T)
R

## (3) Calculate the  transition probability matrix: P(t)

Pctmc <- function(R,t){
  N=dim(R)[1]
  rvals=apply(R,1,sum)
  r=max(rvals)
  Phat=R/r
  diag(Phat) <- 1-rvals/r
  rt=r*t
  K=ceiling(max(c(rt+5*sqrt(rt),20)))
  P=diag(N)*dpois(0,rt)
  Pk=diag(N)
  for(k in 1:K){
    Pk=Pk%*%Phat
    P=P+dpois(k,rt)*Pk
  }
  P
}

## Calculate the expected occupancy times
## m <- ij(t) = Sum from r=0 to r=t of P^(r)

Mctmc <- function(R,t){
  N=dim(R)[1]
  rvals=apply(R,1,sum)
  r=max(rvals)
  Phat=R/r
  diag(Phat) <- 1-rvals/r
  rt=r*t
  K=ceiling(max(c(rt+5*sqrt(rt),20)))
  M=diag(N)*(1-ppois(0,rt))
  Pk=diag(N)
  for(k in 1:K){
    Pk=Pk%*%Phat
    M=M+(1-ppois(k,rt))*Pk
  }
  M/r
}


M3=Mctmc(R,3)
M3
apply(M3,1,sum) 

## Find stationary distribution PI of markov chain when it exists

Pctmc(R,100000) ## set t to infinity

## or the complex way

R
rvals=apply(R,1,sum)
Q=R-diag(rvals)
E=eigen(t(Q))
v.n=Re(E$vectors[,ncol(E$vectors)])
pi=v.n/sum(v.n)
pi

