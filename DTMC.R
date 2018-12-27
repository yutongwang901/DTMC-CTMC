## R functions of discrete-time markov chain(DTMC)

## Enter a probability transition matrix
## commas separate entries in a row
## colons separate rows

P=matrix(c(.1,.9,
           .9,.1),nrow=2,ncol=2,byrow=TRUE)
P 

## Calculate the t-step transition probability matrix：P^(t)

install.packages("expm") 
library(expm)

P2=P%^%2
P2 ## P^(2)
P3=P%^%3
P3 ## P^(3)

## Calculate the expected occupancy times
## m_ij(t) = Sum from r=0 to r=t of P^(r)

## verify that P^0 = I
P%^%0

## start out M at time zeros
n=2 ## number of states 
t=5
Mt=diag(n)  
for(r in 1:t){
  Mt=Mt+P%^%r
}
Mt

## Find stationary distribution PI of markov chain when it exists
## PI = normalized eigenvector of P' corresoponding to the eigenvalue of 1
##

EE=eigen(t(P))
evects=EE$vect
evals=EE$val
## change "index" to be the column index of the "evals" matrix with diagonal = 1
index=which(evals==1)
PI=evects[,index]
## normalize
PI=PI/sum(PI)
PI
## verify by considering P^100
P%^%100


