# Factorize an number using Sieve of Eratosthenes -method
divisors <- function( m ){
  d_list <- c()
  for(d in 1:(m/2)){
    if(m %% d == 0){
      d_list <- c(d_list,d)
    }
  }
  if(m > 1) d_list <- c(d_list,m)
  return(d_list)
}

t_recursion <- function( q , t , K = 1 ){
  t_q <- 0
  for(m in 1:(q-1)){
    t_q <- t_q + t[q-m]*sum(sapply(X = divisors(m),FUN = function(d) d*K^(1-(1==d))*t[d]))
  }
  t_q <- t_q/(q-1)
  return(c(t,t_q))
}

computeT <- function( q_max , K = 1 ){
  t <- 1 # base case t_q(1) == 1
  for(i in 2:max(3,q_max)){
    t <- t_recursion(i,t,K)
  }
  return(t[1:q_max])
}

countNonPlantedTrees <- function( q_max, K = 1 ){
  t <- computeT( q_max+1 , K )
  t_diff <- t[2:(q_max+1)] - t[1:q_max]*K^(1:q_max > 1)
  return(t_diff[2:q_max])
}

generateTable1 <- function( q_max = 15, K = 1 ){
  require(package = "xtable") # package for outputting latex source code
  t <- computeT(q_max,K)
  t_planted <- c(1, t[1:(q_max-1)]*K^(1:(q_max - 1) > 1))
  dt <- data.frame("t" = t,"t_planted" = t_planted, "t_non-planted" = t - t_planted)
  xt <- xtable(dt)
  print.xtable(x = xt)
}