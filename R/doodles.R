# Factorize an number using Sieve of Eratosthenes -method
divisors <- function(m){
  d_list <- c()
  for(d in 1:(m/2)){
    if(m %% d == 0){
      d_list <- c(d_list,d)
    }
  }
  if(m > 1) d_list <- c(d_list,m)
  return(d_list)
}

t_recursion <- function(q,t_lessThanQ,K = 1){
  t_q <- 0
  for(m in 1:(q-1)){
    t_q <- t_q + t_lessThanQ[q-m]*sum(sapply(X = divisors(m),FUN = function(d) d*K^(1-(1==d))*t_lessThanQ[d]))
  }
  t_q <- t_q/(q-1)
  return(c(t_lessThanQ,t_q))
}


computeT <- function(q_max,K=1){
  t <- 1 # base case t_q(1) == 1
  for(i in 2:max(3,q_max)){
    t <- t_recursion(i,t,K)
  }
  return(t[1:q_max])
}

# generateTable1 <- function(q_max,K=1){
#   t <- c(1,1)
#   for(i in 3:max(3,q_max+1)){
#     t <- t_recursion(i,t,K)
#   }
# #  print(t)
#   t[2:(q_max+1)] - t[1:q_max]
# }

generateTable1 <- function(q_max,K=1){
  t <- computeT( q_max+1 , K )
#  t <- c(1,1)
#   for(i in 3:max(3,q_max+1)){
#     t <- t_recursion(i,t,K)
#   }
#  #  print(t)
  t_diff <- t[2:(q_max+1)] - t[1:q_max]*K^(1:q_max > 1)
  return(t_diff[2:q_max])
}