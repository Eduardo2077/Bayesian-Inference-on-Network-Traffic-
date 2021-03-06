# Packages
  library(TruncatedDistributions) #TruncatedDistributions: Use Truncated Distribution functions
  library(readr)
  library(matlib)
  library(pracma)
  library(dplyr)
  library(extraDistr)
  library(networkTomography)
# Defining Y, A, initial X and lambda   
  Y = c(2,28,9,5,34,16,35)
  A = as.matrix(read_table2("C:\\Users\\eduar\\Downloads\\Bayesian-Inference-on-Network-Traffic--master\\matriz A",col_names = FALSE))
# Using QR decomposition
  pivot = qr(A)$pivot
  A1 = A[,pivot[1:length(Y)]]
  A2 = A[,-pivot[1:length(Y)]]
  s.A1 = solve(A1)
# FunÃ§Ã£o para gerar valores iniciais viÃ¡veis para X
  X0 = function(Y,A){
    A1 = A[,pivot[1:length(Y)]]
    A2 = A[,-pivot[1:length(Y)]]
    s.A1 = solve(A1)
    X1 = X2 = vector()
    lambda = r = 0
    while(r!=TRUE){
      lambda = runif(1,1,20)
      X2 = rpois(5,lambda)
      X1 = s.A1%*%(Y-A2%*%X2)
      r = all(X1>=0)
  }
    return(list(X1 = X1,X2 = as.matrix(X2,ncol = 1)))
  }
  
#Solução possível: X1 = [2,0,3,5,5,1,16] e X2 = [6,7,7,10,9]
#                  X1 = [2,0,0,4,5,5,10] e X2 = [9,9,7,10,15]
  
# Gibbs and Metropolis-Hastings Algorithms
  
# Gerando os lambdas utilizando a posteiori condicional a X
  x0 = X0(Y,A)
  lambda = list(X1 = mapply(rtgamma, n = 1, shape = x0[[1]]+1,scale = 1, b =100),
                X2 = mapply(rtgamma, n = 1, shape = x0[[2]]+1,scale = 1, b =100))
  lambda$X1 = as.matrix(lambda$X1,ncol = 1)
  lambda$X2 = as.matrix(lambda$X2,ncol = 1)
# Definindo a posteriori de Xi dado X2,-i, lambda e Y
  p = function(x,lambdai,X1,lambda1){
    prdt = prod(lambda1^X1/factorial(X1))
    res = lambdai^x/factorial(x)*prdt
    return(res)
  }
# Gerando posteiori condicional para Xi dado X2,-i utilizando Metropolis Hasting within Gibbs
# 1 - Gerar Xi* ~ P(lambda_i)
# 2 - Aceitar com probabilidade p.
  
  lista.teste = list(X6 = c(1,9),X8 = c(1,33),X9 = c(7,15),X11 = c(1,27),X12 = c(6,34))
  for(j in 1:6000){
    x1.r = x0$X1[,j]
    x2.r = x0$X2[,j]
  for(i in 1:5){
    r = FALSE
    while(r!=TRUE){
      q = rpois(1,lambda$X2[i,j])
      v.teste = x2.r
      v.teste[i] = q
      X1.teste = s.A1%*%(Y-A2%*%v.teste)
      r = all(X1.teste>=0)
    }
      u = runif(1,0,1)
      rho = min(1,dpois(x2.r[i],lambda$X2[i,j])/dpois(q,lambda$X2[i,j])*p(q,lambda$X2[i,j],X1.teste,lambda$X1[,j])/p(x2.r[i],lambda$X2[i,j],x1.r,lambda$X1[,j]))
      if(u<rho){
        x2.r[i] = q
      }else{x2.r[i] = x2.r[i]}
      x1.r = s.A1%*%(Y-A2%*%x2.r)
    }
  
  x0$X1 = cbind(x0$X1,x1.r)
  x0$X2 = cbind(x0$X2,x2.r)
  
  lambda$X1 = cbind(lambda$X1,mapply(rtgamma, n = 1, shape = x0$X1[,j+1]+1,scale = 1, b =100))
  lambda$X2 = cbind(lambda$X2,mapply(rtgamma, n = 1, shape = x0$X2[,j+1]+1,scale = 1, b =100))
}

  par(mfrow = c(4,3))  
  for(i in 1:7){
  barplot(table(x0$X1[i,-c(1000)]))
  }
  for(i in 1:5){
    barplot(table(x0$X2[i,-c(1000)]))
  }
  
  Xverdadeiro = c(2,2,0,8,5,6,7,7,4,9,17,11) 
    
  par(mfrow = c(4,3))  
  for(i in 1:7){
    plot(density(lambda$X1[i,-c(10000)]),ylim = c(0,.2))
    curve(dgamma(x,shape =  Xverdadeiro[i]+1,scale = 1),add=TRUE)
  }
  for(i in 1:5){
    plot(density(lambda$X2[i,-c(10000)]),ylim = c(0,.2))  
    curve(dgamma(x,shape =  Xverdadeiro[i+7]+1,scale = 1),add=TRUE)
  }
  

