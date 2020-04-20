# Packages
  library(TruncatedDistributions) #TruncatedDistributions: Use Truncated Distribution functions
  library(readr)
  library(matlib)
  library(pracma)
# Defining Y, A, initial X and lambda   
  Y = c(2,28,9,5,34,16,35)
  A = as.matrix(read_table2("Documentos/Seminarios/matriz A",col_names = FALSE))
# Using QR decomposition
  pivot = qr(A)$pivot
  A1 = A[,pivot[1:length(Y)]]
  A2 = A[,-pivot[1:length(Y)]]
  s.A1 = solve(A1)
# Valores iniciais para lambda e X
  X2 = c(0,0,0,0,0)
  X1 = s.A1%*%(Y-A2%*%X2);X1

  
    