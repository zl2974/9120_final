
DiMSC = function(x_list,
                 K,
                 alpha = 0.01,
                 beta  = 40,
                 .max_iter = 10,
                 .tol = 1e-6) {
  
  ##### update function
  .update_z = function(X,g_list,L,alpha,beta){
    A = C = t(X) %*% X
    G = Reduce("+",g_list)
    B = alpha * L + beta * G
    Z = maotai::sylvester(A,B,C)
    
    return(Z)
  }
  
  HSIC = function(Z) {
    
    G = t(Z) %*% Z
    H = diag(nrow(Z)) - matrix(1 / nrow(Z), nrow(Z), nrow(Z))
    G  = H %*% G %*% H
    # G = (nrow(Z)-1)^(-2) * G
    return(G)
    
  }
  
  w_function = function(x) SNFtool::affinityMatrix(
    as.matrix(
      dist(t(x))
    ),
    K = ncol(x)/K,
    sigma = 0.6
  )
  
  .update_l = function(z_list,K){
    w = do.call("+",z_list)
    w = normalize(w)
    w = eigen(w)$vec[,1:K]
    w = w%*%t(w)
    diag(colSums(w)) - w
  }
  
  
  ##### initialize
  x_list = lapply(x_list, as.matrix)
  x_list = lapply(x_list, scale)
  x_list = lapply(x_list, t)
  
  dim_list = lapply(x_list, dim)
  
  g_list = lapply(dim_list, function(d) matrix(0,d[2],d[2]))
  
  z_list = lapply(x_list, w_function)
  z_list[[1]] = cor(x_list[[1]])
  l_list = .update_l(z_list,K)
  
  ##### updating
  iter = 1
  Diff = Inf
  while (iter < .max_iter & Diff > .tol) {
    
    past_Z = z_list[[1]]
    
    g_list = lapply(z_list, HSIC)
    
    for(i in 1:length(x_list)){
      z_list[[i]] = .update_z(x_list[[i]], g_list[-i], l_list, alpha, beta)
      g_list[[i]] = HSIC(z_list[[i]])
      l_list = .update_l(z_list,K)
      message(i)
    }
    
    Diff = norm(past_Z - z_list[[1]],"2")
    message(Diff)
    iter = iter + 1
  }
  
  A = Reduce("+",z_list)
  A = (A +t(A))/2
  Y = RSpectra::eigs(A,K,which = "LR")$vec
  Y = kmeans(Y,3)$cluster
  
  result = list(Z = z_list,
                diff = Diff,
                A = A,
                cluster = Y)
  
  return(result)
}



DiMSC.grid_search =
  function(x_list,
           K,
           len = NA,
           alpha = exp(seq(log(0.001), log(100), len = 5)),
           beta = exp(seq(log(0.001), log(100), len = 5)),
           w_function = NA,
           .max_iter = 100,
           .tol = 1e-6,
           .mc = getOption("mc.cores", 2L)) {
    
    if (!is.na(len))
      alpha = beta  = exp(seq(log(1e-4), log(1e+2), len = len))
    
    result = expand.grid(alpha = alpha,
                         beta = beta)
    
    result$result =
      mcmapply(
        DiMSC,
        list(x_list),
        K,
        result$alpha,
        result$beta,
        w_function,
        .max_iter,
        .tol,
        mc.cores = .mc,
        SIMPLIFY = F
      )
    
    result$cluster = lapply(result$result,function(x) x$cluster)
    
    return(result)
  }
           
