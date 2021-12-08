

PMSC = function(
  x_list,
  K,
  alpha =100,
  beta = 0.001,
  gamma = 0.001,
  .max_iter = 200,
  .tol = 1e-4
){
  ######## update function
  
  update_z = function(X,C,alpha,beta){
    
    H = as.matrix(dist(C))^2
    
    first_term = solve(t(X) %*% X + alpha * diag(ncol(X)))
    
    Z = sapply(1:ncol(X),
               function(x) {
                 first_term %*% (t(X) %*% X[,x] - beta / 4 * H[,x])
               })
    
    return(Z)
  }
  
  
  
  update_c = function(Z,YY,alpha,beta,gamma,weight,K) {
    
    Z = (Z + t(Z))/2
    
    D = diag(colSums(Z))
    
    L = D - Z
    
    M = beta * L - 2 * gamma * weight * YY + gamma * weight * diag(nrow(Z))
    
    C = Re(RSpectra::eigs(M,k = K,which = "SR")$vectors)
    
    return(C)
    
  }
  
  
  update_y = function(c_list, weight,K) {
    P = mapply(function(f, w) {
      w * (diag(nrow(f)) - 2 * f %*% t(f))
    }, c_list, weight,
    SIMPLIFY = F)
    
    P = Reduce("+", P)
    
    P = (P + t(P))/2
    
    Y = Re(RSpectra::eigs(P,k = K,which = "SR")$vectors)
    
    return(Y)
  }
  
  update_weight = function(YY,c_list){
    
    weight = sapply(c_list, function(f){
      min(1e+2,1/(2*norm(YY-f%*%t(f),"F")))
    })
    
    return(weight)
  }
  
  
  vec2mat = function(vec){
    K = unique(vec)
    M = lapply(K, function(x) as.numeric(vec == x))
    M = do.call("cbind",M)
    return(M)
  }
  
  
  ######## initialize 
  # x_list = list(iris[,-5],iris[,-5])
  x_list = lapply(x_list, as.matrix)
  x_list = lapply(x_list, scale)
  x_list = lapply(x_list, t)
  
  dim_list = lapply(x_list, dim)
  
  weight = 1/length(x_list)
  
  # c_list = lapply(x_list, function(x) vec2mat(kmeans(t(x),K)$cluster))
  c_list = lapply(dim_list, function(d) matrix(rbinom(d[2]*K,1,weight),d[2],K))
  
  Y = update_y(c_list,weight,K)
  
  
  ######## updating 
  iter = 1
  Diff = Inf
  
  while (iter<.max_iter & Diff > .tol
         ) {
    
    past_YY = YY = Y%*%t(Y)

    z_list = mapply(update_z, x_list,c_list,alpha,beta,
                    SIMPLIFY = F)
    
    c_list = mapply(update_c, z_list,list(YY),alpha,beta,gamma,weight,K,
                    SIMPLIFY = F)
    
    Y = update_y(c_list,weight,K)
    
    YY = Y%*%t(Y)
    
    # weight = update_weight(YY,c_list)
    
    Diff = abs(norm(past_YY - YY))
    
    iter = iter+1
  }
  
  if(Diff >.tol) message("Early Stop")
  
  Y = kmeans(Y,K)$cluster
  # Y = apply(Y, 1, which.max)
  
  result = list(
    Z = z_list,
    clusters = c_list,
    Y = Y,
    weight = weight,
    diff = Diff
  )
  
  return(result)
}



PMSC.grid_search = function(x_list,
                            K,
                            len = NA,
                            alpha = exp(seq(log(0.001),log(100),len = 5)),
                            beta = exp(seq(log(0.001),log(100),len = 5)),
                            gamma = exp(seq(log(0.001),log(100),len = 5)),
                            .max_iter = 100,
                            .tol = 1e-4,
                            .mc = getOption("mc.cores", 2L)
){
  
  if(!is.na(len)) alpha = beta = gamma = exp(seq(log(0.001),log(100),len = len))
  
  result = expand.grid(alpha = alpha, beta = beta, gamma = gamma)
  result$result = mcmapply(
    PMSC,
    list(x_list),
    K,
    result$alpha,
    result$beta,
    result$gamma,
    .max_iter,
    .tol,
    SIMPLIFY = F,
    mc.cores = .mc
  )
  
  result$cluster = lapply(result$result,function(x) x$Y)
  
  return(result)
}
