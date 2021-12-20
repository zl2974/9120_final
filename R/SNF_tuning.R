
SNF_tuning = function(data_list, K) {
  library(tidyverse)
  neighbour = seq(10, 50, 5)
  sigma = seq(0.3, 0.8, 0.05)
  data_list = lapply(data_list, function(x)
    as.matrix(dist(x)))
  
  tuning = expand.grid(K = neighbour,
                       sigma = sigma) %>%
    
    mutate(affinity = map2(
      K,
      sigma,
      ~ lapply(data_list, function(x)
        SNFtool::affinityMatrix(x, .x, .y))
    ))
  
  tuning$SNF = lapply(tuning$affinity,SNFtool::SNF)
  tuning$cluster = mapply(SNFtool::spectralClustering,tuning$SNF,K,SIMPLIFY = F)
  
  return(tuning[,c("K","sigma","cluster")])
}
