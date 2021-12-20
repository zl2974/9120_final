###### Normalize matrix

normalize = function(.matrix, type = 1) {
    if (type == 1) {
        D = colSums(.matrix)
        D[D == 0] = .Machine$double.eps
        di = diag(1 / D)
        NL = .matrix %*% di
    }
    if (type == 2) {
        D = rowSums(.matrix)
        D[D == 0] = .Machine$double.eps
        di = diag(1 / D)
        NL = di %*% .matrix
    }
    if (type == 3){
      D = rowSums(.matrix) + .Machine$double.eps
      di = diag(sqrt((1 / D)))
      NL = di %*% .matrix %*% di
                
    }
    return(NL)
}
