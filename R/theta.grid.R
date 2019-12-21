#' Build parameter grid in the proper enumeration order
#'  
#' @export
#' 
theta.grid = function(pi1, pi2, lambda1, lambda2, lambda3) {
  
  pi1.len = length(pi1)
  pi2.len = length(pi2)
  lambda1.len = length(lambda1)
  lambda2.len = length(lambda2)
  lambda3.len = length(lambda3)
  
  res = matrix(nrow = prod(pi1.len, pi2.len, lambda1.len, lambda2.len, 
                           lambda3.len), 
               ncol = 5)
  
  colnames(res) = c('pi1', 'pi2', 'lambda1', 'lambda2', 'lambda3')
  
  ind = 1
  for(i in 1:pi1.len) {
    li = pi1[i]
    for(j in 1:pi2.len) {
      lj = pi2[j]
      for(k in 1:lambda1.len) {
        lk = lambda1[k]
        for(l in 1:lambda2.len) {
          ll = lambda2[l]
          for(m in 1:lambda3.len) {
            res[ind,] = c(li, lj, lk, ll, lambda3[m])
            ind = ind + 1
    }}}}}
  
  res
}