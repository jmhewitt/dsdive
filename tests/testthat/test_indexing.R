context("3d array indexing routines")

test_that("(x,y,z) coords to flattened indices", {
  # set up dimensions
  num.bins = 7
  num.stages = 3
  
  # build test arrays
  tx.raw = array(1:(num.bins^2 * num.stages), 
                 dim = c(num.bins, num.bins, num.stages))
  tx.vec = as.numeric(tx.raw)
  
  # check indexing
  for(i in 1:num.bins) {
    for(j in 1:num.bins) {
      for(k in 1:num.stages)
        expect_equal(tx.raw[i,j,k], 
                     tx.vec[toInd(x = i, y = j, z = k, 
                                  x.max = num.bins, y.max = num.bins)])
    }
  }
})

test_that("flattened indices to (x,y,z) coords", {
  # set up dimensions
  num.bins = 60
  num.stages = 3
  
  # build test arrays
  tx.raw = array(1:(num.bins^2 * num.stages), 
                 dim = c(num.bins, num.bins, num.stages))
  tx.vec = as.numeric(tx.raw)
  
  # check indexing
  for(i in 1:num.bins) {
    for(j in 1:num.bins) {
      for(k in 1:num.stages)
        expect_equal(c(i,j,k),
                     fromInd(ind = tx.raw[i,j,k], 
                             x.max = num.bins, y.max = num.bins))
      }
  }
})