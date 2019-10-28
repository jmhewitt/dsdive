beta = matrix(c(1,   0, -1, 
                0, -.5,  0), 
              nrow = 2, byrow = TRUE)

lambda = c(3, 2, 1)

sub.tx = c(60, .05)

surf.tx = .02

dsdive.tx.params(t0 = 0, num.depths = 200, d0 = 0, d0.last = NULL, s0 = 1, 
                 beta = beta, lambda = lambda, sub.tx = sub.tx, 
                 surf.tx = surf.tx)

dsdive.tx.params(t0 = 0, num.depths = 200, d0 = 60, d0.last = NULL, s0 = 2, 
                 beta = beta, lambda = lambda, sub.tx = sub.tx, 
                 surf.tx = surf.tx)

dsdive.tx.params(t0 = 0, num.depths = 200, d0 = 60, d0.last = 61, s0 = 2, 
                 beta = beta, lambda = lambda, sub.tx = sub.tx, 
                 surf.tx = surf.tx)
