data('dive.sim')
attach(dive.sim)
attach(dive.sim$params)

# uniformized transition rate
r.unif = max(outer(lambda, 2 * depth.bins$halfwidth, '/'))

# single-step transition matrix
m = dsdive.tx.matrix.uniformized(depth.bins = depth.bins, beta = beta, 
                                 lambda = lambda, s0 = 2, 
                                 rate.uniformized = r.unif)

# number of states
k = nrow(m)

# number of transitions
N = 10

# list of transition matrices
B = lapply(1:N, function(i) m)

# starting and ending coordinates
x0 = 5
xN = 10

# encode likelihood information
L = matrix(0, nrow = k, ncol = N+1)
L[x0,1] = 1                           # fixed starting location
L[xN,N+1] = 1                         # fixed ending location
L[,-c(1,N+1)] = 1/k                   # free transitions
L[x0,2] = 0                           # force a transition by 2nd step
L = sweep(L, 2, colSums(L), '/')      # restandardize likelihood

# forward filter
a = ff(B = B, L = L)

# backwards sample
set.seed(2019)
y = bs(a = a, B = B, L = L)

# execute forward filtering and backwards sampling via wrapper
set.seed(2019)
x = ffbs(B = B, L = L)

# sampling steps are identical
identical(x,y)

detach(dive.sim$params)
detach(dive.sim)
