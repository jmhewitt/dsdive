devtools::document()

# construct transition matrix
# B = matrix(runif(100), nrow = 10, ncol = 10)
B = matrix(0, nrow = 30000, ncol = 30000)
diag(B) = 4
diag(B[,-1]) = 1
diag(B[-1,]) = 1
B = B / rowSums(B)

library(Matrix)
library(spam)

B = as.dgCMatrix.spam(as.spam(B))


# set initial conditions
init = 1

# set number of transitions to make
t = 20

# forward sample transitions
s = numeric(t)
s[1] = init
for(i in 2:t) {
  s[i] = sample(x = 1:nrow(B), size = 1, prob = B[s[i-1],])
}

# observe process at some indices
na.inds = sample(x = 2:(t-1), size = .5*t)
O = s
O[na.inds] = NA

# wrap transition matrices
Bt = vector('list', t)
for(i in 1:t) {
  Bt[[i]] = B
}

devtools::document()
# draw from posterior
x = t(replicate(1e3, ffbs(Bt = Bt, O = O)))

# verify sampling reproduces observations
all(t(x[,!is.na(O)]) == O[!is.na(O)])


#
# check some posterior distributions of missing values
#

ind = which(is.na(O))[9]
table(x[,ind])/sum(table(x[,ind]))
s[ind]
