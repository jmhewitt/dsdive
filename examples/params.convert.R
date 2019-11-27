data('dive.sim')
attach(dive.sim)
attach(params)

spec = list(beta.sd = 1, beta.absmax = 10, lambda.sd = 1, sub.tx.mean = 0, 
            sub.tx.sd = 1, surf.tx.mean = 0, surf.tx.sd = 1)

v = params.toVec(par = params, spec = spec)

L = params.toList(par = v)

detach(params)
detach(dive.sim)
