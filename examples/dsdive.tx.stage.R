data('dive.sim')
attach(dive.sim)
attach(params)

p = dsdive.tx.stage(t0 = sim$times, d0 = sim$depths, sub.tx = sub.tx, 
                    surf.tx = surf.tx, t0.dive = 0)

detach(params)
detach(dive.sim)