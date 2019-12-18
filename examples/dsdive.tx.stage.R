data('dive.sim')
attach(dive.sim)
attach(params)

p = dsdive.tx.stage(t0 = sim$times[1], sub.tx = sub.tx, 
                    surf.tx = surf.tx, t0.dive = 0, t.stage2 = NA, 
                    rates = params$lambda[1], t.scale = 60, 
                    model = 'conditional', s0 = 1)

detach(params)
detach(dive.sim)
