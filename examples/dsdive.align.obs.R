data('dive.sim')
attach(dive.sim)
attach(params)

t.stages = sim$times[c(FALSE, diff(sim$stages)==1)]
  
dsdive.align.obs(depths = sim.obs$depths, times = sim.obs$times, 
                 t.stages = t.stages, offset = -300, 
                 offset.tf = -300)

detach(params)
detach(dive.sim)
