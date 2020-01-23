data('dive.sim')
attach(dive.sim)

n.impute = 10
tstep = 1

imputed.midpoint = crawl.impute(depth.bins = depth.bins, 
                                depths = sim.obs$depths, 
                                times = sim.obs$times, N = n.impute, 
                                depths.impute = 'midpoint', tstep = tstep)

imputed.unif = crawl.impute(depth.bins = depth.bins, 
                                depths = sim.obs$depths, 
                                times = sim.obs$times, N = n.impute, 
                                depths.impute = 'uniform', tstep = tstep)

imputed.tnorm = crawl.impute(depth.bins = depth.bins, 
                            depths = sim.obs$depths, 
                            times = sim.obs$times, N = n.impute, 
                            depths.impute = 'truncnorm', tstep = tstep, 
                            sd.scale = .125)

pl.m = plot(x = sim.obs, depth.bins = depth.bins, errorbars = TRUE, 
            imputed.list = imputed.midpoint)

pl.u = plot(x = sim.obs, depth.bins = depth.bins, errorbars = TRUE, 
            imputed.list = imputed.unif)

pl.n = plot(x = sim.obs, depth.bins = depth.bins, errorbars = TRUE, 
            imputed.list = imputed.tnorm)

ggpubr::ggarrange(pl.u, pl.n, pl.m, 
                  labels = c('Uniform', 'Trunc. norm.', 'Midpoint'), 
                  ncol = 1, nrow = 3)

detach(dive.sim)
