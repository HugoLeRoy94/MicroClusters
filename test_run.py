import numpy as np
import MC

E1,E2 = 1.,0.
Interactions = [[0.,0.,0.],[0.,E1,0.],[0.,0.,E2]]

step_tot = 50*10**6
check_step = 100

mc = MC.MC(100,100,0,0,interactions=Interactions,temperature=1.)
av_c_size = np.zeros(check_step,dtype=float)
for n_steps in range(check_step):
    mc.monte_carlo_steps(step_tot//check_step)
    av_c_size[n_steps] = np.mean(mc.box.cluster_size())