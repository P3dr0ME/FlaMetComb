#%%
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science'])
import time
start_time = time.time()


gas_mix_incomplete = ct.Solution('gri30.yaml')

#%% Cálculo de velocidad de llama
T_0 = 300 # K

phi = 0.8
p = 10*ct.one_atm


print(f"\033[1;36m### EQUIVALENCE RATIO: {phi} ###\033[0m")
gas_mix_incomplete.set_equivalence_ratio(phi, 'CH4', 'O2: 1.0, N2: 3.76')


print(f"\033[1;36m## PRESSURE (Pa): {p} ##\033[0m")
gas_mix_incomplete.TP = T_0, p


# Llama
flame = ct.FreeFlame(gas=gas_mix_incomplete, width=0.001)
flame.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)


flame.solve(loglevel=0, refine_grid=True, auto=True)

print(f"\033[1;36m# LAMINAR BURNING SPEED (cm/s): {flame.velocity[0]*100} cm/s ###\033[0m")
