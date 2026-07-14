#%%
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science'])
import time
start_time = time.time()


gas = ct.Solution('gri30.yaml')

#%%

phi = 0.3

print(f"\033[1;36m### EQUIVALENCE RATIO: {phi} ###\033[0m")
gas.set_equivalence_ratio(phi, 'CH4', 'O2')



#%% Cálculo de velocidad de llama
T_0 = 300 # K

phi = 0.3
p = 10*ct.one_atm


print(f"\033[1;36m### EQUIVALENCE RATIO: {phi} ###\033[0m")
gas.set_equivalence_ratio(phi, 'CH4', 'O2: 1.0, N2: 3.76')


print(f"\033[1;36m## PRESSURE (Pa): {p} ##\033[0m")
gas.TP = T_0, p

print("Before flame creation")
gas()

# Llama
flame = ct.FreeFlame(gas=gas, width=0.001)
flame.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
print("Before flame solve")
gas()



flame.solve(loglevel=0, refine_grid=True, auto=True)
print("After flame solve")
gas()

print(f"\033[1;36m# LAMINAR BURNING SPEED (cm/s): {flame.velocity[0]*100} cm/s ###\033[0m")

#%%
# phi_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4]
# vel_list = {phi: [] for phi in phi_list}
# vel_list
# %%
