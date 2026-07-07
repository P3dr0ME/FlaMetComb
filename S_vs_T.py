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
p = ct.one_atm # atm
phi = 1.4

T_list = np. linspace(300,650,15)


vel_list = []
# Hay que vaciar vel_list al acabar cada bucle de presión
# para no juntar velocidades de todas las phi.
for j, T in enumerate(T_list):
    print(f"\033[1;36m## Temperature (K): {T} ##\033[0m")

    gas_mix_incomplete.set_equivalence_ratio(phi, 'CH4', 'O2: 1.0, N2: 3.76')
    # Restablecer composición a la de los reactantes antes de resolver la llama de nuevo

    gas_mix_incomplete.TP = T, p

    # Llama
    flame = ct.FreeFlame(gas=gas_mix_incomplete, width=0.03)
    flame.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)

    if j != 0:
        flame.set_initial_guess(data=flame_sol_previa)

    flame.solve(loglevel=0, refine_grid=True, auto=True)
    flame_sol_previa = flame.to_array()

    print(f"\033[1;36m# LAMINAR BURNING SPEED (cm/s): {flame.velocity[0]*100} cm/s ###\033[0m")
    vel_list.append(flame.velocity[0] * 100) # cm/s
    # plt.plot(flame.grid, flame.velocity)
    # plt.show()


    print("Run time %s seconds ---" % (time.time() - start_time))

#%% Plot Speed - phi
plt.figure(figsize=(8,8))
plt.plot(T_list,
        vel_list,
        label="Cantera (GRI 3.0)",
        marker="o")
plt.grid(True, which='both', alpha=0.5)
plt.xlabel("Temperature, T [K] \n \n" + r"$\phi$" + f"= {phi}    p = {p} K")
plt.ylabel("Flame speed [cm/s]")
plt.legend(loc='best', fontsize=10)
plt.savefig("S_vs_T.svg")
plt.show()

print("Complete run time --- %s seconds ---" % (time.time() - start_time))

# %%
