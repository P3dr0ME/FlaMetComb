#%%
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science'])
import time
start_time = time.time()

# INPUTS:
type = "air" # "oxi" or "air"
p = ct.one_atm # atm
phi = 1.4

gas = ct.Solution('gri30.yaml')

#%% Cálculo de velocidad de llama

T_r_list = np. linspace(300,650,15)

vel_list = [{T_r: None for T_r in T_r_list}]
# Hay que vaciar vel_list al acabar cada bucle de presión
# para no juntar velocidades de todas las phi.

oxidizer = "O2" if type == "oxi" else "O2:1, N2:3.76"

for j, T_r in enumerate(T_r_list):
    print(f"\033[1;36m## TEMPERATURA (K): {T_r} ##\033[0m")

    gas.set_equivalence_ratio(phi, 'CH4', f'{oxidizer}')
    # Restablecer composición a la de los reactantes antes de resolver la llama de nuevo

    gas.TP = T_r, p

    # Llama
    flame = ct.FreeFlame(gas=gas, width=0.03)
    flame.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)

    if j != 0:
        flame.set_initial_guess(data=flame_sol_previa)

    flame.solve(loglevel=0, refine_grid=True, auto=True)
    flame_sol_previa = flame.to_array()

    print(f"\033[1;36m# VELOCIDAD DE LLAMA (cm/s): {flame.velocity[0]*100} cm/s ###\033[0m")
    vel_list[T_r].append(flame.velocity[0] * 100) # cm/s
    # plt.plot(flame.grid, flame.velocity)
    # plt.show()


    print("Tiempo de cómputo %s segundos ---" % (time.time() - start_time))

print("TIEMPO DE CÓMPUTO TOTAL --- %s segundos ---" % (time.time() - start_time))

#%% Plot Speed - phi
plt.figure(figsize=(8,8))
plt.plot(T_r_list,
        vel_list,
        label="GRI 3.0",
        marker="o")
plt.grid(True, which='both', alpha=0.5)
plt.xlabel("Temperatura [K] \n \n" + r"$\phi$" + f"= {phi}    p = {p} K")
plt.ylabel("Velocidad de llama [cm/s]")
plt.legend(loc='best', fontsize=10)
plt.savefig(f"./plots/S_vs_T/S_vs_T_{type}.svg")
plt.show()



# %%
