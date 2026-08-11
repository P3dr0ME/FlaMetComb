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
T_r = 300 # K


gas = ct.Solution('gri30.yaml')

#%% Cálculo de velocidad de llama

phi_list = [0.8, 1.0, 1.2, 1.4]
p_list = [np.float64(i*ct.one_atm) for i in np.geomspace(1,30,10)]
# Creo una lista con más densidad de puntos al principio. En Pa.

vel_list = {phi: [] for phi in phi_list}
# Diccionario de listas (de momento vacías) cuyas keys son los valores de phi_list.

oxidizer = "O2" if type == "oxi" else "O2:1, N2:3.76"

for j, phi in enumerate(phi_list):
    print(f"\033[1;36m### RATIO DE EQUIVALENCIA: {phi} ###\033[0m")

    for p in p_list:
        print(f"\033[1;36m   ## PRESIÓN (Pa): {p} ##\033[0m")

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

        print(f"\033[1;36m      # VELOCIDAD DE LLAMA (cm/s): {flame.velocity[0]*100} cm/s ###\033[0m")
        vel_list[phi].append(flame.velocity[0] * 100) # cm/s

        print("      Tiempo de cómputo %s segundos ---" % (time.time() - start_time))

print("TIEMPO DE CÓMPUTO TOTAL --- %s segundos ---" % (time.time() - start_time))

#%% Plot Flame Speed - phi
plt.figure(figsize=(8,8))
for phi in phi_list:
    plt.plot(p_list,
            vel_list[phi],
            label=r"$\phi$" + f"= {phi} (Cantera - GRI3.0)",
            marker=".")
plt.grid(True, which='both', alpha=0.5)
plt.xlabel("Presión, p [Pa] \n \n" + f"$T_{{0}}$ = {T_r} K")
plt.ylabel("Velocidad de llama [cm/s]")
plt.legend(loc='best', fontsize=10)
plt.savefig(f"plots/S_vs_p/S_vs_p_{type}.svg")
plt.show()
