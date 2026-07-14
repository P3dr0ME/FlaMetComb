#%%
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science'])
import time
start_time = time.time()


gas = ct.Solution('gri30.yaml')

#%% Cálculo de velocidad de llama
T_0 = 300 # K

phi_list = [0.8, 1.0, 1.2, 1.4]
p_list = [np.float64(i*ct.one_atm) for i in np.geomspace(1,30,10)]
# Creo una lista con más densidad de puntos al principio. En Pa.

vel_list = {phi: [] for phi in phi_list}
# Diccionario de listas (de momento vacías) cuyas keys son los valores de phi_list.
for j, phi in enumerate(phi_list):
    print(f"\033[1;36m### EQUIVALENCE RATIO: {phi} ###\033[0m")

    for p in p_list:
        print(f"\033[1;36m   ## PRESSURE (Pa): {p} ##\033[0m")

        gas.set_equivalence_ratio(phi, 'CH4', 'O2')
        # Restablecer composición a la de los reactantes antes de resolver la llama de nuevo

        gas.TP = T_0, p

        # Llama
        flame = ct.FreeFlame(gas=gas, width=0.03)
        flame.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)

        if j != 0:
            flame.set_initial_guess(data=flame_sol_previa)

        flame.solve(loglevel=0, refine_grid=True, auto=True)
        flame_sol_previa = flame.to_array()

        print(f"\033[1;36m      # LAMINAR BURNING SPEED (cm/s): {flame.velocity[0]*100} cm/s ###\033[0m")
        vel_list[phi].append(flame.velocity[0] * 100) # cm/s

        print("      Run time %s seconds ---" % (time.time() - start_time))

print("Complete run time --- %s seconds ---" % (time.time() - start_time))

#%% Plot Flame Speed - phi
plt.figure(figsize=(8,8))
for phi in phi_list:
    plt.plot(p_list,
            vel_list[phi],
            label=r"$\phi$" + f"= {phi} (GRI3.0)",
            marker=".")
plt.grid(True, which='both', alpha=0.5)
plt.xlabel("Pressure, p [Pa] \n \n" + f"$T_{{0}}$ = {T_0} K")
plt.ylabel("Flame speed [cm/s]")
plt.legend(loc='best', fontsize=10)
plt.savefig(f"plots/S_vs_p_oxi.svg")
plt.show()
