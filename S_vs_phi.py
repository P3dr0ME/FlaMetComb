import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','ieee'])


gas_mix_incomplete = ct.Solution('gri30.yaml')

#%%
T_0 = 298 # K
p = ct.one_atm

vel_list = []
phi_list = np.linspace(0.5, 1.8, 14)
for j, phi in enumerate(phi_list):
    # Contador j no se usa
    print(f"\033[1;36m### EQUIVALENCE RATIO: {phi} ###\033[0m")

    gas_mix_incomplete.TP = T_0, p

    gas_mix_incomplete.set_equivalence_ratio(phi, 'CH4', 'O2: 1.0, N2: 3.76')

    # Llama
    flame = ct.FreeFlame(gas=gas_mix_incomplete, width=0.03)
        # Clase FreeFlame --> llama de premezcla 1D
        # width declara grid en intervalo [0,width]
        # y que solver determine ptos. intermedios.
    flame.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
        # Criterios que solver seguirá para refinar grid.
        # Por ejemplo, slope dice que si dif. máx de valores en nodos adyacentes
        # supera el 6% de la máx diferencia del perfil, refina.
    # flame.soret_enabled = True  # optionally enable
    # print(flame.grid)
    flame.solve(loglevel=1, refine_grid=True, auto=True)
        # Método que resuelve ecs. (las de fluidos en dif. finitas).
    print(f"\033[1;36m### LAMINAR BURNING SPEED: {flame.velocity[0]*100} cm/s ###\033[0m")
    vel_list.append(flame.velocity[0] * 100) # cm/s
        # flame.velocity[0] = velocidad en primer grid point (inlet).



#%% Plot Speed - phi
plt.figure(figsize=(8,8))

plt.plot(phi_list,
        vel_list,
        label="speed (cm/s)",
        lw=2)

plt.grid(True, which='both', alpha=0.5)

plt.xlabel("Equivalence ratio, $\phi$")
plt.ylabel("Flame speed [cm/s]")

plt.legend(loc='best', fontsize=10)

plt.show()
