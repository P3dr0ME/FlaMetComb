#%%
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science'])
import time
start_time = time.time()

# INPUTS:
type = "oxi" # "oxi" or "air"
T_R = 298 # K
p = ct.one_atm

N = 20
phi_begin = 0.5
phi_end = 1.8
N = 20

# Definición del gas con el modelo GRI3.0
gas = ct.Solution('gri30.yaml')

#%% Cálculo de velocidad de llama

phi_list = [phi_begin + (phi_end-phi_begin) * (n/(N-1) + np.sin(2*np.pi*n/(N-1))/(2*np.pi)) for n in range(N)]
# Puntos que se concentran en el centro del intervalo (donde más varía la curva).

vel_list = {phi: None for phi in phi_list}
T_P_llama = {phi: None for phi in phi_list}
cp_ave = [None for phi in phi_list]

oxidizer = "O2" if type == "oxi" else "O2:1, N2:3.76"

for j, phi in enumerate(phi_list):
    # Contador j vale 0, 1, 2 (las posiciones de la lista iterada, en vez de valores)
    print(f"\033[1;36m### RATIO DE EQUIVALENCIA: {phi} ###\033[0m")

    gas.TP = T_R, p
    gas.set_equivalence_ratio(phi, 'CH4', f'{oxidizer}')

    # Llama
    flame = ct.FreeFlame(gas=gas, width=0.03)
        # Clase FreeFlame --> llama de premezcla 1D
        # width crea grid en intervalo [0,width]
        # y que solver determine ptos. intermedios.

    flame.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
        # Criterios que solver seguirá para refinar grid.
        # Por ejemplo, slope dice que si dif. máx de valores en nodos adyacentes
        # supera el 6% de la máx diferencia del perfil, añade puntos intermedios.
    # print(flame.grid)

    if j != 0:
        flame.set_initial_guess(data=flame_sol_previa)
        # Uso solución de llama con el phi anterior como punto de partida para este
        # Opcional. Tiempo cómputo se reduce 49%.

    flame.solve(loglevel=1, refine_grid=True, auto=True)
        # Método del objeto flame que resuelve ecs. de fluidos en dif. finitas.
    flame_sol_previa = flame.to_array()
        # Guardo la solución de un bucle en esta variable.

    print(f"\033[1;36m### LAMINAR BURNING SPEED (cm/s): {flame.velocity[0]*100} cm/s ###\033[0m")

    vel_list[phi] = flame.velocity[0] * 100 # cm/s
        # flame.velocity[0] = velocidad en primer grid point (inlet).
    T_P_llama[phi] = flame.T[-1] # K. Sirve para aproximar T_ad analítica.

    print("Tiempo de cómputo %s segundos ---" % (time.time() - start_time))

print("TIEMPO DE CÓMPUTO TOTAL --- %s segundos ---" % (time.time() - start_time))

#%% Cálculo analítico de la velocidad de llama
species_dict = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}
ideal_species = [species_dict[S] for S in ("CH4", "O2", "N2", "CO2", "H2O")]
gas_ideal_aux  = ct.Solution(thermo="ideal-gas",
                    species=ideal_species,
                    transport_model='mixture-averaged',
                    kinetics='gas')

M_CH4 = 16.04e-3 # kg/mol
M_aire = 137.33e-3 # kg/mol
M_O2 = 32e-3 # kg/mol
f_s = 1*M_CH4/(2*M_aire) if type == "air" else 1*M_CH4/(2*M_O2) # Dosado estequiométrico
LHV = 50.048e6 # J/kg. Springer, Appendix 1

A_0 = 8.3e5
a = -0.3
b = 1.3
E_a = 30 # kcal/mol
R_u_a = 1.98591e-3 # kcal/(mol·K). Constante universal de los gases.
R_u_SI = 8.31447 # J/(mol·K). Constante universal de los gases.
T_act = E_a/R_u_a # K. Temperatura de activación

T_ig = (T_act - np.sqrt(T_act**2 - 4*T_act*T_R)) / 2

T_ad_analytic = {phi: None for phi in phi_list}
vel_list_analytic = {phi: None for phi in phi_list}

for i in range(len(phi_list)):
    phi=phi_list[i]
    #Cálculo de T_ad analítica
    gas_ideal_aux.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")

    C_CH4 = (gas_ideal_aux["CH4"].X[0]*p)/(R_u_SI*T_R) # Fracción molar de CH4
    C_02 = (gas_ideal_aux["O2"].X[0]*p)/(R_u_SI*T_R) # Fracción molar de O2



    T_ave = (T_P_llama[phi]+T_R)/2
    gas_ideal_aux.TP = T_ave, p
    cp_ave[phi] = gas_ideal_aux.cp # lista de cp de productos como base de datos para calcular T_ad analítica
    if 0<phi<=1: # pobre
        T_ad_analytic[phi] = T_R + ( phi*f_s*LHV ) / ( (1+phi*f_s)*cp_ave[phi] )
    elif phi>1: # rica
        T_ad_analytic[phi] = T_R + ( f_s*LHV ) / ( (1+phi*f_s)*cp_ave[phi] )

    k_ave = gas_ideal_aux.thermal_conductivity # W/m/K. Conductividad térmica
    rho_ave = gas_ideal_aux.density
    alpha_ave = k_ave/(rho_ave*cp_ave[phi]) # m2/s. Difusividad térmica

    r_f_ave = 1.0 * A_0*(C_CH4)^a*(C_02)^b*exp(-T_act/T_ave)
    tau_q = C_CH4/r_f_ave # s. Tiempo químico promedio

    vel_list_analytic[phi] = np.sqrt( (alpha_ave/tau_q) * (T_ad_analytic[phi]-T_ig)/(T_ig-T_R) ) * 100 # cm/s

#%% Plot Speed - phi
plt.figure(figsize=(8,8))

plt.plot(phi_list,
        vel_list.values(),
        label="Cantera (GRI 3.0)",
        marker=".")

plt.plot(phi_list,
        vel_list_analytic.values(),
        label="Cantera (GRI 3.0)",
        marker="")

plt.grid(True, which='both', alpha=0.5)

plt.xlabel("Ratio de equivalencia, " + r"$\phi$" + f"\n \n p = {p} Pa    $T_{{0}}$ = {T_R} K")
plt.ylabel("Velocidad de llama [cm/s]")

plt.legend(loc='best', fontsize=10)

plt.savefig(f"./plots/S_vs_phi/S_vs_phi_{type}.svg")
plt.show()
