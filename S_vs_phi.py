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
T_r = 298 # K
p = ct.one_atm
phi_begin = 0.3; phi_end = 1.8; N = 5

# Definición del gas con el modelo GRI3.0
gas_real = ct.Solution('gri30.yaml')

#%% Cálculo de velocidad de llama

oxidizer = "O2" if type == "oxi" else "O2:1, N2:3.76"

phi_list = [phi_begin + (phi_end-phi_begin) * (n/(N-1) + np.sin(2*np.pi*n/(N-1))/(2*np.pi)) for n in range(N)]
# Puntos que se concentran en el centro del intervalo (donde más varía la curva).

# Inicialización de diccionarios de forma compacta
vel_list, T_P_llama, T_ig_llama, q_p, cp_real = (
    {phi: None for phi in phi_list} for _ in range(4)
)

for j, phi in enumerate(phi_list):
    # Contador j vale 0, 1, 2 (las posiciones de la lista iterada, en vez de valores)
    print(f"\033[1;36m### RATIO DE EQUIVALENCIA: {phi} ###\033[0m")

    gas_real.TP = T_r, p
    gas_real.set_equivalence_ratio(phi, 'CH4', f'{oxidizer}')

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
    T_ig_llama[phi] = flame.T[0] # K. Sirve para comparar T_ig de Cantera y T_ig analítica.

    cp_real[phi] = flame.cp[-1]

    print("Tiempo de cómputo %s segundos ---" % (time.time() - start_time))

print("TIEMPO DE CÓMPUTO TOTAL --- %s segundos ---" % (time.time() - start_time))

#%% Cálculo analítico de la velocidad de llama
species_dict = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}
ideal_species = [species_dict[S] for S in ("CH4", "O2", "N2", "CO2", "H2O")]
gas_ideal_aux  = ct.Solution(thermo="ideal-gas",
                    species=ideal_species,
                    transport_model='mixture-averaged',
                    kinetics='gas')

M_CH4 = 16.04e-3; M_aire = 137.33e-3; M_O2 = 32e-3 # kg/mol
f_s = 1*M_CH4/(2*M_aire) if type == "air" else 1*M_CH4/(2*M_O2) # Dosado estequiométrico
LHV = 50.048e6 # J/kg. Springer, Appendix 1

A_0 = 8.3e5 # (mol/cm^3)^(1-a-b) /s
a = -0.3; b = 1.3
E_a = 30 # kcal/mol
R_u_a = 1.98591e-3 # kcal/(mol·K). Constante universal de los gases.
R_u_SI = 8.31447 # J/(mol·K). Constante universal de los gases.
T_act = E_a/R_u_a # K. Temperatura de activación

T_ig = (T_act - np.sqrt(T_act**2 - 4*T_act*T_r)) / 2

T_ave = {phi: (T_P_llama[phi]+T_ig)/2 for phi in phi_list}

# Inicialización de diccionarios de forma compacta
T_ad_analytic, vel_list_analytic, cp_ave, k_ave, rho_ave_reactants = (
    {phi: None for phi in phi_list} for _ in range(5)
)

for i in range(len(phi_list)):
    phi=phi_list[i]
    #Cálculo de T_ad analítica
    gas_ideal_aux.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
    gas_ideal_aux.TP = T_r, p
    rho_ave_reactants[phi] = gas_ideal_aux.density
    # La densidad se aproxima con la de reactantes. Por eso se calcula aquí a la T_r.

    X_CH4 = gas_ideal_aux.X[gas_ideal_aux.species_index("CH4")]
    X_O2 = gas_ideal_aux.X[gas_ideal_aux.species_index("O2")]

    C_CH4 = (X_CH4*p)/(R_u_SI*T_r) / 1e6 # mol/cm^3
    C_O2 = (X_O2*p)/(R_u_SI*T_r) / 1e6 # mol/cm^3

    gas_ideal_aux.TP = T_ave[phi], p
    cp_ave[phi] = gas_ideal_aux.cp # lista de cp de productos como base de datos para calcular T_ad analítica
    # cp_ave[phi] = 2400
    if 0<phi<=1: # pobre
        T_ad_analytic[phi] = T_r + ( phi*f_s*LHV ) / ( (1+phi*f_s)*cp_ave[phi] )
    elif phi>1: # rica
        T_ad_analytic[phi] = T_r + ( f_s*LHV ) / ( (1+phi*f_s)*cp_ave[phi] )

    k_ave[phi] = gas_ideal_aux.thermal_conductivity # W/m/K. Conductividad térmica

    alpha_ave = k_ave[phi]/(rho_ave_reactants[phi]*cp_ave[phi]) # m^2/s. Difusividad térmica

    r_f_ave = 1.0 * A_0 * C_CH4**a * C_O2**b * np.exp(-T_act/T_ave[phi])
    tau_q = C_CH4/r_f_ave # s. Tiempo químico promedio


    vel_list_analytic[phi] = np.sqrt( (alpha_ave/tau_q) * (T_ad_analytic[phi]-T_ig)/(T_ig-T_r) ) * 100 # cm/s

#%% Plot Speed - phi
plt.figure(figsize=(8,8))

plt.title(
    r"\bf{Velocidad\ de\ llama\ frente\ a\ ratio\ de\ equivalencia}" + "\n"
    f"{'Oxígeno' if type == 'oxi' else 'Aire'} \n"
    fr"$p = {p}$ Pa $\quad T_{{0}} = {T_r}$ K",
    fontsize=11,
    pad=15
)

plt.plot(phi_list,
        vel_list.values(),
        label="Real (Cantera - GRI 3.0)",
        marker=".")

plt.plot(phi_list,
        vel_list_analytic.values(),
        label="Analítica",
        marker="")

plt.grid(True, which='both', alpha=0.5)

plt.xlabel("Ratio de equivalencia, " + r"$\phi$")
plt.ylabel("Velocidad de llama [cm/s]")

plt.legend(loc='best', fontsize=10)

plt.savefig(f"./plots/S_vs_phi/S_vs_phi_{type}.svg")
plt.show()

#%% Plot T vs x
plt.figure(figsize=(8,8))

plt.plot(flame.grid,
        flame.T,

        marker="")

plt.grid(True, which='both', alpha=0.5)

plt.xlabel("x")
plt.ylabel("T [K]")

plt.legend(loc='best', fontsize=10)

plt.show()

# %% Plot T_ad vs phi
plt.figure(figsize=(8,8))

plt.plot(phi_list,
        T_ad_analytic.values(),
        label="Analítica",
        marker=""
        )

plt.grid(True, which='both', alpha=0.5)

plt.xlabel("Ratio de equivalencia, " + r"$\phi$")
plt.ylabel("T_ad [K]")

plt.legend(loc='best', fontsize=10)

plt.show()
