#%% PREÁMBULO
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(["science"])
import pandas as pd
import time
from scipy.optimize import newton

start_time = time.time()



#%% INPUTS
comburente = "air"  # "oxi" or "air"

# Inputs de phi
phi_list = np.geomspace(5, 8, num=10)

# Inputs de p
p_list = ( np.geomspace(1, 1, num=1)*ct.one_atm )

# Inputs T_r
T_r_list = np.linspace( 298, 298, num=1 ) # K



#%% ESPECIES Y GASES
# Oxidante
oxidizer = ( "O2" if comburente == "oxi" else f"O2:1, N2:{79/21}" )

species_dict = { S.name: S for S in ct.Species.list_from_file("gri30.yaml") }

list_ideal_species = ( "CH4", "O2", "N2", "CO2", "H2O" )
ideal_species = [ species_dict[S] for S in list_ideal_species ]

# Gases
gas_real = ct.Solution( "gri30.yaml" )
gas_ideal = ct.Solution(
    thermo="ideal-gas",
    species=ideal_species,
    transport_model="mixture-averaged",
    kinetics="gas"
)

# Masas molares
M_CH4 = 16.04 / 1e3 # kg/mol
M_O2 = 32.00 / 1e3 # kg/mol
M_N2 = 28.014 / 1e3 # kg/mol
M_aire = ( M_O2 + 79/21 * M_N2 ) # kg/mol

# Dosado estequiométrico
f_s = ( M_CH4 / (2 * M_aire) if comburente == "air" else M_CH4 / (2 * M_O2) )

# Poder calorífico
LHV = 50.048e6 # J/kg
# Springer, Appendix 1.

# Parámetros de la reacción global
A_0 = 8.3e5 # (mol/cm3)^(1-a-b) / s
a = -0.3
b = 1.3

# Energía de activación
E_a = 30 # kcal/mol
# Constante universal de gases
R_u_a = 1.98591e-3 # kcal/(mol K)
R_u_SI = 8.31447 # J/(mol K)

# Temperatura de activación
T_act = E_a / R_u_a # K



#%% FUNCIÓN: LLAMA DE CANTERA
def Llama_Cantera( phi, p, T_r, flame_sol_previa=None, loglevel=0 ):
    gas_real.TP = T_r, p
    gas_real.set_equivalence_ratio(phi, "CH4", oxidizer )

    flame = ct.FreeFlame( gas=gas_real, width=0.03 )

    flame.set_refine_criteria( ratio=3, slope=0.06, curve=0.12 )

    # Utilizar solución previa como punto de partida
    if flame_sol_previa is not None:
        flame.set_initial_guess( data=flame_sol_previa )

    flame.solve( loglevel=loglevel, refine_grid=True, auto=True )

    return flame



#%% FUNCIÓN: LLAMA ANALÍTICA
# cp^_i(T) según polinomios NASA7/9 de gri30.yaml
def cp_molar(sp, T):
    return gas_real.species(f"{sp}").thermo.cp(T) / 1e3 # J/mol/K

def Llama_analitica( phi, p, T_r):
    T_ig_analitica = ( T_act - np.sqrt(T_act**2 - 4*T_act*T_r) ) / 2

    # Moles de reactantes
    eps = 0 if phi<=1 else phi-1
    delta = 2/phi-2 if phi<=1 else 0

    n_CH4_r = 1+eps
    n_O2_r = 2+delta
    n_N2_r = 79/21*(2+delta) if comburente=="air" else 0

    # Fracciones molares de reactantes
    X_CH4_r = n_CH4_r/(n_CH4_r+n_O2_r+n_N2_r)
    X_O2_r = n_O2_r/(n_CH4_r+n_O2_r+n_N2_r)
    X_N2_r = n_N2_r/(n_CH4_r+n_O2_r+n_N2_r)

    # Concentraciones de reactantes
    C_CH4_r = X_CH4_r * p / (R_u_SI * T_r) / 1e6 # mol/cm3
    C_O2_r = X_O2_r * p / (R_u_SI * T_r) / 1e6 # mol/cm3

    # Densidad de reactantes
    M_r = X_CH4_r*M_CH4 + X_O2_r*M_O2 + X_N2_r*M_N2
    rho_r = p*M_r/(R_u_SI*T_r)

    # Moles de productos
    n_CH4_p = eps
    n_O2_p = delta
    n_N2_p = 79/21*(2+delta) if comburente == "air" else 0
    n_CO2_p = 1
    n_H2O_p = 2

    # Masa de la mezcla
    m_mezcla = ( n_CH4_r*M_CH4 + n_O2_r*M_O2 + n_N2_r*M_N2 )

    # Calor de combustión por kg de mezcla
    q0_p = (phi*f_s*LHV)/(1+phi*f_s) if phi<=1 else (f_s*LHV)/(1+phi*f_s)

    # T_ad (Método 1 pero con cp de cantera)
    def f(T_p):
        T_ave_Tad = (T_p + T_r) / 2
        cp_p = (
              n_CH4_p * cp_molar("CH4", T_ave_Tad)
            + n_O2_p  * cp_molar("O2",  T_ave_Tad)
            + n_N2_p  * cp_molar("N2",  T_ave_Tad)
            + n_CO2_p * cp_molar("CO2", T_ave_Tad)
            + n_H2O_p * cp_molar("H2O", T_ave_Tad)
        ) / m_mezcla
        # Calculo cp_p igual que Método 1
        return (T_p-T_r)*cp_p - q0_p
    T_ad_m1 = newton(f, x0=1000, tol=1e-12)

    # Propiedades a T_ave
    T_ave = (T_ad_m1 + T_ig_analitica) / 2
    gas_ideal.TP = T_ave, p
    k_ave = gas_ideal.thermal_conductivity # W/(m·K)
    cp_ave = gas_ideal.cp_mass
        # No es cp_p de T_ad, es cp promedio de zona que tiene tanto reactantes como productos y a T_ave diferente.
        # Mejor aproximarlo con función de cantera como k_ave.

    # Difusividad térmica
    alpha_ave = k_ave / (rho_r * cp_ave) # m2/s

    # Ritmo de consumo del fuel promedio
    r_f_ave = A_0 * C_CH4_r**a * C_O2_r**b * np.exp( -T_act / T_ave ) # mol/(cm3 s)

    # Tiempo químico
    tau_q = C_CH4_r / r_f_ave # s

    # Velocidad de llama analítica
    S_L_analitica = 100 * np.sqrt( (alpha_ave/ tau_q) * (T_ad_m1 - T_ig_analitica)/(T_ig_analitica - T_r) )  # cm/s

    # Diccionario de datos y resultados analíticos
    return {
        "T_ig_analitica": T_ig_analitica,
        "T_ave": T_ave,
        "T_ad_m1": T_ad_m1,
        "rho_r": rho_r,
        "C_CH4_r": C_CH4_r,
        "C_O2_r": C_O2_r,
        "cp_ave": cp_ave,
        "k_ave": k_ave,
        "alpha_ave": alpha_ave,
        "r_f_ave": r_f_ave,
        "tau_q": tau_q,
        "S_L_analitica": S_L_analitica
    }


#%% FUNCIÓN DE RESULTADOS
def CALCULO_LLAMA( phi, p, T_r, flame_sol_previa=None, loglevel=0 ):
    print( f"\033[1;36m phi = {phi:.6g} | p = {p/ct.one_atm:.6g} atm | T_r = {T_r:.6g} K \033[0m" )

    # Cantera
    flame = Llama_Cantera( phi=phi, p=p, T_r=T_r, flame_sol_previa=flame_sol_previa, loglevel=loglevel )

    # Analítico
    res_analiticos = Llama_analitica( phi=phi, p=p, T_r=T_r)

    #Unión de resultados
    res = {
        (phi, p, T_r):
        {
        "S_L_Cantera": flame.velocity[0]*100, # cm/s
        "T_p_Cantera": flame.T[-1], # K
        "T_ig_Cantera": flame.T[0], # K
        **res_analiticos
        }
    }
    # la key es una tuple (como una lista pero inmutable) con los datos (lists no pueden ser dict keys)
    # el value es un diccionario con sus resultados

    print( f"    Tiempo total = {time.time() - start_time:.2f} s" )

    return res, flame.to_array() # Esta flame.to_array se realimenta a la propia función en siguiente bucle como flame_sol_previa




#%% CÁLCULO DE LLAMA
RESULTADOS = {}
for phi in phi_list:
    for T_r in T_r_list:
        flame_sol_previa=None
        for p in p_list:
            if (phi, p, T_r) not in RESULTADOS: # Calcular solo si no se ha calculado ya para los mismos datos
                res, flame_sol_previa = CALCULO_LLAMA( phi=phi, p=p, T_r=T_r, flame_sol_previa=flame_sol_previa, loglevel=1 )
                RESULTADOS.update(res)


#%% S vs phi — CURVA Y TABLA
p_ref, T_ref = 1*ct.one_atm, 298

datos = pd.DataFrame.from_dict(
    {phi: r for (phi, p, T), r in RESULTADOS.items() if p == p_ref and T == T_ref},
    orient="index"
).rename_axis("phi")

datos = datos.sort_index()

# Curva
plt.figure(figsize=(8, 8))
plt.plot(datos.index, datos["S_L_Cantera"], ".-", label="Cantera (GRI3.0)")
plt.plot(datos.index, datos["S_L_analitica"], ".-", label="Analítica")

plt.xlabel(r"$\phi$")
plt.ylabel(r"$S_L$ (cm/s)")
plt.title(
    f"{r"\bf{Oxígeno}" if comburente == 'oxi' else r"\bf{Aire}"} \n"
    fr"$p = {p*ct.one_atm}$ Pa $\quad T_r = {T_r}$ K",
    fontsize=11,
    pad=15
)
plt.grid(True, which="both", alpha=0.5)
plt.legend()

# plt.savefig(f"./Res/S_L/S_vs_phi_{comburente}.svg")
plt.show()

# Tabla
Tabla_S_phi = datos[["S_L_Cantera", "S_L_analitica"]].rename(columns={
    "S_L_Cantera": "S_L Cantera (cm/s)",
    "S_L_analitica": "S_L analítica (cm/s)"
})

print(Tabla_S_phi)

Tabla_S_phi.to_csv(
    f"./Res/S_L/Tabla_S_vs_phi_{comburente}.csv"
)
