#%% PREÁMBULO
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(["science"])
import pandas as pd
from scipy.optimize import newton
import os



#%% INPUTS
comburente = "oxi"  # "oxi" or "air"
borrar_csv_antiguo = False

# Seleccionar curva a representar (descomentar 1)
# "1 {dato}" --> dato_list debe contener 1 solo elemento entre corchetes []

# # S vs p | 1 phi, varios T_r
# selector_plot = "S vs p | 1 phi, varios T_r"
# phi_list = [1.00]
# p_list = np.linspace(0.1e6, 6e6, num=5) # Pa
# p_list = [0.2e6, 0.3e6, 0.5e6, 0.7e6, ct.one_atm]
# T_r_list = [300, 373, 500] # K

# S vs phi | 1 T_r, varios p
selector_plot = "S vs phi | 1 T_r, varios p"
phi_list = np.linspace(0.1, 2.0, num=20)
p_list = [ct.one_atm] # Pa
T_r_list = [298] # K

# # S vs phi | 1 p, varios T_r
# selector_plot = "S vs phi | 1 p, varios T_r"
# phi_list =
# p_list = # Pa
# T_r_list = # K

# # S vs T_r | 1 phi, varios p
# selector_plot = "S vs T_r | 1 phi, varios p"
# phi_list = [1.00]
# p_list = [0.1e6, 0.5e6, 1e6] # Pa
# T_r_list = np.linspace(300, 700, num=8) # K

#%% ESPECIES Y GASES
# Oxidante
oxidizer = ( "O2" if comburente == "oxi" else f"O2:1, N2:{79/21}" )

# Gas
gas_real = ct.Solution( "gri30.yaml" )

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
        # Sale casi idéntico a calcular gas_ideal.density

    # Propiedades a T_ave
    gas_real.TP = T_r, p
    gas_real.set_equivalence_ratio(phi, "CH4", oxidizer )

    gas_real.equilibrate("HP")
    T_ad_real = gas_real.T

    T_ave = (T_ad_real + T_ig_analitica) / 2
    gas_real.TP = T_ave, p

    k_ave = gas_real.thermal_conductivity # W/(m·K)
    cp_ave = gas_real.cp_mass
        # En las llamas, calcular sin tener en cuenta cinética o disociación (ideal) no es buena hipótesis
        # Se depende de una buena aproximación de cp, k, T_ad, etc., que se tomará del Método 3 real.

    # Difusividad térmica
    alpha_ave = k_ave / (rho_r * cp_ave) # m2/s

    # Ritmo de consumo del fuel promedio
    r_f_ave = A_0 * C_CH4_r**a * C_O2_r**b * np.exp(-T_act/T_ave) # mol/(cm3 s)

    # Tiempo químico
    tau_q = C_CH4_r / r_f_ave # s

    # Velocidad de llama analítica
    S_L_analitica = 100 * np.sqrt( (alpha_ave/ tau_q) * (T_ad_real - T_ig_analitica)/(T_ig_analitica - T_r) )  # cm/s

    # Diccionario de datos y resultados analíticos
    return {
        "T_ig_analitica": T_ig_analitica,
        "T_ave": T_ave,
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


#%% FUNCIÓN DE CÁLCULO CONJUNTO
def CALCULO_LLAMA( phi, p, T_r, flame_sol_previa=None, loglevel=0 ):
    print( f"\033[1;36m phi = {phi:.6g} | p = {p/1e6} MPa | T_r = {T_r:.6g} K \033[0m" )

    # Cantera
    try:
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
        # key es una tuple con los datos (lists no pueden ser dict keys)
        # value es un diccionario con resultados

        return res, flame.to_array() # Esta flame.to_array se realimenta a la propia función en siguiente bucle como flame_sol_previa

    except Exception as e: # Si falla flame porque no converge, que avise dónde y siga.
        print(f"\033[1;31m FALLO: phi={phi}, p={p/1e6} MPa, T_r={T_r} K \033[0m")
        print(f"    {type(e).__name__}: {e}")

        return {}, None



#%% CÁLCULO DE LLAMA
RESULTADOS = {}
for phi in phi_list:
    for T_r in T_r_list:
        flame_sol_previa=None
        for p in p_list:
            if (phi, p, T_r) not in RESULTADOS: # Calcular solo si no se ha calculado ya para los mismos datos
                res, flame_sol_previa = CALCULO_LLAMA( phi=phi, p=p, T_r=T_r, flame_sol_previa=flame_sol_previa, loglevel=1 )
                RESULTADOS.update(res)



#%% TABLA CSV
ruta = f"./Res/S_L/RESULTADOS_{comburente}.csv"

# Borrar CSV antiguo (si existe)
if borrar_csv_antiguo:
    os.path.exists(ruta) and os.remove(ruta)

# Resultados nuevos
RESULTADOS_CSV_NUEVOS = pd.DataFrame([
    {
        "phi": phi,
        "p": p,
        "T_r": T_r,
        "S_L_Cantera": datos["S_L_Cantera"],
        "S_L_analitica": datos["S_L_analitica"]
    }
    for (phi, p, T_r), datos in RESULTADOS.items()
])

# Leer resultados anteriores y combinar
try:
    RESULTADOS_CSV_ANTIGUOS = pd.read_csv(ruta)
    RESULTADOS_CSV = pd.concat([RESULTADOS_CSV_ANTIGUOS, RESULTADOS_CSV_NUEVOS], ignore_index=True)
except FileNotFoundError:
    RESULTADOS_CSV = RESULTADOS_CSV_NUEVOS

# Eliminar posibles duplicados
RESULTADOS_CSV = RESULTADOS_CSV.drop_duplicates( subset=["phi", "p", "T_r"], keep="last" )

# Guardar
RESULTADOS_CSV.to_csv(ruta, index=False)



#%% PLOT
if selector_plot == "S vs p | 1 phi, varios T_r":
    # S vs p | 1 phi, varios T_r
    plt.figure(figsize=(8,8))
    plt.title(
    f"{r"\bf{Oxígeno}" if comburente == 'oxi' else r"\bf{Aire}"} \n"
    fr"$\phi = {phi_list[0]}$",
    fontsize=11,
    pad=15
    )

    for T_r in T_r_list:
        d = RESULTADOS_CSV[(np.isclose(RESULTADOS_CSV.phi, phi_list)) & (np.isclose(RESULTADOS_CSV.T_r, T_r))].sort_values("p")
        color = plt.gca()._get_lines.get_next_color()
        plt.plot(d.p/1e6, d.S_L_Cantera/100, ".-", label=fr"Cantera, $T_r={T_r}$ K", color=color)
        plt.plot(d.p/1e6, d.S_L_analitica/100, "--", label=fr"Analítica, $T_r={T_r}$ K", color=color)

    plt.xlabel("$p$ (MPa)"); plt.ylabel("$S_L$ (m/s)")
    plt.grid(alpha=.5); plt.legend()

    plt.savefig(f"./Res/S_L/S_vs_p_{comburente}.svg")
    plt.show()

elif selector_plot == "S vs phi | 1 T_r, varios p":
    # S vs phi | 1 T_r, varios p
    plt.figure(figsize=(8,6))
    plt.title(
    f"{r"\bf{Oxígeno}" if comburente == 'oxi' else r"\bf{Aire}"} \n"
    fr"$T_r = {T_r_list[0]}$ K",
    fontsize=11,
    pad=15
    )

    for p in p_list:
        d = RESULTADOS_CSV[(np.isclose(RESULTADOS_CSV.T_r, T_r_list)) & (np.isclose(RESULTADOS_CSV.p, p))].sort_values("phi")
        color = plt.gca()._get_lines.get_next_color()
        plt.plot(d.phi, d.S_L_Cantera, ".-", label=fr"Cantera, $p={p/1e6}$ MPa", color=color)
        plt.plot(d.phi, d.S_L_analitica, "--", label=fr"Analítica, $p={p/1e6}$ MPa", color=color)

    plt.xlabel(r"$\phi$"); plt.ylabel("$S_L$ (cm/s)")
    plt.grid(alpha=.5); plt.legend()

    plt.savefig(f"./Res/S_L/S_vs_phi_varios_p_{comburente}.svg")
    plt.show()

elif selector_plot == "S vs phi | 1 p, varios T_r":
    # S vs phi | 1 p, varios T_r
    plt.figure(figsize=(8,8))
    plt.title(
    f"{r"\bf{Oxígeno}" if comburente == 'oxi' else r"\bf{Aire}"} \n"
    fr"$p = {p_list[0]/1e6}$ MPa",
    fontsize=11,
    pad=15
    )

    for T_r in T_r_list:
        d = RESULTADOS_CSV[(np.isclose(RESULTADOS_CSV.p, p_list)) & (np.isclose(RESULTADOS_CSV.T_r, T_r))].sort_values("phi")
        color = plt.gca()._get_lines.get_next_color()
        plt.plot(d.phi, d.S_L_Cantera, ".-", label=fr"Cantera, $T_r={T_r}$ K", color=color)
        plt.plot(d.phi, d.S_L_analitica, "--", label=fr"Analítica, $T_r={T_r}$ K", color=color)

    plt.xlabel(r"$\phi$"); plt.ylabel("$S_L$ (cm/s)")
    plt.grid(alpha=.5); plt.legend()

    plt.savefig(f"./Res/S_L/S_vs_phi_varios_T_r_{comburente}.svg")
    plt.show()

elif selector_plot == "S vs T_r | 1 phi, varios p":
    # S vs T_r | 1 phi, varios p
    plt.figure(figsize=(8,8))
    plt.title(
    f"{r"\bf{Oxígeno}" if comburente == 'oxi' else r"\bf{Aire}"} \n"
    fr"$\phi = {phi_list[0]}$",
    fontsize=11,
    pad=15
    )

    for p in p_list:
        d = RESULTADOS_CSV[(np.isclose(RESULTADOS_CSV.phi, phi_list)) & (np.isclose(RESULTADOS_CSV.p, p))].sort_values("T_r")
        color = plt.gca()._get_lines.get_next_color()
        plt.plot(d.T_r, d.S_L_Cantera/100, ".-", label=fr"Cantera, $p={p/1e6}$ MPa", color=color)
        plt.plot(d.T_r, d.S_L_analitica/100, "--", label=fr"Analítica, $p={p/1e6}$ MPa", color=color)

    plt.xlabel("$T_r$ (K)"); plt.ylabel("$S_L$ (m/s)")
    plt.grid(alpha=.5); plt.legend()

    plt.savefig(f"./Res/S_L/S_vs_T_r_{comburente}.svg")
    plt.show()

# %%
