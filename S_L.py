#%% PREÁMBULO
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(["science"])
import pandas as pd
from scipy.optimize import newton



#%% INPUTS
comburente = "air"  # "oxi" or "air"

# Inputs de phi
phi_list = [0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8]
p_list_plot = [1, 5, 10, 30]       # atm
T_r_list_plot = [298, 400, 500, 600] # K
# phi_list = np.geomspace( 0.5, 1.8, num=25 )
# phi_list_forced = [1.0, 1.2, 1.4]
# phi_list = np.sort(np.unique(np.concatenate(([phi_list, phi_list_forced]))))

# Inputs de p
p_list = [0.1e6, 0.2e6, 0.5e6, 0.8e6, 1.0e6, 2.027e6, 4.0e6, 6.0e6] # Pa
# p_list = np.geomspace( 0.1e6, 6e6, num=20 ) # Pa
# p_list_forced = [0.2e6, 0.5e6, 1e6, 2e6, 20*ct.one_atm]
# p_list = np.sort(np.unique(np.concatenate(([p_list, p_list_forced]))))


# Inputs T_r
T_r_list = [298, 300, 343, 373, 375, 400, 418, 443, 473, 500, 600, 700]
# T_r_list = np.linspace( 298, 700, num=16 ) # K
# T_r_list_forced = [300, 343, 373, 375, 400, 418, 443, 473, 500]
# T_r_list = np.sort(np.unique(np.concatenate(([T_r_list, T_r_list_forced]))))



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
    print( f"\033[1;36m phi = {phi:.6g} | p = {p/ct.one_atm:.6g} atm | T_r = {T_r:.6g} K \033[0m" )

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
        # la key es una tuple (como una lista pero inmutable) con los datos (lists no pueden ser dict keys)
        # el value es un diccionario con sus resultados

        return res, flame.to_array() # Esta flame.to_array se realimenta a la propia función en siguiente bucle como flame_sol_previa

    except Exception as e: # Si falla flame porque no converge, que avise dónde y siga.
        print(f"\033[1;31m FALLO: phi={phi}, p={p/ct.one_atm} atm, T_r={T_r} K \033[0m")
        print(f"    {type(e).__name__}: {e}")

        return {}, None



#%% CÁLCULO DE LLAMA
RESULTADOS = {}
for phi in phi_list:
    for T_r in T_r_list_plot_phi:
        flame_sol_previa=None
        for p in p_list_plot_phi:
            if (phi, p, T_r) not in RESULTADOS: # Calcular solo si no se ha calculado ya para los mismos datos
                res, flame_sol_previa = CALCULO_LLAMA( phi=phi, p=p, T_r=T_r, flame_sol_previa=flame_sol_previa, loglevel=1 )
                RESULTADOS.update(res)

for T_r in T_r_list:
    for phi in phi_list_plot_T_r:
        flame_sol_previa=None
        for p in p_list_plot_T_r:
            if (phi, p, T_r) not in RESULTADOS: # Calcular solo si no se ha calculado ya para los mismos datos
                res, flame_sol_previa = CALCULO_LLAMA( phi=phi, p=p, T_r=T_r, flame_sol_previa=flame_sol_previa, loglevel=1 )
                RESULTADOS.update(res)

for p in p_list:
    for phi in phi_list_plot_p:
        flame_sol_previa=None
        for T_r in T_r_list_plot_p:
            if (phi, p, T_r) not in RESULTADOS: # Calcular solo si no se ha calculado ya para los mismos datos
                res, flame_sol_previa = CALCULO_LLAMA( phi=phi, p=p, T_r=T_r, flame_sol_previa=flame_sol_previa, loglevel=1 )
                RESULTADOS.update(res)



#%% Guardar en CSV
RESULTADOS_CSV = pd.DataFrame([
    {
        "phi": phi,
        "p": p/ct.one_atm,
        "T": T_r,
        "S_L_Cantera": datos["S_L_Cantera"],
        "S_L_analitica": datos["S_L_analitica"]
    }
    for (phi, p, T_r), datos in RESULTADOS.items()
])

RESULTADOS_CSV.to_csv(f"./Res/S_L/RESULTADOS_{comburente}.csv", index=False)



#%% Leer CSV
RESULTADOS_CSV = pd.read_csv(f"./Res/S_L/RESULTADOS_{comburente}.csv")

#%% S vs p | 1 phi, varios T_r
phi_ref = 1.0
T_list = [298, 400, 500, 600]

plt.figure(figsize=(8,8))
for T in T_list:
    d = RESULTADOS_CSV[(np.isclose(RESULTADOS_CSV.phi, phi_ref)) & (np.isclose(RESULTADOS_CSV.T, T))].sort_values("p")
    plt.plot(d.p, d.S_L_Cantera, ".-", label=fr"Cantera, $T_r$={T} K")
    plt.plot(d.p, d.S_L_analitica, "--", label=fr"Analítica, $T_r$={T} K")

plt.xlabel("$p$ [atm]"); plt.ylabel("$S_L$ [cm/s]")
plt.grid(alpha=.5); plt.legend(); plt.show()



#%% S vs phi | 1 T_r, varios p
T_ref = 298
p_list = [1, 5, 10, 30]

plt.figure(figsize=(8,8))
for p in p_list:
    d = RESULTADOS_CSV[(np.isclose(RESULTADOS_CSV.T, T_ref)) & (np.isclose(RESULTADOS_CSV.p, p))].sort_values("phi")
    plt.plot(d.phi, d.S_L_Cantera, ".-", label=fr"Cantera, $p$={p} atm")
    plt.plot(d.phi, d.S_L_analitica, "--", label=fr"Analítica, $p$={p} atm")

plt.xlabel(r"$\phi$"); plt.ylabel("$S_L$ [cm/s]")
plt.grid(alpha=.5); plt.legend(); plt.show()



#%% S vs phi | 1 p, varios T_r
p_ref = 1
T_list = [298, 400, 500, 600]

plt.figure(figsize=(8,8))
for T in T_list:
    d = RESULTADOS_CSV[(np.isclose(RESULTADOS_CSV.p, p_ref)) & (np.isclose(RESULTADOS_CSV.T, T))].sort_values("phi")
    plt.plot(d.phi, d.S_L_Cantera, ".-", label=fr"Cantera, $T_r$={T} K")
    plt.plot(d.phi, d.S_L_analitica, "--", label=fr"Analítica, $T_r$={T} K")

plt.xlabel(r"$\phi$"); plt.ylabel("$S_L$ [cm/s]")
plt.grid(alpha=.5); plt.legend(); plt.show()



#%% S vs T_r | 1 phi, varios p
phi_ref = 1.0
p_list = [1, 5, 10, 30]

plt.figure(figsize=(8,8))
for p in p_list:
    d = RESULTADOS_CSV[(np.isclose(RESULTADOS_CSV.phi, phi_ref)) & (np.isclose(RESULTADOS_CSV.p, p))].sort_values("T")
    plt.plot(d.T, d.S_L_Cantera, ".-", label=fr"Cantera, $p$={p} atm")
    plt.plot(d.T, d.S_L_analitica, "--", label=fr"Analítica, $p$={p} atm")

plt.xlabel("$T_r$ [K]"); plt.ylabel("$S_L$ [cm/s]")
plt.grid(alpha=.5); plt.legend(); plt.show()
