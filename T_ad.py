#%% PREÁMBULO
# https://www.cantera.org/3.1/userguide/flame-temperature.html
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(["science"])
import pandas as pd
from scipy.optimize import newton
from scipy.integrate import quad
import os



#%% INPUTS
comburente = "air"  # "oxi" or "air"
borrar_csv_antiguo = False


# Seleccionar curva a representar (descomentar 1)
# "1 {dato}" --> dato_list debe contener 1 solo elemento entre corchetes []

# # T_ad vs phi | 1 T_r, varios p
# selector_plot = "T_ad vs phi | 1 T_r, varios p"
# phi_list = np.sort(np.unique(np.concatenate([np.linspace(0.2, 1.8, num=25), [0.6, 1.00]])))
# p_list = [ct.one_atm]  # Pa
# T_r_list = [298]  # K

# # T_ad vs phi | 1 p, varios T_r
# selector_plot = "T_ad vs phi | 1 p, varios T_r"
# phi_list =
# p_list = # Pa
# T_r_list = # K

# T_ad vs p | 1 phi, varios T_r
selector_plot = "T_ad vs p | 1 phi, varios T_r"
phi_list = [1.00]
p_list = np.linspace(0.1*ct.one_atm, 30*ct.one_atm, num=25) # Pa
T_r_list = [300] # K

# # T_ad vs T_r | 1 phi, varios p
# selector_plot = "T_ad vs T_r | 1 phi, varios p"
# phi_list = [1.00]
# p_list = [ct.one_atm] # Pa
# T_r_list = np.linspace(300, 750, num=10) # K



#%% ESPECIES Y GASES
# Oxidante
oxidizer = "O2" if comburente == "oxi" else f"O2:1, N2:{79/21}"

# Especies (para gas ideal reducido)
species_dict = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}

list_ideal_species = "CH4", "O2", "N2", "CO2", "H2O"
ideal_species = [species_dict[S] for S in list_ideal_species]

# Gases
gas_ideal = ct.Solution(thermo="ideal-gas",
                    species=ideal_species,
                    transport_model="mixture-averaged",
                    kinetics="gas")
gas_real = ct.Solution("gri30.yaml")

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

# Entalpías de formación estándar (Método 2)
Dh0_f = {"CH4": -74.87e3, "O2": 0, "N2": 0, "CO2": -393.52e3, "H2O": -241.83e3} # J/mol
# Springer, Table 2.2.



#%% FUNCIONES AUXILIARES
def M(sp):
    return gas_real.species(sp).molecular_weight / 1e3 # kg/mol

# cp^_i(T) según polinomios NASA7/9 de gri30.yaml
def cp_molar(sp, T):
    return gas_real.species(f"{sp}").thermo.cp(T) / 1e3 # J/mol/K

# h^_s,i(T0->T)
def h_s_molar(sp, T0, T):
    return quad( lambda T_: cp_molar(sp, T_), T0, T )[0] # J/mol

def Dh0_f_Cantera(sp):
    return gas_real.species(sp).thermo.h(298.15) / 1e3 # J/mol

def moles_reactantes_productos(phi):
    eps = 0 if phi<=1 else phi-1
    delta = 2/phi-2 if phi<=1 else 0
    n_ir = {"CH4": 1+eps, "O2": 2+delta, "N2": 79/21*(2+delta) if comburente=="air" else 0, "CO2": 0, "H2O": 0}
    n_ip = {"CH4": eps, "O2": delta, "N2": 79/21*(2+delta) if comburente=="air" else 0, "CO2": 1, "H2O": 2}
    m_mezcla = n_ir["CH4"] * M_CH4 * (1 + 1/(phi*f_s))
    return n_ir, n_ip, m_mezcla



#%% FUNCIÓN: T_ad MÉTODO 1 (analítico, cp promedio)
def Metodo1_T_ad( phi, T_r ):
    n_ir, n_ip, m_mezcla = moles_reactantes_productos(phi)
    q0_p_m1 = ( (phi*f_s*LHV)/(1+phi*f_s)/1e6 if phi<=1 else (f_s*LHV)/(1+phi*f_s)/1e6 ) # MJ/kg

    def f1(T_p):
        cp = sum( n_ip[sp]*cp_molar(sp, (T_p+T_r)/2) for sp in list_ideal_species ) / m_mezcla
        return (T_p - T_r)*cp - q0_p_m1*1e6

    T_ad_m1 = newton(f1, x0=1000, tol=1e-12)
    return T_ad_m1, q0_p_m1



#%% FUNCIÓN: T_ad MÉTODO 2 (analítico, balance de entalpías H_p = H_r)
def Metodo2_T_ad( phi, T_r ):
    n_ir, n_ip, m_mezcla = moles_reactantes_productos(phi)

    Q0_p = ( sum(n_ip[sp]*Dh0_f[sp] for sp in list_ideal_species)
           - sum(n_ir[sp]*Dh0_f[sp] for sp in list_ideal_species) )

    def f2(T_p):
        H_p = sum( n_ip[sp] * h_s_molar(sp, 298.15, T_p) for sp in list_ideal_species )
        H_r = sum( n_ir[sp] * h_s_molar(sp, 298.15, T_r) for sp in list_ideal_species )
        return H_p - (-Q0_p + H_r)

    T_ad_m2 = newton(f2, x0=1000, tol=1e-12)
    q0_p_m2 = -Q0_p / m_mezcla / 1e6 # MJ/kg
    return T_ad_m2, q0_p_m2



#%% FUNCIÓN: T_ad MÉTODO 3 (Cantera, equilibrio HP)
def Metodo3_T_ad( phi, p, T_r, gas ):
    gas.TP = T_r, p
    gas.set_equivalence_ratio(phi, "CH4", oxidizer)

    Y_r = gas.mass_fraction_dict()
    h_r_T_r = gas.enthalpy_mass # J/kg

    gas.equilibrate("HP")
    # equilibrate() calcula el estado que minimiza el potencial de Gibbs; como
    # se ha impuesto H cte., su T es la T_ad.

    Y_p = gas.mass_fraction_dict()
    T_ad_m3 = gas.T

    gas.TP = T_r, p
    # El gas conserva las Y de productos pero se lleva a T_r, p para poder calcular h_p(T_r).
    h_p_T_r = gas.enthalpy_mass # J/kg
    q0_p_m3 = -(h_p_T_r - h_r_T_r) / 1e6 # MJ/kg

    cp_p = sum( Y_p[sp]/M(sp) * cp_molar(sp, (T_ad_m3+T_r)/2) for sp in Y_p )

    return {"T_ad_m3": T_ad_m3, "q0_p_m3": q0_p_m3, "cp_p": cp_p, "Y_p": Y_p}



#%% FUNCIÓN DE CÁLCULO CONJUNTO
def CALCULO_T_AD( phi, p, T_r ):
    print( f"\033[1;36m phi = {phi} | p = {p/1e6} MPa | T_r = {T_r:.6g} K \033[0m" )

    try:
        T_ad_m1, q0_p_m1 = Metodo1_T_ad( phi=phi, T_r=T_r )
        T_ad_m2, q0_p_m2 = Metodo2_T_ad( phi=phi, T_r=T_r )
        res_real = Metodo3_T_ad( phi=phi, p=p, T_r=T_r, gas=gas_real )
        res_ideal = Metodo3_T_ad( phi=phi, p=p, T_r=T_r, gas=gas_ideal )

        # Comparación del calor de combustión: pérdida por disociación (real vs ideal)
        Dq0_p = (res_ideal["q0_p_m3"] - res_real["q0_p_m3"]) * 1e6 # J/kg
        q0_p_dis = ( sum( res_real["Y_p"][sp]/M(sp) * Dh0_f_Cantera(sp) for sp in res_real["Y_p"] )
                   - sum( res_ideal["Y_p"][sp]/M(sp) * Dh0_f_Cantera(sp) for sp in res_ideal["Y_p"] ) ) # J/kg
        Error_rel = (Dq0_p - q0_p_dis) / Dq0_p * 100 # %

        res = {
            (phi, p, T_r):
            {
                "T_ad_real": res_real["T_ad_m3"],
                "T_ad_ideal": res_ideal["T_ad_m3"],
                "T_ad_m1": T_ad_m1,
                "T_ad_m2": T_ad_m2,
                "q0_p_real": res_real["q0_p_m3"],
                "q0_p_ideal": res_ideal["q0_p_m3"],
                "q0_p_m1": q0_p_m1,
                "q0_p_m2": q0_p_m2,
                "cp_p_real": res_real["cp_p"],
                "cp_p_ideal": res_ideal["cp_p"],
                "Dq0_p_Jkg": Dq0_p,
                "q0_p_dis_Jkg": q0_p_dis,
                "Error_rel_pct": Error_rel,
            }
        }
        # key es una tuple con los datos (lists no pueden ser dict keys)
        # value es un diccionario con resultados

        return res

    except Exception as e: # Si falla la resolución, que avise dónde y siga.
        print(f"\033[1;31m FALLO: phi={phi}, p={p/1e6} MPa, T_r={T_r} K \033[0m")
        print(f"    {type(e).__name__}: {e}")

        return {}



#%% CÁLCULO DE T_ad
RESULTADOS = {}
for phi in phi_list:
    for T_r in T_r_list:
        for p in p_list:
            if (phi, p, T_r) not in RESULTADOS: # Calcular solo si no se ha calculado ya para los mismos datos
                res = CALCULO_T_AD( phi=phi, p=p, T_r=T_r )
                RESULTADOS.update(res)



#%% TABLA CSV
ruta = f"./Res/T_ad/T_ad_RESULTADOS_{comburente}.csv"

# Borrar CSV antiguo (si existe)
if borrar_csv_antiguo:
    os.path.exists(ruta) and os.remove(ruta)

# Resultados nuevos
RESULTADOS_CSV_NUEVOS = pd.DataFrame([
    {
        "phi": phi,
        "p": p,
        "T_r": T_r,
        **datos
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
os.makedirs(os.path.dirname(ruta), exist_ok=True)
RESULTADOS_CSV.to_csv(ruta, index=False)



#%% CONFIGURACIÓN DE CADA SELECTOR DE GRÁFICO
# x: variable en eje X | grupo: variable con la que se generan varias curvas | fijo: variable(s) que se mantienen constantes (se toma su primer valor de la lista)
CONFIG_SELECTOR = {
    "T_ad vs phi | 1 T_r, varios p": dict(x="phi", x_label=r"$\phi$", grupo="p", grupo_list=p_list, fijo={"T_r": T_r_list[0]}),
    "T_ad vs phi | 1 p, varios T_r": dict(x="phi", x_label=r"$\phi$", grupo="T_r", grupo_list=T_r_list, fijo={"p": p_list[0]}),
    "T_ad vs p | 1 phi, varios T_r":  dict(x="p",   x_label="$p$ (MPa)", grupo="T_r", grupo_list=T_r_list, fijo={"phi": phi_list[0]}),
    "T_ad vs T_r | 1 phi, varios p":  dict(x="T_r", x_label="$T_r$ (K)", grupo="p",  grupo_list=p_list,  fijo={"phi": phi_list[0]}),
}

def _leyenda(grupo, val):
    if grupo == "p":
        return fr"$p={val/1e6}$ MPa"
    if grupo == "T_r":
        return fr"$T_r={val}$ K"
    return f"{grupo}={val}"

def _subtitulo(fijo):
    partes = []
    for k, v in fijo.items():
        if k == "phi":
            partes.append(fr"$\phi = {v}$")
        elif k == "p":
            partes.append(fr"$p = {v/1e6}$ MPa")
        elif k == "T_r":
            partes.append(fr"$T_r = {v}$ K")
    return "  ".join(partes)

def PLOT_GENERICO( columnas, ylabel, nombre_archivo, escala_x=1.0 ):
    """
    columnas: lista de tuplas (etiqueta_metodo, nombre_columna, estilo_linea)
    """
    cfg = CONFIG_SELECTOR[selector_plot]
    x_col, grupo_col, grupo_list, fijo = cfg["x"], cfg["grupo"], cfg["grupo_list"], cfg["fijo"]

    d_fijo = RESULTADOS_CSV.copy()
    for k, v in fijo.items():
        d_fijo = d_fijo[np.isclose(d_fijo[k], v)]

    plt.figure(figsize=(8,6))
    plt.title(
        f"{r'\bf{Oxígeno}' if comburente == 'oxi' else r'\bf{Aire}'} \n"
        f"{_subtitulo(fijo)}",
        fontsize=11,
        pad=15
    )

    for val in grupo_list:
        d = d_fijo[np.isclose(d_fijo[grupo_col], val)].sort_values(x_col)
        for etiqueta, col, estilo in columnas:
            plt.plot(d[x_col]/escala_x, d[col], estilo, label=f"{etiqueta}, {_leyenda(grupo_col, val)}")

    plt.xlabel(cfg["x_label"]); plt.ylabel(ylabel)
    plt.grid(True, which='both', alpha=0.5); plt.legend()

    plt.savefig(f"./Res/T_ad/{nombre_archivo}_{comburente}.svg")
    plt.show()


# GRÁFICO T_ad
PLOT_GENERICO(
    columnas=[
        ("Método 3 (Cantera - real)", "T_ad_real", ".-"),
        ("Método 3 (Cantera - ideal)", "T_ad_ideal", ".-"),
        ("Método 1 ($c_p$ promedio)", "T_ad_m1", "-"),
        ("Método 2 (balance de entalpías)", "T_ad_m2", "-"),
    ],
    ylabel=r"$T_{ad}$ (K)",
    nombre_archivo="T_ad_vs_" + CONFIG_SELECTOR[selector_plot]["x"],
    escala_x=(1e6 if CONFIG_SELECTOR[selector_plot]["x"] == "p" else 1.0),
)


# GRÁFICO q0_p
PLOT_GENERICO(
    columnas=[
        ("Método 3 (Cantera - real)", "q0_p_real", ".-"),
        ("Método 3 (Cantera - ideal)", "q0_p_ideal", ".-"),
        ("Método 1 ($c_p$ promedio)", "q0_p_m1", "-"),
        ("Método 2 (balance de entalpías)", "q0_p_m2", "-"),
    ],
    ylabel=r"$q^0_p$ (MJ/kg)",
    nombre_archivo="q0_p_vs_" + CONFIG_SELECTOR[selector_plot]["x"],
    escala_x=(1e6 if CONFIG_SELECTOR[selector_plot]["x"] == "p" else 1.0),
)


# GRÁFICO cp productos promedio
PLOT_GENERICO(
    columnas=[
        ("Método 3 (Cantera - real)", "cp_p_real", ".-"),
        ("Método 3 (Cantera - ideal)", "cp_p_ideal", ".-"),
    ],
    ylabel=r"$\tilde{c_p}_p$ (J/(kg K))",
    nombre_archivo="cp_p_vs_" + CONFIG_SELECTOR[selector_plot]["x"],
    escala_x=(1e6 if CONFIG_SELECTOR[selector_plot]["x"] == "p" else 1.0),
)

# %%
