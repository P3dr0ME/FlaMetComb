#%% PREÁMBULO
# https://www.cantera.org/3.1/userguide/flame-temperature.html

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
import pandas as pd
from scipy.optimize import newton
from scipy.integrate import quad
plt.style.use(['science'])


#%% INPUTS
comburente = "oxi" # "oxi" or "air"

T_r = 298 # K
p = ct.one_atm

N = 25 # Número de puntos en el vector phi_list
phi_start = 0.2
phi_stop = 1.8


#Crear vector con valores de ratio de equivalencia.
phi_list = np.linspace(phi_start, phi_stop, N)
if not 1.00 in phi_list: # añadir phi=1.00 al la fuerza
    idx = np.searchsorted(phi_list, 1.00)
    phi_list = np.insert(phi_list, idx, 1.00)

#Definir el oxidante
oxidizer = "O2" if comburente == "oxi" else f"O2:1, N2:{79/21}"

# Especies
species_dict = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}
ideal_species = [species_dict[S] for S in ("CH4", "O2", "N2", "CO2", "H2O")]

# Gases
gas_ideal = ct.Solution(thermo="ideal-gas",
                    species=ideal_species,
                    transport_model='mixture-averaged',
                    kinetics='gas')
gas_real = ct.Solution('gri30.yaml')

# Para el Método 1 y 2
M_CH4 = 16.04 / 1e3 # kg/mol
M_O2 = 32.00 / 1e3 # kg/mol
M_N2 = 28.014 / 1e3 # kg/mol
M_aire = M_O2+79/21*M_N2 # kg/mol
f_s = 1*M_CH4/(2*M_aire) if comburente == "air" else 1*M_CH4/(2*M_O2) # Dosado estequiométrico


#%% CÁLCULO DE T_Ad Método 1 (Analítico con cp promedio)
LHV = 50.048 * 1e6 # J/kg. Springer, Appendix 1.

# Inicializar diccionarios
T_ad_m1, cp_m1, q0_p_m1 = (
    {phi: None for phi in phi_list} for _ in range(3)
)

for phi in phi_list:
    q0_p_m1[phi] = ( phi*f_s*LHV ) / (1+phi*f_s) / 1e6 if phi<=1 else ( f_s*LHV ) / (1+phi*f_s) / 1e6 # MJ/kg

    def f1(T_p):
        gas_ideal.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
        gas_ideal.TP = (T_p+T_r)/2, p # T~=(T_ad+T_r)/2
        return (T_p - T_r)*gas_ideal.cp - q0_p_m1[phi]*1e6

    T_ad_m1[phi] = newton(f1, x0=1000, tol=1e-12)
    cp_m1[phi] = (q0_p_m1[phi]*1e6) / (T_ad_m1[phi] - T_r)

#%% CÁLCULO DE T_Ad Método 2 (Analítico H_p=H_r)
Dh0_f = {"CH4": -74.87e3, "O2": 0, "N2": 0, "CO2": -393.52e3, "H2O": -241.83e3} #J/mol
# Springer, Table 2.2.

# Excesos de combustible y de comburente en función de phi
eps = { phi:  0 if phi<=1 else phi-1  for phi in phi_list }
delta = { phi:  2/phi - 2 if phi<=1 else 0  for phi in phi_list }

# Moles de especies en función de phi
n_ir = {
    phi:  {"CH4": 1+eps[phi], "O2": 2+delta[phi], "N2": 79/21*(2+delta[phi]) if comburente=="air" else 0 , "CO2": 0, "H2O": 0}
    for phi in phi_list
    }
n_ip = {
    phi:  {"CH4": eps[phi], "O2": delta[phi], "N2": 79/21*(2+delta[phi]) if comburente=="air" else 0, "CO2": 1, "H2O": 2}
    for phi in phi_list
    }

# c_p^_i(T), según polinomios NASA7/9 de gri30.yaml
def cp_molar(sp, T):
    return gas_ideal.species(f"{sp}").thermo.cp(T) / 1e3 # J/mol/K

# h^_s,i(T)
def h_s_molar(sp, T0, T):
    return quad( lambda T_: cp_molar(sp, T_), T0, T )[0]  # J/mol
    # quad necesita una función de 1 sola var para integrar
    # con lambda x: f1(x) defino esa función anónimamente
    # quad saca lista [valor_integral, error]

# Inicializar diccionarios
T_ad_m2, q0_p_m2 , cp_m2_ave = (
    {phi: None for phi in phi_list} for _ in range(3)
)

for phi in phi_list:
    # Calor de combustión
    Q0_p = (
        sum(n_ip[phi][sp] * Dh0_f[sp] for sp in Dh0_f)
      - sum(n_ir[phi][sp] * Dh0_f[sp] for sp in Dh0_f)
    )

    # Función cuya raíz se busca
    def f2(T_p):
        H_p = sum( n_ip[phi][sp] * h_s_molar(sp, 298.15, T_p) for sp in Dh0_f )
        H_r = sum( n_ir[phi][sp] * h_s_molar(sp, 298.15, T_r) for sp in Dh0_f )

        return H_p - (-Q0_p + H_r)

    # Resolver f(Tad) = 0
    T_ad_m2[phi] = newton(f2, x0=1000, tol=1e-12)

    q0_p_m2[phi] = -Q0_p / ( (n_ir[phi]["CH4"]*M_CH4) * (1+1/(phi*f_s)) )  / 1e6 # MJ/kg
    cp_m2_ave[phi] = (q0_p_m2[phi]*1e6) / (T_ad_m2[phi] - T_r)
    # cp,p(T~) según fórmula Método 1 para poder comparar




#%% CÁLCULO DE T_ad Método 3 (Cantera)
# Inicializar diccionarios
T_ad_ideal, T_ad_real, cp_ideal, cp_real, q0_p_ideal, q0_p_real, h_r_T_r_real, h_r_T_r_ideal, h_p_T_r_real, h_p_T_r_ideal, q0_p_real, q0_p_ideal = ( {phi: None for phi in phi_list} for _ in range(12) )
Y_r_real, Y_p_real, Y_r_ideal, Y_p_ideal, = ( {phi: {} for phi in phi_list} for _ in range(4) )

for phi in phi_list:
    # IDEAL:
    gas_ideal.TP = T_r, p
    gas_ideal.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")

    Y_r_ideal[phi] = gas_ideal.mass_fraction_dict()
    h_r_T_r_ideal[phi] = gas_ideal.enthalpy_mass # J/kg

    gas_ideal.equilibrate("HP")
        # equilibrate() calcula el estado de equilibrio, que minimiza el potencial de Gibbs. Como esta combustión es espontánea,
        # ese estado final es el posterior a la combustión y como hemos impuesto H cte., su T es la T_ad.

    Y_p_ideal[phi] = gas_ideal.mass_fraction_dict()
    # cp_ideal[phi] = gas_ideal.cp

    T_ad_ideal[phi] = gas_ideal.T

    gas_ideal.TP = T_r, p
    # Ahora el gas tiene las Y de productos, pero a T_r y p. Por eso se puede calcular h_p(T_r).
    h_p_T_r_ideal[phi] = gas_ideal.enthalpy_mass # J/kg
    # Calor de combustión por kg de mezcla:
    q0_p_ideal[phi] = -(h_p_T_r_ideal[phi] - h_r_T_r_ideal[phi]) / 1e6 # MJ/kg.

    cp_ideal[phi] = (q0_p_ideal[phi]*1e6) / (T_ad_ideal[phi] - T_r)
    # cp,p(T~) según fórmula Método 1 para poder comparar

    # REAL:
    gas_real.TP = T_r, p
    # Se restablece T y p iniciales en cada bucle para calcular la T_ad
    gas_real.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")

    Y_r_real[phi] = gas_real.mass_fraction_dict()
    h_r_T_r_real[phi] = gas_real.enthalpy_mass # J/kg

    gas_real.equilibrate("HP")

    Y_p_real[phi] = gas_real.mass_fraction_dict()

    T_ad_real[phi] = gas_real.T

    gas_real.TP = T_r, p
    # Ahora el gas tiene las Y de productos, pero a T_r y p. Por eso se puede calcular h_p(T_r).
    h_p_T_r_real[phi] = gas_real.enthalpy_mass # J/kg
    # Calor de combustión por kg de mezcla:
    q0_p_real[phi] = -(h_p_T_r_real[phi] - h_r_T_r_real[phi]) / 1e6 # MJ/kg.

    cp_real[phi] = (q0_p_real[phi]*1e6) / (T_ad_real[phi] - T_r)
    # cp,p(T~) según fórmula Método 1 para poder comparar


#%% TABLA COMPRACIÓN q0_p
def M(sp):
    return gas_real.species(sp).molecular_weight / 1e3 # kg/mol

def Dh0_f_Cantera(sp):
    return gas_real.species(sp).thermo.h(298.15) /1e3 # J/mol
    # CHECK: muy parecidas a las que da Springer, Table 2.2.

# Diferencia de calor de combustión por masa de mezcla
Dq0_p = { phi: (q0_p_ideal[phi] - q0_p_real[phi]) for phi in phi_list } # MJ/kg

# Calor por masa de mezcla que se pierde en formación de productos que real sí tiene pero ideal no.
q0_p_dis = {phi: (
       + sum( (Y_p_real[phi][sp]/M(sp) * Dh0_f_Cantera(sp)) for sp in Y_p_real[phi] )
       - sum( (Y_p_ideal[phi][sp]/M(sp) * Dh0_f_Cantera(sp)) for sp in Y_p_ideal[phi] )
    )
    / 1e6 # MJ/kg
    for phi in phi_list}

Tabla_q0_p = pd.DataFrame({
    "q0_p,ideal - q0_p,real (MJ/kg)": Dq0_p,
    "q0_p disociación (MJ/kg)":  q0_p_dis,
}, index=list(phi_list))

Tabla_q0_p.index.name = "phi" # Cabecero columna primera

print(Tabla_q0_p)
Tabla_q0_p.to_csv(f"./results/T_ad/Tabla_q0_p_{comburente}.csv")


#%% TABLA cp
Tabla_c_p = pd.DataFrame({
    "cp_ideal": cp_ideal,
    "cp_real": cp_real,
    "cp_m1": cp_m1,
    "cp_m2_ave": cp_m2_ave
}, index=list(phi_list))

Tabla_c_p.index.name = "phi" # Cabecero columna primera

print(Tabla_c_p)
Tabla_c_p.to_csv(f"./results/T_ad/Tabla_c_p_{comburente}.csv")

#%% GRÁFICO T - phi
plt.figure(figsize=(8,8))

plt.title(
    #r"\bf{Temperatura\ adiabática\ de\ llama\ frente\ a\ ratio\ de\ equivalencia}" + "\n"
    f"{r"\bf{Oxígeno}" if comburente == 'oxi' else r"\bf{Aire}"} \n"
    fr"$p = {p}$ Pa $\quad T_{{r}} = {T_r}$ K",
    fontsize=11,
    pad=15
)

plt.plot(phi_list, T_ad_real.values(), '.-', label="Método 3 (Cantera - real)")
plt.plot(phi_list, T_ad_ideal.values(), '.-', label="Método 3 (Cantera - ideal)")
plt.plot(phi_list, T_ad_m1.values(), '', label=r"Método 1 ($c_p$ promedio)")
plt.plot(phi_list, T_ad_m2.values(), '', label="Método 2 (balance de entalpías)")

plt.grid(True, which='both', alpha=0.5)

plt.xlabel(r"$\phi$")
plt.ylabel(r"$T_{ad}$ (K)")
plt.legend()

plt.xlim(phi_list[0], phi_list[-1])

plt.savefig(f"./results/T_ad/T_ad_vs_phi_{comburente}.svg")
plt.show() # muestra el gráfico y debe ir después de plt.savefig()



#%% GRÁFICO DE q0_p vs phi
plt.figure(figsize=(8,8))

plt.title(
    #r"\bf{Calor\ de\ combustión\ por\ kg\ de\ mezcla\ frente\ a\ ratio\ de\ equivalencia}" + "\n"
    f"{r"\bf{Oxígeno}" if comburente == 'oxi' else r"\bf{Aire}"} \n"
    fr"$p = {p}$ Pa $\quad T_{{0}} = {T_r}$ K",
    fontsize=11,
    pad=15
)

plt.plot(phi_list, q0_p_real.values(), '.-', label="Método 3 (Cantera - real)")
plt.plot(phi_list, q0_p_ideal.values(), '.-', label="Método 3 (Cantera - ideal)")
plt.plot(phi_list, q0_p_m1.values(), '', label=r"Método 1 ($c_p$ promedio)")
plt.plot(phi_list, q0_p_m2.values(), '', label="Método 2 (balance de entalpías)")

plt.grid(True, which='both', alpha=0.5)

plt.xlabel(r"$\phi$")
plt.ylabel(r"$q^0_p$ (MJ/kg)")
plt.legend()

plt.xlim(phi_list[0], phi_list[-1])

plt.savefig(f"./results/T_ad/q0_p_vs_phi_{comburente}.svg")
plt.show()

# %%
