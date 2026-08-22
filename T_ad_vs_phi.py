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

phi_forced_list = [0.6, 1.00] # valores que añadir a la fuerza.



#Crear vector con valores de ratio de equivalencia.
phi_list = np.sort( np.unique( np.concatenate([np.linspace(phi_start, phi_stop, N), phi_forced_list]) ) )

#Definir el oxidante
oxidizer = "O2" if comburente == "oxi" else f"O2:1, N2:{79/21}"

# Especies
species_dict = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}

list_ideal_species = "CH4", "O2", "N2", "CO2", "H2O"
ideal_species = [species_dict[S] for S in list_ideal_species]

# Gases
gas_ideal = ct.Solution(thermo="ideal-gas",
                    species=ideal_species,
                    transport_model='mixture-averaged',
                    kinetics='gas')
gas_real = ct.Solution('gri30.yaml')

#%% Para el Método 1 y 2
M_CH4 = 16.04 / 1e3 # kg/mol
M_O2 = 32.00 / 1e3 # kg/mol
M_N2 = 28.014 / 1e3 # kg/mol
M_aire = M_O2+79/21*M_N2 # kg/mol
f_s = 1*M_CH4/(2*M_aire) if comburente == "air" else 1*M_CH4/(2*M_O2) # Dosado estequiométrico

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

m_mezcla = {phi: n_ir[phi]["CH4"]*M_CH4*(1+1/(phi*f_s)) for phi in phi_list }

#%% Para Método 2 y 3
# Masa molar
def M(sp):
    return gas_real.species(sp).molecular_weight / 1e3 # kg/mol

# cp^_i(T) según polinomios NASA7/9 de gri30.yaml
def cp_molar(sp, T):
    return gas_real.species(f"{sp}").thermo.cp(T) / 1e3 # J/mol/K



#%% CÁLCULO DE T_Ad Método 1 (Analítico con cp promedio)
LHV = 50.048 * 1e6 # J/kg. Springer, Appendix 1.

# Inicializar diccionarios
T_ad_m1, cp_p_m1, q0_p_m1 = (
    {phi: None for phi in phi_list} for _ in range(3)
)

for phi in phi_list:
    q0_p_m1[phi] = (phi*f_s*LHV) / (1+phi*f_s) / 1e6 if phi<=1 else (f_s*LHV) / (1+phi*f_s) / 1e6 # MJ/kg

    def f1(T_p):
        # cp de productos promedio
        cp = sum( n_ip[phi][sp]*cp_molar(sp, (T_p+T_r)/2) for sp in list_ideal_species) / m_mezcla[phi]
        return (T_p - T_r)*cp - q0_p_m1[phi]*1e6

    T_ad_m1[phi] = newton(f1, x0=1000, tol=1e-12)

    # cp de productos promedio (en m1, def = equivM1)
    # cp_p_m1[phi] = sum( n_ip[phi][sp]*cp_molar(sp, (T_ad_m1[phi]+T_r)/2) for sp in list_ideal_species) / m_mezcla[phi]



#%% CÁLCULO DE T_Ad Método 2 (Analítico H_p=H_r)
Dh0_f = {"CH4": -74.87e3, "O2": 0, "N2": 0, "CO2": -393.52e3, "H2O": -241.83e3} # J/mol
# Springer, Table 2.2.

# h^_s,i(T)
def h_s_molar(sp, T0, T):
    return quad( lambda T_: cp_molar(sp, T_), T0, T )[0]  # J/mol
    # quad necesita una función de 1 sola var para integrar
    # con lambda x: f1(x) defino esa función anónimamente
    # quad saca lista [valor_integral, error]

# Inicializar diccionarios
T_ad_m2, q0_p_m2 , cp_p_m2_def, cp_p_m2_equivM1 = (
    {phi: None for phi in phi_list} for _ in range(4)
)

for phi in phi_list:
    # Calor de combustión
    Q0_p = (
        sum(n_ip[phi][sp] * Dh0_f[sp] for sp in list_ideal_species)
      - sum(n_ir[phi][sp] * Dh0_f[sp] for sp in list_ideal_species)
    )

    # Función cuya raíz se busca
    def f2(T_p):
        H_p = sum( n_ip[phi][sp] * h_s_molar(sp, 298.15, T_p) for sp in list_ideal_species )
        H_r = sum( n_ir[phi][sp] * h_s_molar(sp, 298.15, T_r) for sp in list_ideal_species )

        return H_p - (-Q0_p + H_r)

    # Resolver f(Tad) = 0
    T_ad_m2[phi] = newton(f2, x0=1000, tol=1e-12)

    q0_p_m2[phi] = -Q0_p / m_mezcla[phi] / 1e6 # MJ/kg

    # cp de productos promedio
    # cp_p_m2_def[phi] = sum( n_ip[phi][sp]*cp_molar(sp, (T_ad_m2[phi]+T_r)/2) for sp in list_ideal_species ) / m_mezcla[phi]
    # cp_p_m2_equivM1[phi] = q0_p_m2[phi]*1e6 /( T_ad_m2[phi]-T_r )


#%% CÁLCULO DE T_ad Método 3 (Cantera)
# Inicializar diccionarios
T_ad_ideal, T_ad_real, cp_p_ideal_def, cp_p_real_def, cp_p_ideal_equivM1, cp_p_real_equivM1, q0_p_ideal, q0_p_real, h_r_T_r_real, h_r_T_r_ideal, h_p_T_r_real, h_p_T_r_ideal, q0_p_real, q0_p_ideal = ( {phi: None for phi in phi_list} for _ in range(14) )
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

    T_ad_ideal[phi] = gas_ideal.T

    gas_ideal.TP = T_r, p
    # Ahora el gas tiene las Y de productos, pero a T_r y p. Por eso se puede calcular h_p(T_r).
    h_p_T_r_ideal[phi] = gas_ideal.enthalpy_mass # J/kg
    # Calor de combustión por kg de mezcla:
    q0_p_ideal[phi] = -(h_p_T_r_ideal[phi] - h_r_T_r_ideal[phi]) / 1e6 # MJ/kg.

    # cp de productos promedio
    cp_p_ideal_def[phi] = sum( Y_p_ideal[phi][sp]/M(sp)*cp_molar(sp, (T_ad_ideal[phi]+T_r)/2) for sp in Y_p_ideal[phi] )
    # cp_p_ideal_equivM1[phi] = q0_p_ideal[phi]*1e6 / ( T_ad_ideal[phi]-T_r )

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

    # cp de productos promedio
    cp_p_real_def[phi] = sum( Y_p_real[phi][sp]/M(sp)*cp_molar(sp, (T_ad_real[phi]+T_r)/2) for sp in Y_p_real[phi] )
    # cp_p_real_equivM1[phi] = q0_p_real[phi]*1e6 / ( T_ad_real[phi]-T_r )


#%% TABLA COMPARACIÓN q0_p
def Dh0_f_Cantera(sp):
    return gas_real.species(sp).thermo.h(298.15) /1e3 # J/mol
    # CHECK: muy parecidas a las que da Springer, Table 2.2.

# Diferencia de calor de combustión por masa de mezcla
Dq0_p = { phi: (q0_p_ideal[phi] - q0_p_real[phi])* 1e6 for phi in phi_list } # J/kg

# Calor por masa de mezcla que se pierde en formación de productos que real sí tiene pero ideal no.
q0_p_dis = {phi: (
       + sum( (Y_p_real[phi][sp]/M(sp) * Dh0_f_Cantera(sp)) for sp in Y_p_real[phi] )
       - sum( (Y_p_ideal[phi][sp]/M(sp) * Dh0_f_Cantera(sp)) for sp in Y_p_ideal[phi] )
    ) # J/kg
    for phi in phi_list}

Error_rel = {phi: (Dq0_p[phi] - q0_p_dis[phi]) / Dq0_p[phi] * 100  for phi in phi_list} # %

Tabla_Dif_q0_p = pd.DataFrame({
    "Dq0_p (J/kg)": Dq0_p,
    "q0_p disociación (J/kg)":  q0_p_dis,
    "Error relativo (porcentaje)": Error_rel
}, index=list(phi_list))

Tabla_Dif_q0_p.index.name = "phi" # Cabecero columna primera

print(Tabla_Dif_q0_p)
Tabla_Dif_q0_p.to_csv(f"./Res/T_ad/Tabla_Dif_q0_p_{comburente}.csv")



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

plt.savefig(f"./Res/T_ad/T_ad_vs_phi_{comburente}.svg")
plt.show(block=False) # muestra el gráfico y debe ir después de plt.savefig()



#%% GRÁFICO DE q0_p vs phi
plt.figure(figsize=(8,8))

plt.title(
    #r"\bf{Calor\ de\ combustión\ por\ kg\ de\ mezcla\ frente\ a\ ratio\ de\ equivalencia}" + "\n"
    f"{r"\bf{Oxígeno}" if comburente == 'oxi' else r"\bf{Aire}"} \n"
    fr"$p = {p}$ Pa $\quad T_r = {T_r}$ K",
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

plt.savefig(f"./Res/T_ad/q0_p_vs_phi_{comburente}.svg")
plt.show(block=False)



#%% GRÁFICO DE cp productos promedio vs phi
plt.figure(figsize=(8,8))

plt.title(
    #r"\bf{Calor\ específico\ de\ productos\ promedio\ frente\ a\ ratio\ de\ equivalencia}" + "\n"
    f"{r"\bf{Oxígeno}" if comburente == 'oxi' else r"\bf{Aire}"} \n"
    fr"$p = {p}$ Pa $\quad T_r = {T_r}$ K",
    fontsize=11,
    pad=15
)

plt.plot(phi_list, cp_p_real_def.values(), '.-', label="Método 3 (Cantera - real)")
plt.plot(phi_list, cp_p_ideal_def.values(), '.-', label="Método 3 (Cantera - ideal)")
# plt.plot(phi_list, cp_p_m1.values(), '', label=r"Método 1 ($c_p$ promedio)")
# plt.plot(phi_list, cp_p_m2_def.values(), '', label="Método 2 (balance de entalpías)")

plt.grid(True, which='both', alpha=0.5)

plt.xlabel(r"$\phi$")
plt.ylabel(r"$\tilde{{c_p}}_{p}$ (J/(kgK))")
plt.legend()

plt.xlim(phi_list[0], phi_list[-1])

plt.savefig(f"./Res/T_ad/cp_p_vs_phi_{comburente}.svg")
plt.show(block=False)
