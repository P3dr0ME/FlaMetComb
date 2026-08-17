#%% PREÁMBULO
# https://www.cantera.org/3.1/userguide/flame-temperature.html

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
from scipy.optimize import newton

plt.style.use(['science'])

#%% INPUTS
type = "air" # "oxi" or "air"
T_r = 298 # K
p = ct.one_atm
N=100 # Número de puntos en el vector phi_list

#Crear vector con valores de ratio de equivalencia.
phi_list = np.linspace(0.2, 1.8, N)

#Definir el oxidante
oxidizer = "O2" if type == "oxi" else "O2:1, N2:3.76"

# Especies
species_dict = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}
ideal_species = [species_dict[S] for S in ("CH4", "O2", "N2", "CO2", "H2O")]

# Gases
gas_ideal = ct.Solution(thermo="ideal-gas",
                    species=ideal_species,
                    transport_model='mixture-averaged',
                    kinetics='gas')
gas_real = ct.Solution('gri30.yaml')


#%% CÁLCULO DE T_Ad Método 1 (Analítico con cp promedio)

M_CH4 = 16.04e-3;    M_aire = 137.33e-3;    M_O2 = 32e-3 # kg/mol
f_s = 1*M_CH4/(2*M_aire) if type == "air" else 1*M_CH4/(2*M_O2) # Dosado estequiométrico
LHV = 50.048e6 # J/kg. Springer, Appendix 1.

# Inicializar diccionarios
T_ad_ideal_aux, T_ad_ana_m1, cp_ave, q_p_ana_m1 = (
    {phi: None for phi in phi_list} for _ in range(4)
)

for phi in phi_list:
    gas_ideal.TP = T_r, p
    gas_ideal.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
    gas_ideal.equilibrate("HP")
    T_ad_ideal_aux[phi] = gas_ideal.T
    gas_ideal.TP = (T_ad_ideal_aux[phi]+T_r)/2, p # T=(T_ad+T_r)/2

    cp_ave[phi] = gas_ideal.cp

    if 0<phi<=1: # pobre
        T_ad_ana_m1[phi] = T_r + ( phi*f_s*LHV ) / ( (1+phi*f_s)*cp_ave[phi] )
    elif phi>1: # rica
        T_ad_ana_m1[phi] = T_r + ( f_s*LHV ) / ( (1+phi*f_s)*cp_ave[phi] )

    q_p_ana_m1[phi] = (T_ad_ana_m1[phi]-T_r)*cp_ave[phi] / 1e6 # MJ/kg
        # Idéntica a la q_p_analítica de q_p.py (misma fórmula al fin y al cabo)



#%% CÁLCULO DE T_Ad Método 2 (Analítico H_p=H_r)
Dh0_f = {"CH4": -74.87e3, "O2": 0, "N2": 0, "CO2": -393.52e3, "H2O": -241.83e3} #J/mol
# Springer, Table 2.2.

# Excesos de combustible y de comburente en función de phi
eps = { phi:  0 if phi<=1 else phi-1  for phi in phi_list }
delta = { phi:  2/phi - 2 if phi<=1 else 0  for phi in phi_list }

# Moles de especies en función de phi
n_ir = { phi:  {"CH4": 1+eps[phi], "O2": 2+delta[phi], "N2": 79/21*(2+delta[phi]) , "CO2": 0, "H2O": 0} for phi in phi_list }
n_ip = { phi:  {"CH4": eps[phi], "O2": delta[phi], "N2": 79/21*(2+delta[phi]) , "CO2": 1, "H2O": 2} for phi in phi_list }

# Calor específico molar
cp_molar_r = {sp: gas_ideal.species(f"{sp}").thermo.cp((T_r-298)/2) / 1e3 for sp in Dh0_f} # J/mol/K

# Inicializar diccionarios
T_ad_ana_m2, q_p_ana_m2 = (
    {phi: None for phi in phi_list} for _ in range(2)
)

for phi in phi_list:

    # Calor de combustión
    Q0_p = (
        sum(n_ip[phi][sp] * Dh0_f[sp] for sp in Dh0_f)
      - sum(n_ir[phi][sp] * Dh0_f[sp] for sp in Dh0_f)
    )

    # Función cuya raíz se busca
    def f(T_p):
        cp_molar_p = { sp: gas_ideal.species(sp).thermo.cp((T_p + 298)/2) / 1e3 for sp in Dh0_f }  # J/mol/K
        H_p = sum( n_ip[phi][sp] * cp_molar_p[sp] * (T_p - 298) for sp in Dh0_f )
        H_r = sum( n_ir[phi][sp] * cp_molar_r[sp] * (T_r - 298) for sp in Dh0_f )

        return H_p - (-Q0_p + H_r)

    # Resolver f(Tad) = 0
    T_ad_ana_m2[phi] = newton(f, x0=2000, tol=1e-6, maxiter=100)
    q_p_ana_m2[phi] = -Q0_p / ( (n_ir[phi]["CH4"]*M_CH4)*(1+1/(phi*f_s)) ) / 1e6 # MJ/kg

#%% CÁLCULO DE T_ad Método 3 (Cantera)
# Inicializar diccionarios
T_ad_ideal, T_ad_real, cp_ideal, cp_real, q_p_ideal, q_p_real = (
    {phi: None for phi in phi_list} for _ in range(6)
)

for phi in phi_list:

    gas_ideal.TP = T_r, p
    gas_real.TP = T_r, p
    # Se restablece T y p iniciales en cada bucle para calcular la T_ad
    # para cada phi cuando se parte de estas cond. iniciales.

    gas_ideal.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
    gas_real.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
        # El método set_equivalence_ratio de la clase Solution toma
            # · 1 valor phi,
            # · 1 str con nombres de especies y su X en el fuel (si se sabe)
            # · 1 str con nombres de especies y sus X en el oxidizer (si se saben)
        # CH4 no lleva X porque no hay más especies en fuel (X=1).

    gas_ideal.equilibrate("HP")
    gas_real.equilibrate("HP")
        # La función equilibrate() calcula el estado de equilibrio, a p y T ctes.,
        # que minimiza el potencial de Gibbs. Como esta combustión es espontánea,
        # ese estado final es el posterior a la combustión y como hemos impuesto H cte., su T es la T_ad.

    #cp másicos
    cp_real[phi] = gas_real.cp
    cp_ideal[phi] = gas_ideal.cp

    # Por tanto la T de gas_mix ahora será la Tad
    T_ad_ideal[phi] = gas_ideal.T
    T_ad_real[phi] = gas_real.T

    # Calor de combustión por kg de mezcla:
    q_p_ideal[phi] = (T_ad_ideal[phi]-T_r)*cp_ideal[phi] / 1e6 # MJ/kg
    q_p_real[phi] = (T_ad_real[phi]-T_r)*cp_real[phi] / 1e6 # MJ/kg
        # Aquí se calcula q_p de los gases de Cantera con la fórmula analítica de T_ad.
        # Por eso los q_p no salen idénticos a los de q_p.py.


#%% GRÁFICO T - phi
plt.figure(figsize=(8,8))
plt.title(
    r"\bf{Temperatura\ adiabática\ de\ llama\ frente\ a\ ratio\ de\ equivalencia}" + "\n"
    f"{'Oxígeno' if type == 'oxi' else 'Aire'} \n"
    fr"$p = {p}$ Pa $\quad T_{{0}} = {T_r}$ K",
    fontsize=11,
    pad=15
)

plt.plot(phi_list, T_ad_real.values(), label="Real (GRI-Mech 3.0)", marker="")
plt.plot(phi_list, T_ad_ideal.values(), label="Ideal (GRI-Mech 3.0)", marker="")
plt.plot(phi_list, T_ad_ana_m1.values(), label=r"Analítica ($c_p$ promedio)", marker="")
plt.plot(phi_list, T_ad_ana_m2.values(), label=fr"Analítica (balance de entalpías)", marker="")

plt.grid(True, which='both', alpha=0.5)

plt.xlabel("Ratio de equivalencia, "+ r"$\phi$")
plt.ylabel(r"Temperatura adiabática, $T_{ad}$ [K]")
plt.legend()

plt.xlim(phi_list[0], phi_list[-1])

plt.savefig(f"./plots/T_ad/T_ad_vs_phi_{type}.svg")
plt.show() # muestra el gráfico y debe ir después de plt.savefig()

#%% GRÁFICO DE q_p vs phi
plt.figure(figsize=(8,8))

plt.plot(phi_list, [q_p_real[phi] for phi in phi_list], '', label="Real (Gri-Mech 3.0)")
plt.plot(phi_list, [q_p_ideal[phi] for phi in phi_list], '', label="Ideal (Gri-Mech 3.0)")
plt.plot(phi_list, [q_p_ana_m1[phi] for phi in phi_list], '.', label=r"Analítica ($c_p$ promedio)")
plt.plot(phi_list, [q_p_ana_m2[phi] for phi in phi_list], '-', label="Analítica (balance de entalpías)")


plt.grid(True, which='both', alpha=0.5)

plt.xlabel("Ratio de equivalencia, " + r"$\phi$")
plt.ylabel("Calor de combustión por kg de mezcla, " + r"$q_p$ (MJ/kg)")
plt.legend()

plt.show()

# %%
