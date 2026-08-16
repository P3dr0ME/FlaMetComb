#%% PREÁMBULO
# https://www.cantera.org/3.1/userguide/flame-temperature.html

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science'])

#%% INPUTS
type = "air" # "oxi" or "air"
T_r = 298 # K
p = ct.one_atm
N=100 # Número de puntos en el vector phi_list

#%% CÁLCULO DE T_ad
species_dict = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}
ideal_species = [species_dict[S] for S in ("CH4", "O2", "N2", "CO2", "H2O")]
gas_ideal = ct.Solution(thermo="ideal-gas",
                    species=ideal_species,
                    transport_model='mixture-averaged',
                    kinetics='gas')
gas_ideal_aux  = ct.Solution(thermo="ideal-gas", # Para calcular cp_ave en T=(Tad-T_r)/2
                    species=ideal_species,
                    transport_model='mixture-averaged',
                    kinetics='gas')

gas_real = ct.Solution('gri30.yaml')

#Crear vector con valores de ratio de equivalencia.
phi_list = np.linspace(0.2, 1.8, N)

oxidizer = "O2" if type == "oxi" else "O2:1, N2:3.76"

M_CH4 = 16.04e-3 # kg/mol
M_aire = 137.33e-3 # kg/mol
M_O2 = 32e-3 # kg/mol
f_s = 1*M_CH4/(2*M_aire) if type == "air" else 1*M_CH4/(2*M_O2) # Dosado estequiométrico
LHV = 50.048e6 # J/kg. Springer, Appendix 1.

# Inicializar diccionarios
T_ad_ideal, T_ad_real, T_ad_analitica, cp_ave, cp_ideal, cp_real, cp_ave, q_p_ideal, q_p_real, q_p_analitica = (
    {phi: None for phi in phi_list} for _ in range(10)
)

for i in range(len(phi_list)):
    phi=phi_list[i]
    # Se restablece T y p iniciales en cada bucle para calcular la T_ad
    # para cada phi cuando se parte de estas cond. iniciales.
    gas_ideal.TP = T_r, p
    gas_real.TP = T_r, p
    # El método set_equivalence_ratio de la clase Solution toma
        # · 1 valor phi,
        # · 1 str con nombres de especies y su X en el fuel (si se sabe)
        # · 1 str con nombres de especies y sus X en el oxidizer (si se saben)
    gas_ideal.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
        # CH4 no lleva X porque no hay más especies en fuel (X=1).
    gas_real.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
    # La función equilibrate() calcula el estado de equilibrio, a p y T ctes.,
    # que minimiza el potencial de Gibbs. Como esta combustión es espontánea,
    # ese estado final es el posterior a la combustión y como hemos impuesto H cte., su T es la T_ad.
    gas_ideal.equilibrate("HP");     gas_real.equilibrate("HP")
    #cp másicos
    cp_real[phi] = gas_real.cp;    cp_ideal[phi] = gas_ideal.cp
    # Por tanto la T de gas_mix ahora será la Tad
    T_ad_ideal[phi] = gas_ideal.T
    T_ad_real[phi] = gas_real.T

    gas_ideal_aux.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
    gas_ideal_aux.equilibrate("HP")
    gas_ideal_aux.TP = (T_ad_ideal[phi]+T_r)/2, p # T=(T_ad+T_r)/2
    cp_ave[phi] = gas_ideal_aux.cp # lista de cp (másicos) de productos de reacción completa
    if 0<phi<=1: # pobre
        T_ad_analitica[phi] = T_r + ( phi*f_s*LHV ) / ( (1+phi*f_s)*cp_ave[phi] )
    elif phi>1: # rica
        T_ad_analitica[phi] = T_r + ( f_s*LHV ) / ( (1+phi*f_s)*cp_ave[phi] )

    # Calor de combustión por kg de mezcla:
    q_p_ideal[phi] = (T_ad_ideal[phi]-T_r)*cp_ideal[phi] / 1e6 # MJ/kg
    q_p_real[phi] = (T_ad_real[phi]-T_r)*cp_real[phi] / 1e6 # MJ/kg
        # Aquí se calcula q_p de los gases de Cantera con la fórmula analítica de T_ad.
        # Por eso los q_p no salen idénticos a los de q_p.py.
    q_p_analitica[phi] = (T_ad_analitica[phi]-T_r)*cp_ave[phi] / 1e6 # MJ/kg
        # Idéntica a la q_p_analítica de q_p.py (misma fórmula al fin y al cabo)


#%% GRÁFICO T - phi
plt.figure(figsize=(8,8))
plt.title(
    r"\bf{Temperatura\ adiabática\ de\ llama\ frente\ a\ ratio\ de\ equivalencia}" + "\n"
    f"{'Oxígeno' if type == 'oxi' else 'Aire'} \n"
    fr"$p = {p}$ Pa $\quad T_{{0}} = {T_r}$ K",
    fontsize=11,
    pad=15
)

plt.plot(phi_list,
        T_ad_real.values(),
        label="Real (GRI-Mech 3.0)",
        marker="."
        )

plt.plot(phi_list,
        T_ad_ideal.values(),
        label="Ideal (GRI-Mech 3.0)",
        marker="."
        )

plt.plot(phi_list,
        T_ad_analitica.values(),
        label="Analítica",
        marker=""
        )

ax = plt.gca()
# Ajuste Eje Y (Números cada 200, líneas cada 50)
ax.yaxis.set_major_locator(plt.MultipleLocator(200))
ax.yaxis.set_minor_locator(plt.MultipleLocator(50))
# Ajuste Eje X (Números cada 0.2, líneas cada 0.05)
ax.xaxis.set_major_locator(plt.MultipleLocator(0.2))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))

plt.grid(True, which='both', alpha=0.5)

plt.xlabel("Ratio de equivalencia, "+ r"$\phi$")
plt.ylabel(fr"Temperatura adiabática, T_{{ad}} [K]")

plt.legend(loc='best', fontsize=10)
    # muestra las label definidas en plt.plot

plt.xlim(phi_list[0], phi_list[-1])
# plt.ylim(1400,2400)

plt.savefig(f"./plots/T_ad/T_ad_vs_phi_{type}.svg")
plt.show() # muestra el gráfico y debe ir después de plt.savefig()

# PLOT q_p vs phi
#%% GRÁFICO DE q_p vs phi
plt.figure(figsize=(8,8))

plt.plot(phi_list, [q_p_real[phi] for phi in phi_list], 'o-', label="Real (Gri-Mech 3.0)")
plt.plot(phi_list, [q_p_ideal[phi] for phi in phi_list], 'o-', label="Ideal (Gri-Mech 3.0)")
plt.plot(phi_list, [q_p_analitica[phi] for phi in phi_list], 's--', label="Analítico (con LHV)")


plt.grid(True, which='both', alpha=0.5)

plt.xlabel("Ratio de equivalencia, " + r"$\phi$")
plt.ylabel("Calor de combustión por kg de mezcla, " + r"$q_p$ (MJ/kg)")
plt.legend()

plt.show()
