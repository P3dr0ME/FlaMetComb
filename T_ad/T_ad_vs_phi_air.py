#%% https://www.cantera.org/3.1/userguide/flame-temperature.html

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science'])

# Definir todos los objetos species según el modelo Gri30
species_dict = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}
ideal_species = [species_dict[S] for S in ("CH4", "O2", "N2", "CO2", "H2O")]
gas_ideal = ct.Solution(thermo="ideal-gas",
                    species=ideal_species,
                    transport_model='mixture-averaged',
                    kinetics='gas')
gas_ideal_aux  = ct.Solution(thermo="ideal-gas", # Para calcular cp_ave en T=(Tad-T_0)/2
                    species=ideal_species,
                    transport_model='mixture-averaged',
                    kinetics='gas')

gas = ct.Solution('gri30.yaml')

#%% Cálculo de T_ad
#Crear vector con valores de ratio de equivalencia.
phi_list = np.linspace(0.2, 1.8, 100)

# Inicializar (crear con todo 0) elos vectores.
T_ad_ideal = np.zeros(phi_list.shape)
T_ad_incomplete = np.zeros(phi_list.shape)
    # ERROR: poner T_ad_incomplete = T_ad_ideal.
T_ad_analytic = np.zeros(phi_list.shape)
cp_ave = np.zeros(phi_list.shape)

T_0 = 298 # K
p = ct.one_atm

M_CH4 = 16.04e-3 # kg/mol
M_aire = 137.33e-3 # kg/mol
f_s = 1*M_CH4/(2*M_aire) # Dosado estequiométrico
LHV = 50.048e6 # J/kg. Springer, Appendix 1.

for i in range(len(phi_list)):
    phi=phi_list[i]
    # Se restablece T y p iniciales en cada bucle para calcular la T_ad
    # para cada phi cuando se parte de estas cond. iniciales.
    gas_ideal.TP = T_0, p
    gas.TP = T_0, p
    # El método set_equivalence_ratio de la clase Solution toma
        # · 1 valor phi,
        # · 1 str con nombres de especies y su X en el fuel (si se sabe)
        # · 1 str con nombres de especies y sus X en el oxidizer (si se saben)
    gas_ideal.set_equivalence_ratio(phi, "CH4", "O2:1, N2:3.76")
        # CH4 no lleva X porque no hay más especies en fuel (X=1).
    gas.set_equivalence_ratio(phi, "CH4", "O2:1, N2:3.76")
    # La función equilibrate() calcula el estado de equilibrio, a p y T ctes.,
    # que minimiza el potencial de Gibbs. Como esta combustión es espontánea,
    # ese estado final es el posterior a la combustión y como hemos impuesto H cte., su T es la T_ad.
    gas_ideal.equilibrate("HP")
    gas.equilibrate("HP")
    # Por tanto la T de gas_mix ahora será la Tad
    T_ad_ideal[i] = gas_ideal.T
    T_ad_incomplete[i] = gas.T

    gas_ideal_aux.set_equivalence_ratio(phi, "CH4", "O2:1, N2:3.76")
    gas_ideal_aux.TP = (T_ad_ideal[i]+T_0)/2, p # T=(Tad+T0)/2
    cp_ave[i] = gas_ideal_aux.cp # lista de cp (másicos) de productos de reacción completa
    if 0<phi<=1: # pobre
        T_ad_analytic[i] = T_0 + ( phi*f_s*LHV ) / ( (1+phi*f_s)*cp_ave[i] )
    elif phi>1: # rica
        T_ad_analytic[i] = T_0 + ( f_s*LHV ) / ( (1+phi*f_s)*cp_ave[i] )


#%% Plot T - phi
plt.figure(figsize=(8,8))

plt.plot(phi_list,
        T_ad_ideal,
        label="Ideal (Cantera - GRI3.0)",
        marker="."
        )

plt.plot(phi_list,
        T_ad_incomplete,
        label="Incompleta (Cantera - GRI3.0)",
        marker="."
        )
    # phi es vector de abscisas y T_ad_incomplete, de ordenadas
    # lw: line width

plt.plot(phi_list,
        T_ad_analytic,
        label="Analítica",
        marker="."
        )

ax = plt.gca()
# Ajuste Eje Y (Números cada 200, líneas cada 50)
ax.yaxis.set_major_locator(plt.MultipleLocator(200))
ax.yaxis.set_minor_locator(plt.MultipleLocator(50))
# Ajuste Eje X (Números cada 0.2, líneas cada 0.05)
ax.xaxis.set_major_locator(plt.MultipleLocator(0.2))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))

plt.grid(True, which='both', alpha=0.5)

plt.xlabel("Ratio de equivalencia, "+ r"$\phi$"+ f"\n \n p = {p} Pa    $T_{{0}}$ = {T_0} K")
plt.ylabel("Temperatura adiabática, T_{{ad}} [K]")

plt.legend(loc='best', fontsize=10)
    # muestra las label definidas en plt.plot

plt.xlim(phi_list[0], phi_list[-1])
# plt.ylim(1400,2400)

plt.savefig("./plots/T_ad/T_ad_vs_phi_air.svg")
plt.show()
    # muestra el gráfico
    # plt.show() debe ir después de plt.savefig()
