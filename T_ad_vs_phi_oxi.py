#%% https://www.cantera.org/3.1/userguide/flame-temperature.html

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science'])

# Definir todos los objetos species según el modelo Gri30
species_dict = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}
complete_species = [species_dict[S] for S in ("CH4", "O2", "CO2", "H2O")]
gas_com = ct.Solution(thermo="ideal-gas",
                    species=complete_species,
                    transport_model='mixture-averaged',
                    kinetics='gas')

gas = ct.Solution('gri30.yaml')

#%% Cálculo de T_ad
#Crear vector con valores de ratio de equivalencia.
phi_list = np.linspace(0.2, 1.8, 100)

# Inicializar (crear con todo 0) elos vectores.
T_ad_complete = np.zeros(phi_list.shape)
T_ad_incomplete = np.zeros(phi_list.shape)
    # ERROR: poner T_ad_incomplete = T_ad_complete.

T_0 = 298 # K
p = ct.one_atm

for phi in range(len(phi_list)):
    # Se restablece T y p iniciales en cada bucle para calcular la T_ad
    # para cada phi cuando se parte de estas cond. iniciales.
    gas_com.TP = T_0, p
    gas.TP = T_0, p

    # El método set_equivalence_ratio de la clase Solution toma
        # · 1 valor phi,
        # · 1 str con nombres de especies y su X en el fuel (si se sabe)
        # · 1 str con nombres de especies y sus X en el oxidizer (si se saben)
    gas_com.set_equivalence_ratio(phi_list[phi], "CH4", "O2")
        # CH4 no lleva X porque no hay más especies en fuel (X=1).

    gas.set_equivalence_ratio(phi_list[phi], "CH4", "O2")

    # La función equilibrate() calcula el estado de equilibrio, a p y T ctes.,
    # que minimiza el potencial de Gibbs. Como esta combustión es espontánea,
    # ese estado final es el posterior a la combustión y como hemos impuesto H cte., su T es la T_ad.
    gas_com.equilibrate("HP")
    gas.equilibrate("HP")

    # Por tanto la T de gas_mix ahora será la Tad
    T_ad_complete[phi] = gas_com.T
    T_ad_incomplete[phi] = gas.T

#%% Plot T - phi
plt.figure(figsize=(8,8))

plt.plot(phi_list,
        T_ad_complete,
        label="Complete combustion",
        marker="."
        )
plt.plot(phi_list,
        T_ad_incomplete,
        label="Incomplete (GRI3.0)",
        marker="."
        )
    # phi es vector de abscisas y T_ad_incomplete, de ordenadas
    # lw: line width

ax = plt.gca()
# Ajuste Eje Y (Números cada 200, líneas cada 50)
ax.yaxis.set_major_locator(plt.MultipleLocator(200))
ax.yaxis.set_minor_locator(plt.MultipleLocator(50))
# Ajuste Eje X (Números cada 0.2, líneas cada 0.05)
ax.xaxis.set_major_locator(plt.MultipleLocator(0.2))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))

plt.grid(True, which='both', alpha=0.5)

plt.xlabel("Equivalence ratio, "+ r"$\phi$"+ f"\n \n p = {p} Pa    $T_{{0}}$ = {T_0} K")
plt.ylabel("Temperature [K]")

plt.legend(loc='best', fontsize=10)
    # muestra las label definidas en plt.plot

# plt.xlim(0.6,1.8)
# plt.ylim(1400,2400)

plt.savefig("plots/T_ad_vs_phi_oxi.svg")
plt.show()
    # muestra el gráfico
    # plt.show() debe ir después de plt.savefig()
