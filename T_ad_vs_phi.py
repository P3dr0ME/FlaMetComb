#%% https://www.cantera.org/3.1/userguide/flame-temperature.html

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science','ieee'])

# Definir todos los objetos "Species" según el modelo Gri30
species = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}
# species es un diccionario: {'H2': <Species H2>, 'H': <Species H>, ...}
# keys: nombres de especies
# values: objetos Species


#%%
# Crea una lista a partir de values de species (<Species H2> por ejemplo) pero añadiendo
# solo los values cuyas keys son las especies de la reacción de combustión completa.
complete_species = [species[S] for S in ("CH4", "O2", "N2", "CO2", "H2O")]
# Para que fuera un dict: complete_comb_species =
# {S: species[S] for S in ("CH4", "O2", "N2", "CO2", "H2O")}

# Para reacción incomplete añadimos las 53 especies del GRI3.0
incomplete_species = species.values()

#%%
gas_mix_complete = ct.Solution(thermo="ideal-gas",
                               species=complete_species,
                               transport_model='mixture-averaged',
                               kinetics='gas')
# A diferencia de Tutorial.py, aquí solo damos a la clase Solution "ideal gas"
# (el modelo termo que queremos) y la lista reducida de especies sacadas de gri30.yaml,
# en lugar de sacar todo tal cual de gri30.yaml.
gas_mix_incomplete = ct.Solution(thermo="ideal-gas",
                                 species=incomplete_species,
                                 transport_model='mixture-averaged',
                                 kinetics='gas')
# print(f"Transport model: {gas_mix_incomplete.transport_model}")
    # de gri30.yaml no he sacado transport_model, por eso no está definido



#%%
#Crear vector con valores de ratio de equivalencia.
phi = np.linspace(0.6, 1.8, 100)

# Inicializar (crear con todo 0) elos vectores.
T_ad_complete = np.zeros(phi.shape)
T_ad_incomplete = np.zeros(phi.shape)
    # ERROR: poner T_ad_incomplete = T_ad_complete.

T_0 = 298 # K
p = ct.one_atm


for i in range(len(phi)):
    #print("%"*80)
    #print(f"Equivalence ratio {phi[i]}")

    # Se restablece T y p iniciales en cada bucle para calcular la T_ad
    # para cada phi cuando se parte de estas cond. iniciales.
    gas_mix_complete.TP = T_0, p
    gas_mix_incomplete.TP = T_0, p

    # El método set_equivalence_ratio de la clase Solution toma
        # · 1 valor phi,
        # · 1 str con nombres de especies y su X en el fuel (si se sabe)
        # · 1 str con nombres de especies y sus X en el oxidizer (si se saben)
    gas_mix_complete.set_equivalence_ratio(phi[i], "CH4", "O2:1, N2:3.76")
        # CH4 no lleva X porque no hay más especies en fuel (X=1).

    gas_mix_incomplete.set_equivalence_ratio(phi[i], "CH4", "O2:1, N2:3.76")

    # Calcula el equilibrio químico, manteniendo entalpía y presión ctes., de gas_mix y lo guarda en gas_mix.
    # El equlibrio acaba calculando la combustión porque es la reacción espontánea.
    gas_mix_complete.equilibrate("HP")
    gas_mix_incomplete.equilibrate("HP")

    # Por tanto la T de gas_mix ahora será la Tad
    T_ad_complete[i] = gas_mix_complete.T
    T_ad_incomplete[i] = gas_mix_incomplete.T
    #print(f"T_ad_complete = {T_ad_complete[i]}")
    #print(f"T_ad_incomplete = {T_ad_incomplete[i]}")

#%% Plot T - phi
plt.figure(figsize=(8,8))
# Tamaño de ventana

plt.plot(phi,
        T_ad_complete,
        label="complete combustion",
        lw=2)
plt.plot(phi,
        T_ad_incomplete,
        label="GRI3.0 (Cantera)",
        lw=2)
    # phi es vector de abscisas y T_ad_incomplete, de ordenadas
    # lw: line width

# ax = plt.gca()
# # Ajuste Eje Y (Números cada 200, líneas cada 50)
# ax.yaxis.set_major_locator(plt.MultipleLocator(200))
# ax.yaxis.set_minor_locator(plt.MultipleLocator(50))
# # Ajuste Eje X (Números cada 0.2, líneas cada 0.05)
# ax.xaxis.set_major_locator(plt.MultipleLocator(0.2))
# ax.xaxis.set_minor_locator(plt.MultipleLocator(0.05))

plt.grid(True, which='both', alpha=0.5)

plt.xlabel(f"Equivalence ratio, $\phi$ \n p = {p} Pa    $T_{{0}}$ = {T_0} K")
plt.ylabel("Temperature [K]")

plt.legend(loc='best', fontsize=10)
    # muestra las label definidas en plt.plot

# plt.xlim(0.6,1.8)
# plt.ylim(1400,2400)

plt.show()
    # muestra el gráfico
# plt.savefig("T_ad_vs_phi.svg")
