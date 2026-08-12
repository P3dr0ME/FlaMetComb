#https://cantera.org/3.1/userguide/heating-value.html
import cantera as ct
import numpy as np


# INPUTS:
type = "oxi" # "oxi" or "air"
T_r = 298 # K
p = ct.one_atm
N = 5

#%% Cálculo de q_p
gas = ct.Solution('gri30.yaml')
products = ct.Solution('gri30.yaml') # para sacar X productos.

phi_list = np.linspace(0.2, 1.8, N)

# Inicializar (crear con todo 0) elos vectores.
q_p = {phi: None for phi in phi_list}

oxidizer = "O2" if type == "oxi" else "O2:1, N2:3.76"

for i in range(len(phi_list)):
    phi=phi_list[i]

    products.TP = T_r, p
    products.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
    products.equilibrate("HP")

    #H_R(T_R)
    gas.TP = T_r, p
    gas.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
    Y_CH4 = gas.Y[gas.species_index("CH4")]
    H_1 = gas.enthalpy_mass # J/kg

    #H_P(T_R)
    gas.X = products.X # Se iguala la composición de gas a la de products, que es la de los productos de combustión.
    # gas.X = {"CO2": 1, "H2O": 2} # Con estas X, para phi=1 se obtiene q_p=LHV.
    gas.TP = T_r, p
    H_2 = gas.enthalpy_mass # J/kg

    q_p[phi] = phi * -(H_2 - H_1) / Y_CH4 / 1e6 if phi <= 1 else -(H_2 - H_1) / Y_CH4 / 1e6 # MJ/kg. Dividir por Y_CH4 para obtener q_p por kg de combustible, no por kg de mezcla.
    # q_p es como el LHV, pero para cualquier phi.

    if phi == 1:
        print(gas())


#%% Tabla q_p vs phi
print("q_p values:")
for phi, value in q_p.items():
    print(f"  φ = {phi}: {value:.2f} MJ/kg")
