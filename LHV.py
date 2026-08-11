#%% https://www.cantera.org/3.1/userguide/flame-temperature.html

import cantera as ct
import numpy as np


# INPUTS:
type = "oxi" # "oxi" or "air"
T_r = 298 # K
p = ct.one_atm
N = 5

#%% Cálculo de LHV
gas = ct.Solution('gri30.yaml')
gas_aux = ct.Solution('gri30.yaml') # para sacar X productos.


phi_list = np.linspace(0.2, 1.8, N)

# Inicializar (crear con todo 0) elos vectores.
LHV = {phi: None for phi in phi_list}

oxidizer = "O2" if type == "oxi" else "O2:1, N2:3.76"

for i in range(len(phi_list)):
    phi=phi_list[i]

    print(gas.X)
    gas_aux.TP = T_r, p
    gas_aux.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
    gas_aux.equilibrate("HP")
    print(gas.X)

    #H_R(T_R)
    gas.TP = T_r, p
    gas.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
    H1 = gas.enthalpy_mass # J/kg
    Y_CH4 = gas.Y[gas.species_index("CH4")]
    print(gas.X)

    gas.X = gas_aux.X.copy() # Se iguala la composición de gas a la de gas_aux, que es la de los productos de combustión.
    H2 = gas.enthalpy_mass # J/kg


    LHV[phi] = -(H2 - H1) / Y_CH4 / 1e6 # MJ/kg. LHV va por masa de CH4, mientras que gas.enthalpy_mass va por masa de mezcla. Por eso se divide entre Y_CH4.

#%% Tabla LHV vs phi
print("LHV values:")
for phi, value in LHV.items():
    print(f"  φ = {phi}: {value:.2f} MJ/kg")
