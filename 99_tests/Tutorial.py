#%% https://cantera.org/stable/userguide/python-tutorial.html
import cantera as ct
import numpy as np

# Crear un objeto de fase gas usando el mecanismo GRI 3.0 (metano)
gas = ct.Solution('gri30.yaml')

#%% Para ver estado de mezcla, se llama al objeto como si fuera una función
gas()

#%% Definir estado: Temperatura, Presión y Fracciones molares
# Función gas.TPX lee temperatura, presión y fracciones molares
gas.TPX = 1200, ct.one_atm, 'CH4:1, O2:2, N2:7.52'
gas()

#%% Ver resultado inicial
print(f"Estado inicial: {gas.T} K, {gas.P} Pa")

#%% Realizar un paso de reacción química (equilibrio)
gas.equilibrate('HP')
gas()
print(f"Temperatura tras combustión: {gas.T:.2f} K")
