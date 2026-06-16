import cantera as ct

# 1. Create a pure fluid object for Methane
# Methane.yaml is included in Cantera's data files
fluid = ct.PureFluid('gri30.yaml')

# 2. Define the pressure (1 atm for normal boiling point)
pressure = 101325 # Pa

# 3. Set the state to saturated liquid/vapor
# We set Pressure and a vapor fraction (0 for liquid, 1 for vapor)
fluid.P = pressure

# 4. Print the temperature (in Kelvin)
print(f'Methane boiling temperature at {pressure/1000:.2f} kPa: {fluid.T:.2f} K')
print(f'Methane boiling temperature: {fluid.T - 273.15:.2f} °C')
