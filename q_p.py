##% PREÁMBULO
# https://cantera.org/3.1/userguide/heating-value.html
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science'])

#%% INPUTS
type = "oxi" # "oxi" or "air"
T_r = 298 # K
p = ct.one_atm
N = 10

#%% CÁLCULO DE q_p
gas_real = ct.Solution('gri30.yaml')

species_dict = {S.name: S for S in ct.Species.list_from_file("gri30.yaml")}
ideal_species = [species_dict[S] for S in ("CH4", "O2", "N2", "CO2", "H2O")]
gas_ideal = ct.Solution(thermo="ideal-gas",
                    species=ideal_species,
                    transport_model='mixture-averaged',
                    kinetics='gas')

oxidizer = "O2" if type == "oxi" else "O2:1, N2:3.76"

phi_list = np.linspace(0.2, 1.8, N)

# Inicialización de diccionarios de forma compacta
h_r_T_r_real, h_r_T_r_ideal, h_p_T_r_real, h_p_T_r_ideal, q_p_real, q_p_ideal = ( {phi: None for phi in phi_list} for _ in range(6) )
X_r_real, X_p_real, X_r_ideal, X_p_ideal, = ( {phi: {} for phi in phi_list} for _ in range(4) )

# Bucle
for i in range(len(phi_list)):
    phi=phi_list[i]

    # H_R(T_r):
    gas_real.TP = T_r, p
    gas_real.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
    X_r_real[phi] = gas_real.mole_fraction_dict()
    h_r_T_r_real[phi] = gas_real.enthalpy_mass # J/kg

    gas_ideal.TP = T_r, p
    gas_ideal.set_equivalence_ratio(phi, "CH4", f"{oxidizer}")
    X_r_ideal[phi] = gas_ideal.mole_fraction_dict()
    h_r_T_r_ideal[phi] = gas_ideal.enthalpy_mass # J/kg

    # H_P(T_r):
    gas_real.equilibrate("HP")
    X_p_real[phi] = gas_real.mole_fraction_dict()
    gas_real.TP = T_r, p
    # Ahora el gas tiene las X de productos, pero a T_r y p. Por eso se puede calcular h_p(T_r).
    h_p_T_r_real[phi] = gas_real.enthalpy_mass # J/kg
    # Calor de combustión por unidad de masa de MEZCLA (combustible + comburente):
    q_p_real[phi] = -(h_p_T_r_real[phi] - h_r_T_r_real[phi]) / 1e6 # MJ/kg.
        # T_ad = T_r + q_p (phi) / cp_ave (phi)

    gas_ideal.equilibrate("HP")
    X_p_ideal[phi] = gas_ideal.mole_fraction_dict()
    gas_ideal.TP = T_r, p
    h_p_T_r_ideal[phi] = gas_ideal.enthalpy_mass # J/kg
    q_p_ideal[phi] = -(h_p_T_r_ideal[phi] - h_r_T_r_ideal[phi]) / 1e6 # MJ/kg.


#%% CÁLCULO DE q_p_analítico
M_CH4 = 16.04e-3; M_aire = 137.33e-3; M_O2 = 32e-3 # kg/mol
f_s = 1*M_CH4/(2*M_aire) if type == "air" else 1*M_CH4/(2*M_O2) # Dosado estequiométrico
LHV = 50.048e6 # J/kg. Springer, Appendix 1

q_p_analitica = {phi: phi*f_s/(phi*f_s+1)*LHV/1e6 if phi<=1 else f_s/(phi*f_s+1)*LHV/1e6 for phi in phi_list} # MJ/kg. Springer, Appendix 1


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
