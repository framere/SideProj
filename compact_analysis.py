import numpy as np
from scipy.constants import e, hbar, angstrom, physical_constants
import scipy.constants as const
import help_functions as hf

# ============================ Load data =============================
data = np.loadtxt('potential_data', skiprows=1)
x_data = data[:, 0]  # x values in Angstroms
E_data = data[:, 1]  # Energy values in eV
E_shifted = E_data - E_data.min()  # Shift energies so that minimum is at 0 eV

# ============================= Define constants and parameters =============================
mass_u = physical_constants['atomic mass constant'][0]  # Atomic mass unit in kg
omega = 265.495409 *1e-3 * e / hbar  # Convert meV to J and then to angular frequency
m_eff = 13.585135873543415 * mass_u  # Effective mass in kg

# ============================= Select model to fit =============================
MODEL = 'Morse'  # choose between 'quadratic', 'polynomial', and 'morse'

if MODEL == 'Quadratic':
    potential_func = hf.quadratic_model
    guess_quadratic = [0.0, 100.0, 0.0]  # [offset, curvature, center]
    params_quadratic, _ = hf.fit_model(x_data, hf.quadratic_model, guess_quadratic, E_shifted)
    potential_params = params_quadratic
elif MODEL == 'Polynomial':
    potential_func = hf.polynomial_model
    guess_polynomial = [0.0, 0.0, -1.0, 1000.0]  # [a3, a2, a1, a0]
    params_polynomial, _ = hf.fit_model(x_data, hf.polynomial_model, guess_polynomial, E_shifted)
    potential_params = params_polynomial
elif MODEL == 'Morse':
    potential_func = hf.morse_model
    guess_morse = [E_shifted.max() - E_shifted.min(), 1.0, x_data.mean(), E_shifted.min()]  # Initial guess for [D_e, a, r_e, offset]
    params_morse, _ = hf.fit_model(x_data, hf.morse_model, guess_morse, E_shifted)
    potential_params = params_morse

num_levels = 2  # number of energy levels to compute    


# ============================= Plot results =============================
hf.plot_fit(x_data, E_shifted, potential_func, potential_params, MODEL)

# ============================ Numerical calculations =============================
# --- build grid ----
x_min = -0.1
x_max = 0.15
N = 2000     
# --- solve schrodinger ---
x_A, E_numerical, PSI_numerical = hf.solve_schrodinger(m_eff, x_min, x_max, N, potential_func, potential_params, num_levels, e, hbar)

print("Numerical energy levels (in eV):")
for n, energy in enumerate(E_numerical):
    print(f"n={n}: {energy:.4f} eV")
# --- plot numerical energy levels ---
hf.plot_numerical_energy_levels(x_A, potential_func, potential_params, E_numerical, PSI_numerical, MODEL)


# ============================= Analytical solution Morse =============================
if MODEL != 'Morse':
    raise ValueError("Analytical Morse solution can only be computed if the Morse model is selected.")

# --- compute analytical energy levels ---
energy_levels = hf.analytical_morse_energy_levels(params_morse, m_eff, hbar, e, num_levels)

print("Analytical energy levels (in eV):")
for i, energy in enumerate(energy_levels):
    print(f"n={i}: {energy / e:.4f} eV")

x_vals = np.linspace(x_min, x_max, 1000)
hf.plot_energy_levels(params_morse, energy_levels, x_vals, e)