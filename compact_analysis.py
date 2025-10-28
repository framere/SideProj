import numpy as np
import scipy.constants as const
import help_functions as hf

# ============================ Load data =============================
data = np.loadtxt('potential_data', skiprows=1)
x_data = data[:, 0]  # x values in Angstroms
E_data = data[:, 1]  # Energy values in eV
E_shifted = E_data - E_data.min()  # Shift energies so that minimum is at 0 eV

# ============================= Define constants and parameters =============================
e = const.elementary_charge  # Elementary charge in Coulombs
mass_u = const.physical_constants['atomic mass constant'][0]  # Atomic mass unit in kg
hbar = const.hbar  # Reduced Planck's constant in J·s
angstrom = 1e-10  # 1 Angstrom in meters
omega = 265.495409 *1e-3 * e / hbar  # Convert meV to J and then to angular frequency
m_eff = 13.585135873543415 * mass_u  # Effective mass in kg


# ============================= Fit models to data =============================
guess_quadratic = [0.0, 100.0, 0.0]  # [offset, curvature, center]
params_quadratic, _ = hf.fit_model(x_data, hf.quadratic_model, guess_quadratic, E_shifted)

guess_polynomial = [0.0, 0.0, -1.0, 1000.0]  # [a3, a2, a1, a0]
params_polynomial, _ = hf.fit_model(x_data, hf.polynomial_model, guess_polynomial, E_shifted)

guess_morse = [E_shifted.max() - E_shifted.min(), 1.0, x_data.mean(), E_shifted.min()]  # Initial guess for [D_e, a, r_e, offset]
params_morse, _ = hf.fit_model(x_data, hf.morse_model, guess_morse, E_shifted)

# ============================= Plot results =============================
hf.plot_fits(x_data, E_shifted, params_quadratic, params_polynomial, params_morse)


# ============================= Analytical solution Morse =============================
# convert in SI units for further calculations
D_e = params_morse[0] * e  # in J
a = params_morse[1] / angstrom  # in 1/m
r_e = params_morse[2] * angstrom  # in m
offset = params_morse[3] * e  # in J

# Solving the Schrödinger equation
lambda_param = np.sqrt(2 * m_eff * D_e) / (a * hbar) # see parametrization wikipedia
n_max = int(2 * lambda_param - 0.5)
energy_levels = []
N_energy_levels = 2
for n in range(N_energy_levels):
    epsilon_n = (2 * lambda_param - n - 0.5) * (n + 0.5) * (hbar * a)**2 / (2 * m_eff)
    energy_levels.append(epsilon_n + offset)

print("Calculated energy levels (in eV):")
for i, energy in enumerate(energy_levels):
    print(f"n={i}: {energy / e:.4f} eV")

x_vals = np.linspace(-0.14, 0.23, 1000)
hf.plot_energy_levels(params_morse, energy_levels, x_vals, e)

