import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.linalg import eigh

colors = sns.color_palette("tab10")
folder = 'figures/'

# ================================= Define fits ================================
# --- Quadratic model ---
def quadratic_model(x, a, b, c):
    return a + b * (x - c)**2

# --- Polynomial model (cubic) ---
def polynomial_model(x, a3, a2, a1, a0):
    return a3 * x**3 + a2 * x**2 + a1 * x + a0

# --- Morse potential model ---
def morse_model(x, D_e, a, r_e, offset):
    return D_e * (1 - np.exp(-a * (x - r_e)))**2 + offset

# ============================= Fit models to data =============================
def fit_model(x, model, guess, E_shifted):
    popt, pcov = curve_fit(model, x, E_shifted, p0=guess)
    return popt, pcov


# plot results function
def plot_fit(x_data, E_shifted, model, params, name):
    X_axis = np.linspace(x_data.min()*1.1, x_data.max()*1.1, 200)
    plt.figure(figsize=(6, 4))
    plt.scatter(x_data, E_shifted, label='Data', marker='^', color=colors[3], s=40, zorder=5)
    plt.plot(X_axis, model(X_axis, *params), label=f'{name} Fit', color=colors[0])
    plt.xlabel('Position (Angstrom)')
    plt.ylabel('Energy (eV)')
    plt.title('Potential Energy Surface Fitting')
    plt.legend()
    plt.grid()
    sns.despine()
    plt.savefig(f'{folder}potential_fit_{name}.pdf', bbox_inches='tight')
    plt.close()

# ============================ Numerical solution ============================
def solve_schrodinger(m_eff, x_min, x_max, N, potential_func, potential_params, num_levels, e, hbar):
    # Build grid
    x_A = np.linspace(x_min, x_max, N)
    x = x_A * 1e-10  # Convert Angstrom to meters
    dx = x[1] - x[0]

    # Construct potential energy array
    V_ev = potential_func(x_A, *potential_params)
    V = V_ev * e  # Convert eV to Joules

    # enforce high walls to help localization
    V[0]  = V[-1] = np.max(V)*2

    # Construct Hamiltonian matrix using finite difference method
    # --- kinetic (finite-diff) ---
    coeff = hbar**2 / (2.0 * m_eff * dx**2)
    diag  = 2*coeff + V
    off   = -coeff * np.ones(N-1)

    # Dirichlet boundary conditions 
    H = np.diag(diag) + np.diag(off,1) + np.diag(off,-1)


    # Solve eigenvalue problem
    E_J, PSI = eigh(H, eigvals=(0,num_levels-1))
    E_eV = E_J / e

    # normalize wavefunctions
    for n in range(num_levels):
        norm = np.sqrt(np.sum(np.abs(PSI[:,n])**2)*dx)
        PSI[:,n] /= norm

    return x_A, E_eV, PSI

# plot numerical energy levels function
def plot_numerical_energy_levels(x_A, potential_func, potential_params, E_eV, PSI, model):
    plt.figure(figsize=(6,4))
    V_eV = potential_func(x_A, *potential_params)
    plt.plot(x_A, V_eV, 'k', lw=1.5, label='V(x)')

    for n, energy in enumerate(E_eV):
        min_index = np.where(potential_func(x_A, *potential_params) <= energy)[0][0]
        max_index = np.where(potential_func(x_A, *potential_params) <= energy)[0][-1]
        xi, xj = x_A[min_index], x_A[max_index]
        plt.hlines(energy, xi, xj, colors='r')
        # place level label centered above the line
        xmid = 0.5 * (xi + xj)
        vp = potential_func(x_A, *potential_params)
        y_offset = 0.01 * (vp.max() - vp.min())
        plt.text(xmid, energy + y_offset, f"n={n}", ha='center', va='bottom', color='r', fontsize=8)
        plt.plot(x_A, energy + 0.15*PSI[:,n]/np.max(np.abs(PSI[:,n])),
                label=f'n={n}, E={energy:.4f} eV')

    plt.xlabel('x (Å)')
    plt.ylabel('Energy (eV)')
    plt.legend()
    plt.title('Bound states from finite difference')
    plt.grid(True)
    sns.despine()
    plt.savefig(f'{folder}numerical_solution_{model}.pdf', bbox_inches='tight')
    plt.close()

# ============================= Analytical solution Quadratic =============================
def analytical_quadratic_energy_levels(params_quadratic, m_eff, hbar_ineV, e, N_energy_levels):
    a, b, c = params_quadratic
    b_J_m2 = b * e / (1e-10)**2
    omega = np.sqrt(2 * b_J_m2 / m_eff)
    energy_levels = []
    for n in range(N_energy_levels):
        energy = hbar_ineV * omega * (n + 0.5) + a
        energy_levels.append(energy * e)  # convert to Joules
    return energy_levels


# ============================= Analytical solution Morse =============================
def analytical_morse_energy_levels(params_morse, m_eff, hbar, e, N_energy_levels):
    # convert in SI units for further calculations
    D_e = params_morse[0] * e  # in J
    a = params_morse[1] / 1e-10  # in 1/m
    r_e = params_morse[2] * 1e-10  # in m
    offset = params_morse[3] * e  # in J

    # Solving the Schrödinger equation
    lambda_param = np.sqrt(2 * m_eff * D_e) / (a * hbar) # see parametrization wikipedia
    n_max = int(2 * lambda_param - 0.5)
    energy_levels = []
    for n in range(N_energy_levels):
        epsilon_n = (2 * lambda_param - n - 0.5) * (n + 0.5) * (hbar * a)**2 / (2 * m_eff)
        energy_levels.append(epsilon_n + offset)

    return energy_levels

# plot energy levels function
def plot_energy_levels(model, params, energy_levels, x_min, x_max, e, name):
    plt.figure(figsize=(6, 4))
    x_vals = np.linspace(x_min, x_max, 1000)
    plt.plot(x_vals, model(x_vals, *params), label='Fit', color=colors[2])
    for n, energy in enumerate(energy_levels):
        min_index = np.where(model(x_vals, *params) <= energy / e)[0][0]
        max_index = np.where(model(x_vals, *params) <= energy / e)[0][-1]
        xi, xj = x_vals[min_index], x_vals[max_index]
        plt.hlines(energy / e, xi, xj, colors='r', label='Energy levels' if n == 0 else "")
        # place level label centered above the line
        xmid = 0.5 * (xi + xj)
        vp = model(x_vals, *params)
        y_offset = 0.01 * (vp.max() - vp.min())
        plt.text(xmid, energy / e + y_offset, f"n={n}", ha='center', va='bottom', color='r', fontsize=8)
        
    plt.xlabel('Position (Angstrom)')
    plt.ylabel('Energy (eV)')
    plt.title(f'Bound states from analytical {name} solution')
    plt.legend()
    sns.despine()
    plt.grid()
    plt.savefig(f'{folder}analytical_solution_{name}.pdf', bbox_inches='tight')
    plt.close()