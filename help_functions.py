import numpy as np
import scipy.constants as const
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns

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
def plot_fits(x_data, E_shifted, params_quadratic, params_polynomial, params_morse):
    X_axis = np.linspace(x_data.min()*1.1, x_data.max()*1.1, 200)
    plt.figure(figsize=(6, 4))
    plt.scatter(x_data, E_shifted, label='Data', marker='^', color=colors[3], s=40, zorder=5)
    plt.plot(X_axis, quadratic_model(X_axis, *params_quadratic), label='Quadratic Fit', color=colors[0])
    plt.plot(X_axis, polynomial_model(X_axis, *params_polynomial), label='Polynomial Fit', color=colors[1])
    plt.plot(X_axis, morse_model(X_axis, *params_morse), label='Morse Fit', color=colors[2])
    plt.xlabel('Position (Angstrom)')
    plt.ylabel('Energy (eV)')
    plt.title('Potential Energy Surface Fitting')
    plt.legend()
    plt.grid()
    sns.despine()
    plt.savefig(f'{folder}potential_fits.pdf', bbox_inches='tight')
    plt.show()

# ============================= Analytical solution Morse =============================
# plot energy levels function
def plot_energy_levels(params_morse, energy_levels, x_vals, e):
    plt.plot(x_vals, morse_model(x_vals, *params_morse), label='Morse Fit', color=colors[2])
    for n, energy in enumerate(energy_levels):
        min_index = np.where(morse_model(x_vals, *params_morse) <= energy / e)[0][0]
        max_index = np.where(morse_model(x_vals, *params_morse) <= energy / e)[0][-1]
        xi, xj = x_vals[min_index], x_vals[max_index]
        plt.hlines(energy / e, xi, xj, colors='r', label='Energy levels' if n == 0 else "")
        # place level label centered above the line
        xmid = 0.5 * (xi + xj)
        vp = morse_model(x_vals, *params_morse)
        y_offset = 0.01 * (vp.max() - vp.min())
        plt.text(xmid, energy / e + y_offset, f"n={n}", ha='center', va='bottom', color='r', fontsize=8)
        
    plt.xlabel('Position (Angstrom)')
    plt.ylabel('Energy (eV)')
    plt.title('Morse Potential with Energy Levels')
    plt.legend()
    sns.despine()
    plt.grid()
    plt.savefig(f'{folder}morse_energy_levels.pdf', bbox_inches='tight')
    plt.show()