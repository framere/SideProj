using LinearAlgebra
using QuantumOptics
using LsqFit
using Plots

# read the data from potential_data
data = readdlm("potential_data", skipstart=1)

x = data[:, 1]
E = data[:, 2]

# define the model function and fit
model(x, p) = p[1] .+ p[2] * x .+ p[3] * x.^2 .+ p[4] * x.^3
initial_params = [1000.0, - 1.0, 0.0, 0.0]
fit_result = curve_fit(model, x, E, initial_params)

best_fit_params = fit_result.param

println("Best fit parameters: ", best_fit_params)

# plot the data and the fit
x_plot = range(minimum(x), stop=maximum(x), length=100)
scatter(x, E, label="Data", xlabel="x", ylabel="E [eV]", title="CO Vibrational Potential on Pt(111) Surface")
plot!(x_plot, model(x_plot, best_fit_params), label="Fit", linewidth=2)
