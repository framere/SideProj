using LinearAlgebra
using QuantumOptics
using LsqFit
using Plots

# read the data from potential_data
data = read("potential_data", skipheader=1)

x = data[:, 1]
E = data[:, 2]

# define the model
model(x, p) = p[1] .+ p[2] * x.^2 .+ p[3] * x.^3
initial_params = [0.0, 1.0, 0.0]
fit_result = curve_fit(model, x, E, initial_params)

best_fit_params = fit_result.param

println("Best fit parameters: ", best_fit_params)

# plot the data and the fit
scatter(x, E, label="Data", xlabel="x", ylabel="E", title="Data and Fit")
plot!(x, model(x, best_fit_params), label="Fit", linewidth=2)
