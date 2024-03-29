# Compare the numerical to the analytical solution for a uniform spherical shell
include("../福島.jl")
include("../layer.jl")

import PhysicalConstants.CODATA2018: G
import Unitful: @u_str, uconvert

# Calculate the gravitational potential using both analytical and numerical
# methods and print the results.
function compare_analytical_numerical(layer, evp, δ)
    units = u"m^2/s^2"
    println("Calculating numerical result")
    @time result = 福島2017_compute_potential(layer, evp, δ)
    println("Numerical result at ", evp, ": ", uconvert(units, result))
    result = 福島2017_spherical_layer_potential_analytic(layer, evp)
    println("Analytical result at ", evp, ": ", uconvert(units, result))
end

ρ = 1u"g/cm^3"
R_0 = 6380.0u"km"
R_T, R_B = R_0 + 10.0u"km", R_0 - 40.0u"km"

layer = FiniteBodyLayer(
    G,
    0,
    (n, ϕ, λ) -> ρ,
    (ϕ, λ) -> R_T,
    (ϕ, λ) -> R_B,
    R_0
)

evp = EvaluationPoint(R_0, 0, 0)

δ = 10^-6

println("R_0, zero point")
compare_analytical_numerical(layer, EvaluationPoint(R_0, 0, 0), δ)

println("R_0, another point (should be the same as the last one)")
compare_analytical_numerical(layer, EvaluationPoint(R_0, 0.25π, 0.25π), δ)

println("R_0, zero point, double the precision")
compare_analytical_numerical(layer, EvaluationPoint(R_0, 0, 0), δ^2)

println("center point")
compare_analytical_numerical(layer, EvaluationPoint(0u"m", 0, 0), δ)

println("3000 km")
compare_analytical_numerical(layer, EvaluationPoint(3000.0u"km", 0, 0), δ)

println("8000 km")
compare_analytical_numerical(layer, EvaluationPoint(8000.0u"km", 0, 0), δ)

println("Mass")
@time println("Numerical: ", uconvert(u"kg", layer_mass(layer)))
H = R_T - R_B
analytical_mass = (4π / 3) * ρ * (R_T^2 + R_T * R_B + R_B^2) * H
println("Analytical: ", uconvert(u"kg", analytical_mass))

println("Volume")
@time println("Numerical: ", uconvert(u"m^3", layer_volume(layer)))
outer_volume = (4 / 3) * π * R_T^3
inner_volume = (4 / 3) * π * R_B^3
println("Analytical: ", uconvert(u"m^3", outer_volume - inner_volume))

println("Centre of mass")
@time println("Numerical radius: ", uconvert(u"km", layer_centre_of_mass(layer).R))
println("Analytical radius: ", 0)
