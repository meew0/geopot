# Compare the numerical to the analytical solution
include("../福島.jl")
include("../layer.jl")

import PhysicalConstants.CODATA2018: G
import Unitful: @u_str, uconvert

ρ = 1u"g/cm^3"
R_0 = 6380.0u"km"

layer = FiniteBodyLayer(
    G,
    0,
    (n, ϕ, λ) -> ρ,
    (ϕ, λ) -> (λ <= π) ? R_0 : 0u"km",
    (ϕ, λ) -> 0u"km",
    R_0
)

println("Mass")
@time println("Numerical: ", uconvert(u"kg", layer_mass(layer)))
analytical_mass = (2π / 3) * ρ * R_0^3
println("Analytical: ", uconvert(u"kg", analytical_mass))

println("Volume")
@time println("Numerical: ", uconvert(u"m^3", layer_volume(layer)))
volume = (2π / 3) * R_0^3
println("Analytical: ", uconvert(u"m^3", volume))

println("Centre of mass")
@time println("Numerical radius: ", uconvert(u"km", layer_com(layer).R))
println("Analytical radius: ", 3 * R_0 / 8)
