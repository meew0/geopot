# Compare numerical and analytical solutions of mass, volume, and
# the radius of the centre of mass, for a hemisphere
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

analytical_volume = (2π / 3) * R_0^3

@time println("Numerical mass: ", uconvert(u"kg", layer_mass(layer)))
println("Analytical mass: ", uconvert(u"kg", analytical_volume * ρ))

@time println("Numerical volume: ", uconvert(u"m^3", layer_volume(layer)))
println("Analytical volume: ", uconvert(u"m^3", analytical_volume))

@time println("Numerical radius of CoM: ", uconvert(u"km", layer_centre_of_mass(layer).R))
println("Analytical radius of CoM: ", 3 * R_0 / 8)
