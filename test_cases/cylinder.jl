# Compare numerical and analytical solutions of mass, volume, and
# the radius of the centre of mass, for a cylinder
include("../福島.jl")
include("../layer.jl")

import PhysicalConstants.CODATA2018: G
import Unitful: @u_str, uconvert

ρ = 1u"g/cm^3"
R_0 = 6380.0u"km"

# Dimension of the cylinder (both radius and half-height)
d_cyl = 5000.0u"km"

function R_T(ϕ, λ)
    if ϕ < (-π/4)
        -d_cyl / sin(ϕ)
    elseif ϕ > (π/4)
        d_cyl / sin(ϕ)
    else
        d_cyl / cos(ϕ)
    end
end

layer = FiniteBodyLayer(
    G,
    0,
    (n, ϕ, λ) -> ρ,
    R_T,
    (ϕ, λ) -> 0u"km",
    R_0
)

analytical_volume = 2π * d_cyl^3

@time println("Numerical mass: ", uconvert(u"kg", layer_mass(layer)))
println("Analytical mass: ", uconvert(u"kg", ρ * analytical_volume))

@time println("Numerical volume: ", uconvert(u"m^3", layer_volume(layer)))
println("Analytical volume: ", uconvert(u"m^3", analytical_volume))

@time println("Numerical radius of CoM: ", uconvert(u"km", layer_centre_of_mass(layer).R))
println("Analytical radius of CoM: ", 0u"km")
