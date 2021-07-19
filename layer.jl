# Utility functions for working with finite body layers

include("types.jl")
include("福島.jl")

# Calculate a layer's total mass as a volume integral of the density at each
# point. The first integral of the density (along the radius) is calculated
# analytically.
function layer_mass(layer::FiniteBodyLayer)::Mass
    function area_density(x)
        ϕ, λ = x

        # Antiderivative of the density along a radius (capital rho)
        Ρ(r) = sum(
            (layer.ρ(n, ϕ, λ) * r^(n + 1) * R_0^-n / (n + 1)) for n = 0:layer.N
        )

        top = layer.R_T(ϕ, λ)
        bottom = layer.R_B(ϕ, λ)
        ustrip(
            u"kg",
            cos(ϕ) * (1 / 3) * (Ρ(top) * top^2 - Ρ(bottom) * bottom^2),
        )
    end

    I, E = hcubature(area_density, (-π / 2, 0), (π / 2, 2π))
    I * u"kg"
end

# Calculate a layer's total volume as a surface integral of the thickness at
# each point
function layer_volume(layer::FiniteBodyLayer)::Volume
    function thickness(x)
        ϕ, λ = x
        ustrip(
            u"m^3",
            cos(ϕ) * (1 / 3) * (layer.R_T(ϕ, λ)^3 - layer.R_B(ϕ, λ)^3),
        )
    end

    I, E = hcubature(thickness, (-π / 2, 0), (π / 2, 2π))
    I * u"m^3"
end

function gravitational_potential(layer::FiniteBodyLayer)::GravitationalPotential
    error("NYI")
end
