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
            (layer.ρ(n, ϕ, λ) * r^(n + 1) * layer.R_0^-n / (n + 1)) for
            n = 0:layer.N
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

# Calculate the position of a layer's centre of mass by first calculating its
# moments about the three Cartesian axis planes, then dividing those by the mass
# to get the Cartesian coordinates of the centre of mass, and finally
# transforming those into spherical coordinates.
function layer_com(layer::FiniteBodyLayer)::EvaluationPoint
    function moment_integrand(x, cos_ϕ_offset, cos_λ_offset, cos_λ_exponent)
        ϕ, λ = x

        # Antiderivative of the density along a radius (capital rho)
        Ρ(r) = sum(
            (layer.ρ(n, ϕ, λ) * r^(n + 1) * layer.R_0^-n / (n + 1)) for
            n = 0:layer.N
        )

        ψ = cos(ϕ + cos_ϕ_offset) * cos(λ + cos_λ_offset)^cos_λ_exponent

        top = layer.R_T(ϕ, λ)
        bottom = layer.R_B(ϕ, λ)
        ustrip(
            u"kg*m",
            ψ * cos(ϕ) * (1 / 3) * (Ρ(top) * top^3 - Ρ(bottom) * bottom^3),
        )
    end

    # Lower and upper bounds of integration
    l, u = ((-π / 2, 0), (π / 2, 2π))

    # Calculate moments
    # Moment in xy plane: * -sin ϕ
    M_xy, E = hcubature(x -> moment_integrand(x, π / 2, 0, 0), l, u)

    # Moment in xz plane: * cos ϕ sin λ
    M_xz, E = hcubature(x -> moment_integrand(x, 0, -π / 2, 1), l, u)

    # Moment in yz plane: * cos ϕ cos λ
    M_yz, E = hcubature(x -> moment_integrand(x, 0, 0, 1), l, u)

    # Calculate mass and convert moments to coordinates
    M = layer_mass(layer)
    x, y, z = (M_yz, M_xz, M_xy) .* u"kg*m" ./ M

    # Convert Cartesian to polar coordinates
    R = hypot(x, y, z)
    Φ = asin(z / R)
    Λ = atan(y, x) # (atan2)

    EvaluationPoint(R, Φ, Λ)
end

function gravitational_potential(layer::FiniteBodyLayer)::GravitationalPotential
    error("NYI")
end
