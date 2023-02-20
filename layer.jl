# Utility functions for working with finite body layers

include("types.jl")

using HCubature
import Unitful: @u_str, ustrip

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

# The moment integrand for centre of mass calculation (see next function).
# Implementing it as a macro like this, rather than as a function taking
# coefficients for the trigonometric functions used to compute ψ, may not
# be the most stylish solution but it significantly reduces compilation
# time.
macro moment_integrand(ψ_expr)
    return quote
        function(x)
            ϕ, λ = x

            # Antiderivative of the density along a radius (capital rho)
            Ρ(r) = sum(
                (layer.ρ(n, ϕ, λ) * r^(n + 1) * layer.R_0^-n / (n + 1)) for
                n = 0:layer.N
            )

            ψ = $ψ_expr

            top = layer.R_T(ϕ, λ)
            bottom = layer.R_B(ϕ, λ)
            ustrip(
                u"kg*m",
                ψ * cos(ϕ) * (1 / 4) * (Ρ(top) * top^3 - Ρ(bottom) * bottom^3),
            )
        end
    end
end

# Calculate the position of a layer's centre of mass by first calculating its
# moments about the three Cartesian axis planes, then dividing those by the mass
# to get the Cartesian coordinates of the centre of mass, and finally
# transforming those into spherical coordinates.
function layer_centre_of_mass(layer::FiniteBodyLayer, maxevals = 1000)::EvaluationPoint
    # Lower and upper bounds of integration
    l, u = ((-π / 2, 0), (π / 2, 2π))

    # Calculate moments
    # Moment in xy plane: * -sin ϕ
    xy_moment_integrand = @moment_integrand -sin(ϕ)
    M_xy, E = hcubature(xy_moment_integrand, l, u, maxevals = maxevals)

    # Moment in xz plane: * cos ϕ sin λ
    xz_moment_integrand = @moment_integrand cos(ϕ) * sin(λ)
    M_xz, E = hcubature(xz_moment_integrand, l, u, maxevals = maxevals)

    # Moment in yz plane: * cos ϕ cos λ
    yz_moment_integrand = @moment_integrand cos(ϕ) * cos(λ)
    M_yz, E = hcubature(yz_moment_integrand, l, u, maxevals = maxevals)

    # Calculate mass and convert moments to coordinates
    M = layer_mass(layer)
    x, y, z = (M_yz, M_xz, M_xy) .* u"kg*m" ./ M

    # Convert Cartesian to spherical coordinates
    R = hypot(x, y, z)
    Φ = asin(z / R)
    Λ = atan(y, x) # (atan2)

    EvaluationPoint(R, Φ, Λ)
end
