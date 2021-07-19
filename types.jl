# Definitions of types used in multiple places

import Unitful: @u_str

# Define some shorthands for dimensions we need
Length = typeof(1.0u"m")
Area = typeof(1.0u"m^2")
Volume = typeof(1.0u"m^3")

Mass = typeof(1.0u"kg")
Angle = Real # for now
GravitationalPotential = typeof(1.0u"m^2/s^2")

# Layer within a general finite body, with its top and bottom surfaces modelled
# as arbitrary real-valued functions of spherical surface coordinates, and its
# density distribution modelled as a low-degree polynomial of the radius at each
# set of surface coordinates.
# Also used to specify the value of the gravitational constant that should be
# used for computation.
struct FiniteBodyLayer
    # Gravitational constant
    G::typeof(1.0u"kg^-1 * m^3 * s^-2")

    # Degree of the density polynomial
    N::Integer

    # Function ρ(n, ϕ, λ) returning the n-th degree coefficient of the density
    # polynomial at the given spherical coordinate
    ρ::Function

    # Functions R_T(ϕ, λ) and R_B(ϕ, λ) returning the radius of the top and
    # bottom surface, respectively, of the layer at the given spherical
    # coordinate.
    R_T::Function
    R_B::Function

    # Radius of the layer's reference sphere
    R_0::Length
end

# A set of three-dimensional spherical coordinates representing one point at
# which the surface integral is evaluated
struct EvaluationPoint
    R::Length
    Φ::Angle
    Λ::Angle
end
