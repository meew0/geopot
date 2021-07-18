# Calculation of the gravitational potential of a single layer within a general
# finite-body-like solar system object, according to Fukushima, Toshio (2017):
# Precise and fast computation of gravitational field of general finite body and
# its application to gravitational study of asteroid Eros. Astron J 154:145
# (https://doi.org/10.3847/1538-3881/aa88b8)

using DoubleExponentialFormulas
import Unitful: @u_str, ustrip

# Define some shorthands for dimensions we need
Length = typeof(1.0u"m")
Angle = Float64 # for now
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
    N::Int64

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

# Closed-form expressions for the value of the first coefficient function C at
# degrees n <= 3 (equations 16, 18, 20, 22)
function low_degree_C(n::Int64, B::Float64, α::Float64)::Float64
    if n == 0
        α^2 - 3α * B + (3 / 2) * B^2
    elseif n == 1
        (α - B) * (α^2 - 5α * B + (5 / 2) * B^2)
    elseif n == 2
        α^4 - 10 * α^3 * B + (45 / 2) * α^2 * B^2 - (35 / 2) * α * B^3 +
        (35 / 8) * B^4
    elseif n == 3
        (α - B) * (
            α^4 - 14 * α^3 * B + (77 / 2) * α^2 * B^2 - (63 / 2) * α * B^3 +
            (63 / 8) * B^4
        )
    else
        error(
            "Coefficients higher than 3 are not implemented as closed form expressions",
        )
    end
end

# Closed-form expressions for the value of the second coefficient function D at
# degrees n <= 3 (equations 17, 19, 21, 23)
function low_degree_D(n::Int64, B::Float64, α::Float64, ζ::Float64)::Float64
    if n == 0
        (2α - (3 / 2) * B) + ζ * (1 / 2)
    elseif n == 1
        (3 * α^2 - (35 / 6) * α * B + (5 / 2) * B^2) +
        ζ * ((3 / 2) * α - (5 / 6) * B) +
        ζ^2 * (1 / 3)
    elseif n == 2
        (4 * α^3 - (43 / 3) * α^2 * B + (175 / 12) * α * B^2 - (35 / 8) * B^3) +
        ζ * (3 * α^2 - (49 / 12) * α * B + (35 / 24) * B^2) +
        ζ^2 * ((4 / 3) * α - (7 / 12) * B) +
        ζ^3 * (1 / 4)
    elseif n == 3
        (
            5 * α^4 - (85 / 3) * α^3 * B + (1001 / 20) * α^2 * B^2 -
            (273 / 8) * α * B^3 + (63 / 8) * B^4
        ) +
        ζ * (
            5 * α^3 - (145 / 12) * α^2 * B + (399 / 40) * α * B^2 -
            (21 / 8) * B^3
        ) +
        ζ^2 * ((10 / 3) * α^2 - (69 / 20) * α * B + (21 / 20) * B^2) +
        ζ^3 * ((5 / 4) * α - (9 / 20) * B) # In the paper, (9/20) * B^2 is used here instead of (9/20) * B, which I assume is an error.
        ζ^4 * (1 / 5)
    else
        error(
            "Coefficients higher than 3 are not implemented as closed form expressions",
        )
    end
end

# Equation (62): recurrence relation for the computation of higher order
# integrals, implemented iteratively to avoid problems with the Julia compiler.
function iterative_L(j::Int64, t::Float64, L::Float64, S::Float64, A::Float64)::Float64
    # Find the root value of L at which the recurrence relation would terminate
    j_root = j % 2
    L_j = S * j_root + L * (1 - j_root) # L if j_root == 0, S if 1. Avoids branching

    for j in j_root:2:j
        jL_j = t^(j - 1) * S - (j - 1) * A * L_j
        L_j = jL_j / j
    end

    L_j
end

# Equation (14): the indefinite radial line integral, expressed analytically
function I_n(n::Int64, layer::FiniteBodyLayer, evp::EvaluationPoint, r::Length, ϕ::Angle, λ::Angle)::Float64
    # Calculate non-dimensional quantities from equation (11)
    α = evp.R / layer.R_0
    ξ = ϕ - evp.Φ
    η = λ - evp.Λ
    ζ = (r - evp.R) / layer.R_0

    # Calculate values of auxiliary functions from equations (9) and (10)
    cosine_term = cos(ϕ) * cos(evp.Φ) * sin(η / 2)^2
    B = 2α * (sin(ξ / 2)^2 + cosine_term)
    A = B * (2α - B)

    # Calculate the normalised distance between the internal point and the
    # evaluation point, according to equation (8)
    S = sqrt(A + (B + ζ)^2)

    # Calculate the value of the logarithmic function L, according to equation
    # (15)
    if B + ζ >= 0
        L = log(S + B + ζ)
    else
        L = log(A / (S - B - ζ))
    end

    if false#n <= 3
        # Calculate the coefficient polynomials using low-degree closed form
        # expressions
        C = low_degree_C(n, B, α)
        D = low_degree_D(n, B, α, ζ)

        # Bring everything together
        C * L + D * S
    else
        # Use the compact expression of the radial integral in terms of the
        # recurrence relation for L instead (equation (58))
        sum(
            binomial(n + 2, j) *
            (α - B)^(n + 2 - j) *
            iterative_L(j, ζ + B, L, S, A) for j = 0:(n+2)
        )
    end
end

# Equation (13): Kernel function of degree n, calculated as a difference between
# two line integrals
function K_n(n::Int64, layer::FiniteBodyLayer, evp::EvaluationPoint, ϕ::Angle, λ::Angle)::Float64
    top_radius = layer.R_T(ϕ, λ)
    top_term = I_n(n, layer, evp, top_radius, ϕ, λ)

    bottom_radius = layer.R_B(ϕ, λ)
    bottom_term = I_n(n, layer, evp, bottom_radius, ϕ, λ)

    top_term - bottom_term
end

# Equation (30): the radial term that will be doubly integrated.
function Q_n(
    n::Int64,
    layer::FiniteBodyLayer,
    evp::EvaluationPoint,
    ξ::Angle,
    η::Angle,
)::GravitationalPotential
    ϕ = evp.Φ + ξ
    λ = evp.Λ + η

    constant_term = layer.G * layer.R_0^2
    density_term = layer.ρ(n, ϕ, λ)
    kernel_term = K_n(n, layer, evp, ϕ, λ)

    constant_term * density_term * kernel_term
end

# Equation (29): Single integration of the radial term along one latitude line.
function P_n(n::Int64, layer::FiniteBodyLayer, evp::EvaluationPoint, δ::Float64, ξ::Angle)::GravitationalPotential
    f(η) = ustrip(u"m^2/s^2", Q_n(n, layer, evp, ξ, η))
    I, E = quadde(f, 0, 2π, rtol = δ)
    I * u"m^2/s^2"
end

# Equation (28): Split quadrature integration of terms corresponding to a single
# degree of the density polynomial over the entire surface of the sphere.
function V_n(n::Int64, layer::FiniteBodyLayer, evp::EvaluationPoint, δ::Float64)::GravitationalPotential
    # The cos(ϕ) term is NOT present in the final rewritten expression for V in
    # equation (28). It is however present in the non-rewritten expressions in
    # equations (7) and (12), and without it the result is completely wrong,
    # so I assume this term was just forgotten in the final expression.
    f(ξ) = ustrip(u"m^2/s^2", P_n(n, layer, evp, δ, ξ) * cos(evp.Φ + ξ))

    # Split quadrature integration: the rewritten term contains a singularity
    # at ϕ = Φ (⇔ ξ = 0). The DE integration method can tolerate this but only
    # at the boundary of the integration interval. So the integration must be
    # performed separately for the latitudes above and below Φ.
    I, E = quadde(f, -π / 2 - evp.Φ, 0, π / 2 - evp.Φ, rtol = δ)

    # Add the unit that was lost in integration
    I * u"m^2/s^2"
end

# Equation (6): calculation of the gravitational potential of the given layer
# at a given evaluation point, computed as a finite series of terms
# corresponding to individual degree terms of the polynomial approximation of
# the layer's mass density.
function V(layer::FiniteBodyLayer, evp::EvaluationPoint, δ::Float64)::GravitationalPotential
    sum(V_n(n, layer, evp, δ) for n = 0:layer.N)
end

# Alias for the function V to be used by external code
福島2017_compute_potential = V

# Equation (44): analytically compute the potential of a homogeneous spherical
# layer. Assumes that the three functions specified in the layer are constant!
function 福島2017_spherical_layer_potential_analytic(
    layer::FiniteBodyLayer,
    evp::EvaluationPoint,
)
    R_T = layer.R_T(evp.Φ, evp.Λ)
    R_B = layer.R_B(evp.Φ, evp.Λ)
    ρ = layer.ρ(0, evp.Φ, evp.Λ)
    H = R_T - R_B

    if evp.R < R_B
        2π * layer.G * ρ * (R_T + R_B) * H
    elseif evp.R < R_T
        (2π * layer.G * ρ / 3) * (3R_T^2 - evp.R^2 - 2R_B^3 / evp.R)
    else
        M = (4π / 3) * ρ * (R_T^2 + R_T * R_B + R_B^2) * H
        layer.G * M / evp.R
    end
end
