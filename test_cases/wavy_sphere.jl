# Test case of a wavy sphere with a cubic density distribution
include("../福島.jl")

import PhysicalConstants.CODATA2018: G
import Unitful: @u_str, uconvert

function calculate_and_time(layer, evp, δ)
    units = u"m^2/s^2"
    @time result = 福島2017_compute_potential(layer, evp, δ)
    println("Numerical result at ", evp, ": ", uconvert(units, result))
end

R_0 = 6380.0u"km"

function ρ(n, ϕ, λ)::typeof(1.0u"g/cm^3")
    if n == 0
        ϕ * λ * 1u"g/cm^3"
    elseif n == 1
        4.7975u"g/cm^3"
    elseif n == 2
        6u"g/cm^3"
    elseif n == 3
        2u"g/cm^3"
    end
end

layer = FiniteBodyLayer(
    G,
    3,
    ρ,
    (ϕ, λ) -> R_0 + 20u"km" * sin(ϕ * 5) + 50u"km",
    (ϕ, λ) -> R_0 - 20u"km" * sin(λ * 5) - 50u"km",
    R_0
)

evp = EvaluationPoint(R_0, 0, 0)

δ = 10^-6

println("R_0, zero point")
calculate_and_time(layer, EvaluationPoint(R_0, 0, 0), δ)

println("R_0, another point (should NOT be the same as the last one)")
calculate_and_time(layer, EvaluationPoint(R_0, 0.25π, 0.25π), δ)

println("R_0, zero point, double the precision")
calculate_and_time(layer, EvaluationPoint(R_0, 0, 0), δ^2)

println("center point")
calculate_and_time(layer, EvaluationPoint(0u"m", 0, 0), δ)

println("3000 km")
calculate_and_time(layer, EvaluationPoint(3000.0u"km", 0, 0), δ)

println("8000 km")
calculate_and_time(layer, EvaluationPoint(8000.0u"km", 0, 0), δ)
