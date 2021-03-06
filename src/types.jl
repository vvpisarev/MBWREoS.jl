abstract type AbstractEoSComponent end

abstract type AbstractEoSMixture end

const NothingOrT{T} = Union{Nothing,T}

struct MBWREoSComponent{T<:Number} <: AbstractEoSComponent
    # meta information
    name::String

    # physical parameters
    Pc::T  # critical pressure
    acentric_factor::T
    RTc::T   # R * critical temperature
    Zc::T # critical Z-factor
    molar_mass::T  # [kg mol⁻¹] molar mass
    carbon_number::Int  # [dimless] number of carbons

    function MBWREoSComponent{T}(
        ;
        name::AbstractString="No Name",
        critical_pressure::Number=NaN,
        critical_temperature::Number=NaN,
        acentric_factor::Number=NaN,
        Zc::Number=NaN,
        molar_mass::Number=NaN,
        carbon_number::Integer=0,
        kw...
    ) where {T}
        RTc = GAS_CONSTANT_SI * critical_temperature
        new{T}(
            name,
            critical_pressure,
            acentric_factor,
            RTc,
            Zc,
            molar_mass,
            carbon_number,
        )
    end
end
MBWREoSComponent(; x...) = MBWREoSComponent{Float64}(; x...)

Base.eltype(::MBWREoSComponent{T}) where {T} = T

for func in (:molar_mass, :name, :acentric_factor, :carbon_number)
    expr = :($(func)(c::MBWREoSComponent) = getfield(c, $(QuoteNode(func))))
    eval(expr)
    eval(:(export $func))
end

#=
Mixture
=#

struct MBWREoSMixture{T} <: AbstractEoSMixture
    components::Vector{MBWREoSComponent{T}}

    kij::Matrix{T} # binary interaction coefficient for f
    lij::Matrix{T} # binary interaction coefficient for h

    function MBWREoSMixture(
        ;
        components::AbstractVector{MBWREoSComponent{T}},
        kij::NothingOrT{AbstractMatrix}=nothing,
        lij::NothingOrT{AbstractMatrix}=nothing,
        kw...
    ) where {T}
        nc = length(components)
        kmatr = kij === nothing ? zeros(T, nc, nc) : kij
        lmatr = lij === nothing ? zeros(T, nc, nc) : lij
        new{T}(components, kmatr, lmatr)
    end
end

@inline Base.@propagate_inbounds function Base.getindex(
    mix::MBWREoSMixture,
    i::Integer
)
    return mix.components[i]
end

struct MBWRThermoBuffer{T}
    fij::Matrix{T}
    hij::Matrix{T}
    matr::Matrix{T}
    vec1::Vector{T}
    vec2::Vector{T}
end

function MBWRThermoBuffer{T}(n::Integer) where {T}
    fij = Matrix{T}(undef, n, n)
    hij = similar(fij)
    matr = similar(fij)
    vec1 = similar(fij, (n,))
    vec2 = similar(vec1)
    return MBWRThermoBuffer{T}(fij, hij, matr, vec1, vec2)
end

MBWRThermoBuffer(n::Integer) = MBWRThermoBuffer{Float64}(n)

function MBWRThermoBuffer(mix::MBWREoSMixture{T}) where {T}
    nc = ncomponents(mix)
    return MBWRThermoBuffer{T}(nc)
end

function MBWRThermoBuffer(mix::MBWREoSMixture{Tm}, nmol::AbstractVector{Tn}) where {Tm, Tn}
    nc = ncomponents(mix)
    nc == length(nmol) ||
        throw(DimensionMismatch("Number of mixture components is not equal to the length of moles vector"))
    T = promote_type(Tm, Tn)
    return MBWRThermoBuffer{T}(nc)
end

"""
    thermo_buffer(mix[, nmol])

Create a buffer for intermediate calculations of mixture thermodynamic properties.

See also: [`pressure`](@ref), [`log_c_activity`](@ref), [`log_c_activity!`](@ref),
[`log_c_activity_wj`](@ref), [`log_c_activity_wj!`](@ref)
"""
thermo_buffer(mix::MBWREoSMixture) = MBWRThermoBuffer(mix)
thermo_buffer(mix::MBWREoSMixture, nmol) = MBWRThermoBuffer(mix, nmol)
