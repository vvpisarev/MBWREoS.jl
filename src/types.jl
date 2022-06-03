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

struct MBWREoSMixture{T} <: AbstractEoSMixture{T}
    components::Vector{MBWREoSComponent{T}}

    kij::Matrix{T} # binary interaction coefficient for f
    lij::Matrix{T} # binary interaction coefficient for h

    function MBWREoSMixture(;
        components::AbstractVector{MBWREoSComponent{T}},
        kij::Union{Nothing,AbstractMatrix}=nothing,
        lij::Union{Nothing,AbstractMatrix}=nothing,
        kw...
    ) where {T}
        nc = length(components)
        kmatr = isnothing(kij) ? zeros(T, nc, nc) : kij
        lmatr = isnothing(lij) ? zeros(T, nc, nc) : lij
        new{T}(components, kmatr, lmatr)
    end
end

struct MBWRThermoBuffer{T} <: AbstractEoSThermoBuffer
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
