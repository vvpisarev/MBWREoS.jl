"""
    molar_mass(c::AbstractEoSComponent)

Return the molar mass of a component.
"""
function molar_mass end

"""
    carbon_number(c::AbstractEoSComponent)

Return the number of carbons atoms in the hydrocarbon chain.
"""
function carbon_number end

"""
    name(c::AbstractEoSComponent)

Return the component species name.
"""
function name end

"""
    describe(c::AbstractEoSComponent)

Return `Dict` of parameters. Useful for logging.
"""
function describe end

"""
    load(::Type{T}; name::AbstractString, databases) where {T<:AbstractEoSComponent}

Load component `::T` by its `name` by joining the data on this species in `databases`.
"""
function load end


Base.length(m::AbstractEoSMixture) = m.number_of_components
Base.show(io::IO, x::AbstractEoSMixture) = print(io, "$(typeof(x))($(name(x)))")

components(m::AbstractEoSMixture) = m.components
name(m::AbstractEoSMixture) = join(map(name, components(m)), " + ")
describe(m::AbstractEoSMixture) = Dict{String,Any}("noparameters" => NaN)
load(::Type{T}; names, physdb) where {T<:AbstractEoSMixture} = error("NotImpemented")

function describe(x::MBWREoSComponent)
    return Dict{String,Any}(
        "data structure" => repr(x),
        "name" => name(x),
        "critical pressure [Pa]" => x.Pc,
        "critical temperature [K]" => x.RTc / GAS_CONSTANT_SI,
        "pitzer acentric factor" => x.acentric_factor,
        "molar mass [kg mol⁻¹]" => x.molar_mass,
        "number of carbons atoms" => x.carbon_number,
        "eos" => "MBWR (propane as reference)",
    )
end
