const _DBROOT = joinpath(abspath(@__DIR__), "../data")
const MBWR_PHYSICS_PARAMS = ComponentDatabase(joinpath(_DBROOT, "nist.csv"))

function load(::Type{<:MBWREoSComponent};
    name::AbstractString,
    component_dbs = (MBWR_PHYSICS_PARAMS,)
)
    comp_properties = foldl(
        (dict, db) -> merge!(dict, getentry(db, name)),
        component_dbs;
        init = Dict{Symbol, Any}()
    )

    return __load_mbwr_comp__(; comp_properties...)
end

function __load_mbwr_comp__(;
    name::AbstractString,
    molecular_mass,
    number_carbons,
    critical_temperature,
    critical_pressure,
    critical_compressibility,
    acentric_factor,
    extra_kw...  # unneeded keywords to construct component
)
    __load_mbwr_comp__(
        name,
        molecular_mass,
        number_carbons,
        critical_temperature,
        critical_pressure,
        critical_compressibility,
        acentric_factor,
    )
end

function __load_mbwr_comp__(
    name::AbstractString,
    molecular_mass,
    number_carbons,
    critical_temperature,
    critical_pressure,
    critical_compressibility,
    acentric,
)
    return MBWREoSComponent(;
        name = name,
        critical_pressure = critical_pressure,
        critical_temperature = critical_temperature,
        acentric_factor = acentric,
        Zc = critical_compressibility,
        molar_mass = molecular_mass,
        carbon_number = number_carbons
    )
end

function CubicEoS.load(::Type{<:MBWREoSMixture};
    names,
    component_dbs=(MBWR_PHYSICS_PARAMS, Data.brusilovsky_comp()),
    mix_eos_db::MixtureDatabase=Data.brusilovsky_mix_adjusted(),
)
    components = map(names) do name
        load(MBWREoSComponent; name=name, component_dbs=component_dbs)
    end
    components = collect(components)  # No constructor of mixture for tuple

    corrections = getmatrix(mix_eos_db, names)
    return MBWREoSMixture(; components=components, corrections...)
end
