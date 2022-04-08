module MBWREoS

export MBWREoSComponent, MBWREoSMixture

import CubicEoS
using CubicEoS: AbstractEoSComponent, AbstractEoSMixture
using CubicEoS: name

using CubicEoSDatabase


include("constants.jl")
include("types.jl")
include("interface.jl")
include("dbload.jl")
include("basic_thermo.jl")
include("chempotential.jl")

end # module
