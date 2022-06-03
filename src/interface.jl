acentric_factor(x::MBWREoSComponent) = x.acentric_factor
CubicEoS.carbon_number(x::MBWREoSComponent) = x.carbon_number
CubicEoS.molar_mass(x::MBWREoSComponent) = x.molar_mass
CubicEoS.name(x::MBWREoSComponent) = x.name

CubicEoS.ncomponents(mix::MBWREoSMixture) = length(mix.components)
CubicEoS.components(x::MBWREoSMixture) = x.components
CubicEoS.thermo_buffer(mix::MBWREoSMixture) = MBWRThermoBuffer(mix)
CubicEoS.thermo_buffer(mix::MBWREoSMixture, nmol) = MBWRThermoBuffer(mix, nmol)
