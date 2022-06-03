using MBWREoS
using CubicEoS
using Test


c1c5 = CubicEoS.load(MBWREoSMixture; names=("methane", "n-pentane"))

volume = 1e-6
nmol = 5000 * [0.4, 0.6] * volume
RT = CubicEoS.GAS_CONSTANT_SI * 300

logca = CubicEoS.log_c_activity(c1c5, nmol, volume, RT)
logca = CubicEoS.log_c_activity!(logca, c1c5, nmol, volume, RT)

logca, jacobian = CubicEoS.log_c_activity_wj(c1c5, nmol, volume, RT)
logca, jacobian = CubicEoS.log_c_activity_wj!(logca, jacobian, c1c5, nmol, volume, RT)

vt_stability(c1c5, nmol, volume, RT, CubicEoS.VTStabilityIdealIdentityState)
