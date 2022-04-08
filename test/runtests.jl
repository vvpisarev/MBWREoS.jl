using MBWREoS
using CubicEoS
using Test


c1c5 = MBWREoS.load(MBWREoSMixture; names=("methane", "n-pentane"))

volume = 1e-6
nmol = 5000 * [0.4, 0.6] * volume
RT = CubicEoS.GAS_CONSTANT_SI * 300

MBWREoS.log_c_activity(c1c5, nmol, volume, RT)
