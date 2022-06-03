function exc_helmholtz(
    mix::MBWREoSMixture, nmol, vol, RT;
    buf = thermo_buffer(mix, nmol),
)
    Tx = RT / GAS_CONSTANT_SI
    total_mol = sum(nmol)
    ρ = total_mol / vol
    fx, hx = scaling_coeffs(mix, nmol, Tx; buf)
    Tref = Tx / fx
    ρref = ρ * hx

    ρmoll = ρref * 1e-3
    G1, G2, G3, G4, G5 = (
        -0.2804337729e-3,
        0.1180666107,
        -0.3756325860e+1,
        0.5624374521e3,
        -0.9354759605e5,
    )
    G6, G7, G8, G9 = (
        -0.4557405505e-4,
        0.1530044332,
        -0.1078107476e+3,
        0.2218072099e+5,
    )
    G10, G11, G12 = (
        0.6629473971e-5,
        -0.6199354447e-2,
        0.6754207966e+1,
    )
    G13 = 0.6472837570e-3
    G14, G15 = -0.6804325262e-1, -0.9726162355e+1
    G16 = 0.5097956459e-2
    G17, G18 = -0.1004655900e-3, 0.4363693352e-1
    G19 = -0.1249351947e-2
    G20, G21 = 0.2644755879e+5, -0.7944237270e+7
    G22, G23 = -0.7299920845e+3, 0.5381095003e+8
    G24, G25 = 0.3450217377e+1, 0.9936666689e+3
    G26, G27 = -0.2166699036, -0.1612103424e+5
    G28, G29 = -0.3633126990e-3, 0.1108612343e+1
    G30, G31, G32 = -0.1330932838e-4, -0.3157701101e-2, 0.1423083811
    iTsq = inv(Tref^2)
    iT = inv(Tref)
    a2 = G1 * Tref + G2 * sqrt(Tref) + G3 + (G4 * Tref + G5) * iTsq
    a3 = G6 * Tref + G7 + (G8 * Tref + G9) * iTsq
    a4 = G10 * Tref + G11 + G12 * iT
    a5 = G13
    a6 = (G14 * Tref + G15) * iTsq
    a7 = G16 * iT
    a8 = (G17 * Tref + G18) * iTsq
    a9 = G19 * iTsq

    b3 = (G20 * Tref + G21) * iT * iTsq
    b5 = (G22 * Tref^2 + G23) * iTsq^2
    b7 = (G24 * Tref + G25) * iT * iTsq
    b9 = (G26 * Tref^2 + G27) * iTsq^2
    b11 = (G28 * Tref + G29) * iT * iTsq
    b13 = (G32 + Tref * (G31 + Tref * G30)) * iTsq^2

    a = (a2, a3 / 2, a4 / 3, a5 / 4, a6 / 5, a7 / 6, a8 / 7, a9 / 8)
    Aexc = evalpoly(ρmoll, a) * ρmoll * 1e3
    ρr = C3_VC * ρref
    ρcmoll = 1e-3 / C3_VC
    In = 1 - exp(-ρr^2)
    bfactor = ρcmoll^2
    Ifactor = ρr^2
    Ar = 0.5 * In * b3 * bfactor
    for (n, b) in enumerate((b5, b7, b9, b11, b13))
        bfactor *= ρcmoll^2
        In = n * In - Ifactor * exp(-ρr^2)
        Ar += 0.5 * In * b * bfactor
        Ifactor *= ρr^2
    end
    Aexc += Ar * 1e3
    return Aexc * fx * total_mol
end

#####
##### CubicEoS interface for chemical potential
#####

function CubicEoS.log_c_activity!(
    log_ca::AbstractVector,
    mix::MBWREoSMixture,
    nmol::AbstractVector,
    volume,
    RT;
    buf = thermo_buffer(mix, nmol),
)
    Aexc(nmol) = exc_helmholtz(mix, nmol, volume, RT; buf)
    log_ca .= ForwardDiff.gradient(Aexc, nmol) ./ RT
    return log_ca
end

function CubicEoS.log_c_activity_wj!(
    log_ca::AbstractVector,
    jacobian::AbstractMatrix,
    mix::MBWREoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::MBWRThermoBuffer=thermo_buffer(mix, nmol),
)
    Aexc(nmol) = exc_helmholtz(mix, nmol, volume, RT; buf)
    log_ca .= ForwardDiff.gradient(Aexc, nmol) ./ RT
    jacobian .= ForwardDiff.hessian(Aexc, nmol) ./ RT
    return log_ca, jacobian
end
