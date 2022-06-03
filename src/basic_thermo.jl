function __θ__(Tr, acentric)
    a, b = 0.05202976, -0.7498189
    return 1 + (acentric - C3_OMEGA) * (a + b * log(Tr))
end

function __ϕ__(Tr, Zc, acentric)
    a, b = 0.1435971, -0.2821562
    return (1 - (acentric - C3_OMEGA) * (a + b * log(Tr))) * C3_ZC / Zc
end

function scaling_coeffs(chx::MBWREoSComponent, T)
    acentric = acentric_factor(chx)
    Pc, RTc, Zc = chx.Pc, chx.RTc, chx.Zc
    Vc = Zc * RTc / Pc
    Tc = RTc / GAS_CONSTANT_SI
    Tr = T / Tc
    fx = (Tc / C3_TC) * __θ__(Tr, acentric)
    hx = (Vc / C3_VC) * __ϕ__(Tr, Zc, acentric)
    return (; fx, hx)
end

function __shape_factors!__(mix::MBWREoSMixture, T, fs, hs)
    comps = components(mix)
    @inbounds for i in eachindex(comps, fs, hs)
        subst = comps[i]
        Tc = subst.RTc / GAS_CONSTANT_SI
        Zc = subst.Zc
        Vc = Zc * subst.RTc / subst.Pc
        acentric = acentric_factor(subst)
        fs[i] = (Tc / C3_TC) .* __θ__(T / Tc, acentric)
        hs[i] = (Vc / C3_VC) .* __ϕ__(T / Tc, Zc, acentric)
    end
    return (f = fs, h = hs)
end

function __shape_matrices!__(fs, hs, fmatr, hmatr)
    fmatr .= sqrt.(fs .* fs')
    hmatr .= 0.125 .* (cbrt.(hs) .+ cbrt.(hs')).^3
    return (fm = fmatr, hm = hmatr)
end

"""
    scaling_coeffs(mixture::MBWREoSMixture, nmol, T[; buf])

Return scaling coefficients `fx`, `hx` and `Mx` for the TRAPP EoS for `mixture` with
composition `nmol` and thermal energy `RT`. Allocations may be avoided by passing `buf`.

# Arguments
- `nmol::AbstractVector`: Vector of composition
- `RT::Real`: Thermal energy (J mol⁻¹)

# Keywords
- `buf`: see ?pressure(mixture)
"""
function scaling_coeffs(mix::MBWREoSMixture, nmol::AbstractVector, T;
    buf::MBWRThermoBuffer=thermo_buffer(mix, nmol)
)
    fs, hs = __shape_factors!__(mix, T, buf.vec1, buf.vec2)
    fm, hm = __shape_matrices!__(fs, hs, buf.fij, buf.hij)
    comp = components(mix)
    # TODO: use buf.vec2
    # For some reasons buf.vec2 .= nmol ./ sum(nmol) don't work w ForwardDiff
    x = nmol ./ sum(nmol)

    hx = sum(x[i] * x[j] * hm[i,j] for i in eachindex(comp), j in eachindex(comp))
    fx = sum(
        x[i] * x[j] * hm[i,j] * fm[i,j] for i in eachindex(comp), j in eachindex(comp)
    ) / hx
    return (; fx, hx)
end

function __ref_pressure__(ρ, T)
    ρmoll = ρ * 1e-3
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
    iTsq = inv(T^2)
    iT = inv(T)
    a1 = GAS_CONSTANT_SI * T * 1e-3
    a2 = G1 * T + G2 * sqrt(T) + G3 + (G4 * T + G5) * iTsq
    a3 = G6 * T + G7 + (G8 * T + G9) * iTsq
    a4 = G10 * T + G11 + G12 * iT
    a5 = G13
    a6 = (G14 * T + G15) * iTsq
    a7 = G16 * iT
    a8 = (G17 * T + G18) * iTsq
    a9 = G19 * iTsq

    b3 = (G20 * T + G21) * iT * iTsq
    b5 = (G22 * T^2 + G23) * iTsq^2
    b7 = (G24 * T + G25) * iT * iTsq
    b9 = (G26 * T^2 + G27) * iTsq^2
    b11 = (G28 * T + G29) * iT * iTsq
    b13 = (G32 + T * (G31 + T * G30)) * iTsq^2

    a = (a1, a2, a3, a4, a5, a6, a7, a8, a9)
    b = (b3, b5, b7, b9, b11, b13)
    p_MPa = evalpoly(ρmoll, a) * ρmoll + evalpoly(ρmoll^2, b) * ρmoll^3 * exp(-(C3_VC * ρ)^2)
    return p_MPa * 1e6
end

function CubicEoS.pressure(substance::MBWREoSComponent, nmol::Real, V::Real, RT::Real)
    molvol = V / nmol
    T = RT / GAS_CONSTANT_SI

    fx, hx = scaling_coeffs(substance, T)
    T_ref, v_ref = T / fx, molvol / hx
    ρ_ref = inv(v_ref)
    p_ref = __ref_pressure__(ρ_ref, T_ref)

    return p_ref * fx / hx
end

function CubicEoS.wilson_saturation_pressure(substance::MBWREoSComponent, RT::Real)
    return CubicEoS.wilson_saturation_pressure(
        substance.Pc, substance.RTc, substance.acentric_factor, RT
    )
end

function CubicEoS.pressure(
    mixture::MBWREoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf=thermo_buffer(mixture, nmol),
)
    T = RT / GAS_CONSTANT_SI
    fx, hx = scaling_coeffs(mixture, nmol, T; buf)
    T_ref, v_ref = T / fx, volume / (sum(nmol) * hx)
    ρ_ref = inv(v_ref)
    p_ref = __ref_pressure__(ρ_ref, T_ref)

    return p_ref * fx / hx
end
