# Coherent elastic neutrino–nucleus scattering (CEvNS) and cross-section tables

const G_F = 1.1663788e-11        # MeV⁻²
const ħc_MeV_fm = 197.3269804    # MeV fm
const ħc_MeV_cm = 1.973269804e-11
const AMU_MeV = 931.49410242

"""
    Isotope(Z, N, mass, abundance)

An isotope of a target: its numbers of protons and neutrons, its mass (MeV) and
its abundance (molar fraction).
"""
struct Isotope
    Z::Int
    N::Int
    mass::Float64
    abundance::Float64
end

isotope(Z, A, mass_u, abundance) = Isotope(Z, A - Z, mass_u * AMU_MeV - Z * 0.51099895, abundance)

# Natural abundances (IUPAC) and atomic masses (AME2020, u)
const NUCLEI = Dict(
    :Ar => [
        isotope(18, 36, 35.967545, 0.003336),
        isotope(18, 38, 37.962732, 0.000629),
        isotope(18, 40, 39.962383, 0.996035),
    ],
    :Ge => [
        isotope(32, 70, 69.924249, 0.2057),
        isotope(32, 72, 71.922076, 0.2745),
        isotope(32, 73, 72.923459, 0.0775),
        isotope(32, 74, 73.921178, 0.3650),
        isotope(32, 76, 75.921403, 0.0773),
    ],
    :Xe => [
        isotope(54, 124, 123.905892, 0.000952),
        isotope(54, 126, 125.904298, 0.000890),
        isotope(54, 128, 127.903531, 0.019102),
        isotope(54, 129, 128.904781, 0.264006),
        isotope(54, 130, 129.903509, 0.040710),
        isotope(54, 131, 130.905084, 0.212324),
        isotope(54, 132, 131.904155, 0.269086),
        isotope(54, 134, 133.905395, 0.104357),
        isotope(54, 136, 135.907214, 0.088573),
    ],
)

nucleus(name::Symbol) = get(NUCLEI, name) do
    throw(ArgumentError("unknown nucleus $name; available: $(join(sort(collect(keys(NUCLEI))), ", "))"))
end

j1(x) = x < 1e-3 ? x / 3 : (sin(x) / x - cos(x)) / x   # spherical Bessel function

"""
    helm_form_factor(q, A; s=0.9, radius=nothing) -> Float64

The Helm form factor at momentum transfer `q` (fm⁻¹) for mass number `A`:
F = 3 j₁(q rₙ)/(q rₙ) exp(-(q s)²/2), with the skin thickness `s` (fm) and the
radius rₙ of Lewin and Smith (Astropart. Phys. 6, 87),
rₙ² = c² + 7/3 π² a² - 5 s², c = 1.23 A^(1/3) - 0.6 fm, a = 0.52 fm, unless
`radius` (fm) is given.
"""
function helm_form_factor(q, A; s = 0.9, radius = nothing)
    c = 1.23 * A^(1 / 3) - 0.6
    rn = radius === nothing ? sqrt(c^2 + 7 / 3 * π^2 * 0.52^2 - 5s^2) : radius
    x = q * rn
    return x == 0 ? 1.0 : 3j1(x) / x * exp(-(q * s)^2 / 2)
end

"""
    klein_nystrand_form_factor(q, A; ak=0.7) -> Float64

The Klein–Nystrand form factor, F = 3 j₁(q R)/(q R) / (1 + q² a_k²) with
R = 1.2 A^(1/3) fm and the Yukawa range `ak` (fm).
"""
function klein_nystrand_form_factor(q, A; ak = 0.7)
    x = q * 1.2 * A^(1 / 3)
    return (x == 0 ? 1.0 : 3j1(x) / x) / (1 + (q * ak)^2)
end

form_factor_function(ff::Function) = ff
function form_factor_function(ff::Symbol)
    ff === :helm && return helm_form_factor
    ff === :klein_nystrand && return klein_nystrand_form_factor
    ff === :none && return (q, A) -> 1.0
    throw(
        ArgumentError(
            "unknown form factor $ff; use :helm, :klein_nystrand, :none or a function (q, A) -> F",
        ),
    )
end

"""
    cevns_differential(E, T, target; form_factor=:helm, sin2θW=0.23857) -> Float64

The CEvNS cross section dσ/dT (cm² per MeV, per nucleus) of a neutrino of energy
`E` (MeV) giving the nucleus a recoil energy `T` (MeV):

dσ/dT = G_F² M / 4π · Q_W² (1 - T/E - M T/2E² + T²/2E²) F²(q),  q² = 2 M T,

with Q_W = N - (1 - 4 sin²θ_W) Z. `target` is an `Isotope`, or `:Ar`, `:Ge` or
`:Xe` for the average over the natural abundances. `form_factor` is `:helm`,
`:klein_nystrand`, `:none` or a function `(q, A) -> F` of q in fm⁻¹.
"""
function cevns_differential(E, T, iso::Isotope; form_factor = :helm, sin2θW = 0.23857)
    M = iso.mass
    kin = 1 - T / E - M * T / (2E^2) + T^2 / (2E^2)
    kin > 0 || return 0.0
    QW = iso.N - (1 - 4sin2θW) * iso.Z
    q = sqrt(2M * T) / ħc_MeV_fm
    F = form_factor_function(form_factor)(q, iso.Z + iso.N)
    return G_F^2 * M / (4π) * QW^2 * kin * F^2 * ħc_MeV_cm^2
end
cevns_differential(E, T, target::Symbol; kw...) =
    sum(i.abundance * cevns_differential(E, T, i; kw...) for i in nucleus(target))

"""
    max_recoil(E, M) -> Float64

The largest recoil energy (MeV) of a nucleus of mass `M` hit by a neutrino of
energy `E` (MeV): 2E²/(M + 2E).
"""
max_recoil(E, M) = 2E^2 / (M + 2E)

"""
    cevns_cross_section(E, target; form_factor=:helm, sin2θW=0.23857, recoil_max=Inf) -> Float64

The CEvNS cross section (cm² per nucleus) at neutrino energy `E` (MeV): dσ/dT
integrated over the recoil energies up to `recoil_max` (MeV; all of them by
default). `target` and the keywords as in [`cevns_differential`](@ref).

The tables of SNOwGLoBES are cross sections for recoils below 10 keV, the
window of its CEvNS smearing matrices; for argon above ~15 MeV that is a small
part of the total.
"""
function cevns_cross_section(E, iso::Isotope; recoil_max = Inf, n = 400, kw...)
    Tmax = min(max_recoil(E, iso.mass), recoil_max)
    Tmax > 0 || return 0.0
    # Simpson's rule on n intervals
    h = Tmax / n
    s = cevns_differential(E, 0.0, iso; kw...) + cevns_differential(E, Tmax, iso; kw...)
    for k = 1:(n - 1)
        s += (isodd(k) ? 4 : 2) * cevns_differential(E, k * h, iso; kw...)
    end
    return s * h / 3
end
cevns_cross_section(E, target::Symbol; kw...) =
    sum(i.abundance * cevns_cross_section(E, i; kw...) for i in nucleus(target))

"""
    CrossSection(energies, sigma)

A cross-section table: total cross sections (cm² per target) at `energies`
(MeV), one column of `sigma` per flavour ν_e, ν_μ, ν_τ, ν̄_e, ν̄_μ, ν̄_τ.
"""
struct CrossSection
    energies::Vector{Float64}
    sigma::Matrix{Float64}
    function CrossSection(energies::AbstractVector, sigma::AbstractMatrix)
        size(sigma) == (length(energies), 6) ||
            throw(ArgumentError("sigma must have one row per energy and six columns"))
        return new(collect(Float64, energies), Matrix{Float64}(sigma))
    end
end

Base.show(io::IO, x::CrossSection) = print(
    io,
    "CrossSection($(length(x.energies)) energies, $(fmt(x.energies[1]))–$(fmt(x.energies[end])) MeV)",
)

"""
    cross_section(channel) -> CrossSection

The cross-section table of a SNOwGLoBES channel (by name, e.g. `"nue_Ar40"`).
SNOwGLoBES tabulates σ/E in 10⁻³⁸ cm²/GeV at log₁₀(E/GeV); this is σ in cm².
"""
function cross_section(name::Union{Symbol, AbstractString})
    path = data_file("xscns", "xs_$name.dat")
    isfile(path) || throw(ArgumentError("no cross section for channel $name"))
    t = permutedims(
        reduce(hcat, [parse.(Float64, l[1:7]) for l in data_lines(path) if length(l) >= 7]),
    )
    E_GeV = exp10.(t[:, 1])
    return CrossSection(1000E_GeV, t[:, 2:7] .* E_GeV .* 1e-38)
end

"""
    cevns_cross_sections(target; energies, kw...) -> CrossSection

A CEvNS cross-section table for `target` (`:Ar`, `:Ge`, `:Xe`), the same for all
flavours, for [`SupernovaDetector`](@ref)'s `cross_sections`. The energies
(MeV) are log-spaced from 0.5 to 200 MeV by default; the keywords are those of
[`cevns_cross_section`](@ref).

The CEvNS smearing matrices of SNOwGLoBES spread the events over the recoil
energies of the channel binning only, so give the same `recoil_max`:

```julia
chans = channels(:xenon_coh)
xs = Dict(c.name => cevns_cross_sections(:Xe; recoil_max = c.binning.emax) for c in chans)
SupernovaDetector(:lz, chans; cross_sections = xs)
```
"""
function cevns_cross_sections(
    target::Symbol;
    energies = exp10.(range(log10(0.5), log10(200), 1001)),
    kw...,
)
    σ = [cevns_cross_section(E, target; kw...) for E in energies]
    return CrossSection(energies, repeat(σ, 1, 6))
end

# The table in GLoBES format: log10(E/GeV), σ/E in 10⁻³⁸ cm²/GeV
function write_cross_section(path::AbstractString, x::CrossSection)
    open(path, "w") do io
        for (k, E) in enumerate(x.energies)
            E > 0 || continue
            E_GeV = E / 1000
            vals = x.sigma[k, :] ./ E_GeV ./ 1e-38
            println(
                io,
                @sprintf("%.8f", log10(E_GeV)),
                "  ",
                join((@sprintf("%.8e", v) for v in vals), "  "),
            )
        end
    end
    return path
end
