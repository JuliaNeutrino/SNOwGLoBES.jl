# Supernova neutrino fluences: tables, the pinched (Garching) parameterisation
# and the time-dependent Garching model

const KPC_CM = 3.0856775814913673e21   # 1 kpc in cm
const MEV_PER_ERG = 624150.907446076
const FLAVOUR_NAMES = (:nue, :numu, :nutau, :nuebar, :numubar, :nutaubar)

# SNOwGLoBES and GLoBES tabulate fluences per 0.2 MeV ("per energy bin")
const FILE_BIN_MEV = 0.2

"""
    Fluence(energies, spectra)
    Fluence(energies; nue, numu, nutau, nuebar, numubar, nutaubar)

The time-integrated flux of a supernova at the detector, dΦ/dE in neutrinos per
cm² per MeV, at `energies` (MeV) for the six flavours ν_e, ν_μ, ν_τ, ν̄_e,
ν̄_μ, ν̄_τ (the columns of `spectra`). Read the flavours as `f.nue`, ...,
`f.nutaubar`.

Fluences add and scale: `f1 + f2`, `2f`.
"""
struct Fluence
    energies::Vector{Float64}
    spectra::Matrix{Float64}
    function Fluence(energies::AbstractVector, spectra::AbstractMatrix)
        size(spectra) == (length(energies), 6) ||
            throw(ArgumentError("spectra must have one row per energy and six columns"))
        issorted(energies) || throw(ArgumentError("energies must be increasing"))
        return new(collect(Float64, energies), Matrix{Float64}(spectra))
    end
end

Fluence(energies::AbstractVector; nue, numu, nutau, nuebar, numubar, nutaubar) =
    Fluence(energies, hcat(nue, numu, nutau, nuebar, numubar, nutaubar))

function Base.getproperty(f::Fluence, name::Symbol)
    i = findfirst(==(name), FLAVOUR_NAMES)
    i === nothing || return getfield(f, :spectra)[:, i]
    return getfield(f, name)
end
Base.propertynames(::Fluence) = (:energies, :spectra, FLAVOUR_NAMES...)

function same_grid(a::Fluence, b::Fluence)
    a.energies == b.energies ||
        throw(ArgumentError("fluences on different energy grids; use `regrid`"))
    return a.energies
end
Base.:+(a::Fluence, b::Fluence) = Fluence(same_grid(a, b), a.spectra + b.spectra)
Base.:*(x::Real, f::Fluence) = Fluence(f.energies, x * f.spectra)
Base.:*(f::Fluence, x::Real) = x * f
Base.:(==)(a::Fluence, b::Fluence) = a.energies == b.energies && a.spectra == b.spectra
Base.isapprox(a::Fluence, b::Fluence; kw...) =
    a.energies == b.energies && isapprox(a.spectra, b.spectra; kw...)

function Base.show(io::IO, f::Fluence)
    n = number_of_neutrinos(f)
    print(
        io,
        "Fluence($(length(f.energies)) energies, $(f.energies[1])–$(f.energies[end]) MeV; ",
        "ν_e ",
        fmt(n.nue),
        ", ν̄_e ",
        fmt(n.nuebar),
        ", ν_μ ",
        fmt(n.numu),
        " per cm²)",
    )
end
fmt(x) = @sprintf("%.3g", x)

# Trapezoidal integral of y over x
trapz(x, y) = sum((x[i + 1] - x[i]) * (y[i] + y[i + 1]) / 2 for i = 1:(length(x) - 1); init = 0.0)

"""
    number_of_neutrinos(f) -> NamedTuple
    mean_energy(f) -> NamedTuple

The number of neutrinos per cm² of each flavour (`nue`, ..., `nutaubar`), and
their mean energies (MeV).
"""
number_of_neutrinos(f::Fluence) =
    NamedTuple{FLAVOUR_NAMES}(Tuple(trapz(f.energies, f.spectra[:, i]) for i = 1:6))

function mean_energy(f::Fluence)
    n = number_of_neutrinos(f)
    return NamedTuple{FLAVOUR_NAMES}(
        Tuple(trapz(f.energies, f.energies .* f.spectra[:, i]) / n[i] for i = 1:6),
    )
end
@doc (@doc number_of_neutrinos) mean_energy

"""
    regrid(f, energies) -> Fluence

The fluence interpolated linearly to new energies (MeV), zero outside the
original ones.
"""
function regrid(f::Fluence, energies::AbstractVector)
    out = zeros(length(energies), 6)
    for (k, E) in enumerate(energies)
        f.energies[1] <= E <= f.energies[end] || continue
        i = clamp(searchsortedlast(f.energies, E), 1, length(f.energies) - 1)
        w = (E - f.energies[i]) / (f.energies[i + 1] - f.energies[i])
        out[k, :] .= (1 - w) .* f.spectra[i, :] .+ w .* f.spectra[i + 1, :]
    end
    return Fluence(energies, out)
end

"""
    fluences() -> Vector{String}

The names of the fluence tables that come with SNOwGLoBES, for [`fluence`](@ref):
`livermore` (Totani et al. 1998), `gvkm` (Gava et al. 2009, with collective
oscillations) and others.
"""
fluences() = sort([
    f[1:(end - 4)] for
    f in readdir(data_file("fluxes")) if endswith(f, ".dat") && !occursin("info", f)
])

"""
    fluence(name; distance=10) -> Fluence
    fluence(path; distance=10) -> Fluence

A fluence table of SNOwGLoBES (`:livermore`, `:gvkm`, ..., see [`fluences`](@ref)),
or a file in its format: energy in GeV, then the six flavours in neutrinos per
cm² per 0.2 MeV, for a supernova 10 kpc away. `distance` (kpc) rescales the
fluence to another distance.

```julia
fluence(:livermore)   # the Livermore model at 10 kpc
fluence(:livermore; distance = 1)   # the same supernova at 1 kpc: 100 times more neutrinos
```
"""
function fluence(name::Union{Symbol, AbstractString}; distance::Real = 10)
    path = isfile(String(name)) ? String(name) : data_file("fluxes", "$name.dat")
    isfile(path) || throw(ArgumentError("no fluence $name; available: $(join(fluences(), ", "))"))
    rows = [parse.(Float64, l[1:7]) for l in data_lines(path) if length(l) >= 7]
    table = permutedims(reduce(hcat, rows))
    return (10 / distance)^2 * Fluence(1000table[:, 1], table[:, 2:7] / FILE_BIN_MEV)
end

# The fluence as a GLoBES flux table: energy (GeV) and neutrinos per cm² per 0.2 MeV
flux_table(f::Fluence) = hcat(f.energies / 1000, f.spectra * FILE_BIN_MEV)

"""
    write_fluence(path, f)

Write a fluence in the format of the SNOwGLoBES flux files.
"""
function write_fluence(path::AbstractString, f::Fluence)
    open(path, "w") do io
        for row in eachrow(flux_table(f))
            println(io, join((@sprintf("%.8e", x) for x in row), "  "))
        end
    end
    return path
end

"""
    pinched_spectrum(E, mean_energy, α) -> Float64

The quasi-thermal ("pinched", Garching) spectrum, normalised to one, per MeV:
φ(E) = (α+1)^(α+1) / (⟨E⟩ Γ(α+1)) (E/⟨E⟩)^α exp(-(α+1)E/⟨E⟩).
α = 2 is a Fermi–Dirac-like thermal spectrum; larger α pinches it.
"""
function pinched_spectrum(E, mean_energy, α)
    mean_energy > 0 || return 0.0
    x = E / mean_energy
    return (α + 1)^(α + 1) / gamma(α + 1) * x^α * exp(-(α + 1) * x) / mean_energy
end

"""
    pinched_fluence(; mean_energy, α, total_energy, distance=10, energies=0:0.2:200) -> Fluence

The fluence of a supernova whose spectra are pinched (see [`pinched_spectrum`](@ref)),
from the mean energies (MeV), pinching parameters `α` and `total_energy` emitted
(erg) of ν_e, ν̄_e and each ν_x (ν_μ, ν_τ, ν̄_μ, ν̄_τ), as three-element tuples
in this order. `distance` is in kpc.

This is the `pinched` tool of SNOwGLoBES. No flavour conversion: apply
[`oscillate`](@ref) for it.

```julia
pinched_fluence(mean_energy = (12, 15, 18), α = (2, 3, 4), total_energy = (5e52, 5e52, 5e52))
```
"""
function pinched_fluence(;
    mean_energy,
    α,
    total_energy,
    distance::Real = 10,
    energies::AbstractVector = 0:0.2:200,
)
    length(mean_energy) == length(α) == length(total_energy) == 3 ||
        throw(ArgumentError("give mean_energy, α and total_energy for ν_e, ν̄_e and ν_x"))
    area = 4π * (distance * KPC_CM)^2
    function spectrum(k)
        (mean_energy[k] > 0 && total_energy[k] > 0) || return zeros(length(energies))
        n = total_energy[k] * MEV_PER_ERG / mean_energy[k] / area   # neutrinos per cm²
        return [n * pinched_spectrum(E, mean_energy[k], α[k]) for E in energies]
    end
    e, ebar, x = spectrum(1), spectrum(2), spectrum(3)
    return Fluence(energies; nue = e, numu = x, nutau = x, nuebar = ebar, numubar = x, nutaubar = x)
end

"""
    garching_fluences(; distance=10, energies=0:0.2:200) -> (; time, duration, fluences)

The time-dependent Garching model of SNOwGLoBES (Hüdepohl et al.,
Phys. Rev. Lett. 104, 251101): the start `time` (s, after the bounce) and the
`duration` (s) of each time step, and its pinched fluence.

```julia
g = garching_fluences()
sum(g.fluences)   # the time-integrated fluence
```
"""
function garching_fluences(; distance::Real = 10, energies::AbstractVector = 0:0.2:200)
    dir = data_file("fluxes", "garching")
    info = [parse.(Float64, l) for l in data_lines(joinpath(dir, "garching_pinched_info.dat"))]
    key = Dict(
        Int(t[1]) => (t[2], t[3]) for t in
        (parse.(Float64, l) for l in data_lines(joinpath(dir, "garching_pinched_info_key.dat")))
    )
    steps = [key[Int(r[1])] for r in info]
    fl = [
        pinched_fluence(;
            α = r[2:4],
            mean_energy = r[5:7],
            total_energy = r[8:10],
            distance,
            energies,
        ) for r in info
    ]
    return (time = first.(steps), duration = last.(steps), fluences = fl)
end
