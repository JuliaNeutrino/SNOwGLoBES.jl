# The SNOwGLoBES data tables (a lazy artifact) and the detector and channel lists

"""
    data_directory() -> String

The directory of the SNOwGLoBES data tables: `fluxes/`, `xscns/`, `smear/`,
`effic/`, `channels/`, `backgrounds/` and `detector_configurations.dat`. They
are downloaded (34 MB, checksummed) the first time they are needed.
"""
function data_directory()
    root = artifact"snowglobes"
    return joinpath(root, only(readdir(root)))   # the tarball has one top-level directory
end

data_file(parts...) = joinpath(data_directory(), parts...)

# Lines of a data file without comments (starting with `comment`) and blank lines
data_lines(path; comment = "#") =
    [split(s) for s in (strip(l) for l in eachline(path)) if !isempty(s) && !startswith(s, comment)]

"""
    Detector

A detector configuration of SNOwGLoBES: its `name`, `mass` (kt) and
`normalization`, the number of reference targets per atomic mass unit of the
material (2/18 for the free protons of water, 1/40 for ⁴⁰Ar nuclei, ...).
"""
struct Detector
    name::String
    mass::Float64
    normalization::Float64
end

Base.show(io::IO, d::Detector) =
    print(io, "Detector(", d.name, ": ", d.mass, " kt, ", d.normalization, " targets per amu)")

"""
    detectors() -> Vector{Detector}
    detector(name) -> Detector

The detector configurations that come with SNOwGLoBES (from
`detector_configurations.dat`), or one of them by name.

```julia
detector(:ar40kt)   # 40 kt of liquid argon
[d.name for d in detectors()]
```
"""
detectors() = [
    Detector(String(f[1]), parse(Float64, f[2]), parse(Float64, f[3])) for
    f in data_lines(data_file("detector_configurations.dat")) if length(f) >= 3
]

function detector(name::Union{Symbol, AbstractString})
    i = findfirst(d -> d.name == String(name), detectors())
    i === nothing && throw(
        ArgumentError(
            "unknown detector $name; available: $(join((d.name for d in detectors()), ", "))",
        ),
    )
    return detectors()[i]
end
@doc (@doc detectors) detector

"""
    Binning(sampling_points, sampling_min, sampling_max, bins, emin, emax)

The energy grid of a channel, energies in MeV: `sampling_points` true neutrino
energies between `sampling_min` and `sampling_max`, and `bins` detected-energy
bins between `emin` and `emax`. The default is 200 and 200 from 0.5 to 100 MeV.
"""
struct Binning
    sampling_points::Int
    sampling_min::Float64
    sampling_max::Float64
    bins::Int
    emin::Float64
    emax::Float64
end

const DEFAULT_BINNING = Binning(200, 0.5, 100, 200, 0.5, 100)

# Six numbers in GeV, the order of the channel files
binning_from_gev(f) = Binning(
    parse(Int, f[1]),
    1000parse(Float64, f[2]),
    1000parse(Float64, f[3]),
    parse(Int, f[4]),
    1000parse(Float64, f[5]),
    1000parse(Float64, f[6]),
)

"""
    Channel

An interaction channel: its `name` (which names its cross section, smearing and
efficiency files), the neutrino `flavour` (`:e`, `:μ` or `:τ`), `antineutrino`,
the `weight` (number of targets per reference target, e.g. 18 electrons per
⁴⁰Ar nucleus) and its energy `binning`.
"""
struct Channel
    name::String
    flavour::Symbol
    antineutrino::Bool
    weight::Float64
    binning::Binning
end

function Base.show(io::IO, c::Channel)
    nu = (c.antineutrino ? "ν̄" : "ν") * (c.flavour === :e ? "e" : String(c.flavour))
    print(io, "Channel(", c.name, ": ", nu, ", weight ", c.weight, ")")
end

const FLAVOURS = Dict("e" => :e, "m" => :μ, "t" => :τ)

"""
    channels(name) -> Vector{Channel}

The interaction channels of a SNOwGLoBES channel list,
`channels/channels_<name>.dat`: `:argon`, `:water`, `:scint`, `:lead`,
`:heavywater`, `:argon_he`, `:xenon_coh`, ...

A line `% bins...` before the channels sets their binning (the high-energy
lists, up to 200 MeV); lines with the `SN_nu` prefix carry their own binning.

```julia
channels(:argon)   # ν-e elastic scattering, CC and NC on ⁴⁰Ar
```
"""
function channels(name::Union{Symbol, AbstractString})
    path = data_file("channels", "channels_$name.dat")
    isfile(path) || throw(
        ArgumentError("unknown channel list $name; available: $(join(channel_lists(), ", "))"),
    )
    binning = DEFAULT_BINNING
    list = Channel[]
    for line in eachline(path)
        f = split(line)
        isempty(f) && continue
        if startswith(f[1], "%")
            # SNOwGLoBES skips this line, which leaves the high-energy channels on the 100 MeV grid
            nums = filter(!isnothing, tryparse.(Float64, f[2:end]))
            length(nums) == 6 && (binning = binning_from_gev(f[2:7]))
            continue
        end
        f[1] == "SN_nu" && (f = f[2:end])
        b = length(f) >= 11 ? binning_from_gev(f[6:11]) : binning
        push!(list, Channel(f[1], FLAVOURS[f[4]], f[3] == "-", parse(Float64, f[5]), b))
    end
    return list
end

"""
    channel_lists() -> Vector{String}

The names of the channel lists that come with SNOwGLoBES, for [`channels`](@ref).
"""
channel_lists() =
    sort([f[10:(end - 4)] for f in readdir(data_file("channels")) if startswith(f, "channels_")])

"""
    interaction(channel) -> Symbol

The kind of interaction of a channel, from its name: `:ibd` (inverse beta
decay), `:es` (elastic scattering on electrons), `:nc` (neutral current on
nuclei), `:cevns` (coherent elastic scattering on nuclei) or `:cc` (charged
current on nuclei).
"""
function interaction(name::AbstractString)
    startswith(name, "ibd") && return :ibd
    startswith(name, "coh") && return :cevns
    startswith(name, "nc") && return :nc
    occursin(r"^nu[a-z]*_e(_he)?$", name) && return :es
    return :cc
end
interaction(c::Channel) = interaction(c.name)
