# Event rates: a GLoBES experiment per detector, channel list and binning

# The AEDL settings of SNOwGLoBES (glb/preamble.glb, flux.glb, detector.glb):
# rates = flux (per cm² per 0.2 MeV) × σ × target mass × FLUX_NORM / BASELINE²
const BASELINE = 1.5      # km
const FLUX_NORM = 0.0677

"""
    SupernovaDetector(detector, channels; mass, background, cross_sections)

A detector configuration and a channel list of SNOwGLoBES, as GLoBES
experiments ready to compute event rates for any fluence with
[`event_rates`](@ref). Building it reads the cross sections, smearing matrices
and efficiencies once.

- `detector`: a name of [`detectors`](@ref) (`:ar40kt`, `:wc100kt30prct`, ...).
- `channels`: a name of [`channel_lists`](@ref) (`:argon`, `:water`, ...) or a
  vector of [`Channel`](@ref)s.
- `mass`: the detector mass in kt (default: that of the configuration).
- `background`: `"surf"`, `"0800"` or `"4850"` (depth in feet) for the
  backgrounds of SNOwGLoBES, or a vector of events per sampling bin; for
  detectors with a background smearing matrix (`ar40kt`).
- `cross_sections`: a `Dict` from channel names to [`CrossSection`](@ref)s,
  replacing the tables of SNOwGLoBES (e.g. [`cevns_cross_sections`](@ref)).

```julia
dune = SupernovaDetector(:ar40kt, :argon)
event_rates(dune, fluence(:livermore))
```
"""
struct SupernovaDetector
    detector::Detector
    channels::Vector{Channel}
    setup::Setup
    experiments::Vector{Experiment}
    experiment_of::Vector{Int}          # the experiment of each channel
    background::Union{Nothing, Int}     # the experiment with the background channel
    dir::String
end

function SupernovaDetector(
    det::Union{Symbol, AbstractString},
    chans::Union{Symbol, AbstractString, AbstractVector{Channel}};
    mass::Union{Nothing, Real} = nothing,
    background = nothing,
    cross_sections::AbstractDict = Dict{String, CrossSection}(),
)
    d = detector(det)
    mass === nothing || (d = Detector(d.name, mass, d.normalization))
    list = chans isa AbstractVector ? collect(Channel, chans) : channels(chans)
    isempty(list) && throw(ArgumentError("no channels"))
    xs = Dict(String(k) => v for (k, v) in cross_sections)
    unknown = setdiff(keys(xs), (c.name for c in list))
    isempty(unknown) ||
        throw(ArgumentError("cross sections for channels not in the list: $(join(unknown, ", "))"))
    check_files(d, list, xs)

    dir = mktempdir()
    flux_file = write_fluence(joinpath(dir, "fluence.dat"), Fluence([0.0, 1.0], zeros(2, 6)))
    for (name, x) in xs
        write_cross_section(joinpath(dir, "xs_$name.dat"), x)
    end
    bg = background_rates(d, background)
    setup = Setup()
    experiments = Experiment[]
    # one experiment per binning; the background goes with the default binning
    groups = unique(c.binning for c in list)
    bg === nothing || DEFAULT_BINNING in groups || push!(groups, DEFAULT_BINNING)
    for b in groups
        members = [c for c in list if c.binning == b]
        with_bg = bg !== nothing && b == DEFAULT_BINNING
        file = joinpath(dir, "experiment$(length(experiments) + 1).glb")
        write(file, aedl(d, members, b, flux_file, dir, xs, with_bg ? bg : nothing))
        push!(experiments, Experiment(file; setup))
    end
    experiment_of = [findfirst(==(c.binning), groups) for c in list]
    bg_experiment = bg === nothing ? nothing : findfirst(==(DEFAULT_BINNING), groups)
    return SupernovaDetector(d, list, setup, experiments, experiment_of, bg_experiment, dir)
end

function Base.show(io::IO, s::SupernovaDetector)
    print(io, "SupernovaDetector(", s.detector.name, ", ", s.detector.mass, " kt, ")
    print(
        io,
        length(s.channels),
        " channels",
        s.background === nothing ? "" : " and a background",
        ")",
    )
end

xsec_path(name, dir, xs) =
    haskey(xs, name) ? joinpath(dir, "xs_$name.dat") : data_file("xscns", "xs_$name.dat")
smear_path(name, d) = data_file("smear", d.name, "smear_$(name)_$(d.name).dat")
effic_path(name, d) = data_file("effic", d.name, "effic_$(name)_$(d.name).dat")

function check_files(d::Detector, list, xs)
    missing_files = String[]
    for c in list
        haskey(xs, c.name) ||
            isfile(xsec_path(c.name, "", xs)) ||
            push!(missing_files, "xscns/xs_$(c.name).dat")
        isfile(smear_path(c.name, d)) ||
            push!(missing_files, "smear/$(d.name)/smear_$(c.name)_$(d.name).dat")
    end
    isempty(missing_files) || throw(
        ArgumentError("$(d.name) cannot use these channels; missing $(join(missing_files, ", "))"),
    )
end

function background_rates(d::Detector, background)
    background === nothing && return nothing
    isfile(smear_path("bg_chan", d)) || throw(
        ArgumentError(
            "$(d.name) has no background smearing matrix (smear/$(d.name)/smear_bg_chan_$(d.name).dat)",
        ),
    )
    background isa AbstractVector && return collect(Float64, background)
    path = data_file("backgrounds", "bg_$background.dat")
    isfile(path) ||
        throw(ArgumentError("no background $background; use \"surf\", \"0800\" or \"4850\""))
    return parse.(Float64, split(strip(read(path, String), ['{', '}', ' ', '\n']), ','))
end

aedl_list(v) = "{" * join(v, ", ") * "}"
gev(x) = x / 1000

function aedl(d::Detector, members, b::Binning, flux_file, dir, xs, bg)
    io = IOBuffer()
    print(
        io,
        """
        %!GLoBES
        \$version = "3.0.0"
        \$profiletype = 1
        \$baseline = $BASELINE
        \$bins = $(b.bins)
        \$emin = $(gev(b.emin))
        \$emax = $(gev(b.emax))
        \$sampling_points = $(b.sampling_points)
        \$sampling_min = $(gev(b.sampling_min))
        \$sampling_max = $(gev(b.sampling_max))
        nuflux(#supernova_flux)<
            @flux_file = "$flux_file"
            @time = 1
            @power = 1
            @norm = $FLUX_NORM
        >
        \$target_mass = $(d.mass * d.normalization)
        """,
    )
    names = [c.name for c in members]
    bg === nothing || push!(names, "bg_chan")
    for n in names
        println(io, "include \"", smear_path(n, d), "\"")
    end
    for c in members
        println(io, "cross(#$(c.name))< @cross_file = \"$(xsec_path(c.name, dir, xs))\" >")
    end
    bg === nothing ||
        println(io, "cross(#bg_chan)< @cross_file = \"$(data_file("xscns", "xs_zero.dat"))\" >")
    for c in members
        sign = c.antineutrino ? "-" : "+"
        fl = Dict(:e => "e", :μ => "m", :τ => "t")[c.flavour]
        println(io, "channel(#$(c.name)_signal)<")
        println(io, "    @channel = #supernova_flux: $sign: $fl: $fl: #$(c.name): #$(c.name)_smear")
        # SNOwGLoBES (run mode 0) misses these files: it names them with an unset variable
        eff = effic_path(c.name, d)
        isfile(eff) && println(io, "    @post_smearing_efficiencies = ", strip(read(eff, String)))
        println(io, ">")
    end
    if bg !== nothing
        println(io, "channel(#bg_chan_signal)<")
        println(io, "    @channel = #supernova_flux: -: e: e: #bg_chan: #bg_chan_smear")
        println(io, "    @pre_smearing_background = ", aedl_list(bg))
        println(io, ">")
    end
    print(
        io,
        """
        rule(#supernova_events)<
            @signal = 1.0@#$(members[1].name)_signal
            @signalerror = 0.011 : 0.00005
            @background = 0.0@#$(members[1].name)_signal
            @backgrounderror = 0.011 : 0.00005
            @sys_on_function = "chiSpectrumTilt"
            @sys_off_function = "chiNoSysSpectrum"
            @energy_window = $(gev(b.emin)) : $(gev(b.emax))
        >
        """,
    )
    return String(take!(io))
end

"""
    ChannelRates

The events of one channel: `energies` (MeV, the true neutrino energies) and
`unsmeared` (events per sampling bin, before the detector response), and
`detected_energies` (MeV) and `smeared` (events per detected-energy bin, with
the smearing and the efficiencies). The channel weight is applied.
"""
struct ChannelRates
    name::String
    energies::Vector{Float64}
    unsmeared::Vector{Float64}
    detected_energies::Vector{Float64}
    smeared::Vector{Float64}
end

Base.show(io::IO, r::ChannelRates) = print(
    io,
    "ChannelRates(",
    r.name,
    ": ",
    fmt(sum(r.unsmeared)),
    " events, ",
    fmt(sum(r.smeared)),
    " detected)",
)

"""
    SupernovaRates

The events of all channels of a [`SupernovaDetector`](@ref) for one fluence.
Index it by channel name (`r["nue_Ar40"]` or `r[:nue_Ar40]`) for the
[`ChannelRates`](@ref); `r.background` holds the background, if any.
"""
struct SupernovaRates
    detector::Detector
    channels::Vector{ChannelRates}
    background::Union{Nothing, ChannelRates}
end

function Base.getindex(r::SupernovaRates, name::Union{Symbol, AbstractString})
    i = findfirst(c -> c.name == String(name), r.channels)
    i === nothing && throw(KeyError(name))
    return r.channels[i]
end
Base.keys(r::SupernovaRates) = [c.name for c in r.channels]
Base.iterate(r::SupernovaRates, s...) = iterate(r.channels, s...)
Base.length(r::SupernovaRates) = length(r.channels)

"""
    total_events(r::SupernovaRates; smeared=true, interaction=nothing) -> Float64
    total_events(r::ChannelRates; smeared=true) -> Float64

The number of events of all channels, or of those of one kind of
[`interaction`](@ref) (`:ibd`, `:es`, `:cc`, `:nc`, `:cevns`), detected
(`smeared=true`) or before the detector response. The background is not
included.
"""
GLoBES.total_events(r::ChannelRates; smeared::Bool = true) = sum(smeared ? r.smeared : r.unsmeared)
function GLoBES.total_events(r::SupernovaRates; smeared::Bool = true, interaction = nothing)
    chans =
        interaction === nothing ? r.channels :
        filter(c -> SNOwGLoBES.interaction(c.name) === interaction, r.channels)
    return sum((total_events(c; smeared) for c in chans); init = 0.0)
end

function Base.show(io::IO, ::MIME"text/plain", r::SupernovaRates)
    println(io, "Supernova events in ", r.detector.name, " (", r.detector.mass, " kt):")
    w = maximum(length(c.name) for c in r.channels)
    println(io, "  ", rpad("channel", w), "  ", lpad("events", 10), "  ", lpad("detected", 10))
    for c in r.channels
        println(
            io,
            "  ",
            rpad(c.name, w),
            "  ",
            lpad(@sprintf("%.2f", sum(c.unsmeared)), 10),
            "  ",
            lpad(@sprintf("%.2f", sum(c.smeared)), 10),
        )
    end
    for kind in (:ibd, :es, :cc, :nc, :cevns)
        any(c -> interaction(c.name) === kind, r.channels) || continue
        label = Dict(:ibd => "IBD", :es => "ES", :cc => "CC", :nc => "NC", :cevns => "CEvNS")[kind]
        println(
            io,
            "  ",
            rpad("total $label", w),
            "  ",
            lpad(@sprintf("%.2f", total_events(r; smeared = false, interaction = kind)), 10),
            "  ",
            lpad(@sprintf("%.2f", total_events(r; interaction = kind)), 10),
        )
    end
    print(
        io,
        "  ",
        rpad("total", w),
        "  ",
        lpad(@sprintf("%.2f", total_events(r; smeared = false)), 10),
        "  ",
        lpad(@sprintf("%.2f", total_events(r)), 10),
    )
    r.background === nothing || print(
        io,
        "\n  ",
        rpad("background", w),
        "  ",
        lpad(@sprintf("%.2f", sum(r.background.unsmeared)), 10),
        "  ",
        lpad(@sprintf("%.2f", sum(r.background.smeared)), 10),
    )
end
Base.show(io::IO, r::SupernovaRates) = print(
    io,
    "SupernovaRates(",
    r.detector.name,
    ": ",
    length(r.channels),
    " channels, ",
    fmt(total_events(r)),
    " events)",
)

const NO_OSCILLATIONS = OscParams(θ12 = 0, θ13 = 0, θ23 = 0, Δm21 = 0, Δm31 = 0)

"""
    event_rates(detector::SupernovaDetector, f::Fluence) -> SupernovaRates
    event_rates(detector, fs::AbstractVector{Fluence}) -> Vector{SupernovaRates}
    event_rates(fluence, channels, detector; kw...) -> SupernovaRates

The supernova events of each channel for a fluence (e.g. a time series). The
fluence is the flux at the detector: apply [`oscillate`](@ref) first for
flavour conversion. The last form builds the [`SupernovaDetector`](@ref) (with
the keywords) and accepts a fluence name.

```julia
event_rates(:livermore, :argon, :ar40kt)
```
"""
function GLoBES.event_rates(s::SupernovaDetector, f::Fluence)
    table = flux_table(f)
    for e in s.experiments
        e.data.fluxes[1].storage = table   # swap the flux table of the loaded experiment
    end
    simulate!(NO_OSCILLATIONS; setup = s.setup)
    function rates(e, name, weight)
        c = findfirst(==("$(name)_signal"), channel_names(e))
        return ChannelRates(
            name,
            e.data.smear_data[1].simbincenter * 1000,
            weight * channel_rates(e, c; smeared = false),
            bin_centers(e) * 1000,
            weight * channel_rates(e, c),
        )
    end
    chans = [
        rates(s.experiments[s.experiment_of[k]], c.name, c.weight) for
        (k, c) in enumerate(s.channels)
    ]
    bg = s.background === nothing ? nothing : rates(s.experiments[s.background], "bg_chan", 1.0)
    return SupernovaRates(s.detector, chans, bg)
end

GLoBES.event_rates(s::SupernovaDetector, fs::AbstractVector{Fluence}) =
    [event_rates(s, f) for f in fs]

function GLoBES.event_rates(f::Union{Fluence, Symbol, AbstractString}, chans, det; kw...)
    fl = f isa Fluence ? f : fluence(f)
    return event_rates(SupernovaDetector(det, chans; kw...), fl)
end
