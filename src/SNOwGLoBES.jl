"""
Supernova neutrino event rates with GLoBES.jl: the detectors, channels,
cross sections and smearing of SNOwGLoBES, fluences (tables, pinched spectra,
the time-dependent Garching model), MSW flavour conversion with Earth crossing,
and CEvNS cross sections.
"""
module SNOwGLoBES

using GLoBES
using GLoBES: Setup
import GLoBES: event_rates, total_events, cross_section
using LazyArtifacts
using LinearAlgebra
using Printf
using SpecialFunctions: gamma

export data_directory,
    Detector, detectors, detector, Binning, Channel, channels, channel_lists, interaction
export Fluence,
    fluence,
    fluences,
    write_fluence,
    regrid,
    number_of_neutrinos,
    mean_energy,
    pinched_spectrum,
    pinched_fluence,
    garching_fluences
export oscillate
export Isotope,
    helm_form_factor,
    klein_nystrand_form_factor,
    cevns_differential,
    cevns_cross_section,
    cevns_cross_sections,
    max_recoil,
    CrossSection,
    cross_section
export SupernovaDetector, ChannelRates, SupernovaRates, event_rates, total_events

include("data.jl")
include("fluence.jl")
include("oscillate.jl")
include("cevns.jl")
include("rates.jl")

end
