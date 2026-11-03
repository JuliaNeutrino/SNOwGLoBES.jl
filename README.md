# SNOwGLoBES.jl

Supernova neutrino event rates with [GLoBES.jl](https://github.com/JuliaNeutrino/GLoBES.jl): a Julia port of [SNOwGLoBES](https://github.com/SNOwGLoBES/snowglobes). It uses the detectors, interaction channels, cross sections, smearing matrices and efficiencies of SNOwGLoBES, without its Perl scripts and C driver, and adds:

- fluences from tables, pinched (Garching) spectra and the time-dependent Garching model, in physical units (cm⁻² MeV⁻¹);
- adiabatic MSW flavour conversion in the supernova with full three-flavour mixing, and the Earth crossing;
- CEvNS cross sections with the Helm or Klein–Nystrand form factor;
- one loaded detector for any number of fluences (a 300-step time series takes about a second).

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaNeutrino/SNOwGLoBES.jl")   # not registered yet
using SNOwGLoBES, GLoBES

dune = SupernovaDetector(:ar40kt, :argon)            # 40 kt of liquid argon
f = oscillate(fluence(:livermore), nufit_params("6.1"))   # Livermore model, normal ordering
r = event_rates(dune, f)                              # events per channel
total_events(r; interaction = :cc)
```

The data tables of SNOwGLoBES (34 MB) are a lazy artifact, downloaded the first time they are needed. The notebook [`notebooks/01_supernova_events.ipynb`](notebooks/01_supernova_events.ipynb) goes through a whole workflow ([open in Colab](https://colab.research.google.com/github/JuliaNeutrino/SNOwGLoBES.jl/blob/main/notebooks/01_supernova_events.ipynb)).

## Differences from SNOwGLoBES

- The efficiencies are applied. In the default run mode, `supernova.pl` builds their file names from a variable it never sets.
- The high-energy channel lists (`argon_he`, `water_he`, ...) use the 200 MeV binning of their header, which SNOwGLoBES skips.
- Backgrounds (`background = "surf"`, `"0800"`, `"4850"` for `ar40kt`) are found. `supernova.pl` looks for a file name that does not exist.
- Flavour conversion keeps θ13, uses the ν̄_μ fluence for antineutrinos and computes the Earth crossing. `oscillate_flux` sets θ13 = 0, uses ν_μ for both and ignores the path length it is given.
- `pinched_fluence` takes the distance, and the exact erg-to-MeV factor.
- CEvNS cross sections: the correct speed of light (the SNOwGLoBES script has 2.988×10¹⁰ cm/s, 0.7% on σ), an accurate integral over recoil energies (5% off at 2 MeV), the low-energy sin²θ_W = 0.2386, the Lewin–Smith Helm radius, and a_k = 0.7 fm for Klein–Nystrand. The SNOwGLoBES tables are cross sections for recoils below 10 keV, the window of its CEvNS smearing matrices; `cevns_cross_sections(target; recoil_max)` gives the same window.

## Tests

`julia --project=. test/runtests.jl`. The event rates of every channel of six detectors (argon, water, scintillator, lead, high-energy argon, xenon CEvNS) and of the argon background agree with C GLoBES to 10⁻⁵ (`test/c_reference.txt`, made by `test/reference/make_reference.jl`). The events also agree with N_targets ∫Φσ dE, the flavour conversion with the formulas of SNOwGLoBES for θ13 = 0 and with the probabilities of GLoBES.jl, and the CEvNS cross sections with the SNOwGLoBES tables for its settings.

SNOwGLoBES is by K. Scholberg and collaborators; this package is under the GPL-2.0, like SNOwGLoBES and GLoBES.
