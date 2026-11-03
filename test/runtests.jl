using SNOwGLoBES, GLoBES, Test

# Events of each channel from C GLoBES (test/reference/make_reference.jl)
const C_REFERENCE =
    [split(l) for l in eachline(joinpath(@__DIR__, "c_reference.txt")) if !startswith(l, "#")]

@testset "data" begin
    @test length(detectors()) == 18
    @test detector(:ar40kt).mass == 40 && detector(:ar40kt).normalization == 0.025
    @test_throws ArgumentError detector(:nope)
    argon = channels(:argon)
    @test length(argon) == 14
    @test argon[1].name == "nue_e" && argon[1].weight == 18 && argon[1].flavour === :e
    @test argon[8].antineutrino && argon[8].name == "nuebar_Ar40"
    @test interaction.(argon[[1, 7, 9]]) == [:es, :cc, :nc]
    @test interaction("ibd") === :ibd && interaction("coh_helm_nue_Xe") === :cevns
    # the binning header of the high-energy lists (SNOwGLoBES leaves them at 100 MeV)
    @test channels(:argon_he)[1].binning == Binning(400, 0.5, 200, 400, 0.5, 200)
    @test channels(:xenon_coh)[1].binning.emax ≈ 0.01
    @test "argon" in channel_lists() && "livermore" in SNOwGLoBES.fluences()
end

@testset "fluences" begin
    f = fluence(:livermore)
    @test length(f.energies) == 501 && f.energies[end] == 100
    @test mean_energy(f).nue ≈ 11.2 atol = 0.01
    @test fluence(:livermore; distance = 5) ≈ 4f
    @test number_of_neutrinos(f + f).nuebar ≈ 2number_of_neutrinos(f).nuebar
    path = write_fluence(tempname(), f)
    @test fluence(path) ≈ f

    # pinched spectra: normalised, with the given mean energy
    E = 0:0.01:300
    @test SNOwGLoBES.trapz(E, pinched_spectrum.(E, 12, 3)) ≈ 1 rtol = 1e-6
    @test SNOwGLoBES.trapz(E, E .* pinched_spectrum.(E, 12, 3)) ≈ 12 rtol = 1e-6
    p = pinched_fluence(;
        mean_energy = (10, 14, 16),
        α = (2.5, 3, 2),
        total_energy = (5e52, 4e52, 3e52),
        energies = E,
    )
    n = number_of_neutrinos(p)
    @test n.nue ≈ 5e52 * SNOwGLoBES.MEV_PER_ERG / 10 / (4π * (10 * SNOwGLoBES.KPC_CM)^2) rtol = 1e-6
    @test mean_energy(p).nuebar ≈ 14 rtol = 1e-6
    @test p.numu == p.nutau == p.numubar == p.nutaubar
    @test number_of_neutrinos(
        pinched_fluence(; mean_energy = (10, 0, 0), α = (2, 2, 2), total_energy = (1e52, 0, 0)),
    ).nuebar == 0

    g = garching_fluences()
    @test length(g.fluences) == length(g.time) == length(g.duration) == 300
    @test issorted(g.time) && all(>(0), g.duration)
end

@testset "flavour conversion" begin
    f = fluence(:livermore)
    p = nufit_params("6.1")
    inverted = nufit_params("6.1"; ordering = :inverted)
    total(x) =
        (n = number_of_neutrinos(x); (n.nue + n.numu + n.nutau, n.nuebar + n.numubar + n.nutaubar))
    for q in (p, inverted), L in (0, 8000)
        @test all(total(oscillate(f, q; earth_path = L)) .≈ total(f))
    end
    # θ13 = 0: the formulas of SNOwGLoBES (supernova_mixing.cc)
    s2 = sin(p.θ12)^2
    no = oscillate(f, OscParams(p; θ13 = 0, θ23 = π / 4))   # and maximal θ23
    @test no.nue ≈ f.numu
    @test no.nuebar ≈ (1 - s2) * f.nuebar + s2 * f.numubar
    @test no.numu ≈ (f.nue + f.numu) / 2
    io = oscillate(f, OscParams(inverted; θ13 = 0, θ23 = π / 4))
    @test io.nue ≈ s2 * f.nue + (1 - s2) * f.numu
    @test io.nuebar ≈ f.numubar
    # the Earth: vacuum for a short path, regeneration of ν̄_e through the mantle
    @test oscillate(f, p; earth_path = 1e-3).spectra ≈ oscillate(f, p).spectra rtol = 1e-6
    U = pmns_matrix(p)
    P = SNOwGLoBES.earth_probabilities(p, U, 30.0, prem_profile(8000, 20), true)
    @test sum(P; dims = 1) ≈ ones(1, 3) && sum(P; dims = 2) ≈ ones(3, 1)
    @test abs(P[1, 1] - abs2(U[1, 1])) > 1e-3
    # the evolution agrees with the flavour probabilities of GLoBES.jl
    flavours = (:e, :μ, :τ)
    for anti in (false, true), E in (10.0, 60.0)
        S = SNOwGLoBES.earth_evolution(p, U, E, ([8000.0], [4.5]), anti)
        @test [abs2(S[a, b]) for a = 1:3, b = 1:3] ≈ [
            probability(
                p,
                flavours[b],
                flavours[a];
                E = E / 1000,
                L = 8000,
                ρ = 4.5,
                antineutrino = anti,
            ) for a = 1:3, b = 1:3
        ] atol = 1e-10
    end
    @test_throws ArgumentError oscillate(f, p; earth_path = 20000)
end

@testset "CEvNS" begin
    # coherent limit: σ = G_F² Q_W² E² / 4π
    iso = SNOwGLoBES.nucleus(:Ar)[3]
    QW = iso.N - (1 - 4 * 0.23857) * iso.Z
    @test cevns_cross_section(1.0, iso; form_factor = :none) ≈
          SNOwGLoBES.G_F^2 * QW^2 / 4π * SNOwGLoBES.ħc_MeV_cm^2 rtol = 1e-3
    @test helm_form_factor(0, 40) == 1 && 0 < helm_form_factor(0.5, 40) < 1
    @test klein_nystrand_form_factor(0, 40) == 1
    @test max_recoil(50, iso.mass) ≈ 2 * 50^2 / (iso.mass + 100)
    @test cevns_differential(10.0, max_recoil(10.0, iso.mass) * 1.01, iso) == 0

    # the settings of SNOwGLoBES reproduce its table, up to its 0.33% error in c
    # and its coarse Riemann sum over the recoil energies
    table = cross_section("coh_helm_nue_Ar")
    old = [
        Isotope(18, 18, 3.350e4, 0.003336),
        Isotope(18, 20, 3.536e4, 0.000629),
        Isotope(18, 22, 3.722e4, 0.996035),
    ]
    ff(q, A) = helm_form_factor(q, A; radius = 1.1425 * A^(1 / 3))
    for E in (10.0, 30.0, 100.0)
        i = argmin(abs.(table.energies .- E))
        σ = sum(
            j.abundance * cevns_cross_section(
                table.energies[i],
                j;
                form_factor = ff,
                sin2θW = 0.231,
                recoil_max = 0.01,
            ) for j in old
        )
        @test table.sigma[i, 1] / σ ≈ (2.988 / 2.998)^2 rtol = 3e-3
    end
    xs = cevns_cross_sections(:Xe; energies = [5.0, 20.0])
    @test xs.sigma[:, 1] == xs.sigma[:, 6] ≈ cevns_cross_section.([5.0, 20.0], :Xe)
end

# Within 1e-5: C GLoBES interpolates cross sections on its fixed 1001-point grid,
# GLoBES.jl on the grid of the file (the tables for water, scintillator, lead)
@testset "event rates against C GLoBES" begin
    groups = unique(r[1:3] for r in C_REFERENCE if r[4] != "bg_surf")
    for (f, chans, det) in groups
        r = event_rates(Symbol(f), Symbol(chans), Symbol(det))
        for ref in filter(x -> x[1:3] == [f, chans, det] && x[4] != "bg_surf", C_REFERENCE)
            c = r[ref[4]]
            @test total_events(c; smeared = false) ≈ parse(Float64, ref[5]) rtol = 1e-5
            @test total_events(c) ≈ parse(Float64, ref[6]) rtol = 1e-5
        end
    end
    ref = only(filter(x -> x[4] == "bg_surf", C_REFERENCE))
    r = event_rates(:livermore, :argon, :ar40kt; background = "surf")
    @test sum(r.background.unsmeared) ≈ parse(Float64, ref[5]) rtol = 1e-6
    @test sum(r.background.smeared) ≈ parse(Float64, ref[6]) rtol = 1e-6
    @test total_events(r) ≈ total_events(event_rates(:livermore, :argon, :ar40kt))
end

@testset "event rates" begin
    f = fluence(:livermore)
    dune = SupernovaDetector(:ar40kt, :argon)
    r = event_rates(dune, f)
    @test length(r) == 14 && keys(r)[7] == "nue_Ar40" && r[:nue_Ar40] === r["nue_Ar40"]
    @test total_events(r) ≈ sum(total_events(c) for c in r)
    @test total_events(r; interaction = :nc) +
          total_events(r; interaction = :cc) +
          total_events(r; interaction = :es) ≈ total_events(r)
    c = r["nue_Ar40"]
    @test length(c.energies) == length(c.unsmeared) == 200 &&
          c.detected_energies[end] ≈ 100 - 99.5 / 400
    # the events are N_targets ∫ Φ σ dE
    x = cross_section("nue_Ar40")
    Es = 0.5:0.01:100
    σ = [SNOwGLoBES.regrid(Fluence(x.energies, x.sigma), [E]).nue[1] for E in Es]
    targets = 40 * 0.025 * 1e9 * 6.02214076e23
    @test total_events(c; smeared = false) ≈
          targets * SNOwGLoBES.trapz(Es, SNOwGLoBES.regrid(f, Es).nue .* σ) rtol = 2e-3
    # linear in the fluence, the mass and the cross sections
    @test total_events(event_rates(dune, 2f)) ≈ 2total_events(r)
    @test total_events(event_rates(f, :argon, :ar40kt; mass = 10)) ≈ total_events(r) / 4
    y = CrossSection(x.energies, 2x.sigma)
    r2 = event_rates(f, :argon, :ar40kt; cross_sections = Dict("nue_Ar40" => y))
    @test total_events(r2["nue_Ar40"]) ≈ 2total_events(c) rtol = 1e-6
    @test total_events(r2["nue_e"]) ≈ total_events(r["nue_e"])
    # a time series: the sum of the steps is the time-integrated fluence
    g = garching_fluences()
    steps = event_rates(dune, g.fluences[1:50])
    @test sum(total_events, steps) ≈ total_events(event_rates(dune, sum(g.fluences[1:50])))
    @test_throws ArgumentError SupernovaDetector(:ar40kt, :water)
    @test_throws ArgumentError SupernovaDetector(:wc100kt30prct, :water; background = "surf")
    @test occursin("total CC", sprint(show, MIME"text/plain"(), r))
end
