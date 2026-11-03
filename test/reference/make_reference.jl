# Writes test/c_reference.txt: the events of each channel of a few
# configurations, computed by C GLoBES from the experiments SNOwGLoBES.jl writes.
#
#   julia --project=. test/reference/make_reference.jl /path/to/rates_reference
#
# with rates_reference.c built against C GLoBES 3.0.11 patched with
# test/reference/c-fixes.patch of GLoBES.jl. C GLoBES cannot open long paths,
# so the experiments are rewritten with paths relative to a directory holding
# links to the data. It also loses definitions in files with many large
# smearing matrices, so each channel gets its own experiment (the rates of a
# channel do not depend on the others).
using SNOwGLoBES

const CONFIGS = [
    (:livermore, :argon, :ar40kt),
    (:gvkm, :water, :wc100kt30prct),
    (:livermore, :scint, :scint20kt),
    (:livermore, :lead, :halo1),
    (:livermore, :argon_he, :ar40kt_he),
    (:livermore, :xenon_coh, :lz),
]

driver = abspath(ARGS[1])
work = mktempdir()
for d in ("smear", "xscns", "effic")
    symlink(joinpath(data_directory(), d), joinpath(work, d))
end

# Events and detected events of the GLoBES channels of a detector, from C
function c_rates(s::SupernovaDetector)
    text = read(joinpath(s.dir, "experiment1.glb"), String)
    text =
        replace(text, data_directory() * "/" => "", joinpath(s.dir, "fluence.dat") => "fluence.dat")
    write(joinpath(work, "experiment.glb"), text)
    lines = cd(() -> readlines(`$driver experiment.glb`), work)
    return [parse.(Float64, split(l)[3:4]) for l in lines]
end

open(joinpath(@__DIR__, "..", "c_reference.txt"), "w") do out
    println(
        out,
        "# fluence channels detector channel events detected (C GLoBES 3.0.11 + c-fixes.patch)",
    )
    for (f, chans, det) in CONFIGS
        write_fluence(joinpath(work, "fluence.dat"), fluence(f))
        for c in channels(chans)
            events, detected = only(c_rates(SupernovaDetector(det, [c])))
            println(out, "$f $chans $det $(c.name) $(c.weight * events) $(c.weight * detected)")
        end
    end
    # the background, with one channel of ar40kt
    write_fluence(joinpath(work, "fluence.dat"), fluence(:livermore))
    s = SupernovaDetector(:ar40kt, channels(:argon)[1:1]; background = "surf")
    events, detected = c_rates(s)[2]
    println(out, "livermore argon ar40kt bg_surf $events $detected")
end
