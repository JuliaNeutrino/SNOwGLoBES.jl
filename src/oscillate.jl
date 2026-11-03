# Flavour conversion of supernova neutrinos: adiabatic MSW conversion in the
# star, then optionally the crossing of the Earth

"""
    oscillate(f, params; earth_path=0, layers=20) -> Fluence

The fluence at the detector after adiabatic MSW conversion in the supernova,
for the oscillation parameters `params` (an `OscParams`); the sign of `Δm31`
sets the mass ordering. The neutrinos leave the star as mass eigenstates:

| ordering | ν_e leaves as | ν̄_e leaves as |
|:---------|:--------------|:--------------|
| normal   | ν_3           | ν̄_1           |
| inverted | ν_2           | ν̄_3           |

and the ν_x (the mean of ν_μ and ν_τ) fill the other two. At the Earth the
flux of flavour α is Σᵢ P(νᵢ → ν_α) Fᵢ, with P = |U_αi|² or, for a path of
`earth_path` km through the Earth (PREM, `layers` layers), the probability of
the matter evolution (Earth regeneration).

There are no collective oscillations, shock effects or non-adiabatic
transitions. Unlike the `oscillate_flux` tool of SNOwGLoBES, θ13 is not
neglected and the Earth crossing is computed.

```julia
p = nufit_params("6.1")   # normal ordering
oscillate(fluence(:livermore), p)
oscillate(fluence(:livermore), p; earth_path = 8000)   # detected through 8000 km of Earth
```
"""
function oscillate(f::Fluence, p::OscParams; earth_path::Real = 0, layers::Integer = 20)
    normal = p.Δm31 > 0
    U = pmns_matrix(p)
    earth_path <= 2GLoBES.GLB_EARTH_RADIUS ||
        throw(ArgumentError("earth_path must be at most the Earth's diameter"))
    path = earth_path > 0 ? prem_profile(earth_path, layers) : nothing
    out = similar(f.spectra)
    for (k, E) in enumerate(f.energies)
        F = view(f.spectra, k, :)
        for anti in (false, true)
            Fe, Fx = anti ? (F[4], (F[5] + F[6]) / 2) : (F[1], (F[2] + F[3]) / 2)
            # the mass eigenstate that ν_e (ν̄_e) leaves the star as
            i_e = anti ? (normal ? 1 : 3) : (normal ? 3 : 2)
            Fi = [i == i_e ? Fe : Fx for i = 1:3]
            P = path === nothing || E <= 0 ? abs2.(U) : earth_probabilities(p, U, E, path, anti)
            out[k, (anti ? 4 : 1):(anti ? 6 : 3)] .= P * Fi
        end
    end
    return Fluence(f.energies, out)
end

"""
    earth_probabilities(params, U, E, (lengths, densities), anti) -> Matrix

P(νᵢ → ν_α) for mass eigenstates entering the Earth, row α and column i, at
energy `E` (MeV) through layers of the given lengths (km) and densities (g/cm³).
"""
earth_probabilities(p::OscParams, U, E, path, anti::Bool) =
    abs2.(earth_evolution(p, U, E, path, anti) * (anti ? conj(U) : U))

# The flavour evolution matrix S (ν_β → ν_α is S[α, β]) through the layers
function earth_evolution(p::OscParams, U, E, (lengths, densities), anti::Bool)
    Um = anti ? conj(U) : U
    H0 = Um * Diagonal([0, p.Δm21, p.Δm31]) * Um' / (2E * 1e6)   # eV
    S = Matrix{ComplexF64}(I, 3, 3)
    for (L, ρ) in zip(lengths, densities)
        V = (anti ? -1 : 1) * ρ * GLoBES.GLB_V_FACTOR * GLoBES.GLB_Ne_MANTLE
        λ, Q = eigen(Hermitian(H0 + Diagonal([V, 0, 0])))
        S = Q * Diagonal(cis.(-λ * L / GLoBES.GLB_EV_TO_KM_FACTOR)) * Q' * S
    end
    return S
end
