# Connections: Puzzles ↔ Dissections

Two projects, one algebraic structure.

## The objects

| | KTW Puzzles (`~/git/puzzles`) | Stomachion (`~/git/wooden`) |
|---|---|---|
| **Space** | Λ = ring of symmetric functions | Diss = dissections of rectangles |
| **Basis element** | s_λ (Schur function, indexed by partition λ) | D (dissection, indexed by piece shapes) |
| **Product** | s_μ · s_ν = Σ c^λ_{μ,ν} s_λ (induction) | D₁ · D₂ = juxtapose along shared edge |
| **Coproduct** | Δ(s_λ) = Σ c^λ_{μ,ν} s_μ ⊗ s_ν (restriction) | Δ(D) = Σ_{fault lines} D\|_A ⊗ D\|_B |
| **Witness** | KTW puzzle filling the triangle | tiling of the dissected region |
| **Multiplicity** | c^λ_{μ,ν} = puzzle count | piece-assignment count (pieces fit either half via rotation) |
| **Primitive** | p_n (power sum) | irreducible dissection (no fault lines) |
| **The special one** | s_{(3,2,1)}: 25 coproduct terms, c=2 at ((2,1),(2,1)) | Stomachion: primitive, 536 tilings, unexplained by factoring |

## The coproduct

### Symmetric functions

```
Δ(s_λ) = Σ_{μ,ν} c^λ_{μ,ν} s_μ ⊗ s_ν
```

Computed in `transfer_operators.py:lr_via_transfer()` via:
1. T_free chains → skew Kostka numbers K̃_{λ/μ,ρ}
2. Kostka matrix inversion (unitriangular in dominance order) → c^λ_{μ,ν}

Cross-validated against puzzle enumeration in `lr_coefficients.py`.

### Dissections

```
Δ(D) = D ⊗ 1 + 1 ⊗ D + Σ_{fault lines ℓ} D|_A(ℓ) ⊗ D|_B(ℓ)
```

Computed in `solver/dissection_coalgebra.py` via:
1. `fault_lines()` — detect axis-aligned lines no piece crosses
2. `restrict()` — split pieces and shift coordinates
3. `factorize()` — recurse to irreducible components

Tiling count formula:
```
|Sol(D)| = Σ_α |Sol(D_A^α)| × |Sol(D_B^α)|
```
where α ranges over valid piece-to-region assignments (the "LR coefficient" of dissections).

## The failed approach that taught us something

The first attempt at `lr_via_transfer` evaluated s_{λ/μ}(1,...,1) at multiple n and solved a linear system. This fails because hook-content polynomials are linearly dependent:

```
dim_n(3,1) − 3·dim_n(2,2) + dim_n(2,1,1) ≡ 0   for all n
```

The principal specialization (all variables equal) loses information about Schur functions. The Kostka matrix — which tracks SSYT by content, not just by total — has no such degeneracy.

Lesson: the content (distribution across rows) matters, not just the total count.

## The Connes-Kreimer connection

The dissection coalgebra mirrors the **Connes-Kreimer Hopf algebra of Feynman graphs**:

| Feynman graphs | Dissections |
|---|---|
| 1PI (one-particle-irreducible) | prime (no fault lines) |
| subdivergence | fault line |
| renormalization = recursive subtraction | factorize along fault lines |
| forest formula | Sol = ∐ ∏ Sol(Rⱼ) |

The Stomachion is a "1PI graph" — it has no subdivergences (fault lines), so its tiling count is not explained by renormalization (factoring). All 536 tilings are "genuinely interacting."

## Open questions

1. **Antipode**: does the dissection coalgebra have an antipode, making it a full Hopf algebra? The Connes-Kreimer antipode comes from recursive forest subtraction. For dissections, the antipode would "invert" juxtaposition.

2. **Exchange moves as algebra**: exchange moves generate multiplicity *within* irreducible dissections (not across fault lines). They might correspond to a braid group action or crystal basis — a different algebraic structure from the coproduct.

3. **Transfer matrix for dissections**: could row-by-row transfer matrices (as in the KTW puzzle solver) be adapted to count tilings of a dissection? The puzzle's secret trick (△ determinism) might have a dissection analogue.

4. **The Stomachion as algebraic outlier**: in Λ, certain partitions give large LR coefficients (s_{(3,2,1)} has c=2). Is the Stomachion analogous — a primitive element whose "representation dimension" (536) is anomalously large among primes of comparable size?

## Files

| File | What it does |
|---|---|
| `~/git/puzzles/transfer_operators.py` | Fock space operators, `lr_via_transfer()` coproduct |
| `~/git/puzzles/lr_coefficients.py` | KTW puzzle enumeration, `lr_count()`, `lr_puzzles()` |
| `~/git/wooden/solver/dissection_coalgebra.py` | Dissection coproduct, `fault_lines()`, `factorize()` |
| `~/git/wooden/solver/stomachion.py` | Canonical 14-piece Stomachion coordinates |
| `~/git/wooden/solver/primes.py` | `is_irreducible()` — prime dissection detection |
| `~/git/wooden/solver/explore.py` | `random_dissection()`, `count_tilings()` |
