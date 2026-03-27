#!/usr/bin/env python3
"""
Coalgebra of dissections.

Formalizes the coproduct structure on rectangle dissections:

  Δ(D) = D ⊗ 1 + 1 ⊗ D + Σ_{fault lines ℓ} D|_L(ℓ) ⊗ D|_R(ℓ)

Irreducible (prime) dissections are PRIMITIVE elements: Δ(D) = D ⊗ 1 + 1 ⊗ D.
Reducible dissections factor along fault lines, and tiling counts multiply:

  |Sol(D)| = ∏_j |Sol(D_j)|

This is the dissection analogue of the Hopf coproduct on symmetric functions,
where Δ(s_λ) = Σ c^λ_{μ,ν} s_μ ⊗ s_ν and the puzzle count equals c^λ_{μ,ν}.

Analogy table:
  Symmetric functions          Dissections
  ─────────────────────        ──────────────────────────
  s_λ (Schur function)         D (dissection)
  Δ(s_λ) (Hopf coproduct)     Σ D|_A ⊗ D|_B (fault-line splitting)
  c^λ_{μ,ν} (LR coefficient)  piece-assignment multiplicity
  p_n (primitive = power sum)  irreducible dissection (no fault lines)
  KTW puzzle witness            tiling of the dissected region
"""

from .geometry import Polygon, area2
from .stomachion import PIECE_VERTS, PIECE_NAMES


# ================================================================
# Fault line detection & restriction
# ================================================================

def fault_lines(pieces, width, height):
    """Find all axis-aligned fault lines of a dissection.

    A fault line y=k (or x=k) exists when NO piece straddles it.
    Returns list of ('h', k) or ('v', k).
    """
    y_ranges = [(min(v[1] for v in p), max(v[1] for v in p)) for p in pieces]
    x_ranges = [(min(v[0] for v in p), max(v[0] for v in p)) for p in pieces]

    faults = []
    for k in range(1, height):
        if not any(y_min < k < y_max for y_min, y_max in y_ranges):
            faults.append(('h', k))
    for k in range(1, width):
        if not any(x_min < k < x_max for x_min, x_max in x_ranges):
            faults.append(('v', k))
    return faults


def restrict(pieces, fault, width, height):
    """Split a dissection along a fault line into two sub-dissections.

    Returns ((pieces_A, w_A, h_A), (pieces_B, w_B, h_B)).
    Coordinates shifted so each sub-region starts at (0, 0).
    """
    direction, k = fault

    if direction == 'h':
        bot, top = [], []
        for p in pieces:
            if max(v[1] for v in p) <= k:
                bot.append(p)
            else:
                top.append(tuple((x, y - k) for x, y in p))
        return (bot, width, k), (top, width, height - k)
    else:
        left, right = [], []
        for p in pieces:
            if max(v[0] for v in p) <= k:
                left.append(p)
            else:
                right.append(tuple((x - k, y) for x, y in p))
        return (left, k, height), (right, width - k, height)


# ================================================================
# Coproduct & factorization
# ================================================================

def coproduct(pieces, width, height):
    """Reduced coproduct Δ̃(D) = Σ D|_A ⊗ D|_B over fault lines.

    Returns list of ((pieces_A, w_A, h_A), (pieces_B, w_B, h_B)).
    Empty list ⟺ D is primitive (irreducible).
    """
    return [restrict(pieces, f, width, height)
            for f in fault_lines(pieces, width, height)]


def factorize(pieces, width, height):
    """Recursively split into irreducible (primitive) components.

    Returns list of (pieces, width, height) — the irreducible factors.
    Tiling count of the original = product of factor tiling counts.
    """
    faults = fault_lines(pieces, width, height)
    if not faults:
        return [(pieces, width, height)]
    a, b = restrict(pieces, faults[0], width, height)
    return factorize(*a) + factorize(*b)


# ================================================================
# Demo
# ================================================================

def demo():
    print("=" * 60)
    print("  Coalgebra of Dissections")
    print("=" * 60)

    # --- The Stomachion: a primitive element ---
    print("\n[1] The Stomachion — primitive element")
    print("-" * 50)

    faults = fault_lines(PIECE_VERTS, 12, 12)
    n_pieces = len(PIECE_VERTS)
    total_a2 = sum(area2(p) for p in PIECE_VERTS)

    print(f"  {n_pieces} pieces, 2×area = {total_a2} (12×12 = {2*144})")
    print(f"  Fault lines: {faults}")
    print(f"  Primitive? {len(faults) == 0}")
    if not faults:
        print(f"  Δ̃(Stomachion) = 0")
        print(f"  → The Stomachion is a generator of the coalgebra.")
        print(f"  → Its 536 tilings are NOT explained by sub-factoring.")

    # --- Construct a reducible dissection by stacking ---
    print(f"\n[2] Constructed reducible dissection")
    print("-" * 50)

    from .explore import random_dissection, count_tilings

    # Use asymmetric halves (12×4 + 12×8) so pieces can't swap sides
    bot = top = None
    for seed in range(300):
        if bot is None:
            bot = random_dissection(12, 4, 3, seed=seed)
        elif top is None:
            top = random_dissection(12, 8, 5, seed=seed)
        if top and bot:
            break

    if not (top and bot):
        print("  Failed to generate sub-dissections")
        return

    # Stack: bottom at y∈[0,4], top at y∈[4,12]
    top_shifted = [tuple((x, y + 4) for x, y in p) for p in top]
    D = bot + top_shifted

    print(f"  {len(D)} pieces = {len(bot)} (bottom 12×4) + {len(top)} (top 12×8)")
    print(f"  2×area = {sum(area2(p) for p in D)}")

    # --- Coproduct ---
    print(f"\n[3] Coproduct Δ̃(D)")
    print("-" * 50)

    faults_D = fault_lines(D, 12, 12)
    print(f"  Fault lines: {faults_D}")

    terms = coproduct(D, 12, 12)
    for i, (a, b) in enumerate(terms):
        pa, wa, ha = a
        pb, wb, hb = b
        print(f"  Term {i+1}: ({len(pa)}pc in {wa}×{ha}) ⊗ ({len(pb)}pc in {wb}×{hb})")

    # --- Full factorization ---
    print(f"\n[4] Irreducible factorization")
    print("-" * 50)

    factors = factorize(D, 12, 12)
    for i, (fp, fw, fh) in enumerate(factors):
        fl = fault_lines(fp, fw, fh)
        print(f"  Factor {i+1}: {len(fp)} pieces in {fw}×{fh}"
              f"  (primitive: {len(fl)==0})")

    # --- Tiling count and assignment multiplicity ---
    print(f"\n[5] Tiling counts & piece-assignment multiplicity")
    print("-" * 50)

    total = count_tilings(D, 12, 12)
    print(f"  |Sol(D)| = {total}  (all placements, all orientations)")

    if terms:
        (pa, wa, ha), (pb, wb, hb) = terms[0]
        ca = count_tilings(pa, wa, ha)
        cb = count_tilings(pb, wb, hb)
        fixed = ca * cb
        print(f"  |Sol(bot {wa}×{ha})| × |Sol(top {wb}×{hb})| = {ca} × {cb} = {fixed}")
        print(f"  Ratio = {total} / {fixed} = {total // fixed}")
        print()
        print(f"  The ratio counts piece-assignment multiplicities:")
        print(f"  rotated pieces may fit in either half, giving extra tilings.")
        print(f"  This is the dissection analogue of LR coefficients c > 1.")
        print(f"  |Sol(D)| = Σ_α |Sol(D_A^α)| × |Sol(D_B^α)|")

    # --- Summary ---
    print(f"\n[6] The analogy")
    print("-" * 50)
    print(f"  In Λ:    Δ(s_λ) = Σ c^λ_{{μ,ν}} s_μ ⊗ s_ν")
    print(f"           puzzle count = c^λ_{{μ,ν}}")
    print(f"           primitive ⟺ s_λ = p_n (power sum)")
    print(f"")
    print(f"  In Diss: Δ(D)   = Σ_{{fault ℓ}} D|_A ⊗ D|_B")
    print(f"           tiling count = |Sol(D)| = ∏|Sol(Dⱼ)|")
    print(f"           primitive ⟺ D irreducible (no fault lines)")
    print(f"           The Stomachion: primitive with 536 tilings.")
    print()


if __name__ == '__main__':
    demo()
