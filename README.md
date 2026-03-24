# Stomachion

An interactive web puzzle and computational investigation of Archimedes' Stomachion — the oldest known combinatorial puzzle.

**[Play it here](https://raggedr.github.io/stomachion)**

The Stomachion dissects a square into 14 pieces (11 triangles, 2 quadrilaterals, and a pentagon) on an integer lattice. It has exactly **536 distinct tilings** (or 2,144 if you distinguish congruent pieces) — a remarkably large number for a prime dissection.

## Interactive puzzle

`index.html` is a self-contained drag-and-drop puzzle. Grab pieces and place them on the board to tile the square.

- **Drag** to move pieces
- **R** = rotate, **F** = flip (select a piece first)
- Scroll wheel or right-click also rotate
- **Scatter** to shuffle, **Solve** to see the canonical tiling

## Solver

A Python DLX (Dancing Links) exact-cover solver that enumerates all tilings of the 12x12 square.

```
cd solver
python stomachion.py    # verify the 17,152 raw solutions
python experiment.py    # random dissection experiments
python primes.py        # prime dissection analysis
```

Key files:
- `dlx.py` — Algorithm X with dancing links
- `geometry.py` — lattice polygon operations
- `pieces.py` — Stomachion piece definitions and transformations
- `solver.py` — exact-cover formulation for lattice dissections
- `primes.py` — prime (irreducible) dissection detection and enumeration
- `experiment.py` — random dissection generation and tiling count distributions
- `explore.py` — search for high-count dissections
- `EXPERIMENT.md` — full experimental results and analysis

## Paper

`paper/primes.tex` — *Prime Dissections and the Optimality of Archimedes' Stomachion*

Introduces the notion of a **prime dissection** (one with no fault lines) and shows that any lattice dissection factors uniquely into prime sub-dissections — analogous to the fundamental theorem of arithmetic. Presents computational evidence that the Stomachion is maximal among prime dissections of comparable size, and that the distribution of tiling counts for random dissections follows a Poisson-like law in the number of independent exchange moves.
