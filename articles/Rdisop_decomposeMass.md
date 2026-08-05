# Rdisop mass decomposition (\`decomposeMass\`)

Algorithms: Böcker et al. (DECOMP / Money Changing Problem). See
Bioconductor vignette *Mass decomposition with the Rdisop package*.

------------------------------------------------------------------------

## 1. High-level workflow (formula calculator)

``` text
m/z + polarity/adduct + ppm + element set
        ↓
1. Convert ion m/z → neutral (or charged) mass consistently
        ↓
2. Enumerate formulas (Rdisop::decomposeMass / decomposeIsotopes)
        ↓
3. Filter chemically (parity, H/C, heteroatom rules, DBE / Senior)
        ↓
4. Rank (ppm, RDBe, isotope fit if MS1 pattern available)
        ↓
5. Return tidy table of candidates
```

| Input           | Role                             | Example            |
|-----------------|----------------------------------|--------------------|
| `mz` / `mass`   | Measured ion m/z or neutral mass | `203.0526`         |
| `ppm` / `mzabs` | Mass tolerance                   | `5` ppm            |
| `z`             | Charge (sign = polarity)         | `+1` / `-1`        |
| `adduct`        | Optional ion type                | `[M+H]+`, `[M-H]-` |
| `elements`      | Allowed CHNOPS…                  | `C,H,N,O`          |

Decide explicitly whether the input is **ion m/z** or **neutral mass**.
Do not mix “neutral mass + `z = 1`” unless that is intentional.

### Important: `z` does **not** convert m/z → formula mass

Rdisop’s `z` is **not** a full ESI charge-state converter.

| Call | What `z` actually does |
|----|----|
| `decomposeMass(..., z = k)` | Does **not** multiply the input by `\|k\|`. Same `mass` with `z = 1` or `z = 2` yields the **same** formulas (mass window still centered on the value you passed). `z` mainly affects charge-related bookkeeping (e.g. parity / validity / how masses are reported). |
| `getMolecule(formula, z = k)` | Adjusts by $`k \times m_e`$, then for $`\|k\| > 1`$**divides** so the returned value looks like an **m/z**. Example: `C4H7N6O4` → ~203.05 (`z = 0/1`) vs ~101.53 (`z = 2`), i.e. $`(M - 2 m_e)/2`$. |

Empirical check (`mass = 203.0526`, CHNO, 5 ppm):

``` r

getFormula(decomposeMass(203.0526, ppm = 5, z = 1, elements = elements))
# C4H7N6O4, C3H11N2O8, ...
getFormula(decomposeMass(203.0526, ppm = 5, z = 2, elements = elements))
# same list — still ~203 Da formulas, not ~406 Da
getMass(getMolecule("C4H7N6O4", z = 2))
# ~101.5259  ← display m/z only, not what decomposeMass searches
```

**Correct pre-conversion** before calling `decomposeMass` (bare charge /
electron correction; adducts still separate):

``` math
M_{\text{search}} = m/z \times |z| \;\pm\; z \cdot m_e
```

For `m/z = 203.0526`, `z = 2` you must decompose **~406.1 Da**, not 203
and not 101.5.

- **101.5** is what `getMolecule(..., z = 2)` *reports* for a ~203 Da
  formula (output m/z).
- **406** is the mass scale of a 2+ ion whose *observed* m/z is 203.

[`MassTools::calcMF`](https://rdrr.io/pkg/MassTools/man/calcMF.html)
only does `mz + z * m_e` then calls `decomposeMass` — it also does
**not** multiply by `|z|`, so it has the same gap for $`|z| > 1`$.

**Rule for a formula calculator:** convert m/z → search mass yourself
(`mz * |z| ± z*me` and/or adduct correction via MSCC), then call
`decomposeMass(..., z = 0)` (or `z = ±1` only if you intentionally want
ion-parity semantics on an already-converted mass).

In the MSdev / MSCC split:

| Layer | Responsibility |
|----|----|
| **Rdisop** | Enumerate candidates from mass |
| **MSCC** | Adduct mass, formula format, isotope pattern |
| **MSdev** | Call calculator on features (`mzmed`), attach `candidate.formula` |

------------------------------------------------------------------------

## 2. What `decomposeMass` does

`decomposeMass` is a thin wrapper:

``` text
decomposeMass(mass, ppm, ...)
  → decomposeIsotopes(c(mass), intensity = 1, ...)
  → C++ Money-Changing / mass-decomposition solver
```

Core problem (Money Changing Problem / MCP):

> Find non-negative integers $`(n_C, n_H, n_N, \ldots)`$ such that  
> $`\sum n_e \cdot m_e`$ lies within the mass tolerance window.

Typical R call:

``` r

library(Rdisop)

elements <- initializeElements(c("C", "H", "N", "O"))
results <- decomposeMass(
  mass = 203.0526,
  ppm = 5,
  z = 1,
  elements = elements
)
getFormula(results)
```

Useful accessors:
[`getFormula()`](https://rdrr.io/pkg/Rdisop/man/getMolecule.html),
[`getMass()`](https://rdrr.io/pkg/Rdisop/man/getMolecule.html),
[`getScore()`](https://rdrr.io/pkg/Rdisop/man/getMolecule.html),
[`getValid()`](https://rdrr.io/pkg/Rdisop/man/getMolecule.html),
[`getIsotope()`](https://rdrr.io/pkg/Rdisop/man/getMolecule.html).

With measured isotope peaks, prefer
`decomposeIsotopes(masses, intensities, ...)` — same enumeration core,
better ranking.

------------------------------------------------------------------------

## 3. Internal steps

``` text
target mass + ppm (+ mzabs)
        │
        ▼
build element alphabet (monoisotopic masses), sorted by mass
        │
        ▼
MCP / mass decomposition
  → all integer vectors n with |Σ n·m − mass| ≤ tol
        │
        ▼
optional chemical validity (DBE, nitrogen rule, …)
        │
        ▼
molecule objects (formula, exactmass, score, isotopes, …)
```

------------------------------------------------------------------------

## 4. Brute force vs Rdisop

**Naive enumeration** (nested loops over atom counts) scales roughly as
$`\prod_e (n_e^{\max}+1)`$ — exponential in the number of element types
when bounds grow with mass.

**Rdisop / DECOMP** does **not** do that. It uses efficient MCP
algorithms (dynamic programming + smart backtracking; Böcker & Lipták):

- Finds **all** compositions whose mass falls in the tolerance window
- Runtime is typically closer to **proportional to the number of
  solutions**, not to a full nested-loop lattice walk
- Far better memory profile than classical full DP tables for this
  problem

So: **yes**, it enumerates all mass-matching compositions for the
allowed elements; **no**, it is not a dumb “try every combination”
search.

------------------------------------------------------------------------

## 5. Does adding elements explode the cost?

**Yes in practice** — mainly because the **solution space** grows, not
because the search is naive.

| Factor | Effect |
|----|----|
| More element species | Larger alphabet → more ways to hit the same mass → more candidates |
| Higher mass | More atoms possible → more solutions |
| Wider `ppm` / `mzabs` | Wider window → more solutions |
| Algorithm | MCP-style; pay per solution + DP overhead, not $`k^{n}`$ nested loops |

Rdisop’s vignette notes that the result list grows with mass, allowed
ppm, and the allowed elements list.

------------------------------------------------------------------------

## 6. Practical controls

1.  Keep the alphabet small (CHNO or CHNOPS first).
2.  Tighten `ppm` / `mzabs`.
3.  Use `minElements` / `maxElements` to cap atom counts.
4.  Prefer `decomposeIsotopes` when M+1 / M+2 intensities exist.
5.  Post-filter (H/C ranges, RDBE, nitrogen rule) and keep top-$`N`$ by
    ppm / score.
6.  Expect a sharp jump in candidates when adding Cl, Br, metals, or
    many heteroatoms.

------------------------------------------------------------------------

## 7. Suggested return table

After decomposition, tidy to something like:

| Column          | Meaning                                       |
|-----------------|-----------------------------------------------|
| `formula`       | Sum formula string                            |
| `exactmass`     | Calculated exact mass                         |
| `ppm`           | $`(exactmass - target) / target \times 10^6`$ |
| `score`         | Rdisop / isotope score                        |
| `DBE` / `valid` | Chemical plausibility hints                   |

------------------------------------------------------------------------

## 8. Bottom line

- `decomposeMass` returns essentially **all** compositions for the
  allowed elements that fit the mass window.
- The engine is an **MCP algorithm**, not brute-force nested loops.
- Adding element types still hurts a lot: the **number of valid
  formulas** grows combinatorially, so CPU and memory grow with output
  size.
- **`z` is not a substitute for m/z → mass conversion.** For
  $`|z| > 1`$, multiply (and correct electrons/adducts) yourself before
  decomposing.
- A usable formula calculator = charge/adduct mass handling + Rdisop
  enumeration + chemical filters + ranking (isotope / DB).
