# self_referential_fixed_points

Verification code for the paper

> Mesut Ismail, **Self-Referential Fixed Points and the Physical Horizon of NP-Complete Computation:
> A Structural Analog of Gödel–Löb Incompleteness in Resource-Bounded Computation**, October 2026.

The paper studies a machine whose task size is set by the budget that the computation itself will
leave behind. With initial budget B₀, cost function c(n) and specification φ(B), self-consistency is a
fixed point of the self-consumption operator

    T(B) = B₀ − c(φ(B)).

The paper proves elementary theorems about this operator (Section 4), formulates resource counterparts
of Gödel's and Löb's theorems (Section 5), and evaluates the fixed point for the observable universe
(Sections 6–8). **Every number printed in the paper, and every claim of Sections 4–8 that can be tested
by computation, is recomputed by `verify_selfref.py` as a PASS/FAIL check.** The script is the only
file needed; `runlog_verify_selfref.txt` is its output as run by the author.

## Running

```
pip install mpmath
python3 verify_selfref.py
```

Requirements: Python ≥ 3.9 and `mpmath` (tested with Python 3.13, mpmath 1.3.0). No data files, no
network. Runtime about 30 s. All arithmetic is done at 50 decimal digits; every printed value is
compared with the precision at which it appears in the paper. The last line must read

```
49 checks, 49 PASS, 0 FAIL
```

and the exit status is 1 if any check fails.

## Inputs (Section 6.1 of the paper)

| Quantity | Value | Source |
|---|---|---|
| k_B | 1.380649 × 10⁻²³ J/K | SI (exact) |
| ħ, c, G | CODATA 2018 | |
| H₀ | 67.4 km s⁻¹ Mpc⁻¹ | Planck 2018 |
| R_H | 4.4 × 10²⁶ m | comoving radius |
| t_U | 4.35 × 10¹⁷ s | age of the universe |
| T_CMB | 2.7255 K | |
| E_U | 4 × 10⁶⁹ J | reference energy of the paper |

The budget is B_U = E_U / (k_B T_CMB ln 2); the canonical fixed point is B* = B_U/2 with
n* = log₂ B_U − 1.

## What is checked, section by section

The script prints its checks in groups A–I. Each group corresponds to a part of the paper.

### A — Section 6.1, cosmological parameters
V = (4/3)πR_H³ ≈ 3.6×10⁸⁰ m³; ρ_c = 3H₀²/(8πG) ≈ 8.53×10⁻²⁷ kg/m³; the three energy estimates
ρ_c c²(c t_U)³ ≈ 1.7×10⁶⁹ J, ρ_c c²·(4π/3)(c t_U)³ ≈ 7.1×10⁶⁹ J, ρ_c c²V ≈ 2.7×10⁷¹ J; that the reference
E_U = 4×10⁶⁹ J lies between the first two; and that the third would raise n* by 6.1.

### B — Section 6.2, fixed point at CMB temperature
ε_CMB = 2.608×10⁻²³ J; B_U = 1.534×10⁹²; log₂B_U = 306.234; n* = 305.234 (equation (5), boxed);
{log₂B_U} = 0.234 < log₂(3/2) = 0.585, so the discrete fixed point of Theorem 4.7 exists, with
n*_d = 305 and B*_d = B_U − 2³⁰⁵ = 0.575·B_U.

### C — Table 1 and Section 6.3, thermal regimes and sensitivity
All five rows of Table 1 (ε, B_U and n* at 300 K, 77 K, 4 K, 2.7255 K, 0.01 K); the derivatives
∂n*/∂ln E_U = 1/ln 2 = 1.443 and log₂10 = 3.32 per order of magnitude.

### D — Theorem 4.3, stability
The counterexample stated in Theorem 4.3(b): for δ(B) = −1 + 0.1B/B₀ one has δ′ > 0 but
sup|T′(B)| ≤ 2^(−0.9)(1 + 0.1 ln 2) < 0.58, and the fixed point is attracting (iteration converges).
The sufficient conditions of Theorem 4.3(b) — d ≥ 0, δ′ > 0 ⇒ repelling; d ≤ 0, −1 < s ≤ 0 ⇒ attracting —
on 20 000 random pairs (d, s).

### E — Proposition 4.6, power-law consumption T(B) = B₀ − κB^α
On 3000 random (κ, α, B₀): T′(B*) = −α(B₀ − B*)/B*, and |T′(B*)| < 1 ⇔ B* > αB₀/(1+α).

### F — Theorem 4.7, discrete model φ_d(B) = ⌊log₂B⌋
By brute force on about 6000 budgets (all integers 2…2999 and 3000 random reals up to 2⁶⁰):
a fixed point exists iff B₀ < 3·2^(m−1) with m = ⌊log₂B₀⌋ (equivalently {log₂B₀} < log₂(3/2)), it is
unique, n*_d = m − 1 and B*_d = B₀ − 2^(m−1); otherwise {B₀ − 2^m, B₀ − 2^(m−1)} is a two-cycle with
sizes m − 1 and m. The example B₀ = 7 ({log₂7} = 0.807, two-cycle {3, 5}); exactly half of the integers
in every [2^m, 2^(m+1)) admit a fixed point (m = 1…12); and the referee's table
B₀ = 4, 5, 8, 11, 16, 23 → B*_d = 2, 3, 4, 7, 8, 15 and B₀ = 6, 7, 12, 15, 24 → none.

### F2 — Section 6.3, range of scenarios and robustness of the discrete fixed point
n* ranges from 304.0 to 311.3 over E_U ∈ [1.7×10⁶⁹, 2.7×10⁷¹] J; with 1.5·E_U the fractional part is
0.819 > 0.585 and no discrete fixed point exists.

### G — Section 5, elementary facts used in the proofs
2^⌊log₂B⌋ ≤ B < 2^(⌊log₂B⌋+1) for B ≥ 1 (Theorem 5.4); n_G = ⌊log₂B₀⌋ + 1 gives 2^n_G > B₀, i.e.
the instance of size n_G is infeasible, ¬□_B Π_x, for every B ≤ B₀ (Theorems 5.6, 5.7); and the closure
condition of Theorem 5.9, φ(T(B)) = φ(B) ⇔ B = B₀/2, with the overshoot (B > B*) and undershoot
(B < B*) cases. The provability predicate ⊩_B of Definition 5.1 is not computed: ⊩_B ψ_x requires the
truth value of a SAT instance, which the script does not evaluate.

### H — Sections 6.4–6.5, information bounds, reversible computation, faster algorithms
ℓ_P = 1.616×10⁻³⁵ m; Bekenstein bound 2πk_B R E/(ħc) with R = c t_U, E = E_U: 1.5×10¹²² bits;
holographic bound A/(4ℓ_P²) of the Hubble sphere: of order 10¹²² bits; both exceed B_U by about
thirty orders of magnitude. Margolus–Levitin: N_ML = 2E_U t_U/(πħ) = 1.05×10¹²¹, n*_ML = 401.
PPSZ: 1.30704ⁿ ≤ B_U/2 ⇔ n ≤ 790.

### I — Sections 7 and 8, quantum search and other cost functions
Corollary 7.2: n*_Q = 2 × 305.234 = 610.47. Section 8.1: matched polynomial sizes (B_U/2)^(1/k) for
k = 3, 10, 30, 100, 1000 (4.2×10³⁰, 1.5×10⁹, 1156, 8.3, 1.2) and the exponent k ≈ 37 at which the size
equals 305; the mismatched operator T(B) = B_U − (log₂B)^k has fixed point B_U − (log₂B_U)^k to
relative precision 10⁻⁸⁰ (k = 1, 3, 10). Section 8.4: 2^(n²) → n* = 17.5; 2^(2ⁿ) → n* = 8.25.

## What is not checked

Theorems 4.1, 4.3, 4.6, 4.7 and the results of Section 5 are proved in the paper; the script tests
their statements numerically on many instances, which supports but does not replace the proofs.
The physical assumptions (Landauer accounting, Assumption 2.2; the reference energy E_U) are inputs,
not outputs, of the computation.

## Citation

If you use this code, cite the paper above. Author: Mesut Ismail, Department of Programming and
Computer Technologies, Technical University of Sofia, ORCID 0009-0001-0496-964X.

## License

MIT.
