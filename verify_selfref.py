#!/usr/bin/env python3
"""
verify_selfref.py -- recomputes the numerical values and checks by computation
the elementary claims of

  M. Ismail, "Self-Referential Fixed Points and the Physical Horizon of
  NP-Complete Computation: A Structural Analog of Goedel--Loeb Incompleteness
  in Resource-Bounded Computation" (October 2026).
Repository: https://github.com/mesito/self_referential_fixed_points

Each check prints PASS or FAIL; the exit status is 1 if any check fails.
Requires Python >= 3.9 and mpmath.   Usage:  python3 verify_selfref.py
"""
import random
import sys

from mpmath import mp, mpf, log, pi, floor, sqrt, findroot

mp.dps = 50
RESULTS = []


def check(name, cond, info=""):
    RESULTS.append(bool(cond))
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f"   ({info})" if info else ""))


def close(x, printed, abs_=None, rel=None):
    x, printed = mpf(x), mpf(printed)
    return abs(x - printed) <= (mpf(abs_) if abs_ is not None else mpf(rel) * abs(printed))


def lg(x, b=2):
    return log(x) / log(b)


def sig(x, d):
    return mp.nstr(mpf(x), d, min_fixed=0, max_fixed=0).replace("e+", "e")


# physical constants (SI) and cosmological inputs of Section 6.1
kB, hbar, c, G = mpf("1.380649e-23"), mpf("1.054571817e-34"), mpf("299792458"), mpf("6.67430e-11")
ln2 = log(2)
E_U, T_CMB, t_U, R_H = mpf("4e69"), mpf("2.7255"), mpf("4.35e17"), mpf("4.4e26")
H0 = mpf("67.4") * 1000 / mpf("3.0856775814913673e22")

print("A. Section 6.1 -- cosmological parameters")
V = 4 * pi / 3 * R_H ** 3
rho_c = 3 * H0 ** 2 / (8 * pi * G)
check("V = (4/3) pi R_H^3 = 3.6e80 m^3", sig(V, 2) == "3.6e80", sig(V, 4))
check("rho_c = 3 H0^2/(8 pi G) = 8.53e-27 kg/m^3", sig(rho_c, 3) == "8.53e-27", sig(rho_c, 4))
E1 = rho_c * c ** 2 * (c * t_U) ** 3
E2 = rho_c * c ** 2 * 4 * pi / 3 * (c * t_U) ** 3
E3 = rho_c * c ** 2 * V
check("rho_c c^2 (c t_U)^3 = 1.7e69 J", sig(E1, 2) == "1.7e69", sig(E1, 4))
check("rho_c c^2 (4pi/3)(c t_U)^3 = 7.1e69 J", sig(E2, 2) == "7.1e69", sig(E2, 4))
check("E_U = 4e69 J lies between the two", E1 < E_U < E2)
check("rho_c c^2 V = 2.7e71 J", sig(E3, 2) == "2.7e71", sig(E3, 4))
check("... which raises n* by log2(2.7e71/4e69) = 6.1", close(lg(E3 / E_U), "6.1", abs_="0.05"), mp.nstr(lg(E3 / E_U), 4))

print("B. Section 6.2 -- fixed point at CMB temperature")
eps = kB * T_CMB * ln2
B_U = E_U / eps
check("eps_CMB = k_B T ln 2 = 2.608e-23 J", sig(eps, 4) == "2.608e-23", sig(eps, 6))
check("B_U = 1.534e92", sig(B_U, 4) == "1.534e92", sig(B_U, 6))
check("log2 B_U = 306.234", close(lg(B_U), "306.234", abs_="0.0005"), mp.nstr(lg(B_U), 9))
check("n* = log2 B_U - 1 = 305.234 (abstract: 305.23)", close(lg(B_U) - 1, "305.234", abs_="0.0005"))
th_U = lg(B_U) - floor(lg(B_U))
check("{log2 B_U} = 0.234 < log2(3/2) = 0.585: discrete fixed point exists",
      close(th_U, "0.234", abs_="0.0005") and th_U < lg(mpf(3) / 2) and close(lg(mpf(3) / 2), "0.585", abs_="0.0005"))
check("n*_d = 305, B*_d = B_U - 2^305 = 0.575 B_U",
      int(floor(lg(B_U))) - 1 == 305 and close((B_U - mpf(2) ** 305) / B_U, "0.575", abs_="0.0005"))

print("C. Table 1 and Section 6.3 -- thermal regimes and sensitivity")
for T, e_p, b_p, n_p in [("300", "2.87e-21", "1.39e90", 298), ("77", "7.37e-22", "5.43e90", 300),
                         ("4", "3.83e-23", "1.04e92", 305), ("2.7255", "2.61e-23", "1.53e92", 305),
                         ("0.01", "9.57e-26", "4.18e94", 313)]:
    e = kB * mpf(T) * ln2
    b = E_U / e
    n = lg(b) - 1
    check(f"T = {T} K: eps = {e_p}, B_U = {b_p}, n* ~ {n_p}",
          sig(e, 3) == e_p and sig(b, 3) == b_p and int(floor(n + mpf("0.5"))) == n_p,
          f"{sig(e, 4)}, {sig(b, 4)}, {mp.nstr(n, 6)}")
check("d n*/d ln E_U = 1/ln 2 = 1.443", close(1 / ln2, "1.443", abs_="0.0005"))
check("one order of magnitude moves n* by log2 10 = 3.32", close(lg(10), "3.32", abs_="0.005"))

print("D. Theorem 4.3 -- stability")
B0 = mpf(1)
Td = lambda B: B0 - B * 2 ** (-1 + mpf("0.1") * B / B0)
dTd = lambda B: -2 ** (-1 + mpf("0.1") * B / B0) * (1 + B * ln2 * mpf("0.1") / B0)
mx = max(abs(dTd(B0 * k / 4000)) for k in range(1, 4001))
check("delta = -1 + 0.1 B/B0: delta' > 0 but sup |T'| <= 2^-0.9 (1 + 0.1 ln 2) < 0.58",
      mx <= 2 ** mpf("-0.9") * (1 + mpf("0.1") * ln2) + mpf(10) ** -40 and mx < mpf("0.58"), mp.nstr(mx, 5))
Bs = findroot(lambda B: Td(B) - B, mpf("0.6"))
x = mpf("0.05")
for _ in range(300):
    x = Td(x)
check("... and its fixed point is attracting", abs(x - Bs) < mpf(10) ** -30)
random.seed(1)
ok = True
for _ in range(20000):
    d, s = mpf(random.uniform(-3, 3)), mpf(random.uniform(-0.999, 4))
    t = abs(2 ** d * (1 + s))
    if d >= 0 and s > 0:
        ok &= t > 1
    if d <= 0 and -1 < s <= 0 and not (d == 0 and s == 0):
        ok &= t < 1
check("sufficient conditions of Theorem 4.3(b) (20000 random cases)", ok)

print("E. Proposition 4.6 -- power-law consumption")
ok = True
for _ in range(3000):
    a = mpf(random.uniform(0.05, 4))
    B0 = mpf(10) ** mpf(random.uniform(-2, 6))
    k = mpf(10) ** mpf(random.uniform(-3, 3)) * B0 ** (1 - a)
    lo, hi = mpf(0), B0
    for _ in range(300):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if mid + k * mid ** a > B0 else (mid, hi)
    Bst = (lo + hi) / 2
    der = -a * k * Bst ** (a - 1)
    ok &= abs(der + a * (B0 - Bst) / Bst) < mpf(10) ** -20 * (1 + abs(der))
    if abs(abs(der) - 1) > mpf(10) ** -20:
        ok &= ((abs(der) < 1) == (Bst > a * B0 / (1 + a)))
check("T'(B*) = -alpha(B0-B*)/B*, attracting iff B* > alpha B0/(1+alpha) (3000 cases)", ok)

print("F. Theorem 4.7 -- discrete model")


def phi_d(B):
    n = int(floor(lg(B)))
    while 2 ** (n + 1) <= B:
        n += 1
    while n > 0 and 2 ** n > B:
        n -= 1
    return n


def fixed_points(B0):
    out, n = [], 0
    while 2 ** n <= B0:
        B = B0 - 2 ** n
        if B >= 1 and phi_d(B) == n:
            out.append((B, n))
        n += 1
    return out


ok, okc = True, True
samples = [mpf(k) for k in range(2, 3000)] + [mpf(2) ** mpf(random.uniform(1, 60)) for _ in range(3000)]
for B0 in samples:
    m = phi_d(B0)
    fp = fixed_points(B0)
    if B0 < 3 * mpf(2) ** (m - 1):
        ok &= (len(fp) == 1 and fp[0][1] == m - 1 and fp[0][0] == B0 - mpf(2) ** (m - 1))
    else:
        ok &= (fp == [])
        a, a2 = B0 - mpf(2) ** m, B0 - mpf(2) ** (m - 1)
        okc &= (B0 - mpf(2) ** phi_d(a) == a2 and B0 - mpf(2) ** phi_d(a2) == a
                and phi_d(a) == m - 1 and phi_d(a2) == m)
check(f"fixed point iff B0 < 3*2^(m-1), unique, n*_d = m-1 (brute force, {len(samples)} budgets)", ok)
check("otherwise {B0 - 2^m, B0 - 2^(m-1)} is a two-cycle with sizes m-1, m", okc)
check("B0 = 7: no fixed point, {log2 7} = 0.807, two-cycle {3, 5}",
      fixed_points(mpf(7)) == [] and close(lg(7) - 2, "0.807", abs_="0.0005")
      and 7 - 2 ** phi_d(mpf(3)) == 5 and 7 - 2 ** phi_d(mpf(5)) == 3)
check("exactly half of the integers in [2^m, 2^(m+1)) admit a fixed point (m = 1..12)",
      all(sum(1 for B in range(2 ** m, 2 ** (m + 1)) if fixed_points(mpf(B))) == 2 ** (m - 1) for m in range(1, 13)))

check("critic's table: B0 = 4,5,8,11,16,23 -> B*_d = 2,3,4,7,8,15; 6,7,12,15,24 -> none",
      [fixed_points(mpf(b)) for b in (4, 5, 8, 11, 16, 23)] == [[(2, 1)], [(3, 1)], [(4, 2)], [(7, 2)], [(8, 3)], [(15, 3)]]
      and all(fixed_points(mpf(b)) == [] for b in (6, 7, 12, 15, 24)))
print("F2. Section 6.3 -- range of scenarios and robustness of the discrete fixed point")
rng = [lg(E / eps) - 1 for E in (E1, E3)]
check("n* ranges from 304.0 to 311.3 over E_U in [1.7e69, 2.7e71]",
      close(rng[0], "304.0", abs_="0.05") and close(rng[1], "311.3", abs_="0.05"), f"{mp.nstr(rng[0],5)}, {mp.nstr(rng[1],5)}")
L15 = lg(mpf("1.5") * B_U)
check("1.5 E_U: {log2 B_U'} = 0.819 > 0.585, no discrete fixed point",
      close(L15 - floor(L15), "0.819", abs_="0.0005") and L15 - floor(L15) > lg(mpf(3) / 2))
print("F3. Theorem 4.8 -- periodic orbits of the discrete operator (integer budgets)")


def Td_int(B, B0):
    return B0 - (1 << (B.bit_length() - 1))


cnt = {"fp_only": 0, "both": 0, "cyc_only": 0}
bad = 0
for B0 in range(2, 100001):
    m = B0.bit_length() - 1
    r = B0 - (1 << m)
    cands = {B0 - (1 << k) for k in range(m + 1)} | {B0}   # only image points can be periodic
    cands = {x for x in cands if 1 <= x <= B0}
    fps = sorted(x for x in cands if Td_int(x, B0) == x)
    cyc = set()
    for a in cands:
        b = Td_int(a, B0)
        if 1 <= b <= B0 and b != a and Td_int(b, B0) == a:
            cyc.add(frozenset((a, b)))
    pfp = [(1 << (m - 1)) + r] if r < (1 << (m - 1)) else []
    pcyc = {frozenset((r, B0 - (1 << (r.bit_length() - 1))))} if r > 0 else set()
    if fps != pfp or cyc != pcyc:
        bad += 1
    if r == 0:
        cnt["fp_only"] += 1
    elif pfp:
        cnt["both"] += 1
    else:
        cnt["cyc_only"] += 1
check("fixed points and two-cycles as stated in Theorem 4.8, all integer B0 in [2, 100000]", bad == 0)
check("regimes: 16 power-of-two budgets, 65519 with both, 34464 with a two-cycle only",
      cnt == {"fp_only": 16, "both": 65519, "cyc_only": 34464}, str(cnt))
check("B0 = 5: fixed point 3 and two-cycle {1, 4}", Td_int(3, 5) == 3 and Td_int(1, 5) == 4 and Td_int(4, 5) == 1)
bb = 0
for B0 in range(2, 4001):
    m = B0.bit_length() - 1
    r = B0 - (1 << m)
    if r == 0:
        fp = 1 << (m - 1)
        for B in range(1, B0):
            x, st = B, 0
            while x != fp and st < 3:
                x, st = Td_int(x, B0), st + 1
            bb += (x != fp or st > 2)
        continue
    i = r.bit_length() - 1
    fp = (1 << (m - 1)) + r if r < (1 << (m - 1)) else None
    cyc = {r, B0 - (1 << i)}
    for B in range(1, B0 + 1):
        k = B.bit_length() - 1
        x, st = B, 0
        if i < k <= m - 1:
            if fp is None:
                bb += 1
                continue
            while x != fp and st < 3:
                x, st = Td_int(x, B0), st + 1
            bb += (x != fp or st > 2)
        else:
            while x not in cyc and st < 3:
                x, st = Td_int(x, B0), st + 1
            bb += (x not in cyc or st > 2)
check("Theorem 4.8(c): basins and 'at most two steps', all B for all B0 <= 4000", bb == 0)
B0i = int(floor(B_U))
mi = B0i.bit_length() - 1
ri = B0i - (1 << mi)
ii = ri.bit_length() - 1
check("B_U: m = 306, r/2^m = 0.176, i = 303: two-cycle between sizes 303 and 306; fixed-point basin sizes 304-305",
      mi == 306 and close(mpf(ri) / mpf(2) ** mi, "0.176", abs_="0.0005") and ii == 303 and ri < (1 << (mi - 1)),
      f"{mi}, {mp.nstr(mpf(ri)/mpf(2)**mi, 4)}, {ii}")
print("G. Section 5 -- facts used in the proofs")
ok = True
for _ in range(5000):
    B = mpf(2) ** mpf(random.uniform(0, 80))
    n = phi_d(B)
    ok &= 2 ** n <= B < 2 ** (n + 1)
check("2^floor(log2 B) <= B < 2^(floor(log2 B)+1) for B >= 1", ok)
check("n_G = floor(log2 B0) + 1 gives 2^n_G > B0 (B0 = 2, 7, 1000, B_U)",
      all(2 ** (phi_d(B) + 1) > B for B in (mpf(2), mpf(7), mpf(1000), B_U)))
ok = True
for _ in range(3000):
    B0 = mpf(10) ** mpf(random.uniform(0.5, 60))
    B = B0 * mpf(random.uniform(0.001, 0.999))
    closes = abs(lg(B0 - B) - lg(B)) < mpf(10) ** -30
    ok &= closes == (abs(B - B0 / 2) < mpf(10) ** -30 * B0)
    if B > B0 / 2 * (1 + mpf(10) ** -20):
        ok &= B > B0 - B
    if B < B0 / 2 * (1 - mpf(10) ** -20):
        ok &= B < B0 - B
check("Theorem 5.9: closure phi(T(B)) = phi(B) iff B = B0/2; overshoot/undershoot cases", ok)

print("H. Sections 6.4-6.5 -- information bounds, reversible computation, PPSZ")
S_bek = 2 * pi * (c * t_U) * E_U / (hbar * c) / ln2
lP = sqrt(hbar * G / c ** 3)
S_holo = 4 * pi * (c / H0) ** 2 / (4 * lP ** 2) / ln2
check("l_P = 1.616e-35 m", sig(lP, 4) == "1.616e-35")
check("Bekenstein bound, R = c t_U, E = E_U: 1.5e122 bits", sig(S_bek, 2) == "1.5e122", sig(S_bek, 4))
check("holographic bound of the Hubble sphere: of order 1e122 bits", 122 <= log(S_holo, 10) < 123, sig(S_holo, 3))
check("both exceed B_U by about thirty orders of magnitude",
      29 < log(S_bek / B_U, 10) < 31 and 29 < log(S_holo / B_U, 10) < 31)
N_ML = 2 * E_U * t_U / (pi * hbar)
check("N_ML = 2 E_U t_U/(pi hbar) = 1.05e121", sig(N_ML, 3) == "1.05e121", sig(N_ML, 4))
check("n*_ML = log2 N_ML - 1 = 401", int(floor(lg(N_ML) - 1 + mpf("0.5"))) == 401, mp.nstr(lg(N_ML) - 1, 6))
n_ppsz = lg(B_U / 2) / lg(mpf("1.30704"))
check("1.30704^n <= B_U/2  <=>  n <= 790", int(floor(n_ppsz)) == 790, mp.nstr(n_ppsz, 6))
nm3, nm6 = lg(B_U / 2) - lg(10 ** 3), lg(B_U / 2) - lg(10 ** 6)
check("bit-cost model c = m 2^n: n* = log2(B_U/2) - log2 m = 295 (m = 1e3), 285 (m = 1e6)",
      int(floor(nm3 + mpf("0.5"))) == 295 and int(floor(nm6 + mpf("0.5"))) == 285, f"{mp.nstr(nm3,5)}, {mp.nstr(nm6,5)}")

print("I. Sections 7-8 -- quantum and other cost functions")
check("n*_Q = 2 x 305.234 = 610.47", close(2 * (lg(B_U) - 1), "610.47", abs_="0.005"))
check("matched n^3: (B_U/2)^(1/3) = 4e30 (4.2e30)", sig((B_U / 2) ** (mpf(1) / 3), 2) == "4.2e30")
check("matched n^10: (B_U/2)^(1/10) = 1.5e9", sig((B_U / 2) ** (mpf(1) / 10), 2) == "1.5e9")
v30, v100, v1000 = [(B_U / 2) ** (mpf(1) / k) for k in (30, 100, 1000)]
check("matched n^k, k = 30, 100, 1000: 1156, 8.3, 1.2",
      int(floor(v30 + mpf("0.5"))) == 1156 and close(v100, "8.3", abs_="0.05") and close(v1000, "1.2", abs_="0.05"),
      f"{mp.nstr(v30,5)}, {mp.nstr(v100,3)}, {mp.nstr(v1000,3)}")
k305 = log(B_U / 2) / log(305)
check("exponent with n* = 305: k = 37", close(k305, "37", abs_="0.5"), mp.nstr(k305, 4))
ok = True
for kk in (1, 3, 10):
    lo, hi = mpf(1), B_U
    for _ in range(400):
        mid = (lo + hi) / 2
        lo, hi = (lo, mid) if mid + lg(mid) ** kk > B_U else (mid, hi)
    ok &= abs((lo + hi) / 2 - (B_U - lg(B_U) ** kk)) / B_U < mpf(10) ** -80
check("mismatched phi = log2 B, c = n^k: B* = B_U - (log2 B_U)^k (k = 1, 3, 10)", ok)
check("EXPTIME-type 2^(n^2): n* = sqrt(log2(B_U/2)) = 17.5", close(sqrt(lg(B_U / 2)), "17.5", abs_="0.05"))
check("doubly exponential 2^(2^n): n* = log2 log2(B_U/2) = 8 (8.25)", close(lg(lg(B_U / 2)), "8.25", abs_="0.005"))

nfail = RESULTS.count(False)
print(f"{len(RESULTS)} checks, {len(RESULTS) - nfail} PASS, {nfail} FAIL")
sys.exit(1 if nfail else 0)
