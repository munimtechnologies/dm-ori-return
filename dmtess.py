# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""DM Ori TESS arm, rules frozen in AMENDMENT-2026-11-25-tess.md (stdlib only, tested).

TESS saw DM Ori in three sectors only: 6 (Dec 2018), 98 (Dec 2025) and 110 (Nov 2026). Sector 6
caught a ~1 mag dip that ASAS-SN g also recorded, so it calibrates TESS flux against ASAS-SN g
flux (relative to quiescence) for this star. Comparison stars put sectors 98 and 110 on the
sector-6 flux scale. The rules below turn that into an ASAS-SN-equivalent depth per sector and
into a map of which ingress dates of a 2013-14-like occultation the snapshots allow.
"""
import math

VERSION = "0.2.0-tess"
RA = 85.9125
DEC = -5.067453
REF_SECTOR = 6
TARGET_SECTORS = (98, 110)

# Occultation geometry from Rodriguez et al. (2016) as frozen in PREREGISTRATION.md
T2_START = 2456320.5          # 2013 ingress start (end of January 2013)
P_ING_LO = 4452               # ingress-to-ingress separations allowed by the 2000 bounds (days)
P_ING_HI = 4692
# Deep phases after ingress start, from the 2013-14 shape: ~80 d ingress, ~30 d partial recovery
# from JD 2456511 (ingress + ~191 d), duration > 420 d.
DEEP_WINDOWS = ((80.0, 185.0), (225.0, 420.0))

DEEP_DMAG = 0.75              # half the published depth (as in the v0.1.0 rules)
BRIGHT_DMAG = 0.30
PAIR_HALF_WINDOW = 0.25       # days: TESS points matched to one ASAS-SN point
MIN_TESS_PER_PAIR = 3
MIN_PAIRS = 6
MIN_CORR = 0.90
MIN_SPAN = 0.40               # range of ASAS-SN g relative flux among the pairs
MIN_COMPS = 4
MAX_SYS = 0.10
SYS_FLOOR = 0.02
DMAG_CAP = 3.0

# ASAS-SN Sky Patrol v2 g points around TESS sector 6, copied from the v0.1.0 output
# (results/v0.1.0/asassn_points.csv: jd, relative flux r, camera). r = flux / quiescent reference.
ASASSN_S6 = (
    (2458462.88481, 0.8520, "br"),
    (2458464.57656, 0.8655, "bF"),
    (2458465.69078, 0.8783, "br"),
    (2458466.59188, 0.8620, "bF"),
    (2458470.88926, 0.8703, "bB"),
    (2458476.56659, 0.5987, "bj"),
    (2458476.90147, 0.7620, "br"),
    (2458480.61614, 0.7338, "bj"),
    (2458483.97896, 0.3866, "bB"),
    (2458486.53128, 0.4047, "bj"),
    (2458486.93895, 0.3647, "bB"),
    (2458487.69468, 0.3530, "bj"),
    (2458488.69119, 0.4040, "bj"),
    (2458489.70249, 0.3314, "bj"),
    (2458490.70148, 0.3482, "bj"),
    (2458491.81794, 0.2962, "br"),
    (2458493.77778, 0.3866, "bB"),
    (2458495.85663, 0.4966, "bB"),
    (2458497.81459, 0.5735, "bB"),
)


def median(xs):
    v = sorted(x for x in xs if x is not None and not (isinstance(x, float) and math.isnan(x)))
    n = len(v)
    if n == 0:
        return None
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


def mad(xs):
    m = median(xs)
    if m is None:
        return None
    return median([abs(x - m) for x in xs])


def pearson(x, y):
    n = len(x)
    if n < 2:
        return None
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)


def ols(x, y):
    """Least squares y = a + b x. Returns (a, b)."""
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    b = sxy / sxx
    return my - b * mx, b


def pairs(asassn, t, f, half=PAIR_HALF_WINDOW, min_n=MIN_TESS_PER_PAIR):
    """For each ASAS-SN point (jd, r, cam) take the median TESS flux within +-half days.
    t must be sorted. Returns list of dicts (jd, r, f, n) for points with >= min_n TESS points."""
    out = []
    for jd, r, cam in asassn:
        sel = [fv for tv, fv in zip(t, f) if abs(tv - jd) <= half]
        if len(sel) >= min_n:
            out.append({"jd": jd, "r": r, "f": median(sel), "n": len(sel), "camera": cam})
    return out


def calibrate(prs):
    """Linear calibration TESS flux = a + b * r_g. Gate: n >= MIN_PAIRS, b > 0,
    Pearson >= MIN_CORR, span of r >= MIN_SPAN."""
    n = len(prs)
    res = {"n_pairs": n, "ok": False}
    if n < MIN_PAIRS:
        res["reason"] = "too few pairs"
        return res
    x = [p["r"] for p in prs]
    y = [p["f"] for p in prs]
    span = max(x) - min(x)
    corr = pearson(x, y)
    res.update({"span": span, "corr": corr, "r_min": min(x), "r_max": max(x)})
    if span < MIN_SPAN:
        res["reason"] = "span too small"
        return res
    a, b = ols(x, y)
    res.update({"a": a, "b": b, "f_quiescent": a + b,
                "nonvarying_share": a / (a + b) if (a + b) != 0 else None})
    if b <= 0:
        res["reason"] = "slope not positive"
        return res
    if corr is None or corr < MIN_CORR:
        res["reason"] = "correlation below %.2f" % MIN_CORR
        return res
    res["ok"] = True
    return res


def equiv_r(F, a, b):
    return (F - a) / b


def equiv_dmag(F, a, b):
    r = equiv_r(F, a, b)
    if r <= 10 ** (-0.4 * DMAG_CAP):
        return DMAG_CAP
    return -2.5 * math.log10(r)


def comp_scale(ratios):
    """Sector scale factor k (median comparison-star ratio to the reference sector) and its
    systematic fraction sys = max(SYS_FLOOR, 1.4826 * MAD / k)."""
    k = median(ratios)
    if k is None or k <= 0:
        return None, None
    m = mad(ratios)
    return k, max(SYS_FLOOR, 1.4826 * m / k)


def classify(dmag_lo, dmag_hi):
    """dmag_lo: depth with the brighter-side scale; dmag_hi: with the fainter-side scale."""
    if dmag_lo >= DEEP_DMAG:
        return "DEEP"
    if dmag_hi <= BRIGHT_DMAG:
        return "BRIGHT"
    if dmag_lo > BRIGHT_DMAG and dmag_hi < DEEP_DMAG:
        return "MID"
    return "UNCERTAIN"


def sector_level(m, k, sys, a, b):
    """m: target median flux in the sector (e-/s); k: comparison scale to sector 6; sys: its
    systematic fraction. Returns dict with depths and the label."""
    F = m / k
    d_mid = equiv_dmag(F, a, b)
    d_lo = equiv_dmag(m / (k * (1.0 - sys)), a, b)
    d_hi = equiv_dmag(m / (k * (1.0 + sys)), a, b)
    return {"F_on_ref_scale": F, "r_equiv": equiv_r(F, a, b), "dmag": d_mid,
            "dmag_lo": d_lo, "dmag_hi": d_hi, "label": classify(d_lo, d_hi),
            "tess_fraction_of_quiescent": F / (a + b)}


def in_deep_window(t, ti):
    d = t - ti
    return any(lo <= d <= hi for lo, hi in DEEP_WINDOWS)


def ingress_map(snapshots):
    """snapshots: list of (t_mid_jd, label). For each allowed ingress start
    ti = T2_START + p (p = P_ING_LO..P_ING_HI, 1-day steps): 'excluded' if a BRIGHT snapshot
    falls in a deep window of ti; 'consistent' if not excluded and a DEEP snapshot falls in one;
    otherwise 'open'."""
    rows = []
    for p in range(P_ING_LO, P_ING_HI + 1):
        ti = T2_START + p
        bright = [t for t, lab in snapshots if lab == "BRIGHT" and in_deep_window(t, ti)]
        deep = [t for t, lab in snapshots if lab == "DEEP" and in_deep_window(t, ti)]
        state = "excluded" if bright else ("consistent" if deep else "open")
        rows.append({"p": p, "ti": ti, "state": state})
    n = len(rows)
    frac = {s: sum(1 for r in rows if r["state"] == s) / n for s in ("excluded", "consistent", "open")}
    return rows, frac


def runs(rows, state):
    out = []
    cur = None
    for r in rows:
        if r["state"] == state:
            if cur is None:
                cur = [r["ti"], r["ti"]]
            else:
                cur[1] = r["ti"]
        elif cur is not None:
            out.append(cur)
            cur = None
    if cur is not None:
        out.append(cur)
    return out


def headline(s98, s110):
    """Pre-registered headline from the two sector labels (None if the sector failed)."""
    a = s98 or "MISSING"
    b = s110 or "MISSING"
    if a == "DEEP" and b in ("DEEP", "MID"):
        return "DIM_DEC2025_NOT_RECOVERED_NOV2026"
    if a == "DEEP" and b == "BRIGHT":
        return "DIM_DEC2025_RECOVERED_NOV2026"
    if a == "BRIGHT" and b == "BRIGHT":
        return "BRIGHT_BOTH"
    if a == "BRIGHT" and b in ("DEEP", "MID"):
        return "BRIGHT_DEC2025_DIM_NOV2026"
    return "MIXED_%s_%s" % (a, b)


def jd_to_date(jd):
    """Calendar date (UTC, proleptic Gregorian) for a Julian Date; Fliegel-Van Flandern."""
    z = int(math.floor(jd + 0.5))
    a = z + 32044
    b = (4 * a + 3) // 146097
    c = a - 146097 * b // 4
    d = (4 * c + 3) // 1461
    e = c - 1461 * d // 4
    m = (5 * e + 2) // 153
    day = e - (153 * m + 2) // 5 + 1
    month = m + 3 - 12 * (m // 10)
    year = 100 * b + d - 4800 + m // 10
    return "%04d-%02d-%02d" % (year, month, day)
