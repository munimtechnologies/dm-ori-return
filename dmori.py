# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""DM Ori return test: frozen rules (PREREGISTRATION.md, 2026-11-24). Stdlib only, tested.

Question: did the young star DM Ori start the deep (~1.5 mag), >1-year dimming that Rodriguez et al.
(2016, ApJ 831, 74; arXiv:1607.07455) predicted for July 2026, 12.5 years after the 2013-2014 event?
"""
import csv
import io
import math
from datetime import date, timedelta

VERSION = "0.1.0"

RA, DEC = 85.912500, -5.067453  # Tycho-2 J2000 05:43:39.0 -05:04:02.83 (Rodriguez et al. 2016)

_J2000_0 = 2451544.5  # JD of 2000-01-01 00:00


def jd(y, m, d):
    return _J2000_0 + (date(y, m, d) - date(2000, 1, 1)).days


def date_str(t):
    return (date(2000, 1, 1) + timedelta(days=int(math.floor(t - _J2000_0)))).isoformat()


def decimal_year(t):
    y = int(date_str(t)[:4])
    a, b = jd(y, 1, 1), jd(y + 1, 1, 1)
    return y + (t - a) / (b - a)


# --- the published events (Rodriguez et al. 2016) ---------------------------------------------
T2_START = jd(2013, 1, 28)   # ingress began "near the end of UT 2013 January"
T2_END = jd(2014, 9, 1)      # egress "during the 2014 seasonal observing gap" ("2014 September")
T2_CENTER = 0.5 * (T2_START + T2_END)
PRED = jd(2026, 7, 1)        # "the next predicted dimming in July of 2026"
DAYS_PER_YEAR = 365.25
# 2000-2002 event: last bright NSVS point 2000-03-25, first dim ASAS point 2000-11-20.
E1_LAST_BRIGHT = jd(2000, 3, 25)
E1_FIRST_DIM = jd(2000, 11, 20)
P_ING_LO = T2_START - E1_FIRST_DIM     # 4452 d = 12.19 yr (ingress-to-ingress, shortest)
P_ING_HI = T2_START - E1_LAST_BRIGHT   # 4692 d = 12.85 yr (longest)
INGRESS_DAYS = 80.0                    # "an ~80 day period" (2013 ingress)
MIN_DURATION = 420.0                   # 2013-14 event lasted > 420 d
# Every ingress-to-ingress period in [P_ING_LO, P_ING_HI] with the 2013 shape puts DM Ori deep in
# the dimming during this window (2026-02-21 .. 2026-06-01):
CORE_LO = T2_START + P_ING_HI + INGRESS_DAYS
CORE_HI = T2_START + P_ING_LO + MIN_DURATION

# --- frozen thresholds -------------------------------------------------------------------------
BIN_DAYS = 30
BIN_ORIGIN = 2450000.5
MIN_PER_BIN = 3
MIN_SNR = 5.0
DIM_R = 10 ** (-0.4 * 0.75)     # >= 0.75 mag below quiescence (half the published ~1.5 mag depth)
BRIGHT_R = 10 ** (-0.4 * 0.30)  # within 0.30 mag of quiescence
MAX_GAP = 150.0                  # days of no data allowed inside one episode (seasonal gap ~110 d)
MIN_DIM_BINS = 2
REF_WINDOWS = {"V": (jd(2014, 12, 1), jd(2019, 1, 1)), "g": (jd(2017, 12, 1), jd(2025, 1, 1))}
MIN_REF = 20
HALF_DURATION = 210.0            # half of the 2013-14 lower limit (>420 d)
P_GRID = [round(11.0 + 0.01 * k, 2) for k in range(301)]  # 11.00 .. 14.00 yr
P_STATED = (11.5, 13.5)          # 12.5 +- 1 yr

J_2024_07 = jd(2024, 7, 1)
J_2025_07 = jd(2025, 7, 1)
J_2025_10 = jd(2025, 10, 1)
J_2026_01 = jd(2026, 1, 1)
J_2026_02 = jd(2026, 2, 1)
J_2026_05 = jd(2026, 5, 1)
J_2026_08 = jd(2026, 8, 1)
J_2027_01 = jd(2027, 1, 1)


def median(xs):
    s = sorted(xs)
    n = len(s)
    if n == 0:
        return float("nan")
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def robust_scatter(xs):
    if len(xs) < 2:
        return float("nan")
    m = median(xs)
    return 1.4826 * median([abs(x - m) for x in xs])


def dmag(r):
    return -2.5 * math.log10(r) if r > 0 else float("inf")


# --- normalisation ----------------------------------------------------------------------------
def clean(by_band):
    """{'V'|'g': [(t, flux, err, cam)]} -> same, keeping flux > 0 and flux/err >= MIN_SNR."""
    out = {}
    for band, rows in by_band.items():
        out[band] = [(t, f, e, c) for t, f, e, c in rows if f > 0 and e > 0 and f / e >= MIN_SNR]
    return out


def references(by_band):
    """Quiescent reference flux per (band, camera) from REF_WINDOWS: the camera's own median if it
    has >= MIN_REF points in the window, else the band median over the window (all cameras).
    Bands with < MIN_REF points in the window get no reference (unused)."""
    refs = {}
    meta = {}
    for band, rows in by_band.items():
        if band not in REF_WINDOWS:
            continue
        a, b = REF_WINDOWS[band]
        inw = [r for r in rows if a <= r[0] < b]
        if len(inw) < MIN_REF:
            meta[band] = {"n_ref": len(inw), "used": False}
            continue
        band_ref = median([r[1] for r in inw])
        meta[band] = {"n_ref": len(inw), "used": True, "band_ref_flux": band_ref, "cameras": {}}
        for cam in sorted(set(r[3] for r in rows)):
            own = [r[1] for r in inw if r[3] == cam]
            if len(own) >= MIN_REF:
                refs[(band, cam)] = median(own)
                meta[band]["cameras"][cam] = {"ref": "own", "n": len(own), "flux": refs[(band, cam)]}
            else:
                refs[(band, cam)] = band_ref
                meta[band]["cameras"][cam] = {"ref": "band", "n_own": len(own), "flux": band_ref}
    return refs, meta


def normalise(by_band, refs):
    out = []
    for band, rows in by_band.items():
        for t, f, e, cam in rows:
            ref = refs.get((band, cam))
            if ref:
                out.append((t, f / ref, band, cam))
    out.sort()
    return out


# --- binning and states -----------------------------------------------------------------------
def state_of(r):
    if r <= DIM_R:
        return "DIM"
    if r >= BRIGHT_R:
        return "BRIGHT"
    return "MID"


def bin_points(norm, width=BIN_DAYS, origin=BIN_ORIGIN, min_n=MIN_PER_BIN):
    groups = {}
    for t, r, band, cam in norm:
        k = int(math.floor((t - origin) / width))
        groups.setdefault(k, []).append((t, r, band))
    bins = []
    for k in sorted(groups):
        pts = groups[k]
        if len(pts) < min_n:
            continue
        mr = median([p[1] for p in pts])
        lo = origin + k * width
        bins.append({"k": k, "t_lo": lo, "t_hi": lo + width, "t_mid": lo + width / 2.0, "n": len(pts),
                     "med_r": mr, "dmag": dmag(mr), "state": state_of(mr),
                     "bands": "+".join(sorted(set(p[2] for p in pts)))})
    return bins


def episodes(bins, max_gap=MAX_GAP, min_dim=MIN_DIM_BINS):
    """Maximal runs of non-BRIGHT bins (DIM or MID) in which consecutive bins are separated by at
    most max_gap days without data. An episode qualifies if it holds >= min_dim DIM bins."""
    eps = []
    cur = None
    for b in bins:
        if b["state"] == "BRIGHT":
            if cur:
                eps.append(cur)
            cur = None
            continue
        if cur and b["t_lo"] - cur["t_end"] <= max_gap:
            cur["bins"].append(b)
            cur["t_end"] = b["t_hi"]
        else:
            if cur:
                eps.append(cur)
            cur = {"t_start": b["t_lo"], "t_end": b["t_hi"], "bins": [b]}
    if cur:
        eps.append(cur)
    out = []
    for e in eps:
        nd = sum(1 for b in e["bins"] if b["state"] == "DIM")
        mr = min(b["med_r"] for b in e["bins"])
        first_dim = next((b["t_lo"] for b in e["bins"] if b["state"] == "DIM"), None)
        out.append({"t_start": e["t_start"], "t_end": e["t_end"], "date_start": date_str(e["t_start"]),
                    "date_end": date_str(e["t_end"]), "n_bins": len(e["bins"]), "n_dim": nd,
                    "n_mid": len(e["bins"]) - nd, "min_r": mr, "max_dmag": dmag(mr),
                    "first_dim": first_dim, "qualifies": nd >= min_dim})
    return out


def overlaps(ep, a, b):
    return ep["t_start"] < b and ep["t_end"] > a


def control(eps):
    """Known answer: a qualifying episode must overlap the 2013-2014 event."""
    hit = [e for e in eps if e["qualifies"] and overlaps(e, T2_START, T2_END)]
    return ("PASS" if hit else "FAIL"), hit


def count_state(bins, state, a, b):
    return sum(1 for x in bins if x["state"] == state and a <= x["t_mid"] < b)


def main_outcome(bins, eps):
    post = [e for e in eps if e["qualifies"] and e["t_end"] >= J_2024_07]
    recent = [e for e in post if e["t_end"] > J_2025_07]
    info = {"n_post_episodes": len(post)}
    if recent:
        on = [e for e in recent if overlaps(e, J_2026_01, J_2027_01)]
        label = "RETURNED_ON_SCHEDULE" if on else "RETURNED_EARLY"
        ep = (on or recent)[0]
        info.update({"episode": ep, "onset": date_str(ep["first_dim"]) if ep["first_dim"] else None})
        return label, info
    if post:
        info["episode"] = post[0]
        return "DIMMED_BEFORE_2025_07", info
    shallow = [e for e in eps if e["t_end"] > J_2025_07 and e["n_bins"] >= 3]
    a = count_state(bins, "BRIGHT", J_2026_02, J_2026_05)
    b = count_state(bins, "BRIGHT", J_2026_08, float("inf"))
    c = count_state(bins, "BRIGHT", J_2025_10, float("inf"))
    d = count_state(bins, "DIM", J_2025_07, float("inf"))
    info.update({"bright_feb_apr_2026": a, "bright_from_aug_2026": b, "bright_from_oct_2025": c,
                 "dim_bins_from_jul_2025": d, "shallow_episodes": shallow})
    if shallow:
        return "SHALLOW_ONLY", info
    if a >= 1 and b >= 1 and c >= 4:
        return "NOT_RETURNED", info
    return "INCONCLUSIVE_COVERAGE", info


def bracket_core(bins):
    """Bins in the window every bracketed period predicts to be dim (CORE_LO..CORE_HI)."""
    sel = [b for b in bins if CORE_LO <= b["t_mid"] < CORE_HI]
    n = len(sel)
    nd = sum(1 for b in sel if b["state"] == "DIM")
    nb = sum(1 for b in sel if b["state"] == "BRIGHT")
    if n >= 1 and nd / n >= 0.5:
        lab = "BRACKET_DIM"
    elif n >= 2 and nb == n:
        lab = "BRACKET_BRIGHT"
    else:
        lab = "BRACKET_UNDETERMINED"
    return lab, {"date_lo": date_str(CORE_LO), "date_hi": date_str(CORE_HI), "n_bins": n, "n_dim": nd,
                 "n_bright": nb, "dmag": [round(b["dmag"], 3) for b in sel]}


def period_map(bins, grid=P_GRID, half=HALF_DURATION):
    rows = []
    for P in grid:
        tc = T2_CENTER + P * DAYS_PER_YEAR
        inb = [b for b in bins if tc - half <= b["t_mid"] <= tc + half]
        n = len(inb)
        nb = sum(1 for b in inb if b["state"] == "BRIGHT")
        nd = sum(1 for b in inb if b["state"] == "DIM")
        if n < 3:
            st = "open"
        elif nd / n >= 0.5:
            st = "supported"
        elif nb >= 3 and nd / n < 0.25:
            st = "excluded"
        else:
            st = "open"
        rows.append({"P": P, "t_center": tc, "date_center": date_str(tc), "n": n, "n_bright": nb,
                     "n_dim": nd, "state": st})
    return rows


def fraction(rows, state, lo=P_STATED[0], hi=P_STATED[1]):
    sel = [r for r in rows if lo - 1e-9 <= r["P"] <= hi + 1e-9]
    return sum(1 for r in sel if r["state"] == state) / len(sel) if sel else float("nan")


def runs(rows, state):
    out = []
    cur = None
    for r in rows:
        if r["state"] == state:
            if cur is None:
                cur = [r["P"], r["P"]]
            else:
                cur[1] = r["P"]
        elif cur is not None:
            out.append(tuple(cur))
            cur = None
    if cur is not None:
        out.append(tuple(cur))
    return out


def chart_period(rows):
    """0.25-yr bins from 11.0 to 14.0: fractions excluded / supported / open."""
    out = []
    for k in range(12):
        lo = 11.0 + 0.25 * k
        hi = lo + 0.25
        sel = [r for r in rows if lo - 1e-9 <= r["P"] < hi - 1e-9]
        n = len(sel)
        out.append((lo, hi, n,
                    sum(1 for r in sel if r["state"] == "excluded") / n,
                    sum(1 for r in sel if r["state"] == "supported") / n))
    return out


# --- ZTF (secondary, not in the verdict: DM Ori is near ZTF saturation when bright) ------------
def ztf_bins(text, width=BIN_DAYS, origin=BIN_ORIGIN):
    """Per (filter, 30-d bin): all rows, unflagged rows, median mag of unflagged rows."""
    rd = csv.DictReader(io.StringIO(text))
    g = {}
    for row in rd:
        try:
            t = float(row["hjd"])
            mag = float(row["mag"])
            flags = int(float(row["catflags"]))
        except (KeyError, ValueError, TypeError):
            continue
        t = t if t > 2400000 else t + 2400000.5
        k = int(math.floor((t - origin) / width))
        e = g.setdefault((row.get("filtercode", "?"), k), {"all": 0, "mags": []})
        e["all"] += 1
        if flags == 0 and mag == mag:
            e["mags"].append(mag)
    out = []
    for (fc, k) in sorted(g, key=lambda x: (x[1], x[0])):
        e = g[(fc, k)]
        lo = origin + k * width
        out.append({"filter": fc, "t_mid": lo + width / 2.0, "date": date_str(lo + width / 2.0),
                    "n_all": e["all"], "n_clean": len(e["mags"]),
                    "median_mag_clean": median(e["mags"]) if e["mags"] else None})
    return out
