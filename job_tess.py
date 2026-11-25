# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""DM Ori TESS arm (AMENDMENT-2026-11-25-tess.md).

python job_tess.py --data-dir /data --out-dir /out
Needs numpy and astropy. Fetches 21x21 TESScut cutouts for sectors 6, 98 and 110, measures DM Ori
and comparison stars matched by sky position, calibrates sector-6 TESS flux against the
simultaneous ASAS-SN g points (dmtess.ASASSN_S6) and grades sectors 98 and 110 with dmtess rules.
"""
import argparse
import csv
import glob
import io
import json
import os
import sys
import time
import traceback
import urllib.request
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import dmtess  # noqa: E402

START = time.time()
CUT_URL = ("https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra={ra:.6f}&dec={dec:.6f}"
           "&y={n}&x={n}&units=px&sector={s}")
SECTOR_URL = "https://mast.stsci.edu/tesscut/api/v0.1/sector?ra={ra:.6f}&dec={dec:.6f}&radius=0m"
BAD_BITS = 1 | 2 | 4 | 8 | 32 | 128
NPIX = 21
MAX_COMPS = 12


def log(msg):
    print("[%7.1fs] %s" % (time.time() - START, msg), flush=True)


def write_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, default=str)
    os.replace(tmp, path)


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def get(url, tries=4, timeout=900):
    last = None
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:
            last = repr(exc)
            log("GET try %d failed (%s): %s" % (k + 1, url[:90], last))
            time.sleep(30 * (k + 1))
    raise RuntimeError(last)


def fetch_cutout(sector, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    have = glob.glob(os.path.join(out_dir, "*s%04d*_%dx%d_*.fits" % (sector, NPIX, NPIX)))
    if have:
        return have[0]
    blob = get(CUT_URL.format(ra=dmtess.RA, dec=dmtess.DEC, n=NPIX, s=sector))
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        names = [x for x in z.namelist() if x.endswith(".fits")]
        z.extractall(out_dir, members=names)
    log("TESScut sector %d: %s" % (sector, ", ".join(names)))
    return os.path.join(out_dir, names[0])


def load(path):
    import numpy as np
    from astropy.io import fits
    from astropy.wcs import WCS
    with fits.open(path, memmap=False) as h:
        d = h[1].data
        t = np.array(d["TIME"], dtype=float)
        flux = np.array(d["FLUX"], dtype=np.float32)
        q = np.array(d["QUALITY"], dtype=int)
        wcs = WCS(h[2].header)
    nt, ny, nx = flux.shape
    good = np.isfinite(t) & ((q & BAD_BITS) == 0) & np.isfinite(flux[:, ny // 2, nx // 2])
    t, flux = t[good], flux[good]
    med = np.nanmedian(flux, axis=0)
    cy, cx = ny // 2, nx // 2
    sub = med[cy - 1:cy + 2, cx - 1:cx + 2]
    dy, dx = np.unravel_index(int(np.nanargmax(sub)), sub.shape)
    py, px = cy - 1 + dy, cx - 1 + dx
    yy, xx = np.mgrid[0:ny, 0:nx]
    dist = np.hypot(yy - py, xx - px)
    thr = np.nanpercentile(med, 60)
    bmask = (dist >= 4) & (med < thr) & np.isfinite(med)
    bkg = np.nanmedian(flux[:, bmask], axis=1)
    return {"t": t + 2457000.0, "flux": flux, "bkg": bkg, "med": med, "bkg_med": float(np.nanmedian(bkg)),
            "target_px": (int(py), int(px)), "wcs": wcs, "shape": (ny, nx), "n_good": int(len(t)),
            "n_all": int(nt), "n_bkg_pixels": int(bmask.sum())}


def ap_flux(S, py, px):
    import numpy as np
    blk = S["flux"][:, py - 1:py + 2, px - 1:px + 2].reshape(len(S["t"]), 9)
    return np.nansum(blk, axis=1) - 9.0 * S["bkg"]


def local_max(med, y, x):
    v = med[y, x]
    return all(v > med[y + a, x + b] for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0))


def find_comps(S):
    """Local maxima in the reference median image, 0.2-5x the target peak (above background),
    >= 4 px from the target and >= 2 px from the edge; brightest MAX_COMPS kept."""
    import numpy as np
    med = S["med"] - S["bkg_med"]
    ny, nx = S["shape"]
    ty, tx = S["target_px"]
    tpk = med[ty, tx]
    cands = []
    for y in range(2, ny - 2):
        for x in range(2, nx - 2):
            if not np.isfinite(med[y, x]) or (abs(y - ty) < 4 and abs(x - tx) < 4):
                continue
            if 0.2 * tpk <= med[y, x] <= 5.0 * tpk and local_max(med, y, x):
                cands.append((float(med[y, x]), y, x))
    cands.sort(reverse=True)
    out = []
    for v, y, x in cands[:MAX_COMPS]:
        ra, dec = S["wcs"].all_pix2world([[x, y]], 0)[0]
        out.append({"id": len(out), "ra": float(ra), "dec": float(dec), "peak_ref": v})
    return out


def match(S, comp):
    """Pixel of the comparison star in sector S: the brightest pixel within +-1 of the WCS position;
    None if off the usable area or not a local maximum."""
    import numpy as np
    x, y = S["wcs"].all_world2pix([[comp["ra"], comp["dec"]]], 0)[0]
    xi, yi = int(round(float(x))), int(round(float(y)))
    ny, nx = S["shape"]
    best = None
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            yy, xx = yi + a, xi + b
            if 2 <= yy <= ny - 3 and 2 <= xx <= nx - 3 and np.isfinite(S["med"][yy, xx]):
                if best is None or S["med"][yy, xx] > S["med"][best[0], best[1]]:
                    best = (yy, xx)
    if best is None or not local_max(S["med"], best[0], best[1]):
        return None
    ty, tx = S["target_px"]
    if abs(best[0] - ty) < 4 and abs(best[1] - tx) < 4:
        return None
    return best


def main():
    import numpy as np
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    rpath = os.path.join(a.out_dir, "results_tess.json")
    res = {"version": dmtess.VERSION, "target": "DM Ori", "rules": "AMENDMENT-2026-11-25-tess.md",
           "stage": "start", "sectors": {}}
    write_json(rpath, res)

    try:
        secs = sorted(int(r["sector"]) for r in json.loads(
            get(SECTOR_URL.format(ra=dmtess.RA, dec=dmtess.DEC)).decode()).get("results", []))
    except Exception as exc:
        secs = [dmtess.REF_SECTOR] + list(dmtess.TARGET_SECTORS)
        res["sector_list_error"] = repr(exc)
    res["sector_list"] = secs
    log("sectors %s" % secs)

    data = {}
    for s in secs:
        try:
            S = load(fetch_cutout(s, os.path.join(a.data_dir, "tesscut21")))
            data[s] = S
            x, y = S["wcs"].all_world2pix([[dmtess.RA, dmtess.DEC]], 0)[0]
            S["wcs_check_px"] = [float(y), float(x)]
            log("sector %d: %d good cadences, target pixel %s, WCS puts DM Ori at row %.2f col %.2f"
                % (s, S["n_good"], S["target_px"], float(y), float(x)))
        except Exception as exc:
            res["sectors"][str(s)] = {"error": repr(exc), "traceback": traceback.format_exc()[-1500:]}
            log("sector %d failed: %r" % (s, exc))
    if dmtess.REF_SECTOR not in data:
        res["outcome"] = "NO_REFERENCE_SECTOR"
        write_json(rpath, res)
        log("OUTCOME NO_REFERENCE_SECTOR")
        return

    ref = data[dmtess.REF_SECTOR]
    comps = find_comps(ref)
    res["n_comps_found"] = len(comps)
    comp_rows = []
    comp_med = {}
    for s, S in sorted(data.items()):
        for c in comps:
            pos = match(S, c)
            if pos is None:
                continue
            f = ap_flux(S, pos[0], pos[1])
            m = float(np.nanmedian(f))
            if m > 0:
                comp_med[(s, c["id"])] = m
                comp_rows.append((s, c["id"], "%.6f" % c["ra"], "%.6f" % c["dec"], pos[0], pos[1], "%.2f" % m))
    write_csv(os.path.join(a.out_dir, "tess_comps.csv"),
              ["sector", "comp_id", "ra", "dec", "row", "col", "median_flux"], comp_rows)

    # target series and sector levels
    tgt = {}
    level_rows = []
    for s, S in sorted(data.items()):
        f = ap_flux(S, *S["target_px"])
        ok = np.isfinite(f)
        tgt[s] = (S["t"][ok], f[ok])
        ratios = [comp_med[(s, c["id"])] / comp_med[(dmtess.REF_SECTOR, c["id"])] for c in comps
                  if (s, c["id"]) in comp_med and (dmtess.REF_SECTOR, c["id"]) in comp_med]
        k, sysf = dmtess.comp_scale(ratios) if ratios else (None, None)
        m = float(np.nanmedian(f[ok]))
        info = {"date_mid": dmtess.jd_to_date(float(np.median(S["t"][ok]))),
                "t_mid": float(np.median(S["t"][ok])), "t_first": float(S["t"][ok][0]),
                "t_last": float(S["t"][ok][-1]), "n": int(ok.sum()), "target_median_flux": m,
                "n_comps": len(ratios), "comp_ratios": ratios, "k": k, "sys": sysf,
                "raw_ratio_to_ref": None, "target_px": list(S["target_px"]),
                "wcs_check_px": S.get("wcs_check_px"), "n_bkg_pixels": S["n_bkg_pixels"]}
        res["sectors"][str(s)] = info
    mref = res["sectors"][str(dmtess.REF_SECTOR)]["target_median_flux"]
    for s in res["sectors"]:
        if "target_median_flux" in res["sectors"][s]:
            res["sectors"][s]["raw_ratio_to_ref"] = res["sectors"][s]["target_median_flux"] / mref
    res["stage"] = "photometry"
    write_json(rpath, res)

    # calibration against ASAS-SN g in sector 6
    t6, f6 = tgt[dmtess.REF_SECTOR]
    prs = dmtess.pairs(dmtess.ASASSN_S6, t6.tolist(), f6.tolist())
    cal = dmtess.calibrate(prs)
    res["calibration"] = cal
    write_csv(os.path.join(a.out_dir, "calibration_pairs.csv"), ["jd", "date", "asassn_g_rel_flux", "tess_flux", "n_tess", "camera"],
              [("%.5f" % p["jd"], dmtess.jd_to_date(p["jd"]), "%.4f" % p["r"], "%.2f" % p["f"], p["n"], p["camera"])
               for p in prs])
    log("calibration: %s" % json.dumps({k2: v for k2, v in cal.items()}, default=str))

    comps_ok = all(res["sectors"].get(str(s), {}).get("n_comps", 0) >= dmtess.MIN_COMPS and
                   (res["sectors"][str(s)].get("sys") or 1.0) <= dmtess.MAX_SYS
                   for s in dmtess.TARGET_SECTORS if "target_median_flux" in res["sectors"].get(str(s), {}))
    res["comps_ok"] = comps_ok

    labels = {}
    daily_rows = []
    if cal["ok"]:
        A, B = cal["a"], cal["b"]
        for s in sorted(tgt):
            info = res["sectors"][str(s)]
            k = info["k"] if s != dmtess.REF_SECTOR else 1.0
            sysf = info["sys"] if s != dmtess.REF_SECTOR else 0.0
            if k is None:
                continue
            lev = dmtess.sector_level(info["target_median_flux"], k, sysf, A, B)
            lev["extrapolated"] = bool(lev["r_equiv"] < cal["r_min"] or lev["r_equiv"] > cal["r_max"])
            info["level"] = lev
            if s in dmtess.TARGET_SECTORS:
                labels[s] = lev["label"]
            # no-comparison robustness (reported only)
            info["level_no_comps"] = dmtess.sector_level(info["target_median_flux"], 1.0, dmtess.SYS_FLOOR, A, B)
            t, f = tgt[s]
            day = np.floor(t - 0.5)
            for dd in np.unique(day):
                sel = day == dd
                if sel.sum() < 10:
                    continue
                F = float(np.median(f[sel])) / k
                daily_rows.append((s, "%.1f" % (dd + 0.5), dmtess.jd_to_date(dd + 0.5), int(sel.sum()), "%.2f" % F,
                                   "%.4f" % dmtess.equiv_r(F, A, B), "%.3f" % dmtess.equiv_dmag(F, A, B)))
    write_csv(os.path.join(a.out_dir, "tess_daily.csv"),
              ["sector", "jd_day_start", "date", "n", "flux_on_s6_scale", "asassn_g_equiv_rel_flux", "dmag_equiv"],
              daily_rows)
    write_csv(os.path.join(a.out_dir, "tess_sector_levels.csv"),
              ["sector", "date_mid", "target_median_flux", "n_comps", "k", "sys", "raw_ratio_to_ref",
               "dmag_equiv", "dmag_lo", "dmag_hi", "label", "dmag_equiv_no_comps"],
              [(s, v.get("date_mid"), "%.2f" % v["target_median_flux"], v.get("n_comps"),
                "" if v.get("k") is None else "%.4f" % v["k"], "" if v.get("sys") is None else "%.4f" % v["sys"],
                "%.4f" % v["raw_ratio_to_ref"],
                "%.3f" % v["level"]["dmag"] if "level" in v else "", "%.3f" % v["level"]["dmag_lo"] if "level" in v else "",
                "%.3f" % v["level"]["dmag_hi"] if "level" in v else "", v["level"]["label"] if "level" in v else "",
                "%.3f" % v["level_no_comps"]["dmag"] if "level_no_comps" in v else "")
               for s, v in sorted(res["sectors"].items(), key=lambda kv: int(kv[0])) if "target_median_flux" in v])

    if not cal["ok"]:
        res["outcome"] = "CALIBRATION_FAILED"
    elif not comps_ok:
        res["outcome"] = "COMPS_FAILED"
    else:
        snaps = [(res["sectors"][str(s)]["t_mid"], labels[s]) for s in dmtess.TARGET_SECTORS if s in labels]
        rows, frac = dmtess.ingress_map(snaps)
        write_csv(os.path.join(a.out_dir, "ingress_map.csv"), ["p_days", "ingress_jd", "ingress_date", "state"],
                  [(r["p"], "%.1f" % r["ti"], dmtess.jd_to_date(r["ti"]), r["state"]) for r in rows])
        res["ingress_map"] = {"fractions": frac,
                              "excluded_runs": [[dmtess.jd_to_date(x), dmtess.jd_to_date(y)] for x, y in dmtess.runs(rows, "excluded")],
                              "consistent_runs": [[dmtess.jd_to_date(x), dmtess.jd_to_date(y)] for x, y in dmtess.runs(rows, "consistent")],
                              "open_runs": [[dmtess.jd_to_date(x), dmtess.jd_to_date(y)] for x, y in dmtess.runs(rows, "open")]}
        res["labels"] = {str(s): v for s, v in labels.items()}
        res["outcome"] = dmtess.headline(labels.get(98), labels.get(110))
    res["stage"] = "done"
    res["runtime_s"] = time.time() - START
    write_json(rpath, res)
    with open(os.path.join(a.out_dir, "summary_tess.md"), "w") as fh:
        fh.write("# DM Ori TESS arm %s\n\nOutcome: %s\n\n" % (dmtess.VERSION, res["outcome"]))
        fh.write("Calibration: %s\n\n" % json.dumps(cal, default=str))
        for s, v in sorted(res["sectors"].items(), key=lambda kv: int(kv[0])):
            fh.write("Sector %s: %s\n\n" % (s, json.dumps({k2: v.get(k2) for k2 in
                                                          ("date_mid", "n_comps", "k", "sys", "raw_ratio_to_ref",
                                                           "level", "level_no_comps")}, default=str)))
        if "ingress_map" in res:
            fh.write("Ingress map: %s\n" % json.dumps(res["ingress_map"]))
    log("OUTCOME %s" % res["outcome"])


if __name__ == "__main__":
    main()
