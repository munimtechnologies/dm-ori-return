# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""DM Ori return test (pre-registered 2026-11-24, PREREGISTRATION.md).

python job.py --data-dir /data --out-dir /out
Primary data: ASAS-SN Sky Patrol v2 (V 2012-2018, g 2017-2026). Secondary, not in the verdict:
ZTF (saturates when DM Ori is bright) and TESS full-frame-image cutouts (sector flux levels).
Writes /out/results.json (rewritten after every stage, key "stage"), CSVs and summary.md.
"""
import argparse
import csv
import json
import os
import sys
import time
import traceback
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import dmori  # noqa: E402

START = time.time()
ZTF_URL = ("https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?POS=CIRCLE%20{ra:.6f}%20{dec:.6f}"
           "%20{rad:.6f}&BANDNAME=g,r,i&FORMAT=csv")


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


def get_asassn(data_dir, res):
    import asassn_io
    out = {}
    res["data"]["asassn"] = out
    try:
        out["client"] = asassn_io.ensure_client()
        path = os.path.join(data_dir, "asassn_lc.csv")
        info = None
        for k in range(3):
            try:
                info = asassn_io.fetch(dmori.RA, dmori.DEC, path, log=log)
                if info.get("status") == "ok":
                    break
            except Exception as exc:
                info = {"status": "error", "error": repr(exc)}
                log("ASAS-SN try %d: %r" % (k + 1, exc))
            time.sleep(60)
        out["fetch"] = info
        if not info or info.get("status") != "ok":
            return None
        by_band, lmeta = asassn_io.load_points(path)
        out["load"] = lmeta
        return by_band
    except Exception as exc:
        out["error"] = repr(exc)
        out["traceback"] = traceback.format_exc()[-2000:]
        log("ASAS-SN failed: %r" % exc)
        return None


def get_ztf(data_dir, out_dir, res):
    info = {}
    res["data"]["ztf"] = info
    url = ZTF_URL.format(ra=dmori.RA, dec=dmori.DEC, rad=0.0006)
    info["url"] = url
    text = None
    for k in range(4):
        try:
            with urllib.request.urlopen(url, timeout=300) as resp:
                text = resp.read().decode("utf-8", "replace")
            break
        except Exception as exc:
            info["error_try_%d" % (k + 1)] = repr(exc)
            time.sleep(30 * (k + 1))
    if text is None:
        info["status"] = "failed"
        return
    with open(os.path.join(data_dir, "ztf_lc.csv"), "w") as fh:
        fh.write(text)
    rows = dmori.ztf_bins(text)
    info.update({"status": "ok", "n_bins": len(rows), "n_rows": sum(r["n_all"] for r in rows),
                 "n_clean": sum(r["n_clean"] for r in rows)})
    write_csv(os.path.join(out_dir, "ztf_bins.csv"),
              ["filter", "bin_mid_jd", "date", "n_all", "n_clean", "median_mag_clean"],
              [(r["filter"], "%.1f" % r["t_mid"], r["date"], r["n_all"], r["n_clean"],
                "" if r["median_mag_clean"] is None else "%.3f" % r["median_mag_clean"]) for r in rows])
    log("ZTF: %d rows, %d clean, %d bins" % (info["n_rows"], info["n_clean"], len(rows)))


def get_tess(data_dir, out_dir, res):
    import io_tess
    info = {"sectors": {}}
    res["data"]["tess"] = info
    try:
        secs = io_tess.sectors(dmori.RA, dmori.DEC, log=log)
    except Exception as exc:
        info["error"] = "sector list: %r" % exc
        return
    info["sector_list"] = secs
    rows = []
    for s in secs:
        try:
            path = io_tess.fetch_cutout(dmori.RA, dmori.DEC, s, os.path.join(data_dir, "tesscut"), log=log)
            ph = io_tess.photometry(path)
        except Exception as exc:
            info["sectors"][str(s)] = {"error": repr(exc)}
            log("TESS sector %d failed: %r" % (s, exc))
            continue
        if ph["n_good"] < 50:
            info["sectors"][str(s)] = {"note": "too few points", "n_good": ph["n_good"]}
            continue
        t = ph["bjd"]
        m3 = dmori.median(ph["f3"])
        m1 = dmori.median(ph["f1"])
        rows.append({"sector": s, "t_mid": dmori.median(t), "date": dmori.date_str(dmori.median(t)),
                     "median_f3": m3, "median_f1": m1, "n": ph["n_good"],
                     "frac_aperture_in_peak": ph["frac_aperture_in_peak"]})
        info["sectors"][str(s)] = rows[-1]
    pre = [r for r in rows if r["t_mid"] < dmori.J_2025_07]
    ref3 = dmori.median([r["median_f3"] for r in pre]) if pre else None
    ref1 = dmori.median([r["median_f1"] for r in pre]) if pre else None
    for r in rows:
        r["ratio_f3"] = r["median_f3"] / ref3 if ref3 else None
        r["ratio_f1"] = r["median_f1"] / ref1 if ref1 else None
    info["reference_sectors"] = [r["sector"] for r in pre]
    write_csv(os.path.join(out_dir, "tess_sectors.csv"),
              ["sector", "date", "n", "median_f3", "median_f1", "ratio_f3_to_pre2025_07", "ratio_f1_to_pre2025_07"],
              [(r["sector"], r["date"], r["n"], "%.2f" % r["median_f3"], "%.2f" % r["median_f1"],
                "" if r["ratio_f3"] is None else "%.4f" % r["ratio_f3"],
                "" if r["ratio_f1"] is None else "%.4f" % r["ratio_f1"]) for r in rows])
    log("TESS sectors: %s" % ", ".join("%d %s r3=%s" % (r["sector"], r["date"], r["ratio_f3"]) for r in rows))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="/data")
    ap.add_argument("--out-dir", default="/out")
    a = ap.parse_args()
    os.makedirs(a.data_dir, exist_ok=True)
    os.makedirs(a.out_dir, exist_ok=True)
    rpath = os.path.join(a.out_dir, "results.json")
    res = {"version": dmori.VERSION, "target": "DM Ori", "ra": dmori.RA, "dec": dmori.DEC,
           "stage": "start", "data": {}, "preregistration": "PREREGISTRATION.md (2026-11-24)",
           "constants": {"T2_START": dmori.T2_START, "T2_END": dmori.T2_END, "T2_CENTER": dmori.T2_CENTER,
                         "PRED": dmori.PRED, "DIM_R": dmori.DIM_R, "BRIGHT_R": dmori.BRIGHT_R,
                         "REF_WINDOWS": dmori.REF_WINDOWS}}
    write_json(rpath, res)

    by_band = get_asassn(a.data_dir, res)
    res["stage"] = "asassn"
    write_json(rpath, res)
    if not by_band:
        res["outcome"] = "NO_DATA"
    else:
        cl = dmori.clean(by_band)
        refs, rmeta = dmori.references(cl)
        norm = dmori.normalise(cl, refs)
        res["data"]["asassn"]["n_clean"] = {b: len(v) for b, v in cl.items()}
        res["data"]["asassn"]["references"] = rmeta
        res["data"]["asassn"]["n_normalised"] = len(norm)
        if norm:
            res["data"]["asassn"]["first"] = dmori.date_str(norm[0][0])
            res["data"]["asassn"]["last"] = dmori.date_str(norm[-1][0])
        write_csv(os.path.join(a.out_dir, "asassn_points.csv"), ["jd", "date", "r", "dmag", "band", "camera"],
                  [("%.5f" % t, dmori.date_str(t), "%.4f" % r, "%.3f" % dmori.dmag(r), b, c) for t, r, b, c in norm])
        bins = dmori.bin_points(norm)
        eps = dmori.episodes(bins)
        write_csv(os.path.join(a.out_dir, "bins.csv"),
                  ["bin_start_jd", "date_start", "decimal_year_mid", "n", "median_r", "dmag", "state", "bands"],
                  [("%.1f" % b["t_lo"], dmori.date_str(b["t_lo"]), "%.3f" % dmori.decimal_year(b["t_mid"]),
                    b["n"], "%.4f" % b["med_r"], "%.3f" % b["dmag"], b["state"], b["bands"]) for b in bins])
        write_csv(os.path.join(a.out_dir, "chart_lightcurve.csv"), ["year", "dmag_below_quiescence"],
                  [("%.3f" % dmori.decimal_year(b["t_mid"]), "%.3f" % b["dmag"]) for b in bins])
        write_csv(os.path.join(a.out_dir, "episodes.csv"),
                  ["date_start", "date_end", "n_bins", "n_dim", "n_mid", "max_dmag", "qualifies"],
                  [(e["date_start"], e["date_end"], e["n_bins"], e["n_dim"], e["n_mid"], "%.3f" % e["max_dmag"],
                    int(e["qualifies"])) for e in eps])
        ctrl, hit = dmori.control(eps)
        label, info = dmori.main_outcome(bins, eps)
        blab, binfo = dmori.bracket_core(bins)
        res["bracket_core"] = {"label": blab, "info": binfo,
                               "P_ingress_days": [dmori.P_ING_LO, dmori.P_ING_HI]}
        pm = dmori.period_map(bins)
        write_csv(os.path.join(a.out_dir, "period_map.csv"),
                  ["P_yr", "center_date", "n_bins", "n_bright", "n_dim", "state"],
                  [("%.2f" % r["P"], r["date_center"], r["n"], r["n_bright"], r["n_dim"], r["state"]) for r in pm])
        write_csv(os.path.join(a.out_dir, "chart_period_map.csv"),
                  ["P_lo_yr", "P_hi_yr", "n_trials", "frac_excluded", "frac_supported"],
                  [("%.2f" % x[0], "%.2f" % x[1], x[2], "%.3f" % x[3], "%.3f" % x[4]) for x in dmori.chart_period(pm)])
        res.update({"n_bins": len(bins), "n_bins_by_state": {s: sum(1 for b in bins if b["state"] == s)
                                                             for s in ("DIM", "MID", "BRIGHT")},
                    "episodes": eps, "control": {"verdict": ctrl, "episodes": hit},
                    "main": {"label": label, "info": info},
                    "period_map": {"excluded_frac_11_5_13_5": dmori.fraction(pm, "excluded"),
                                   "supported_frac_11_5_13_5": dmori.fraction(pm, "supported"),
                                   "excluded_runs": dmori.runs(pm, "excluded"),
                                   "supported_runs": dmori.runs(pm, "supported"),
                                   "open_runs": dmori.runs(pm, "open")}})
        res["outcome"] = label if ctrl == "PASS" else "CONTROL_FAILED_" + label
        log("control %s; main %s; excluded frac %.3f" % (ctrl, label, res["period_map"]["excluded_frac_11_5_13_5"]))
    res["stage"] = "verdict"
    write_json(rpath, res)

    for name, fn in (("ztf", lambda: get_ztf(a.data_dir, a.out_dir, res)),
                     ("tess", lambda: get_tess(a.data_dir, a.out_dir, res))):
        try:
            fn()
        except Exception as exc:
            res["data"].setdefault(name, {})["error"] = repr(exc)
            log("%s failed: %r" % (name, exc))
        res["stage"] = name
        write_json(rpath, res)

    res["stage"] = "done"
    res["runtime_s"] = time.time() - START
    write_json(rpath, res)
    with open(os.path.join(a.out_dir, "summary.md"), "w") as fh:
        fh.write("# DM Ori return test v%s\n\nOutcome: %s\n\n" % (dmori.VERSION, res.get("outcome")))
        if "main" in res:
            fh.write("Control (2013-14 event recovered): %s\n\n" % res["control"]["verdict"])
            fh.write("Main: %s\n\n%s\n\n" % (res["main"]["label"], json.dumps(res["main"]["info"], default=str)))
            fh.write("Bracket core (2026-02-21..2026-06-01): %s %s\n\n"
                     % (res["bracket_core"]["label"], json.dumps(res["bracket_core"]["info"])))
            fh.write("Period map 11.5-13.5 yr: excluded %.3f, supported %.3f\n"
                     % (res["period_map"]["excluded_frac_11_5_13_5"], res["period_map"]["supported_frac_11_5_13_5"]))
    log("OUTCOME %s" % res.get("outcome"))


if __name__ == "__main__":
    main()
