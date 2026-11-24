# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
"""Tests for the frozen DM Ori rules (stdlib only). Synthetic times sit mid-bin, far from edges."""
import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import dmori  # noqa: E402
from dmori import jd  # noqa: E402

O, W = dmori.BIN_ORIGIN, dmori.BIN_DAYS
SEASON = (1, 2, 3, 4, 8, 9, 10, 11, 12)  # Orion is observable roughly August to April


def build(r_func, t_stop, start=None):
    start = start or jd(2012, 1, 1)
    by = {"V": [], "g": []}
    for k in range(int(math.floor((start - O) / W)), int(math.floor((t_stop - O) / W))):
        mid = O + W * k + W / 2.0
        if int(dmori.date_str(mid)[5:7]) not in SEASON:
            continue
        r = r_func(mid)
        band = "V" if mid < jd(2018, 6, 1) else "g"
        for off in (12.0, 14.0, 16.0, 18.0):
            by[band].append((O + W * k + off, 1000.0 * r, 10.0, "bA"))
    return by


def in2013(t):
    return dmori.T2_START + 40 <= t <= dmori.T2_END - 40


def run(by):
    cl = dmori.clean(by)
    refs, _ = dmori.references(cl)
    norm = dmori.normalise(cl, refs)
    bins = dmori.bin_points(norm)
    eps = dmori.episodes(bins)
    ctrl, _ = dmori.control(eps)
    label, info = dmori.main_outcome(bins, eps)
    return bins, eps, ctrl, label, info


STOP = jd(2026, 11, 20)


class TestDates(unittest.TestCase):
    def test_constants(self):
        self.assertEqual(jd(2000, 1, 1), 2451544.5)
        self.assertEqual(dmori.T2_START, 2456320.5)
        self.assertEqual(dmori.T2_END, 2456901.5)
        self.assertEqual(dmori.J_2026_01, 2461041.5)
        self.assertEqual(dmori.date_str(dmori.PRED), "2026-07-01")
        self.assertEqual(dmori.P_ING_LO, 4452.0)
        self.assertEqual(dmori.P_ING_HI, 4692.0)
        self.assertEqual(dmori.date_str(dmori.CORE_LO), "2026-02-21")
        self.assertEqual(dmori.date_str(dmori.CORE_HI), "2026-06-01")
        self.assertEqual(dmori.date_str(jd(2013, 1, 28) + 0.4), "2013-01-28")
        self.assertAlmostEqual(dmori.decimal_year(jd(2026, 7, 2)), 2026.5, places=2)

    def test_thresholds(self):
        self.assertAlmostEqual(dmori.dmag(dmori.DIM_R), 0.75, places=9)
        self.assertEqual(dmori.state_of(0.25), "DIM")
        self.assertEqual(dmori.state_of(0.65), "MID")
        self.assertEqual(dmori.state_of(1.0), "BRIGHT")


class TestNormalise(unittest.TestCase):
    def test_clean_snr(self):
        out = dmori.clean({"g": [(1.0, 100.0, 30.0, "a"), (2.0, 100.0, 10.0, "a"), (3.0, -5.0, 1.0, "a")]})
        self.assertEqual([p[0] for p in out["g"]], [2.0])

    def test_references(self):
        t0 = jd(2020, 3, 10) + 0.3
        rows = [(t0 + i, 100.0, 1.0, "bA") for i in range(25)]
        rows += [(t0 + 40 + i, 80.0, 1.0, "bB") for i in range(5)]
        rows += [(jd(2026, 3, 10) + 0.3 + i, 30.0, 1.0, "bC") for i in range(5)]
        v = [(jd(2016, 2, 10) + 0.3 + i, 50.0, 1.0, "bA") for i in range(5)]
        refs, meta = dmori.references({"g": rows, "V": v})
        self.assertEqual(refs[("g", "bA")], 100.0)
        self.assertEqual(refs[("g", "bB")], 100.0)  # too few own points: band median
        self.assertEqual(refs[("g", "bC")], 100.0)
        self.assertNotIn(("V", "bA"), refs)
        self.assertFalse(meta["V"]["used"])
        norm = dmori.normalise({"g": rows}, refs)
        self.assertEqual(sorted(set(round(p[1], 6) for p in norm)), [0.3, 0.8, 1.0])


class TestScenarios(unittest.TestCase):
    def test_on_schedule(self):
        f = lambda t: 0.25 if (in2013(t) or jd(2025, 11, 1) <= t <= jd(2027, 3, 1)) else 1.0  # noqa: E731
        bins, eps, ctrl, label, info = run(build(f, STOP))
        self.assertEqual(ctrl, "PASS")
        self.assertEqual(label, "RETURNED_ON_SCHEDULE")
        self.assertIn(info["onset"][:7], ("2025-10", "2025-11"))
        self.assertGreater(info["episode"]["max_dmag"], 1.4)
        self.assertEqual(dmori.bracket_core(bins)[0], "BRACKET_DIM")

    def test_not_returned_and_period_map(self):
        f = lambda t: 0.25 if in2013(t) else 1.0  # noqa: E731
        bins, eps, ctrl, label, info = run(build(f, STOP))
        self.assertEqual(ctrl, "PASS")
        self.assertEqual(label, "NOT_RETURNED")
        self.assertGreaterEqual(info["bright_feb_apr_2026"], 1)
        self.assertGreaterEqual(info["bright_from_aug_2026"], 1)
        self.assertEqual(dmori.bracket_core(bins)[0], "BRACKET_BRIGHT")
        self.assertEqual(dmori.bracket_core([])[0], "BRACKET_UNDETERMINED")
        pm = {r["P"]: r for r in dmori.period_map(bins)}
        self.assertEqual(pm[12.5]["state"], "excluded")
        self.assertEqual(pm[11.0]["state"], "excluded")
        self.assertEqual(pm[14.0]["state"], "open")
        rows = dmori.period_map(bins)
        self.assertGreater(dmori.fraction(rows, "excluded"), 0.5)
        self.assertEqual(dmori.fraction(rows, "supported"), 0.0)
        ch = dmori.chart_period(rows)
        self.assertEqual(len(ch), 12)
        self.assertEqual(sum(x[2] for x in ch), 300)

    def test_supported_period(self):
        f = lambda t: 0.25 if (in2013(t) or jd(2025, 11, 1) <= t <= jd(2027, 3, 1)) else 1.0  # noqa: E731
        bins, _, _, _, _ = run(build(f, STOP))
        pm = {r["P"]: r for r in dmori.period_map(bins)}
        self.assertEqual(pm[12.5]["state"], "supported")

    def test_early(self):
        f = lambda t: 0.25 if (in2013(t) or jd(2024, 9, 1) <= t <= jd(2025, 9, 30)) else 1.0  # noqa: E731
        self.assertEqual(run(build(f, STOP))[3], "RETURNED_EARLY")

    def test_dimmed_before_2025_07(self):
        f = lambda t: 0.25 if (in2013(t) or jd(2024, 8, 15) <= t <= jd(2025, 4, 30)) else 1.0  # noqa: E731
        self.assertEqual(run(build(f, STOP))[3], "DIMMED_BEFORE_2025_07")

    def test_shallow(self):
        f = lambda t: 0.25 if in2013(t) else (0.65 if t >= jd(2025, 11, 1) else 1.0)  # noqa: E731
        self.assertEqual(run(build(f, STOP))[3], "SHALLOW_ONLY")

    def test_coverage(self):
        f = lambda t: 0.25 if in2013(t) else 1.0  # noqa: E731
        self.assertEqual(run(build(f, jd(2025, 6, 1)))[3], "INCONCLUSIVE_COVERAGE")

    def test_control_fails_without_2013_event(self):
        _, _, ctrl, label, _ = run(build(lambda t: 1.0, STOP))
        self.assertEqual(ctrl, "FAIL")
        self.assertEqual(label, "NOT_RETURNED")


class TestEpisodes(unittest.TestCase):
    def _bin(self, t_lo, state, r):
        return {"t_lo": t_lo, "t_hi": t_lo + 30, "t_mid": t_lo + 15, "med_r": r, "state": state, "n": 4}

    def test_gap_splits(self):
        b = [self._bin(2460000.5, "DIM", 0.3), self._bin(2460230.5, "DIM", 0.3)]
        eps = dmori.episodes(b)
        self.assertEqual(len(eps), 2)
        self.assertFalse(any(e["qualifies"] for e in eps))

    def test_bright_splits_and_mid_joins(self):
        b = [self._bin(2460000.5, "DIM", 0.3), self._bin(2460030.5, "MID", 0.6),
             self._bin(2460060.5, "DIM", 0.3), self._bin(2460090.5, "BRIGHT", 1.0),
             self._bin(2460120.5, "DIM", 0.3)]
        eps = dmori.episodes(b)
        self.assertEqual([e["n_bins"] for e in eps], [3, 1])
        self.assertTrue(eps[0]["qualifies"])
        self.assertEqual(eps[0]["n_mid"], 1)

    def test_runs(self):
        rows = [{"P": 11.0, "state": "open"}, {"P": 11.01, "state": "excluded"},
                {"P": 11.02, "state": "excluded"}, {"P": 11.03, "state": "open"},
                {"P": 11.04, "state": "excluded"}]
        self.assertEqual(dmori.runs(rows, "excluded"), [(11.01, 11.02), (11.04, 11.04)])


class TestZtf(unittest.TestCase):
    def test_bins(self):
        t = jd(2026, 9, 1)
        k = int(math.floor((t - O) / W))
        base = O + W * k
        text = ("oid,hjd,mag,magerr,catflags,filtercode\n"
                "1,%.5f,13.10,0.02,0,zr\n1,%.5f,13.30,0.02,0,zr\n1,%.5f,11.20,0.02,32768,zr\n"
                % (base + 10.0, base + 14.0, base + 18.0))
        rows = dmori.ztf_bins(text)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["n_all"], 3)
        self.assertEqual(rows[0]["n_clean"], 2)
        self.assertAlmostEqual(rows[0]["median_mag_clean"], 13.2)


if __name__ == "__main__":
    unittest.main()
