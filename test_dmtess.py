# Copyright 2026 Munim, Inc.
# SPDX-License-Identifier: Apache-2.0
import math
import unittest

import dmtess as d


class TestStats(unittest.TestCase):
    def test_median_mad(self):
        self.assertEqual(d.median([3, 1, 2]), 2)
        self.assertEqual(d.median([4, 1, 2, 3]), 2.5)
        self.assertIsNone(d.median([]))
        self.assertEqual(d.mad([1, 2, 3, 4, 100]), 1)

    def test_pearson_ols(self):
        x = [0.0, 1.0, 2.0, 3.0]
        y = [1.0, 3.0, 5.0, 7.0]
        self.assertAlmostEqual(d.pearson(x, y), 1.0)
        a, b = d.ols(x, y)
        self.assertAlmostEqual(a, 1.0)
        self.assertAlmostEqual(b, 2.0)
        self.assertIsNone(d.pearson([1.0, 1.0], [2.0, 3.0]))


class TestCalibration(unittest.TestCase):
    def _series(self, a=2000.0, b=4000.0):
        # TESS flux follows a + b*r with r a step from 0.85 to 0.35 at t = 10.
        t = [i * 0.05 for i in range(0, 401)]
        f = [a + b * (0.85 if tv < 10 else 0.35) for tv in t]
        return t, f

    def test_pairs_and_fit(self):
        t, f = self._series()
        asassn = [(2.0, 0.85, "x"), (4.0, 0.85, "x"), (6.0, 0.85, "x"), (12.0, 0.35, "x"),
                  (14.0, 0.35, "x"), (16.0, 0.35, "x"), (50.0, 0.5, "x")]
        prs = d.pairs(asassn, t, f)
        self.assertEqual(len(prs), 6)  # the point at t = 50 has no TESS data
        self.assertTrue(all(p["n"] >= 3 for p in prs))
        cal = d.calibrate(prs)
        self.assertTrue(cal["ok"], cal)
        self.assertAlmostEqual(cal["a"], 2000.0, places=6)
        self.assertAlmostEqual(cal["b"], 4000.0, places=6)
        self.assertAlmostEqual(cal["f_quiescent"], 6000.0, places=6)
        self.assertAlmostEqual(cal["nonvarying_share"], 1.0 / 3.0, places=6)

    def test_gates(self):
        prs = [{"r": 0.8, "f": 10.0}] * 3
        self.assertEqual(d.calibrate(prs)["reason"], "too few pairs")
        prs = [{"r": 0.8 + 0.01 * i, "f": 10.0 + i} for i in range(8)]
        self.assertEqual(d.calibrate(prs)["reason"], "span too small")
        prs = [{"r": 0.3 + 0.1 * i, "f": 10.0 - i} for i in range(8)]
        self.assertEqual(d.calibrate(prs)["reason"], "slope not positive")
        f = [10, 30, 12, 31, 9, 33, 11, 35]
        prs = [{"r": 0.3 + 0.1 * i, "f": float(f[i])} for i in range(8)]
        self.assertFalse(d.calibrate(prs)["ok"])

    def test_embedded_asassn(self):
        self.assertEqual(len(d.ASASSN_S6), 19)
        jds = [p[0] for p in d.ASASSN_S6]
        self.assertEqual(jds, sorted(jds))
        self.assertTrue(all(0.2 < p[1] < 1.0 for p in d.ASASSN_S6))


class TestLevels(unittest.TestCase):
    def test_equiv_dmag(self):
        a, b = 2000.0, 4000.0
        self.assertAlmostEqual(d.equiv_dmag(6000.0, a, b), 0.0)
        self.assertAlmostEqual(d.equiv_dmag(2000.0 + 4000.0 * 0.5, a, b), 2.5 * math.log10(2.0))
        self.assertEqual(d.equiv_dmag(1999.0, a, b), d.DMAG_CAP)

    def test_comp_scale(self):
        k, sys = d.comp_scale([1.0, 1.0, 1.0, 1.0])
        self.assertEqual(k, 1.0)
        self.assertEqual(sys, d.SYS_FLOOR)
        k, sys = d.comp_scale([0.9, 1.0, 1.1, 1.0, 1.0])
        self.assertEqual(k, 1.0)
        self.assertAlmostEqual(sys, d.SYS_FLOOR)  # MAD is 0
        k, sys = d.comp_scale([0.8, 0.9, 1.0, 1.1, 1.2])
        self.assertAlmostEqual(sys, 1.4826 * 0.1, places=9)

    def test_classify(self):
        self.assertEqual(d.classify(0.8, 1.0), "DEEP")
        self.assertEqual(d.classify(-0.1, 0.2), "BRIGHT")
        self.assertEqual(d.classify(0.4, 0.6), "MID")
        self.assertEqual(d.classify(0.25, 0.5), "UNCERTAIN")
        self.assertEqual(d.classify(0.6, 0.9), "UNCERTAIN")

    def test_sector_level(self):
        a, b = 2000.0, 4000.0
        # true r = 0.25 (1.505 mag) seen through a scale k = 0.8
        m = (a + b * 0.25) * 0.8
        lev = d.sector_level(m, 0.8, 0.02, a, b)
        self.assertAlmostEqual(lev["r_equiv"], 0.25, places=9)
        self.assertEqual(lev["label"], "DEEP")
        self.assertLess(lev["dmag_lo"], lev["dmag"])
        self.assertGreater(lev["dmag_hi"], lev["dmag"])
        lev = d.sector_level(6000.0 * 1.1, 1.1, 0.02, a, b)
        self.assertEqual(lev["label"], "BRIGHT")
        self.assertAlmostEqual(lev["tess_fraction_of_quiescent"], 1.0)


class TestIngressMap(unittest.TestCase):
    def test_windows(self):
        self.assertTrue(d.in_deep_window(100.0, 0.0))
        self.assertFalse(d.in_deep_window(200.0, 0.0))  # partial recovery
        self.assertTrue(d.in_deep_window(300.0, 0.0))
        self.assertFalse(d.in_deep_window(500.0, 0.0))
        self.assertFalse(d.in_deep_window(50.0, 0.0))

    def test_map(self):
        n = d.P_ING_HI - d.P_ING_LO + 1
        self.assertEqual(n, 241)
        t = d.T2_START + 4552.0
        rows, frac = d.ingress_map([(t, "BRIGHT")])
        exc = [r["p"] for r in rows if r["state"] == "excluded"]
        self.assertEqual(exc, list(range(4452, 4473)))
        self.assertAlmostEqual(frac["excluded"], 21.0 / 241.0)
        rows, frac = d.ingress_map([(t, "DEEP")])
        self.assertAlmostEqual(frac["consistent"], 21.0 / 241.0)
        self.assertAlmostEqual(frac["open"], 220.0 / 241.0)
        rows, frac = d.ingress_map([(t, "BRIGHT"), (t, "DEEP")])
        self.assertAlmostEqual(frac["consistent"], 0.0)
        self.assertEqual(d.runs(rows, "excluded"), [[d.T2_START + 4452, d.T2_START + 4472]])
        rows, frac = d.ingress_map([(t, "UNCERTAIN")])
        self.assertAlmostEqual(frac["open"], 1.0)

    def test_headline(self):
        self.assertEqual(d.headline("DEEP", "MID"), "DIM_DEC2025_NOT_RECOVERED_NOV2026")
        self.assertEqual(d.headline("DEEP", "BRIGHT"), "DIM_DEC2025_RECOVERED_NOV2026")
        self.assertEqual(d.headline("BRIGHT", "BRIGHT"), "BRIGHT_BOTH")
        self.assertEqual(d.headline(None, "BRIGHT"), "MIXED_MISSING_BRIGHT")

    def test_dates(self):
        self.assertEqual(d.jd_to_date(2451544.5), "2000-01-01")
        self.assertEqual(d.jd_to_date(d.T2_START), "2013-01-28")
        self.assertEqual(d.jd_to_date(d.T2_START + d.P_ING_LO), "2025-04-07")
        self.assertEqual(d.jd_to_date(d.T2_START + d.P_ING_HI), "2025-12-03")


if __name__ == "__main__":
    unittest.main()
