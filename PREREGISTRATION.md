# Pre-registration: did DM Ori's predicted 2026 occultation happen?

Fixed 2026-11-24 (UTC), before any photometry of DM Ori after 2016 was retrieved or looked at by us.
Code: `dmori.py` (rules), `job.py` (run), `test_dmori.py` (tests). Version 0.1.0.

## 1. Setup (credit)

Rodriguez, Stassun, Cargile, Shappee, Siverd, Pepper, Lund, Kochanek et al. (2016), "DM Ori: A Young
Star Occulted by a Disturbance in its Protoplanetary Disk", ApJ 831, 74,
[arXiv:1607.07455](https://arxiv.org/abs/1607.07455), found two long, deep dimmings of the young
star DM Ori (Tycho-2 J2000 05:43:39.0 -05:04:02.83; quiescent V ~ 11.7):

- 2000-2002: ~1.5 mag; last bright point 2000-03-25 (NSVS), first dim point 2000-11-20 (ASAS),
  first bright point after the seasonal gap 2002-08-02; duration > 396 d and < 860 d.
- 2013-2014: ~1.5-1.7 mag; ingress began near the end of January 2013 and took ~80 d; a ~30-day
  partial recovery from JD 2456511; egress in the 2014 seasonal gap; duration > 420 d and < 585 d.

They proposed a periodic occultation by a warp or other disturbance more than 6 au out in the
disk, possibly caused by a forming protoplanet, with event centres ~12.5 yr apart. They warned that two
events cannot prove strict periodicity, and called for monitoring "especially during the next
predicted dimming in July of 2026". We found no published test (scoop check in
`notes/scoop-check-2026-11-24.md`).

## 2. Question

Between mid-2025 and the end of the available data (expected November 2026), was DM Ori in a deep
dimming of the kind seen in 2000-2002 and 2013-2014?

## 3. Data

- **Primary (only source used in the verdict):** ASAS-SN Sky Patrol v2 (pyasassn client), nearest
  catalogued source within 10 arcsec, quality "G" only, flux/flux_err >= 5. V band (2012-2018) and
  g band (2017 on).
- **Secondary (reported, not in the verdict):** ZTF public light curves (IRSA API). DM Ori is near
  or above ZTF saturation when bright, so unflagged ZTF points are biased towards faint states.
  TESS full-frame-image cutouts (MAST TESScut): median 3x3-aperture flux per sector relative to the
  sectors before 2025-07-01. Both are corroboration only.

## 4. Frozen rules (dmori.py)

- **Quiescent reference** per (band, camera): median flux in V 2014-12-01..2019-01-01 and
  g 2017-12-01..2025-01-01. A camera with fewer than 20 points in its band's window uses the band
  median. (This is deliberately not a whole-light-curve median: a camera that saw the star only
  while dim must not hide a dimming.)
- **Bins:** 30 days, origin JD 2450000.5, at least 3 points; bin value = median of flux/reference.
  DIM if >= 0.75 mag below quiescence (half the published depth), BRIGHT if within 0.30 mag,
  otherwise MID.
- **Episode:** a maximal run of non-BRIGHT bins with gaps of at most 150 days without data.
  An episode qualifies if it holds >= 2 DIM bins.
- **Control (known answer):** a qualifying episode must overlap 2013-01-28..2014-09-01. If the
  control fails, the outcome is prefixed CONTROL_FAILED_ and is not graded as a result.
- **Main outcome:**
  - RETURNED_ON_SCHEDULE: a qualifying episode ending after 2025-07-01 overlaps calendar 2026.
  - RETURNED_EARLY: such an episode exists but ends before 2026-01-01.
  - DIMMED_BEFORE_2025_07: the only qualifying episode after 2024-07-01 ended before 2025-07-01.
  - SHALLOW_ONLY: no qualifying episode, but a non-qualifying episode of >= 3 bins after 2025-07-01.
  - NOT_RETURNED: none of the above, with >= 1 BRIGHT bin in 2026-02-01..2026-05-01, >= 1 BRIGHT
    bin from 2026-08-01, and >= 4 BRIGHT bins from 2025-10-01.
  - INCONCLUSIVE_COVERAGE: otherwise.
- **Bracket core (sharper, period-free check):** ingress-to-ingress separations allowed by the
  2000 bounds are 4452-4692 d (12.19-12.85 yr). With the 2013 shape (ingress ~80 d, duration > 420 d)
  every one of them puts DM Ori deep in the dimming during 2026-02-21..2026-06-01. Label
  BRACKET_DIM if >= 50% of the bins there (at least 1) are DIM; BRACKET_BRIGHT if >= 2 bins there
  and all are BRIGHT; else BRACKET_UNDETERMINED.
- **Period map:** for each period P = 11.00..14.00 yr (0.01 yr steps), the predicted dim interval
  is the 2013-14 event centre (JD 2456611.0) + P, +-210 d (half the minimum duration). P is
  excluded if the interval holds >= 3 bins, >= 3 of them BRIGHT and < 25% DIM; supported if >= 3
  bins and >= 50% DIM; otherwise open. Reported: excluded and supported fractions of 11.5-13.5 yr.

## 5. What the outcomes would mean

- RETURNED_ON_SCHEDULE plus BRACKET_DIM: the third occultation came when predicted, which strongly
  supports a periodic structure in the disk; the onset date pins the period. Grade as positive.
- NOT_RETURNED plus BRACKET_BRIGHT: a strictly periodic occulter with the 2013-14 shape and a period
  of 12.2-12.85 yr is ruled out. A longer period, a changing structure or a non-periodic cause stays
  possible; the period map says which periods remain. Grade as negative.
- Anything else is reported as found, with its label.

## 6. Deviations

Any change after this date will be logged here with its date and reason before the result is
graded. A failed run caused by infrastructure (archive outage, network) is re-run unchanged as -v2.

## 7. Limits

ASAS-SN pixels are 8 arcsec, so nearby stars can dilute a dimming; the 2013-14 control checks that
the same pipeline sees a known event. Seasonal gaps (Orion is behind the Sun around May-July) hide
some months, which is why the verdict uses bins on both sides of the gap. This work was carried out
by AI agents working for MV2 (Munim, Inc.); no human expert has checked it yet.
