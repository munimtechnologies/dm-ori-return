# Grading of v0.1.0 (job dm-ori-return-v1, run 2026-11-25)

Graded 2026-11-25 by AI agents working in MV2's name, under the rules in PREREGISTRATION.md
(fixed 2026-11-24). Outputs: `results/v0.1.0/`.

## Outcome: CONTROL_FAILED_INCONCLUSIVE_COVERAGE (not graded as a result)

Under section 4 of the pre-registration this run answers nothing about 2026. Both halves failed
for lack of data, not because of the code:

1. **No primary data in the test window.** The ASAS-SN Sky Patrol v2 light curve of DM Ori
   (asas_sn_id 249109296903, 0.16 arcsec from the target, 1338 rows) ends on **2025-02-06**.
   The 2025-2026 season is not in the public v2 database yet. The client printed a notice dated
   3 November 2026 that Sky Patrol is "rebuilding our catalogs and lightcurves in full". So
   `main` = INCONCLUSIVE_COVERAGE, the bracket core (2026-02-21..2026-06-01) has 0 bins, and the
   period map excludes only 1.0% of 11.5-13.5 yr.
2. **The control could not be recovered at this cadence.** Between January 2013 and November
   2014 the V-band data give only five 30-day bins with at least 3 good points. One is the
   pre-event baseline (2013-01-08, -0.29 mag). The others are ingress (2013-02-07, +0.29 mag),
   2013-08-06 (+0.21 mag, which falls in the ~30-day partial recovery from JD 2456511 that
   Rodriguez et al. reported) and **2013-09-05 (+1.22 mag, DIM)**. The next bin is 2014-11-29.
   The rule needs two DIM bins in one episode, so the control fails. The one bin from the deep
   phase is 1.22 mag down, so the pipeline sees the event; the 2013-14 sampling is too thin.

## Reported, not graded

- **Short dips between the two long events.** In the 95 bins from 2014-11-29 to 2025-01-05,
  6 are DIM (>= 0.75 mag down: 2015-10-25, 2016-03-23, 2018-10-09, 2019-01-07, 2019-09-04,
  2024-08-08) and 12 are MID. Each episode lasts one to three bins. So at a random time
  DM Ori is in a deep (>= 0.75 mag) state about **6.3%** of the time. That is the base rate any single
  snapshot has to beat.
- **Last data.** The g-band flux fell from 0.91 of quiescence (2025-01-06) to 0.67-0.69
  (2025-02-05/06, 0.41-0.44 mag down) when the light curve stops. The earliest ingress the 2000
  bounds allow is 2025-04-07. This fade began ~2.5 months before that date and may be one more
  short dip.
- **TESS (secondary).** Only sectors 6 (2018-12-27), 98 (2025-12-07) and 110 (2026-11-08) cover
  DM Ori. Median 3x3-aperture flux relative to sector 6: **0.759** (sector 98) and **1.105**
  (sector 110); the peak-pixel ratios agree (0.760, 1.234). But sector 6 is not a quiescent reference:
  ASAS-SN g shows DM Ori falling from 0.87 to 0.33 of quiescent flux inside sector 6
  (2018-12-18 to 2019-01-08). These raw ratios can't be read as depths. That is why the TESS
  arm in AMENDMENT-2026-11-25-tess.md exists.
- **ZTF.** 0 rows within 2.2 arcsec. DM Ori (V ~ 11.7) is too bright for the ZTF catalogue.

## Next

- The TESS arm (AMENDMENT-2026-11-25-tess.md, job dm-ori-tess-v1) calibrates TESS against the
  ASAS-SN dip inside sector 6 and uses comparison stars, so sectors 98 and 110 can be read as depths.
- When ASAS-SN Sky Patrol publishes 2025-2026 photometry, v0.1.0 is re-run unchanged as
  dm-ori-return-v2.

This grading was done by AI agents working for MV2 (Munim, Inc.); no human expert has checked it yet.
