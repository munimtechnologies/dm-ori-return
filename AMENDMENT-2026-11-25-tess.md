# Amendment 2026-11-25: a graded TESS arm

Fixed 2026-11-25 (UTC), before job dm-ori-tess-v1 ran. Rules: `dmtess.py` (tested in
`test_dmtess.py`), run by `job_tess.py`. Version 0.2.0-tess. The rules of PREREGISTRATION.md are not
changed; this adds a second, separately graded arm.

## Why, and what we had already seen

v0.1.0 failed for coverage (GRADING-v0.1.0.md): the public ASAS-SN light curve ends on 2025-02-06.
The only public photometry of DM Ori from the predicted window that we can reach is TESS:
sector 98 (late Nov to Dec 2025) and sector 110 (Oct-Nov 2026). **This amendment was written
after we saw** the v0.1.0 outputs. Those include the raw TESS sector ratios (sector 98 / sector 6 = 0.759,
sector 110 / sector 6 = 1.105) and the ASAS-SN light curve to 2025-02-06, which shows a ~1 mag
dip inside sector 6. We therefore fix in advance the calibration and comparison-star rules that turn
those ratios into depths. We do not adjust them after this run. The depth thresholds (0.75 and
0.30 mag) and the occultation geometry are the ones frozen on 2026-11-24.

## Rules (dmtess.py)

1. **Photometry.** 21x21-pixel TESScut cutouts of sectors 6, 98 and 110. Target: 3x3 aperture on the
   brightest of the central 3x3 pixels; per-cadence background = median of pixels >= 4 px away
   whose median value is below the 60th percentile (as in v0.1.0). Quality bits 1, 2, 4, 8, 32 and 128
   are dropped.
2. **Comparison stars.** Local maxima in the sector-6 median image, 0.2-5x the target's peak
   (background-subtracted), >= 4 px from the target and >= 2 px from the cutout edge. The brightest
   12 are kept and matched in sectors 98 and 110 by sky position (cutout WCS; brightest pixel
   within +-1 px, which must be a local maximum). Sector scale k = median over matched stars of
   (median flux in the sector / median flux in sector 6). Systematic fraction sys =
   max(0.02, 1.4826 x MAD / k). **Gate COMPS:** >= 4 matched stars and sys <= 0.10 in each target sector.
3. **Calibration.** Each ASAS-SN g point in `dmtess.ASASSN_S6` (copied from
   results/v0.1.0/asassn_points.csv, relative to quiescence) is paired with the median sector-6
   TESS target flux within +-0.25 d, if >= 3 TESS points fall there. Fit TESS flux = a + b x r_g.
   **Gate CALIBRATION:** >= 6 pairs, span of r_g >= 0.40, b > 0, Pearson correlation >= 0.90.
   The intercept a takes up flux that does not dim (neighbours in the 63-arcsec aperture). The
   slope b assumes the TESS-to-g depth ratio of the December 2018 dip.
4. **Depth per sector.** F = (target median flux) / k; ASAS-SN-equivalent depth
   dmag = -2.5 log10((F - a) / b) (capped at 3). Recomputed with k(1 - sys) and k(1 + sys):
   **DEEP** if even the shallower depth is >= 0.75 mag; **BRIGHT** if even the deeper depth is
   <= 0.30 mag; **MID** if both lie between; otherwise **UNCERTAIN**.
5. **Ingress map.** Allowed ingress starts (the 2000 bounds, as in v0.1.0): 2025-04-07 to
   2025-12-03, 1-day steps (241 dates). With the 2013-14 shape, an ingress at ti puts DM Ori in a
   deep phase during ti + 80..185 d and ti + 225..420 d (the gap is the ~30-day partial recovery
   seen in August 2013). A BRIGHT sector midpoint inside a deep phase **excludes** ti. A DEEP one
   makes ti **consistent** (unless another sector excludes it). Otherwise ti stays **open**.
6. **Headline** (from the sector-98 and sector-110 labels):
   DIM_DEC2025_NOT_RECOVERED_NOV2026 (DEEP, DEEP/MID), DIM_DEC2025_RECOVERED_NOV2026 (DEEP, BRIGHT),
   BRIGHT_BOTH, BRIGHT_DEC2025_DIM_NOV2026 (BRIGHT, DEEP/MID), otherwise MIXED_<98>_<110>.
   Failed gates give CALIBRATION_FAILED or COMPS_FAILED and no grade.

## What the outcomes would mean (fixed now)

- **BRIGHT_BOTH:** a 2013-14-like occultation is ruled out for ingress starts 2025-04-07..04-26
  and 2025-06-05..12-03 (~84% of the allowed range). Only 2025-04-27..06-04 stays open, where both
  snapshots fall in non-deep phases. Graded as a strong negative for the strict 12.5-yr period.
- **DIM_DEC2025_...:** DM Ori was >= 0.75 mag down in December 2025, inside the predicted
  occultation window. This is necessary for the prediction, but a single 27-day snapshot is not
  proof, because short dips put DM Ori in such a state ~6.3% of the time (GRADING-v0.1.0.md). The
  daily light curve in sector 98 is reported. A sector-long deep state favours the occultation; a
  dip that starts and ends inside the sector favours a short dip. That reading is descriptive,
  not graded. RECOVERED_NOV2026 is consistent with ingress starts in 2025-04-07..04-26 or
  2025-06-05..09-13.
- **BRIGHT_DEC2025_DIM_NOV2026:** consistent only with a late ingress (2025-09-19..12-03).
- **MID / UNCERTAIN** labels are reported as found.

## Limits

TESS pixels are 21 arcsec, so the aperture holds neighbours. The calibration absorbs their constant
light but assumes one dimming colour. Depths beyond the calibrated range (r_g 0.30-0.88) are
extrapolations and are flagged. Comparison stars in Orion can be variable; the median and the MAD
systematic limit this. Two snapshots cannot measure a period. This work was carried out by AI
agents working for MV2 (Munim, Inc.); no human expert has checked it yet.
