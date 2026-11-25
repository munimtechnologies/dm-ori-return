# dm-ori-return

Did the young star DM Ori start the deep, year-long occultation that was predicted for July 2026?

DM Ori, a young star near the Orion star-forming region, faded by about 1.5 magnitudes for more
than a year in 2000-2002 and again in 2013-2014. Rodriguez et al. (2016,
[arXiv:1607.07455](https://arxiv.org/abs/1607.07455)) proposed that a warp or other disturbance in
its planet-forming disk, possibly raised by a young protoplanet, passes in front of the star every
~12.5 years, and predicted the next dimming for July 2026. That date has passed and public survey
photometry covers it, but we found no published test.

This repository holds a pre-registered test with ASAS-SN Sky Patrol (primary), ZTF and TESS
(secondary). The rules are in [PREREGISTRATION.md](PREREGISTRATION.md), fixed on 2026-11-24 before
any post-2016 photometry was retrieved. The same pipeline must first recover the known 2013-14 event
(control).

## Status (2026-11-25)

- **v0.1.0: CONTROL_FAILED_INCONCLUSIVE_COVERAGE**, not a result. The public ASAS-SN Sky Patrol v2
  light curve ends on 2025-02-06, before the predicted window, and its 2013-14 sampling is too thin
  for the control. See [GRADING-v0.1.0.md](GRADING-v0.1.0.md) and `results/v0.1.0/`. It will be
  re-run unchanged when Sky Patrol publishes 2025-2026 photometry.
- **TESS arm** ([AMENDMENT-2026-11-25-tess.md](AMENDMENT-2026-11-25-tess.md)): TESS sectors 98
  (Dec 2025) and 110 (Nov 2026) are read as depths. The calibration uses the ASAS-SN dip that TESS
  sector 6 caught in December 2018, and comparison stars set each sector's scale. Job dm-ori-tess-v1
  is queued.

## Files

- `dmori.py`: references, binning, episodes, control, outcome, bracket check, period map (stdlib, tested)
- `asassn_io.py`, `io_tess.py`: public data access (ASAS-SN Sky Patrol v2, MAST TESScut)
- `job.py`: the v0.1.0 run: `python job.py --data-dir /data --out-dir /out`
- `test_dmori.py`: `python -m unittest -v test_dmori`
- `dmtess.py`, `test_dmtess.py`, `job_tess.py`, `requirements-tess.txt`: the TESS arm
  (`python job_tess.py --data-dir /data --out-dir /out`)

Credit: the events, the disk-warp interpretation and the July 2026 prediction are from J. E.
Rodriguez, K. G. Stassun, P. Cargile, B. J. Shappee et al. (2016, ApJ 831, 74). Data: ASAS-SN Sky
Patrol v2, ZTF (IRSA), TESS (NASA, MAST).

By MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working in his name; no human expert
has checked it yet. Code Apache-2.0; text and results CC BY 4.0.
