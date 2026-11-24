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
(control). Results will appear under `results/`.

- `dmori.py`: references, binning, episodes, control, outcome, bracket check, period map (stdlib, tested)
- `asassn_io.py`, `io_tess.py`: public data access (ASAS-SN Sky Patrol v2, MAST TESScut)
- `job.py`: the full run: `python job.py --data-dir /data --out-dir /out`
- `test_dmori.py`: `python -m unittest -v`

Credit: the events, the disk-warp interpretation and the July 2026 prediction are from J. E.
Rodriguez, K. G. Stassun, P. Cargile, B. J. Shappee et al. (2016, ApJ 831, 74). Data: ASAS-SN Sky
Patrol v2, ZTF (IRSA), TESS (NASA, MAST).

By MV2, Munim, Inc. (mv2@munimtech.com). Made with AI agents working in his name; no human expert
has checked it yet. Code Apache-2.0; text and results CC BY 4.0.
