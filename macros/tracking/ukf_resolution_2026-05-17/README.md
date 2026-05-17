# UKF momentum-resolution study — HYDRA Prototype

Reproduces the σ_p/p measurement on the HYDRA Prototype with the
realistic `good_evt` (³He + π⁻) generator. Final accepted result:

| Stage | σ_p/p (Gauss core, long chord ≥16 cm) | bias |
|-------|--------------------------------------:|-----:|
| baseline Pratt+GN (no vertex) | 22.6 % | +17.2 % |
| **+ vertex constraint** | **3.1 %** (σ_R/R), 3.0 % (σ_p/p seed) | +0.2 % |
| **+ UKF, helix seed direction** | **4.7 %** | +0.2 % |
| **+ Huber GN (k=0.1 cm)** | **4.3 %** | +0.4 % |

ATTPCROOT 4 % reference: effectively met at UKF level; comfortably
under at seed level.

## Pipeline

The full chain is:

```
sim   →  lang   →   reco       →  tracking            →  ukf
sim_X    lang_X     output_reco_X   output_tracking_X     output_ukf_X
       (Langevin)   (Cal2Hit)       (TripClust+Pratt+GN+  (FitterUKF)
                                     vertex+Huber)
```

Each stage's macro is in its canonical directory (`macros/sim`,
`macros/proj`, `macros/reco`, `macros/tracking`). The `reproduce.sh`
script below runs the full chain with the settings that produced the
results above.

## Reproduce

```bash
cd $VMCWORKDIR/glad-tpc/macros/tracking/ukf_resolution_2026-05-17
./reproduce.sh                  # full chain + analysis, ~10 min
./reproduce.sh analysis_only    # if you already have sim_goodevt2k.root etc.
```

This regenerates `plots/*.png` and prints σ_p/p numbers to stdout.

## Key knobs (run_tracking.C env vars)

| Env var | Default | Effect |
|---------|---------|--------|
| `USE_VERTEX` | `0` | Set to `1` to enable the vertex pseudo-hit constraint |
| `VERTEX_SIGMA_CM` | `0.05` | Vertex pseudo-hit Gaussian σ |
| `HUBER_K_CM` | `0.10` | Huber threshold (set to `0` for pure L2) |
| `USE_RIEMANN` | `0` | Use Riemann RANSAC instead of TripClust |
| `SUFFIX` | `""` | Filename suffix (`_goodevt2k` for the 2k run) |

## Macros (in `macros/`)

- `vertex_refit.C` — standalone Pratt+GN refit with vertex constraint,
  reads tracking output, prints/plots σ_R/R and σ_p/p with and without
  vertex. **Main resolution measurement.**
- `measure_ukf.C` — reads UKF output (`output_ukf_X.root`), computes
  σ_p_total/p_total and σ_p_T/p_T Gauss cores per chord bin.
- `gluckstern_floor.C` — per-event theoretical floor σ_R from
  σ_xy·sqrt(720/N)/L² for sanity-check.
- `profile_outliers.C` — splits the long-chord events into good vs
  catastrophic bins and compares predictors.
- `evtdisp_z.C` / `evtdisp_z_multi.C` / `evtdisp_z_he3.C` — event
  displays in (z, x) chamber-local frame with reco hits, MC truth,
  vertex, fitted circle, truth circle, pad-plane outline.

## Plots (in `plots/`)

- `evtdisp_z_horizontal.png` — 9 single-cluster events showing the
  full kinematic picture (vertex outside chamber → diagonal track
  inside pad plane → fitted vs MC-truth circles).
- `evtdisp_z_multi.png` — 9 multi-cluster events (real ³He + π⁻
  separations, δ-electron showers, one nuclear interaction).
- `evtdisp_z_he3.png` — the 6 events where ³He reaches the active gas.
- `vertex_refit.png` — R_fit/R_truth - 1 histogram, baseline (red)
  vs with-vertex (blue), long + mid chord.
- `gluckstern_floor.png` — theoretical σ_R/R floor vs chord length
  vs measured deviation.
- `outlier_profile.png` — feature distributions for good vs
  catastrophic fits (10 % long-chord outliers).
- `ukf_sigp_goodevt2k_hk0.1.png` — final UKF σ_p/p distribution.

## Geometry recap

HYDRA Prototype: pad plane is the (x, z) plane, [0, 8.8] × [0, 25.6]
cm in chamber-local frame. B = (0, 2 T, 0) (GLAD horizontal dipole).
Drift along ±y. Target at world (−2.7, 0, 227) cm → chamber-local
x = −6.9 cm (the vertex pseudo-hit x).

The good_evt ASCII generator produces 2 primaries per event: π⁻ and
³He. π⁻ reaches the active gas in 317/500 events; ³He only in 6/500.
Most multi-cluster events are physical (³He recoils, δ-electron
showers), not algorithm over-splits.

## Commits

```
3ba70e7 feat(reco): Huber-weighted Gauss-Newton in Pratt+GN circle fit
d6546ad diag: outlier profile macro for σ_R fit tails
81e42fa fix(reco): seed UKF direction from PRA helix (R, θ_y)
3626cb2 feat(reco): vertex-constrained Pratt+GN seed (σ_R/R 23%→3%)
f7e1b51 feat(macros): event-display diagnostics for HYDRA Prototype
e459275 feat(reco): port Riemann RANSAC + expose TripClust params
```
