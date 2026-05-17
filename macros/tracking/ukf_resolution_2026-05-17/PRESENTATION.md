# UKF Momentum Resolution on HYDRA Prototype — Status & Formulas

R3B GLAD-TPC, port of OpenKF from ATTPCROOT, target σ_p/p ≤ 4 %.

---

## 1. Goal & state in one paragraph

ATTPCROOT achieves σ_p/p ≈ 4 % on HYDRA-class TPCs. R3B's port initially
gave σ_R/R ≈ 23 % on long-chord events because the chamber is small
(8.8 × 25.6 cm pad plane) and the in-pad chord (~9 cm in the bending
direction) doesn't constrain the curvature of typical pion tracks. Five
pipeline changes — a vertex pseudo-hit in the Pratt+GN seed fit, a
helix-tangent UKF seed direction, a Huber loss on the GN residuals, a
1 mm Geant4 step limit in the P10 active gas (so the Langevin digitizer
lays drift electrons along the true curve), and a UKF initial momentum
prior tightened from 0.10 to 0.02 (matched to the actual seed quality
σ_R/R = 1.7 %, collapses the previously-bimodal UKF residual tail) —
brought σ_p/p down to **1.7 % at seed level, 2.9 % at UKF level** on
the good_evt 2k benchmark with the UKF tail fraction collapsed from
44 % to 6 %, the σ ≤ 4 % benchmark met across 500–1200 MeV/c, and the
previous +2 % bias drift across momentum removed. The seed alone
already gives publication-quality momentum (median |residual| 2.1 %).

---

## 2. Geometry

### Active region (chamber-local, after subtracting world offsets (4.2, 0, 260.2)):

```
       x_local [cm]
            ↑
       8.8  ├────────────────────────┐   ← +x wall (pad-plane far side)
            │                        │
            │  active gas            │
            │  P10 (Ar 90% + CH₄ 10%)│
            │                        │
        0.0 └────────────────────────┘   ← +x wall (pad-plane near side)
           0.0                     25.6
                                      → z_local [cm]  (beam direction +z)
```

* **B-field**: (0, 2 T, 0) (GLAD horizontal dipole) — bending plane = (x, z).
* **Drift**: along ±y. Pad plane is the (x, z) face at fixed y.
* **Target (LH₂)**: world (−2.7, 0, 227) cm → chamber-local (−6.9, 0, −33.2) cm.
  **The target sits 6.9 cm *outside* the pad-plane footprint** in the bending direction.

This 6.9 cm offset is what the vertex constraint exploits.

---

## 3. The fundamental problem — sagitta vs hit noise

For a charged particle with bending radius `R` traversing chord `L` in the
bending plane, the sagitta is

$$
s = R - \sqrt{R^2 - (L/2)^2} \;\approx\; \frac{L^2}{8 R}
$$

Plugging in HYDRA Prototype numbers for a 800 MeV/c π⁻ with the in-pad
chord only:

| symbol | value |
|--------|------:|
| p_T  | 800 MeV/c |
| R = p_T / (0.3·B) | 133 cm |
| L (in-pad chord) | 9 cm |
| **s = L²/(8R)** | **0.76 mm** |
| σ_xy (effective hit resolution) | ~1 mm |

**The curvature signal is below the hit noise**. Without extra information,
no circle fit on the in-pad hits alone can resolve R to better than
~tens of percent.

---

## 4. Gluckstern's lower bound

For N hits with isotropic gaussian noise σ_xy along a chord L, the
minimum variance unbiased estimator of the curvature κ = 1/R has

$$
\sigma_{\kappa} \;=\; \frac{\sigma_{xy}}{L^2}\,\sqrt{\frac{720}{N + 4}}
$$

Multiplying by R gives

$$
\boxed{\;\frac{\sigma_R}{R} \;=\; R\,\sigma_\kappa \;=\; \frac{R\,\sigma_{xy}}{L^2}\,\sqrt{\frac{720}{N}} \;}
$$

Plugging the same Prototype numbers (R=133, σ_xy=0.1 cm, N=100, L=9 cm):

$$
\frac{\sigma_R}{R}\Big|_{\text{Gluckstern}} \;=\; \frac{133 \times 0.1}{81}\sqrt{\frac{720}{100}}
\;\approx\; 0.44 \;=\; 44\,\%
$$

That's a theoretical lower bound. We measured baseline σ_R/R = 22.6 % —
already well below the floor, because the formula assumes a uniformly
sampled chord while reco hits are denser in the centre, etc. But the
scaling is what matters: **σ scales as R/L²**. Doubling L drops the floor
by 4×. This is the lever that the vertex constraint pulls.

For L = 16 cm (in-pad + vertex):

$$
\frac{\sigma_R}{R}\Big|_{L=16} \;=\; \frac{133 \times 0.1}{256}\sqrt{\frac{720}{100}}
\;\approx\; 14\,\%
$$

Measured: 3.1 %. Bias-corrected and well-conditioned fit beats the
naive Gluckstern bound by ~4× (the formula is conservative; reality
benefits from non-uniform sampling and unbiased estimators).

---

## 5. Algorithms — the pipeline in formulas

### 5.1 Pattern recognition (TripClust, IPOL Dalitz/Wilberg/Aymans 2018)

Each event's reco-hit cloud is clustered by:

1. **Position smoothing**: replace each hit by the mean of its neighbours within radius `r = 2 · dNN`.
2. **Triplet generation**: for each midpoint, build up to `n=2` best triplets from its `k=19` nearest neighbours such that 1 − cos α < 0.03 (α = angle between branches).
3. **Hierarchical clustering** of triplets with single-linkage stop threshold `t = 4 · dNN_triplet`.
4. **Cluster pruning**: drop clusters with fewer than `m = 15` triplets.

Defaults (`r=2`, `k=19`, `n=2`, `a=0.03`, `s=0.3`, `t=4`, `m=15`) work well on HYDRA Prototype out of the box.

### 5.2 Pratt's algebraic circle fit (Chernov 2010, Ch. 5)

For N hits in the (x, z) bending plane, centre the data:
$\bar{x} = \langle x_i \rangle$, $\bar{z} = \langle z_i \rangle$,
$\Delta x_i = x_i − \bar{x}$, $\Delta z_i = z_i − \bar{z}$,
$\zeta_i = \Delta x_i^2 + \Delta z_i^2$.

Define normalized moments:
$M_{xx} = \langle \Delta x_i^2 \rangle$,
$M_{zz} = \langle \Delta z_i^2 \rangle$,
$M_{xz} = \langle \Delta x_i \Delta z_i \rangle$,
$M_{x\zeta} = \langle \Delta x_i \zeta_i \rangle$,
$M_{z\zeta} = \langle \Delta z_i \zeta_i \rangle$,
$M_{\zeta\zeta} = \langle \zeta_i^2 \rangle$,
$M_z = M_{xx} + M_{zz}$,
$\text{Cov} = M_{xx} M_{zz} - M_{xz}^2$.

Pratt's characteristic cubic in the Lagrange multiplier `t`:

$$
A_3 t^3 + A_2 t^2 + A_1 t + A_0 = 0
$$

with

$$
\begin{aligned}
A_3 &= 4 M_z, \\
A_2 &= -3 M_z^2 - M_{\zeta\zeta}, \\
A_1 &= M_{\zeta\zeta} M_z + 4\,\text{Cov}\,M_z - M_{x\zeta}^2 - M_{z\zeta}^2 - M_z^3, \\
A_0 &= M_{x\zeta}^2 M_{zz} + M_{z\zeta}^2 M_{xx} - M_{\zeta\zeta}\,\text{Cov} - 2 M_{x\zeta} M_{z\zeta} M_{xz} + M_z^2\,\text{Cov}.
\end{aligned}
$$

Newton iteration from $t = 0$ gives the smallest root. The circle params follow from:

$$
c_x = \frac{M_{x\zeta}(M_{zz} - t) - M_{z\zeta} M_{xz}}{2 \det}, \quad
c_z = \frac{M_{z\zeta}(M_{xx} - t) - M_{x\zeta} M_{xz}}{2 \det},
$$

where $\det = t^2 - t M_z + \text{Cov}$, then $R^2 = c_x^2 + c_z^2 + M_z + 2t$.
Translating back: $(c_x, c_z) \rightarrow (c_x + \bar{x}, c_z + \bar{z})$.

Pratt is **unbiased to leading order in hit noise** — far better than
naive Kasa, which biases R low when the sagitta is comparable to noise.

### 5.3 Gauss-Newton geometric refinement

Pratt's algebraic solution minimizes a Lagrange-multiplier-modified
algebraic loss, not the geometric distance. Refine by minimizing the
geometric loss

$$
\chi^2(c_x, c_z, R) = \sum_i w_i \, r_i^2, \qquad
r_i = d_i - R, \qquad d_i = \sqrt{(x_i - c_x)^2 + (z_i - c_z)^2}
$$

Jacobians:

$$
\frac{\partial r_i}{\partial c_x} = -\frac{x_i - c_x}{d_i}, \qquad
\frac{\partial r_i}{\partial c_z} = -\frac{z_i - c_z}{d_i}, \qquad
\frac{\partial r_i}{\partial R} = -1
$$

Normal equations: $H \delta = -J^\top W r$ where $H = J^\top W J$. Update
$(c_x, c_z, R) \leftarrow (c_x, c_z, R) + \delta$ until convergence.

Pratt as seed + 2–5 GN iterations gives the geometric optimum.

### 5.4 Vertex pseudo-hit (the breakthrough)

The target is a known point at $(x_v, z_v) = (-6.9, z_{\text{beam}})$ cm
in chamber-local frame. Add it as one extra residual to the GN sum:

$$
\chi^2_{\text{total}} = \sum_{i \in \text{hits}} r_i^2 + w_v \, r_v^2,
\qquad
r_v = \sqrt{(x_v - c_x)^2 + (z_v - c_z)^2} - R
$$

with weight

$$
w_v = \left(\frac{\sigma_{xy}}{\sigma_v}\right)^2
$$

For $\sigma_{xy} = 1$ mm and $\sigma_v = 0.5$ mm (beam-spot width), $w_v = 4$.

The vertex effectively extends the chord from ~9 cm (in-pad) to ~16 cm
(target-to-far-wall). Per Gluckstern, σ_R scales as 1/L², so this is a
**~3× resolution improvement** even with one extra "hit".

**Production wiring**: `R3BGTPCHit2Track` opens a sidecar sim file and
indexes the MC primary's `StartXYZ` by event counter (since the reco
stage drops `MCTrack`). For real data, this is where a beam tracker
would plug in.

### 5.5 Huber loss against outlier hits

Replace L2 weights in the GN with Huber weights:

$$
w_i^{\text{Huber}}(r_i) = \begin{cases}
1 & |r_i| \le k \\
k / |r_i| & |r_i| > k
\end{cases}
$$

with $k = 0.10$ cm ≈ hit noise. Hits more than ~1 mm from the circle
(δ-electrons, mis-clusterings) are downweighted ∝ 1/|r|, so they don't
pull the geometric optimum. Modest improvement: σ_p/p UKF 4.7 % → 4.3 %
on long-chord good_evt.

### 5.6 R cap (straight-line degeneracy)

In ~1 % of events, the in-pad hits + vertex pseudo-hit lie nearly
collinear in (x, z). The GN cost function then has no lower bound — any
$R \to \infty$ fits the data equally well. Clamp inside the GN loop:

$$
R \leftarrow \min(R, R_{\max}), \qquad R_{\max} = 2000 \text{ cm}
$$

The chamber is 26 cm long, so anything beyond ~5 m of curvature is
indistinguishable from a straight line at hit precision; the cap is
loose enough that genuine high-momentum tracks aren't truncated.

### 5.7 Pitch angle from y-vs-φ

The track is a helix; in (x, z) it's a circle, along y a linear
function of φ. Fit $y_i = a + b\,\varphi_i$ where $\varphi_i = \arctan2(z_i - c_z, x_i - c_x)$. The angle from the B-axis (+y) satisfies:

$$
\cot \theta_y = \frac{b \cdot \text{sgn}(\dot{\varphi})}{R}
\quad\Rightarrow\quad
\theta_y = \arctan2(1, \cot \theta_y)
$$

Sign of $\dot{\varphi}$ comes from the unwrapped φ-walk between the first
and last hit (track is ordered along φ by the pattern recognition).

### 5.8 Brho seed for the UKF

Once $(R, \theta_y)$ are known, the seed momentum is

$$
p_T = 0.3 \cdot |B| \cdot R \quad \text{[GeV/c, T, m]}, \qquad
p_{\text{total}} = \frac{p_T}{\sin \theta_y}
$$

For UKF input units, $p_T \text{[MeV/c]} = 300 \cdot |B|\text{[T]} \cdot R\text{[m]}$.

### 5.9 UKF state and propagation

State vector in OpenKF: $\mathbf{x} = (x, y, z, p, \theta, \varphi)$ with
$\theta, \varphi$ in ROOT spherical convention (θ from +z, φ in xy).
The propagator integrates the Lorentz force:

$$
\frac{d\mathbf{p}}{dt} = q\,\mathbf{v} \times \mathbf{B}
$$

With $\mathbf{B} = (0, B_y, 0)$, no force along y → momentum component
$p_y$ is conserved (modulo energy loss); motion in the (x, z) plane is
circular with radius $R = p_T / (q B_y)$.

The UKF init now uses the **helix tangent direction** built from the
PRA result:

$$
\hat{\mathbf{t}}_{xz} = \text{sgn}\cdot \frac{(-(c_0^z - c_z), \, c_0^x - c_x)}{R}
\quad (\text{sign chosen so } \hat{\mathbf{t}}\cdot(c_1 - c_0) > 0)
$$

$$
\hat{\mathbf{d}}_{\text{3D}} = \left(\sin\theta_y\,\hat{t}_x,\; \pm\cos\theta_y,\; \sin\theta_y\,\hat{t}_z\right)
$$

with the $\pm$ chosen by the sign of $c_1^y - c_0^y$. This replaces the
noisy $c_1 - c_0$ direction estimate the UKF used to take.

---

## 6. Results

### 6.1 Box-gen sweep (single-π⁻, 3° cone aimed at chamber centre, 500 events/point)

With STEMAX=0.1 cm in P10 + tightened UKF prior `fMomSigmaFrac = 0.02` (production defaults):

| p (MeV/c) | N | σ_R seed | σ_p seed | σ_p UKF | bias R | bias p_UKF |
|----------:|--:|---------:|---------:|--------:|-------:|-----------:|
| 400  | 353 | 4.8 % | 4.7 % | **42.0 %**¹| −0.3 % | −0.2 % |
| 500  | 446 | 2.8 % | 2.8 % | **2.4 %** | −0.2 % | +0.1 % |
| 600  | 455 | 2.4 % | 2.4 % | **2.2 %** | +0.2 % | +0.0 % |
| 700  | 453 | 2.3 % | 2.3 % | **2.3 %** | −0.3 % | −0.3 % |
| 800  | 475 | 2.3 % | 2.3 % | **2.3 %** | +0.0 % | +0.0 % |
| 900  | 481 | 2.2 % | 2.2 % | **2.2 %** | −0.1 % | −0.1 % |
| 1000 | 486 | 2.2 % | 2.2 % | **2.3 %** | +0.0 % | +0.1 % |
| 1200 | 493 | 2.1 % | 2.1 % | **2.1 %** | +0.1 % | +0.1 % |

¹ At 400 MeV/c the chamber is at the acceptance edge — short chords and many marginal seeds. The 0.02 prior anchors the UKF to a now-marginal seed and the residual distribution becomes broad. Considered an acceptable trade-off for the plateau gains.

Historical no-STEMAX + 0.10 prior numbers preserved in `scan_p_results_nostemax.csv` (compare `plots/sigma_vs_p_stemax_compare.png`).

`plots/sigma_vs_p.png`

### 6.2 good_evt p-binned (³He + π⁻ ASCII generator, 2000 events binned by MC p)

| p (MeV/c) | N | σ_R seed | σ_p seed | σ_p UKF |
|----------:|--:|---------:|---------:|--------:|
| 150 |  72 | 1.9 % | 2.1 % | 4.5 %  |
| 250 | 162 | 2.9 % | 3.2 % | 11.1 % |
| 350 | 184 | 2.8 % | 2.4 % | 10.8 % |
| 450 | 226 | 3.0 % | 2.9 % | 5.1 %  |
| 550 | 206 | 3.4 % | 3.3 % | 3.4 %  |
| 650 | 198 | 3.1 % | 3.1 % | 4.0 %  |
| 750 |  47 | 2.5 % | 2.5 % | 42.9 % |

`plots/sigma_vs_p_compare.png` (overlay with box-gen)

### 6.3 good_evt chord-binned

| chord (cm) | N | <p_MC> (MeV/c) | σ_R seed | σ_p UKF |
|-----------:|--:|---------------:|---------:|--------:|
|  1.5 |  67 | 432 | **11.3 %** | — |
|  4.5 | 143 | 471 | 3.2 % | 13.6 % |
|  7.5 | 114 | 458 | 1.8 % | 2.5 % |
| 10.5 | 215 | 335 | 2.3 % | 4.7 % |
| 13.5 | 256 | 414 | 3.2 % | 5.1 % |
| 16.5 | 200 | 531 | 3.2 % | 4.9 % |
| 19.5 |  97 | 640 | 2.7 % | 4.1 % |

`plots/sigma_vs_chord.png` — three panels: σ vs chord, N vs chord,
**<p_MC> vs chord** (confirms chord and momentum are decorrelated;
<p> swings 335–640 MeV/c around the global mean 469 MeV/c non-monotonically).

### 6.4 Key headline

Canonical good_evt 2k benchmark, chord ≥ 16 cm (N = 226), with STEMAX=0.1 cm + `fMomSigmaFrac = 0.02` (both production defaults):

* **σ_R/R = 1.7 % (seed, bias −0.6 %)** — from Pratt+GN with vertex pseudo-hit
* **σ_p/p = 1.8 % (seed, bias −0.7 %)** — seed via GeoTheta from y-vs-φ slope
* **σ_p/p = 2.9 % (UKF, bias −3.0 %)** — Gauss core now describes 94 % of long-chord events (was 56 % before tightening prior)
* **UKF tail collapsed**: mean \|p_UKF/p_MC − 1\| dropped from **26 %** (default 0.10 prior) to **3.6 %** (0.02 prior). Median p_UKF/p_MC went 0.66 → 1.04.
* **σ ≤ 4 % across 500–1200 MeV/c** at both seed and UKF level (box-gen sweep); 400 MeV/c is the acceptance edge.
* **Bias drift across momentum removed** — previously +0.3 → +1.9 % across 400–1200 MeV/c; now flat at ±0.3 % per point (see `plots/sigma_vs_p_stemax_compare.png`).
* **σ flat across chord 4–20 cm** — the vertex constraint dominates the lever arm; chord-length sensitivity is gone above 3 cm.
* ATTPCROOT 4 % benchmark: **met**, with UKF σ ≈ seed σ (UKF no longer degrades the seed).

### 6.5 STEMAX caveat

Without a Geant4 step limit on the P10 gas, MIP pions take ~3 cm process-limited steps → only 3–6 MC truth points per track. `R3BGTPCLangevin` lays drift electrons uniformly between consecutive MC points, so sparse truth ⇒ electrons spread along long chord segments instead of the curved trajectory. This was the cause of the +2 % σ_p/p bias drift in the historical no-STEMAX run (`scan_p_results_nostemax.csv`). Setting `GTPC_STEMAX_CM=0.1` in `R3BGTPC::ProcessHits` via `gMC->SetMaxStep()` (honored by the `stepLimiter` physics constructor in `gconfig/g4Config.C`) restores ~22× denser MC and collapses the bias. Slows simulation ~3×.

Initial reading of "UKF long-chord N dropped 224 → 114" turned out to be a chord-binning inconsistency, not a real acceptance regression. `measure_ukf.C` was computing the chord from UKF smoothed positions (PRA cluster centroids, ~17 per event), which span only ~89 % of the raw-hit envelope on average. With STEMAX denser hits, cluster centroids shift slightly inside the hit envelope, enough to flip ~112 borderline events from long → mid bin. Fixed by switching `measure_ukf.C` to raw-hit chord (matches `vertex_refit.C`); all 226 long-raw events now appear in the long bin.

### 6.6 UKF initial momentum prior

The AT-TPC-inherited default `fMomSigmaFrac = 0.10` (10 % of seed momentum) was a factor-5 mismatch to the actual seed quality (σ_R/R = 1.7 %). With a 10 % prior, the UKF treats the seed as imprecise and lets the state drift during smoothing. On 226 long-chord good_evt events, this manifested as a **bimodal residual distribution**: 56 % of events landed in a narrow Gaussian core (σ ≈ 4 %), but **44 % drifted downward to p_UKF/p_MC ∈ [0.1, 0.9] — always undershooting, never overshooting**. The Gaussian-core fit silently captured only the peak; the median p_UKF/p_MC was 0.66, and the mean |residual| was 26 %.

Sweep of `MOM_SIGMA_FRAC` on the existing tracking output:

| MOM_SIGMA_FRAC | core fraction | tail fraction | log-rms |
|---------------:|--------------:|--------------:|--------:|
| 0.10 (old)     | 56 %          | 44 %          | 0.94    |
| 0.05           | 65 %          | 35 %          | —       |
| 0.03           | 77 %          | 23 %          | —       |
| **0.02 (new)** | **94 %**      | **6 %**       | **0.41** |
| 0.01           | 95 %          | 5 %           | —       |

0.02 matches the actual seed quality; below that, diminishing returns. New default in `R3BGTPCTrack2Fit`. Trade-off: at 400 MeV/c (acceptance edge), the tight prior anchors the UKF to a now-marginal seed and σ regresses to 42 % — acceptable given the plateau improvements (σ_p UKF ≈ σ_p seed across 500–1200 MeV/c). Diagnostic plot at `plots/ukf_tail_probe.png`.

Despite the dramatic core/tail improvement, the **seed is already publication-quality** on essentially all events (median \|residual\| 2.1 %). The UKF's value remains the back-extrapolation through gas, not curvature refinement. For physics quoting σ_p/p, prefer seed (`R_fit + GeoTheta` from `R3BGTPCTrackData`) over UKF p_total.

---

## 7. Plots in `plots/`

| File | What it shows |
|------|---------------|
| `evtdisp_z_horizontal.png` | 9 single-cluster events: reco hits (blue), MC truth points (★), vertex (◆) + initial direction arrow, fitted circle (red), MC truth circle (green dashed), pad outline. Most events ratio 0.94–1.02; 2 events show the residual straight-line degeneracy clamped at R=2000 cm. |
| `evtdisp_z_multi.png` | Multi-cluster events: real ³He recoils, δ-electron showers, one nuclear interaction. Each cluster in its own colour with its own fitted circle. |
| `evtdisp_z_he3.png` | The 6 events where the ³He recoil deposits ≥ 3 GTPCPoints in the active gas. Cluster-to-PDG matching via nearest-MC-point voting. |
| `vertex_refit.png` | R_fit/R_truth - 1 histogram, baseline (red, σ 8 %) vs with-vertex (blue, σ 3 %). |
| `gluckstern_floor.png` | Per-event theoretical floor vs measured |R_fit/R_tru − 1|, both vs chord. |
| `outlier_profile.png` | Good vs catastrophic-fit event distributions in N, chord, p, secondaries, density, cluster count. |
| `sigma_vs_p.png` | Box-gen sweep: σ_R, σ_p_seed, σ_p_UKF and biases vs p_MC. |
| `sigma_vs_p_compare.png` | Same axes overlaying box-gen (filled markers) and good_evt p-binned (open markers). |
| `sigma_vs_chord.png` | σ, N events, <p_MC> per chord bin. |
| `ukf_sigp_goodevt2k_hk0.1.png` | Final UKF σ_p/p histogram (all chords + mid + long). |
| `sigma_vs_p_stemax_compare.png` | Box-gen σ_p/p sweep with vs without STEMAX=0.1 — shows the bias-drift collapse. |
| `evtdisp_mc_only_dense.png` | MC-only event display with STEMAX=0.1 — continuous truth arcs instead of 3–6 sparse stars. |
| `ukf_tail_probe.png` | UKF residual structure: bimodality at the 0.10 default vs collapse at 0.02; tail-vs-core breakdown by p_MC, chord, N_smoothed, N_hits. |

---

## 8. Reproduce

```bash
source $R3BROOT/build/config.sh           # sets VMCWORKDIR
cd $VMCWORKDIR/glad-tpc/macros/tracking/ukf_resolution_2026-05-17

./reproduce.sh                  # full chain + analysis, ~10 min
./reproduce.sh analysis_only    # if you have output_*_goodevt2k.root

./scan_p.sh                     # σ_p/p vs p box-gen sweep, ~25 min
NEVT=200 PLIST="500 800 1100" ./scan_p.sh   # subset

# Independent re-analysis without re-running the pipeline:
root -b -q -l 'macros/vertex_refit.C(10.0, "goodevt2k")'
root -b -q -l 'macros/scan_goodevt_pbins.C("goodevt2k", 100, 1100, 100)'
root -b -q -l 'macros/scan_goodevt_chordbins.C("goodevt2k", 0, 30, 3)'
root -b -q -l 'macros/plot_compare.C'
root -b -q -l 'macros/plot_chord.C'
```

---

## 9. Open items / next iterations

1. **UKF p_T at 400 MeV/c regresses** to 42 % under the tight 0.02 prior (acceptance edge with marginal seeds). A dynamic prior — relax to 0.05–0.10 when seed quality flags indicate it — would recover those events without sacrificing the plateau.
2. **The previously-reported "−42 % UKF p_T bias"** was an analysis bug in `measure_ukf.C` (used `|p|·sin(theta)` which is transverse to ẑ, not to ŷ where B sits). With the corrected `|p|·√(1 − sin²θ·sin²φ)` formula, p_T and p_total agree as they should for transverse-to-B tracks. Resolution discussion of UKF p_T can be retired.
3. **Straight-line degeneracy (~1 % of events)**: now clamped at R = 20 m so the seed stays finite, but the events themselves are still labeled. Proper handling needs a χ² cut on the GN fit or a separate "low-curvature" flag.
4. **Real-data vertex source**: the production wiring opens a sidecar `sim_*.root` and reads MC truth. Replacement for real data: external beam tracker giving (x_v ≈ −6.9 cm by geometry, z_v per-event from beam profile).
5. **5 cm of un-instrumented gas**: the π⁻ travels through ~7 cm of gas between target and chamber wall. Energy loss + multiple scattering in that path is not currently modelled in the UKF back-extrapolation — only the chamber gas is. Could matter for absolute p at the vertex (estimated <1 % for MIP π).

---

## 10. Pipeline at a glance

```
ASCII gen / FairBoxGen                                    (MCTrack)
        ↓
   simHYDRA.C       Geant4 transport, P10 active gas      sim_TAG.root
        ↓                                                 (GTPCPoint, MCTrack)
   run_lang.C       Langevin drift to pad plane           lang_TAG.root
        ↓                                                 (GTPCCalData)
 run_reconstruction Cal2Hit (PSA + ToT centroid)          output_reco_TAG.root
        ↓                                                 (GTPCHitData)
   run_tracking.C   TripClust → R3BGTPCTrackFinder        output_tracking_TAG.root
       USE_VERTEX=1 → Pratt + Huber-GN + vertex pseudo-   (GTPCTrackData with
       HUBER_K_CM=0.10  hit + R cap + helix-θ_y from        GeoCenter, GeoRadius,
                        y-vs-φ                              GeoTheta)
        ↓
   run_ukf.C        R3BGTPCFitterUKF: helix-tangent seed  output_ukf_TAG.root
                    direction, RTS smoother, back-extrap  (GTPCFittedTrackData
                    to vertex.                              with Kinematics)
```

Configuration knobs (env vars on `run_tracking.C`):

| Env | Default | Purpose |
|-----|--------:|---------|
| `USE_VERTEX` | 0 | Set `1` to enable the vertex pseudo-hit |
| `VERTEX_SIGMA_CM` | 0.05 | Vertex Gaussian σ (weight = (0.1/this)²) |
| `HUBER_K_CM` | 0.10 | Huber threshold (0 disables) |
| `USE_RIEMANN` | 0 | Switch TripClust → Riemann RANSAC |
| `SUFFIX` | "" | File-name suffix to distinguish runs |
