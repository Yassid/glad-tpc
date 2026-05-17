#!/usr/bin/env bash
# Reproduce the σ_p/p study on HYDRA Prototype, 2000 good_evt events.
#
# Usage:
#   ./reproduce.sh                  full chain + analysis
#   ./reproduce.sh analysis_only    skip sim/lang/reco/tracking/ukf,
#                                   only re-run the analysis macros
#                                   (requires output_*_goodevt2k.root)
#
# Requires VMCWORKDIR set (source build/config.sh in R3BRoot first).

set -uo pipefail   # NOT -e: ROOT exits non-zero on shutdown heap warnings
                   # even when the FairRunAna output is written correctly.
                   # We check file existence after each stage instead.
MODE="${1:-full}"
[ -z "${VMCWORKDIR:-}" ] && { echo "ERROR: VMCWORKDIR not set (source build/config.sh)"; exit 1; }

require_file() {
    if [ ! -s "$1" ]; then
        echo "ERROR: expected output $1 not produced (or empty)"
        exit 1
    fi
}

TAG=goodevt2k
NEVT=2000
HERE=$(cd "$(dirname "$0")" && pwd)
TRACKDIR="$VMCWORKDIR/glad-tpc/macros/tracking"
PLOTS="$HERE/plots"
mkdir -p "$PLOTS"

if [ "$MODE" != "analysis_only" ]; then
  echo "===== [1/5] sim ($NEVT events, good_evt) ====="
  cd "$VMCWORKDIR/glad-tpc/macros/sim"
  SUFFIX="_${TAG}" root -b -q -l "simHYDRA.C($NEVT, \"Prototype\", \"good_evt\")" || true
  require_file "Prototype/sim_${TAG}.root"

  echo "===== [2/5] Langevin drift ====="
  cd "$VMCWORKDIR/glad-tpc/macros/proj"
  SUFFIX="_${TAG}" root -b -q -l run_lang.C || true
  require_file "Prototype/lang_${TAG}.root"

  echo "===== [3/5] Cal2Hit reco ====="
  cd "$VMCWORKDIR/glad-tpc/macros/reco"
  SUFFIX="_${TAG}" root -b -q -l "run_reconstruction.C(\"lang_${TAG}.root\")" || true
  require_file "output_reco_${TAG}.root"

  echo "===== [4/5] tracking (TripClust + Pratt+GN + vertex + Huber) ====="
  cd "$TRACKDIR"
  USE_VERTEX=1 HUBER_K_CM=0.10 SUFFIX="_${TAG}" \
    root -b -q -l "run_tracking.C(\"output_reco_${TAG}.root\")" || true
  require_file "output_tracking_${TAG}.root"

  echo "===== [5/5] UKF fit (helix-direction seed) ====="
  SUFFIX="_${TAG}" root -b -q -l run_ukf.C || true
  require_file "output_ukf_${TAG}.root"
fi

echo "===== analysis: σ_p/p numbers ====="
cd "$TRACKDIR"

echo "--- seed-level (vertex_refit.C, vertex weight = 10) ---"
root -b -q -l "$HERE/macros/vertex_refit.C(10.0, \"${TAG}\")"
mv vertex_refit.png "$PLOTS/" 2>/dev/null || true

echo "--- UKF-level (measure_ukf.C) ---"
root -b -q -l "$HERE/macros/measure_ukf.C(\"output_ukf_${TAG}.root\", \"Prototype/sim_${TAG}.root\", \"${TAG}\", -211, true)"
mv "ukf_sigp_${TAG}.png" "$PLOTS/ukf_sigp_goodevt2k_hk0.1.png" 2>/dev/null || true

echo "--- Gluckstern floor (σ_xy = 1 mm) ---"
ln -sf "$VMCWORKDIR/glad-tpc/macros/sim/Prototype/sim_${TAG}.root" sim_goodevt.root 2>/dev/null || true
ln -sf "output_tracking_${TAG}.root" output_tracking_goodevt.root 2>/dev/null || true
root -b -q -l "$HERE/macros/gluckstern_floor.C(1.0)"
mv gluckstern_floor.png "$PLOTS/" 2>/dev/null || true
rm -f sim_goodevt.root output_tracking_goodevt.root

echo "--- outlier profile ---"
root -b -q -l "$HERE/macros/profile_outliers.C(\"${TAG}\")"
mv outlier_profile.png "$PLOTS/" 2>/dev/null || true

echo "===== event displays ====="
root -b -q -l "$HERE/macros/evtdisp_z.C"
mv evtdisp_z_horizontal.png "$PLOTS/" 2>/dev/null || true
root -b -q -l "$HERE/macros/evtdisp_z_multi.C"
mv evtdisp_z_multi.png "$PLOTS/" 2>/dev/null || true
root -b -q -l "$HERE/macros/evtdisp_z_he3.C"
mv evtdisp_z_he3.png "$PLOTS/" 2>/dev/null || true

echo
echo "Done. Plots in $PLOTS/"
ls -1 "$PLOTS/"
