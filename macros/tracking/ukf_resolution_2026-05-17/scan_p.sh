#!/usr/bin/env bash
# σ_p/p vs pion momentum sweep on HYDRA Prototype.
#
# For each P in PLIST: 500-event box-gen sim of π⁻ at p=P, full chain
# (sim → lang → reco → tracking [vertex+Huber] → UKF), per-point σ_p/p
# extraction via analyze_scan_point.C. Writes CSV + summary plot.
#
# Usage:
#   source build/config.sh                # set VMCWORKDIR
#   cd glad-tpc/macros/tracking/ukf_resolution_2026-05-17
#   ./scan_p.sh                            # default PLIST
#   PLIST="300 600 900" NEVT=200 ./scan_p.sh
#
# Per-point output files live in their canonical macros/{sim,proj,reco,
# tracking}/ folders with suffix _pX so multiple points coexist on disk.

set -uo pipefail
[ -z "${VMCWORKDIR:-}" ] && { echo "ERROR: VMCWORKDIR not set"; exit 1; }

NEVT=${NEVT:-500}
PLIST=${PLIST:-"200 300 400 500 600 700 800 900 1000 1200"}
# STEMAX=0.1 cm is the production default: it densifies MC truth (~22× more
# GTPCPoints) so Langevin lays drift electrons along the true curve, removing
# the +2 % bias drift seen with sparse MC. Override to "" (empty) for the
# historical no-limit baseline (see scan_p_results_nostemax.csv).
export GTPC_STEMAX_CM=${GTPC_STEMAX_CM:-0.1}
# Optional TAG distinguishes parallel runs (e.g. TAG=_nostemax). Empty = default.
TAG=${TAG:-}

HERE=$(cd "$(dirname "$0")" && pwd)
TRACKDIR="$VMCWORKDIR/glad-tpc/macros/tracking"
OUT_CSV="$HERE/scan_p_results${TAG}.csv"
LOG_DIR=/tmp/r3b_scan_p${TAG}
mkdir -p "$LOG_DIR" "$HERE/plots"

require_file() { [ ! -s "$1" ] && { echo "ERROR: $1 missing"; return 1; } || return 0; }

echo "p_MeV,N_fit,sigma_R_pct,sigma_p_seed_pct,sigma_p_UKF_pct,bias_R_pct,bias_p_seed_pct,bias_p_UKF_pct" >"$OUT_CSV"

for P in $PLIST; do
    export P_MEV=$P
    SFX="${TAG}_p${P}"
    export SUFFIX=$SFX
    echo "===== p = $P MeV/c   (SUFFIX=$SFX) ====="

    pushd "$VMCWORKDIR/glad-tpc/macros/sim" >/dev/null
    SUFFIX="$SFX" P_MEV=$P root -b -q -l "simHYDRA.C($NEVT, \"Prototype\", \"box\")" \
        >"$LOG_DIR/sim${SFX}.log" 2>&1 || true
    require_file "Prototype/sim${SFX}.root" || { popd >/dev/null; continue; }
    popd >/dev/null

    pushd "$VMCWORKDIR/glad-tpc/macros/proj" >/dev/null
    SUFFIX="$SFX" root -b -q -l run_lang.C >"$LOG_DIR/lang${SFX}.log" 2>&1 || true
    require_file "Prototype/lang${SFX}.root" || { popd >/dev/null; continue; }
    popd >/dev/null

    pushd "$VMCWORKDIR/glad-tpc/macros/reco" >/dev/null
    SUFFIX="$SFX" root -b -q -l "run_reconstruction.C(\"lang${SFX}.root\")" \
        >"$LOG_DIR/reco${SFX}.log" 2>&1 || true
    require_file "output_reco${SFX}.root" || { popd >/dev/null; continue; }
    popd >/dev/null

    pushd "$TRACKDIR" >/dev/null
    USE_VERTEX=1 HUBER_K_CM=0.10 SUFFIX="$SFX" \
        root -b -q -l "run_tracking.C(\"output_reco${SFX}.root\")" \
        >"$LOG_DIR/tracking${SFX}.log" 2>&1 || true
    require_file "output_tracking${SFX}.root" || { popd >/dev/null; continue; }

    SUFFIX="$SFX" root -b -q -l run_ukf.C >"$LOG_DIR/ukf${SFX}.log" 2>&1 || true
    require_file "output_ukf${SFX}.root" || { popd >/dev/null; continue; }

    LINE=$(root -b -q -l "$HERE/macros/analyze_scan_point.C($P,\"$TAG\")" 2>/dev/null \
           | awk '/^RESULT/{print}')
    echo "  $LINE"
    echo "$LINE" | awk '{print $2","$3","$4","$5","$6","$7","$8","$9}' >>"$OUT_CSV"
    popd >/dev/null
done

echo
echo "===== Summary ====="
column -t -s, "$OUT_CSV"

# Generate the sweep plot only for the baseline (tag-less) run; comparison
# plots are handled separately so a STEMAX rerun doesn't overwrite them.
if [ -z "$TAG" ]; then
    root -b -q -l "$HERE/macros/plot_scan_p.C" 2>/dev/null
    echo
    echo "Outputs:"
    echo "  $OUT_CSV"
    echo "  $HERE/plots/sigma_vs_p.png"
else
    echo
    echo "Outputs:"
    echo "  $OUT_CSV"
fi
