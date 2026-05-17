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

HERE=$(cd "$(dirname "$0")" && pwd)
TRACKDIR="$VMCWORKDIR/glad-tpc/macros/tracking"
OUT_CSV="$HERE/scan_p_results.csv"
LOG_DIR=/tmp/r3b_scan_p
mkdir -p "$LOG_DIR" "$HERE/plots"

require_file() { [ ! -s "$1" ] && { echo "ERROR: $1 missing"; return 1; } || return 0; }

echo "p_MeV,N_fit,sigma_R_pct,sigma_p_seed_pct,sigma_p_UKF_pct,bias_R_pct,bias_p_seed_pct,bias_p_UKF_pct" >"$OUT_CSV"

for P in $PLIST; do
    export P_MEV=$P
    export SUFFIX=_p${P}
    echo "===== p = $P MeV/c ====="

    pushd "$VMCWORKDIR/glad-tpc/macros/sim" >/dev/null
    SUFFIX="_p$P" P_MEV=$P root -b -q -l "simHYDRA.C($NEVT, \"Prototype\", \"box\")" \
        >"$LOG_DIR/sim_p${P}.log" 2>&1 || true
    require_file "Prototype/sim_p${P}.root" || { popd >/dev/null; continue; }
    popd >/dev/null

    pushd "$VMCWORKDIR/glad-tpc/macros/proj" >/dev/null
    SUFFIX="_p$P" root -b -q -l run_lang.C >"$LOG_DIR/lang_p${P}.log" 2>&1 || true
    require_file "Prototype/lang_p${P}.root" || { popd >/dev/null; continue; }
    popd >/dev/null

    pushd "$VMCWORKDIR/glad-tpc/macros/reco" >/dev/null
    SUFFIX="_p$P" root -b -q -l "run_reconstruction.C(\"lang_p${P}.root\")" \
        >"$LOG_DIR/reco_p${P}.log" 2>&1 || true
    require_file "output_reco_p${P}.root" || { popd >/dev/null; continue; }
    popd >/dev/null

    pushd "$TRACKDIR" >/dev/null
    USE_VERTEX=1 HUBER_K_CM=0.10 SUFFIX="_p$P" \
        root -b -q -l "run_tracking.C(\"output_reco_p${P}.root\")" \
        >"$LOG_DIR/tracking_p${P}.log" 2>&1 || true
    require_file "output_tracking_p${P}.root" || { popd >/dev/null; continue; }

    SUFFIX="_p$P" root -b -q -l run_ukf.C >"$LOG_DIR/ukf_p${P}.log" 2>&1 || true
    require_file "output_ukf_p${P}.root" || { popd >/dev/null; continue; }

    LINE=$(root -b -q -l "$HERE/macros/analyze_scan_point.C($P)" 2>/dev/null \
           | awk '/^RESULT/{print}')
    echo "  $LINE"
    echo "$LINE" | awk '{print $2","$3","$4","$5","$6","$7","$8","$9}' >>"$OUT_CSV"
    popd >/dev/null
done

echo
echo "===== Summary ====="
column -t -s, "$OUT_CSV"

# Generate the sweep plot
root -b -q -l "$HERE/macros/plot_scan_p.C" 2>/dev/null

echo
echo "Outputs:"
echo "  $OUT_CSV"
echo "  $HERE/plots/sigma_vs_p.png"
