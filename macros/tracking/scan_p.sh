#!/usr/bin/env bash
#
# σ_p/p scan for the R3B HYDRA UKF pipeline, mirroring the ATTPCROOT
# scan_zfix_widefield_aim_2k driver. For each momentum point P (MeV/c):
#   1. sim   (FairBoxGenerator @ p=P, narrow theta/phi cone aimed at the
#             Prototype chamber)
#   2. lang  (Langevin drift)
#   3. reco  (Cal2Hit)
#   4. tracking (TripletClust → R3BGTPCTrack2Find with Kasa geo fit)
#   5. ukf   (R3BGTPCTrack2Fit; default = back-extrap on, helix POCA, no
#             straight tail since R3B's field map covers upstream too)
#
# Outputs go under glad-tpc/macros/{sim,proj,reco,tracking}/output_*_pX.root
# via the SUFFIX=_pX env var the stage macros now honour.
#
# Usage:
#   cd glad-tpc/macros/tracking
#   ./scan_p.sh

set -eo pipefail
: "${LD_LIBRARY_PATH:=}"
: "${ROOT_INCLUDE_PATH:=}"

NEVT=${NEVT:-1000}
PLIST=${PLIST:-"200 400 600 800 1000 1200"}

# Resolve install paths (the R3BRoot build dir holds the rootmaps; FairRoot
# and FairSoft are needed for ROOT + Geant4 + VMC).
: "${VMCWORKDIR:=$HOME/fair_install/R3BRoot}"
: "${FAIRROOTPATH:=$HOME/fair_install/FairRootInstall}"
: "${SIMPATH:=$HOME/fair_install/FairSoft/install}"
# ATTPCROOT's config.sh sets the Geant4 data env vars we need; source it
# first, then re-export VMCWORKDIR so it points at R3BRoot (config.sh
# clobbers it) and prepend R3BRoot's build/lib to LD_LIBRARY_PATH.
source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh"
VMCWORKDIR=$HOME/fair_install/R3BRoot
export VMCWORKDIR FAIRROOTPATH SIMPATH
export LD_LIBRARY_PATH="$VMCWORKDIR/build/lib:$LD_LIBRARY_PATH"

OUT_CSV="$VMCWORKDIR/glad-tpc/macros/tracking/scan_p_results.csv"
echo "p_MeV,n_fit,n_thr,sigma_p_over_p_core,bias_core" >"$OUT_CSV"

mkdir -p "$VMCWORKDIR/glad-tpc/macros/sim/Prototype"
mkdir -p "$VMCWORKDIR/glad-tpc/macros/proj/Prototype"

for P in $PLIST; do
   export P_MEV=$P
   export SUFFIX=_p${P}
   # Use the Kasa-derived seed (no truth override). With B = ŷ the bending
   # plane is (x, z) and the Kasa fit can now extract a meaningful radius.
   unset SEED_P_MEV
   echo "===== p=$P MeV/c ====="

   pushd "$VMCWORKDIR/glad-tpc/macros/sim" >/dev/null
   root -b -q -e ".L simHYDRA.C" -e "simHYDRA($NEVT, \"Prototype\", \"box\")" \
        >"/tmp/r3b_scan_sim${SUFFIX}.log" 2>&1
   popd >/dev/null

   pushd "$VMCWORKDIR/glad-tpc/macros/proj" >/dev/null
   root -b -q "run_lang.C(\"Prototype\")" \
        >"/tmp/r3b_scan_lang${SUFFIX}.log" 2>&1
   popd >/dev/null

   pushd "$VMCWORKDIR/glad-tpc/macros/reco" >/dev/null
   root -b -q "run_reconstruction.C(\"lang${SUFFIX}.root\")" \
        >"/tmp/r3b_scan_reco${SUFFIX}.log" 2>&1
   popd >/dev/null

   pushd "$VMCWORKDIR/glad-tpc/macros/tracking" >/dev/null
   root -b -q "run_tracking.C(\"output_reco${SUFFIX}.root\")" \
        >"/tmp/r3b_scan_track${SUFFIX}.log" 2>&1
   root -b -q "run_ukf.C(\"output_tracking${SUFFIX}.root\", \"output_ukf${SUFFIX}.root\")" \
        >"/tmp/r3b_scan_ukf${SUFFIX}.log" 2>&1

   # Tally σ_p/p directly from the per-point output.
   LINE=$(root -b -q "analyze_scan_point.C(${P})" 2>&1 \
          | awk '/^RESULT/{print}')
   echo "$LINE"
   echo "$LINE" | awk '{print $2","$3","$4","$5","$6}' >>"$OUT_CSV"
   popd >/dev/null

done

echo "===== Done. Summary: ====="
column -t -s, "$OUT_CSV"
