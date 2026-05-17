/******************************************************************************
 *   Copyright (C) 2018-2026 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 ******************************************************************************/

#pragma once

#include "AtELossModel.h"

#include <Math/Vector3D.h>
#include <Rtypes.h>
#include <limits>
#include <memory>

class R3BGTPCTrackData;
class R3BGTPCFittedTrackData;

namespace kf
{
class TrackFitterUKF;
}

/// Unscented-Kalman track fitter for the GLAD-TPC.
///
/// Companion to the existing GenFit-based R3BGTPCFitter; this one uses the
/// OpenKF UKF (vendored under glad-tpc/openkf/). State vector is the same
/// 6-tuple (x, y, z, |p|, theta, phi) ATTPCROOT's AtFitterUKF uses, so the
/// experience there ports across with minor adaptations.
///
/// Typical usage:
///     auto eloss = std::make_unique<AtTools::AtELossCATIMA>(rho);
///     eloss->SetProjectile(A, Z, mAmu);
///     eloss->SetMaterial({{1, 1, 1}});
///     R3BGTPCFitterUKF fitter(chargeC, mass_MeV, std::move(eloss));
///     fitter.SetBField({0, 0, 2.0});  // Tesla
///     auto fitted = fitter.FitTrack(track);
class R3BGTPCFitterUKF
{
  public:
    R3BGTPCFitterUKF(double charge, double mass_MeV, std::unique_ptr<AtTools::AtELossModel> elossModel);
    ~R3BGTPCFitterUKF();

    // -- configuration
    void SetBField(const ROOT::Math::XYZVector& bField) { fBField = bField; }
    void SetMomentumSeed(double p_MeV) { fMomentumSeed = p_MeV; }
    void SetMeasurementSigma(double sigma_mm) { fMeasSigma_mm = sigma_mm; }
    void SetMomentumSigmaFrac(double frac) { fMomSigmaFrac = frac; }
    void SetMinClusters(int n) { fMinClusters = n; }
    /// Maximum seed circle radius (cm) accepted by the UKF. The Pratt+GN
    /// fit clamps R at 20 m to avoid straight-line catastrophes, but seeds
    /// near the cap have no real curvature information — feeding them to a
    /// tight-prior UKF lock the fit to a wildly-wrong p. Default 500 cm
    /// covers all physical pion radii (R = 67 cm at 400 MeV/c, 200 cm at
    /// 1.2 GeV/c). Set <= 0 to disable.
    void SetMaxSeedRadius_cm(double r) { fMaxSeedRadius_cm = r; }
    void SetEnableEnergyStraggling(bool enable) { fEnableEnStraggling = enable; }
    void SetELossScaleFactor(double f) { fELossScaleFactor = f; }

    /// Length-unit conversion: multiply R3BGTPCTrackData positions by this
    /// factor before feeding them to the UKF (and divide output positions by
    /// it before storing). R3BRoot stores hits in cm (FairRoot convention)
    /// so default = 10. ATTPCROOT users would pass 1.0.
    void SetInputUnit_mm(double factor) { fInputUnit_mm = factor; }

    // -- back-extrapolation knobs (all distances in mm internally)
    /// Closed-form POCA on the PRA circle + helix-pitch z. When off, falls
    /// back to a linear step along the initial momentum direction.
    void SetUseHelixBackExtrap(bool on) { fUseHelixBackExtrap = on; }
    /// Cap on the back-extrapolated arc length (mm). Prevents the helix
    /// from wrapping around for low-momentum tracks.
    void SetBackExtrapMaxPath(double mm) { fBackExtrapMaxPath = mm; }
    /// Force the back-extrapolated vertex onto the beam axis by appending a
    /// chord from POCA to (0,0). Useful when the PRA circle has a few-mm
    /// offset from the true vertex.
    void SetForceVertexOnBeamAxis(bool on) { fForceVertexOnBeamAxis = on; }
    /// Rotate Kinematics.phi by arc/R during back-extrap so the stored
    /// azimuth is at the vertex, not at the first cluster.
    void SetUpdateAnglesOnBackExtrap(bool on) { fUpdateAnglesOnBackExtrap = on; }
    /// Optional straight-line tail from POCA to this x-plane (mm, input
    /// frame). For setups where the production vertex sits in a field-free
    /// region upstream of the chamber. NaN (default) disables.
    void SetBackExtrapTargetX(double x_in) { fBackExtrapTargetX = x_in; }
    void SetUKFParameters(double alpha, double beta, double kappa)
    {
        fAlpha = alpha;
        fBeta = beta;
        fKappa = kappa;
    }

    /// Fit a single PRA track. Returns nullptr if the track is rejected (too
    /// few clusters, bad geometry, etc.). Otherwise returns a fitted-track
    /// object with kinematics, vertex, smoothed positions, and chi^2/ndf.
    std::unique_ptr<R3BGTPCFittedTrackData> FitTrack(R3BGTPCTrackData* track);

  private:
    void InitUKF();
    double GetSeedMomentum(R3BGTPCTrackData* track) const;

    double fCharge;     ///< elementary charges (Coulombs)
    double fMass_MeV;   ///< rest mass in MeV/c^2
    std::unique_ptr<AtTools::AtELossModel> fELossModel;
    std::unique_ptr<kf::TrackFitterUKF> fUKF;

    ROOT::Math::XYZVector fBField{ 0., 0., 2.0 };
    double fMomentumSeed{ -1.0 }; ///< -1: derive from Brho (GeoRadius)
    double fMeasSigma_mm{ 1.0 };
    double fMomSigmaFrac{ 0.1 };     ///< matches AtFitterUKF default
    bool fEnableEnStraggling{ true };
    double fELossScaleFactor{ 1.0 };
    // Minimum input clusters to attempt UKF. At 5 the chamber's short-chord
    // tracks (400-600 MeV/c, chord 4-7 cm) are silently rejected because
    // TripletClust forms <5 clusters. Lowering to 3 recovers σ_p/p UKF at
    // 400 MeV/c from 42 % (small-N Gauss artefact) to 3.8 % at a cost of
    // +0.2-0.4 pp at 500-600 MeV/c.
    int fMinClusters{ 3 };
    double fMaxSeedRadius_cm{ 500.0 };
    double fAlpha{ 1e-3 };
    double fBeta{ 2.0 };
    double fKappa{ 0.0 };
    double fInputUnit_mm{ 10.0 }; ///< R3B default: input positions are cm → ×10 to mm

    // -- back-extrap
    bool fUseHelixBackExtrap{ true };
    bool fForceVertexOnBeamAxis{ false };
    bool fUpdateAnglesOnBackExtrap{ true };
    double fBackExtrapMaxPath{ 450.0 };
    double fBackExtrapTargetX{ std::numeric_limits<double>::quiet_NaN() };
};
