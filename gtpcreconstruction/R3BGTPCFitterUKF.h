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
    void SetEnableEnergyStraggling(bool enable) { fEnableEnStraggling = enable; }
    void SetELossScaleFactor(double f) { fELossScaleFactor = f; }
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
    int fMinClusters{ 5 };
    double fAlpha{ 1e-3 };
    double fBeta{ 2.0 };
    double fKappa{ 0.0 };
};
