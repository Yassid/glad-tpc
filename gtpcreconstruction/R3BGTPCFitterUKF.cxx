/******************************************************************************
 *   Copyright (C) 2018-2026 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 ******************************************************************************/

#include "R3BGTPCFitterUKF.h"

#include "AtKinematics.h"
#include "AtPropagator.h"
#include "R3BGTPCFittedTrackData.h"
#include "R3BGTPCHitClusterData.h"
#include "R3BGTPCTrackData.h"
#include "TrackFitterUKF.h"

#include <FairLogger.h>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <TMatrixD.h>

#include <algorithm>
#include <cmath>

namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr double kInitAngSigma_deg = 1.0;
}

R3BGTPCFitterUKF::R3BGTPCFitterUKF(double charge,
                                   double mass_MeV,
                                   std::unique_ptr<AtTools::AtELossModel> elossModel)
    : fCharge(charge)
    , fMass_MeV(mass_MeV)
    , fELossModel(std::move(elossModel))
{
}

R3BGTPCFitterUKF::~R3BGTPCFitterUKF() = default;

void R3BGTPCFitterUKF::InitUKF()
{
    AtTools::AtPropagator propagator(fCharge, fMass_MeV, std::move(fELossModel));
    propagator.SetBField(fBField);

    auto stepper = std::make_unique<AtTools::AtRK4Stepper>();
    fUKF = std::make_unique<kf::TrackFitterUKF>(std::move(propagator), std::move(stepper));
    fUKF->setParameters(static_cast<float>(fAlpha), static_cast<float>(fBeta), static_cast<float>(fKappa));
    fUKF->fEnableEnStraggling = fEnableEnStraggling;
    fUKF->fELossScaleFactor = fELossScaleFactor;
}

double R3BGTPCFitterUKF::GetSeedMomentum(R3BGTPCTrackData* track) const
{
    if (fMomentumSeed > 0)
        return fMomentumSeed;

    // Brho seed from the PRA circle: p_T = 0.3 · B · R, then divide by sin(θ).
    // GeoRadius is in mm, B field in T, theta in rad. Returns p in MeV/c.
    const double R_mm = track->GetGeoRadius();
    if (!std::isfinite(R_mm) || R_mm <= 0)
    {
        LOG(warn) << "R3BGTPCFitterUKF: GeoRadius not set or invalid (" << R_mm
                  << "), falling back to 100 MeV/c seed";
        return 100.0;
    }

    double theta = track->GetGeoTheta();
    if (!std::isfinite(theta))
        theta = kPi / 2.0;
    double sinTheta = std::sin(theta);
    if (std::abs(sinTheta) < 0.1)
        sinTheta = (sinTheta < 0 ? -0.1 : 0.1);

    const double Brho = std::abs(fBField.Z() * R_mm / 1000.0 / sinTheta); // T·m
    return Brho * 0.3 * 1000.0;                                           // → MeV/c
}

std::unique_ptr<R3BGTPCFittedTrackData> R3BGTPCFitterUKF::FitTrack(R3BGTPCTrackData* track)
{
    if (track == nullptr)
        return nullptr;
    if (fUKF == nullptr)
        InitUKF();

    auto* clusters = track->GetHitClusterArray();
    if (clusters == nullptr || static_cast<int>(clusters->size()) < fMinClusters)
    {
        LOG(info) << "R3BGTPCFitterUKF: skipping track " << track->GetTrackId() << " — "
                  << (clusters ? clusters->size() : 0) << " clusters < " << fMinClusters;
        return nullptr;
    }

    // -- 1. Seed momentum + initial pose from first two clusters
    const double p_seed = GetSeedMomentum(track);
    const auto& c0 = clusters->at(0);
    ROOT::Math::XYZPoint initialPos(c0.GetX(), c0.GetY(), c0.GetZ());

    ROOT::Math::XYZVector initialMom;
    const auto& c1 = clusters->at(1);
    ROOT::Math::XYZVector dir(c1.GetX() - c0.GetX(), c1.GetY() - c0.GetY(), c1.GetZ() - c0.GetZ());
    if (dir.R() > 1e-6)
        initialMom = p_seed * dir.Unit();
    else
        initialMom = ROOT::Math::XYZVector(p_seed, 0, 0);

    // -- 2. Initial covariance: σ_pos = fMeasSigma_mm, σ_p = fMomSigmaFrac·p,
    //       σ_angles = 1° (matches AtFitterUKF::GetInitialCovariance)
    TMatrixD P0(6, 6);
    P0.Zero();
    const double sig_pos2 = fMeasSigma_mm * fMeasSigma_mm;
    P0(0, 0) = P0(1, 1) = P0(2, 2) = sig_pos2;
    const double sig_mom = fMomSigmaFrac * p_seed;
    P0(3, 3) = sig_mom * sig_mom;
    const double angSig = kInitAngSigma_deg * kPi / 180.0;
    P0(4, 4) = P0(5, 5) = angSig * angSig;

    fUKF->SetInitialState(initialPos, initialMom, P0);

    TMatrixD R0(3, 3);
    R0.Zero();
    R0(0, 0) = R0(1, 1) = R0(2, 2) = fMeasSigma_mm * fMeasSigma_mm;
    fUKF->SetMeasCov(R0);

    // -- 3. Forward filter pass
    bool converged = true;
    try
    {
        for (size_t i = 1; i < clusters->size(); ++i)
        {
            const auto& ci = clusters->at(i);
            ROOT::Math::XYZPoint meas(ci.GetX(), ci.GetY(), ci.GetZ());

            fUKF->predictUKF(meas);

            // Sanity check the predicted state before correcting; UKF can
            // produce non-finite components when sigma-points blow up.
            const auto& state = fUKF->vecX();
            bool stateValid = true;
            for (int d = 0; d < 6; ++d)
            {
                if (!std::isfinite(state[d]) || std::abs(state[d]) > 1e6)
                {
                    stateValid = false;
                    break;
                }
            }
            if (!stateValid)
            {
                LOG(debug) << "R3BGTPCFitterUKF: invalid state at cluster " << i << "; aborting fit";
                converged = false;
                break;
            }

            fUKF->correctUKF(meas);
        }
    }
    catch (const std::exception& e)
    {
        LOG(warn) << "R3BGTPCFitterUKF: forward pass exception: " << e.what();
        converged = false;
    }

    // -- 4. RTS smoother
    if (converged)
    {
        try
        {
            fUKF->smoothUKF();
        }
        catch (const std::exception& e)
        {
            LOG(warn) << "R3BGTPCFitterUKF: smoother exception: " << e.what();
            converged = false;
        }
    }

    auto fitted = std::make_unique<R3BGTPCFittedTrackData>();
    fitted->SetTrackId(track->GetTrackId());
    fitted->SetConverged(converged);

    if (!converged)
        return fitted;

    // -- 5. Extract kinematics from the smoothed state at first cluster.
    //       Back-extrapolation to the production vertex is a TODO; for now
    //       fKinematics == fKinematicsXtr (both at first cluster).
    const auto& smoothed = fUKF->GetSmoothedStates();
    if (smoothed.empty())
        return fitted;

    const auto& s0 = smoothed.front();
    const double vx = s0[0];
    const double vy = s0[1];
    const double vz = s0[2];
    const double p_first = s0[3];
    const double theta_first = s0[4];
    const double phi_first = s0[5];
    const double KE_first = std::sqrt(p_first * p_first + fMass_MeV * fMass_MeV) - fMass_MeV;

    fitted->SetVertex(ROOT::Math::XYZVector(vx, vy, vz));
    fitted->SetKinematicsXtr(KE_first, theta_first, phi_first);
    fitted->SetKinematics(KE_first, theta_first, phi_first);

    std::vector<ROOT::Math::XYZPoint> smoothedPositions;
    smoothedPositions.reserve(smoothed.size());
    for (const auto& s : smoothed)
        smoothedPositions.emplace_back(s[0], s[1], s[2]);
    fitted->SetSmoothedPositions(std::move(smoothedPositions));

    // -- 6. chi^2 / ndf (rough: sum of (measurement - smoothed)^2 / sigma^2)
    double chi2 = 0;
    int ndf = 0;
    const double sig2 = fMeasSigma_mm * fMeasSigma_mm;
    for (size_t i = 0; i < clusters->size() && i < smoothed.size(); ++i)
    {
        const auto& ci = clusters->at(i);
        const auto& si = smoothed[i];
        const double dx = ci.GetX() - si[0];
        const double dy = ci.GetY() - si[1];
        const double dz = ci.GetZ() - si[2];
        chi2 += (dx * dx + dy * dy + dz * dz) / sig2;
        ndf += 3;
    }
    ndf -= 6; // subtract free parameters
    fitted->SetChi2(chi2);
    fitted->SetNdf(std::max(ndf, 1));

    return fitted;
}
