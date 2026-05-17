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
    // GeoRadius is in input units (cm in R3B default, mm in ATTPCROOT) and is
    // scaled to mm via fInputUnit_mm. B field in T, theta in rad. Returns p in MeV/c.
    const double R_in = track->GetGeoRadius();
    if (!std::isfinite(R_in) || R_in <= 0)
    {
        LOG(warn) << "R3BGTPCFitterUKF: GeoRadius not set or invalid (" << R_in
                  << "), falling back to 100 MeV/c seed";
        return 100.0;
    }
    const double R_mm = R_in * fInputUnit_mm;

    double theta = track->GetGeoTheta();
    if (!std::isfinite(theta))
        theta = kPi / 2.0;
    double sinTheta = std::sin(theta);
    if (std::abs(sinTheta) < 0.1)
        sinTheta = (sinTheta < 0 ? -0.1 : 0.1);

    // R3B GLAD: B = (0, B_y, 0). The Kasa fit in R3BGTPCTrackFinder gives
    // GeoRadius in the bending plane (x,z); GeoTheta is the angle from the
    // field direction (+y). Brho uses |B| (magnitude) for generality.
    const double Bmag = std::sqrt(fBField.X() * fBField.X() + fBField.Y() * fBField.Y() + fBField.Z() * fBField.Z());
    const double Brho = std::abs(Bmag * R_mm / 1000.0 / sinTheta); // T·m
    return Brho * 0.3 * 1000.0;                                    // → MeV/c
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

    // Seed-quality cut: with the tight 0.02 momentum prior, the UKF locks
    // to the seed. Seeds near the Pratt+GN R cap (20 m) carry no curvature
    // information and produce wildly wrong p; better to reject them here
    // than let them poison the σ_p/p histogram. Cf. plots/probe_400.png.
    const double R_seed_cm = track->GetGeoRadius();
    if (fMaxSeedRadius_cm > 0 && std::isfinite(R_seed_cm) && R_seed_cm > fMaxSeedRadius_cm)
    {
        LOG(info) << "R3BGTPCFitterUKF: skipping track " << track->GetTrackId()
                  << " — seed R_fit = " << R_seed_cm << " cm > " << fMaxSeedRadius_cm
                  << " cm (degenerate, near GN R cap)";
        return nullptr;
    }

    // -- 1. Seed momentum + initial pose from first two clusters. All
    //       positions are scaled to mm via fInputUnit_mm so the UKF runs in
    //       the same unit system as ATTPCROOT (where AtPropagator was
    //       written and validated).
    const double p_seed = GetSeedMomentum(track);
    const auto& c0 = clusters->at(0);
    ROOT::Math::XYZPoint initialPos(c0.GetX() * fInputUnit_mm,
                                    c0.GetY() * fInputUnit_mm,
                                    c0.GetZ() * fInputUnit_mm);

    // Build the seed momentum direction from the PRA helix params (R,
    // GeoCenter, θ_y) rather than the noisy c1−c0 difference. The tangent
    // in the (x,z) bending plane at first cluster is perpendicular to the
    // radial direction from (cx,cz)→c0; the y-component is governed by
    // θ_y (the angle from +y of the helix). Sign of the in-plane tangent
    // and of py is fixed by the second cluster.
    ROOT::Math::XYZVector initialMom;
    auto cen = track->GetGeoCenter();
    const double cx = cen.first, cz = cen.second;
    const double theta_y = track->GetGeoTheta();
    const auto& c1 = clusters->at(1);
    bool seedFromHelix = std::isfinite(cx) && std::isfinite(cz) && std::isfinite(theta_y);
    if (seedFromHelix)
    {
        const double rx = c0.GetX() - cx;
        const double rz = c0.GetZ() - cz;
        const double rn = std::hypot(rx, rz);
        if (rn < 1e-6) seedFromHelix = false;
        if (seedFromHelix)
        {
            // Two tangent candidates; pick the one pointing toward c1.
            const double t1x = -rz / rn, t1z = rx / rn;
            const double dx = c1.GetX() - c0.GetX(), dz = c1.GetZ() - c0.GetZ();
            const double sgn = (t1x * dx + t1z * dz >= 0) ? 1.0 : -1.0;
            const double tx_xz = sgn * t1x, tz_xz = sgn * t1z;
            const double sinTh = std::sin(theta_y);
            const double cosTh = std::cos(theta_y);
            const double dy_sgn = (c1.GetY() - c0.GetY() >= 0) ? 1.0 : -1.0;
            // 3D unit direction: in-plane component magnitude sin(θ_y), y comp cos(θ_y) with the
            // sign matching the observed y drift between the first two clusters.
            const double dx3 = sinTh * tx_xz;
            const double dy3 = dy_sgn * cosTh;
            const double dz3 = sinTh * tz_xz;
            initialMom = p_seed * ROOT::Math::XYZVector(dx3, dy3, dz3);
        }
    }
    if (!seedFromHelix)
    {
        ROOT::Math::XYZVector dir(c1.GetX() - c0.GetX(), c1.GetY() - c0.GetY(), c1.GetZ() - c0.GetZ());
        if (dir.R() > 1e-6)
            initialMom = p_seed * dir.Unit();
        else
            initialMom = ROOT::Math::XYZVector(p_seed, 0, 0);
    }

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
            ROOT::Math::XYZPoint meas(ci.GetX() * fInputUnit_mm,
                                      ci.GetY() * fInputUnit_mm,
                                      ci.GetZ() * fInputUnit_mm);

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

    // -- 5. Extract kinematics from the smoothed state at first cluster +
    //       back-extrapolate to the production vertex along the PRA helix.
    //       Two outputs:
    //         KinematicsXtr — at the first cluster (pre-back-extrap, purely
    //                         measurement-driven; useful sanity check).
    //         Kinematics    — at the back-extrapolated vertex (canonical:
    //                         this is what physics analyses read).
    const auto& smoothed = fUKF->GetSmoothedStates();
    if (smoothed.empty())
        return fitted;

    const auto& s0 = smoothed.front();
    double vx = s0[0]; // mm
    double vy = s0[1];
    double vz = s0[2];
    double p_s = s0[3];
    double theta_s = s0[4];
    double phi_s = s0[5];
    const double KE_first = std::sqrt(p_s * p_s + fMass_MeV * fMass_MeV) - fMass_MeV;
    fitted->SetKinematicsXtr(KE_first, theta_s, phi_s);

    // --- Back-extrapolation to beam axis along the PRA helix ---
    // Ported from AtFitterUKF (commit cbc97ade); see there for the
    // derivation of the POCA-on-circle closed form and the chord append for
    // ForceVertexOnBeamAxis. Distances are mm; input frame conversion is
    // done at the boundary via fInputUnit_mm.
    if (p_s > 0 && fBackExtrapMaxPath > 0.)
    {
        double pathLength = 0.;

        // PRA circle parameters from the upstream Kasa fit, scaled to mm.
        const auto geoCenter = track->GetGeoCenter();
        const double cx = geoCenter.first * fInputUnit_mm;
        const double cy = geoCenter.second * fInputUnit_mm;
        const double R = track->GetGeoRadius() * fInputUnit_mm;
        const double dCenter = std::sqrt(cx * cx + cy * cy);
        const bool circleValid = fUseHelixBackExtrap && std::isfinite(R) && (R > 1.0) && (dCenter > 1.0);

        if (circleValid)
        {
            // POCA of the PRA circle to (0, 0). The closer of the two
            // extrema is (cx, cy) · (1 − R / d) for any d > 0.
            const double f = 1.0 - R / dCenter;
            const double pocaX = cx * f;
            const double pocaY = cy * f;

            // Use the raw hit closest to the beam axis as the back-extrap
            // start (instead of the cluster-centroid smoothed state), which
            // removes the half-cluster-spacing · cot(θ) z-bias on upward
            // tracks.
            double rxFirst = vx;
            double ryFirst = vy;
            double rzFirst = vz;
            {
                const auto& hits = track->GetHitArray();
                double bestR2 = rxFirst * rxFirst + ryFirst * ryFirst;
                for (const auto& h : hits)
                {
                    const double hx = h.GetX() * fInputUnit_mm;
                    const double hy = h.GetY() * fInputUnit_mm;
                    const double hz = h.GetZ() * fInputUnit_mm;
                    const double r2 = hx * hx + hy * hy;
                    if (r2 < bestR2)
                    {
                        bestR2 = r2;
                        rxFirst = hx;
                        ryFirst = hy;
                        rzFirst = hz;
                    }
                }
            }

            // Arc length along the PRA circle from raw first hit to POCA;
            // pick the shorter wrap.
            const double phi1 = std::atan2(ryFirst - cy, rxFirst - cx);
            const double phi0 = std::atan2(pocaY - cy, pocaX - cx);
            double dPhi = std::abs(phi1 - phi0);
            if (dPhi > kPi)
                dPhi = 2 * kPi - dPhi;
            double arc = R * dPhi;

            double endX = pocaX;
            double endY = pocaY;
            if (fForceVertexOnBeamAxis)
            {
                arc += std::sqrt(pocaX * pocaX + pocaY * pocaY);
                endX = 0.0;
                endY = 0.0;
            }
            arc = std::min(arc, fBackExtrapMaxPath);

            const double sinTheta = std::max(std::sin(theta_s), 0.1);
            const double cotTheta = std::cos(theta_s) / sinTheta;

            vx = endX;
            vy = endY;
            vz = rzFirst - arc * cotTheta;
            pathLength = arc;

            if (fUpdateAnglesOnBackExtrap)
            {
                const double dphi_mag = arc / R;
                const double signFactor = (fCharge < 0) ? -1.0 : +1.0;
                phi_s += signFactor * dphi_mag;
                while (phi_s > kPi)
                    phi_s -= 2 * kPi;
                while (phi_s <= -kPi)
                    phi_s += 2 * kPi;
            }
        }
        else
        {
            // Linear fallback (used when the PRA circle is degenerate).
            const double rXY = std::sqrt(vx * vx + vy * vy);
            const double sinTheta = std::sin(theta_s);
            pathLength = (sinTheta > 0.1) ? rXY / sinTheta : rXY;
            pathLength = std::min(pathLength, fBackExtrapMaxPath);
            ROOT::Math::Polar3DVector momDir(1.0, theta_s, phi_s);
            ROOT::Math::XYZVector dir(momDir);
            vx -= dir.X() * pathLength;
            vy -= dir.Y() * pathLength;
            vz -= dir.Z() * pathLength;
        }

        // Optional straight-line tail from POCA to fBackExtrapTargetX (in
        // input units). For wide-field setups where the production vertex
        // sits inside the field region, leave fBackExtrapTargetX at NaN.
        if (!std::isnan(fBackExtrapTargetX))
        {
            const double targetX_mm = fBackExtrapTargetX * fInputUnit_mm;
            if (vx > targetX_mm)
            {
                ROOT::Math::Polar3DVector momDir(1.0, theta_s, phi_s);
                ROOT::Math::XYZVector dir(momDir);
                if (std::abs(dir.X()) > 0.01)
                {
                    const double tailPath = (vx - targetX_mm) / dir.X();
                    vx -= dir.X() * tailPath;
                    vy -= dir.Y() * tailPath;
                    vz -= dir.Z() * tailPath;
                }
            }
        }

        // Energy correction along the back-extrapolated path (signed: we
        // add eLost since we're going *backward* in time, so the particle
        // had MORE energy at the vertex).
        double KE_at_cluster = std::sqrt(p_s * p_s + fMass_MeV * fMass_MeV) - fMass_MeV;
        if (auto* elossModel = fUKF->GetPropagator().GetELossModel())
        {
            const double dEdx = elossModel->GetdEdx(KE_at_cluster);
            const double eLost = dEdx * pathLength;
            const double KE_at_vertex = KE_at_cluster + eLost;
            if (KE_at_vertex > 0)
                p_s = std::sqrt(KE_at_vertex * KE_at_vertex + 2 * KE_at_vertex * fMass_MeV);
        }
    }

    const double KE_vtx = std::sqrt(p_s * p_s + fMass_MeV * fMass_MeV) - fMass_MeV;

    // Vertex + smoothed positions are stored in the input unit system so
    // they line up with the upstream R3BGTPCHitData / R3BGTPCTrackData.
    const double outScale = 1.0 / fInputUnit_mm;
    fitted->SetVertex(ROOT::Math::XYZVector(vx * outScale, vy * outScale, vz * outScale));
    fitted->SetKinematics(KE_vtx, theta_s, phi_s);

    std::vector<ROOT::Math::XYZPoint> smoothedPositions;
    smoothedPositions.reserve(smoothed.size());
    for (const auto& s : smoothed)
        smoothedPositions.emplace_back(s[0] * outScale, s[1] * outScale, s[2] * outScale);
    fitted->SetSmoothedPositions(std::move(smoothedPositions));

    // -- 6. chi^2 / ndf (rough: sum of (measurement - smoothed)^2 / sigma^2).
    //       Both terms are in the input unit system (mm cancellation done
    //       implicitly via outScale on smoothed; cluster positions stay in
    //       input units).
    double chi2 = 0;
    int ndf = 0;
    const double sig_in = fMeasSigma_mm * outScale;
    const double sig2 = sig_in * sig_in;
    for (size_t i = 0; i < clusters->size() && i < smoothed.size(); ++i)
    {
        const auto& ci = clusters->at(i);
        const double dx = ci.GetX() - smoothed[i][0] * outScale;
        const double dy = ci.GetY() - smoothed[i][1] * outScale;
        const double dz = ci.GetZ() - smoothed[i][2] * outScale;
        chi2 += (dx * dx + dy * dy + dz * dz) / sig2;
        ndf += 3;
    }
    ndf -= 6; // subtract free parameters
    fitted->SetChi2(chi2);
    fitted->SetNdf(std::max(ndf, 1));

    return fitted;
}
