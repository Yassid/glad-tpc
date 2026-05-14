/******************************************************************************
 *   Copyright (C) 2018-2026 GSI Helmholtzzentrum für Schwerionenforschung   *
 *   Copyright (C) 2018-2026 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 ******************************************************************************/

#pragma once

#include "TObject.h"
#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <limits>
#include <vector>

/// Output of R3BGTPCFitterUKF: post-UKF kinematics, vertex (back-extrapolated
/// to the closest approach to the beam axis), smoothed cluster positions, and
/// the standard chi^2/ndf quality metric. Two kinematics blocks are kept:
///   - fKinematicsXtr — at the first cluster (pre-back-extrap, measurement-
///     dominated; useful as a sanity check decoupled from the back-extrap
///     energy-loss correction).
///   - fKinematics    — at the back-extrapolated vertex (canonical: this is
///     what physics analyses should use).
class R3BGTPCFittedTrackData : public TObject
{
  public:
    struct Kinematics
    {
        Double_t kineticEnergy{ std::numeric_limits<Double_t>::quiet_NaN() };
        Double_t theta{ std::numeric_limits<Double_t>::quiet_NaN() };
        Double_t phi{ std::numeric_limits<Double_t>::quiet_NaN() };
    };

    R3BGTPCFittedTrackData() = default;
    ~R3BGTPCFittedTrackData() override = default;

    // -- setters
    void SetTrackId(Int_t id) { fTrackId = id; }
    void SetKinematics(const Kinematics& k) { fKinematics = k; }
    void SetKinematicsXtr(const Kinematics& k) { fKinematicsXtr = k; }
    void SetKinematics(Double_t ke, Double_t theta, Double_t phi) { fKinematics = { ke, theta, phi }; }
    void SetKinematicsXtr(Double_t ke, Double_t theta, Double_t phi) { fKinematicsXtr = { ke, theta, phi }; }
    void SetVertex(const ROOT::Math::XYZVector& v) { fVertex = v; }
    void SetSmoothedPositions(std::vector<ROOT::Math::XYZPoint> pos) { fSmoothedPositions = std::move(pos); }
    void SetChi2(Double_t chi2) { fChi2 = chi2; }
    void SetNdf(Int_t ndf) { fNdf = ndf; }
    void SetConverged(Bool_t conv) { fConverged = conv; }

    // -- getters
    Int_t GetTrackId() const { return fTrackId; }
    const Kinematics& GetKinematics() const { return fKinematics; }
    const Kinematics& GetKinematicsXtr() const { return fKinematicsXtr; }
    const ROOT::Math::XYZVector& GetVertex() const { return fVertex; }
    const std::vector<ROOT::Math::XYZPoint>& GetSmoothedPositions() const { return fSmoothedPositions; }
    Double_t GetChi2() const { return fChi2; }
    Int_t GetNdf() const { return fNdf; }
    Double_t GetChi2OverNdf() const { return fNdf > 0 ? fChi2 / fNdf : std::numeric_limits<Double_t>::quiet_NaN(); }
    Bool_t IsConverged() const { return fConverged; }

  protected:
    Int_t fTrackId{ -1 };
    Kinematics fKinematics{};    ///< at back-extrapolated vertex
    Kinematics fKinematicsXtr{}; ///< at first cluster (pre-back-extrap)
    ROOT::Math::XYZVector fVertex{ std::numeric_limits<Double_t>::quiet_NaN(),
                                   std::numeric_limits<Double_t>::quiet_NaN(),
                                   std::numeric_limits<Double_t>::quiet_NaN() };
    std::vector<ROOT::Math::XYZPoint> fSmoothedPositions{};
    Double_t fChi2{ std::numeric_limits<Double_t>::quiet_NaN() };
    Int_t fNdf{ 0 };
    Bool_t fConverged{ kFALSE };

    ClassDefOverride(R3BGTPCFittedTrackData, 1);
};
