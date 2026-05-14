/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#pragma once

#include "R3BGTPCHitClusterData.h"
#include "R3BGTPCHitData.h"
#include "TObject.h"
#include <limits>
#include <memory>
#include <stdint.h>
#include <utility>

class R3BGTPCTrackData : public TObject
{

  public:
    // Default Constructor
    R3BGTPCTrackData();

    /** Standard Constructor

    **/
    R3BGTPCTrackData(std::size_t trackId,
                     std::vector<R3BGTPCHitData> hitArray,
                     std::vector<R3BGTPCHitClusterData> hitClusterArray);

    // Destructor
    virtual ~R3BGTPCTrackData() {}

    // Getters
    Int_t GetTrackId() { return fTrackId; }
    std::vector<R3BGTPCHitData>& GetHitArray() { return fHitArray; }
    std::vector<R3BGTPCHitClusterData>* GetHitClusterArray() { return &fHitClusterArray; }

    // Geometry from the pattern-recognition circle/helix fit (Riemann/RANSAC,
    // Kasa, ...). Populated by the track-finder; consumed by the Kalman-style
    // fitters (e.g. R3BGTPCFitterUKF) for the helix-POCA back-extrapolation
    // and the initial-momentum seed. Defaults are NaN so a missing population
    // is obvious rather than silently propagating zero.
    std::pair<Double_t, Double_t> GetGeoCenter() const { return fGeoCenter; }
    Double_t GetGeoRadius() const { return fGeoRadius; }
    Double_t GetGeoTheta() const { return fGeoTheta; }

    // Setters
    void SetTrackId(Int_t val) { fTrackId = val; }
    void AddHit(R3BGTPCHitData& hit) { fHitArray.push_back(hit); }
    void AddClusterHit(std::shared_ptr<R3BGTPCHitClusterData> hitCluster)
    {
        fHitClusterArray.push_back(std::move(*hitCluster));
    }
    void SetGeoCenter(std::pair<Double_t, Double_t> center) { fGeoCenter = center; }
    void SetGeoRadius(Double_t radius) { fGeoRadius = radius; }
    void SetGeoTheta(Double_t theta) { fGeoTheta = theta; }

  protected:
    Int_t fTrackId{ -1 };                  // Track Id
    std::vector<R3BGTPCHitData> fHitArray; // Track Hit Array
    std::vector<R3BGTPCHitClusterData> fHitClusterArray;

    // Pattern-recognition geometry (set by track-finder; NaN until populated)
    std::pair<Double_t, Double_t> fGeoCenter{ std::numeric_limits<Double_t>::quiet_NaN(),
                                              std::numeric_limits<Double_t>::quiet_NaN() };
    Double_t fGeoRadius{ std::numeric_limits<Double_t>::quiet_NaN() };
    Double_t fGeoTheta{ std::numeric_limits<Double_t>::quiet_NaN() };

    ClassDef(R3BGTPCTrackData, 2)
};
