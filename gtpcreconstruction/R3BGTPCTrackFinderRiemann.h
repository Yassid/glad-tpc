/******************************************************************************
 * Copyright (C) 2018-2026 GSI Helmholtzzentrum für Schwerionenforschung GmbH *
 *         Copyright (C) 2018-2026 Members of R3B Collaboration               *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *              GNU Lesser General Public Licence (LGPL) version 3,           *
 *                     copied verbatim in the file "LICENSE".                 *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

/********************************************************************
 * R3BGTPCTrackFinderRiemann                                          *
 *                                                                    *
 * Smart-seed RANSAC pattern recognition for the GLAD-TPC, ported     *
 * verbatim from ATTPCROOT's AtTrackFinderRiemann. Operates on a flat *
 * std::vector<R3BGTPCHitData> and returns a vector of                *
 * R3BGTPCTrackData with the inlier hit set and the refined circle    *
 * geometry (GeoCenter, GeoRadius, GeoTheta) populated so the         *
 * downstream UKF stage can read its seed straight from the track.    *
 *                                                                    *
 * Algorithm (per RANSAC pass):                                       *
 *   - seed1 = random unassigned hit;                                 *
 *   - seed2/seed3 = picked from seed1's k nearest neighbors in 3D;   *
 *   - analytic circle through 3 seeds + z = a + b·phi line anchored  *
 *     on the two seeds with the largest |Δphi|;                      *
 *   - 3D inlier cut: |xy circle residual| < fInlierDist AND          *
 *                    |z - (a + b·phi)|     < fZInlierDist;           *
 *   - Kasa LSQ refit + LSQ z-line refit + 3D re-collect;             *
 *   - largest-contiguous-arc filter + (optional) sliding-window arc  *
 *     walk to recover track tails through small gaps.                *
 *                                                                    *
 * The bending plane in R3B GLAD is (x, z) since B = (0, B_y, 0); the *
 * adapter below feeds (x, z) to the algorithm as the 2D plane and    *
 * uses y as the "z" axis (the helix axis along B).                   *
 ********************************************************************/

#pragma once

#include "R3BGTPCHitData.h"
#include "R3BGTPCTrackData.h"

#include <Rtypes.h>

#include <vector>

class R3BGTPCTrackFinderRiemann
{
  public:
    R3BGTPCTrackFinderRiemann() = default;
    virtual ~R3BGTPCTrackFinderRiemann() = default;

    /// Run RANSAC on the full hit list and return the assembled tracks.
    /// Each returned track has its hit array filled with the inliers and
    /// GeoCenter/GeoRadius/GeoTheta set from the refined circle.
    std::vector<R3BGTPCTrackData> FindTracks(const std::vector<R3BGTPCHitData>& hits);

    // RANSAC + inlier knobs
    void SetInlierDist(double d) { fInlierDist = d; }
    void SetZInlierDist(double d) { fZInlierDist = d; }
    void SetMinHitsPerTrack(int n) { fMinHitsPerTrack = n; }
    void SetMaxTracks(int n) { fMaxTracks = n; }
    void SetMaxIterations(int n) { fMaxIterations = n; }
    void SetSeed(unsigned int s) { fSeed = s; }
    void SetKNN(int k) { fK_NN = k; }
    void SetMaxPhiGap(double deg) { fMaxPhiGap = deg; }
    void SetUseArcWalkExtend(bool use) { fUseArcWalkExtend = use; }
    void SetArcWalkWindow(int n) { fArcWalkWindow = n; }
    void SetArcWalkMaxMiss(int n) { fArcWalkMaxMiss = n; }

    /// Minimum mutual distance between the 3 RANSAC seed hits, in cm. Required
    /// so the seed triple spans enough chord to anchor R against hit noise.
    /// Set to 0 to disable. Useful for small chambers (HYDRA Prototype) where
    /// the kNN-locality heuristic from AT-TPC over-clusters the seeds.
    void SetMinSeedSpread(double cm) { fMinSeedSpread = cm; }

    /// Minimum candidate circle radius, in cm. Reject the seed if its analytic
    /// circle is tighter than this — tiny noise-circles through 3 nearby
    /// points are the main RANSAC failure mode on near-straight tracks.
    void SetMinCircleR(double cm) { fMinCircleR = cm; }

  private:
    double fInlierDist{ 0.5 };       ///< Max |(x,z) circle residual| in cm (~5 mm)
    double fZInlierDist{ 1.5 };      ///< Max |y - (a + b·phi)| in cm (~15 mm along helix)
    int fMinHitsPerTrack{ 10 };      ///< Reject candidates below this
    int fMaxTracks{ 6 };             ///< Cap iterative extraction
    int fMaxIterations{ 400 };       ///< RANSAC iterations per pass
    unsigned int fSeed{ 12345 };     ///< RNG seed
    int fK_NN{ 10 };                 ///< kNN seeds per RANSAC iteration
    double fMaxPhiGap{ 60.0 };       ///< Max φ gap (deg) before arc split
    bool fUseArcWalkExtend{ true };  ///< Sliding-window tail-recovery walk
    int fArcWalkWindow{ 10 };
    int fArcWalkMaxMiss{ 5 };
    double fMinSeedSpread{ 0.0 };    ///< Min seed-triple mutual distance (cm); 0 = off
    double fMinCircleR{ 0.0 };       ///< Min candidate circle radius (cm); 0 = off
};
