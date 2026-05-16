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

#pragma once

#include "TClonesArray.h"
#include <Rtypes.h> // for THashConsistencyHolder, ClassDef

#include "cluster.h" // for Cluster
#include <stdio.h>   // for size_t

#include <memory> // for unique_ptr
#include <vector> // for vector

#include "R3BGTPCHitClusterData.h"
#include "R3BGTPCTrackData.h"

struct tc_params
{
    float s;
    size_t k;
    size_t n;
    size_t m;
    float r;
    float a;
    float t;
    float _padding;
};

class R3BGTPCTrackFinder
{
  private:
    tc_params inputParams{ .s = 0.3, .k = 19, .n = 2, .m = 15, .r = 2, .a = 0.03, .t = 4.0 };

  public:
    R3BGTPCTrackFinder();
    virtual ~R3BGTPCTrackFinder() = default;
    void Clusterize(R3BGTPCTrackData& track, Float_t distance, Float_t radius);

    /// Kasa LSQ circle fit on the track's xy hit positions + LSQ refit of
    /// z = a + b·phi on (phi, z) around the circle centre. Sets fGeoCenter,
    /// fGeoRadius, fGeoTheta on the track. Required by R3BGTPCFitterUKF's
    /// Brho seed and back-extrapolation.
    void SetTrackInitialParameters(R3BGTPCTrackData& track);

    void eventToClusters(TClonesArray* hitCA, PointCloud& cloud);
    std::unique_ptr<R3BGTPCTrackData> clustersToTrack(PointCloud& cloud,
                                                      const std::vector<cluster_t>& clusters,
                                                      TClonesArray* trackCA,
                                                      TClonesArray* hitCA);

    void SetScluster(float s) { inputParams.s = s; }
    void SetKtriplet(size_t k) { inputParams.k = k; }
    void SetNtriplet(size_t n) { inputParams.n = n; }
    void SetMcluster(size_t m) { inputParams.m = m; }
    void SetRsmooth(float r) { inputParams.r = r; }
    void SetAtriplet(float a) { inputParams.a = a; }
    void SetTcluster(float t) { inputParams.t = t; }

    /// Live tc_params view for R3BGTPCHit2Track::Exec — keeps the legacy
    /// triplclust Opt-defaults code path while letting macros override
    /// individual knobs via the SetXxx setters above.
    const tc_params& GetInputParams() const { return inputParams; }

    /// When true, R3BGTPCHit2Track multiplies r and s by the data's dNN
    /// (matches the upstream triplclust Opt defaults: r=2·dNN, s=0.3·dNN).
    /// When false, r and s are used literally — useful if the macro pins
    /// absolute mm-scale values.
    void SetUseDnnScaling(bool use) { fUseDnnScaling = use; }
    bool UseDnnScaling() const { return fUseDnnScaling; }

    /// Add a vertex constraint to the Gauss-Newton circle fit, in the (x, z)
    /// bending plane and in chamber-local cm. The constraint enters as one
    /// extra residual with effective weight = (σ_xy / fVertexSigma)². Set
    /// fVertexSigma <= 0 (or call SetVertexConstraint with neg) to disable.
    /// This is the highest-leverage knob on small chambers: extending the
    /// effective chord from ~9 cm (in-pad-plane) to ~16 cm (target-to-far
    /// wall) drops the Gluckstern floor by (L_new/L_old)² ≈ 3×.
    void SetVertexConstraint(double xL, double zL, double sigma_cm = 0.1)
    {
        fVertexX = xL;
        fVertexZ = zL;
        fVertexSigma = sigma_cm;
    }
    void DisableVertexConstraint() { fVertexSigma = -1.0; }

    /// Huber-weighted Gauss-Newton refit. When fHuberK_cm > 0, each
    /// residual r in the GN inner loop is multiplied by w(r) = 1 if
    /// |r| ≤ fHuberK_cm, else fHuberK_cm / |r|. This downweights outlier
    /// hits (δ-rays, misclustering) so they don't pull the circle. k≈1.5×σ_hit
    /// (i.e. ~1.5 mm for σ_xy=1 mm) gives 95% efficiency at the normal.
    void SetHuberK(double k_cm) { fHuberK_cm = k_cm; }

  private:
    bool fUseDnnScaling{ true };
    double fVertexX{ 0.0 };
    double fVertexZ{ 0.0 };
    double fVertexSigma{ -1.0 }; // <0 disables
    double fHuberK_cm{ -1.0 };    // <=0 disables (pure L2)

  public:

    ClassDef(R3BGTPCTrackFinder, 1);
};
