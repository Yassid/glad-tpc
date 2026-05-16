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

#include "R3BGTPCTrackFinder.h"

#include "TMath.h"
#include <Math/Point3D.h> // for PositionVector3D

#include <boost/core/checked_delete.hpp>  // for checked_delete
#include <boost/smart_ptr/shared_ptr.hpp> // for shared_ptr

#include <algorithm>
#include <cmath>    // for sqrt
#include <iostream> // for cout, cerr
#include <memory>   // for allocator_traits<>::value_...
#include <utility>  // for move

#include "dnn.h"
#include "graph.h"
#include "option.h"
#include "output.h"
#include "pointcloud.h"

#include <Math/Point3D.h>
#include <Math/Point3Dfwd.h>
#include <Math/Vector3D.h>
#include <Math/Vector3Dfwd.h>
#include <Rtypes.h>
#include <TObject.h>

constexpr auto cRED = "\033[1;31m";
constexpr auto cYELLOW = "\033[1;33m";
constexpr auto cNORMAL = "\033[0m";
constexpr auto cGREEN = "\033[1;32m";

R3BGTPCTrackFinder::R3BGTPCTrackFinder() {}

void R3BGTPCTrackFinder::eventToClusters(TClonesArray* hitCA, PointCloud& cloud)
{

    Int_t nHits = hitCA->GetEntries();
    R3BGTPCHitData** hitData;
    hitData = new R3BGTPCHitData*[nHits];

    if (hitData)
    {

        for (Int_t iHit = 0; iHit < nHits; iHit++)
        {
            Point point;
            hitData[iHit] = (R3BGTPCHitData*)(hitCA->At(iHit));
            point.x = hitData[iHit]->GetX();
            point.y = hitData[iHit]->GetY();
            point.z = hitData[iHit]->GetZ();
            point.SetID(iHit);
            cloud.push_back(point);
        }

        delete hitData;
    }
}

std::unique_ptr<R3BGTPCTrackData> R3BGTPCTrackFinder::clustersToTrack(PointCloud& cloud,
                                                                      const std::vector<cluster_t>& clusters,
                                                                      TClonesArray* trackCA,
                                                                      TClonesArray* hitCA)
{

    std::vector<R3BGTPCTrackData> tracks;
    std::vector<Point> points = cloud;

    for (size_t cluster_index = 0; cluster_index < clusters.size(); ++cluster_index)
    {

        R3BGTPCTrackData track; // One track per cluster

        const std::vector<size_t>& point_indices = clusters[cluster_index];
        if (point_indices.size() == 0)
            continue;

        // add points
        for (std::vector<size_t>::const_iterator it = point_indices.begin(); it != point_indices.end(); ++it)
        {

            const Point& point = cloud[*it];

            Int_t nHits = hitCA->GetEntries();

            R3BGTPCHitData** hitData;
            hitData = new R3BGTPCHitData*[nHits];

            if (hitData)
            {

                hitData[point.GetID()] = (R3BGTPCHitData*)(hitCA->At(point.GetID()));
                track.AddHit(*hitData[point.GetID()]);
                delete hitData;
            }

            // remove current point from vector points
            for (std::vector<Point>::iterator p = points.begin(); p != points.end(); p++)
            {
                if (*p == point)
                {
                    points.erase(p);
                    break;
                }
            }

        } // Point indices

        track.SetTrackId(cluster_index);
        Clusterize(track, 0.70, 1.5);

        tracks.push_back(track);

    } // Clusters loop

    std::cout << cRED << " Tracks found " << tracks.size() << cNORMAL << "\n";

    // TODO
    // Dump noise into pattern event
    // auto retEvent = std::make_unique<AtPatternEvent>();
    // for (const auto &point : points)
    // retEvent->AddNoise(event.GetHit(point.intensity));

    for (auto& track : tracks)
    {
        if (track.GetHitArray().size() >= 3)
            SetTrackInitialParameters(track);

        TClonesArray& clref = *trackCA;
        Int_t size = clref.GetEntriesFast();
        auto* persisted = new (clref[size]) R3BGTPCTrackData(
            track.GetTrackId(), std::move(track.GetHitArray()), std::move(*track.GetHitClusterArray()));
        persisted->SetGeoCenter(track.GetGeoCenter());
        persisted->SetGeoRadius(track.GetGeoRadius());
        persisted->SetGeoTheta(track.GetGeoTheta());
    }

    return NULL;
}

void R3BGTPCTrackFinder::Clusterize(R3BGTPCTrackData& track, Float_t distance, Float_t radius)
{

    std::vector<R3BGTPCHitData> hitArray = track.GetHitArray();
    std::vector<R3BGTPCHitData> hitTBArray;
    int clusterID = 0;

    std::cout << " Number of hits per track : " << hitArray.size() << "\n";

    if (hitArray.size() > 0)
    {

        ROOT::Math::XYZVector refPos{ hitArray.at(0).GetX(), hitArray.at(0).GetY(), hitArray.at(0).GetZ() };

        for (auto iHit = 0; iHit < hitArray.size(); ++iHit)
        {

            auto hit = hitArray.at(iHit);
            ROOT::Math::XYZVector hitPos{ hit.GetX(), hit.GetY(), hit.GetZ() };

            // Check distance with respect to reference Hit
            Double_t distRef = TMath::Sqrt((hitPos - refPos).Mag2());

            if (distRef < distance)
            {

                continue;
            }
            else
            {

                Double_t clusterQ = 0.0;
                hitTBArray.clear();
                std::copy_if(hitArray.begin(),
                             hitArray.end(),
                             std::back_inserter(hitTBArray),
                             [&refPos, radius](R3BGTPCHitData& hitIn)
                             {
                                 ROOT::Math::XYZVector hitInPos{ hitIn.GetX(), hitIn.GetY(), hitIn.GetZ() };
                                 return TMath::Sqrt((hitInPos - refPos).Mag2()) < radius;
                             });

                if (hitTBArray.size() > 0)
                {
                    double x = 0, y = 0, z = 0;

                    int timeStamp = 0;
                    std::shared_ptr<R3BGTPCHitClusterData> hitCluster = std::make_shared<R3BGTPCHitClusterData>();
                    hitCluster->SetClusterID(clusterID);
                    Double_t hitQ = 0.0;
                    std::for_each(hitTBArray.begin(),
                                  hitTBArray.end(),
                                  [&x, &y, &z, &hitQ, &timeStamp](R3BGTPCHitData& hitInQ)
                                  {
                                      ROOT::Math::XYZPoint pos{ hitInQ.GetX(), hitInQ.GetY(), hitInQ.GetZ() };
                                      x += pos.X() * hitInQ.GetEnergy();
                                      y += pos.Y() * hitInQ.GetEnergy();
                                      z += pos.Z();
                                      hitQ += hitInQ.GetEnergy();
                                      // TODO
                                      // timeStamp += hitInQ.GetTimeStamp();
                                  });
                    x /= hitQ;
                    y /= hitQ;
                    z /= hitTBArray.size();
                    // timeStamp /= hitTBArray.size();

                    ROOT::Math::XYZPoint clustPos(x, y, z);
                    Bool_t checkDistance = kTRUE;

                    // Check distance with respect to existing clusters
                    for (auto iClusterHit : *track.GetHitClusterArray())
                    {
                        ROOT::Math::XYZPoint iclusterHitPos{ iClusterHit.GetX(),
                                                             iClusterHit.GetY(),
                                                             iClusterHit.GetZ() };
                        if (TMath::Sqrt((iclusterHitPos - clustPos).Mag2()) < distance)
                        {
                            // std::cout<<" Cluster with less than  : "<<distance<<" found
                            // "<<"\n";
                            checkDistance = kFALSE;
                            continue;
                        }
                    }

                    if (checkDistance)
                    {

                        hitCluster->SetEnergy(hitQ);
                        hitCluster->SetX(x);
                        hitCluster->SetY(y);
                        hitCluster->SetZ(z);
                        // hitCluster->SetTime(timeStamp);
                        ++clusterID;
                        track.AddClusterHit(hitCluster);
                    }

                } // if hitTBArray>0

            } // if distance

            ROOT::Math::XYZVector refPosBuff{ hitArray.at(iHit).GetX(),
                                              hitArray.at(iHit).GetY(),
                                              hitArray.at(iHit).GetZ() };

            refPos = refPosBuff;

        } // for Hit array

    } // if array size
}

void R3BGTPCTrackFinder::SetTrackInitialParameters(R3BGTPCTrackData& track)
{
    // Circle fit in the bending plane (perpendicular to B). R3B GLAD has
    // B = (0, B_y, 0) so the bending plane is (x, z). We use Pratt's
    // method (N. Chernov, "Circular and Linear Regression", Ch. 5) which
    // is a closed-form algebraic fit unbiased to leading order in the
    // hit noise — substantially better than naive Kasa, which biases R
    // low when sagitta is comparable to noise. A short Gauss-Newton
    // refinement on the perpendicular residual then gives the geometric
    // optimum.
    //
    // GeoCenter stores (cx, cz); GeoTheta stores the angle from +ŷ.
    const auto& hits = track.GetHitArray();
    const std::size_t n = hits.size();
    if (n < 3)
        return;

    // Pratt fit (centred data, Newton on the characteristic cubic).
    double xb = 0, zb = 0;
    for (const auto& h : hits) { xb += h.GetX(); zb += h.GetZ(); }
    xb /= n;
    zb /= n;
    double Mxx = 0, Myy = 0, Mxy = 0, Mxz = 0, Myz = 0, Mzz = 0;
    for (const auto& h : hits)
    {
        const double dx = h.GetX() - xb;
        const double dy = h.GetZ() - zb;
        const double z = dx * dx + dy * dy;
        Mxx += dx * dx;
        Myy += dy * dy;
        Mxy += dx * dy;
        Mxz += dx * z;
        Myz += dy * z;
        Mzz += z * z;
    }
    const double nd = static_cast<double>(n);
    Mxx /= nd; Myy /= nd; Mxy /= nd; Mxz /= nd; Myz /= nd; Mzz /= nd;
    const double Mz = Mxx + Myy;
    const double Cov_xy = Mxx * Myy - Mxy * Mxy;
    const double Mxz2 = Mxz * Mxz;
    const double Myz2 = Myz * Myz;
    // Cubic A3·t³ + A2·t² + A1·t + A0 = 0 (Pratt's characteristic equation)
    const double A3 = 4.0 * Mz;
    const double A2 = -3.0 * Mz * Mz - Mzz;
    const double A1 = Mzz * Mz + 4.0 * Cov_xy * Mz - Mxz2 - Myz2 - Mz * Mz * Mz;
    const double A0 = Mxz2 * Myy + Myz2 * Mxx - Mzz * Cov_xy - 2.0 * Mxz * Myz * Mxy + Mz * Mz * Cov_xy;
    const double A22 = A2 + A2;
    const double A33 = A3 + A3 + A3;
    // Newton's method on the cubic, starting at t=0 (smallest root)
    double tnew = 0, ynew = 1e20;
    for (int iter = 0; iter < 100; ++iter)
    {
        const double yold = ynew;
        ynew = A0 + tnew * (A1 + tnew * (A2 + tnew * A3));
        if (std::abs(ynew) > std::abs(yold)) break;
        const double Dy = A1 + tnew * (A22 + tnew * A33);
        if (std::abs(Dy) < 1e-30) break;
        const double told = tnew;
        tnew = told - ynew / Dy;
        if (std::abs(tnew - told) < 1e-12 * std::abs(tnew == 0 ? 1.0 : tnew)) break;
        if (tnew < 0) { tnew = 0; break; }
    }
    const double DET = tnew * tnew - tnew * Mz + Cov_xy;
    if (std::abs(DET) < 1e-30) return;
    const double Xcenter = (Mxz * (Myy - tnew) - Myz * Mxy) / (DET * 2.0);
    const double Ycenter = (Myz * (Mxx - tnew) - Mxz * Mxy) / (DET * 2.0);
    double cx = Xcenter + xb;
    double cz = Ycenter + zb;
    const double R2_pratt = Xcenter * Xcenter + Ycenter * Ycenter + Mz + 2.0 * tnew;
    if (!(R2_pratt > 0)) return;
    double R = std::sqrt(R2_pratt);

    // Geometric refinement (Gauss-Newton on Σ(d − R)²) from the Pratt seed.
    // Pratt is already unbiased to leading order; this nails the geometric
    // optimum in 2-5 iterations.
    //
    // Vertex constraint: when fVertexSigma > 0, append one extra residual
    // for a virtual point at (fVertexX, fVertexZ) with effective weight
    // w = (sigma_hit / fVertexSigma)² ≈ 1 (we assume sigma_hit ~ 1 mm and
    // sigma_vertex ~ 1 mm by default). The lever arm is the geometric
    // distance from the vertex to the centroid of the in-pad hits.
    const double w_vtx = (fVertexSigma > 0) ? (0.1 / fVertexSigma) * (0.1 / fVertexSigma) : 0.0;
    const bool useHuber = (fHuberK_cm > 0.0);
    for (int iter = 0; iter < 30; ++iter)
    {
        double Jxx = 0, Jxy = 0, Jxr = 0, Jyy = 0, Jyr = 0, Jrr = 0;
        double bx = 0, by = 0, br = 0;
        for (const auto& h : hits)
        {
            const double dx = h.GetX() - cx;
            const double dy = h.GetZ() - cz;
            const double d = std::sqrt(dx * dx + dy * dy);
            if (d < 1e-9) continue;
            const double r = d - R;
            // Huber re-weighting: w(r) = 1 if |r| ≤ k, else k/|r|. Effectively
            // shrinks far-from-circle hits (δ-rays, misclusterings) so they
            // don't pull the geometric optimum.
            const double w = (useHuber && std::abs(r) > fHuberK_cm) ? (fHuberK_cm / std::abs(r)) : 1.0;
            const double Jx = -dx / d, Jy = -dy / d, JR = -1.0;
            Jxx += w * Jx * Jx; Jxy += w * Jx * Jy; Jxr += w * Jx * JR;
            Jyy += w * Jy * Jy; Jyr += w * Jy * JR; Jrr += w * JR * JR;
            bx  += w * Jx * r;  by  += w * Jy * r;  br  += w * JR * r;
        }
        if (w_vtx > 0.0)
        {
            const double dx = fVertexX - cx;
            const double dy = fVertexZ - cz;
            const double d = std::sqrt(dx * dx + dy * dy);
            if (d > 1e-9)
            {
                const double r = d - R;
                const double Jx = -dx / d, Jy = -dy / d, JR = -1.0;
                Jxx += w_vtx * Jx * Jx; Jxy += w_vtx * Jx * Jy; Jxr += w_vtx * Jx * JR;
                Jyy += w_vtx * Jy * Jy; Jyr += w_vtx * Jy * JR; Jrr += w_vtx * JR * JR;
                bx  += w_vtx * Jx * r;  by  += w_vtx * Jy * r;  br  += w_vtx * JR * r;
            }
        }
        const double dt = Jxx * (Jyy * Jrr - Jyr * Jyr) - Jxy * (Jxy * Jrr - Jyr * Jxr)
                          + Jxr * (Jxy * Jyr - Jyy * Jxr);
        if (std::abs(dt) < 1e-12) break;
        const double dxc = -(bx * (Jyy * Jrr - Jyr * Jyr) - Jxy * (by * Jrr - Jyr * br)
                             + Jxr * (by * Jyr - Jyy * br)) / dt;
        const double dyc = -(Jxx * (by * Jrr - Jyr * br) - bx * (Jxy * Jrr - Jyr * Jxr)
                             + Jxr * (Jxy * br - by * Jxr)) / dt;
        const double dRc = -(Jxx * (Jyy * br - by * Jyr) - Jxy * (Jxy * br - by * Jxr)
                             + bx * (Jxy * Jyr - Jyy * Jxr)) / dt;
        cx += dxc;
        cz += dyc;
        R += dRc;
        if (std::abs(dxc) + std::abs(dyc) + std::abs(dRc) < 1e-6) break;
    }
    if (!(R > 0) || !std::isfinite(R))
        return;
    track.SetGeoCenter({ cx, cz });
    track.SetGeoRadius(R);

    // LSQ refit of y = a + b·phi on (phi, y), with phi computed around the
    // fitted circle centre in (x,z). The angle from the B direction satisfies
    //     cot(θ_y) = sign(ψ̇)·b/R,
    // where sign(ψ̇) is the rotation direction along the track in (x,z).
    if (n < 2)
    {
        track.SetGeoTheta(TMath::Pi() / 2.0);
        return;
    }
    double sphi = 0., sy_p = 0., sphi2 = 0., sphiy = 0.;
    for (const auto& h : hits)
    {
        const double phi = std::atan2(h.GetZ() - cz, h.GetX() - cx);
        const double y = h.GetY();
        sphi += phi;
        sy_p += y;
        sphi2 += phi * phi;
        sphiy += phi * y;
    }
    const double detPhi = nd * sphi2 - sphi * sphi;
    if (std::abs(detPhi) < 1e-9)
    {
        track.SetGeoTheta(TMath::Pi() / 2.0);
        return;
    }
    const double b = (nd * sphiy - sphi * sy_p) / detPhi;

    // Resolve rotation sign from the φ-walk in (x,z) between first and last
    // hit (pattern recognition returns hits in track order). Unwrap dφ to
    // (-π, π] so a single half-turn track keeps its direction.
    const double phi_first = std::atan2(hits.front().GetZ() - cz, hits.front().GetX() - cx);
    const double phi_last = std::atan2(hits.back().GetZ() - cz, hits.back().GetX() - cx);
    double dphi = phi_last - phi_first;
    while (dphi > TMath::Pi())
        dphi -= 2.0 * TMath::Pi();
    while (dphi <= -TMath::Pi())
        dphi += 2.0 * TMath::Pi();
    const double sign_psi_dot = (dphi >= 0) ? 1.0 : -1.0;

    const double cot_th = b * sign_psi_dot / R;
    track.SetGeoTheta(std::atan2(1.0, cot_th));
}
