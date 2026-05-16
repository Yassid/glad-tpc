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

#include "R3BGTPCTrackFinderRiemann.h"

#include <TMath.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <set>
#include <utility>
#include <vector>

namespace {

constexpr auto cRED = "\033[1;31m";
constexpr auto cNORMAL = "\033[0m";

struct Circle
{
    double cx{ 0 }, cy{ 0 }, r{ 0 };
    bool valid{ false };
};

struct ZLine
{
    double a{ 0 }, b{ 0 };
    bool valid{ false };
};

/// Analytic circle through 3 (x,y) points.
Circle circleThroughThree(const std::array<double, 2>& a,
                          const std::array<double, 2>& b,
                          const std::array<double, 2>& c)
{
    const double ax = a[0], ay = a[1];
    const double bx = b[0], by = b[1];
    const double cx = c[0], cy = c[1];

    const double d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if (std::abs(d) < 1e-9)
        return {};

    const double a2 = ax * ax + ay * ay;
    const double b2 = bx * bx + by * by;
    const double c2 = cx * cx + cy * cy;
    const double ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d;
    const double uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d;
    const double r = std::hypot(ax - ux, ay - uy);

    // Centimetre units: reject sub-mm radii and anything beyond 1 m of curvature.
    if (!(r > 0.05) || !(r < 1.0e3))
        return {};

    return { ux, uy, r, true };
}

inline double pointCircleDist(const Circle& c, double x, double y)
{
    return std::abs(std::hypot(x - c.cx, y - c.cy) - c.r);
}

/// y = a + b*phi line through the two seeds with the largest |Δphi|. (In R3B
/// GLAD the helix-axis coordinate is y, not z; the algorithm is otherwise the
/// same as the ATTPCROOT (x,y) plane + z-along-B case.)
ZLine zLineFromSeeds(double phi0, double y0, double phi1, double y1, double phi2, double y2)
{
    const double dp01 = std::abs(phi1 - phi0);
    const double dp02 = std::abs(phi2 - phi0);
    const double dp12 = std::abs(phi2 - phi1);
    double pa, ya, pb, yb, dp;
    if (dp01 >= dp02 && dp01 >= dp12)
    {
        pa = phi0; ya = y0; pb = phi1; yb = y1; dp = dp01;
    }
    else if (dp02 >= dp12)
    {
        pa = phi0; ya = y0; pb = phi2; yb = y2; dp = dp02;
    }
    else
    {
        pa = phi1; ya = y1; pb = phi2; yb = y2; dp = dp12;
    }
    if (dp < 0.05)
        return {};
    const double b = (yb - ya) / (pb - pa);
    const double a = ya - b * pa;
    return { a, b, true };
}

/// LSQ refit of y = a + b*phi on a hit set, with phi around (cx, cy).
ZLine zLineLSQ(const std::vector<std::array<double, 2>>& xy,
               const std::vector<double>& yAxis,
               const std::vector<std::size_t>& idx,
               const Circle& c)
{
    if (idx.size() < 2)
        return {};
    double sphi = 0, sy = 0, sphi2 = 0, sphiy = 0;
    const double n = static_cast<double>(idx.size());
    for (std::size_t i : idx)
    {
        const double phi = std::atan2(xy[i][1] - c.cy, xy[i][0] - c.cx);
        sphi += phi;
        sy += yAxis[i];
        sphi2 += phi * phi;
        sphiy += phi * yAxis[i];
    }
    const double det = n * sphi2 - sphi * sphi;
    if (std::abs(det) < 1e-9)
        return { 0.0, 0.0, false };
    const double b = (n * sphiy - sphi * sy) / det;
    const double a = (sy - b * sphi) / n;
    return { a, b, true };
}

/// Kasa linear-LSQ circle refit on selected indices (Cramer's rule on 3x3 normals).
Circle kasaRefit(const std::vector<std::array<double, 2>>& xy, const std::vector<std::size_t>& idx)
{
    if (idx.size() < 3)
        return {};
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, sxr = 0, syr = 0, sr = 0;
    for (std::size_t k : idx)
    {
        const double x = xy[k][0], y = xy[k][1];
        const double r2 = x * x + y * y;
        sx += x; sy += y;
        sxx += x * x; syy += y * y; sxy += x * y;
        sxr += x * r2; syr += y * r2; sr += r2;
    }
    const double n = static_cast<double>(idx.size());
    const double det = sxx * (syy * n - sy * sy) - sxy * (sxy * n - sy * sx) + sx * (sxy * sy - syy * sx);
    if (std::abs(det) < 1e-9)
        return {};
    const double dA = sxr * (syy * n - sy * sy) - sxy * (syr * n - sy * sr) + sx * (syr * sy - syy * sr);
    const double dB = sxx * (syr * n - sy * sr) - sxr * (sxy * n - sy * sx) + sx * (sxy * sr - syr * sx);
    const double dC = sxx * (syy * sr - syr * sy) - sxy * (sxy * sr - syr * sx) + sxr * (sxy * sy - syy * sx);
    const double A = dA / det, B = dB / det, C = dC / det;
    const double cx = A / 2.0, cy = B / 2.0;
    const double r2 = C + cx * cx + cy * cy;
    if (!(r2 > 0))
        return {};
    return { cx, cy, std::sqrt(r2), true };
}

/// Arc-walk tail recovery — sliding-window local Kasa refit walking outward
/// in φ from both ends of the inlier set. See the ATTPCROOT class comment.
std::vector<std::size_t> arcWalkExtend(const std::vector<std::array<double, 2>>& xy,
                                       const std::vector<double>& yVals,
                                       const std::vector<std::size_t>& inliers,
                                       const std::vector<std::size_t>& liveAll,
                                       const Circle& globalCircle,
                                       double xyTol,
                                       double yTol,
                                       int windowSize = 10,
                                       int maxMiss = 5)
{
    if (inliers.size() < static_cast<std::size_t>(std::max(3, windowSize / 2)))
        return inliers;

    auto phiOf = [&](std::size_t i)
    { return std::atan2(xy[i][1] - globalCircle.cy, xy[i][0] - globalCircle.cx); };

    std::set<std::size_t> inlierSet(inliers.begin(), inliers.end());
    std::vector<std::pair<double, std::size_t>> phiInl;
    phiInl.reserve(inliers.size());
    for (std::size_t i : inliers)
        phiInl.emplace_back(phiOf(i), i);
    std::sort(phiInl.begin(), phiInl.end());

    std::vector<std::pair<double, std::size_t>> phiCand;
    phiCand.reserve(liveAll.size());
    for (std::size_t i : liveAll)
    {
        if (inlierSet.count(i))
            continue;
        phiCand.emplace_back(phiOf(i), i);
    }
    std::sort(phiCand.begin(), phiCand.end());

    auto checkAccept = [&](const std::vector<std::size_t>& win, std::size_t cand) -> bool
    {
        const Circle lc = kasaRefit(xy, win);
        if (!lc.valid)
            return false;
        if (pointCircleDist(lc, xy[cand][0], xy[cand][1]) >= xyTol)
            return false;
        const ZLine zl = zLineLSQ(xy, yVals, win, lc);
        if (zl.valid)
        {
            const double phi = std::atan2(xy[cand][1] - lc.cy, xy[cand][0] - lc.cx);
            if (std::abs(yVals[cand] - (zl.a + zl.b * phi)) >= yTol)
                return false;
        }
        return true;
    };

    // Forward walk (+phi).
    {
        auto it = std::upper_bound(phiCand.begin(), phiCand.end(), phiInl.back(),
                                   [](const std::pair<double, std::size_t>& a,
                                      const std::pair<double, std::size_t>& b)
                                   { return a.first < b.first; });
        int miss = 0;
        for (; it != phiCand.end() && miss <= maxMiss; ++it)
        {
            const std::size_t k = std::min(static_cast<std::size_t>(windowSize), phiInl.size());
            std::vector<std::size_t> win;
            win.reserve(k);
            for (auto rit = phiInl.end() - k; rit != phiInl.end(); ++rit)
                win.push_back(rit->second);
            if (checkAccept(win, it->second))
            {
                phiInl.emplace_back(it->first, it->second);
                inlierSet.insert(it->second);
                miss = 0;
            }
            else
            {
                ++miss;
            }
        }
    }

    // Backward walk (−phi).
    {
        auto it = std::lower_bound(phiCand.begin(), phiCand.end(), phiInl.front(),
                                   [](const std::pair<double, std::size_t>& a,
                                      const std::pair<double, std::size_t>& b)
                                   { return a.first < b.first; });
        int miss = 0;
        while (it != phiCand.begin() && miss <= maxMiss)
        {
            --it;
            const std::size_t k = std::min(static_cast<std::size_t>(windowSize), phiInl.size());
            std::vector<std::size_t> win;
            win.reserve(k);
            for (auto fit = phiInl.begin(); fit != phiInl.begin() + k; ++fit)
                win.push_back(fit->second);
            if (checkAccept(win, it->second))
            {
                phiInl.insert(phiInl.begin(), { it->first, it->second });
                inlierSet.insert(it->second);
                miss = 0;
            }
            else
            {
                ++miss;
            }
        }
    }

    std::vector<std::size_t> out;
    out.reserve(phiInl.size());
    for (auto& p : phiInl)
        out.push_back(p.second);
    return out;
}

/// Largest contiguous arc: consecutive sorted phi gaps must stay below maxGap.
std::vector<std::size_t> largestContiguousArc(const std::vector<std::array<double, 2>>& xy,
                                              const std::vector<std::size_t>& inliers,
                                              const Circle& c,
                                              double maxGapDeg)
{
    if (inliers.size() < 2)
        return inliers;

    const double maxGap = maxGapDeg * M_PI / 180.0;
    std::vector<std::pair<double, std::size_t>> phiHits;
    phiHits.reserve(inliers.size());
    for (std::size_t i : inliers)
    {
        const double phi = std::atan2(xy[i][1] - c.cy, xy[i][0] - c.cx);
        phiHits.emplace_back(phi, i);
    }
    std::sort(phiHits.begin(), phiHits.end());

    const std::size_t n = phiHits.size();
    std::size_t gapStart = 0;
    double maxObservedGap = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        const double next = (i + 1 < n) ? phiHits[i + 1].first : (phiHits[0].first + 2 * M_PI);
        const double g = next - phiHits[i].first;
        if (g > maxObservedGap)
        {
            maxObservedGap = g;
            gapStart = i;
        }
    }
    if (maxObservedGap <= maxGap)
    {
        std::vector<std::size_t> out;
        out.reserve(inliers.size());
        for (auto& p : phiHits)
            out.push_back(p.second);
        return out;
    }

    std::vector<std::size_t> bestRun;
    for (std::size_t s = 0; s < n; ++s)
    {
        const std::size_t start = (gapStart + 1 + s) % n;
        std::vector<std::size_t> run{ phiHits[start].second };
        for (std::size_t step = 1; step < n; ++step)
        {
            const std::size_t cur = (start + step) % n;
            const std::size_t prev = (start + step - 1) % n;
            double g = phiHits[cur].first - phiHits[prev].first;
            if (cur < prev)
                g += 2 * M_PI;
            if (g > maxGap)
                break;
            run.push_back(phiHits[cur].second);
        }
        if (run.size() > bestRun.size())
            bestRun = std::move(run);
    }
    return bestRun;
}

/// Pratt's algebraic circle fit + Gauss-Newton geometric refinement.
/// Project memory documents this as the cure for Kasa's noise-induced low-R
/// bias on near-straight tracks (Chernov, "Circular and Linear Regression",
/// Ch. 5). Run as a FINAL refit on the chosen inlier set; the cheaper
/// algebraic Kasa is still fine for the inlier-classifier step.
Circle prattGNRefit(const std::vector<std::array<double, 2>>& xy, const std::vector<std::size_t>& idx)
{
    if (idx.size() < 3)
        return {};
    double xb = 0, zb = 0;
    for (std::size_t k : idx) { xb += xy[k][0]; zb += xy[k][1]; }
    const double nd = static_cast<double>(idx.size());
    xb /= nd; zb /= nd;
    double Mxx = 0, Myy = 0, Mxy = 0, Mxz = 0, Myz = 0, Mzz = 0;
    for (std::size_t k : idx) {
        const double dx = xy[k][0] - xb;
        const double dy = xy[k][1] - zb;
        const double z2 = dx * dx + dy * dy;
        Mxx += dx * dx; Myy += dy * dy; Mxy += dx * dy;
        Mxz += dx * z2; Myz += dy * z2; Mzz += z2 * z2;
    }
    Mxx /= nd; Myy /= nd; Mxy /= nd; Mxz /= nd; Myz /= nd; Mzz /= nd;
    const double Mz = Mxx + Myy;
    const double Cov = Mxx * Myy - Mxy * Mxy;
    const double Mxz2 = Mxz * Mxz, Myz2 = Myz * Myz;
    const double A3 = 4.0 * Mz;
    const double A2 = -3.0 * Mz * Mz - Mzz;
    const double A1 = Mzz * Mz + 4.0 * Cov * Mz - Mxz2 - Myz2 - Mz * Mz * Mz;
    const double A0 = Mxz2 * Myy + Myz2 * Mxx - Mzz * Cov - 2.0 * Mxz * Myz * Mxy + Mz * Mz * Cov;
    const double A22 = A2 + A2, A33 = A3 + A3 + A3;
    double t = 0, y = 1e20;
    for (int it = 0; it < 100; ++it) {
        const double y_old = y;
        y = A0 + t * (A1 + t * (A2 + t * A3));
        if (std::abs(y) > std::abs(y_old)) break;
        const double Dy = A1 + t * (A22 + t * A33);
        if (std::abs(Dy) < 1e-30) break;
        const double t_old = t;
        t = t_old - y / Dy;
        if (std::abs(t - t_old) < 1e-12 * std::abs(t == 0 ? 1.0 : t)) break;
        if (t < 0) { t = 0; break; }
    }
    const double DET = t * t - t * Mz + Cov;
    if (std::abs(DET) < 1e-30) return {};
    const double Xc = (Mxz * (Myy - t) - Myz * Mxy) / (DET * 2.0);
    const double Yc = (Myz * (Mxx - t) - Mxz * Mxy) / (DET * 2.0);
    double cx = Xc + xb, cy = Yc + zb;
    const double R2 = Xc * Xc + Yc * Yc + Mz + 2.0 * t;
    if (!(R2 > 0)) return {};
    double R = std::sqrt(R2);
    // Gauss-Newton on Σ(d − R)² — Pratt is already unbiased to leading order;
    // this nails the geometric optimum in 2-5 iterations.
    for (int it = 0; it < 20; ++it) {
        double Jxx=0, Jxy=0, Jxr=0, Jyy=0, Jyr=0, Jrr=0;
        double bx=0, by=0, br=0;
        for (std::size_t k : idx) {
            const double dx = xy[k][0] - cx;
            const double dy = xy[k][1] - cy;
            const double d = std::sqrt(dx * dx + dy * dy);
            if (d < 1e-9) continue;
            const double r = d - R;
            const double Jx = -dx / d, Jy = -dy / d, JR = -1.0;
            Jxx += Jx * Jx; Jxy += Jx * Jy; Jxr += Jx * JR;
            Jyy += Jy * Jy; Jyr += Jy * JR; Jrr += JR * JR;
            bx += Jx * r; by += Jy * r; br += JR * r;
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
        cx += dxc; cy += dyc; R += dRc;
        if (std::abs(dxc) + std::abs(dyc) + std::abs(dRc) < 1e-6) break;
    }
    if (!(R > 0) || !std::isfinite(R)) return {};
    return { cx, cy, R, true };
}

/// Brute-force 3D kNN: for each hit return up to k nearest-neighbour indices.
std::vector<std::vector<std::size_t>>
precomputeKNN(const std::vector<std::array<double, 2>>& xy, const std::vector<double>& y, int k)
{
    const std::size_t N = xy.size();
    std::vector<std::vector<std::size_t>> nb(N);
    for (std::size_t i = 0; i < N; ++i)
    {
        std::vector<std::pair<double, std::size_t>> dists;
        dists.reserve(N - 1);
        for (std::size_t j = 0; j < N; ++j)
        {
            if (j == i)
                continue;
            const double dx = xy[j][0] - xy[i][0];
            const double dy = xy[j][1] - xy[i][1];
            const double dz = y[j] - y[i];
            dists.emplace_back(dx * dx + dy * dy + dz * dz, j);
        }
        const std::size_t kk = std::min(static_cast<std::size_t>(k), dists.size());
        if (kk == 0)
            continue;
        std::nth_element(dists.begin(), dists.begin() + kk, dists.end());
        nb[i].reserve(kk);
        for (std::size_t m = 0; m < kk; ++m)
            nb[i].push_back(dists[m].second);
    }
    return nb;
}

} // namespace

std::vector<R3BGTPCTrackData>
R3BGTPCTrackFinderRiemann::FindTracks(const std::vector<R3BGTPCHitData>& hits)
{
    std::vector<R3BGTPCTrackData> outTracks;

    const int nHits = static_cast<int>(hits.size());
    if (nHits < fMinHitsPerTrack)
        return outTracks;

    // R3B GLAD has B = (0, B_y, 0), so the bending plane is (x, z) and the
    // helix axis is y. Feed (x, z) to the algorithm as the "xy plane" and y
    // as the "z" coordinate (i.e. the helix-axis variable).
    std::vector<std::array<double, 2>> xy;
    std::vector<double> yVals;
    xy.reserve(nHits);
    yVals.reserve(nHits);
    for (const auto& h : hits)
    {
        xy.push_back({ h.GetX(), h.GetZ() });
        yVals.push_back(h.GetY());
    }

    const auto neighbors = precomputeKNN(xy, yVals, fK_NN);

    std::vector<bool> assigned(nHits, false);
    std::mt19937 rng(fSeed);

    int tracksFound = 0;
    for (int trk = 0; trk < fMaxTracks; ++trk)
    {
        std::vector<std::size_t> live;
        for (std::size_t i = 0; i < xy.size(); ++i)
            if (!assigned[i])
                live.push_back(i);
        if (static_cast<int>(live.size()) < fMinHitsPerTrack)
            break;

        Circle bestCircle;
        ZLine bestLine;
        std::vector<std::size_t> bestInliers;

        std::uniform_int_distribution<std::size_t> pickLive(0, live.size() - 1);

        for (int it = 0; it < fMaxIterations; ++it)
        {
            const std::size_t s0 = live[pickLive(rng)];

            std::vector<std::size_t> liveNeighbors;
            liveNeighbors.reserve(neighbors[s0].size());
            for (std::size_t nb : neighbors[s0])
                if (!assigned[nb])
                    liveNeighbors.push_back(nb);
            if (liveNeighbors.size() < 2)
                continue;

            std::uniform_int_distribution<std::size_t> pickNb(0, liveNeighbors.size() - 1);
            const std::size_t s1 = liveNeighbors[pickNb(rng)];
            const std::size_t s2 = liveNeighbors[pickNb(rng)];
            if (s1 == s2 || s0 == s1 || s0 == s2)
                continue;

            // Optional: require the seed triple to span at least fMinSeedSpread
            // cm pairwise. Anchors R against hit noise on near-straight tracks.
            if (fMinSeedSpread > 0.0)
            {
                auto dist2 = [&](std::size_t a, std::size_t b)
                {
                    const double dx = xy[a][0] - xy[b][0];
                    const double dy = xy[a][1] - xy[b][1];
                    const double dz = yVals[a] - yVals[b];
                    return dx * dx + dy * dy + dz * dz;
                };
                const double s2thr = fMinSeedSpread * fMinSeedSpread;
                if (dist2(s0, s1) < s2thr || dist2(s0, s2) < s2thr || dist2(s1, s2) < s2thr)
                    continue;
            }

            const Circle cand = circleThroughThree(xy[s0], xy[s1], xy[s2]);
            if (!cand.valid)
                continue;
            if (fMinCircleR > 0.0 && cand.r < fMinCircleR)
                continue;

            const double phi0 = std::atan2(xy[s0][1] - cand.cy, xy[s0][0] - cand.cx);
            const double phi1 = std::atan2(xy[s1][1] - cand.cy, xy[s1][0] - cand.cx);
            const double phi2 = std::atan2(xy[s2][1] - cand.cy, xy[s2][0] - cand.cx);
            const ZLine zline = zLineFromSeeds(phi0, yVals[s0], phi1, yVals[s1], phi2, yVals[s2]);

            std::vector<std::size_t> inl;
            inl.reserve(live.size());
            for (std::size_t i : live)
            {
                if (pointCircleDist(cand, xy[i][0], xy[i][1]) >= fInlierDist)
                    continue;
                if (zline.valid)
                {
                    const double phi = std::atan2(xy[i][1] - cand.cy, xy[i][0] - cand.cx);
                    if (std::abs(yVals[i] - (zline.a + zline.b * phi)) >= fZInlierDist)
                        continue;
                }
                inl.push_back(i);
            }

            if (inl.size() > bestInliers.size())
            {
                bestInliers = std::move(inl);
                bestCircle = cand;
                bestLine = zline;
            }
        }

        if (static_cast<int>(bestInliers.size()) < fMinHitsPerTrack)
            break;

        const std::vector<std::size_t> ransacCluster = bestInliers;
        auto rejectAndContinue = [&]()
        {
            for (std::size_t i : ransacCluster)
                assigned[i] = true;
        };

        const Circle refinedCircle = kasaRefit(xy, bestInliers);
        if (!refinedCircle.valid)
        {
            rejectAndContinue();
            continue;
        }
        const ZLine refinedLine = zLineLSQ(xy, yVals, bestInliers, refinedCircle);

        std::vector<std::size_t> finalInliers;
        finalInliers.reserve(bestInliers.size());
        for (std::size_t i : live)
        {
            if (pointCircleDist(refinedCircle, xy[i][0], xy[i][1]) >= fInlierDist)
                continue;
            if (refinedLine.valid)
            {
                const double phi = std::atan2(xy[i][1] - refinedCircle.cy, xy[i][0] - refinedCircle.cx);
                if (std::abs(yVals[i] - (refinedLine.a + refinedLine.b * phi)) >= fZInlierDist)
                    continue;
            }
            finalInliers.push_back(i);
        }
        if (static_cast<int>(finalInliers.size()) < fMinHitsPerTrack)
        {
            rejectAndContinue();
            continue;
        }
        bestCircle = refinedCircle;
        bestInliers.swap(finalInliers);

        bestInliers = largestContiguousArc(xy, bestInliers, bestCircle, fMaxPhiGap);
        if (static_cast<int>(bestInliers.size()) < fMinHitsPerTrack)
        {
            rejectAndContinue();
            continue;
        }

        if (fUseArcWalkExtend)
        {
            std::vector<std::size_t> liveAll;
            liveAll.reserve(live.size());
            for (std::size_t i : live)
                liveAll.push_back(i);
            bestInliers = arcWalkExtend(xy, yVals, bestInliers, liveAll, bestCircle,
                                        fInlierDist, fZInlierDist,
                                        fArcWalkWindow, fArcWalkMaxMiss);
        }

        // Final Pratt + Gauss-Newton refit on the chosen inlier set. Algebraic
        // Kasa biases R low when sagitta is comparable to hit noise (~1 mm in
        // R3B Prototype), so use it only as inlier classifier and let Pratt
        // produce the reported geometric circle.
        if (const Circle gnFinal = prattGNRefit(xy, bestInliers); gnFinal.valid)
            bestCircle = gnFinal;

        // Sort inliers along φ around the refined circle so the hit list is
        // ordered along the track — downstream consumers expect this.
        std::sort(bestInliers.begin(), bestInliers.end(),
                  [&](std::size_t a, std::size_t b)
                  {
                      const double pa = std::atan2(xy[a][1] - bestCircle.cy, xy[a][0] - bestCircle.cx);
                      const double pb = std::atan2(xy[b][1] - bestCircle.cy, xy[b][0] - bestCircle.cx);
                      return pa < pb;
                  });

        R3BGTPCTrackData track;
        track.SetTrackId(tracksFound++);
        for (std::size_t i : bestInliers)
        {
            R3BGTPCHitData h = hits[i];
            track.AddHit(h);
            assigned[i] = true;
        }

        // Algorithm-internal circle + z-line are essentially free output; expose
        // them on the track so a downstream stage can read the PR geometry
        // without re-fitting. cot(θ_y) = sign(ψ̇) · b / R (b = dy/dphi).
        track.SetGeoCenter({ bestCircle.cx, bestCircle.cy });
        track.SetGeoRadius(bestCircle.r);
        if (refinedLine.valid && bestInliers.size() >= 2)
        {
            const double phi_first = std::atan2(xy[bestInliers.front()][1] - bestCircle.cy,
                                                xy[bestInliers.front()][0] - bestCircle.cx);
            const double phi_last = std::atan2(xy[bestInliers.back()][1] - bestCircle.cy,
                                               xy[bestInliers.back()][0] - bestCircle.cx);
            double dphi = phi_last - phi_first;
            while (dphi > TMath::Pi())
                dphi -= 2.0 * TMath::Pi();
            while (dphi <= -TMath::Pi())
                dphi += 2.0 * TMath::Pi();
            const double sign_psi_dot = (dphi >= 0) ? 1.0 : -1.0;
            const double cot_th = refinedLine.b * sign_psi_dot / bestCircle.r;
            track.SetGeoTheta(std::atan2(1.0, cot_th));
        }
        else
        {
            track.SetGeoTheta(TMath::Pi() / 2.0);
        }

        outTracks.push_back(std::move(track));
    }

    std::cout << cRED << " R3BGTPCTrackFinderRiemann: tracks found " << outTracks.size() << cNORMAL << "\n";
    return outTracks;
}
