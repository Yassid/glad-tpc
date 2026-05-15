/// @file smoke_geo.C
/// @brief Smoke test for R3BGTPCTrackFinder::SetTrackInitialParameters.
///
/// Builds a synthetic R3BGTPCTrackData with hits sampled along a known
/// helix, calls SetTrackInitialParameters, and checks that the recovered
/// GeoCenter / GeoRadius / GeoTheta match the truth within tolerance.
///
/// Run: root -b -q smoke_geo.C

#include <Math/Vector3D.h>
#include <iostream>

void smoke_geo()
{
    gSystem->Load("libR3BGTPCData");
    gSystem->Load("libR3BOpenKF");
    gSystem->Load("libR3BGTPCReconstruction");

    // ---- truth ----
    const double R_true = 333.33;             // mm
    const double cx_true = 0.0;
    const double cy_true = -R_true;
    const double theta_true = 75.0 * M_PI / 180.0; // upward-going (75 deg, dz/ds > 0)
    const double cot_th = std::cos(theta_true) / std::sin(theta_true);
    const int N = 30;

    R3BGTPCTrackData track;
    for (int i = 0; i < N; ++i) {
        // CW motion: psi DECREASES with arc s
        const double s = (i * 250.0) / (N - 1); // 250 mm total arc
        const double psi = M_PI / 2.0 - s / R_true;
        const double x = cx_true + R_true * std::cos(psi);
        const double y = cy_true + R_true * std::sin(psi);
        const double z = s * cot_th;

        R3BGTPCHitData h(x, y, z, 0., 1.);
        track.AddHit(h);
    }
    track.SetTrackId(0);

    R3BGTPCTrackFinder finder;
    finder.SetTrackInitialParameters(track);

    const auto c = track.GetGeoCenter();
    const double R_fit = track.GetGeoRadius();
    const double theta_fit = track.GetGeoTheta();

    std::cout << "\n=== Kasa fit ===\n";
    std::cout << "Truth   center (" << cx_true << ", " << cy_true << ")  R = " << R_true
              << "  theta = " << theta_true * 180. / M_PI << " deg\n";
    std::cout << "Fit     center (" << c.first << ", " << c.second << ")  R = " << R_fit
              << "  theta = " << theta_fit * 180. / M_PI << " deg\n";

    const double dC = std::sqrt((c.first - cx_true) * (c.first - cx_true)
                                + (c.second - cy_true) * (c.second - cy_true));
    const double dR = R_fit - R_true;
    const double dTh = (theta_fit - theta_true) * 180. / M_PI;

    bool ok = (dC < 1.0) && (std::abs(dR) < 1.0) && (std::abs(dTh) < 1.0);
    std::cout << "Δcenter = " << dC << " mm,  ΔR = " << dR << " mm,  Δθ = " << dTh
              << " deg  " << (ok ? "[PASS]" : "[FAIL]") << "\n";
}
