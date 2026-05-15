/// @file smoke_ukf.C
/// @brief Standalone smoke test for R3BGTPCFitterUKF.
///
/// Builds a synthetic R3BGTPCTrackData containing N hit clusters sampled
/// along a circular arc in (x,y) at fixed z (= helix with pz=0). Pushes it
/// through R3BGTPCFitterUKF::FitTrack and prints the recovered kinematics.
///
/// Purpose: verify that the vendored OpenKF + AtPropagator + catima chain
/// links and runs in R3BRoot's build, independently of the full sim →
/// lang → reco → tracking pipeline.
///
/// Run: root -b -q smoke_ukf.C

#include <Math/Vector3D.h>
#include <TMatrixD.h>
#include <iostream>

void smoke_ukf()
{
    gSystem->Load("libR3BGTPCData");
    gSystem->Load("libR3BOpenKF");
    gSystem->Load("libR3BGTPCReconstruction");

    // ---- parameters ----
    const double mass_proton_MeV = 938.272;
    const double charge_C = 1.602176634e-19;
    const double B_T = 2.0;
    const double p_true_MeV = 200.0;
    const double theta_true = M_PI / 2.0; // perpendicular to B
    const double phi0_true = 0.0;
    const int N_clusters = 20;
    const double sigma_xy_mm = 1.0;

    // Helix radius: R[mm] = p[MeV] / (0.3 · B[T]). For a positive charge in
    // +B_z, F = qv×B points toward -y at v=+x, so the centre is at (0, -R)
    // (CW motion in (x,y)). Earlier mistake: centring at (0, +R) gave a CCW
    // trajectory — the UKF then inferred ~3× the seed momentum trying to
    // reconcile the wrong-sign curvature with a positive-charge propagator.
    const double R_mm = p_true_MeV / (0.3 * B_T);
    const double cx = 0.0;
    const double cy = -R_mm;

    // ---- build synthetic R3BGTPCTrackData ----
    auto track = std::make_unique<R3BGTPCTrackData>();
    track->SetTrackId(0);
    // Geo fields will be populated below by SetTrackInitialParameters once
    // the hits have been added — same path the real pipeline takes.

    // Sample clusters at uniform arc length, walking counter-clockwise from
    // (0, 0, 0) along the lower half of the circle.
    const double arc_total = 200.0; // mm
    const double arc_step = arc_total / (N_clusters - 1);
    std::cout << "Synthetic helix: p = " << p_true_MeV << " MeV/c, R = " << R_mm
              << " mm, sampling " << N_clusters << " clusters over " << arc_total
              << " mm arc.\n";

    for (int i = 0; i < N_clusters; ++i) {
        const double s = i * arc_step;
        // CW motion: psi DECREASES with arc length. Start at psi=+pi/2
        // (position (cx, cy + R) = (0, 0)) and sweep clockwise.
        const double psi = M_PI / 2.0 - s / R_mm;
        const double x = cx + R_mm * std::cos(psi);
        const double y = cy + R_mm * std::sin(psi);

        R3BGTPCHitClusterData hc;
        hc.SetX(x);
        hc.SetY(y);
        hc.SetZ(0.0);
        hc.SetEnergy(1.0);
        TMatrixD cov(3, 3);
        cov.Zero();
        cov(0, 0) = cov(1, 1) = cov(2, 2) = sigma_xy_mm * sigma_xy_mm;
        hc.SetCovMatrix(cov);
        hc.SetClusterID(i);

        auto hcp = std::make_shared<R3BGTPCHitClusterData>(hc);
        track->AddClusterHit(hcp);

        // Also push as raw hits, since SetTrackInitialParameters works on
        // the hit array (representative of the production pipeline).
        R3BGTPCHitData rawHit(x, y, 0.0, 0.0, 1.0);
        track->AddHit(rawHit);
    }

    std::cout << "Track has " << track->GetHitClusterArray()->size() << " clusters.\n";

    // Populate fGeoCenter / fGeoRadius / fGeoTheta from the hits, the way
    // R3BGTPCTrackFinder::clustersToTrack does it in the real pipeline.
    R3BGTPCTrackFinder finder;
    finder.SetTrackInitialParameters(*track);
    std::cout << "Kasa: center (" << track->GetGeoCenter().first << ", " << track->GetGeoCenter().second
              << ")  R = " << track->GetGeoRadius() << " mm  theta = "
              << track->GetGeoTheta() * 180. / M_PI << " deg\n";

    // ---- configure energy loss (hydrogen gas at HYDRA-like density) ----
    auto eloss = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5);
    eloss->SetProjectile(1, 1, 1);
    eloss->SetMaterial({ std::make_tuple(1, 1, 1) });

    // ---- run the fitter ----
    R3BGTPCFitterUKF fitter(charge_C, mass_proton_MeV, std::move(eloss));
    fitter.SetBField({ 0., 0., B_T });
    // No SetMomentumSeed: Brho is derived from the geo fields populated
    // by SetTrackInitialParameters above. This is the end-to-end flow that
    // the production FairTask wrapper will use.
    fitter.SetMeasurementSigma(sigma_xy_mm);
    fitter.SetMinClusters(5);
    // Smoke test generates positions directly in mm — override the R3B
    // default (cm → mm conversion factor = 10).
    fitter.SetInputUnit_mm(1.0);

    auto fitted = fitter.FitTrack(track.get());

    if (!fitted) {
        std::cerr << "[FAIL] FitTrack returned nullptr.\n";
        return;
    }
    if (!fitted->IsConverged()) {
        std::cerr << "[FAIL] Fit did not converge.\n";
        return;
    }

    const auto& kin = fitted->GetKinematics();
    const double KE_true = std::sqrt(p_true_MeV * p_true_MeV + mass_proton_MeV * mass_proton_MeV) - mass_proton_MeV;

    std::cout << "\n=== UKF result ===\n";
    std::cout << "Truth     KE = " << KE_true << " MeV  theta = " << theta_true * 180. / M_PI
              << " deg  phi = " << phi0_true * 180. / M_PI << " deg\n";
    std::cout << "Fit       KE = " << kin.kineticEnergy << " MeV  theta = " << kin.theta * 180. / M_PI
              << " deg  phi = " << kin.phi * 180. / M_PI << " deg\n";
    std::cout << "Vertex (mm) = (" << fitted->GetVertex().X() << ", " << fitted->GetVertex().Y() << ", "
              << fitted->GetVertex().Z() << ")\n";
    std::cout << "chi2 / ndf = " << fitted->GetChi2() << " / " << fitted->GetNdf()
              << " = " << fitted->GetChi2OverNdf() << "\n";
    std::cout << "Smoothed positions: " << fitted->GetSmoothedPositions().size() << "\n";

    const double dKE = (kin.kineticEnergy - KE_true) / KE_true;
    if (std::abs(dKE) < 0.1)
        std::cout << "[PASS] |ΔKE/KE| = " << dKE << " (< 0.1)\n";
    else
        std::cerr << "[WARN] |ΔKE/KE| = " << dKE << " out of expected tolerance\n";
}
