/// @file run_ukf.C
/// @brief Stage 5 of the HYDRA reco chain: UKF fit on R3BGTPCTrackData.
///
/// Reads ./output_tracking.root, runs R3BGTPCTrack2Fit, writes
/// ./output_ukf.root with the GTPCFittedTrackData branch alongside the
/// upstream branches.
///
/// Run after run_tracking.C:
///   root -b -q run_ukf.C
///
/// For a non-default particle / gas / field, set them before calling Init().

void run_ukf(TString fileName = "output_tracking.root", TString outName = "output_ukf.root")
{
    TString workDir = gSystem->Getenv("VMCWORKDIR");
    const char* suffix_env = gSystem->Getenv("SUFFIX");
    TString suffix = suffix_env ? suffix_env : "";
    // If caller passed the default filenames AND a SUFFIX is set, apply it
    // to keep file naming consistent across the chain.
    if (suffix.Length() > 0 && fileName == "output_tracking.root")
        fileName = "output_tracking" + suffix + ".root";
    if (suffix.Length() > 0 && outName == "output_ukf.root")
        outName = "output_ukf" + suffix + ".root";
    TString inFile = workDir + "/glad-tpc/macros/tracking/" + fileName;
    TString outFile = workDir + "/glad-tpc/macros/tracking/" + outName;

    auto* fRun = new FairRunAna();
    fRun->SetSource(new FairFileSource(inFile));
    fRun->SetOutputFile(outFile.Data());

    auto* fitTask = new R3BGTPCTrack2Fit();
    // pion default — match the HYDRA scan in ATTPCROOT
    fitTask->SetParticle(1.602176634e-19, 139.57039);
    fitTask->SetProjectile(1, 1, 0.1395); // π± in u
    fitTask->SetBField({ 0., 2.0, 0. });  // T — R3B GLAD: horizontal dipole along +y
    fitTask->SetMeasurementSigma(1.0);
    fitTask->SetMomentumSigmaFrac(0.1);
    // Seed override: SEED_P_MEV env var sets a fixed initial momentum
    // (overrides Brho from PRA). Useful for diagnosing whether σ_p/p
    // is seed-limited vs measurement-limited in small chambers.
    const char* seed_env = gSystem->Getenv("SEED_P_MEV");
    if (seed_env && std::atof(seed_env) > 0) {
        const double p_seed = std::atof(seed_env);
        fitTask->SetMomentumSeed(p_seed);
        std::cout << "[run_ukf] Using seed override p = " << p_seed << " MeV/c\n";
    }
    fitTask->SetMinClusters(5);
    fitTask->SetEnableEnergyStraggling(true);

    fRun->AddTask(fitTask);
    fRun->Init();
    fRun->Run(0, 0);

    std::cout << "UKF stage finished. Output: " << outFile << std::endl;
}
