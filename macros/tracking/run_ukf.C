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
    TString inFile = workDir + "/glad-tpc/macros/tracking/" + fileName;
    TString outFile = workDir + "/glad-tpc/macros/tracking/" + outName;

    auto* fRun = new FairRunAna();
    fRun->SetSource(new FairFileSource(inFile));
    fRun->SetOutputFile(outFile.Data());

    auto* fitTask = new R3BGTPCTrack2Fit();
    // pion default — match the HYDRA scan in ATTPCROOT
    fitTask->SetParticle(1.602176634e-19, 139.57039);
    fitTask->SetProjectile(1, 1, 0.1395); // π± in u
    fitTask->SetBField({ 0., 0., 2.0 });  // T
    fitTask->SetMeasurementSigma(1.0);
    fitTask->SetMomentumSigmaFrac(0.1);
    fitTask->SetMinClusters(5);
    fitTask->SetEnableEnergyStraggling(true);

    fRun->AddTask(fitTask);
    fRun->Init();
    fRun->Run(0, 0);

    std::cout << "UKF stage finished. Output: " << outFile << std::endl;
}
