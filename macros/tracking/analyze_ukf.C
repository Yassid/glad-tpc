/// @file analyze_ukf.C
/// @brief Compare UKF kinematics vs MC truth on the HYDRA pipeline output.
///
/// Reads:
///   ../sim/Prototype/sim.root        — MC tracks (R3BMCTrack)
///   output_ukf.root                  — UKF fits (R3BGTPCFittedTrackData)
///
/// Prints (KE_MC, KE_fit, ΔKE/KE) for each fitted event.

void analyze_ukf()
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gSystem->Load("libR3BOpenKF");

    TString workDir = gSystem->Getenv("VMCWORKDIR");
    TString mcPath = workDir + "/glad-tpc/macros/sim/Prototype/sim.root";
    TString ukfPath = workDir + "/glad-tpc/macros/tracking/output_ukf.root";
    std::cout << "MC file:  " << mcPath << "\n";
    std::cout << "UKF file: " << ukfPath << "\n";
    TFile fMC(mcPath.Data());
    TFile fUKF(ukfPath.Data());
    if (fMC.IsZombie() || fUKF.IsZombie()) {
        std::cerr << "Failed to open input files.\n";
        return;
    }

    auto* tMC = (TTree*)fMC.Get("evt");
    auto* tUKF = (TTree*)fUKF.Get("evt");
    if (!tMC || !tUKF) {
        std::cerr << "Missing trees.\n";
        return;
    }

    // List branches so we know what's available
    std::cout << "MC branches:\n";
    auto* mcBranches = tMC->GetListOfBranches();
    for (int i = 0; i < mcBranches->GetEntries(); ++i)
        std::cout << "  " << mcBranches->At(i)->GetName() << "\n";

    std::cout << "\nUKF-stage branches:\n";
    auto* ukfBranches = tUKF->GetListOfBranches();
    for (int i = 0; i < ukfBranches->GetEntries(); ++i)
        std::cout << "  " << ukfBranches->At(i)->GetName() << "\n";

    auto* mcTracks = new TClonesArray("R3BMCTrack");
    auto* fitTracks = new TClonesArray("R3BGTPCFittedTrackData");
    tMC->SetBranchAddress("MCTrack", &mcTracks);
    tUKF->SetBranchAddress("GTPCFittedTrackData", &fitTracks);

    const double mass_pi = 0.13957;
    const Long64_t N = std::min(tMC->GetEntries(), tUKF->GetEntries());
    std::cout << "\nEntries: " << N << "\n\n";
    std::cout << "evt  MC_KE(MeV)  Fit_KE(MeV)  dKE/KE  chi2/ndf\n";
    for (Long64_t i = 0; i < N; ++i) {
        tMC->GetEntry(i);
        tUKF->GetEntry(i);
        if (fitTracks->GetEntries() == 0)
            continue;
        // Pick the first non-beam MC particle (pion).
        double keMC = 0.;
        for (int j = 0; j < mcTracks->GetEntries(); ++j) {
            auto* mc = (R3BMCTrack*)mcTracks->At(j);
            if (mc->GetMotherId() == -1)
                continue;
            const double p = std::sqrt(mc->GetPx() * mc->GetPx() + mc->GetPy() * mc->GetPy() + mc->GetPz() * mc->GetPz());
            keMC = std::sqrt(p * p + mass_pi * mass_pi) - mass_pi;
            break;
        }
        auto* fitted = (R3BGTPCFittedTrackData*)fitTracks->At(0);
        const double keFit_MeV = fitted->GetKinematics().kineticEnergy;
        const double keMC_MeV = keMC * 1000.;
        const double dke = (keFit_MeV - keMC_MeV) / std::max(keMC_MeV, 1e-9);
        std::printf("%3lld  %8.2f  %10.2f  %8.3f  %8.3f\n", i, keMC_MeV, keFit_MeV, dke, fitted->GetChi2OverNdf());
    }
}
