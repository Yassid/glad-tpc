/// @file analyze_scan_point.C
/// @brief Compute σ_p/p (Gaussian core) for a single scan point.
///
/// Prints a single RESULT line to stdout that scan_p.sh appends to the
/// CSV summary:
///   RESULT <p_MeV> <n_fit> <n_thr> <sigma_core> <bias_core>

void analyze_scan_point(int p_MeV)
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gSystem->Load("libR3BOpenKF");

    TString workDir = gSystem->Getenv("VMCWORKDIR");
    TString mcPath = TString::Format("%s/glad-tpc/macros/sim/Prototype/sim_p%d.root", workDir.Data(), p_MeV);
    TString ukfPath = TString::Format("%s/glad-tpc/macros/tracking/output_ukf_p%d.root", workDir.Data(), p_MeV);

    TFile fMC(mcPath);
    TFile fUKF(ukfPath);
    if (fMC.IsZombie() || fUKF.IsZombie()) {
        std::cerr << "[analyze_scan_point] missing input files at p=" << p_MeV << "\n";
        return;
    }
    auto* tMC = (TTree*)fMC.Get("evt");
    auto* tUKF = (TTree*)fUKF.Get("evt");
    if (!tMC || !tUKF)
        return;

    auto* mcTracks = new TClonesArray("R3BMCTrack");
    auto* fitTracks = new TClonesArray("R3BGTPCFittedTrackData");
    tMC->SetBranchAddress("MCTrack", &mcTracks);
    tUKF->SetBranchAddress("GTPCFittedTrackData", &fitTracks);

    const double mass_pi = 0.13957;
    const Long64_t N = std::min(tMC->GetEntries(), tUKF->GetEntries());
    // Wide window: at high p in a small chamber the sagitta can be below hit
    // resolution and individual residuals reach -1 (p_fit → 0). Use [-1, 1] and
    // fit a Gaussian to the core within ±0.3 of the histogram median.
    auto* h = new TH1F("h", "", 200, -1.0, 1.0);
    int nFit = 0;
    for (Long64_t i = 0; i < N; ++i) {
        tMC->GetEntry(i);
        tUKF->GetEntry(i);
        if (fitTracks->GetEntries() == 0)
            continue;
        // The box-generated pion is the primary (MotherId == -1). Pick it.
        double pmc_GeV = 0;
        for (int j = 0; j < mcTracks->GetEntries(); ++j) {
            auto* mc = (R3BMCTrack*)mcTracks->At(j);
            if (mc->GetMotherId() != -1)
                continue;
            if (std::abs(mc->GetPdgCode()) != 211)
                continue; // require pion
            pmc_GeV = std::sqrt(mc->GetPx() * mc->GetPx() + mc->GetPy() * mc->GetPy() + mc->GetPz() * mc->GetPz());
            break;
        }
        if (pmc_GeV <= 0)
            continue;
        const double pmc_MeV = pmc_GeV * 1000;
        auto* fitted = (R3BGTPCFittedTrackData*)fitTracks->At(0);
        const double ke_fit = fitted->GetKinematics().kineticEnergy; // MeV
        if (!std::isfinite(ke_fit) || ke_fit <= 0)
            continue;
        const double p_fit = std::sqrt(ke_fit * ke_fit + 2 * ke_fit * mass_pi * 1000);
        h->Fill((p_fit - pmc_MeV) / pmc_MeV);
        ++nFit;
    }

    double mu = 0., sigma = 0.;
    if (nFit > 5) {
        // Core fit: ±0.3 around peak
        const double xpk = h->GetBinCenter(h->GetMaximumBin());
        h->Fit("gaus", "Q0", "", xpk - 0.3, xpk + 0.3);
        auto* g = h->GetFunction("gaus");
        if (g) {
            mu = g->GetParameter(1);
            sigma = g->GetParameter(2);
        }
    }
    const double mean = h->GetMean();
    const double rms = h->GetRMS();
    std::printf("RESULT %d %d %lld %.4f %.4f  mean=%.3f rms=%.3f\n",
                p_MeV, nFit, N, sigma, mu, mean, rms);
}
