/// Dump per-event kinematics for one scan point.

void dump_scan_point(int p_MeV = 800, int nDump = 20)
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gSystem->Load("libR3BOpenKF");

    TString workDir = gSystem->Getenv("VMCWORKDIR");
    TFile fMC(TString::Format("%s/glad-tpc/macros/sim/Prototype/sim_p%d.root", workDir.Data(), p_MeV));
    TFile fUKF(TString::Format("%s/glad-tpc/macros/tracking/output_ukf_p%d.root", workDir.Data(), p_MeV));
    auto* tMC = (TTree*)fMC.Get("evt");
    auto* tUKF = (TTree*)fUKF.Get("evt");

    auto* mcTracks = new TClonesArray("R3BMCTrack");
    auto* fitTracks = new TClonesArray("R3BGTPCFittedTrackData");
    auto* gtpcPoints = new TClonesArray("R3BGTPCPoint");
    tMC->SetBranchAddress("MCTrack", &mcTracks);
    tMC->SetBranchAddress("GTPCPoint", &gtpcPoints);
    tUKF->SetBranchAddress("GTPCFittedTrackData", &fitTracks);

    const double mass_pi = 0.13957;
    const Long64_t N = std::min<Long64_t>(nDump, std::min(tMC->GetEntries(), tUKF->GetEntries()));
    std::cout << "evt nPts MC_px MC_py MC_pz |pMC| KE_MC  fitConverged  KE_fit\n";
    for (Long64_t i = 0; i < N; ++i) {
        tMC->GetEntry(i);
        tUKF->GetEntry(i);
        // primary pion (box generator output: MotherId == -1, |PDG| == 211)
        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mcTracks->GetEntries(); ++j) {
            auto* mc = (R3BMCTrack*)mcTracks->At(j);
            if (mc->GetMotherId() != -1)
                continue;
            if (std::abs(mc->GetPdgCode()) != 211)
                continue;
            pi = mc;
            break;
        }
        if (!pi)
            continue;
        const double p = std::sqrt(pi->GetPx() * pi->GetPx() + pi->GetPy() * pi->GetPy() + pi->GetPz() * pi->GetPz());
        const double keMC = std::sqrt(p * p + mass_pi * mass_pi) - mass_pi;
        const int nPts = gtpcPoints->GetEntries();
        if (fitTracks->GetEntries() == 0) {
            std::printf("%3lld  %3d  %7.3f %7.3f %7.3f  %7.3f  %7.3f  ---\n",
                        i, nPts, pi->GetPx(), pi->GetPy(), pi->GetPz(), p, keMC * 1000);
            continue;
        }
        auto* fitted = (R3BGTPCFittedTrackData*)fitTracks->At(0);
        std::printf("%3lld  %3d  %7.3f %7.3f %7.3f  %7.3f  %7.3f  %d  %.3f MeV\n",
                    i, nPts, pi->GetPx(), pi->GetPy(), pi->GetPz(), p, keMC * 1000,
                    fitted->IsConverged(), fitted->GetKinematics().kineticEnergy);
    }
}
