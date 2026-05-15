/// @file analyze_good_evt.C
/// @brief σ_p/p vs p_MC using the realistic good_evt generator (He-3 + π⁻).

void analyze_good_evt()
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gSystem->Load("libR3BOpenKF");
    gStyle->SetOptStat(0);

    TString workDir = gSystem->Getenv("VMCWORKDIR");
    TFile fSim(workDir + "/glad-tpc/macros/sim/Prototype/sim.root");
    TFile fUKF(workDir + "/glad-tpc/macros/tracking/output_ukf.root");
    auto* tSim = (TTree*)fSim.Get("evt");
    auto* tUKF = (TTree*)fUKF.Get("evt");
    auto* mc = new TClonesArray("R3BMCTrack");
    auto* fits = new TClonesArray("R3BGTPCFittedTrackData");
    tSim->SetBranchAddress("MCTrack", &mc);
    tUKF->SetBranchAddress("GTPCFittedTrackData", &fits);

    auto* gRes = new TGraph();
    gRes->SetMarkerStyle(20);
    gRes->SetMarkerColor(kBlue + 1);
    auto* hRes = new TH1F("hRes", "(p_fit-p_MC)/p_MC;residual;events", 100, -1.5, 1.5);
    int n = 0, nFit = 0;
    for (Long64_t i = 0; i < std::min(tSim->GetEntries(), tUKF->GetEntries()); ++i) {
        tSim->GetEntry(i);
        tUKF->GetEntry(i);
        if (fits->GetEntries() == 0) continue;
        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mc->At(j);
            if (m->GetMotherId() == -1 && std::abs(m->GetPdgCode()) == 211) { pi = m; break; }
        }
        if (!pi) continue;
        const double pmc = std::sqrt(pi->GetPx() * pi->GetPx() + pi->GetPy() * pi->GetPy() + pi->GetPz() * pi->GetPz()) * 1000;
        auto* fit = (R3BGTPCFittedTrackData*)fits->At(0);
        const double ke = fit->GetKinematics().kineticEnergy;
        if (!std::isfinite(ke) || ke <= 0) continue;
        const double pfit = std::sqrt(ke * ke + 2 * ke * 139.57);
        gRes->SetPoint(gRes->GetN(), pmc, (pfit - pmc) / pmc);
        hRes->Fill((pfit - pmc) / pmc);
        ++nFit;
    }
    std::cout << "good_evt: fits=" << nFit << "/" << std::min(tSim->GetEntries(), tUKF->GetEntries()) << "\n";

    auto* c = new TCanvas("c", "good_evt analysis", 1500, 600);
    c->Divide(2, 1);
    c->cd(1);
    gRes->SetTitle(";p_{MC} (MeV/c);(p_{fit}-p_{MC})/p_{MC}");
    gRes->Draw("AP");
    auto* line0 = new TLine(0, 0, 1000, 0);
    line0->SetLineStyle(2);
    line0->Draw();
    c->cd(2);
    hRes->Fit("gaus", "Q", "", -0.5, 0.5);
    hRes->Draw();
    auto* g = hRes->GetFunction("gaus");
    if (g) std::cout << "Core σ_p/p (core ±0.5): " << g->GetParameter(2) << " bias " << g->GetParameter(1) << "  histRMS=" << hRes->GetRMS() << "\n";
    c->SaveAs("good_evt_resolution.png");
}
