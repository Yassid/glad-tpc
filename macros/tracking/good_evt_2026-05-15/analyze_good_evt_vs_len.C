/// @file analyze_good_evt_vs_len.C
/// @brief σ_p/p and bias binned by reconstructed track chord length, on
///        good_evt (realistic He-3 + π⁻ kinematics).
///
/// Chord length = √(Δz² + Δx²) of the cluster centroid bounding box (cm).

void analyze_good_evt_vs_len()
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gSystem->Load("libR3BOpenKF");
    gStyle->SetOptStat(0);

    TString workDir = gSystem->Getenv("VMCWORKDIR");
    TFile fSim(workDir + "/glad-tpc/macros/sim/Prototype/sim.root");
    TFile fTrk(workDir + "/glad-tpc/macros/tracking/output_tracking.root");
    TFile fUKF(workDir + "/glad-tpc/macros/tracking/output_ukf.root");
    auto* tSim = (TTree*)fSim.Get("evt");
    auto* tTrk = (TTree*)fTrk.Get("evt");
    auto* tUKF = (TTree*)fUKF.Get("evt");
    auto* mc = new TClonesArray("R3BMCTrack");
    auto* trks = new TClonesArray("R3BGTPCTrackData");
    auto* fits = new TClonesArray("R3BGTPCFittedTrackData");
    tSim->SetBranchAddress("MCTrack", &mc);
    tTrk->SetBranchAddress("GTPCTrackData", &trks);
    tUKF->SetBranchAddress("GTPCFittedTrackData", &fits);

    // 4 bins by chord length (cm): [3,8), [8,12), [12,16), [16,∞)
    const std::vector<std::pair<double, double>> bins = {
        { 3, 8 }, { 8, 12 }, { 12, 16 }, { 16, 30 }
    };
    std::vector<TH1F*> hs;
    std::vector<TH2F*> h2s;
    std::vector<int> nFit(bins.size(), 0);
    for (size_t b = 0; b < bins.size(); ++b) {
        hs.push_back(new TH1F(Form("h_%zu", b),
                              Form("chord %.0f-%.0f cm;(p_{fit}-p_{MC})/p_{MC};events", bins[b].first, bins[b].second),
                              80, -1.2, 0.8));
        h2s.push_back(new TH2F(Form("h2_%zu", b),
                               Form("chord %.0f-%.0f cm;p_{MC} (MeV/c);(p_{fit}-p_{MC})/p_{MC}", bins[b].first, bins[b].second),
                               30, 100, 900, 60, -1.2, 0.8));
    }

    int total = 0;
    for (Long64_t i = 0; i < std::min({ tSim->GetEntries(), tTrk->GetEntries(), tUKF->GetEntries() }); ++i) {
        tSim->GetEntry(i);
        tTrk->GetEntry(i);
        tUKF->GetEntry(i);
        if (fits->GetEntries() == 0 || trks->GetEntries() == 0) continue;
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

        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        auto& cls = *tr->GetHitClusterArray();
        if (cls.size() < 3) continue;
        double zmn = 1e9, zmx = -1e9, xmn = 1e9, xmx = -1e9;
        for (auto& c : cls) {
            zmn = std::min(zmn, (double)c.GetZ()); zmx = std::max(zmx, (double)c.GetZ());
            xmn = std::min(xmn, (double)c.GetX()); xmx = std::max(xmx, (double)c.GetX());
        }
        const double chord_cm = std::sqrt((zmx - zmn) * (zmx - zmn) + (xmx - xmn) * (xmx - xmn));
        const double res = (pfit - pmc) / pmc;

        ++total;
        for (size_t b = 0; b < bins.size(); ++b) {
            if (chord_cm >= bins[b].first && chord_cm < bins[b].second) {
                hs[b]->Fill(res);
                h2s[b]->Fill(pmc, res);
                ++nFit[b];
                break;
            }
        }
    }

    auto* c = new TCanvas("c", "σ vs chord", 1600, 900);
    c->Divide(4, 2);
    std::cout << "\nchord (cm)   N   core σ   bias    hist RMS\n";
    std::cout << "----------------------------------------------\n";
    for (size_t b = 0; b < bins.size(); ++b) {
        c->cd(b + 1);
        hs[b]->Draw();
        if (nFit[b] >= 5) {
            hs[b]->Fit("gaus", "Q", "", -0.5, 0.5);
            auto* g = hs[b]->GetFunction("gaus");
            std::printf("%4.0f-%-4.0f  %3d   %5.3f    %+6.3f    %5.3f\n",
                        bins[b].first, bins[b].second, nFit[b],
                        g ? g->GetParameter(2) : 0.0,
                        g ? g->GetParameter(1) : 0.0,
                        hs[b]->GetRMS());
        } else {
            std::printf("%4.0f-%-4.0f  %3d   (too few)\n", bins[b].first, bins[b].second, nFit[b]);
        }
        c->cd(b + 5);
        h2s[b]->Draw("COLZ");
    }
    std::cout << "----------------------------------------------\n";
    std::cout << "Total fits: " << total << "\n";
    c->SaveAs("good_evt_sigma_vs_chord.png");
}
