// Per-event seed vs UKF momentum tuple for box-gen 400 MeV/c output.
// We want to know: at 400 MeV/c the seed σ is 4.7 % (vs 2.1-2.8 % plateau)
// and UKF σ exploded to 42 %. Is the seed actually bimodal — a good core
// + catastrophic outliers — or just uniformly broader? If bimodal, can we
// flag the bad-seed events on a per-event quantity before the UKF runs?
void probe_400()
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    TString wd = gSystem->Getenv("VMCWORKDIR");
    TString suf = "_p400";
    TFile fS(wd + "/glad-tpc/macros/sim/Prototype/sim" + suf + ".root");
    TFile fT(wd + "/glad-tpc/macros/tracking/output_tracking" + suf + ".root");
    TFile fU(wd + "/glad-tpc/macros/tracking/output_ukf" + suf + ".root");
    auto* tS = (TTree*)fS.Get("evt");
    auto* tT = (TTree*)fT.Get("evt");
    auto* tU = (TTree*)fU.Get("evt");
    auto* mc     = new TClonesArray("R3BMCTrack");
    auto* trks   = new TClonesArray("R3BGTPCTrackData");
    auto* fitted = new TClonesArray("R3BGTPCFittedTrackData");
    tS->SetBranchAddress("MCTrack", &mc);
    tT->SetBranchAddress("GTPCTrackData", &trks);
    tU->SetBranchAddress("GTPCFittedTrackData", &fitted);
    const double m_pi = 139.57039;
    const double B = 2.0;

    int N = 0;
    auto* hSeed = new TH1F("hSeed",";p_seed/p_MC - 1;", 80, -1, 1);
    auto* hUKF  = new TH1F("hUKF", ";p_UKF /p_MC - 1;", 80, -1, 1);
    auto* hSeedVsUKF = new TH2F("h2","p_seed vs p_UKF;p_seed/p_MC - 1;p_UKF/p_MC - 1", 60,-0.5,0.5, 60,-1.0,1.0);

    // Per-event diagnostic: identify what predicts a bad seed
    auto* hChord_g  = new TH1F("hCg","good seed (<10%) ; chord [cm];", 20, 0, 25);
    auto* hChord_b  = new TH1F("hCb","bad seed (>10%); chord [cm];", 20, 0, 25);
    auto* hNhits_g  = new TH1F("hNg","good; N hits;", 20, 0, 200);
    auto* hNhits_b  = new TH1F("hNb","bad; N hits;", 20, 0, 200);
    auto* hSinTh_g  = new TH1F("hSg","good; sin#theta_y;", 20, 0, 1.05);
    auto* hSinTh_b  = new TH1F("hSb","bad; sin#theta_y;", 20, 0, 1.05);
    auto* hRfit_g   = new TH1F("hRg","good; R_fit [cm];", 30, 0, 400);
    auto* hRfit_b   = new TH1F("hRb","bad; R_fit [cm];", 30, 0, 400);

    int nGoodSeed = 0, nBadSeed = 0, nUKFconverged = 0, nUKFcore = 0;
    int nBadSeedButUKFcore = 0, nGoodSeedButUKFtail = 0;
    int dumped = 0;
    for (Long64_t i = 0; i < std::min({tS->GetEntries(), tT->GetEntries(), tU->GetEntries()}); ++i) {
        tS->GetEntry(i); tT->GetEntry(i); tU->GetEntry(i);
        if (trks->GetEntries() == 0) continue;
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        auto& hits = tr->GetHitArray();
        if (hits.size() < 10) continue;

        double R_fit = tr->GetGeoRadius();
        double th_seed = tr->GetGeoTheta();
        if (!std::isfinite(R_fit) || R_fit <= 0 || !std::isfinite(th_seed)) continue;
        double sinTh = std::sin(th_seed);
        if (std::abs(sinTh) < 0.1) continue;
        double pT_seed = 0.3 * B * R_fit / 100.0 * 1000.0;
        double p_seed  = pT_seed / sinTh;

        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mc->At(j);
            if (m->GetMotherId()==-1 && std::abs(m->GetPdgCode())==211) { pi = m; break; }
        }
        if (!pi) continue;
        double pMC = std::hypot(std::hypot(pi->GetPx(), pi->GetPy()), pi->GetPz()) * 1000;

        double xmn=1e9,xmx=-1e9,zmn=1e9,zmx=-1e9;
        for (auto& h : hits) { xmn=std::min(xmn,(double)h.GetX()); xmx=std::max(xmx,(double)h.GetX());
                                zmn=std::min(zmn,(double)h.GetZ()); zmx=std::max(zmx,(double)h.GetZ()); }
        double chord = std::hypot(xmx-xmn, zmx-zmn);

        double r_seed = p_seed / pMC - 1;
        bool seedGood = (std::abs(r_seed) < 0.10);
        if (seedGood) ++nGoodSeed; else ++nBadSeed;
        hSeed->Fill(r_seed);
        TH1F* hC = seedGood ? hChord_g : hChord_b;
        TH1F* hN = seedGood ? hNhits_g : hNhits_b;
        TH1F* hS = seedGood ? hSinTh_g : hSinTh_b;
        TH1F* hR = seedGood ? hRfit_g  : hRfit_b;
        hC->Fill(chord); hN->Fill(hits.size()); hS->Fill(std::abs(sinTh)); hR->Fill(R_fit);

        if (fitted->GetEntries() > 0) {
            auto* ft = (R3BGTPCFittedTrackData*)fitted->At(0);
            if (ft->IsConverged()) {
                ++nUKFconverged;
                const auto& kin = ft->GetKinematicsXtr();
                double pUKF = std::sqrt((kin.kineticEnergy + m_pi)*(kin.kineticEnergy + m_pi) - m_pi*m_pi);
                double r_ukf = pUKF / pMC - 1;
                hUKF->Fill(r_ukf);
                hSeedVsUKF->Fill(r_seed, r_ukf);
                bool ukfCore = std::abs(r_ukf) < 0.10;
                if (ukfCore) ++nUKFcore;
                if (!seedGood && ukfCore) ++nBadSeedButUKFcore;
                if (seedGood && !ukfCore) ++nGoodSeedButUKFtail;
                if (dumped < 10 && (!seedGood || !ukfCore)) {
                    ++dumped;
                    printf("evt %4lld: chord=%4.1f N=%3zu sinTh=%.2f R_fit=%5.1f  p_MC=%4.0f  p_seed=%4.0f (%+5.1f%%)  p_UKF=%4.0f (%+5.1f%%)  %s%s\n",
                           i, chord, hits.size(), std::abs(sinTh), R_fit, pMC, p_seed, 100*r_seed, pUKF, 100*r_ukf,
                           seedGood ? "seedOK" : "BADseed",
                           ukfCore ? " ukfOK" : " UKFtail");
                }
            }
        }
        ++N;
    }
    printf("\n=== 400 MeV/c box-gen, %d total seed-converged events ===\n", N);
    printf("  good seed (|r_seed| < 0.10): %d  (%.0f %%)\n", nGoodSeed, 100.*nGoodSeed/N);
    printf("  bad  seed (|r_seed| > 0.10): %d  (%.0f %%)\n", nBadSeed, 100.*nBadSeed/N);
    printf("  UKF converged: %d   (UKF core %d, %.0f %%)\n", nUKFconverged, nUKFcore, 100.*nUKFcore/std::max(1,nUKFconverged));
    printf("  bad-seed but UKF core: %d  (UKF rescued seed)\n", nBadSeedButUKFcore);
    printf("  good-seed but UKF tail: %d  (UKF broke a good seed)\n", nGoodSeedButUKFtail);

    TCanvas c("c","probe 400", 1600, 1000);
    c.Divide(3, 2, 0.006, 0.025);
    c.cd(1); gPad->SetLogy(); hSeed->SetLineColor(kBlue+1); hSeed->SetLineWidth(2); hSeed->Draw();
             hUKF->SetLineColor(kRed+1); hUKF->SetLineWidth(2); hUKF->Draw("SAME");
             auto* lg = new TLegend(0.55,0.7,0.88,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
             lg->AddEntry(hSeed,"seed","l"); lg->AddEntry(hUKF,"UKF","l"); lg->Draw();
    c.cd(2); hSeedVsUKF->Draw("colz");
             auto* l0a = new TLine(-0.1,-1,-0.1,1); l0a->SetLineStyle(2); l0a->Draw();
             auto* l0b = new TLine( 0.1,-1, 0.1,1); l0b->SetLineStyle(2); l0b->Draw();
             auto* l0c = new TLine(-0.5,-0.1,0.5,-0.1); l0c->SetLineStyle(2); l0c->Draw();
             auto* l0d = new TLine(-0.5, 0.1,0.5, 0.1); l0d->SetLineStyle(2); l0d->Draw();
    auto drawPair = [](TH1F* g, TH1F* b) {
        g->SetLineColor(kBlue+1); g->SetLineWidth(2);
        b->SetLineColor(kRed+1);  b->SetLineWidth(2);
        double mx = std::max(g->GetMaximum(), b->GetMaximum()) * 1.1;
        g->SetMaximum(mx); g->Draw(); b->Draw("SAME");
        auto* lg = new TLegend(0.55,0.75,0.88,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.035);
        lg->AddEntry(g,"good seed","l"); lg->AddEntry(b,"bad seed","l"); lg->Draw();
    };
    c.cd(3); drawPair(hChord_g, hChord_b);
    c.cd(4); drawPair(hNhits_g, hNhits_b);
    c.cd(5); drawPair(hSinTh_g, hSinTh_b);
    c.cd(6); drawPair(hRfit_g,  hRfit_b);
    c.SaveAs(wd + "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17/plots/probe_400.png");
}
