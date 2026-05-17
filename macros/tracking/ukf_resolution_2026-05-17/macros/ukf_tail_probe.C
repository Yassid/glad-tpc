// Diagnose UKF |p|-residual bimodality on long-chord good_evt events.
// Output: per-event p_UKF/p_MC histogram (1D), |residual| seed vs UKF, and
// four core-vs-tail population breakdowns (p_MC, chord, N_smoothed, N_hits).
// Reveals that the UKF's Gaussian-core σ_p/p number is misleading: it
// describes only the narrow peak at p_UKF/p_MC = 1, while ~40 % of events
// sit in a broad shelf at p_UKF/p_MC ∈ [0.1, 0.9]. The seed gives accurate
// momentum on essentially all events; the UKF drops half of them.
void ukf_tail_probe(const char* tag = "goodevt2k", double chordCutCm = 16.0)
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    TString wd = gSystem->Getenv("VMCWORKDIR");
    TFile fS(wd + "/glad-tpc/macros/sim/Prototype/sim_" + tag + ".root");
    TFile fT(wd + "/glad-tpc/macros/tracking/output_tracking_" + tag + ".root");
    TFile fU(wd + "/glad-tpc/macros/tracking/output_ukf_"      + tag + ".root");
    auto* tS = (TTree*)fS.Get("evt");
    auto* tT = (TTree*)fT.Get("evt");
    auto* tU = (TTree*)fU.Get("evt");
    if (!tS || !tT || !tU) { std::cerr << "missing input file\n"; return; }
    auto* mc     = new TClonesArray("R3BMCTrack");
    auto* trks   = new TClonesArray("R3BGTPCTrackData");
    auto* fitted = new TClonesArray("R3BGTPCFittedTrackData");
    tS->SetBranchAddress("MCTrack", &mc);
    tT->SetBranchAddress("GTPCTrackData", &trks);
    tU->SetBranchAddress("GTPCFittedTrackData", &fitted);
    const double m_pi = 139.57039;

    int nCore = 0, nTail = 0, nLongChord = 0;
    auto* hAll  = new TH1F("hAll", ";p_{UKF}/p_{MC};", 80, 0, 2);
    auto* hRat  = new TH1F("hRat", ";|p_{UKF}/p_{MC} - 1|;", 80, 0, 1);
    auto* hSeedRat = new TH1F("hSr", ";|p_{seed}/p_{MC} - 1|;", 80, 0, 1);
    auto* hCoreP   = new TH1F("hCp", "core;p_{MC} [MeV/c];",   12, 0, 1200);
    auto* hTailP   = new TH1F("htp", "tail;p_{MC} [MeV/c];",   12, 0, 1200);
    auto* hCoreCh  = new TH1F("hCc", "core;chord [cm];",       14, 14, 28);
    auto* hTailCh  = new TH1F("htc", "tail;chord [cm];",       14, 14, 28);
    auto* hCoreSP  = new TH1F("hCs", "core;N_{smoothed};",     30, 0, 60);
    auto* hTailSP  = new TH1F("hts", "tail;N_{smoothed};",     30, 0, 60);
    auto* hCoreHit = new TH1F("hCh", "core;N_{hits};",         25, 0, 250);
    auto* hTailHit = new TH1F("hth", "tail;N_{hits};",         25, 0, 250);

    for (Long64_t i = 0; i < std::min({tS->GetEntries(), tT->GetEntries(), tU->GetEntries()}); ++i) {
        tS->GetEntry(i); tT->GetEntry(i); tU->GetEntry(i);
        if (trks->GetEntries() == 0) continue;
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        auto& hits = tr->GetHitArray();
        if (hits.size() < 10) continue;
        double xmn=1e9,xmx=-1e9,zmn=1e9,zmx=-1e9;
        for (auto& h : hits) { xmn=std::min(xmn,(double)h.GetX()); xmx=std::max(xmx,(double)h.GetX());
                                zmn=std::min(zmn,(double)h.GetZ()); zmx=std::max(zmx,(double)h.GetZ()); }
        double chord = std::hypot(xmx-xmn, zmx-zmn);
        if (chord < chordCutCm) continue;
        ++nLongChord;
        if (fitted->GetEntries() == 0) continue;
        auto* ft = (R3BGTPCFittedTrackData*)fitted->At(0);
        if (!ft->IsConverged()) continue;
        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mc->At(j);
            if (m->GetMotherId()==-1 && std::abs(m->GetPdgCode())==211) { pi = m; break; }
        }
        if (!pi) continue;
        double pMC = std::hypot(std::hypot(pi->GetPx(), pi->GetPy()), pi->GetPz()) * 1000;
        const auto& kin = ft->GetKinematicsXtr();
        double pUKF = std::sqrt((kin.kineticEnergy + m_pi)*(kin.kineticEnergy + m_pi) - m_pi*m_pi);
        double r = pUKF / pMC;
        hAll->Fill(r);
        hRat->Fill(std::abs(r-1));

        double R_fit = tr->GetGeoRadius();
        double th_seed = tr->GetGeoTheta();
        if (std::isfinite(R_fit) && R_fit > 0 && std::isfinite(th_seed)) {
            double sinTh = std::max(0.1, std::abs(std::sin(th_seed)));
            double p_seed = 0.3 * 2.0 * R_fit / 100.0 * 1000.0 / sinTh;
            hSeedRat->Fill(std::abs(p_seed/pMC - 1));
        }

        bool core = (std::abs(r - 1) < 0.10);
        if (core) ++nCore; else ++nTail;
        TH1F* hP = core ? hCoreP : hTailP;
        TH1F* hC = core ? hCoreCh : hTailCh;
        TH1F* hS = core ? hCoreSP : hTailSP;
        TH1F* hH = core ? hCoreHit : hTailHit;
        hP->Fill(pMC);
        hC->Fill(chord);
        hS->Fill(ft->GetSmoothedPositions().size());
        hH->Fill(hits.size());
    }
    printf("=== Long-chord (>=%.1f cm) UKF outcome ===\n", chordCutCm);
    printf("  long-chord events total:                 %d\n", nLongChord);
    printf("  core (|p_UKF/p_MC - 1| < 0.10):          %d  (%.0f %%)\n", nCore, 100.*nCore/std::max(1,nCore+nTail));
    printf("  tail (|p_UKF/p_MC - 1| > 0.10):          %d  (%.0f %%)\n", nTail, 100.*nTail/std::max(1,nCore+nTail));
    printf("  seed mean |p_seed/p_MC - 1| = %.3f\n", hSeedRat->GetMean());
    printf("  UKF  mean |p_UKF /p_MC - 1| = %.3f\n", hRat->GetMean());

    TCanvas c("c", "ukf tail", 1600, 1100);
    c.Divide(3, 2, 0.006, 0.025);
    c.cd(1); gPad->SetLogy(); hAll->Draw();
    auto* lT = new TLine(1, 0.5, 1, hAll->GetMaximum()); lT->SetLineColor(kRed); lT->SetLineStyle(2); lT->Draw();
    c.cd(2); hRat->SetLineColor(kRed+1); hRat->SetLineWidth(2); hRat->Draw();
             hSeedRat->SetLineColor(kBlue+1); hSeedRat->SetLineWidth(2); hSeedRat->Draw("SAME");
             auto* lg = new TLegend(0.5, 0.7, 0.88, 0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
             lg->AddEntry(hSeedRat, "seed", "l"); lg->AddEntry(hRat, "UKF", "l"); lg->Draw();
    auto drawPair = [](TH1F* core, TH1F* tail) {
        core->SetLineColor(kBlue+1); core->SetLineWidth(2);
        tail->SetLineColor(kRed+1);  tail->SetLineWidth(2);
        double mx = std::max(core->GetMaximum(), tail->GetMaximum()) * 1.1;
        core->SetMaximum(mx); core->Draw(); tail->Draw("SAME");
        auto* lg = new TLegend(0.6, 0.78, 0.88, 0.88);
        lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.04);
        lg->AddEntry(core, "core", "l"); lg->AddEntry(tail, "tail", "l"); lg->Draw();
    };
    c.cd(3); drawPair(hCoreP, hTailP);
    c.cd(4); drawPair(hCoreCh, hTailCh);
    c.cd(5); drawPair(hCoreSP, hTailSP);
    c.cd(6); drawPair(hCoreHit, hTailHit);
    c.SaveAs(wd + "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17/plots/ukf_tail_probe.png");
}
