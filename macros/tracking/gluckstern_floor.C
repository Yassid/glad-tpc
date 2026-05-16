// Per-event Gluckstern floor for σ(p)/p on the circle fit.
//
//   σ(1/R) = σ_xy · sqrt(720 / (N+4)) / L²        (Gluckstern 1963, ideal)
//   σ_R / R = R · σ(1/R) = σ_xy · sqrt(720/N) · R / L²
//
// Compare against the measured σ_R/R on each chord bin to see how close we
// run to the theoretical floor. If we're 1-2× above, the fit is fine and
// improvements need MORE data (longer chord, vertex constraint, etc.). If
// we're 5×+ above, the fit / clustering is leaving info on the table.
void gluckstern_floor(double sigma_xy_mm = 1.0)
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    TString wd = gSystem->Getenv("VMCWORKDIR");
    TFile fSim(wd + "/glad-tpc/macros/sim/Prototype/sim_goodevt.root");
    TFile fTrk("output_tracking_goodevt.root");
    auto* tS = (TTree*)fSim.Get("evt");
    auto* tT = (TTree*)fTrk.Get("evt");
    auto* mc = new TClonesArray("R3BMCTrack");
    auto* trks = new TClonesArray("R3BGTPCTrackData");
    tS->SetBranchAddress("MCTrack", &mc);
    tT->SetBranchAddress("GTPCTrackData", &trks);
    const double B = 2.0;

    struct Stats { int n=0; double sigSum=0, sigSum2=0, deltaSum=0; };
    Stats binsLong, binsMid, binsShort;
    auto* gFloorVsChord = new TGraph();
    auto* gMeasVsChord  = new TGraph();
    int npts = 0;

    for (Long64_t i = 0; i < std::min(tS->GetEntries(), tT->GetEntries()); ++i) {
        tS->GetEntry(i); tT->GetEntry(i);
        if (trks->GetEntries() == 0) continue;
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        auto& hits = tr->GetHitArray();
        if (hits.size() < 10) continue;
        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mc->At(j);
            if (m->GetMotherId()==-1 && std::abs(m->GetPdgCode())==211) { pi = m; break; }
        }
        if (!pi) continue;
        double pT = std::sqrt(pi->GetPx()*pi->GetPx() + pi->GetPz()*pi->GetPz())*1000;
        double R_truth_cm = pT / (3.0 * B);
        double R_fit = tr->GetGeoRadius();
        if (!std::isfinite(R_fit) || R_fit <= 0) continue;
        double xmn=1e9,xmx=-1e9,zmn=1e9,zmx=-1e9;
        for (auto& h : hits) { xmn=std::min(xmn,(double)h.GetX()); xmx=std::max(xmx,(double)h.GetX());
                               zmn=std::min(zmn,(double)h.GetZ()); zmx=std::max(zmx,(double)h.GetZ()); }
        double L_cm = std::sqrt((xmx-xmn)*(xmx-xmn)+(zmx-zmn)*(zmx-zmn));
        int N = hits.size();

        // Gluckstern σ(1/R) and σ_R/R
        double sig_xy_cm = sigma_xy_mm / 10.0;
        double sig_inv_R = sig_xy_cm * std::sqrt(720.0 / (N + 4)) / (L_cm * L_cm);
        double sigR_floor = R_truth_cm * R_truth_cm * sig_inv_R;
        double floor_frac = sigR_floor / R_truth_cm;

        double deviation = std::abs(R_fit / R_truth_cm - 1.0);

        gFloorVsChord->SetPoint(npts, L_cm, 100*floor_frac);
        gMeasVsChord->SetPoint(npts, L_cm, 100*deviation);
        ++npts;

        Stats* b = nullptr;
        if (L_cm >= 16) b = &binsLong;
        else if (L_cm >= 12) b = &binsMid;
        else b = &binsShort;
        b->n++;
        b->sigSum  += floor_frac;
        b->sigSum2 += floor_frac*floor_frac;
        b->deltaSum += deviation;
    }
    auto report = [](const char* lbl, const Stats& s) {
        if (s.n == 0) return;
        double mF = s.sigSum/s.n;
        double sF = std::sqrt(std::max(0., s.sigSum2/s.n - mF*mF));
        double mD = s.deltaSum/s.n;
        printf("  %-12s  N=%3d   Gluckstern floor median %5.2f%% (rms %.1f%%)   |R_fit-R_tru|/R_tru mean = %5.1f%%\n",
               lbl, s.n, 100*mF, 100*sF, 100*mD);
    };
    printf("=== Gluckstern floor (σ_xy = %.1f mm) ===\n", sigma_xy_mm);
    report("short (<12)", binsShort);
    report("mid (12-16)", binsMid);
    report("long (>=16)", binsLong);

    TCanvas c("c","",1100,500);
    c.cd();
    gFloorVsChord->SetMarkerStyle(20); gFloorVsChord->SetMarkerSize(0.5); gFloorVsChord->SetMarkerColor(kBlue+1);
    gMeasVsChord ->SetMarkerStyle(24); gMeasVsChord ->SetMarkerSize(0.5); gMeasVsChord ->SetMarkerColor(kRed+1);
    auto* h = new TH1F("h",";chord L [cm];|R_{fit}-R_{tru}|/R_{tru} (red), Gluckstern floor (blue) [%]", 30, 0, 30);
    h->SetMinimum(0); h->SetMaximum(200);
    h->Draw();
    gFloorVsChord->Draw("P SAME");
    gMeasVsChord->Draw("P SAME");
    auto* lg = new TLegend(0.55, 0.7, 0.95, 0.88);
    lg->AddEntry(gFloorVsChord, "Gluckstern floor", "p");
    lg->AddEntry(gMeasVsChord, "measured |R_{fit}/R_{tru}-1|", "p");
    lg->Draw();
    c.SaveAs("gluckstern_floor.png");
}
