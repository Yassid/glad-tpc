// Side-by-side robustness of three resolution estimators on the same
// long-chord (>=16 cm) good_evt 2k UKF residual distribution:
//   sigma_core  -- Gaussian fit to the core in [-0.4, 0.4]
//   sigma_RMS3  -- RMS within +/- 3*sigma_core of the mean (truncated)
//   sigma_q     -- (q84 - q16) / 2  (quantile half-width, robust)
// Produces plots/robustness_check.png with the three numbers annotated.
void robustness_check(TString ukfFile = "output_ukf_goodevt2k.root",
                      TString simFile = "Prototype/sim_goodevt2k.root",
                      TString trkFile = "output_tracking_goodevt2k.root")
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

    TString wd = gSystem->Getenv("VMCWORKDIR");
    TFile fS(wd + "/glad-tpc/macros/sim/" + simFile);
    TFile fU(wd + "/glad-tpc/macros/tracking/" + ukfFile);
    TFile fT(wd + "/glad-tpc/macros/tracking/" + trkFile);

    auto* tS = (TTree*)fS.Get("evt");
    auto* tU = (TTree*)fU.Get("evt");
    auto* tT = (TTree*)fT.Get("evt");
    auto* mc     = new TClonesArray("R3BMCTrack");
    auto* fitted = new TClonesArray("R3BGTPCFittedTrackData");
    auto* trks   = new TClonesArray("R3BGTPCTrackData");
    tS->SetBranchAddress("MCTrack", &mc);
    tU->SetBranchAddress("GTPCFittedTrackData", &fitted);
    tT->SetBranchAddress("GTPCTrackData", &trks);

    const double m_pi = 139.57039;
    auto p_from_KE = [&](double ke){ return std::sqrt((ke+m_pi)*(ke+m_pi)-m_pi*m_pi); };

    auto* h = new TH1F("hLong",
        ";p_{UKF}/p_{MC} - 1;events / 0.025",
        80, -1.0, 1.0);

    Long64_t nE = std::min({tS->GetEntries(), tU->GetEntries(), tT->GetEntries()});
    int nL = 0;
    for (Long64_t i = 0; i < nE; ++i) {
        tS->GetEntry(i); tU->GetEntry(i); tT->GetEntry(i);
        if (fitted->GetEntries()==0 || trks->GetEntries()==0) continue;
        auto* ft = (R3BGTPCFittedTrackData*)fitted->At(0);
        if (!ft->IsConverged()) continue;
        const auto& kin = ft->GetKinematicsXtr();
        if (!std::isfinite(kin.kineticEnergy) || kin.kineticEnergy <= 0) continue;
        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* mt = (R3BMCTrack*)mc->At(j);
            if (mt->GetMotherId()==-1 && mt->GetPdgCode()==-211) { pi = mt; break; }
        }
        if (!pi) continue;
        double pMC = std::sqrt(pi->GetPx()*pi->GetPx()+pi->GetPy()*pi->GetPy()+pi->GetPz()*pi->GetPz())*1000;
        if (pMC <= 0) continue;
        auto& hits = ((R3BGTPCTrackData*)trks->At(0))->GetHitArray();
        if (hits.size() < 10) continue;
        double xmn=1e9,xmx=-1e9,zmn=1e9,zmx=-1e9;
        for (auto& q : hits) {
            xmn=std::min(xmn,(double)q.GetX()); xmx=std::max(xmx,(double)q.GetX());
            zmn=std::min(zmn,(double)q.GetZ()); zmx=std::max(zmx,(double)q.GetZ());
        }
        double chord = std::sqrt((xmx-xmn)*(xmx-xmn)+(zmx-zmn)*(zmx-zmn));
        if (chord < 16) continue;
        double pFit = p_from_KE(kin.kineticEnergy);
        h->Fill(pFit/pMC - 1);
        ++nL;
    }
    printf("Long-chord entries: %d\n", nL);

    // 1. Gaussian core sigma in [-0.4, 0.4]
    h->Fit("gaus", "Q", "", -0.4, 0.4);
    auto* fG = h->GetFunction("gaus");
    double mu_g = fG->GetParameter(1);
    double s_g  = fG->GetParameter(2);

    // 2. Truncated RMS within +/- 3 sigma_core of the Gaussian mean
    h->GetXaxis()->SetRangeUser(mu_g - 3*s_g, mu_g + 3*s_g);
    double s_trunc = h->GetRMS();
    double mu_trunc = h->GetMean();
    h->GetXaxis()->SetRange(0, 0);  // restore full range

    // 3. Quantile half-width (q84 - q16) / 2 on the full distribution
    double probs[3] = { 0.16, 0.50, 0.84 };
    double q[3];
    h->GetQuantiles(3, q, probs);
    double s_q = 0.5 * (q[2] - q[0]);

    printf("=== Resolution estimators on long-chord good_evt 2k (N=%d) ===\n", nL);
    printf("  Gaussian core (fit [-0.4,0.4]):  mu = %+.4f   sigma = %.4f   ==> %.2f%%\n", mu_g, s_g, 100*s_g);
    printf("  Truncated RMS (+/- 3 sigma):     mu = %+.4f   sigma = %.4f   ==> %.2f%%\n", mu_trunc, s_trunc, 100*s_trunc);
    printf("  Quantile half-width (q84-q16)/2: median = %+.4f sigma = %.4f  ==> %.2f%%\n", q[1], s_q, 100*s_q);

    // Draw
    auto* c = new TCanvas("c_robust","robustness", 1200, 700);
    c->SetLeftMargin(0.12); c->SetRightMargin(0.04);
    c->SetBottomMargin(0.13); c->SetTopMargin(0.05);
    h->SetLineColor(kBlack); h->SetLineWidth(2);
    h->GetXaxis()->SetRangeUser(-0.4, 0.4);
    h->GetXaxis()->SetTitleSize(0.05);
    h->GetYaxis()->SetTitleSize(0.05);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->Draw("HIST");
    if (auto* f = h->GetFunction("gaus")) {
        f->SetLineColor(kRed); f->SetLineWidth(3);
        f->Draw("SAME");
    }

    // Vertical bands at +/- each estimator
    auto draw_band = [&](double sigma, int color, double yfrac){
        double y = yfrac * h->GetMaximum();
        auto* lP = new TLine(mu_g + sigma, 0, mu_g + sigma, y); lP->SetLineColor(color); lP->SetLineWidth(2); lP->SetLineStyle(2); lP->Draw();
        auto* lM = new TLine(mu_g - sigma, 0, mu_g - sigma, y); lM->SetLineColor(color); lM->SetLineWidth(2); lM->SetLineStyle(2); lM->Draw();
    };
    draw_band(s_g,     kRed,      1.00);
    draw_band(s_trunc, kBlue+1,   0.70);
    draw_band(s_q,     kGreen+2,  0.45);

    auto* leg = new TLegend(0.55, 0.62, 0.96, 0.93);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.040);
    leg->AddEntry((TObject*)nullptr, Form("N = %d  (chord #geq 16 cm)", nL), "");
    leg->AddEntry((TObject*)nullptr, Form("Gaussian core    #sigma = %.2f %%", 100*s_g), "");
    leg->AddEntry((TObject*)nullptr, Form("Truncated RMS    #sigma = %.2f %%", 100*s_trunc), "");
    leg->AddEntry((TObject*)nullptr, Form("(q_{84}#minusq_{16})/2  = %.2f %%", 100*s_q), "");
    leg->Draw();

    TString outDir = wd + "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17/plots/";
    c->SaveAs(outDir + "robustness_check.png");
    printf("Wrote %srobustness_check.png\n", outDir.Data());
}
