// 400 MeV/c acceptance-edge diagnostic.
// Three panels showing why the box-gen chamber acceptance falls off at low p:
//   1. Pipeline funnel: events surviving each stage (gen / reco / seed / UKF),
//      per momentum point.
//   2. Chord distribution per momentum — short chord at low p is the cause.
//   3. R_seed distribution per momentum — R-cap pileup only at the edge.
void plot_accept_400()
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    TString wd = gSystem->Getenv("VMCWORKDIR");

    // Momentum points to overlay (subset of the scan)
    std::vector<int> P = { 400, 600, 800, 1000, 1200 };
    int colors[] = { kRed+1, kOrange+1, kGreen+2, kAzure+2, kViolet+1 };

    // Acceptance funnel: counts at each stage per momentum
    int nGen = 500;  // events generated per point (NEVT from scan_p.sh)
    std::vector<int> nReco(P.size(), 0), nSeed(P.size(), 0), nUKF(P.size(), 0);

    // Chord and R_seed distributions
    std::vector<TH1F*> hChord(P.size()), hRseed(P.size()), hResid(P.size());
    for (size_t k = 0; k < P.size(); ++k) {
        hChord[k] = new TH1F(Form("hC_%d", P[k]), ";chord (x,z) [cm];events", 25, 0, 25);
        hRseed[k] = new TH1F(Form("hR_%d", P[k]), ";R_{fit} [cm];events", 40, 0, 2050);
        hResid[k] = new TH1F(Form("hRes_%d", P[k]), ";p_{UKF}/p_{MC} - 1;events", 80, -1, 1);
        hChord[k]->SetLineColor(colors[k]); hChord[k]->SetLineWidth(2);
        hRseed[k]->SetLineColor(colors[k]); hRseed[k]->SetLineWidth(2);
        hResid[k]->SetLineColor(colors[k]); hResid[k]->SetLineWidth(2);
    }

    const double m_pi = 139.57039;
    for (size_t k = 0; k < P.size(); ++k) {
        TString suf = Form("_p%d", P[k]);
        TFile fS(wd + "/glad-tpc/macros/sim/Prototype/sim" + suf + ".root");
        TFile fR(wd + "/glad-tpc/macros/reco/output_reco" + suf + ".root");
        TFile fT(wd + "/glad-tpc/macros/tracking/output_tracking" + suf + ".root");
        TFile fU(wd + "/glad-tpc/macros/tracking/output_ukf"      + suf + ".root");
        auto* tS = (TTree*)fS.Get("evt");
        auto* tR = (TTree*)fR.Get("evt");
        auto* tT = (TTree*)fT.Get("evt");
        auto* tU = (TTree*)fU.Get("evt");
        if (!tS || !tT || !tU) { printf("missing files for p=%d\n", P[k]); continue; }

        auto* mc     = new TClonesArray("R3BMCTrack");
        auto* recoHits = new TClonesArray("R3BGTPCHitData");
        auto* trks   = new TClonesArray("R3BGTPCTrackData");
        auto* fitted = new TClonesArray("R3BGTPCFittedTrackData");
        tS->SetBranchAddress("MCTrack", &mc);
        if (tR) tR->SetBranchAddress("GTPCHitData", &recoHits);
        tT->SetBranchAddress("GTPCTrackData", &trks);
        tU->SetBranchAddress("GTPCFittedTrackData", &fitted);

        Long64_t nE = std::min(tS->GetEntries(), tT->GetEntries());
        if (tU) nE = std::min(nE, tU->GetEntries());
        for (Long64_t i = 0; i < nE; ++i) {
            tS->GetEntry(i); tT->GetEntry(i);
            if (tR) tR->GetEntry(i);
            if (tU) tU->GetEntry(i);

            R3BMCTrack* pi = nullptr;
            for (int j = 0; j < mc->GetEntries(); ++j) {
                auto* m = (R3BMCTrack*)mc->At(j);
                if (m->GetMotherId()==-1 && std::abs(m->GetPdgCode())==211) { pi = m; break; }
            }
            if (!pi) continue;
            double pMC = std::hypot(std::hypot(pi->GetPx(), pi->GetPy()), pi->GetPz()) * 1000;
            if (tR && recoHits->GetEntries() > 5) ++nReco[k];
            if (trks->GetEntries() == 0) continue;
            auto* tr = (R3BGTPCTrackData*)trks->At(0);
            auto& hits = tr->GetHitArray();
            if (hits.size() < 10) continue;
            double R_fit = tr->GetGeoRadius();
            if (!std::isfinite(R_fit) || R_fit <= 0) continue;
            ++nSeed[k];

            double xmn=1e9,xmx=-1e9,zmn=1e9,zmx=-1e9;
            for (auto& h : hits) { xmn=std::min(xmn,(double)h.GetX()); xmx=std::max(xmx,(double)h.GetX());
                                    zmn=std::min(zmn,(double)h.GetZ()); zmx=std::max(zmx,(double)h.GetZ()); }
            double chord = std::hypot(xmx-xmn, zmx-zmn);
            hChord[k]->Fill(chord);
            hRseed[k]->Fill(std::min(2049.0, R_fit));

            if (fitted->GetEntries() > 0) {
                auto* ft = (R3BGTPCFittedTrackData*)fitted->At(0);
                if (ft->IsConverged()) {
                    ++nUKF[k];
                    const auto& kin = ft->GetKinematicsXtr();
                    if (std::isfinite(kin.kineticEnergy) && kin.kineticEnergy > 0) {
                        double pUKF = std::sqrt((kin.kineticEnergy+m_pi)*(kin.kineticEnergy+m_pi) - m_pi*m_pi);
                        hResid[k]->Fill(pUKF/pMC - 1);
                    }
                }
            }
        }
    }

    printf("\n=== Acceptance funnel ===\n");
    printf("p (MeV/c)  gen   reco   seed   UKF   (reco%%  seed%%  UKF%%)\n");
    for (size_t k = 0; k < P.size(); ++k) {
        printf("  %4d     %3d   %3d    %3d    %3d   ( %4.1f  %4.1f  %4.1f )\n",
               P[k], nGen, nReco[k], nSeed[k], nUKF[k],
               100.*nReco[k]/nGen, 100.*nSeed[k]/nGen, 100.*nUKF[k]/nGen);
    }

    TCanvas c("c", "400 MeV/c acceptance edge", 1800, 600);
    c.Divide(3, 1, 0.008, 0.020);

    // Panel 1: acceptance funnel as fractions
    c.cd(1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    auto* hReco = new TH1F("hRecoF","",P.size(), 0, P.size());
    auto* hSeedF= new TH1F("hSeedF","",P.size(), 0, P.size());
    auto* hUKFF = new TH1F("hUKFF", "",P.size(), 0, P.size());
    for (size_t k = 0; k < P.size(); ++k) {
        hReco ->SetBinContent(k+1, 100.0*nReco[k]/nGen);
        hSeedF->SetBinContent(k+1, 100.0*nSeed[k]/nGen);
        hUKFF ->SetBinContent(k+1, 100.0*nUKF[k]/nGen);
        hReco ->GetXaxis()->SetBinLabel(k+1, Form("%d", P[k]));
        hSeedF->GetXaxis()->SetBinLabel(k+1, Form("%d", P[k]));
        hUKFF ->GetXaxis()->SetBinLabel(k+1, Form("%d", P[k]));
    }
    hReco ->SetLineColor(kGray+2); hReco ->SetLineWidth(3); hReco ->SetMarkerStyle(24); hReco ->SetMarkerSize(1.5);
    hSeedF->SetLineColor(kBlue+1); hSeedF->SetLineWidth(3); hSeedF->SetMarkerStyle(20); hSeedF->SetMarkerSize(1.5); hSeedF->SetMarkerColor(kBlue+1);
    hUKFF ->SetLineColor(kRed+1);  hUKFF ->SetLineWidth(3); hUKFF ->SetMarkerStyle(21); hUKFF ->SetMarkerSize(1.5); hUKFF ->SetMarkerColor(kRed+1);
    hReco->SetMaximum(110); hReco->SetMinimum(0);
    hReco->GetXaxis()->SetTitle("p_{MC} [MeV/c]"); hReco->GetYaxis()->SetTitle("survival fraction [%]");
    hReco->GetXaxis()->SetLabelSize(0.05);
    hReco->Draw("P"); hSeedF->Draw("P SAME"); hUKFF->Draw("P SAME");
    auto* lg = new TLegend(0.18, 0.18, 0.58, 0.38);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.035);
    lg->AddEntry(hReco,  "reco hits (>=5)", "p");
    lg->AddEntry(hSeedF, "seed fit (R>0, hits>=10)", "p");
    lg->AddEntry(hUKFF,  "UKF converged", "p");
    lg->Draw();
    TLatex t; t.SetTextSize(0.028); t.SetTextColor(kGray+2);
    t.DrawLatex(0.1, 5, "Acceptance falls at 400: short chord + R-cap seeds");

    // Panel 2: chord distribution per p
    c.cd(2); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14);
    double maxC = 0;
    for (auto* h : hChord) maxC = std::max(maxC, h->GetMaximum());
    hChord[0]->SetMaximum(maxC * 1.15);
    hChord[0]->Draw("HIST");
    for (size_t k = 1; k < P.size(); ++k) hChord[k]->Draw("HIST SAME");
    auto* l16 = new TLine(16, 0, 16, maxC*1.15);
    l16->SetLineStyle(2); l16->SetLineColor(kGray+2); l16->Draw();
    t.SetTextColor(kGray+2); t.DrawLatex(16.4, maxC*1.0, "long-chord cut");
    auto* lg2 = new TLegend(0.62, 0.55, 0.92, 0.88);
    lg2->SetBorderSize(0); lg2->SetFillStyle(0); lg2->SetTextSize(0.035);
    for (size_t k = 0; k < P.size(); ++k) lg2->AddEntry(hChord[k], Form("%d MeV/c", P[k]), "l");
    lg2->Draw();

    // Panel 3: R_seed distribution per p
    c.cd(3); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetLogy();
    double maxR = 0;
    for (auto* h : hRseed) maxR = std::max(maxR, h->GetMaximum());
    hRseed[0]->SetMaximum(maxR * 5); hRseed[0]->SetMinimum(0.5);
    hRseed[0]->Draw("HIST");
    for (size_t k = 1; k < P.size(); ++k) hRseed[k]->Draw("HIST SAME");
    auto* l500 = new TLine(500, 0.5, 500, maxR*5);
    l500->SetLineStyle(2); l500->SetLineColor(kGray+2); l500->Draw();
    t.DrawLatex(560, maxR*1.5, "SetMaxSeedRadius");
    auto* l2000 = new TLine(2000, 0.5, 2000, maxR*5);
    l2000->SetLineStyle(3); l2000->SetLineColor(kBlack); l2000->Draw();
    t.SetTextColor(kBlack); t.DrawLatex(1550, maxR*0.05, "GN R cap (20 m)");
    auto* lg3 = new TLegend(0.62, 0.55, 0.92, 0.88);
    lg3->SetBorderSize(0); lg3->SetFillStyle(0); lg3->SetTextSize(0.035);
    for (size_t k = 0; k < P.size(); ++k) lg3->AddEntry(hRseed[k], Form("%d MeV/c", P[k]), "l");
    lg3->Draw();

    c.SaveAs(wd + "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17/plots/accept_400.png");

    // Bonus: residual overlay
    TCanvas c2("c2", "residual per p", 900, 600);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetLogy();
    double maxRes = 0;
    for (auto* h : hResid) maxRes = std::max(maxRes, h->GetMaximum());
    if (maxRes > 0) {
        hResid[0]->SetMaximum(maxRes * 5); hResid[0]->SetMinimum(0.5);
        hResid[0]->Draw("HIST");
        for (size_t k = 1; k < P.size(); ++k) hResid[k]->Draw("HIST SAME");
        auto* l0 = new TLine(0, 0.5, 0, maxRes*5); l0->SetLineStyle(2); l0->SetLineColor(kGray+2); l0->Draw();
        auto* lg4 = new TLegend(0.62, 0.55, 0.92, 0.88);
        lg4->SetBorderSize(0); lg4->SetFillStyle(0); lg4->SetTextSize(0.035);
        for (size_t k = 0; k < P.size(); ++k) lg4->AddEntry(hResid[k], Form("%d MeV/c", P[k]), "l");
        lg4->Draw();
    }
    c2.SaveAs(wd + "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17/plots/residual_overlay_p.png");
}
