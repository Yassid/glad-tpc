/// @file draw_padplane_grid.C
/// @brief Grid of pad-plane events with MC truth + reco hits + clusters +
///        UKF smoothed track, each panel zoomed to its own active region.

void draw_padplane_grid(int p_MeV = 600, int nEvents = 6, int nCols = 3)
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gSystem->Load("libR3BOpenKF");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    const double offX_cm = 4.2;
    const double offZ_cm = 260.2;

    TString workDir = gSystem->Getenv("VMCWORKDIR");
    TFile fSim(TString::Format("%s/glad-tpc/macros/sim/Prototype/sim_p%d.root", workDir.Data(), p_MeV));
    TFile fLang(TString::Format("%s/glad-tpc/macros/proj/Prototype/lang_p%d.root", workDir.Data(), p_MeV));
    TFile fTrk(TString::Format("%s/glad-tpc/macros/tracking/output_tracking_p%d.root", workDir.Data(), p_MeV));
    TFile fUKF(TString::Format("%s/glad-tpc/macros/tracking/output_ukf_p%d.root", workDir.Data(), p_MeV));
    auto* tSim = (TTree*)fSim.Get("evt");
    auto* tLang = (TTree*)fLang.Get("evt");
    auto* tTrk = (TTree*)fTrk.Get("evt");
    auto* tUKF = (TTree*)fUKF.Get("evt");

    auto* gtpcPts = new TClonesArray("R3BGTPCPoint");
    auto* mcTracks = new TClonesArray("R3BMCTrack");
    auto* cal = new TClonesArray("R3BGTPCCalData");
    auto* trks = new TClonesArray("R3BGTPCTrackData");
    auto* fits = new TClonesArray("R3BGTPCFittedTrackData");
    tSim->SetBranchAddress("GTPCPoint", &gtpcPts);
    tSim->SetBranchAddress("MCTrack", &mcTracks);
    tLang->SetBranchAddress("GTPCCalData", &cal);
    tTrk->SetBranchAddress("GTPCTrackData", &trks);
    tUKF->SetBranchAddress("GTPCFittedTrackData", &fits);

    const int nColZ = 128, nRowX = 44;
    const double padSize_mm = 2.0;
    const double padPlaneZ_cm = nColZ * padSize_mm / 10.0;
    const double padPlaneX_cm = nRowX * padSize_mm / 10.0;

    // Pick events: converged fit, 30..200 pads, ≥2 MC primary pts
    std::vector<int> picks;
    const Long64_t N = std::min({ tSim->GetEntries(), tLang->GetEntries(), tTrk->GetEntries(), tUKF->GetEntries() });
    for (Long64_t i = 0; i < N && (int)picks.size() < nEvents; ++i) {
        tLang->GetEntry(i);
        tTrk->GetEntry(i);
        tUKF->GetEntry(i);
        tSim->GetEntry(i);
        if (cal->GetEntries() < 30 || cal->GetEntries() > 200) continue;
        if (fits->GetEntries() == 0 || trks->GetEntries() == 0) continue;
        int nPri = 0;
        for (int j = 0; j < gtpcPts->GetEntries(); ++j)
            if (((R3BGTPCPoint*)gtpcPts->At(j))->GetTrackID() == 0) ++nPri;
        if (nPri < 2) continue;
        picks.push_back(i);
    }
    if (picks.empty()) { std::cout << "[grid] no events\n"; return; }

    const int nRows = (picks.size() + nCols - 1) / nCols;
    auto* c = new TCanvas("cgr", Form("Pad plane grid p=%d MeV/c", p_MeV), 480 * nCols, 380 * nRows);
    c->Divide(nCols, nRows, 0.005, 0.02);

    for (size_t k = 0; k < picks.size(); ++k) {
        int i = picks[k];
        tSim->GetEntry(i);
        tLang->GetEntry(i);
        tTrk->GetEntry(i);
        tUKF->GetEntry(i);

        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mcTracks->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mcTracks->At(j);
            if (m->GetMotherId() == -1 && std::abs(m->GetPdgCode()) == 211) { pi = m; break; }
        }
        const double pmc = pi ? std::sqrt(pi->GetPx() * pi->GetPx() + pi->GetPy() * pi->GetPy() + pi->GetPz() * pi->GetPz()) * 1000 : 0;
        auto* fit = (R3BGTPCFittedTrackData*)fits->At(0);
        const double ke = fit->GetKinematics().kineticEnergy;
        const double pfit = std::sqrt(ke * ke + 2 * ke * 139.57);

        auto* hPad = new TH2F(Form("hPad_%zu", k),
                              Form("evt %d  p_fit=%.0f MeV/c  N_pads=%d;z (cm);x (cm)",
                                   i, pfit, cal->GetEntries()),
                              nColZ, 0., padPlaneZ_cm,
                              nRowX, 0., padPlaneX_cm);
        double zlo = 1e9, zhi = -1e9, xlo = 1e9, xhi = -1e9;
        for (int j = 0; j < cal->GetEntries(); ++j) {
            auto* d = (R3BGTPCCalData*)cal->At(j);
            const int pid = d->GetPadId();
            const int icol = pid / nRowX, irow = pid % nRowX;
            const double z_cm = (icol + 0.5) * padSize_mm / 10.0;
            const double x_cm = (irow + 0.5) * padSize_mm / 10.0;
            double q = 0.;
            for (auto a : d->GetADC()) q += a;
            hPad->SetBinContent(icol + 1, irow + 1, q);
            zlo = std::min(zlo, z_cm); zhi = std::max(zhi, z_cm);
            xlo = std::min(xlo, x_cm); xhi = std::max(xhi, x_cm);
        }
        hPad->GetXaxis()->SetRangeUser(std::max(0., zlo - 0.6), std::min(padPlaneZ_cm, zhi + 0.6));
        hPad->GetYaxis()->SetRangeUser(std::max(0., xlo - 0.6), std::min(padPlaneX_cm, xhi + 0.6));

        auto* gMC = new TGraph();   gMC->SetMarkerStyle(29); gMC->SetMarkerColor(kGreen + 2); gMC->SetMarkerSize(2.0);
        auto* gSec = new TGraph();  gSec->SetMarkerStyle(25); gSec->SetMarkerColor(kOrange + 1); gSec->SetMarkerSize(0.9);
        auto* gHit = new TGraph();  gHit->SetMarkerStyle(20); gHit->SetMarkerColor(kBlack); gHit->SetMarkerSize(0.45);
        auto* gCl = new TGraph();   gCl->SetMarkerStyle(24); gCl->SetMarkerColor(kRed); gCl->SetMarkerSize(1.4);
        auto* gFit = new TGraph();  gFit->SetMarkerStyle(34); gFit->SetMarkerColor(kMagenta + 1); gFit->SetMarkerSize(0.9);
        gFit->SetLineColor(kMagenta + 1); gFit->SetLineWidth(2);

        for (int j = 0; j < gtpcPts->GetEntries(); ++j) {
            auto* p = (R3BGTPCPoint*)gtpcPts->At(j);
            auto* g = (p->GetTrackID() == 0) ? gMC : gSec;
            g->SetPoint(g->GetN(), p->GetZ() - offZ_cm, p->GetX() - offX_cm);
        }
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        for (auto& h : tr->GetHitArray()) gHit->SetPoint(gHit->GetN(), h.GetZ(), h.GetX());
        for (auto& cc : *tr->GetHitClusterArray()) gCl->SetPoint(gCl->GetN(), cc.GetZ(), cc.GetX());
        for (auto& s : fit->GetSmoothedPositions()) gFit->SetPoint(gFit->GetN(), s.Z(), s.X());

        c->cd(k + 1);
        gPad->SetRightMargin(0.13);
        hPad->Draw("COLZ");
        gMC->Draw("P SAME");
        gSec->Draw("P SAME");
        gHit->Draw("P SAME");
        gCl->Draw("P SAME");
        gFit->Draw("LP SAME");
    }

    // Single shared legend in top-left of first pad
    c->cd(1);
    auto* leg = new TLegend(0.10, 0.62, 0.45, 0.88);
    leg->SetBorderSize(1); leg->SetFillStyle(1001); leg->SetFillColor(kWhite); leg->SetTextSize(0.035);
    auto* lh = new TGraph(1); lh->SetMarkerStyle(20); lh->SetMarkerColor(kBlack); lh->SetMarkerSize(0.8);
    auto* lm = new TGraph(1); lm->SetMarkerStyle(29); lm->SetMarkerColor(kGreen + 2); lm->SetMarkerSize(1.6);
    auto* lc = new TGraph(1); lc->SetMarkerStyle(24); lc->SetMarkerColor(kRed); lc->SetMarkerSize(1.4);
    auto* lf = new TGraph(1); lf->SetMarkerStyle(34); lf->SetMarkerColor(kMagenta + 1); lf->SetLineColor(kMagenta + 1); lf->SetLineWidth(2);
    leg->AddEntry(lh, "Reco hits", "p");
    leg->AddEntry(lm, "MC primary", "p");
    leg->AddEntry(lc, "Clusters", "p");
    leg->AddEntry(lf, "UKF smoothed", "lp");
    leg->Draw();

    c->SaveAs(Form("padplane_grid_p%d.png", p_MeV));
    std::cout << "Wrote padplane_grid_p" << p_MeV << ".png — " << picks.size() << " events\n";
}
