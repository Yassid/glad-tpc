/// @file draw_tracks.C
/// @brief Display N representative tracks: MC truth points, reconstructed
///        clusters, UKF smoothed positions — 3D + three 2D projections.

void draw_one(int evt, TTree* tMC, TTree* tUKF, TClonesArray* mcTracks,
              TClonesArray* gtpcPoints, TClonesArray* fittedTracks,
              TClonesArray* fitTrackData, TCanvas* c, int padBase);

void draw_tracks(int p_MeV = 800, int nTracks = 6)
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gSystem->Load("libR3BOpenKF");
    gStyle->SetOptStat(0);

    TString workDir = gSystem->Getenv("VMCWORKDIR");
    TString mcPath = TString::Format("%s/glad-tpc/macros/sim/Prototype/sim_p%d.root", workDir.Data(), p_MeV);
    TString trkPath = TString::Format("%s/glad-tpc/macros/tracking/output_tracking_p%d.root", workDir.Data(), p_MeV);
    TString ukfPath = TString::Format("%s/glad-tpc/macros/tracking/output_ukf_p%d.root", workDir.Data(), p_MeV);

    TFile fMC(mcPath);
    TFile fTrk(trkPath);
    TFile fUKF(ukfPath);
    auto* tMC = (TTree*)fMC.Get("evt");
    auto* tTrk = (TTree*)fTrk.Get("evt");
    auto* tUKF = (TTree*)fUKF.Get("evt");

    auto* mcTracks = new TClonesArray("R3BMCTrack");
    auto* gtpcPoints = new TClonesArray("R3BGTPCPoint");
    auto* gtpcTracks = new TClonesArray("R3BGTPCTrackData");
    auto* fittedTracks = new TClonesArray("R3BGTPCFittedTrackData");
    tMC->SetBranchAddress("MCTrack", &mcTracks);
    tMC->SetBranchAddress("GTPCPoint", &gtpcPoints);
    tTrk->SetBranchAddress("GTPCTrackData", &gtpcTracks);
    tUKF->SetBranchAddress("GTPCFittedTrackData", &fittedTracks);

    // Pick events that have ≥40 clusters and a converged fit, and that are
    // visually distinct (spread the choice across the run).
    std::vector<int> picks;
    const Long64_t N = std::min({ tMC->GetEntries(), tTrk->GetEntries(), tUKF->GetEntries() });
    const int stride = std::max<int>(1, N / (nTracks * 3));
    for (Long64_t i = 0; i < N && (int)picks.size() < nTracks; i += stride) {
        tTrk->GetEntry(i);
        tUKF->GetEntry(i);
        if (gtpcTracks->GetEntries() == 0 || fittedTracks->GetEntries() == 0)
            continue;
        auto* tr = (R3BGTPCTrackData*)gtpcTracks->At(0);
        if (tr->GetHitArray().size() < 50)
            continue;
        picks.push_back(i);
    }
    if (picks.empty()) {
        std::cout << "[draw_tracks] no suitable events\n";
        return;
    }

    auto* c3d = new TCanvas("c3d", Form("Tracks p=%d MeV/c — 3D", p_MeV), 1200, 800);
    c3d->Divide(3, 2);
    auto* cxy = new TCanvas("cxy", Form("Tracks p=%d MeV/c — projections", p_MeV), 1600, 900);
    cxy->Divide(3, picks.size());

    for (size_t k = 0; k < picks.size(); ++k) {
        int i = picks[k];
        tMC->GetEntry(i);
        tTrk->GetEntry(i);
        tUKF->GetEntry(i);

        auto* tr = (R3BGTPCTrackData*)gtpcTracks->At(0);
        auto& clusters = *tr->GetHitClusterArray();
        auto& rawHits = tr->GetHitArray();
        auto* fit = (R3BGTPCFittedTrackData*)fittedTracks->At(0);
        const auto& smoothed = fit->GetSmoothedPositions();

        // Graphs:
        //   gMC*  = MC truth GTPCPoint  (gray, small)
        //   gHT*  = reconstructed hits  (light blue, medium)
        //   gCL*  = TripletClust clusters (blue, larger)
        //   gFT*  = UKF smoothed track  (red line)
        auto* gMC3 = new TGraph2D(); gMC3->SetMarkerStyle(29); gMC3->SetMarkerColor(kGreen + 2); gMC3->SetMarkerSize(1.2);
        auto* gHT3 = new TGraph2D(); gHT3->SetMarkerStyle(20); gHT3->SetMarkerColor(kAzure + 1); gHT3->SetMarkerSize(0.5);
        auto* gCL3 = new TGraph2D(); gCL3->SetMarkerStyle(24); gCL3->SetMarkerColor(kBlue + 2); gCL3->SetMarkerSize(1.2);
        auto* gFT3 = new TGraph2D(); gFT3->SetMarkerStyle(34); gFT3->SetMarkerColor(kRed + 1); gFT3->SetMarkerSize(1.4);

        auto* gMCxy = new TGraph(); auto* gHTxy = new TGraph(); auto* gCLxy = new TGraph(); auto* gFTxy = new TGraph();
        auto* gMCxz = new TGraph(); auto* gHTxz = new TGraph(); auto* gCLxz = new TGraph(); auto* gFTxz = new TGraph();
        auto* gMCyz = new TGraph(); auto* gHTyz = new TGraph(); auto* gCLyz = new TGraph(); auto* gFTyz = new TGraph();
        for (auto* g : { gMCxy, gMCxz, gMCyz }) { g->SetMarkerStyle(29); g->SetMarkerColor(kGreen + 2); g->SetMarkerSize(1.0); }
        for (auto* g : { gHTxy, gHTxz, gHTyz }) { g->SetMarkerStyle(20); g->SetMarkerColor(kAzure + 1); g->SetMarkerSize(0.5); }
        for (auto* g : { gCLxy, gCLxz, gCLyz }) { g->SetMarkerStyle(24); g->SetMarkerColor(kBlue + 2); g->SetMarkerSize(1.0); }
        for (auto* g : { gFTxy, gFTxz, gFTyz }) { g->SetMarkerStyle(34); g->SetMarkerColor(kRed + 1); g->SetMarkerSize(1.2); g->SetLineColor(kRed + 1); g->SetLineWidth(2); }

        for (int j = 0; j < gtpcPoints->GetEntries(); ++j) {
            auto* p = (R3BGTPCPoint*)gtpcPoints->At(j);
            const double x = p->GetX(), y = p->GetY(), z = p->GetZ();
            gMC3->SetPoint(gMC3->GetN(), x, y, z);
            gMCxy->SetPoint(gMCxy->GetN(), x, y);
            gMCxz->SetPoint(gMCxz->GetN(), x, z);
            gMCyz->SetPoint(gMCyz->GetN(), z, y);
        }
        for (auto& h : rawHits) {
            gHT3->SetPoint(gHT3->GetN(), h.GetX(), h.GetY(), h.GetZ());
            gHTxy->SetPoint(gHTxy->GetN(), h.GetX(), h.GetY());
            gHTxz->SetPoint(gHTxz->GetN(), h.GetX(), h.GetZ());
            gHTyz->SetPoint(gHTyz->GetN(), h.GetZ(), h.GetY());
        }
        for (auto& h : clusters) {
            gCL3->SetPoint(gCL3->GetN(), h.GetX(), h.GetY(), h.GetZ());
            gCLxy->SetPoint(gCLxy->GetN(), h.GetX(), h.GetY());
            gCLxz->SetPoint(gCLxz->GetN(), h.GetX(), h.GetZ());
            gCLyz->SetPoint(gCLyz->GetN(), h.GetZ(), h.GetY());
        }
        // UKF smoothed positions are stored in input unit (cm in R3B).
        for (auto& s : smoothed) {
            gFT3->SetPoint(gFT3->GetN(), s.X(), s.Y(), s.Z());
            gFTxy->SetPoint(gFTxy->GetN(), s.X(), s.Y());
            gFTxz->SetPoint(gFTxz->GetN(), s.X(), s.Z());
            gFTyz->SetPoint(gFTyz->GetN(), s.Z(), s.Y());
        }

        const double ke_fit = fit->GetKinematics().kineticEnergy;
        const double p_fit = std::sqrt(ke_fit * ke_fit + 2 * ke_fit * 139.57);
        TString hdr = Form("evt %d  N_hit=%lu  N_cl=%lu  p_fit=%.0f MeV/c",
                           i, rawHits.size(), clusters.size(), p_fit);

        c3d->cd(k + 1);
        gHT3->SetTitle(hdr + ";x (cm);y (cm);z (cm)");
        gHT3->Draw("P");
        gMC3->Draw("P SAME");
        gCL3->Draw("P SAME");
        gFT3->Draw("LINE SAME");

        cxy->cd(3 * k + 1);
        gHTxy->SetTitle(hdr + ";x (cm);y (cm)");
        gHTxy->Draw("AP");
        gMCxy->Draw("P SAME");
        gCLxy->Draw("P SAME");
        gFTxy->Draw("LP SAME");

        cxy->cd(3 * k + 2);
        gHTxz->SetTitle("x vs z;x (cm);z (cm)");
        gHTxz->Draw("AP");
        gMCxz->Draw("P SAME");
        gCLxz->Draw("P SAME");
        gFTxz->Draw("LP SAME");

        cxy->cd(3 * k + 3);
        gHTyz->SetTitle("z vs y;z (cm);y (cm)");
        gHTyz->Draw("AP");
        gMCyz->Draw("P SAME");
        gCLyz->Draw("P SAME");
        gFTyz->Draw("LP SAME");
    }

    // Legend
    auto* leg = new TLegend(0.05, 0.95, 0.95, 0.99);
    leg->SetNColumns(4);
    leg->SetBorderSize(0);
    auto* lh = new TGraph(1); lh->SetMarkerStyle(20); lh->SetMarkerColor(kAzure + 1); lh->SetMarkerSize(1.0);
    auto* lm = new TGraph(1); lm->SetMarkerStyle(29); lm->SetMarkerColor(kGreen + 2); lm->SetMarkerSize(1.4);
    auto* lc = new TGraph(1); lc->SetMarkerStyle(24); lc->SetMarkerColor(kBlue + 2); lc->SetMarkerSize(1.4);
    auto* lf = new TGraph(1); lf->SetMarkerStyle(34); lf->SetMarkerColor(kRed + 1); lf->SetLineColor(kRed + 1); lf->SetLineWidth(2);
    leg->AddEntry(lh, "Reco hits", "p");
    leg->AddEntry(lm, "MC truth", "p");
    leg->AddEntry(lc, "Clusters", "p");
    leg->AddEntry(lf, "UKF smoothed", "lp");
    cxy->cd(0);
    leg->Draw();

    c3d->SaveAs(Form("tracks_3d_p%d.png", p_MeV));
    cxy->SaveAs(Form("tracks_proj_p%d.png", p_MeV));
    std::cout << "Wrote tracks_3d_p" << p_MeV << ".png and tracks_proj_p" << p_MeV << ".png\n";
}
