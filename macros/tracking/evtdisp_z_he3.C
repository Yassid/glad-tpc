void evtdisp_z_he3()
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TString wd = gSystem->Getenv("VMCWORKDIR");
    TFile fReco(wd + "/glad-tpc/macros/reco/output_reco_goodevt.root");
    TFile fTrk("output_tracking_goodevt.root");
    TFile fSim(wd + "/glad-tpc/macros/sim/Prototype/sim_goodevt.root");
    auto* tR = (TTree*)fReco.Get("evt");
    auto* tT = (TTree*)fTrk.Get("evt");
    auto* tS = (TTree*)fSim.Get("evt");
    auto* recoHits = new TClonesArray("R3BGTPCHitData");
    auto* trks     = new TClonesArray("R3BGTPCTrackData");
    auto* mc       = new TClonesArray("R3BMCTrack");
    auto* pts      = new TClonesArray("R3BGTPCPoint");
    tR->SetBranchAddress("GTPCHitData", &recoHits);
    tT->SetBranchAddress("GTPCTrackData", &trks);
    tS->SetBranchAddress("MCTrack", &mc);
    tS->SetBranchAddress("GTPCPoint", &pts);

    const double offX = 4.2, offZ = 260.2;
    const double padXmin = 0.0, padXmax = 8.8;
    const double padZmin = 0.0, padZmax = 25.6;
    const int PDG_He3 = 1000020030;
    const int PDG_pim = -211;

    // Find events where He3 has at least 3 GTPCPoints in active gas
    struct EvHit { Long64_t evt; int nHe3pts; };
    std::vector<EvHit> picks;
    for (Long64_t i = 0; i < tS->GetEntries(); ++i) {
        tS->GetEntry(i);
        int nHe = 0;
        for (int j = 0; j < pts->GetEntries(); ++j)
            if (((R3BGTPCPoint*)pts->At(j))->GetPDGCode() == PDG_He3) ++nHe;
        if (nHe >= 3) picks.push_back({i, nHe});
    }
    printf("events with >=3 He3 MC points in gas: %lu\n", picks.size());
    for (auto& p : picks) printf("  evt %lld   %d He3 MC pts\n", p.evt, p.nHe3pts);
    if (picks.empty()) return;

    // Match each cluster to {pi-, He3, other} by which PDG dominates a 3D
    // proximity check against MC GTPCPoints.
    auto matchClusterPDG = [&](const std::vector<R3BGTPCHitData>& hits,
                               const std::vector<std::tuple<double,double,double,int>>& mcPoints) -> int
    {
        std::map<int,int> votes;
        for (auto& h : hits) {
            double bestD2 = 1e18; int bestPdg = 0;
            for (auto& [px,py,pz,pdg] : mcPoints) {
                double d2 = (h.GetX()-px)*(h.GetX()-px) +
                            (h.GetY()-py)*(h.GetY()-py) +
                            (h.GetZ()-pz)*(h.GetZ()-pz);
                if (d2 < bestD2) { bestD2 = d2; bestPdg = pdg; }
            }
            if (bestD2 < 4.0) votes[bestPdg]++;  // <2 cm match
        }
        int winner = 0, best = 0;
        for (auto& [pdg, n] : votes) if (n > best) { best = n; winner = pdg; }
        return winner;
    };

    int N = std::min((int)picks.size(), 6);
    int cols = 3, rows = (N + cols - 1) / cols;
    TCanvas* c = new TCanvas("c","", 1800, 400 * rows);
    c->Divide(cols, rows, 0.006, 0.025);

    int ip = 0;
    for (int q = 0; q < N; ++q) {
        Long64_t i = picks[q].evt;
        ++ip; c->cd(ip);
        gPad->SetLeftMargin(0.10);
        gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.10);
        gPad->SetBottomMargin(0.14);
        tR->GetEntry(i); tT->GetEntry(i); tS->GetEntry(i);

        // Collect MC points in local frame, tagged by PDG
        std::vector<std::tuple<double,double,double,int>> mcPts; // (xl,y,zl,pdg)
        for (int j = 0; j < pts->GetEntries(); ++j) {
            auto* p = (R3BGTPCPoint*)pts->At(j);
            mcPts.emplace_back(p->GetX()-offX, p->GetY(), p->GetZ()-offZ, p->GetPDGCode());
        }

        // Identify cluster-to-PDG mapping
        std::vector<int> clusterPdg(trks->GetEntries(), 0);
        for (int k = 0; k < trks->GetEntries(); ++k) {
            auto* tr = (R3BGTPCTrackData*)trks->At(k);
            clusterPdg[k] = matchClusterPDG(tr->GetHitArray(), mcPts);
        }

        // Frame bounds: pad + clusters + vertices
        R3BMCTrack* piMC = nullptr, *heMC = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mc->At(j);
            if (m->GetMotherId() == -1 && m->GetPdgCode() == PDG_pim) piMC = m;
            if (m->GetMotherId() == -1 && m->GetPdgCode() == PDG_He3) heMC = m;
        }
        double zlo = padZmin, zhi = padZmax, xlo = padXmin, xhi = padXmax;
        auto expandXZ = [&](double zl, double xl){
            zlo = std::min(zlo, zl); zhi = std::max(zhi, zl);
            xlo = std::min(xlo, xl); xhi = std::max(xhi, xl);
        };
        if (piMC) expandXZ(piMC->GetStartZ()-offZ, piMC->GetStartX()-offX);
        if (heMC) expandXZ(heMC->GetStartZ()-offZ, heMC->GetStartX()-offX);
        for (int k = 0; k < trks->GetEntries(); ++k) {
            auto* tr = (R3BGTPCTrackData*)trks->At(k);
            for (auto& h : tr->GetHitArray()) expandXZ(h.GetZ(), h.GetX());
        }
        double margin = 1.0;
        auto* frame = gPad->DrawFrame(zlo-margin, xlo-margin, zhi+margin, xhi+margin);

        TString tinfo = TString::Format("evt %lld", (long long)i);
        for (int k = 0; k < trks->GetEntries(); ++k) {
            auto* tr = (R3BGTPCTrackData*)trks->At(k);
            const char* lab = clusterPdg[k] == PDG_pim ? "#pi^{-}"
                            : clusterPdg[k] == PDG_He3 ? "^{3}He"
                            : clusterPdg[k] == 11      ? "e^{-}"
                            : clusterPdg[k] == 0       ? "?"
                            : Form("PDG%d", clusterPdg[k]);
            tinfo += TString::Format("   c%d=%s(%lu)", k, lab, tr->GetHitArray().size());
        }
        frame->SetTitle(tinfo);
        frame->GetXaxis()->SetTitle("z_{local} (beam) [cm]");
        frame->GetYaxis()->SetTitle("x_{local} [cm]");

        // Pad outline
        auto* box = new TBox(padZmin, padXmin, padZmax, padXmax);
        box->SetLineColor(kGray+2); box->SetLineWidth(2); box->SetLineStyle(7); box->SetFillStyle(0);
        box->Draw("L SAME");

        // MC points coloured by particle (pi-/He3/other)
        auto* gMpi = new TGraph(); int mpi = 0;
        auto* gMhe = new TGraph(); int mhe = 0;
        auto* gMot = new TGraph(); int mot = 0;
        for (auto& [xl,y,zl,pdg] : mcPts) {
            if (pdg == PDG_pim) gMpi->SetPoint(mpi++, zl, xl);
            else if (pdg == PDG_He3) gMhe->SetPoint(mhe++, zl, xl);
            else gMot->SetPoint(mot++, zl, xl);
        }
        if (mot > 0) { gMot->SetMarkerStyle(28); gMot->SetMarkerSize(1.0); gMot->SetMarkerColor(kGray+1); gMot->Draw("P SAME"); }
        if (mpi > 0) { gMpi->SetMarkerStyle(29); gMpi->SetMarkerSize(1.6); gMpi->SetMarkerColor(kBlack); gMpi->Draw("P SAME"); }
        if (mhe > 0) { gMhe->SetMarkerStyle(34); gMhe->SetMarkerSize(1.6); gMhe->SetMarkerColor(kAzure+2); gMhe->Draw("P SAME"); }

        // Clusters coloured by matched PDG
        for (int k = 0; k < trks->GetEntries(); ++k) {
            auto* tr = (R3BGTPCTrackData*)trks->At(k);
            auto& hits = tr->GetHitArray();
            int col = clusterPdg[k] == PDG_pim ? kBlue+1
                    : clusterPdg[k] == PDG_He3 ? kRed+1
                    : clusterPdg[k] == 11      ? kOrange+7
                    : kGray+1;
            auto* g = new TGraph(); int gi = 0;
            for (auto& h : hits) g->SetPoint(gi++, h.GetZ(), h.GetX());
            g->SetMarkerStyle(20); g->SetMarkerSize(0.55); g->SetMarkerColor(col);
            g->Draw("P SAME");
        }

        // Vertices (pi-=magenta, He3=teal) + arrows
        auto drawVertex = [&](R3BMCTrack* m, int col, const char* lab) {
            if (!m) return;
            double vzL = m->GetStartZ() - offZ;
            double vxL = m->GetStartX() - offX;
            auto* g = new TGraph(); g->SetPoint(0, vzL, vxL);
            g->SetMarkerStyle(33); g->SetMarkerSize(2.0); g->SetMarkerColor(col);
            g->Draw("P SAME");
            double mag = std::hypot(m->GetPx(), m->GetPz());
            if (mag > 1e-6) {
                double az = vzL + 2.5*m->GetPz()/mag, ax = vxL + 2.5*m->GetPx()/mag;
                auto* arr = new TArrow(vzL, vxL, az, ax, 0.018, "|>");
                arr->SetLineColor(col); arr->SetLineWidth(2); arr->Draw();
            }
        };
        drawVertex(piMC, kMagenta+2, "vtx pi");
        drawVertex(heMC, kAzure+2, "vtx He");
    }
    c->SaveAs("evtdisp_z_he3.png");
}
