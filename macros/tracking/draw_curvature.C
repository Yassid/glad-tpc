/// @file draw_curvature.C
/// @brief Three views of the actual track curvature for a single event:
///   (a) pad-plane (z, x) with cluster centroids + analytic-truth helix
///   (b) cluster (z, x) residual to a least-squares line — amplifies curvature
///   (c) numeric summary: expected sagitta, observed RMS

void draw_curvature(int p_MeV = 600, int evt = -1)
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gSystem->Load("libR3BOpenKF");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    const double offX_cm = 4.2;
    const double offZ_cm = 260.2;
    const double B_T = 2.0;

    TString workDir = gSystem->Getenv("VMCWORKDIR");
    TFile fSim(TString::Format("%s/glad-tpc/macros/sim/Prototype/sim_p%d.root", workDir.Data(), p_MeV));
    TFile fLang(TString::Format("%s/glad-tpc/macros/proj/Prototype/lang_p%d.root", workDir.Data(), p_MeV));
    TFile fTrk(TString::Format("%s/glad-tpc/macros/tracking/output_tracking_p%d.root", workDir.Data(), p_MeV));
    auto* tSim = (TTree*)fSim.Get("evt");
    auto* tLang = (TTree*)fLang.Get("evt");
    auto* tTrk = (TTree*)fTrk.Get("evt");
    auto* gtpcPts = new TClonesArray("R3BGTPCPoint");
    auto* mcTracks = new TClonesArray("R3BMCTrack");
    auto* cal = new TClonesArray("R3BGTPCCalData");
    auto* trks = new TClonesArray("R3BGTPCTrackData");
    tSim->SetBranchAddress("GTPCPoint", &gtpcPts);
    tSim->SetBranchAddress("MCTrack", &mcTracks);
    tLang->SetBranchAddress("GTPCCalData", &cal);
    tTrk->SetBranchAddress("GTPCTrackData", &trks);

    Long64_t pick = evt;
    if (pick < 0) {
        for (Long64_t i = 0; i < std::min({ tSim->GetEntries(), tLang->GetEntries(), tTrk->GetEntries() }); ++i) {
            tLang->GetEntry(i); tTrk->GetEntry(i);
            if (cal->GetEntries() < 50 || cal->GetEntries() > 200) continue;
            if (trks->GetEntries() == 0) continue;
            if (((R3BGTPCTrackData*)trks->At(0))->GetHitArray().size() < 50) continue;
            pick = i; break;
        }
    }
    tSim->GetEntry(pick); tLang->GetEntry(pick); tTrk->GetEntry(pick);

    R3BMCTrack* pi = nullptr;
    for (int j = 0; j < mcTracks->GetEntries(); ++j) {
        auto* m = (R3BMCTrack*)mcTracks->At(j);
        if (m->GetMotherId() == -1 && std::abs(m->GetPdgCode()) == 211) { pi = m; break; }
    }
    const double px = pi->GetPx(), py = pi->GetPy(), pz = pi->GetPz();
    const double pt = std::sqrt(px * px + pz * pz);
    const double R_cm = pt / (0.3 * B_T) * 100.0;
    const double pTot = std::sqrt(px * px + py * py + pz * pz);

    // First MC primary point as anchor for the analytic helix
    double vx0 = 0, vy0 = 0, vz0 = 0;
    for (int j = 0; j < gtpcPts->GetEntries(); ++j) {
        auto* p = (R3BGTPCPoint*)gtpcPts->At(j);
        if (p->GetTrackID() == 0) { vx0 = p->GetX(); vy0 = p->GetY(); vz0 = p->GetZ(); break; }
    }
    // Centre of curvature in world (x,z): for q<0 in +B_y, F_xz ∝ (vz, -vx)
    const double nhx = pz / pt;
    const double nhz = -px / pt;
    const double cx = vx0 + R_cm * nhx;
    const double cz = vz0 + R_cm * nhz;
    const double phi0 = std::atan2(vz0 - cz, vx0 - cx);
    const double dsxz_per_ds = pt / pTot;
    const double dphi_per_ds = -dsxz_per_ds / R_cm; // CW from +y view (q<0, B>0)

    auto* tr = (R3BGTPCTrackData*)trks->At(0);
    auto& hits = tr->GetHitArray();
    // Find z range of clusters for helix extent
    double zMin = 1e9, zMax = -1e9;
    for (auto& h : hits) { zMin = std::min(zMin, (double)h.GetZ()); zMax = std::max(zMax, (double)h.GetZ()); }

    // ---- Pad plane with helix and clusters ----
    const int nColZ = 128, nRowX = 44;
    const double padSize = 0.2; // cm
    const double Zmax = nColZ * padSize;
    const double Xmax = nRowX * padSize;
    auto* hPad = new TH2F("hPad",
                          Form("p=%d MeV/c (evt %lld), R = %.0f cm, expected sagitta over %.1f cm path = %.2f mm",
                               p_MeV, pick, R_cm, zMax - zMin,
                               (zMax - zMin) * (zMax - zMin) / (8 * R_cm) * 10.0),
                          nColZ, 0., Zmax, nRowX, 0., Xmax);
    double zlo = 1e9, zhi = -1e9, xlo = 1e9, xhi = -1e9;
    for (int j = 0; j < cal->GetEntries(); ++j) {
        auto* d = (R3BGTPCCalData*)cal->At(j);
        const int pid = d->GetPadId();
        const int icol = pid / nRowX, irow = pid % nRowX;
        const double z_cm = (icol + 0.5) * padSize;
        const double x_cm = (irow + 0.5) * padSize;
        double q = 0; for (auto a : d->GetADC()) q += a;
        hPad->SetBinContent(icol + 1, irow + 1, q);
        zlo = std::min(zlo, z_cm); zhi = std::max(zhi, z_cm);
        xlo = std::min(xlo, x_cm); xhi = std::max(xhi, x_cm);
    }
    hPad->GetXaxis()->SetRangeUser(std::max(0., zlo - 0.6), std::min(Zmax, zhi + 0.6));
    hPad->GetYaxis()->SetRangeUser(std::max(0., xlo - 0.6), std::min(Xmax, xhi + 0.6));

    // Analytic helix (the *correct* trajectory for the MC momentum)
    auto* gHelix = new TGraph();
    gHelix->SetLineColor(kRed); gHelix->SetLineWidth(2); gHelix->SetLineStyle(2);
    for (int k = 0; k <= 4000; ++k) {
        const double s = k * 0.02;
        const double phi = phi0 + dphi_per_ds * s;
        const double xw = cx + R_cm * std::cos(phi);
        const double zw = cz + R_cm * std::sin(phi);
        const double z_pad = zw - offZ_cm;
        const double x_pad = xw - offX_cm;
        if (z_pad > zMax + 0.5) break;
        if (z_pad < 0 || x_pad < 0 || z_pad > Zmax || x_pad > Xmax) continue;
        gHelix->SetPoint(gHelix->GetN(), z_pad, x_pad);
    }
    // Cluster centroids
    auto* gCl = new TGraph();
    gCl->SetMarkerStyle(24); gCl->SetMarkerColor(kRed); gCl->SetMarkerSize(1.6);
    for (auto& c : *tr->GetHitClusterArray()) gCl->SetPoint(gCl->GetN(), c.GetZ(), c.GetX());

    // ---- Residual to a LS line on (z, x) using hits ----
    int nHit = hits.size();
    double sz = 0, sx = 0, szz = 0, szx = 0;
    for (auto& h : hits) { sz += h.GetZ(); sx += h.GetX(); szz += h.GetZ() * h.GetZ(); szx += h.GetZ() * h.GetX(); }
    const double det = nHit * szz - sz * sz;
    const double slope = (nHit * szx - sz * sx) / det;
    const double intercept = (sx - slope * sz) / nHit;
    auto* gRes = new TGraph();
    gRes->SetMarkerStyle(20); gRes->SetMarkerColor(kBlack); gRes->SetMarkerSize(0.6);
    double resMin = 1e9, resMax = -1e9;
    for (auto& h : hits) {
        const double pred = intercept + slope * h.GetZ();
        const double res_mm = (h.GetX() - pred) * 10.0;
        gRes->SetPoint(gRes->GetN(), h.GetZ(), res_mm);
        resMin = std::min(resMin, res_mm); resMax = std::max(resMax, res_mm);
    }
    // Analytic-helix residual to the same LS line (sampled at same z grid)
    auto* gHelixRes = new TGraph();
    gHelixRes->SetLineColor(kRed); gHelixRes->SetLineWidth(2); gHelixRes->SetLineStyle(2);
    for (int i = 0; i < gHelix->GetN(); ++i) {
        double z, x;
        gHelix->GetPoint(i, z, x);
        const double pred = intercept + slope * z;
        gHelixRes->SetPoint(gHelixRes->GetN(), z, (x - pred) * 10.0);
    }

    auto* c = new TCanvas("cc", "curvature", 1500, 1100);
    c->Divide(1, 2);
    c->cd(1);
    gPad->SetRightMargin(0.13);
    hPad->Draw("COLZ");
    gHelix->Draw("L SAME");
    gCl->Draw("P SAME");
    auto* leg = new TLegend(0.62, 0.13, 0.85, 0.32);
    leg->SetFillStyle(1001); leg->SetFillColor(kWhite); leg->SetTextSize(0.030);
    leg->AddEntry(hPad, "Digitized pad ADC", "f");
    leg->AddEntry(gCl, "TripletClust clusters", "p");
    leg->AddEntry(gHelix, Form("Analytic helix (p=%d MeV/c)", p_MeV), "l");
    leg->Draw();

    c->cd(2);
    auto* hAxis = new TH2F("axis",
                          Form("Residual to LS straight line   (vertical exaggeration shows curvature);z (cm, pad);residual x [mm]"),
                          100, std::max(0., zlo - 0.6), std::min(Zmax, zhi + 0.6),
                          100, std::min(resMin, -1.) - 0.5, std::max(resMax, 1.) + 0.5);
    hAxis->Draw();
    auto* line0 = new TLine(std::max(0., zlo - 0.6), 0, std::min(Zmax, zhi + 0.6), 0);
    line0->SetLineColor(kGray + 2); line0->SetLineStyle(3); line0->Draw();
    gRes->Draw("P SAME");
    gHelixRes->Draw("L SAME");
    auto* leg2 = new TLegend(0.13, 0.78, 0.45, 0.92);
    leg2->SetFillStyle(1001); leg2->SetFillColor(kWhite); leg2->SetTextSize(0.030);
    leg2->AddEntry(gRes, "Hits − LS straight line", "p");
    leg2->AddEntry(gHelixRes, "Analytic helix residual to same line", "l");
    leg2->Draw();

    c->SaveAs(Form("curvature_p%d_evt%lld.png", p_MeV, pick));
    std::cout << "Wrote curvature_p" << p_MeV << "_evt" << pick << ".png  "
              << "(R=" << R_cm << " cm, sagitta over "
              << (zMax - zMin) << " cm path = "
              << (zMax - zMin) * (zMax - zMin) / (8 * R_cm) * 10. << " mm)\n";
}
