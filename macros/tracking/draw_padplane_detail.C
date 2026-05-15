/// @file draw_padplane_detail.C
/// @brief Single-event pad-plane view with MC truth, reco hits, TripletClust
///        clusters, and UKF smoothed track all overlaid in (z, x) pad-local.
///
/// World → pad transform: z_pad = z_world − 260.2 cm, x_pad = x_world − 4.2 cm
/// (HYDRAprototype_FileSetup_v2_02082022.par: GladOffsetX/Z).

void draw_padplane_detail(int p_MeV = 600, int evt = -1)
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

    // Auto-pick a representative event if none requested: want a converged fit
    // with 30 < n_pads < 200 and ≥3 MC primary points so the helix sample is
    // meaningful.
    Long64_t pickEvt = evt;
    if (pickEvt < 0) {
        const Long64_t N = std::min({ tSim->GetEntries(), tLang->GetEntries(),
                                       tTrk->GetEntries(), tUKF->GetEntries() });
        for (Long64_t i = 0; i < N; ++i) {
            tLang->GetEntry(i);
            tTrk->GetEntry(i);
            tUKF->GetEntry(i);
            tSim->GetEntry(i);
            if (cal->GetEntries() < 30 || cal->GetEntries() > 200) continue;
            if (fits->GetEntries() == 0 || trks->GetEntries() == 0) continue;
            int nPri = 0;
            for (int j = 0; j < gtpcPts->GetEntries(); ++j)
                if (((R3BGTPCPoint*)gtpcPts->At(j))->GetTrackID() == 0) ++nPri;
            if (nPri < 3) continue;
            pickEvt = i;
            break;
        }
    }
    if (pickEvt < 0) { std::cout << "[detail] no suitable event\n"; return; }

    tSim->GetEntry(pickEvt);
    tLang->GetEntry(pickEvt);
    tTrk->GetEntry(pickEvt);
    tUKF->GetEntry(pickEvt);

    // MC primary pion
    R3BMCTrack* pi = nullptr;
    for (int j = 0; j < mcTracks->GetEntries(); ++j) {
        auto* m = (R3BMCTrack*)mcTracks->At(j);
        if (m->GetMotherId() == -1 && std::abs(m->GetPdgCode()) == 211) { pi = m; break; }
    }
    const double pmc_MeV = pi ? std::sqrt(pi->GetPx() * pi->GetPx() + pi->GetPy() * pi->GetPy() + pi->GetPz() * pi->GetPz()) * 1000 : 0;
    auto* fit = fits->GetEntries() > 0 ? (R3BGTPCFittedTrackData*)fits->At(0) : nullptr;
    const double ke_fit = fit ? fit->GetKinematics().kineticEnergy : 0;
    const double p_fit = std::sqrt(ke_fit * ke_fit + 2 * ke_fit * 139.57);

    // ---- Pad-plane 2D ADC ----
    const int nColZ = 128, nRowX = 44;
    const double padSize_mm = 2.0;
    const double padPlaneZ_cm = nColZ * padSize_mm / 10.0;
    const double padPlaneX_cm = nRowX * padSize_mm / 10.0;
    auto* hPad = new TH2F("hPad",
                          Form("p_MC=%.0f MeV/c, p_fit=%.0f MeV/c, N_pads=%d;z (cm, pad-local);x (cm, pad-local)",
                               pmc_MeV, p_fit, cal->GetEntries()),
                          nColZ, 0., padPlaneZ_cm,
                          nRowX, 0., padPlaneX_cm);
    for (int j = 0; j < cal->GetEntries(); ++j) {
        auto* d = (R3BGTPCCalData*)cal->At(j);
        const int pid = d->GetPadId();
        const int icol = pid / nRowX, irow = pid % nRowX;
        const double z_cm = (icol + 0.5) * padSize_mm / 10.0;
        const double x_cm = (irow + 0.5) * padSize_mm / 10.0;
        double q = 0.;
        for (auto a : d->GetADC()) q += a;
        hPad->SetBinContent(icol + 1, irow + 1, q);
    }

    // ---- MC truth primary (trackID==0) ----
    auto* gMC = new TGraph();
    gMC->SetMarkerStyle(29);
    gMC->SetMarkerColor(kGreen + 2);
    gMC->SetMarkerSize(2.5);
    for (int j = 0; j < gtpcPts->GetEntries(); ++j) {
        auto* p = (R3BGTPCPoint*)gtpcPts->At(j);
        if (p->GetTrackID() != 0) continue;
        gMC->SetPoint(gMC->GetN(), p->GetZ() - offZ_cm, p->GetX() - offX_cm);
    }
    // ---- MC secondaries ----
    auto* gSec = new TGraph();
    gSec->SetMarkerStyle(25);
    gSec->SetMarkerColor(kOrange + 1);
    gSec->SetMarkerSize(1.0);
    for (int j = 0; j < gtpcPts->GetEntries(); ++j) {
        auto* p = (R3BGTPCPoint*)gtpcPts->At(j);
        if (p->GetTrackID() == 0) continue;
        gSec->SetPoint(gSec->GetN(), p->GetZ() - offZ_cm, p->GetX() - offX_cm);
    }

    // ---- Reconstructed hits ----
    auto* tr = trks->GetEntries() > 0 ? (R3BGTPCTrackData*)trks->At(0) : nullptr;
    auto* gHit = new TGraph();
    gHit->SetMarkerStyle(20);
    gHit->SetMarkerColor(kBlack);
    gHit->SetMarkerSize(0.5);
    if (tr) {
        for (auto& h : tr->GetHitArray())
            gHit->SetPoint(gHit->GetN(), h.GetZ(), h.GetX());
    }
    // ---- TripletClust clusters ----
    auto* gCl = new TGraph();
    gCl->SetMarkerStyle(24);
    gCl->SetMarkerColor(kRed);
    gCl->SetMarkerSize(1.6);
    if (tr) {
        for (auto& c : *tr->GetHitClusterArray())
            gCl->SetPoint(gCl->GetN(), c.GetZ(), c.GetX());
    }
    // ---- UKF smoothed track ----
    auto* gFit = new TGraph();
    gFit->SetMarkerStyle(34);
    gFit->SetMarkerColor(kMagenta + 1);
    gFit->SetMarkerSize(1.0);
    gFit->SetLineColor(kMagenta + 1);
    gFit->SetLineWidth(2);
    if (fit) {
        for (auto& s : fit->GetSmoothedPositions())
            gFit->SetPoint(gFit->GetN(), s.Z(), s.X());
    }

    // Zoom into the actually-occupied region (track footprint + a 1 cm pad)
    double zlo = 1e9, zhi = -1e9, xlo = 1e9, xhi = -1e9;
    for (int j = 0; j < cal->GetEntries(); ++j) {
        auto* d = (R3BGTPCCalData*)cal->At(j);
        const int pid = d->GetPadId();
        const double z_cm = (pid / nRowX + 0.5) * padSize_mm / 10.0;
        const double x_cm = (pid % nRowX + 0.5) * padSize_mm / 10.0;
        zlo = std::min(zlo, z_cm); zhi = std::max(zhi, z_cm);
        xlo = std::min(xlo, x_cm); xhi = std::max(xhi, x_cm);
    }
    hPad->GetXaxis()->SetRangeUser(std::max(0., zlo - 1.0), std::min(padPlaneZ_cm, zhi + 1.0));
    hPad->GetYaxis()->SetRangeUser(std::max(0., xlo - 1.0), std::min(padPlaneX_cm, xhi + 1.0));

    auto* c = new TCanvas("cdet", "pad plane detail", 1200, 800);
    c->SetRightMargin(0.13);
    c->SetLeftMargin(0.10);
    hPad->Draw("COLZ");
    gMC->Draw("P SAME");
    gSec->Draw("P SAME");
    gHit->Draw("P SAME");
    gCl->Draw("P SAME");
    gFit->Draw("LP SAME");

    auto* leg = new TLegend(0.62, 0.13, 0.86, 0.42);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->AddEntry(hPad, "Digitized pad ADC", "f");
    leg->AddEntry(gMC, "MC primary (#pi^{-})", "p");
    leg->AddEntry(gSec, "MC secondaries", "p");
    leg->AddEntry(gHit, "Reco hits (Cal2Hit)", "p");
    leg->AddEntry(gCl, "TripletClust clusters", "p");
    leg->AddEntry(gFit, "UKF smoothed", "lp");
    leg->Draw();

    c->SaveAs(Form("padplane_detail_p%d_evt%lld.png", p_MeV, pickEvt));
    std::cout << "Wrote padplane_detail_p" << p_MeV << "_evt" << pickEvt << ".png\n";
}
