/// @file check_data_curvature.C
/// @brief Fit a Kasa circle in (x, z) to (a) MC primary points and (b) reco
///        hits, for the same long-chord events. Compare R_MC vs R_reco vs the
///        truth R = p_MC/(0.3·B). If R_reco > R_MC the bias is in digi/reco
///        (data side); if R_MC > truth_R the bias is in the simulation/MC.
///
/// Both fits in cm; convert to MeV/c via Brho.

#include <array>

static bool kasa_xz(const std::vector<std::array<double, 2>>& pts, double& R, double& cx, double& cz)
{
    const std::size_t n = pts.size();
    if (n < 3) return false;
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, sxr = 0, syr = 0, sr = 0;
    for (auto& p : pts) {
        const double x = p[0], y = p[1];
        const double r2 = x * x + y * y;
        sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
        sxr += x * r2; syr += y * r2; sr += r2;
    }
    const double nd = (double)n;
    const double det = sxx * (syy * nd - sy * sy) - sxy * (sxy * nd - sy * sx) + sx * (sxy * sy - syy * sx);
    if (std::abs(det) < 1e-9) return false;
    const double dA = sxr * (syy * nd - sy * sy) - sxy * (syr * nd - sy * sr) + sx * (syr * sy - syy * sr);
    const double dB = sxx * (syr * nd - sy * sr) - sxr * (sxy * nd - sy * sx) + sx * (sxy * sr - syr * sx);
    const double dC = sxx * (syy * sr - syr * sy) - sxy * (sxy * sr - syr * sx) + sxr * (sxy * sy - syy * sx);
    const double A = dA / det, B = dB / det, C = dC / det;
    cx = A / 2; cz = B / 2;
    const double r2 = C + cx * cx + cz * cz;
    if (r2 <= 0) return false;
    R = std::sqrt(r2);
    return true;
}

void check_data_curvature()
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);

    const double offX = 4.2, offZ = 260.2;
    const double B_T = 2.0;
    const double mass_pi = 139.57;

    TString workDir = gSystem->Getenv("VMCWORKDIR");
    TFile fSim(workDir + "/glad-tpc/macros/sim/Prototype/sim.root");
    TFile fTrk(workDir + "/glad-tpc/macros/tracking/output_tracking.root");
    auto* tSim = (TTree*)fSim.Get("evt");
    auto* tTrk = (TTree*)fTrk.Get("evt");
    auto* pts = new TClonesArray("R3BGTPCPoint");
    auto* mc = new TClonesArray("R3BMCTrack");
    auto* trks = new TClonesArray("R3BGTPCTrackData");
    tSim->SetBranchAddress("GTPCPoint", &pts);
    tSim->SetBranchAddress("MCTrack", &mc);
    tTrk->SetBranchAddress("GTPCTrackData", &trks);

    auto* gMCvsTruth = new TGraph();
    auto* gRecoVsMC = new TGraph();
    auto* gRecoVsTruth = new TGraph();
    auto* hMCvsTruth = new TH1F("hMCvsT", ";R_{MC fit}/R_{truth}-1;events", 60, -0.5, 0.5);
    auto* hRecoVsMC = new TH1F("hRecoMC", ";R_{reco hits}/R_{MC}-1;events", 60, -0.5, 0.5);
    auto* hRecoVsTruth = new TH1F("hRecoT", ";R_{reco hits}/R_{truth}-1;events", 60, -0.5, 0.5);

    int n = 0;
    for (Long64_t i = 0; i < std::min(tSim->GetEntries(), tTrk->GetEntries()); ++i) {
        tSim->GetEntry(i); tTrk->GetEntry(i);
        if (trks->GetEntries() == 0) continue;
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        auto& hits = tr->GetHitArray();
        if (hits.size() < 30) continue;

        // chord cut
        double zmn = 1e9, zmx = -1e9, xmn = 1e9, xmx = -1e9;
        for (auto& h : hits) { zmn = std::min(zmn, (double)h.GetZ()); zmx = std::max(zmx, (double)h.GetZ());
                               xmn = std::min(xmn, (double)h.GetX()); xmx = std::max(xmx, (double)h.GetX()); }
        if (std::sqrt((zmx - zmn) * (zmx - zmn) + (xmx - xmn) * (xmx - xmn)) < 16.0) continue;

        // MC primary points in pad-local (x, z)
        std::vector<std::array<double, 2>> mcPts;
        for (int j = 0; j < pts->GetEntries(); ++j) {
            auto* p = (R3BGTPCPoint*)pts->At(j);
            if (std::abs(p->GetPDGCode()) != 211) continue;
            mcPts.push_back({ p->GetX() - offX, p->GetZ() - offZ });
        }
        if (mcPts.size() < 3) continue;

        // Reco hit positions (already pad-local, cm)
        std::vector<std::array<double, 2>> recoPts;
        for (auto& h : hits) recoPts.push_back({ (double)h.GetX(), (double)h.GetZ() });

        // Truth R from MC primary momentum
        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mc->At(j);
            if (m->GetMotherId() == -1 && std::abs(m->GetPdgCode()) == 211) { pi = m; break; }
        }
        if (!pi) continue;
        const double pT_truth_GeV = std::sqrt(pi->GetPx() * pi->GetPx() + pi->GetPz() * pi->GetPz());
        const double R_truth_cm = pT_truth_GeV / (0.3 * B_T) * 100.0;
        const double p_truth_MeV = std::sqrt(pi->GetPx() * pi->GetPx() + pi->GetPy() * pi->GetPy() + pi->GetPz() * pi->GetPz()) * 1000;

        // Kasa on MC primary points
        double R_mc, cx_mc, cz_mc;
        if (!kasa_xz(mcPts, R_mc, cx_mc, cz_mc)) continue;
        // Kasa on reco hits
        double R_reco, cx_re, cz_re;
        if (!kasa_xz(recoPts, R_reco, cx_re, cz_re)) continue;

        gMCvsTruth->SetPoint(n, p_truth_MeV, R_mc / R_truth_cm - 1);
        gRecoVsMC->SetPoint(n, p_truth_MeV, R_reco / R_mc - 1);
        gRecoVsTruth->SetPoint(n, p_truth_MeV, R_reco / R_truth_cm - 1);
        hMCvsTruth->Fill(R_mc / R_truth_cm - 1);
        hRecoVsMC->Fill(R_reco / R_mc - 1);
        hRecoVsTruth->Fill(R_reco / R_truth_cm - 1);
        ++n;
    }
    std::cout << "events: " << n << "\n";
    auto report = [](TH1F* h, const char* l) {
        std::printf("  %-32s mean=%+6.3f RMS=%5.3f  (median %+6.3f)\n",
                    l, h->GetMean(), h->GetRMS(),
                    h->GetBinCenter(h->GetMaximumBin()));
    };
    report(hMCvsTruth,   "R_MC_kasa / R_truth - 1");
    report(hRecoVsMC,    "R_reco / R_MC_kasa - 1");
    report(hRecoVsTruth, "R_reco / R_truth - 1");

    auto* c = new TCanvas("c", "data curvature", 1500, 900);
    c->Divide(3, 2);
    c->cd(1); hMCvsTruth->Draw();
    c->cd(2); hRecoVsMC->Draw();
    c->cd(3); hRecoVsTruth->Draw();
    auto* line = new TLine(); line->SetLineStyle(2); line->SetLineColor(kRed);
    for (int k = 1; k <= 3; ++k) { c->cd(k); line->DrawLine(0, 0, 0, 30); }
    c->cd(4); gMCvsTruth->SetTitle("R_{MC kasa}/R_{truth}-1 vs p_{truth};p (MeV/c);ratio-1"); gMCvsTruth->SetMarkerStyle(7); gMCvsTruth->Draw("AP");
    c->cd(5); gRecoVsMC->SetTitle("R_{reco}/R_{MC kasa}-1 vs p_{truth};p (MeV/c);ratio-1"); gRecoVsMC->SetMarkerStyle(7); gRecoVsMC->Draw("AP");
    c->cd(6); gRecoVsTruth->SetTitle("R_{reco}/R_{truth}-1 vs p_{truth};p (MeV/c);ratio-1"); gRecoVsTruth->SetMarkerStyle(7); gRecoVsTruth->Draw("AP");
    c->SaveAs("data_curvature.png");
}
