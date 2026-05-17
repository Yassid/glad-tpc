// Per-momentum-point σ_p/p tally for the box-generator scan.
// Reads output_tracking_pX.root + output_ukf_pX.root + sim_pX.root and
// returns one CSV row:
//   p_MeV, N_fit, σ_R/R_seed_pct, σ_p/p_seed_pct, σ_p/p_UKF_pct,
//   bias_R_seed_pct, bias_p_seed_pct, bias_p_UKF_pct
// (where σ is the Gaussian core from a fit of the central peak in
//  p/p_MC - 1, and bias is the Gaussian mean of the same fit).
//
// Usage:
//   root -b -q 'analyze_scan_point.C(800)'      → prints CSV line, returns
void analyze_scan_point(int p_mev = 800, const char* tag = "")
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    TString wd = gSystem->Getenv("VMCWORKDIR");
    TString suf = TString(tag) + Form("_p%d", p_mev);
    TFile fS(wd + "/glad-tpc/macros/sim/Prototype/sim" + suf + ".root");
    TFile fT("output_tracking" + suf + ".root");
    TFile fU("output_ukf" + suf + ".root");
    auto* tS = (TTree*)fS.Get("evt");
    auto* tT = (TTree*)fT.Get("evt");
    auto* tU = (TTree*)fU.Get("evt");
    if (!tS || !tT || !tU) {
        printf("RESULT %d 0 NA NA NA NA NA NA\n", p_mev);
        return;
    }
    auto* mc     = new TClonesArray("R3BMCTrack");
    auto* trks   = new TClonesArray("R3BGTPCTrackData");
    auto* fitted = new TClonesArray("R3BGTPCFittedTrackData");
    tS->SetBranchAddress("MCTrack", &mc);
    tT->SetBranchAddress("GTPCTrackData", &trks);
    tU->SetBranchAddress("GTPCFittedTrackData", &fitted);
    const double m_pi = 139.57039;
    const double B = 2.0;

    auto* hR = new TH1F("hR",";R_{fit}/R_{tru}-1;", 80, -1, 1);
    auto* hPSeed = new TH1F("hPS",";p_{seed}/p_{MC}-1;", 80, -1, 1);
    auto* hPUKF  = new TH1F("hPU",";p_{UKF}/p_{MC}-1;",  80, -1, 1);
    int nFit = 0;

    for (Long64_t i = 0; i < std::min({tS->GetEntries(), tT->GetEntries(), tU->GetEntries()}); ++i) {
        tS->GetEntry(i); tT->GetEntry(i); tU->GetEntry(i);
        if (trks->GetEntries() == 0) continue;
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        if (tr->GetHitArray().size() < 10) continue;
        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mc->At(j);
            if (m->GetMotherId()==-1 && m->GetPdgCode()==-211) { pi = m; break; }
        }
        if (!pi) continue;
        double pT_MC = std::sqrt(pi->GetPx()*pi->GetPx() + pi->GetPz()*pi->GetPz()) * 1000;
        double p_MC  = std::sqrt(pi->GetPx()*pi->GetPx() + pi->GetPy()*pi->GetPy() + pi->GetPz()*pi->GetPz()) * 1000;
        double R_truth = pT_MC / (3.0 * B);
        double R_fit = tr->GetGeoRadius();
        double theta = tr->GetGeoTheta();
        if (!std::isfinite(R_fit) || R_fit <= 0) continue;
        if (!std::isfinite(theta)) continue;
        double sinTh = std::sin(theta);
        if (std::abs(sinTh) < 0.1) continue;
        double pT_seed = 0.3 * B * (R_fit / 100.0) * 1000.0;
        double p_seed  = pT_seed / sinTh;
        hR->Fill(R_fit/R_truth - 1);
        hPSeed->Fill(p_seed/p_MC - 1);
        // UKF
        if (fitted->GetEntries() > 0) {
            auto* ft = (R3BGTPCFittedTrackData*)fitted->At(0);
            if (ft->IsConverged()) {
                const auto& kin = ft->GetKinematicsXtr();
                if (std::isfinite(kin.kineticEnergy) && kin.kineticEnergy > 0) {
                    double p_UKF = std::sqrt((kin.kineticEnergy + m_pi)*(kin.kineticEnergy + m_pi) - m_pi*m_pi);
                    hPUKF->Fill(p_UKF/p_MC - 1);
                }
            }
        }
        ++nFit;
    }

    auto gauss = [](TH1F* h) -> std::pair<double,double> {
        if (h->GetEntries() < 5) return {NAN, NAN};
        h->Fit("gaus", "Q", "", -0.5, 0.5);
        auto* f = h->GetFunction("gaus");
        if (!f) return {NAN, NAN};
        return { f->GetParameter(1), f->GetParameter(2) };
    };
    auto [biasR, sigR] = gauss(hR);
    auto [biasPS, sigPS] = gauss(hPSeed);
    auto [biasPU, sigPU] = gauss(hPUKF);

    printf("RESULT %d %d %.3f %.3f %.3f %.3f %.3f %.3f\n",
           p_mev, nFit,
           100*sigR, 100*sigPS, 100*sigPU,
           100*biasR, 100*biasPS, 100*biasPU);
}
