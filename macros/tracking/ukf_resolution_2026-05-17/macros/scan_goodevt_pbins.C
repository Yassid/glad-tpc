// σ_p/p vs pion momentum on the good_evt sim, binned by MC truth momentum.
//
// Reads sim_goodevt2k + output_tracking_goodevt2k + output_ukf_goodevt2k,
// puts each event into a 100 MeV/c bin in p_MC, and Gauss-fits the central
// peak in (R_fit/R_truth - 1), (p_seed/p_MC - 1), (p_UKF/p_MC - 1) per bin.
// Writes the same CSV columns as scan_p_results.csv so the two sweeps can
// be plotted on the same axes.
void scan_goodevt_pbins(TString tag = "goodevt2k",
                        double pMin = 100, double pMax = 1100, double pStep = 100)
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    TString wd = gSystem->Getenv("VMCWORKDIR");
    TFile fS(wd + "/glad-tpc/macros/sim/Prototype/sim_" + tag + ".root");
    TFile fT(wd + "/glad-tpc/macros/tracking/output_tracking_" + tag + ".root");
    TFile fU(wd + "/glad-tpc/macros/tracking/output_ukf_"      + tag + ".root");
    auto* tS = (TTree*)fS.Get("evt");
    auto* tT = (TTree*)fT.Get("evt");
    auto* tU = (TTree*)fU.Get("evt");
    if (!tS || !tT || !tU) { std::cerr << "missing input file\n"; return; }
    auto* mc     = new TClonesArray("R3BMCTrack");
    auto* trks   = new TClonesArray("R3BGTPCTrackData");
    auto* fitted = new TClonesArray("R3BGTPCFittedTrackData");
    tS->SetBranchAddress("MCTrack", &mc);
    tT->SetBranchAddress("GTPCTrackData", &trks);
    tU->SetBranchAddress("GTPCFittedTrackData", &fitted);
    const double m_pi = 139.57039;
    const double B = 2.0;

    const int nBins = (int)std::round((pMax - pMin) / pStep);
    struct Bin { TH1F* hR; TH1F* hPS; TH1F* hPU; int n = 0; };
    std::vector<Bin> bins(nBins);
    for (int i = 0; i < nBins; ++i) {
        double lo = pMin + i*pStep, hi = lo + pStep;
        bins[i].hR  = new TH1F(Form("hR_%d",  i), Form("p in [%g,%g];R_{fit}/R_{tru}-1;", lo, hi), 80, -1, 1);
        bins[i].hPS = new TH1F(Form("hPS_%d", i), Form("p in [%g,%g];p_{seed}/p_{MC}-1;", lo, hi), 80, -1, 1);
        bins[i].hPU = new TH1F(Form("hPU_%d", i), Form("p in [%g,%g];p_{UKF}/p_{MC}-1;",  lo, hi), 80, -1, 1);
    }

    for (Long64_t i = 0; i < std::min({tS->GetEntries(), tT->GetEntries(), tU->GetEntries()}); ++i) {
        tS->GetEntry(i); tT->GetEntry(i); tU->GetEntry(i);
        if (trks->GetEntries() == 0) continue;
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        if (tr->GetHitArray().size() < 10) continue;
        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* mt = (R3BMCTrack*)mc->At(j);
            if (mt->GetMotherId()==-1 && mt->GetPdgCode()==-211) { pi = mt; break; }
        }
        if (!pi) continue;
        double pT_MC = std::sqrt(pi->GetPx()*pi->GetPx() + pi->GetPz()*pi->GetPz()) * 1000;
        double p_MC  = std::sqrt(pi->GetPx()*pi->GetPx() + pi->GetPy()*pi->GetPy() + pi->GetPz()*pi->GetPz()) * 1000;
        if (p_MC < pMin || p_MC >= pMax) continue;
        int b = std::min(nBins-1, std::max(0, (int)((p_MC - pMin) / pStep)));

        double R_truth = pT_MC / (3.0 * B);
        double R_fit = tr->GetGeoRadius();
        double theta = tr->GetGeoTheta();
        if (!std::isfinite(R_fit) || R_fit <= 0 || !std::isfinite(theta)) continue;
        double sinTh = std::sin(theta);
        if (std::abs(sinTh) < 0.1) continue;
        double pT_seed = 0.3 * B * (R_fit / 100.0) * 1000.0;
        double p_seed  = pT_seed / sinTh;

        bins[b].hR ->Fill(R_fit/R_truth - 1);
        bins[b].hPS->Fill(p_seed/p_MC   - 1);
        if (fitted->GetEntries() > 0) {
            auto* ft = (R3BGTPCFittedTrackData*)fitted->At(0);
            if (ft->IsConverged()) {
                const auto& kin = ft->GetKinematicsXtr();
                if (std::isfinite(kin.kineticEnergy) && kin.kineticEnergy > 0) {
                    double p_UKF = std::sqrt((kin.kineticEnergy + m_pi)*(kin.kineticEnergy + m_pi) - m_pi*m_pi);
                    bins[b].hPU->Fill(p_UKF/p_MC - 1);
                }
            }
        }
        bins[b].n++;
    }

    TString csv = wd + "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17/scan_goodevt_pbins.csv";
    std::ofstream of(csv.Data());
    of << "p_MeV,N_fit,sigma_R_pct,sigma_p_seed_pct,sigma_p_UKF_pct,bias_R_pct,bias_p_seed_pct,bias_p_UKF_pct\n";
    printf("p_MeV  N_fit  sigma_R    sigma_p_seed  sigma_p_UKF  bias_R   bias_p_seed  bias_p_UKF\n");
    auto gauss = [](TH1F* h) -> std::pair<double,double> {
        if (h->GetEntries() < 5) return {std::nan(""), std::nan("")};
        h->Fit("gaus", "Q", "", -0.4, 0.4);
        auto* f = h->GetFunction("gaus");
        if (!f) return {std::nan(""), std::nan("")};
        return { f->GetParameter(1), f->GetParameter(2) };
    };
    for (int i = 0; i < nBins; ++i) {
        double pMid = pMin + (i + 0.5)*pStep;
        auto [bR, sR]   = gauss(bins[i].hR);
        auto [bPS, sPS] = gauss(bins[i].hPS);
        auto [bPU, sPU] = gauss(bins[i].hPU);
        printf("%4.0f  %4d  %7.2f%%  %8.2f%%  %8.2f%%   %+6.2f%%  %+8.2f%%  %+6.2f%%\n",
               pMid, bins[i].n, 100*sR, 100*sPS, 100*sPU, 100*bR, 100*bPS, 100*bPU);
        of << pMid << "," << bins[i].n << ","
           << 100*sR  << "," << 100*sPS << "," << 100*sPU << ","
           << 100*bR  << "," << 100*bPS << "," << 100*bPU << "\n";
    }
    of.close();
    printf("\nCSV: %s\n", csv.Data());
}
