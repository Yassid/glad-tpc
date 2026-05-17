// σ_p/p vs (x,z) chord length on the good_evt sim.
// Same per-bin structure as scan_goodevt_pbins.C but binned by chord
// (bounding-box hypotenuse in xz) instead of MC momentum. Tells you
// where the lever arm in the bending plane starts to matter for the
// fit precision.
void scan_goodevt_chordbins(TString tag = "goodevt2k",
                            double cMin = 0, double cMax = 30, double cStep = 3)
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

    const int nBins = (int)std::round((cMax - cMin) / cStep);
    struct Bin { TH1F* hR; TH1F* hPS; TH1F* hPU; int n = 0; double pSum = 0; };
    std::vector<Bin> bins(nBins);
    for (int i = 0; i < nBins; ++i) {
        double lo = cMin + i*cStep, hi = lo + cStep;
        bins[i].hR  = new TH1F(Form("hR_%d",  i), Form("c in [%g,%g];R/R_{tru}-1;", lo, hi), 80, -1, 1);
        bins[i].hPS = new TH1F(Form("hPS_%d", i), Form("c in [%g,%g];p_{seed}/p_{MC}-1;", lo, hi), 80, -1, 1);
        bins[i].hPU = new TH1F(Form("hPU_%d", i), Form("c in [%g,%g];p_{UKF}/p_{MC}-1;",  lo, hi), 80, -1, 1);
    }

    for (Long64_t i = 0; i < std::min({tS->GetEntries(), tT->GetEntries(), tU->GetEntries()}); ++i) {
        tS->GetEntry(i); tT->GetEntry(i); tU->GetEntry(i);
        if (trks->GetEntries() == 0) continue;
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        auto& hits = tr->GetHitArray();
        if (hits.size() < 10) continue;
        double xmn=1e9,xmx=-1e9,zmn=1e9,zmx=-1e9;
        for (auto& h : hits) { xmn=std::min(xmn,(double)h.GetX()); xmx=std::max(xmx,(double)h.GetX());
                               zmn=std::min(zmn,(double)h.GetZ()); zmx=std::max(zmx,(double)h.GetZ()); }
        double chord = std::sqrt((xmx-xmn)*(xmx-xmn)+(zmx-zmn)*(zmx-zmn));
        if (chord < cMin || chord >= cMax) continue;
        int b = std::min(nBins-1, std::max(0, (int)((chord - cMin) / cStep)));

        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* mt = (R3BMCTrack*)mc->At(j);
            if (mt->GetMotherId()==-1 && mt->GetPdgCode()==-211) { pi = mt; break; }
        }
        if (!pi) continue;
        double pT_MC = std::sqrt(pi->GetPx()*pi->GetPx() + pi->GetPz()*pi->GetPz()) * 1000;
        double p_MC  = std::sqrt(pi->GetPx()*pi->GetPx() + pi->GetPy()*pi->GetPy() + pi->GetPz()*pi->GetPz()) * 1000;
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
        bins[b].pSum += p_MC;
    }

    TString csv = wd + "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17/scan_goodevt_chordbins.csv";
    std::ofstream of(csv.Data());
    of << "chord_cm,N_fit,p_avg_MeV,sigma_R_pct,sigma_p_seed_pct,sigma_p_UKF_pct,bias_R_pct,bias_p_seed_pct,bias_p_UKF_pct\n";
    printf("chord  N_fit  <p_MC>   sigma_R    sigma_p_seed  sigma_p_UKF  bias_R   bias_p_seed  bias_p_UKF\n");
    auto gauss = [](TH1F* h) -> std::pair<double,double> {
        if (h->GetEntries() < 5) return {std::nan(""), std::nan("")};
        h->Fit("gaus", "Q", "", -0.4, 0.4);
        auto* f = h->GetFunction("gaus");
        if (!f) return {std::nan(""), std::nan("")};
        return { f->GetParameter(1), f->GetParameter(2) };
    };
    for (int i = 0; i < nBins; ++i) {
        double cMid = cMin + (i + 0.5)*cStep;
        double pAvg = (bins[i].n > 0) ? bins[i].pSum / bins[i].n : 0;
        auto [bR, sR]   = gauss(bins[i].hR);
        auto [bPS, sPS] = gauss(bins[i].hPS);
        auto [bPU, sPU] = gauss(bins[i].hPU);
        printf("%4.1f  %4d  %6.1f   %7.2f%%   %8.2f%%   %8.2f%%   %+6.2f%%  %+8.2f%%  %+6.2f%%\n",
               cMid, bins[i].n, pAvg, 100*sR, 100*sPS, 100*sPU, 100*bR, 100*bPS, 100*bPU);
        of << cMid << "," << bins[i].n << "," << pAvg << ","
           << 100*sR  << "," << 100*sPS << "," << 100*sPU << ","
           << 100*bR  << "," << 100*bPS << "," << 100*bPU << "\n";
    }
    of.close();
    printf("\nCSV: %s\n", csv.Data());
}
