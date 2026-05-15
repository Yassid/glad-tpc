/// @file analyze_bias_origin.C
/// @brief Find where the +19 % bias enters: at the first-cluster fit, or
///        during the helix-POCA back-extrapolation to the vertex.
///
/// For each long-chord event (chord ≥ 16 cm):
///   res_xtr = (p_fit_xtr - p_MC_at_first_cluster) / p_MC_at_first_cluster
///   res_vtx = (p_fit_vtx - p_MC_at_vertex) / p_MC_at_vertex
///
/// p_MC_at_vertex is the box generator vertex momentum from R3BMCTrack.
/// p_MC_at_first_cluster: best we can do is the MC primary point closest to
/// the first cluster; report both with and without that correction.

void analyze_bias_origin()
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gSystem->Load("libR3BOpenKF");
    gStyle->SetOptStat(0);

    TString workDir = gSystem->Getenv("VMCWORKDIR");
    TFile fSim(workDir + "/glad-tpc/macros/sim/Prototype/sim.root");
    TFile fTrk(workDir + "/glad-tpc/macros/tracking/output_tracking.root");
    TFile fUKF(workDir + "/glad-tpc/macros/tracking/output_ukf.root");
    auto* tSim = (TTree*)fSim.Get("evt");
    auto* tTrk = (TTree*)fTrk.Get("evt");
    auto* tUKF = (TTree*)fUKF.Get("evt");
    auto* mc = new TClonesArray("R3BMCTrack");
    auto* pts = new TClonesArray("R3BGTPCPoint");
    auto* trks = new TClonesArray("R3BGTPCTrackData");
    auto* fits = new TClonesArray("R3BGTPCFittedTrackData");
    tSim->SetBranchAddress("MCTrack", &mc);
    tSim->SetBranchAddress("GTPCPoint", &pts);
    tTrk->SetBranchAddress("GTPCTrackData", &trks);
    tUKF->SetBranchAddress("GTPCFittedTrackData", &fits);

    auto* hXtr = new TH1F("hXtr", "p_fit_xtr vs p_MC_first_pt;(p_{fit,xtr}-p_{MC,first pt})/p_{MC,first pt};events", 100, -1.0, 1.5);
    auto* hVtx = new TH1F("hVtx", "p_fit_vtx vs p_MC_vertex;(p_{fit,vtx}-p_{MC,vtx})/p_{MC,vtx};events", 100, -1.0, 1.5);
    auto* hDeltaP_MC = new TH1F("hDeltaPMC", "MC: p at first GTPC point vs at vertex;(p_{first}-p_{vtx})/p_{vtx};events", 100, -1.0, 1.5);
    auto* hDeltaP_fit = new TH1F("hDeltaPFit", "Fit: back-extrap correction;(p_{vtx,fit}-p_{xtr,fit})/p_{xtr,fit};events", 100, -1.0, 1.5);

    int nLong = 0;
    double sumVtx = 0, sumXtr = 0, sumMC = 0, sumFit = 0;
    double sumKasaR = 0, sumPMC = 0, sumPxtr = 0;
    int kasaTooHigh = 0, kasaTooLow = 0;
    for (Long64_t i = 0; i < std::min({ tSim->GetEntries(), tTrk->GetEntries(), tUKF->GetEntries() }); ++i) {
        tSim->GetEntry(i); tTrk->GetEntry(i); tUKF->GetEntry(i);
        if (fits->GetEntries() == 0 || trks->GetEntries() == 0) continue;
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        auto& cls = *tr->GetHitClusterArray();
        if (cls.size() < 3) continue;
        // chord length filter — focus on the working population
        double zmn = 1e9, zmx = -1e9, xmn = 1e9, xmx = -1e9;
        for (auto& c : cls) { zmn = std::min(zmn, (double)c.GetZ()); zmx = std::max(zmx, (double)c.GetZ());
                              xmn = std::min(xmn, (double)c.GetX()); xmx = std::max(xmx, (double)c.GetX()); }
        const double chord_cm = std::sqrt((zmx - zmn) * (zmx - zmn) + (xmx - xmn) * (xmx - xmn));
        if (chord_cm < 16.0) continue;

        // MC primary pion at vertex
        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mc->At(j);
            if (m->GetMotherId() == -1 && std::abs(m->GetPdgCode()) == 211) { pi = m; break; }
        }
        if (!pi) continue;
        const double pmc_vtx = std::sqrt(pi->GetPx() * pi->GetPx() + pi->GetPy() * pi->GetPy() + pi->GetPz() * pi->GetPz()) * 1000;

        // First MC π point inside chamber (lowest z, |PDG|==211).
        double pmc_first = 0; double zmin_pt = 1e9;
        const double mass_pi_MeV = 139.57;
        for (int j = 0; j < pts->GetEntries(); ++j) {
            auto* p = (R3BGTPCPoint*)pts->At(j);
            if (std::abs(p->GetPDGCode()) != 211) continue;
            if (p->GetZ() < zmin_pt) {
                zmin_pt = p->GetZ();
                const double ke_MeV = p->GetKineticEnergy() * 1000.0; // GeV → MeV
                pmc_first = std::sqrt(ke_MeV * ke_MeV + 2 * ke_MeV * mass_pi_MeV);
            }
        }
        if (pmc_first <= 0) continue;

        auto* fit = (R3BGTPCFittedTrackData*)fits->At(0);
        const double ke_vtx = fit->GetKinematics().kineticEnergy;
        const double ke_xtr = fit->GetKinematicsXtr().kineticEnergy;
        if (ke_vtx <= 0 || ke_xtr <= 0) continue;
        const double p_vtx = std::sqrt(ke_vtx * ke_vtx + 2 * ke_vtx * 139.57);
        const double p_xtr = std::sqrt(ke_xtr * ke_xtr + 2 * ke_xtr * 139.57);

        hVtx->Fill((p_vtx - pmc_vtx) / pmc_vtx);
        hXtr->Fill((p_xtr - pmc_first) / pmc_first);
        hDeltaP_MC->Fill((pmc_first - pmc_vtx) / pmc_vtx);
        hDeltaP_fit->Fill((p_vtx - p_xtr) / p_xtr);
        sumVtx += (p_vtx - pmc_vtx) / pmc_vtx;
        sumXtr += (p_xtr - pmc_first) / pmc_first;
        sumMC += (pmc_first - pmc_vtx) / pmc_vtx;
        sumFit += (p_vtx - p_xtr) / p_xtr;
        // Compare Kasa-implied p_T to MC p_T_at_first_cluster (both in x-z plane,
        // since B=ŷ). Kasa: p_T_seed = 0.3·|B|·R. MC: p_T = √(px²+pz²) at first hit.
        const double R_cm = tr->GetGeoRadius();
        const double pT_seed = 0.3 * 2.0 * (R_cm * 0.01) * 1000.0; // MeV/c
        sumKasaR += pT_seed;
        sumPMC += pmc_vtx;
        sumPxtr += p_xtr;
        if (pT_seed > pmc_first) ++kasaTooHigh; else ++kasaTooLow;
        ++nLong;
    }
    std::cout << "long-chord events: " << nLong << "\n";
    if (nLong > 0) {
        std::cout << "MEAN over events:\n";
        std::printf("  (p_fit_vtx - p_MC_vtx)/p_MC_vtx         = %+.3f\n", sumVtx / nLong);
        std::printf("  (p_fit_xtr - p_MC_first)/p_MC_first     = %+.3f\n", sumXtr / nLong);
        std::printf("  (p_MC_first - p_MC_vtx)/p_MC_vtx (Eloss)= %+.3f\n", sumMC / nLong);
        std::printf("  (p_fit_vtx - p_fit_xtr)/p_fit_xtr       = %+.3f\n", sumFit / nLong);
        std::printf("  mean Kasa-implied p_T seed              = %.0f MeV/c\n", sumKasaR / nLong);
        std::printf("  mean p_MC_vtx                           = %.0f MeV/c\n", sumPMC / nLong);
        std::printf("  mean p_fit_xtr                          = %.0f MeV/c\n", sumPxtr / nLong);
        std::printf("  Kasa seed > MC truth in %d/%d events    \n", kasaTooHigh, nLong);
    }

    auto* c = new TCanvas("c", "bias origin", 1500, 900);
    c->Divide(2, 2);
    auto report = [](TH1F* h, const char* label) {
        h->Fit("gaus", "Q", "", -0.3, 0.3);
        auto* g = h->GetFunction("gaus");
        std::printf("  %-40s  N=%d  mean=%+6.3f  rms=%5.3f  gauss(σ=%5.3f μ=%+6.3f)\n",
                    label, (int)h->GetEntries(), h->GetMean(), h->GetRMS(),
                    g ? g->GetParameter(2) : 0, g ? g->GetParameter(1) : 0);
    };
    c->cd(1); hVtx->Draw(); report(hVtx,    "p_fit_vtx vs p_MC_vtx           ");
    c->cd(2); hXtr->Draw(); report(hXtr,    "p_fit_xtr vs p_MC_first_clust   ");
    c->cd(3); hDeltaP_MC->Draw(); report(hDeltaP_MC, "MC: p_first - p_vtx (energy loss)");
    c->cd(4); hDeltaP_fit->Draw(); report(hDeltaP_fit, "Fit: p_vtx - p_xtr (back-extrap) ");
    c->SaveAs("bias_origin.png");
}
