// Measure σ_p/p from R3BGTPCFittedTrackData (UKF output) against MC truth.
// Optional: pass the trk file (default goodevt) and sim file as args.
// useXtr: when true use Kinematics at first cluster (pre-back-extrap); when
// false use the vertex Kinematics. The Xtr value is the UKF's actual fit
// observable; the vertex one adds a back-extrapolation through 7 cm of gas.
void measure_ukf(TString ukfFile = "output_ukf_goodevt.root",
                 TString simFile = "Prototype/sim_goodevt.root",
                 TString tag = "goodevt",
                 int particlePDG = -211,
                 bool useXtr = true)
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    TString wd = gSystem->Getenv("VMCWORKDIR");
    TFile fSim(wd + "/glad-tpc/macros/sim/" + simFile);
    TFile fUkf(ukfFile);
    // Open the tracking file too so the chord cut uses raw-hit span
    // (consistent with vertex_refit.C). UKF smoothed-position chord is
    // ~89 % of raw-hit chord because PRA cluster centroids sit inside the
    // hit envelope — using SP-chord would silently misbin ~half the
    // borderline long-chord events.
    TString trkFile = ukfFile; trkFile.ReplaceAll("output_ukf_", "output_tracking_");
    TFile fTrk(trkFile);
    auto* tS = (TTree*)fSim.Get("evt");
    auto* tU = (TTree*)fUkf.Get("evt");
    auto* tT = (TTree*)fTrk.Get("evt");
    auto* mc     = new TClonesArray("R3BMCTrack");
    auto* fitted = new TClonesArray("R3BGTPCFittedTrackData");
    auto* trks   = new TClonesArray("R3BGTPCTrackData");
    tS->SetBranchAddress("MCTrack", &mc);
    tU->SetBranchAddress("GTPCFittedTrackData", &fitted);
    if (tT) tT->SetBranchAddress("GTPCTrackData", &trks);

    const double m_pi = 139.57039; // MeV
    auto p_from_KE = [&](double ke) { return std::sqrt((ke + m_pi)*(ke + m_pi) - m_pi*m_pi); };

    int nL=0, nM=0, nS=0, nTot=0;
    double sLL=0, sLL2=0, sLM=0, sLM2=0;
    auto* hL = new TH1F("hL", ";p_{fit}/p_{MC}-1;", 80, -1, 1);
    auto* hM = new TH1F("hM", ";p_{fit}/p_{MC}-1;", 80, -1, 1);
    auto* hAll = new TH1F("hAll", ";p_{fit}/p_{MC}-1;", 80, -1, 1);
    Long64_t nEvt = std::min(tS->GetEntries(), tU->GetEntries());
    if (tT) nEvt = std::min(nEvt, tT->GetEntries());
    for (Long64_t i = 0; i < nEvt; ++i) {
        tS->GetEntry(i); tU->GetEntry(i);
        if (tT) tT->GetEntry(i);
        if (fitted->GetEntries() == 0) continue;
        auto* ft = (R3BGTPCFittedTrackData*)fitted->At(0);
        if (!ft->IsConverged()) continue;
        const auto& kin = useXtr ? ft->GetKinematicsXtr() : ft->GetKinematics();
        if (!std::isfinite(kin.kineticEnergy) || kin.kineticEnergy <= 0) continue;
        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* mt = (R3BMCTrack*)mc->At(j);
            if (mt->GetMotherId()==-1 && mt->GetPdgCode()==particlePDG) { pi = mt; break; }
        }
        if (!pi) continue;
        double pMC = std::sqrt(pi->GetPx()*pi->GetPx() + pi->GetPy()*pi->GetPy() + pi->GetPz()*pi->GetPz()) * 1000;
        double pTMC = std::sqrt(pi->GetPx()*pi->GetPx() + pi->GetPz()*pi->GetPz()) * 1000;
        double pFit  = p_from_KE(kin.kineticEnergy);
        // kin.theta and kin.phi are ROOT spherical (theta from z-axis, phi
        // in xy-plane), inherited from XYZVector::Theta()/Phi() in the
        // UKF state. R3B B-field is along y, so pT must be transverse to y:
        // pT = sqrt(px^2 + pz^2) = |p| * sqrt(1 - sin^2(theta)*sin^2(phi)).
        const double sTh = std::sin(kin.theta);
        const double sPhi = std::sin(kin.phi);
        double pTFit = pFit * std::sqrt(std::max(0.0, 1.0 - sTh*sTh*sPhi*sPhi));
        if (pMC <= 0) continue;
        // Print first 5 entries in both metrics for diagnosis
        static int dumped = 0;
        if (dumped < 5) {
            ++dumped;
            printf("  evt %lld: p_MC=%.0f p_fit=%.0f  pT_MC=%.0f pT_fit=%.0f  th_MC=%.2f th_fit=%.2f\n",
                   i, pMC, pFit, pTMC, pTFit, std::acos(pi->GetPy()/(pMC/1000)), kin.theta);
        }
        // Chord from raw hits in R3BGTPCTrackData (consistent with
        // vertex_refit.C). Fall back to UKF smoothed positions if the
        // tracking file isn't available.
        double xmn=1e9,xmx=-1e9,zmn=1e9,zmx=-1e9;
        if (tT && trks->GetEntries() > 0) {
            auto& hits = ((R3BGTPCTrackData*)trks->At(0))->GetHitArray();
            if (hits.size() < 10) continue;
            for (auto& h : hits) { xmn=std::min(xmn,(double)h.GetX()); xmx=std::max(xmx,(double)h.GetX());
                                    zmn=std::min(zmn,(double)h.GetZ()); zmx=std::max(zmx,(double)h.GetZ()); }
        } else {
            const auto& sp = ft->GetSmoothedPositions();
            if (sp.size() < 10) continue;
            for (auto& q : sp) { xmn=std::min(xmn,q.X()); xmx=std::max(xmx,q.X());
                                  zmn=std::min(zmn,q.Z()); zmx=std::max(zmx,q.Z()); }
        }
        double chord = std::sqrt((xmx-xmn)*(xmx-xmn)+(zmx-zmn)*(zmx-zmn));
        double r  = pFit/pMC - 1;
        double rT = pTFit/pTMC - 1;
        nTot++;
        hAll->Fill(r);
        if (chord >= 16) { ++nL; sLL += std::log(pFit/pMC); sLL2 += std::log(pFit/pMC)*std::log(pFit/pMC); hL->Fill(r); }
        else if (chord >= 12) { ++nM; sLM += std::log(pFit/pMC); sLM2 += std::log(pFit/pMC)*std::log(pFit/pMC); hM->Fill(r); }
        else ++nS;
        static auto* hLpT = new TH1F("hLpT",";p_{T,fit}/p_{T,MC}-1 (long);",80,-1,1);
        if (chord >= 16) hLpT->Fill(rT);
        static auto* hMpT = new TH1F("hMpT",";p_{T,fit}/p_{T,MC}-1 (mid);",80,-1,1);
        if (chord >= 12 && chord < 16) hMpT->Fill(rT);
    }
    printf("=== UKF σ_p/p  [%s] ===\n", tag.Data());
    printf("  fitted+converged: %d   (long=%d  mid=%d  short=%d)\n", nTot, nL, nM, nS);
    if (nL > 0) {
        double m=sLL/nL, s=std::sqrt(sLL2/nL-m*m);
        printf("  chord >=16  N=%d  median=%.3f  log-rms=%.3f\n", nL, std::exp(m), s);
        hL->Fit("gaus","Q","",-0.4,0.4);
        if (auto* f = hL->GetFunction("gaus"))
            printf("    p_total Gauss core:  mean=%+.3f  sigma=%.3f\n", f->GetParameter(1), f->GetParameter(2));
        auto* hLpT = (TH1F*)gROOT->FindObject("hLpT");
        if (hLpT) { hLpT->Fit("gaus","Q","",-0.4,0.4);
            if (auto* f = hLpT->GetFunction("gaus"))
                printf("    p_T     Gauss core:  mean=%+.3f  sigma=%.3f\n", f->GetParameter(1), f->GetParameter(2)); }
    }
    if (nM > 0) {
        double m=sLM/nM, s=std::sqrt(sLM2/nM-m*m);
        printf("  chord 12-16 N=%d  median=%.3f  log-rms=%.3f\n", nM, std::exp(m), s);
        hM->Fit("gaus","Q","",-0.4,0.4);
        if (auto* f = hM->GetFunction("gaus"))
            printf("    p_total Gauss core:  mean=%+.3f  sigma=%.3f\n", f->GetParameter(1), f->GetParameter(2));
        auto* hMpT = (TH1F*)gROOT->FindObject("hMpT");
        if (hMpT) { hMpT->Fit("gaus","Q","",-0.4,0.4);
            if (auto* f = hMpT->GetFunction("gaus"))
                printf("    p_T     Gauss core:  mean=%+.3f  sigma=%.3f\n", f->GetParameter(1), f->GetParameter(2)); }
    }
    TCanvas c("c","",1500,500); c.Divide(3,1);
    c.cd(1); hAll->Draw(); c.cd(2); hM->Draw(); c.cd(3); hL->Draw();
    c.SaveAs(TString("ukf_sigp_") + tag + ".png");
}
