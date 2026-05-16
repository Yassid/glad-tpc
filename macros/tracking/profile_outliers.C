// Profile what makes a seed-fit outlier: examine R_fit/R_truth tail
// (|ratio - 1| > 0.2) and compare to good fits to find what predictors
// (chord, hit count, MC angle, presence of 3He/electrons in event, etc.)
// flag the catastrophes.
void profile_outliers(TString tag = "goodevt2k")
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    TString wd = gSystem->Getenv("VMCWORKDIR");
    TFile fS(wd + "/glad-tpc/macros/sim/Prototype/sim_" + tag + ".root");
    TFile fT("output_tracking_" + tag + ".root");
    auto* tS = (TTree*)fS.Get("evt"); auto* tT = (TTree*)fT.Get("evt");
    auto* mc = new TClonesArray("R3BMCTrack");
    auto* trks = new TClonesArray("R3BGTPCTrackData");
    auto* pts  = new TClonesArray("R3BGTPCPoint");
    tS->SetBranchAddress("MCTrack", &mc);
    tS->SetBranchAddress("GTPCPoint", &pts);
    tT->SetBranchAddress("GTPCTrackData", &trks);

    int nGood = 0, nBad = 0;
    auto* hGoodN  = new TH1F("hGoodN" ,";N hits; ", 40, 0, 200);
    auto* hBadN   = new TH1F("hBadN"  ,";N hits; ", 40, 0, 200);
    auto* hGoodCh = new TH1F("hGoodCh",";chord (cm); ", 30, 0, 30);
    auto* hBadCh  = new TH1F("hBadCh" ,";chord (cm); ", 30, 0, 30);
    auto* hGoodP  = new TH1F("hGoodP" ,";p_{MC} (MeV/c);", 40, 0, 1000);
    auto* hBadP   = new TH1F("hBadP"  ,";p_{MC} (MeV/c);", 40, 0, 1000);
    auto* hGoodPdg = new TH1F("hGoodPdg",";n other PDG in gas;",10,0,10);
    auto* hBadPdg  = new TH1F("hBadPdg" ,";n other PDG in gas;",10,0,10);
    auto* hGoodH2 = new TH1F("hGoodH2",";N hits / chord;",40,0,30);
    auto* hBadH2  = new TH1F("hBadH2" ,";N hits / chord;",40,0,30);
    auto* hGoodNcl = new TH1F("hGoodNcl",";n clusters in event;",6,0.5,6.5);
    auto* hBadNcl  = new TH1F("hBadNcl" ,";n clusters in event;",6,0.5,6.5);

    for (Long64_t i = 0; i < std::min(tS->GetEntries(), tT->GetEntries()); ++i) {
        tS->GetEntry(i); tT->GetEntry(i);
        if (trks->GetEntries() == 0) continue;
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        auto& hits = tr->GetHitArray();
        if (hits.size() < 10) continue;
        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mc->At(j);
            if (m->GetMotherId()==-1 && m->GetPdgCode()==-211) { pi = m; break; }
        }
        if (!pi) continue;
        double pT = std::sqrt(pi->GetPx()*pi->GetPx() + pi->GetPz()*pi->GetPz())*1000;
        double pMC = std::sqrt(pi->GetPx()*pi->GetPx()+pi->GetPy()*pi->GetPy()+pi->GetPz()*pi->GetPz())*1000;
        double R_truth = pT / 6.0;
        double R_fit = tr->GetGeoRadius();
        if (!std::isfinite(R_fit) || R_fit <= 0) continue;
        double xmn=1e9,xmx=-1e9,zmn=1e9,zmx=-1e9;
        for (auto& h : hits) { xmn=std::min(xmn,(double)h.GetX()); xmx=std::max(xmx,(double)h.GetX());
                               zmn=std::min(zmn,(double)h.GetZ()); zmx=std::max(zmx,(double)h.GetZ()); }
        double chord = std::sqrt((xmx-xmn)*(xmx-xmn)+(zmx-zmn)*(zmx-zmn));
        if (chord < 12) continue;  // restrict to relevant region

        std::set<int> otherPdgs;
        for (int j = 0; j < pts->GetEntries(); ++j) {
            auto* p = (R3BGTPCPoint*)pts->At(j);
            if (p->GetPDGCode() != -211) otherPdgs.insert(p->GetPDGCode());
        }
        int nOther = otherPdgs.size();
        int nCl = trks->GetEntries();

        double ratio = R_fit / R_truth;
        bool isBad = std::abs(ratio - 1.0) > 0.15;
        if (isBad) {
            ++nBad;
            hBadN->Fill(hits.size()); hBadCh->Fill(chord); hBadP->Fill(pMC);
            hBadPdg->Fill(nOther); hBadH2->Fill(hits.size()/chord); hBadNcl->Fill(nCl);
        } else {
            ++nGood;
            hGoodN->Fill(hits.size()); hGoodCh->Fill(chord); hGoodP->Fill(pMC);
            hGoodPdg->Fill(nOther); hGoodH2->Fill(hits.size()/chord); hGoodNcl->Fill(nCl);
        }
    }
    printf("Good: %d   Bad: %d  (%.0f%% catastrophic fits, |ratio-1|>0.15)\n", nGood, nBad, 100.0*nBad/(nGood+nBad));
    auto cmp = [](const char* lbl, TH1F* g, TH1F* b){
        printf("  %-20s good mean=%.2f rms=%.2f  | bad mean=%.2f rms=%.2f\n",
               lbl, g->GetMean(), g->GetRMS(), b->GetMean(), b->GetRMS());
    };
    cmp("N hits",         hGoodN , hBadN );
    cmp("chord (cm)",     hGoodCh, hBadCh);
    cmp("p_MC (MeV/c)",   hGoodP , hBadP );
    cmp("n other PDG",    hGoodPdg, hBadPdg);
    cmp("hit density",    hGoodH2, hBadH2);
    cmp("n clusters",     hGoodNcl, hBadNcl);

    TCanvas c("c","",1500,900); c.Divide(3,2);
    auto draw2 = [&](int idx, TH1F* g, TH1F* b){
        c.cd(idx);
        g->SetLineColor(kGreen+2); g->Scale(1.0/std::max(1.0,g->Integral()));
        b->SetLineColor(kRed+1);   b->Scale(1.0/std::max(1.0,b->Integral()));
        double mx = std::max(g->GetMaximum(), b->GetMaximum());
        g->SetMaximum(mx*1.2);
        g->Draw(); b->Draw("SAME");
        auto* l = new TLegend(0.55,0.7,0.95,0.88);
        l->AddEntry(g, Form("good (%d)", nGood), "l");
        l->AddEntry(b, Form("bad (%d)", nBad), "l");
        l->Draw();
    };
    draw2(1, hGoodN, hBadN);
    draw2(2, hGoodCh, hBadCh);
    draw2(3, hGoodP, hBadP);
    draw2(4, hGoodPdg, hBadPdg);
    draw2(5, hGoodH2, hBadH2);
    draw2(6, hGoodNcl, hBadNcl);
    c.SaveAs("outlier_profile.png");
}
