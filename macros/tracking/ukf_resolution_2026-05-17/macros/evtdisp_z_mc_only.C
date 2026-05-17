// MC-only event display: shows just the Geant4 GTPCPoints (π±) + vertex +
// initial momentum arrow over the pad-plane outline. No reco hits, no fitted
// circle, no truth-arc overlay — the underlying truth trajectory in isolation.
void evtdisp_z_mc_only()
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TString wd = gSystem->Getenv("VMCWORKDIR");
    TFile fReco(wd + "/glad-tpc/macros/reco/output_reco_goodevt2k.root");
    TFile fTrk("output_tracking_goodevt2k.root");
    TFile fSim(wd + "/glad-tpc/macros/sim/Prototype/sim_goodevt2k.root");
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

    // Same event selection as evtdisp_z.C so panels line up between the two
    // figures and the reader can compare "MC only" vs "MC + reco + fit".
    std::vector<Long64_t> evts;
    for (Long64_t i = 0; i < tR->GetEntries() && (int)evts.size() < 9; ++i) {
        tR->GetEntry(i); tT->GetEntry(i);
        if (trks->GetEntries() != 1) continue;
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        auto& hits = tr->GetHitArray();
        if (hits.size() < 40) continue;
        double xmn=1e9,xmx=-1e9,zmn=1e9,zmx=-1e9;
        for (auto& h : hits) {
            xmn=std::min(xmn,(double)h.GetX()); xmx=std::max(xmx,(double)h.GetX());
            zmn=std::min(zmn,(double)h.GetZ()); zmx=std::max(zmx,(double)h.GetZ());
        }
        if (std::hypot(xmx-xmn, zmx-zmn) < 10.0) continue;
        evts.push_back(i);
    }
    printf("selected %lu events\n", evts.size());

    int cols=3, rows=(evts.size()+cols-1)/cols;
    TCanvas* c = new TCanvas("c","", 1800, 400*rows);
    c->Divide(cols, rows, 0.006, 0.025);
    int ip = 0;
    for (auto i : evts) {
        ++ip;
        c->cd(ip);
        gPad->SetLeftMargin(0.10);
        gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.10);
        gPad->SetBottomMargin(0.14);
        tS->GetEntry(i);

        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mc->At(j);
            if (m->GetMotherId()==-1 && std::abs(m->GetPdgCode())==211) { pi = m; break; }
        }
        if (!pi) continue;
        const double vxL = pi->GetStartX() - offX;
        const double vzL = pi->GetStartZ() - offZ;
        const double px = pi->GetPx(), pz = pi->GetPz();
        const double mag = std::hypot(px, pz);
        const double pTot = std::sqrt(pi->GetPx()*pi->GetPx() + pi->GetPy()*pi->GetPy() + pi->GetPz()*pi->GetPz()) * 1000;

        double zlo = std::min(padZmin, vzL), zhi = std::max(padZmax, vzL);
        double xlo = std::min(padXmin, vxL), xhi = std::max(padXmax, vxL);
        for (int j = 0; j < pts->GetEntries(); ++j) {
            auto* p = (R3BGTPCPoint*)pts->At(j);
            if (std::abs(p->GetPDGCode()) != 211) continue;
            double xl = p->GetX() - offX, zl = p->GetZ() - offZ;
            zlo = std::min(zlo, zl); zhi = std::max(zhi, zl);
            xlo = std::min(xlo, xl); xhi = std::max(xhi, xl);
        }
        double padMargin = 1.0;
        auto* frame = gPad->DrawFrame(zlo-padMargin, xlo-padMargin, zhi+padMargin, xhi+padMargin);
        frame->SetTitle(TString::Format("evt %lld    p_{MC}=%.0f MeV/c", (long long)i, pTot));
        frame->GetXaxis()->SetTitle("z_{local} (beam) [cm]");
        frame->GetYaxis()->SetTitle("x_{local} [cm]");

        auto* box = new TBox(padZmin, padXmin, padZmax, padXmax);
        box->SetLineColor(kGray+2); box->SetLineWidth(2); box->SetLineStyle(7); box->SetFillStyle(0);
        box->Draw("L SAME");

        auto* gM = new TGraph(); int miN = 0;
        for (int j = 0; j < pts->GetEntries(); ++j) {
            auto* p = (R3BGTPCPoint*)pts->At(j);
            if (std::abs(p->GetPDGCode()) != 211) continue;
            gM->SetPoint(miN++, p->GetZ() - offZ, p->GetX() - offX);
        }
        if (miN > 0) {
            gM->SetMarkerStyle(29); gM->SetMarkerSize(1.2); gM->SetMarkerColor(kBlack);
            gM->Draw("P SAME");
        }

        auto* gV = new TGraph(); gV->SetPoint(0, vzL, vxL);
        gV->SetMarkerStyle(33); gV->SetMarkerSize(2.0); gV->SetMarkerColor(kMagenta+2);
        gV->Draw("P SAME");
        double az = vzL + 2.5*pz/mag, ax = vxL + 2.5*px/mag;
        auto* arr = new TArrow(vzL, vxL, az, ax, 0.018, "|>");
        arr->SetLineColor(kMagenta+2); arr->SetLineWidth(2); arr->Draw();
    }
    c->SaveAs("evtdisp_z_mc_only.png");
}
