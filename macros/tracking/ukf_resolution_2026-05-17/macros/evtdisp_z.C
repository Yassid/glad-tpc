void evtdisp_z()
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
        tR->GetEntry(i); tT->GetEntry(i); tS->GetEntry(i);
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        auto& hits = tr->GetHitArray();

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
        const double pT = mag * 1000;
        const double R_truth = pT / (3.0 * 2.0);
        const double R_fit = tr->GetGeoRadius();
        auto cen = tr->GetGeoCenter();

        double nx = -pz/mag, nz = px/mag;
        double cxA = vxL + R_truth*nx, czA = vzL + R_truth*nz;
        double cxB = vxL - R_truth*nx, czB = vzL - R_truth*nz;
        double sumA = 0, sumB = 0;
        for (int j = 0; j < pts->GetEntries(); ++j) {
            auto* p = (R3BGTPCPoint*)pts->At(j);
            if (std::abs(p->GetPDGCode()) != 211) continue;
            double xl = p->GetX() - offX, zl = p->GetZ() - offZ;
            sumA += std::abs(std::hypot(xl-cxA, zl-czA) - R_truth);
            sumB += std::abs(std::hypot(xl-cxB, zl-czB) - R_truth);
        }
        double cxT = (sumA <= sumB) ? cxA : cxB;
        double czT = (sumA <= sumB) ? czA : czB;

        double zlo = std::min(padZmin, vzL), zhi = std::max(padZmax, vzL);
        double xlo = std::min(padXmin, vxL), xhi = std::max(padXmax, vxL);
        for (auto& h : hits) {
            zlo = std::min(zlo,(double)h.GetZ()); zhi = std::max(zhi,(double)h.GetZ());
            xlo = std::min(xlo,(double)h.GetX()); xhi = std::max(xhi,(double)h.GetX());
        }
        double padMargin = 1.0;
        auto* frame = gPad->DrawFrame(zlo-padMargin, xlo-padMargin, zhi+padMargin, xhi+padMargin);
        frame->SetTitle(TString::Format("evt %lld    R_{fit}=%.0f cm    R_{tru}=%.0f cm    ratio=%.2f",
                                        (long long)i, R_fit, R_truth, R_truth>0?R_fit/R_truth:0));
        frame->GetXaxis()->SetTitle("z_{local} (beam) [cm]");
        frame->GetYaxis()->SetTitle("x_{local} [cm]");

        auto* box = new TBox(padZmin, padXmin, padZmax, padXmax);
        box->SetLineColor(kGray+2); box->SetLineWidth(2); box->SetLineStyle(7); box->SetFillStyle(0);
        box->Draw("L SAME");

        // Truth arc through vertex and hits
        auto phiOf = [&](double zl, double xl){ return std::atan2(xl-cxT, zl-czT); };
        double phiVx = phiOf(vzL, vxL);
        double phiH1 = phiOf(hits.front().GetZ(), hits.front().GetX());
        double phiHN = phiOf(hits.back().GetZ(), hits.back().GetX());
        double pmin = std::min({phiVx, phiH1, phiHN});
        double pmax = std::max({phiVx, phiH1, phiHN});
        if (pmax - pmin > TMath::Pi()) {
            auto sh = [](double p){ return p < 0 ? p + 2*TMath::Pi() : p; };
            double a = sh(phiVx), b = sh(phiH1), c2 = sh(phiHN);
            pmin = std::min({a,b,c2}); pmax = std::max({a,b,c2});
        }
        double dphi = pmax-pmin; pmin -= 0.04*dphi; pmax += 0.04*dphi;
        auto* gTraj = new TGraph();
        for (int k = 0; k <= 300; ++k) {
            double p = pmin + (pmax-pmin)*k/300;
            gTraj->SetPoint(k, czT + R_truth*std::cos(p), cxT + R_truth*std::sin(p));
        }
        gTraj->SetLineColor(kGreen+2); gTraj->SetLineWidth(2); gTraj->SetLineStyle(2);
        gTraj->Draw("L SAME");

        auto* gH = new TGraph(); int gi = 0;
        for (auto& h : hits) gH->SetPoint(gi++, h.GetZ(), h.GetX());
        gH->SetMarkerStyle(20); gH->SetMarkerSize(0.55); gH->SetMarkerColor(kBlue+1);
        gH->Draw("P SAME");

        auto* gM = new TGraph(); int miN = 0;
        for (int j = 0; j < pts->GetEntries(); ++j) {
            auto* p = (R3BGTPCPoint*)pts->At(j);
            if (std::abs(p->GetPDGCode()) != 211) continue;
            gM->SetPoint(miN++, p->GetZ() - offZ, p->GetX() - offX);
        }
        if (miN > 0) {
            gM->SetMarkerStyle(29); gM->SetMarkerSize(1.5); gM->SetMarkerColor(kBlack);
            gM->Draw("P SAME");
        }

        auto* gV = new TGraph(); gV->SetPoint(0, vzL, vxL);
        gV->SetMarkerStyle(33); gV->SetMarkerSize(2.0); gV->SetMarkerColor(kMagenta+2);
        gV->Draw("P SAME");
        double az = vzL + 2.5*pz/mag, ax = vxL + 2.5*px/mag;
        auto* arr = new TArrow(vzL, vxL, az, ax, 0.018, "|>");
        arr->SetLineColor(kMagenta+2); arr->SetLineWidth(2); arr->Draw();

        if (std::isfinite(R_fit) && R_fit > 0) {
            double pmnF=1e9, pmxF=-1e9;
            for (auto& h : hits) {
                double phi = std::atan2(h.GetX()-cen.first, h.GetZ()-cen.second);
                pmnF = std::min(pmnF, phi); pmxF = std::max(pmxF, phi);
            }
            double dpF = pmxF-pmnF; pmnF -= 0.05*dpF; pmxF += 0.05*dpF;
            auto* gc = new TGraph();
            for (int k = 0; k <= 200; ++k) {
                double p = pmnF + (pmxF-pmnF)*k/200;
                gc->SetPoint(k, cen.second + R_fit*std::cos(p), cen.first + R_fit*std::sin(p));
            }
            gc->SetLineColor(kRed+1); gc->SetLineWidth(2);
            gc->Draw("L SAME");
        }
    }
    c->SaveAs("evtdisp_z_horizontal.png");
}
