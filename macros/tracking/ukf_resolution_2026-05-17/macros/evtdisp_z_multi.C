void evtdisp_z_multi()
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
    for (Long64_t i = 0; i < tR->GetEntries(); ++i) {
        tT->GetEntry(i);
        if (trks->GetEntries() < 2) continue;
        evts.push_back(i);
    }
    printf("found %lu multi-cluster events\n", evts.size());
    if ((int)evts.size() > 9) evts.resize(9);

    int cols=3, rows=(evts.size()+cols-1)/cols;
    TCanvas* c = new TCanvas("c","", 1800, 400*rows);
    c->Divide(cols, rows, 0.006, 0.025);

    const int clusterColors[] = {kBlue+1, kRed+1, kGreen+2, kOrange+7, kViolet+1, kCyan+2};

    int ip = 0;
    for (auto i : evts) {
        ++ip; c->cd(ip);
        gPad->SetLeftMargin(0.10);
        gPad->SetRightMargin(0.03);
        gPad->SetTopMargin(0.10);
        gPad->SetBottomMargin(0.14);
        tR->GetEntry(i); tT->GetEntry(i); tS->GetEntry(i);

        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mc->At(j);
            if (m->GetMotherId()==-1 && std::abs(m->GetPdgCode())==211) { pi = m; break; }
        }
        double vxL = pi ? pi->GetStartX() - offX : -999;
        double vzL = pi ? pi->GetStartZ() - offZ : -999;

        // Frame includes vertex + all cluster hits + pad outline
        double zlo = std::min(padZmin, vzL), zhi = std::max(padZmax, vzL);
        double xlo = std::min(padXmin, vxL), xhi = std::max(padXmax, vxL);
        const int nTrk = trks->GetEntries();
        for (int k = 0; k < nTrk; ++k) {
            auto* tr = (R3BGTPCTrackData*)trks->At(k);
            for (auto& h : tr->GetHitArray()) {
                zlo = std::min(zlo,(double)h.GetZ()); zhi = std::max(zhi,(double)h.GetZ());
                xlo = std::min(xlo,(double)h.GetX()); xhi = std::max(xhi,(double)h.GetX());
            }
        }
        double margin = 1.0;
        auto* frame = gPad->DrawFrame(zlo-margin, xlo-margin, zhi+margin, xhi+margin);
        // total hits and cluster sizes
        int totalH = 0;
        TString sizes;
        for (int k = 0; k < nTrk; ++k) {
            auto* tr = (R3BGTPCTrackData*)trks->At(k);
            int n = tr->GetHitArray().size();
            totalH += n;
            sizes += TString::Format("%d", n);
            if (k < nTrk-1) sizes += "+";
        }
        frame->SetTitle(TString::Format("evt %lld   %d clusters   hits: %s = %d",
                                        (long long)i, nTrk, sizes.Data(), totalH));
        frame->GetXaxis()->SetTitle("z_{local} (beam) [cm]");
        frame->GetYaxis()->SetTitle("x_{local} [cm]");

        // Pad outline
        auto* box = new TBox(padZmin, padXmin, padZmax, padXmax);
        box->SetLineColor(kGray+2); box->SetLineWidth(2); box->SetLineStyle(7); box->SetFillStyle(0);
        box->Draw("L SAME");

        // Each cluster in its own colour
        for (int k = 0; k < nTrk; ++k) {
            auto* tr = (R3BGTPCTrackData*)trks->At(k);
            auto& hits = tr->GetHitArray();
            auto* gH = new TGraph(); int gi = 0;
            for (auto& h : hits) gH->SetPoint(gi++, h.GetZ(), h.GetX());
            gH->SetMarkerStyle(20); gH->SetMarkerSize(0.55);
            gH->SetMarkerColor(clusterColors[k % 6]);
            gH->Draw("P SAME");

            // fitted circle in same colour
            double R_fit = tr->GetGeoRadius();
            auto cen = tr->GetGeoCenter();
            if (std::isfinite(R_fit) && R_fit > 0 && hits.size() >= 3) {
                double pmnF=1e9, pmxF=-1e9;
                for (auto& h : hits) {
                    double phi = std::atan2(h.GetX()-cen.first, h.GetZ()-cen.second);
                    pmnF = std::min(pmnF, phi); pmxF = std::max(pmxF, phi);
                }
                double dpF = pmxF-pmnF; pmnF -= 0.05*dpF; pmxF += 0.05*dpF;
                auto* gc = new TGraph();
                for (int kk = 0; kk <= 200; ++kk) {
                    double p = pmnF + (pmxF-pmnF)*kk/200;
                    gc->SetPoint(kk, cen.second + R_fit*std::cos(p), cen.first + R_fit*std::sin(p));
                }
                gc->SetLineColor(clusterColors[k % 6]); gc->SetLineWidth(2);
                gc->Draw("L SAME");
            }
        }

        // MC truth points (black stars)
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

        // Vertex + arrow
        if (pi) {
            auto* gV = new TGraph(); gV->SetPoint(0, vzL, vxL);
            gV->SetMarkerStyle(33); gV->SetMarkerSize(2.0); gV->SetMarkerColor(kMagenta+2);
            gV->Draw("P SAME");
            double px = pi->GetPx(), pz = pi->GetPz();
            double mag = std::hypot(px, pz);
            if (mag > 1e-6) {
                double az = vzL + 2.5*pz/mag, ax = vxL + 2.5*px/mag;
                auto* arr = new TArrow(vzL, vxL, az, ax, 0.018, "|>");
                arr->SetLineColor(kMagenta+2); arr->SetLineWidth(2); arr->Draw();
            }
        }
    }
    c->SaveAs("evtdisp_z_multi.png");
}
