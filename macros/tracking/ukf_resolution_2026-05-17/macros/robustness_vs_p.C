// Per-momentum robustness of the three sigma estimators on the box-gen
// scan files. Tests the hypothesis that the apparently flat sigma_p/p
// plateau across 400-1200 MeV/c is a Gauss-core-window artefact: the
// core sigma stays flat because [-0.4, 0.4] only catches the peak,
// while the wings widen with p as Gluckstern predicts.
//
// For each p, fills the p_UKF/p_MC - 1 distribution and computes:
//   sigma_core   -- Gaussian fit to the core in [-0.4, 0.4]
//   sigma_trunc  -- RMS within +/- 3*sigma_core of the Gaussian mean
//   sigma_q      -- (q84 - q16) / 2  (quantile half-width)
//
// Writes:
//   plots/robustness_vs_p.png       -- sigma vs p (3 estimators)
//   plots/robustness_vs_p.csv       -- per-p table
//
// Reads files from glad-tpc/macros/{sim/Prototype,tracking}/.
// Default PLIST matches scan_p.sh.

#include <vector>

void robustness_vs_p(TString plist = "200 300 400 500 600 700 800 900 1000 1200")
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

    const double m_pi = 139.57039;
    auto p_from_KE = [&](double ke){ return std::sqrt((ke+m_pi)*(ke+m_pi) - m_pi*m_pi); };

    TString wd = gSystem->Getenv("VMCWORKDIR");
    TString simDir = wd + "/glad-tpc/macros/sim/Prototype/";
    TString trkDir = wd + "/glad-tpc/macros/tracking/";

    std::vector<int>    pv;
    std::vector<double> sCore, sTrunc, sQ;
    std::vector<double> esCore, esTrunc, esQ;
    std::vector<int>    nv;

    // Tokenize PLIST
    TObjArray* toks = plist.Tokenize(" ");
    printf("\n%-6s | %-6s | %-8s | %-8s | %-8s | %-6s\n",
           "p", "N", "core", "trunc", "quant", "ratio");
    printf("%s\n", TString('-', 60).Data());

    for (int it = 0; it < toks->GetEntries(); ++it) {
        int p = ((TObjString*)toks->At(it))->GetString().Atoi();
        TString simFile = TString::Format("%ssim_p%d.root", simDir.Data(), p);
        TString ukfFile = TString::Format("%soutput_ukf_p%d.root", trkDir.Data(), p);

        if (gSystem->AccessPathName(simFile)) { printf("MISSING %s\n", simFile.Data()); continue; }
        if (gSystem->AccessPathName(ukfFile)) { printf("MISSING %s\n", ukfFile.Data()); continue; }

        TFile fS(simFile);
        TFile fU(ukfFile);
        auto* tS = (TTree*)fS.Get("evt");
        auto* tU = (TTree*)fU.Get("evt");
        if (!tS || !tU) { printf("p=%d: bad trees\n", p); continue; }

        auto* mc     = new TClonesArray("R3BMCTrack");
        auto* fitted = new TClonesArray("R3BGTPCFittedTrackData");
        tS->SetBranchAddress("MCTrack", &mc);
        tU->SetBranchAddress("GTPCFittedTrackData", &fitted);

        // Wider binning at low p (wide wings); fixed [-1, 1] window is
        // enough -- residuals beyond +/- 1 are pathological.
        auto* h = new TH1F(Form("h_p%d", p),
                           ";p_{UKF}/p_{MC} - 1;events / 0.025",
                           80, -1.0, 1.0);

        Long64_t nE = std::min(tS->GetEntries(), tU->GetEntries());
        int nFit = 0;
        for (Long64_t i = 0; i < nE; ++i) {
            tS->GetEntry(i); tU->GetEntry(i);
            if (fitted->GetEntries() == 0) continue;
            auto* ft = (R3BGTPCFittedTrackData*)fitted->At(0);
            if (!ft->IsConverged()) continue;
            const auto& kin = ft->GetKinematicsXtr();
            if (!std::isfinite(kin.kineticEnergy) || kin.kineticEnergy <= 0) continue;
            R3BMCTrack* pi = nullptr;
            for (int j = 0; j < mc->GetEntries(); ++j) {
                auto* mt = (R3BMCTrack*)mc->At(j);
                if (mt->GetMotherId() == -1 && mt->GetPdgCode() == -211) { pi = mt; break; }
            }
            if (!pi) continue;
            double pMC = std::sqrt(pi->GetPx()*pi->GetPx() + pi->GetPy()*pi->GetPy() + pi->GetPz()*pi->GetPz()) * 1000;
            if (pMC <= 0) continue;
            double pFit = p_from_KE(kin.kineticEnergy);
            h->Fill(pFit/pMC - 1);
            ++nFit;
        }

        if (nFit < 20) {
            printf("p=%d: only %d fits, skipping\n", p, nFit);
            delete h; delete mc; delete fitted;
            continue;
        }

        // 1. Gaussian core
        h->Fit("gaus", "Q0", "", -0.4, 0.4);
        auto* fG = h->GetFunction("gaus");
        double mu_g = fG->GetParameter(1);
        double s_g  = fG->GetParameter(2);
        double es_g = fG->GetParError(2);

        // 2. Truncated RMS within +/- 3 sigma_core of the Gaussian mean
        int bLo = h->GetXaxis()->FindBin(mu_g - 3*s_g);
        int bHi = h->GetXaxis()->FindBin(mu_g + 3*s_g);
        h->GetXaxis()->SetRange(bLo, bHi);
        double s_trunc = h->GetRMS();
        double es_trunc = h->GetRMSError();
        h->GetXaxis()->SetRange(0, 0);

        // 3. Quantile half-width
        double probs[3] = { 0.16, 0.50, 0.84 };
        double q[3];
        h->GetQuantiles(3, q, probs);
        double s_q = 0.5 * (q[2] - q[0]);
        // Quantile uncertainty (rough): sigma_q ~ sqrt(p(1-p)/N) / pdf(q).
        // Approximate with sigma_q / sqrt(N) which is the usual scale.
        double es_q = s_q / std::sqrt((double)nFit);

        printf("%-6d | %-6d | %6.2f %% | %6.2f %% | %6.2f %% | %5.2f\n",
               p, nFit, 100*s_g, 100*s_trunc, 100*s_q, s_q/s_g);

        pv.push_back(p);
        nv.push_back(nFit);
        sCore.push_back(s_g);    esCore.push_back(es_g);
        sTrunc.push_back(s_trunc); esTrunc.push_back(es_trunc);
        sQ.push_back(s_q);       esQ.push_back(es_q);

        delete h; delete mc; delete fitted;
    }
    delete toks;

    if (pv.empty()) { printf("No data\n"); return; }

    // ---- Plot ----
    int np = pv.size();
    std::vector<double> px(np), epx(np, 0);
    std::vector<double> yC(np), yT(np), yQ(np);
    std::vector<double> eyC(np), eyT(np), eyQ(np);
    for (int i = 0; i < np; ++i) {
        px[i] = pv[i];
        yC[i]  = 100*sCore[i];  eyC[i]  = 100*esCore[i];
        yT[i]  = 100*sTrunc[i]; eyT[i]  = 100*esTrunc[i];
        yQ[i]  = 100*sQ[i];     eyQ[i]  = 100*esQ[i];
    }
    auto gC = new TGraphErrors(np, px.data(), yC.data(), epx.data(), eyC.data());
    auto gT = new TGraphErrors(np, px.data(), yT.data(), epx.data(), eyT.data());
    auto gQ = new TGraphErrors(np, px.data(), yQ.data(), epx.data(), eyQ.data());
    gC->SetMarkerStyle(20); gC->SetMarkerColor(kRed);    gC->SetLineColor(kRed);    gC->SetMarkerSize(1.2);
    gT->SetMarkerStyle(21); gT->SetMarkerColor(kBlue+1); gT->SetLineColor(kBlue+1); gT->SetMarkerSize(1.1);
    gQ->SetMarkerStyle(22); gQ->SetMarkerColor(kGreen+2);gQ->SetLineColor(kGreen+2);gQ->SetMarkerSize(1.3);

    auto c = new TCanvas("c_rob_vs_p", "robustness vs p", 1100, 700);
    c->SetLeftMargin(0.11); c->SetRightMargin(0.04);
    c->SetBottomMargin(0.13); c->SetTopMargin(0.06);

    double ymax = 0;
    for (double v : yC) ymax = std::max(ymax, v);
    for (double v : yT) ymax = std::max(ymax, v);
    for (double v : yQ) ymax = std::max(ymax, v);
    auto* fr = c->DrawFrame(pv.front()-50, 0, pv.back()+50, 1.4*ymax);
    fr->GetXaxis()->SetTitle("p_{MC} [MeV/c]");
    fr->GetYaxis()->SetTitle("#sigma (p_{UKF}/p_{MC} #minus 1) [%]");
    fr->GetXaxis()->SetTitleSize(0.05);
    fr->GetYaxis()->SetTitleSize(0.05);
    fr->GetXaxis()->SetLabelSize(0.045);
    fr->GetYaxis()->SetLabelSize(0.045);

    // Reference: 4% target line
    auto* ln4 = new TLine(pv.front()-50, 4.0, pv.back()+50, 4.0);
    ln4->SetLineStyle(2); ln4->SetLineColor(kGray+1); ln4->Draw();

    // Gluckstern-shape reference: scale Gauss-core at lowest p
    // (skip if p_low has too few stats)
    int iref = -1;
    for (int i = 0; i < np; ++i) if (nv[i] >= 100) { iref = i; break; }
    if (iref >= 0) {
        auto* gG = new TGraph(np);
        for (int i = 0; i < np; ++i) gG->SetPoint(i, px[i], yC[iref] * px[i] / px[iref]);
        gG->SetLineColor(kGray+2); gG->SetLineStyle(7); gG->SetLineWidth(2);
        gG->Draw("L SAME");
    }

    gC->Draw("PL SAME");
    gT->Draw("PL SAME");
    gQ->Draw("PL SAME");

    auto* lg = new TLegend(0.50, 0.62, 0.95, 0.92);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.038);
    lg->AddEntry(gC, "Gaussian core  [-0.4, 0.4]", "lp");
    lg->AddEntry(gT, "Truncated RMS  (#pm 3#sigma_{core})", "lp");
    lg->AddEntry(gQ, "Quantile (q_{84} - q_{16}) / 2", "lp");
    if (iref >= 0)
        lg->AddEntry((TObject*)nullptr,
            Form("Gluckstern #propto p (anchored at %d MeV/c)", pv[iref]), "");
    lg->AddEntry((TObject*)nullptr, "ATTPCROOT 4 % target", "");
    lg->Draw();

    TString outDir = wd + "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17/plots/";
    c->SaveAs(outDir + "robustness_vs_p.png");
    printf("Wrote %srobustness_vs_p.png\n", outDir.Data());

    // CSV
    FILE* fcsv = fopen((outDir + "robustness_vs_p.csv").Data(), "w");
    fprintf(fcsv, "p_MeV,N,sigma_core_pct,sigma_trunc_pct,sigma_q_pct,ratio_q_over_core\n");
    for (int i = 0; i < np; ++i)
        fprintf(fcsv, "%d,%d,%.3f,%.3f,%.3f,%.3f\n",
                pv[i], nv[i], 100*sCore[i], 100*sTrunc[i], 100*sQ[i], sQ[i]/sCore[i]);
    fclose(fcsv);
    printf("Wrote %srobustness_vs_p.csv\n", outDir.Data());
}
