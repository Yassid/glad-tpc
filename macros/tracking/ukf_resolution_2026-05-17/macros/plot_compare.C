// Side-by-side σ_p/p vs p plot: box-gen sweep and good_evt p-binned
// analysis on the same axes.
void plot_compare()
{
    gStyle->SetOptStat(0);
    TString here = TString(gSystem->Getenv("VMCWORKDIR")) +
                   "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17";

    auto loadCsv = [](const TString& path,
                      std::vector<double>& p,
                      std::vector<double>& sR,
                      std::vector<double>& sPs,
                      std::vector<double>& sPu) {
        std::ifstream f(path.Data());
        if (!f.is_open()) { std::cerr << "missing " << path << "\n"; return; }
        std::string line; std::getline(f, line); // header
        while (std::getline(f, line)) {
            double pV, N, srV, spsV, spuV, br, bps, bpu;
            if (sscanf(line.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                       &pV, &N, &srV, &spsV, &spuV, &br, &bps, &bpu) == 8) {
                if (!std::isfinite(srV)) continue;
                p.push_back(pV); sR.push_back(srV);
                sPs.push_back(spsV); sPu.push_back(spuV);
            }
        }
    };

    std::vector<double> pBG, sR_BG, sPs_BG, sPu_BG;
    std::vector<double> pGE, sR_GE, sPs_GE, sPu_GE;
    loadCsv(here + "/scan_p_results.csv",       pBG, sR_BG, sPs_BG, sPu_BG);
    loadCsv(here + "/scan_goodevt_pbins.csv",   pGE, sR_GE, sPs_GE, sPu_GE);

    auto mk = [](const std::vector<double>& x, const std::vector<double>& y,
                 int color, int style, double sz) {
        auto* g = new TGraph();
        int n = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            if (std::isfinite(y[i])) { g->SetPoint(n++, x[i], std::min(15.0, y[i])); }
        }
        g->SetMarkerColor(color); g->SetMarkerStyle(style); g->SetMarkerSize(sz);
        g->SetLineColor(color); g->SetLineStyle(2); g->SetLineWidth(1);
        return g;
    };

    TCanvas c("c", "resolution comparison", 1600, 600);
    c.Divide(2, 1, 0.005, 0.005);

    auto* hL = new TH1F("hL", ";p_{MC} [MeV/c]; #sigma_{R}/R, #sigma_{p}/p [%]", 14, 100, 1300);
    auto* hR = new TH1F("hR", ";p_{MC} [MeV/c]; #sigma_{p,tot}/p_{tot}  (UKF) [%]", 14, 100, 1300);

    c.cd(1);
    gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.13);
    hL->SetMinimum(0); hL->SetMaximum(8); hL->Draw();
    auto* gR_BG  = mk(pBG, sR_BG,  kBlue+1,  20, 1.3);
    auto* gPs_BG = mk(pBG, sPs_BG, kGreen+2, 21, 1.3);
    auto* gR_GE  = mk(pGE, sR_GE,  kBlue+1,  24, 1.2);
    auto* gPs_GE = mk(pGE, sPs_GE, kGreen+2, 25, 1.2);
    gR_BG ->Draw("PL SAME"); gPs_BG->Draw("PL SAME");
    gR_GE ->Draw("PL SAME"); gPs_GE->Draw("PL SAME");
    auto* l4 = new TLine(100, 4, 1300, 4); l4->SetLineStyle(2); l4->SetLineColor(kGray+2); l4->Draw();
    TLatex t; t.SetTextSize(0.030); t.SetTextColor(kGray+2);
    t.DrawLatex(950, 4.2, "4% target");
    auto* lg1 = new TLegend(0.35, 0.6, 0.88, 0.88);
    lg1->SetBorderSize(0); lg1->SetFillStyle(0); lg1->SetTextSize(0.030);
    lg1->AddEntry(gR_BG,  "#sigma_{R}/R  box-gen (single-#pi^{-})",       "p");
    lg1->AddEntry(gPs_BG, "#sigma_{p}/p  box-gen (seed via #theta_{y})",  "p");
    lg1->AddEntry(gR_GE,  "#sigma_{R}/R  good_evt (binned)",              "p");
    lg1->AddEntry(gPs_GE, "#sigma_{p}/p  good_evt (binned)",              "p");
    lg1->Draw();

    c.cd(2);
    gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.13);
    hR->SetMinimum(0); hR->SetMaximum(15); hR->Draw();
    auto* gPu_BG = mk(pBG, sPu_BG, kRed+1, 22, 1.4);
    auto* gPu_GE = mk(pGE, sPu_GE, kRed+1, 26, 1.3);
    gPu_BG->Draw("PL SAME"); gPu_GE->Draw("PL SAME");
    auto* l4b = new TLine(100, 4, 1300, 4); l4b->SetLineStyle(2); l4b->SetLineColor(kGray+2); l4b->Draw();
    auto* lg2 = new TLegend(0.40, 0.7, 0.88, 0.88);
    lg2->SetBorderSize(0); lg2->SetFillStyle(0); lg2->SetTextSize(0.030);
    lg2->AddEntry(gPu_BG, "UKF  box-gen (single-#pi^{-})",   "p");
    lg2->AddEntry(gPu_GE, "UKF  good_evt (binned)",          "p");
    lg2->Draw();

    c.SaveAs(here + "/plots/sigma_vs_p_compare.png");
}
