// Overlay σ_p/p sweep with and without STEMAX=0.1 cm.
// Reads scan_p_results.csv (default, STEMAX on) and
// scan_p_results_nostemax.csv (historical, no step limit).
void plot_stemax_compare()
{
    gStyle->SetOptStat(0);
    TString here = TString(gSystem->Getenv("VMCWORKDIR")) +
                   "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17";

    auto load = [&](const char* fn, std::vector<double>& p,
                    std::vector<double>& sR, std::vector<double>& sUKF,
                    std::vector<double>& bUKF) {
        std::ifstream f((here + "/" + fn).Data());
        std::string line; std::getline(f, line);
        while (std::getline(f, line)) {
            double pm, N, sr, sps, spu, br, bps, bpu;
            if (sscanf(line.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                       &pm, &N, &sr, &sps, &spu, &br, &bps, &bpu) == 8) {
                if (N >= 100 && std::isfinite(sr) && std::isfinite(spu)) {
                    p.push_back(pm); sR.push_back(sr); sUKF.push_back(spu); bUKF.push_back(bpu);
                }
            }
        }
    };

    std::vector<double> p1, sR1, sU1, b1, p2, sR2, sU2, b2;
    load("scan_p_results.csv",          p1, sR1, sU1, b1);  // STEMAX=0.1
    load("scan_p_results_nostemax.csv", p2, sR2, sU2, b2);  // no limit

    auto mkG = [](std::vector<double>& x, std::vector<double>& y, int col, int mk) {
        auto* g = new TGraph(x.size(), x.data(), y.data());
        g->SetMarkerStyle(mk); g->SetMarkerColor(col); g->SetMarkerSize(1.4);
        g->SetLineColor(col); g->SetLineStyle(2); g->SetLineWidth(2);
        return g;
    };
    auto* gR1 = mkG(p1, sR1, kBlue+1,   20);
    auto* gR2 = mkG(p2, sR2, kBlue+1,   24);
    auto* gU1 = mkG(p1, sU1, kRed+1,    21);
    auto* gU2 = mkG(p2, sU2, kRed+1,    25);
    auto* gB1 = mkG(p1, b1,  kBlack,    20);
    auto* gB2 = mkG(p2, b2,  kBlack,    24);

    TCanvas c("c","STEMAX comparison", 1400, 550);
    c.Divide(2, 1, 0.008, 0.008);

    c.cd(1);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    auto* h1 = new TH1F("h1", ";p_{MC} [MeV/c]; #sigma [%]", 14, 200, 1300);
    h1->SetMaximum(7); h1->SetMinimum(0); h1->Draw();
    gR2->Draw("PL SAME"); gR1->Draw("PL SAME");
    gU2->Draw("PL SAME"); gU1->Draw("PL SAME");
    auto* l4 = new TLine(200, 4, 1300, 4); l4->SetLineStyle(2); l4->SetLineColor(kGray+2); l4->Draw();
    TLatex t; t.SetTextSize(0.030); t.SetTextColor(kGray+2);
    t.DrawLatex(1050, 4.15, "4 % target");
    auto* lg = new TLegend(0.45, 0.62, 0.92, 0.88);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.030);
    lg->AddEntry(gR1, "#sigma_{R}/R  seed (STEMAX=0.1 cm)", "p");
    lg->AddEntry(gR2, "#sigma_{R}/R  seed (no STEMAX)",      "p");
    lg->AddEntry(gU1, "#sigma_{p}/p  UKF  (STEMAX=0.1 cm)",  "p");
    lg->AddEntry(gU2, "#sigma_{p}/p  UKF  (no STEMAX)",      "p");
    lg->Draw();

    c.cd(2);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    auto* h2 = new TH1F("h2", ";p_{MC} [MeV/c]; bias #sigma_{p}/p UKF [%]", 14, 200, 1300);
    h2->SetMaximum(3); h2->SetMinimum(-1); h2->Draw();
    auto* l0 = new TLine(200, 0, 1300, 0); l0->SetLineStyle(1); l0->SetLineColor(kGray+1); l0->Draw();
    gB2->Draw("PL SAME"); gB1->Draw("PL SAME");
    auto* lg2 = new TLegend(0.45, 0.74, 0.92, 0.88);
    lg2->SetBorderSize(0); lg2->SetFillStyle(0); lg2->SetTextSize(0.030);
    lg2->AddEntry(gB1, "bias UKF (STEMAX=0.1 cm)", "p");
    lg2->AddEntry(gB2, "bias UKF (no STEMAX)",     "p");
    lg2->Draw();

    c.SaveAs(here + "/plots/sigma_vs_p_stemax_compare.png");
}
