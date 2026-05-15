/// Plot σ_p/p vs p for the R3B HYDRA scan.

void plot_scan(const char* csv = "scan_p_results.csv")
{
    auto* g = new TGraphErrors();
    auto* gb = new TGraphErrors();
    int n = 0;
    std::ifstream f(csv);
    std::string line;
    std::getline(f, line); // header
    while (std::getline(f, line)) {
        double p, nfit, nthr, sig, bias;
        if (std::sscanf(line.c_str(), "%lf,%lf,%lf,%lf,%lf", &p, &nfit, &nthr, &sig, &bias) == 5) {
            if (sig <= 0)
                continue; // skip points with no fit
            g->SetPoint(n, p, sig * 100);
            g->SetPointError(n, 0, sig / std::sqrt(2.0 * std::max(1.0, nfit)) * 100);
            gb->SetPoint(n, p, bias * 100);
            gb->SetPointError(n, 0, sig / std::sqrt(std::max(1.0, nfit)) * 100);
            ++n;
        }
    }

    auto* c = new TCanvas("c", "σ_p/p scan", 800, 600);
    c->Divide(1, 2);
    c->cd(1);
    g->SetTitle("R3B GLAD-TPC Prototype, truth-seeded UKF;p [MeV/c];#sigma_{p}/p [%]");
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue + 1);
    g->SetLineColor(kBlue + 1);
    g->Draw("AP");
    auto* line4 = new TLine(g->GetX()[0] - 50, 4, g->GetX()[g->GetN() - 1] + 50, 4);
    line4->SetLineStyle(2);
    line4->SetLineColor(kRed);
    line4->Draw();

    c->cd(2);
    gb->SetTitle(";p [MeV/c];bias [%]");
    gb->SetMarkerStyle(20);
    gb->SetMarkerColor(kBlack);
    gb->Draw("AP");

    c->SaveAs("scan_p_results.png");
    std::cout << "Wrote scan_p_results.png\n";
}
