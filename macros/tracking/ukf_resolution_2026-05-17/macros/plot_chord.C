// Plot σ_R/R, σ_p/p (seed and UKF) vs chord length from
// scan_goodevt_chordbins.csv.
void plot_chord()
{
    gStyle->SetOptStat(0);
    TString here = TString(gSystem->Getenv("VMCWORKDIR")) +
                   "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17";
    std::ifstream f((here + "/scan_goodevt_chordbins.csv").Data());
    if (!f.is_open()) { std::cerr << "no csv\n"; return; }
    std::string line; std::getline(f, line);
    auto* gN  = new TGraph();
    auto* gP  = new TGraph();
    auto* gSR = new TGraph(); auto* gSPs = new TGraph(); auto* gSPu = new TGraph();
    int n = 0;
    while (std::getline(f, line)) {
        double c, N, pAvg, sR, sPs, sPu, br, bps, bpu;
        if (sscanf(line.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                   &c, &N, &pAvg, &sR, &sPs, &sPu, &br, &bps, &bpu) == 9) {
            gN->SetPoint(n, c, N);
            if (N >= 10) gP->SetPoint(gP->GetN(), c, pAvg);
            if (std::isfinite(sR) && N >= 10) {
                gSR ->SetPoint(gSR->GetN(),  c, std::min(15.0, sR));
                gSPs->SetPoint(gSPs->GetN(), c, std::min(15.0, sPs));
                if (std::isfinite(sPu))
                    gSPu->SetPoint(gSPu->GetN(), c, std::min(15.0, sPu));
            }
            ++n;
        }
    }
    TCanvas c("c", "σ vs chord", 1800, 500);
    c.Divide(3, 1, 0.005, 0.005);

    c.cd(1);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    auto* hL = new TH1F("hL", ";chord (x,z) [cm]; resolution [%]", 11, 0, 22);
    hL->SetMinimum(0); hL->SetMaximum(15); hL->Draw();
    gSR ->SetMarkerStyle(20); gSR ->SetMarkerColor(kBlue+1);  gSR ->SetMarkerSize(1.5); gSR ->SetLineStyle(2); gSR ->SetLineColor(kBlue+1);
    gSPs->SetMarkerStyle(21); gSPs->SetMarkerColor(kGreen+2); gSPs->SetMarkerSize(1.4); gSPs->SetLineStyle(2); gSPs->SetLineColor(kGreen+2);
    gSPu->SetMarkerStyle(22); gSPu->SetMarkerColor(kRed+1);   gSPu->SetMarkerSize(1.5); gSPu->SetLineStyle(2); gSPu->SetLineColor(kRed+1);
    gSR->Draw("PL SAME"); gSPs->Draw("PL SAME"); gSPu->Draw("PL SAME");
    auto* l4 = new TLine(0, 4, 22, 4); l4->SetLineStyle(2); l4->SetLineColor(kGray+2); l4->Draw();
    TLatex t; t.SetTextSize(0.030); t.SetTextColor(kGray+2);
    t.DrawLatex(15, 4.3, "4% target");
    auto* lg = new TLegend(0.42, 0.65, 0.92, 0.88);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
    lg->AddEntry(gSR,  "#sigma_{R}/R  (seed Pratt+GN+vertex)", "p");
    lg->AddEntry(gSPs, "#sigma_{p}/p  (seed via #theta_{y})",  "p");
    lg->AddEntry(gSPu, "#sigma_{p_{tot}}/p_{tot}  (UKF)",      "p");
    lg->Draw();

    c.cd(2);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    auto* hR = new TH1F("hR", ";chord (x,z) [cm]; N events", 11, 0, 22);
    int Nmax = 0;
    for (int i = 0; i < gN->GetN(); ++i) {
        double x, y; gN->GetPoint(i, x, y);
        if (y > Nmax) Nmax = y;
    }
    hR->SetMaximum(Nmax * 1.2); hR->SetMinimum(0); hR->Draw();
    gN->SetMarkerStyle(20); gN->SetMarkerColor(kBlack); gN->SetMarkerSize(1.4);
    gN->SetLineColor(kBlack); gN->SetLineStyle(2);
    gN->Draw("PL SAME");

    c.cd(3);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    auto* h3 = new TH1F("h3", ";chord (x,z) [cm]; <p_{MC}> per bin [MeV/c]", 11, 0, 22);
    h3->SetMinimum(0); h3->SetMaximum(800); h3->Draw();
    gP->SetMarkerStyle(20); gP->SetMarkerColor(kMagenta+2); gP->SetMarkerSize(1.4);
    gP->SetLineColor(kMagenta+2); gP->SetLineStyle(2);
    gP->Draw("PL SAME");
    // Mean of <p_MC> as horizontal reference line — confirms chord ≠ p proxy.
    double pAvg = 0; int np = 0;
    for (int i = 0; i < gP->GetN(); ++i) { double x, y; gP->GetPoint(i, x, y); pAvg += y; ++np; }
    if (np > 0) pAvg /= np;
    auto* lAvg = new TLine(0, pAvg, 22, pAvg);
    lAvg->SetLineStyle(2); lAvg->SetLineColor(kGray+2); lAvg->Draw();
    TLatex tx; tx.SetTextSize(0.030); tx.SetTextColor(kGray+2);
    tx.DrawLatex(13, pAvg + 30, Form("global <p> = %.0f MeV/c", pAvg));

    c.SaveAs(here + "/plots/sigma_vs_chord.png");
}
