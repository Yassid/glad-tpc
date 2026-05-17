// Plot σ_p/p (seed + UKF) and biases vs pion momentum from scan_p_results.csv
void plot_scan_p()
{
    gStyle->SetOptStat(0);
    TString here = gSystem->Getenv("PWD");
    // The script chdirs into TRACKDIR before calling us; jump back to the
    // resolution folder to find the CSV.
    TString csv = TString(gSystem->Getenv("VMCWORKDIR")) +
                  "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17/scan_p_results.csv";
    std::ifstream f(csv.Data());
    if (!f.is_open()) { std::cerr << "Cannot open " << csv << std::endl; return; }
    std::string line;
    std::getline(f, line); // header
    auto* gSR  = new TGraph();
    auto* gSP  = new TGraph();
    auto* gSU  = new TGraph();
    auto* gBR  = new TGraph();
    auto* gBP  = new TGraph();
    auto* gBU  = new TGraph();
    int n = 0;
    while (std::getline(f, line)) {
        double p, N, sR, sPs, sPu, bR, bPs, bPu;
        if (sscanf(line.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                   &p, &N, &sR, &sPs, &sPu, &bR, &bPs, &bPu) == 8) {
            gSR->SetPoint(n, p, sR);
            gSP->SetPoint(n, p, sPs);
            gSU->SetPoint(n, p, sPu);
            gBR->SetPoint(n, p, bR);
            gBP->SetPoint(n, p, bPs);
            gBU->SetPoint(n, p, bPu);
            ++n;
        }
    }
    if (n == 0) { std::cerr << "No data in CSV" << std::endl; return; }

    TCanvas c("c", "σ_p/p scan", 1500, 600);
    c.Divide(2, 1, 0.005, 0.005);

    c.cd(1);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    auto* hL = new TH1F("hL", ";p_{MC} [MeV/c]; resolution [%]", 14, 100, 1300);
    hL->SetMinimum(0); hL->SetMaximum(15);
    hL->Draw();
    gSR ->SetMarkerStyle(20); gSR ->SetMarkerColor(kBlue+1);  gSR ->SetMarkerSize(1.4);
    gSP ->SetMarkerStyle(21); gSP ->SetMarkerColor(kGreen+2); gSP ->SetMarkerSize(1.4);
    gSU ->SetMarkerStyle(22); gSU ->SetMarkerColor(kRed+1);   gSU ->SetMarkerSize(1.5);
    gSR->Draw("P SAME"); gSP->Draw("P SAME"); gSU->Draw("P SAME");
    auto* line4 = new TLine(100, 4, 1300, 4); line4->SetLineStyle(2); line4->SetLineColor(kGray+2); line4->Draw();
    TLatex t; t.SetTextSize(0.035); t.SetTextColor(kGray+2);
    t.DrawLatex(900, 4.3, "ATTPCROOT 4% target");
    auto* lg1 = new TLegend(0.55, 0.65, 0.95, 0.88);
    lg1->SetBorderSize(0); lg1->SetFillStyle(0);
    lg1->AddEntry(gSR, "#sigma_{R}/R  (seed)", "p");
    lg1->AddEntry(gSP, "#sigma_{p}/p  (seed, via #theta_{y})", "p");
    lg1->AddEntry(gSU, "#sigma_{p_{tot}}/p_{tot}  (UKF)", "p");
    lg1->Draw();

    c.cd(2);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    auto* hR = new TH1F("hR", ";p_{MC} [MeV/c]; bias [%]", 14, 100, 1300);
    hR->SetMinimum(-15); hR->SetMaximum(15);
    hR->Draw();
    gBR ->SetMarkerStyle(20); gBR ->SetMarkerColor(kBlue+1);  gBR ->SetMarkerSize(1.4);
    gBP ->SetMarkerStyle(21); gBP ->SetMarkerColor(kGreen+2); gBP ->SetMarkerSize(1.4);
    gBU ->SetMarkerStyle(22); gBU ->SetMarkerColor(kRed+1);   gBU ->SetMarkerSize(1.5);
    gBR->Draw("P SAME"); gBP->Draw("P SAME"); gBU->Draw("P SAME");
    auto* line0 = new TLine(100, 0, 1300, 0); line0->SetLineStyle(2); line0->SetLineColor(kGray+1); line0->Draw();
    auto* lg2 = new TLegend(0.55, 0.7, 0.95, 0.88);
    lg2->SetBorderSize(0); lg2->SetFillStyle(0);
    lg2->AddEntry(gBR, "bias R (seed)", "p");
    lg2->AddEntry(gBP, "bias p (seed)", "p");
    lg2->AddEntry(gBU, "bias p (UKF)", "p");
    lg2->Draw();

    TString out = TString(gSystem->Getenv("VMCWORKDIR")) +
                  "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17/plots/sigma_vs_p.png";
    c.SaveAs(out);
}
