// Stage-by-stage reconstruction efficiency vs pion momentum.
// For each p in {200..1200}:
//   eps_reco = events with >=1 reco hit / events generated (500)
//   eps_seed = events with >=1 R3BGTPCTrackData track   / generated
//   eps_UKF  = events with converged UKF fit            / generated
//
// Reads files from the canonical chain dirs (no separate CSV needed).
// Produces plots/efficiency_vs_p.png.
void plot_efficiency()
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

    TString wd = gSystem->Getenv("VMCWORKDIR");
    std::vector<int> P = {200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200};
    const int NGEN = 500;

    TGraph *gReco = new TGraph();
    TGraph *gSeed = new TGraph();
    TGraph *gUKF  = new TGraph();

    for (size_t k = 0; k < P.size(); ++k) {
        TString suf = Form("_p%d", P[k]);
        TFile fR((wd + "/glad-tpc/macros/reco/output_reco" + suf + ".root").Data(), "READ");
        TFile fT((wd + "/glad-tpc/macros/tracking/output_tracking" + suf + ".root").Data(), "READ");
        TFile fU((wd + "/glad-tpc/macros/tracking/output_ukf" + suf + ".root").Data(), "READ");
        auto* tR = (TTree*)fR.Get("evt");
        auto* tT = (TTree*)fT.Get("evt");
        auto* tU = (TTree*)fU.Get("evt");
        if (!tR || !tT || !tU) { printf("missing for p=%d\n", P[k]); continue; }

        auto* recoHits = new TClonesArray("R3BGTPCHitData");
        auto* trks     = new TClonesArray("R3BGTPCTrackData");
        auto* fitted   = new TClonesArray("R3BGTPCFittedTrackData");
        tR->SetBranchAddress("GTPCHitData", &recoHits);
        tT->SetBranchAddress("GTPCTrackData", &trks);
        tU->SetBranchAddress("GTPCFittedTrackData", &fitted);

        int nReco = 0, nSeed = 0, nUKF = 0;
        Long64_t nE = std::min({tR->GetEntries(), tT->GetEntries(), tU->GetEntries()});
        for (Long64_t i = 0; i < nE; ++i) {
            tR->GetEntry(i); tT->GetEntry(i); tU->GetEntry(i);
            if (recoHits->GetEntries() >= 1) ++nReco;
            if (trks->GetEntries() >= 1) ++nSeed;
            if (fitted->GetEntries() >= 1) {
                auto* ft = (R3BGTPCFittedTrackData*)fitted->At(0);
                if (ft->IsConverged()) ++nUKF;
            }
        }
        printf("p=%d  reco=%d/%d  seed=%d/%d  UKF=%d/%d\n",
               P[k], nReco, NGEN, nSeed, NGEN, nUKF, NGEN);
        gReco->SetPoint(k, P[k], 100.0 * nReco / NGEN);
        gSeed->SetPoint(k, P[k], 100.0 * nSeed / NGEN);
        gUKF ->SetPoint(k, P[k], 100.0 * nUKF  / NGEN);
    }

    auto* c = new TCanvas("c_eff","efficiency", 1100, 650);
    c->SetLeftMargin(0.12); c->SetRightMargin(0.04);
    c->SetBottomMargin(0.13); c->SetTopMargin(0.05);

    gReco->SetMarkerStyle(20); gReco->SetMarkerColor(kAzure+2); gReco->SetLineColor(kAzure+2); gReco->SetLineWidth(3); gReco->SetMarkerSize(1.4);
    gSeed->SetMarkerStyle(21); gSeed->SetMarkerColor(kOrange+1); gSeed->SetLineColor(kOrange+1); gSeed->SetLineWidth(3); gSeed->SetMarkerSize(1.4);
    gUKF ->SetMarkerStyle(22); gUKF ->SetMarkerColor(kGreen+2);  gUKF ->SetLineColor(kGreen+2);  gUKF ->SetLineWidth(3); gUKF ->SetMarkerSize(1.6);

    auto* mg = new TMultiGraph();
    mg->Add(gReco, "LP"); mg->Add(gSeed, "LP"); mg->Add(gUKF, "LP");
    mg->Draw("A");
    mg->GetXaxis()->SetTitle("p_{MC} [MeV/c]");
    mg->GetYaxis()->SetTitle("efficiency [%]");
    mg->GetXaxis()->SetTitleSize(0.05);
    mg->GetYaxis()->SetTitleSize(0.05);
    mg->GetXaxis()->SetLabelSize(0.045);
    mg->GetYaxis()->SetLabelSize(0.045);
    mg->GetYaxis()->SetRangeUser(0, 110);
    mg->GetXaxis()->SetLimits(150, 1300);

    // Reference horizontal at 100 %
    auto* l100 = new TLine(150, 100, 1300, 100);
    l100->SetLineStyle(2); l100->SetLineColor(kGray+2); l100->Draw();

    auto* leg = new TLegend(0.55, 0.18, 0.94, 0.42);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.045);
    leg->AddEntry(gReco, "reco hits (#geq 1)", "lp");
    leg->AddEntry(gSeed, "seed track (TripletClust + Pratt)", "lp");
    leg->AddEntry(gUKF,  "UKF converged",     "lp");
    leg->Draw();

    TString outDir = wd + "/glad-tpc/macros/tracking/ukf_resolution_2026-05-17/plots/";
    c->SaveAs(outDir + "efficiency_vs_p.png");
    printf("Wrote %sefficiency_vs_p.png\n", outDir.Data());
}
