/// @file draw_padplane.C
/// @brief Digitized track display on the GLAD-TPC pad plane.
///
/// Three panels per event:
///   - Pad plane (z, x) colored by integrated ADC
///   - Side view (pad column z) × time bucket — drift pattern in z
///   - Top  view (pad row    x) × time bucket — drift pattern in x

void draw_padplane(int p_MeV = 600, int nEvents = 4)
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);

    TString workDir = gSystem->Getenv("VMCWORKDIR");
    TFile fLang(TString::Format("%s/glad-tpc/macros/proj/Prototype/lang_p%d.root", workDir.Data(), p_MeV));
    auto* t = (TTree*)fLang.Get("evt");
    if (!t) {
        std::cout << "[draw_padplane] no 'evt' tree at " << fLang.GetName() << "\n";
        return;
    }
    auto* cal = new TClonesArray("R3BGTPCCalData");
    t->SetBranchAddress("GTPCCalData", &cal);

    // Pad plane: 128 columns (z) × 44 rows (x), 2 mm pitch, no offset.
    const int nColZ = 128;
    const int nRowX = 44;
    const double padSize_mm = 2.0;
    const double padPlaneZ_cm = nColZ * padSize_mm / 10.0; // 25.6
    const double padPlaneX_cm = nRowX * padSize_mm / 10.0; // 8.8

    // Pick events with a reasonable number of hit pads (skip empties / mega-deltas)
    std::vector<int> picks;
    for (Long64_t i = 0; i < t->GetEntries() && (int)picks.size() < nEvents; ++i) {
        t->GetEntry(i);
        const int nPads = cal->GetEntries();
        if (nPads < 20 || nPads > 400)
            continue;
        picks.push_back(i);
    }
    if (picks.empty()) {
        std::cout << "[draw_padplane] no suitable events\n";
        return;
    }

    auto* c = new TCanvas("cpp", Form("Pad-plane p=%d MeV/c", p_MeV), 1700, 320 * picks.size());
    c->Divide(3, picks.size());

    for (size_t k = 0; k < picks.size(); ++k) {
        int i = picks[k];
        t->GetEntry(i);

        auto* hPad = new TH2F(Form("hPad_%zu", k),
                              Form("evt %d  N_pads=%d ;z (cm);x (cm)", i, cal->GetEntries()),
                              nColZ, 0., padPlaneZ_cm,
                              nRowX, 0., padPlaneX_cm);
        // Side: z (col) vs time bucket, charge weighted
        const int nTB = 512; // AGET 512 TBs typical; we'll auto-detect length
        auto* hZT = new TH2F(Form("hZT_%zu", k),
                             Form("evt %d  z (col) vs time bucket;z (cm);TB", i),
                             nColZ, 0., padPlaneZ_cm, nTB, 0., nTB);
        auto* hXT = new TH2F(Form("hXT_%zu", k),
                             Form("evt %d  x (row) vs time bucket;x (cm);TB", i),
                             nRowX, 0., padPlaneX_cm, nTB, 0., nTB);

        for (int j = 0; j < cal->GetEntries(); ++j) {
            auto* d = (R3BGTPCCalData*)cal->At(j);
            const int padId = d->GetPadId();
            const int icol = padId / nRowX;
            const int irow = padId % nRowX;
            const double z_cm = (icol + 0.5) * padSize_mm / 10.0;
            const double x_cm = (irow + 0.5) * padSize_mm / 10.0;
            const auto& adc = d->GetADC();
            double padQ = 0.;
            for (size_t b = 0; b < adc.size(); ++b) {
                padQ += adc[b];
                hZT->Fill(z_cm, b, adc[b]);
                hXT->Fill(x_cm, b, adc[b]);
            }
            hPad->Fill(z_cm, x_cm, padQ);
        }

        c->cd(3 * k + 1);
        gPad->SetLogz(0);
        hPad->Draw("COLZ");

        c->cd(3 * k + 2);
        hZT->Draw("COLZ");

        c->cd(3 * k + 3);
        hXT->Draw("COLZ");
    }

    c->SaveAs(Form("padplane_p%d.png", p_MeV));
    std::cout << "Wrote padplane_p" << p_MeV << ".png\n";
}
