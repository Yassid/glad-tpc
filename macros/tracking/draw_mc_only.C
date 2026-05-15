/// @file draw_mc_only.C
/// @brief Visualize MC truth in world coordinates: vertex + primary pion
///        trajectory points + secondaries. Includes the analytical pion
///        path from the box generator (straight + helix in B=2T).

void draw_mc_only(int p_MeV = 800, int nTracks = 6)
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);

    TString workDir = gSystem->Getenv("VMCWORKDIR");
    TFile fMC(TString::Format("%s/glad-tpc/macros/sim/Prototype/sim_p%d.root", workDir.Data(), p_MeV));
    auto* tMC = (TTree*)fMC.Get("evt");
    auto* mcTracks = new TClonesArray("R3BMCTrack");
    auto* gtpcPoints = new TClonesArray("R3BGTPCPoint");
    tMC->SetBranchAddress("MCTrack", &mcTracks);
    tMC->SetBranchAddress("GTPCPoint", &gtpcPoints);

    const double mass_pi_GeV = 0.13957;
    const double B_T = 2.0; // R3B GLAD: B = (0, B_T, 0) — along +ŷ
    const double q = -1.0;  // π⁻

    // Pick the first N events
    auto* c = new TCanvas("cMC", Form("MC truth p=%d MeV/c", p_MeV), 1600, 900);
    c->Divide(3, nTracks);

    int drawn = 0;
    for (Long64_t i = 0; i < tMC->GetEntries() && drawn < nTracks; ++i) {
        tMC->GetEntry(i);
        // Find primary pion
        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mcTracks->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mcTracks->At(j);
            if (m->GetMotherId() == -1 && std::abs(m->GetPdgCode()) == 211) {
                pi = m;
                break;
            }
        }
        if (!pi)
            continue;

        // MC points: split into primary (trackID == 0) and secondaries
        auto* gPriXY = new TGraph();
        auto* gPriXZ = new TGraph();
        auto* gPriYZ = new TGraph();
        auto* gSecXY = new TGraph();
        auto* gSecXZ = new TGraph();
        auto* gSecYZ = new TGraph();
        for (auto* g : { gPriXY, gPriXZ, gPriYZ }) { g->SetMarkerStyle(20); g->SetMarkerColor(kBlue + 2); g->SetMarkerSize(1.0); }
        for (auto* g : { gSecXY, gSecXZ, gSecYZ }) { g->SetMarkerStyle(25); g->SetMarkerColor(kOrange + 1); g->SetMarkerSize(0.7); }

        double xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9, zmin = 1e9, zmax = -1e9;
        for (int j = 0; j < gtpcPoints->GetEntries(); ++j) {
            auto* p = (R3BGTPCPoint*)gtpcPoints->At(j);
            const double x = p->GetX(), y = p->GetY(), z = p->GetZ();
            auto* gXY = (p->GetTrackID() == 0) ? gPriXY : gSecXY;
            auto* gXZ = (p->GetTrackID() == 0) ? gPriXZ : gSecXZ;
            auto* gYZ = (p->GetTrackID() == 0) ? gPriYZ : gSecYZ;
            gXY->SetPoint(gXY->GetN(), x, y);
            gXZ->SetPoint(gXZ->GetN(), x, z);
            gYZ->SetPoint(gYZ->GetN(), z, y);
            xmin = std::min(xmin, x); xmax = std::max(xmax, x);
            ymin = std::min(ymin, y); ymax = std::max(ymax, y);
            zmin = std::min(zmin, z); zmax = std::max(zmax, z);
        }

        // Analytical helix starting from the FIRST primary MC point (entry
        // into the active volume). Vertex from FairBoxGenerator is upstream
        // and the simulation may apply a target-frame transform we don't
        // mirror here, so anchoring at the first MC hit is robust.
        double vx0 = 0, vy0 = 0, vz0 = 0;
        bool haveEntry = false;
        for (int j = 0; j < gtpcPoints->GetEntries(); ++j) {
            auto* p = (R3BGTPCPoint*)gtpcPoints->At(j);
            if (p->GetTrackID() != 0)
                continue;
            vx0 = p->GetX(); vy0 = p->GetY(); vz0 = p->GetZ();
            haveEntry = true;
            break;
        }
        if (!haveEntry) {
            ++drawn;
            continue;
        }
        const double px0 = pi->GetPx(), py0 = pi->GetPy(), pz0 = pi->GetPz();
        const double pt0 = std::sqrt(px0 * px0 + pz0 * pz0); // transverse to B = ŷ
        const double R_m = pt0 / (0.3 * B_T);                // m
        const double R_cm = R_m * 100.0;
        // Center in (x,z): for q<0 in +B_y the Lorentz force is along F_xz =
        // (vz, -vx)·B_y (positive multiplier), so the centre lies in the (vz,
        // -vx) direction from the particle. For q>0, opposite side.
        const double sign = (q < 0) ? +1.0 : -1.0;
        const double nhx = sign * pz0 / pt0;
        const double nhz = -sign * px0 / pt0;
        const double cx = vx0 + R_cm * nhx;
        const double cz = vz0 + R_cm * nhz;
        const double phi0 = std::atan2(vz0 - cz, vx0 - cx);
        auto* gHxXY = new TGraph();
        auto* gHxXZ = new TGraph();
        auto* gHxYZ = new TGraph();
        for (auto* g : { gHxXY, gHxXZ, gHxYZ }) { g->SetLineColor(kRed); g->SetLineWidth(2); g->SetLineStyle(2); }
        // Arclength parametrization in (x,z): dphi/ds_xz = +1/R for q<0 in +B
        // (CCW viewed from +y, with x→right, z→up). Total path: ds_xz/ds =
        // sin(θ_y) = pt0/|p|, dy/ds = py/|p|.
        const double pTot = std::sqrt(px0 * px0 + py0 * py0 + pz0 * pz0);
        const double dy_per_ds = py0 / pTot;
        const double dsxz_per_ds = pt0 / pTot;
        // ω_y = q·B_y/m → for q<0, B_y>0, rotation is CW in (x,z) viewed from +y
        // (dphi/dt < 0). The `sign` variable encodes the centre side; rotation
        // direction is its negative.
        const double dphi_per_ds = (sign > 0 ? -1.0 : +1.0) * dsxz_per_ds / R_cm;
        // Step along arc-length, but only retain segment inside the MC bbox
        // (with a small pad). Otherwise the helix from vertex (-2.7,0,237) is
        // far off-canvas.
        // Anchor: start at the entry point (s=0) and extend forward.
        for (int k = 0; k <= 4000; ++k) {
            const double s = k * 0.05;
            const double phi = phi0 + dphi_per_ds * s;
            const double x = cx + R_cm * std::cos(phi);
            const double z = cz + R_cm * std::sin(phi);
            const double y = vy0 + dy_per_ds * s;
            if (z > zmax + 0.5)
                break;
            gHxXY->SetPoint(gHxXY->GetN(), x, y);
            gHxXZ->SetPoint(gHxXZ->GetN(), x, z);
            gHxYZ->SetPoint(gHxYZ->GetN(), z, y);
        }

        TString hdr = Form("evt %lld  p_MC=(%.3f,%.3f,%.3f) GeV  N_GTPC=%d", i, px0, py0, pz0, gtpcPoints->GetEntries());

        c->cd(3 * drawn + 1);
        gPriXY->SetTitle(hdr + ";x (cm, world);y (cm, world)");
        gPriXY->Draw("AP");
        gSecXY->Draw("P SAME");
        gHxXY->Draw("L SAME");

        c->cd(3 * drawn + 2);
        gPriXZ->SetTitle("x vs z;x (cm);z (cm)");
        gPriXZ->Draw("AP");
        gSecXZ->Draw("P SAME");
        gHxXZ->Draw("L SAME");

        c->cd(3 * drawn + 3);
        gPriYZ->SetTitle("z vs y;z (cm);y (cm)");
        gPriYZ->Draw("AP");
        gSecYZ->Draw("P SAME");
        gHxYZ->Draw("L SAME");

        ++drawn;
    }

    // Legend
    auto* leg = new TLegend(0.05, 0.96, 0.95, 0.99);
    leg->SetNColumns(3);
    leg->SetBorderSize(0);
    auto* lp = new TGraph(1); lp->SetMarkerStyle(20); lp->SetMarkerColor(kBlue + 2); lp->SetMarkerSize(1.4);
    auto* ls = new TGraph(1); ls->SetMarkerStyle(25); ls->SetMarkerColor(kOrange + 1); ls->SetMarkerSize(1.2);
    auto* lh = new TGraph(1); lh->SetLineColor(kRed); lh->SetLineWidth(2); lh->SetLineStyle(2);
    leg->AddEntry(lp, "MC primary (trackID==0)", "p");
    leg->AddEntry(ls, "MC secondaries", "p");
    leg->AddEntry(lh, "Analytic helix (B=2T, π MC p)", "l");
    c->cd(0);
    leg->Draw();

    c->SaveAs(Form("mc_only_p%d.png", p_MeV));
    std::cout << "Wrote mc_only_p" << p_MeV << ".png\n";
}
