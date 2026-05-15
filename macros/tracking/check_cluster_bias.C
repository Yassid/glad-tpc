/// @file check_cluster_bias.C
/// @brief Test the cluster-chord-pull hypothesis. For each cluster, find the
///        nearest MC primary point, then measure the cluster's perpendicular
///        offset from the MC-defined helix segment. If clusters lie
///        systematically toward the *concave* side of the arc, that proves
///        the centroiding bias and explains the +17 % UKF momentum bias.

void check_cluster_bias()
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);

    const double offX = 4.2, offZ = 260.2;

    TString workDir = gSystem->Getenv("VMCWORKDIR");
    TFile fSim(workDir + "/glad-tpc/macros/sim/Prototype/sim.root");
    TFile fTrk(workDir + "/glad-tpc/macros/tracking/output_tracking.root");
    auto* tSim = (TTree*)fSim.Get("evt");
    auto* tTrk = (TTree*)fTrk.Get("evt");
    auto* pts = new TClonesArray("R3BGTPCPoint");
    auto* mc = new TClonesArray("R3BMCTrack");
    auto* trks = new TClonesArray("R3BGTPCTrackData");
    tSim->SetBranchAddress("GTPCPoint", &pts);
    tSim->SetBranchAddress("MCTrack", &mc);
    tTrk->SetBranchAddress("GTPCTrackData", &trks);

    // For each cluster on a long-chord track: compute residual to a 2-point
    // MC chord (nearest-z primary point pair), measured perpendicular to that
    // chord. Sign convention: + means cluster is on the "outside" of the arc
    // (i.e. opposite from the curvature centre); - means "inside" (chord).
    //
    // The Kasa centre (cx, cz) tells us which side is "inside". Sign of
    // residual = sign( (cluster - chord_midpoint) · (chord_midpoint - centre) ).
    auto* hRes = new TH1F("hRes", "cluster residual to MC chord (+ = away from centre);residual [mm];clusters", 100, -5, 5);
    auto* hVsLen = new TH2F("hVsLen", ";chord length (cm);residual [mm]", 30, 0, 30, 80, -5, 5);
    int total = 0;
    for (Long64_t i = 0; i < std::min(tSim->GetEntries(), tTrk->GetEntries()); ++i) {
        tSim->GetEntry(i);
        tTrk->GetEntry(i);
        if (trks->GetEntries() == 0) continue;
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        auto& cls = *tr->GetHitClusterArray();
        if (cls.size() < 5) continue;

        // MC π primary points, sorted by z
        std::vector<std::array<double, 3>> mcPts;
        for (int j = 0; j < pts->GetEntries(); ++j) {
            auto* p = (R3BGTPCPoint*)pts->At(j);
            if (std::abs(p->GetPDGCode()) != 211) continue;
            mcPts.push_back({ p->GetX() - offX, p->GetY(), p->GetZ() - offZ });
        }
        if (mcPts.size() < 2) continue;
        std::sort(mcPts.begin(), mcPts.end(), [](auto& a, auto& b) { return a[2] < b[2]; });

        // chord length
        double zmn = 1e9, zmx = -1e9, xmn = 1e9, xmx = -1e9;
        for (auto& c : cls) { zmn = std::min(zmn, (double)c.GetZ()); zmx = std::max(zmx, (double)c.GetZ());
                              xmn = std::min(xmn, (double)c.GetX()); xmx = std::max(xmx, (double)c.GetX()); }
        const double chord_cm = std::sqrt((zmx - zmn) * (zmx - zmn) + (xmx - xmn) * (xmx - xmn));

        // Curvature centre side (cx, cz) in (x, z), pad-local
        auto gc = tr->GetGeoCenter();
        const double cx = gc.first;
        const double cz = gc.second;

        for (auto& c : cls) {
            // Find two bracketing MC points by z
            double cz_clust = c.GetZ();
            int ia = -1;
            for (size_t k = 0; k + 1 < mcPts.size(); ++k) {
                if (mcPts[k][2] <= cz_clust && cz_clust <= mcPts[k + 1][2]) { ia = k; break; }
            }
            if (ia < 0) continue; // cluster outside MC z range
            const double x1 = mcPts[ia][0], z1 = mcPts[ia][2];
            const double x2 = mcPts[ia + 1][0], z2 = mcPts[ia + 1][2];
            // Perp distance from cluster to chord (x2-x1, z2-z1)
            const double dx_chord = x2 - x1, dz_chord = z2 - z1;
            const double len = std::sqrt(dx_chord * dx_chord + dz_chord * dz_chord);
            if (len < 1e-6) continue;
            // perp vector = (-dz_chord, dx_chord)/len (rotated 90°)
            // residual = (cluster - mid) · perp
            const double midx = 0.5 * (x1 + x2), midz = 0.5 * (z1 + z2);
            const double dx_c = c.GetX() - midx, dz_c = c.GetZ() - midz;
            const double res_cm = (dx_c * -dz_chord + dz_c * dx_chord) / len;
            // Sign: + if cluster is on the OUTSIDE (away from centre)
            // Direction from chord midpoint to centre: (cx - midx, cz - midz)
            const double dirCx = cx - midx, dirCz = cz - midz;
            // Project this onto perp direction; if perp · centre_dir > 0,
            // perp points TOWARD centre; we want to flip the sign so + means away.
            const double sign = (((-dz_chord) * dirCx + dx_chord * dirCz) > 0) ? -1.0 : +1.0;
            const double res_mm_signed = sign * res_cm * 10.0;
            hRes->Fill(res_mm_signed);
            hVsLen->Fill(chord_cm, res_mm_signed);
            ++total;
        }
    }
    std::cout << "clusters analyzed: " << total << "\n";
    std::cout << "  mean residual: " << hRes->GetMean() << " mm\n";
    std::cout << "  RMS:           " << hRes->GetRMS() << " mm\n";
    auto* c = new TCanvas("c", "cluster residual", 1400, 500);
    c->Divide(2, 1);
    c->cd(1); hRes->Draw();
    c->cd(2); hVsLen->Draw("COLZ");
    c->SaveAs("cluster_chord_residual.png");
}
