// Standalone analysis: read tracking output, refit each track's circle in
// (x, z) using Pratt + Gauss-Newton with an ADDITIONAL vertex constraint at
// (x_v=-6.9, z_v=MC vertex z). Compare R_new/R_truth to baseline.
//
// This is the feasibility study for adding a beam-tracker vertex pseudo-hit
// to the seed circle fit. The vertex x is fixed by target geometry; z is the
// per-event beam spot. For real data, z would come from a beam tracker; here
// we use MC truth as the cleanest test.
struct Circle { double cx, cz, R; bool valid; };

static Circle pratt_gn(const std::vector<std::pair<double,double>>& xy,
                       bool useVertex, double vx, double vz, double w_vtx)
{
    const std::size_t n = xy.size();
    if (n < 3) return {0,0,0,false};
    double xb = 0, zb = 0;
    for (auto& p : xy) { xb += p.first; zb += p.second; }
    xb /= n; zb /= n;
    double Mxx=0, Myy=0, Mxy=0, Mxz=0, Myz=0, Mzz=0;
    for (auto& p : xy) {
        double dx = p.first - xb, dy = p.second - zb;
        double z2 = dx*dx + dy*dy;
        Mxx += dx*dx; Myy += dy*dy; Mxy += dx*dy;
        Mxz += dx*z2; Myz += dy*z2; Mzz += z2*z2;
    }
    double nd = (double)n;
    Mxx/=nd; Myy/=nd; Mxy/=nd; Mxz/=nd; Myz/=nd; Mzz/=nd;
    double Mz = Mxx + Myy;
    double Cov = Mxx*Myy - Mxy*Mxy;
    double Mxz2 = Mxz*Mxz, Myz2 = Myz*Myz;
    double A3 = 4.0*Mz;
    double A2 = -3.0*Mz*Mz - Mzz;
    double A1 = Mzz*Mz + 4.0*Cov*Mz - Mxz2 - Myz2 - Mz*Mz*Mz;
    double A0 = Mxz2*Myy + Myz2*Mxx - Mzz*Cov - 2*Mxz*Myz*Mxy + Mz*Mz*Cov;
    double A22 = A2+A2, A33 = A3+A3+A3;
    double t = 0, y = 1e20;
    for (int it = 0; it < 100; ++it) {
        double y_old = y;
        y = A0 + t*(A1 + t*(A2 + t*A3));
        if (std::abs(y) > std::abs(y_old)) break;
        double Dy = A1 + t*(A22 + t*A33);
        if (std::abs(Dy) < 1e-30) break;
        double t_old = t;
        t = t_old - y/Dy;
        if (std::abs(t-t_old) < 1e-12*std::abs(t==0?1.0:t)) break;
        if (t < 0) { t = 0; break; }
    }
    double DET = t*t - t*Mz + Cov;
    if (std::abs(DET) < 1e-30) return {0,0,0,false};
    double Xc = (Mxz*(Myy-t) - Myz*Mxy)/(DET*2.0);
    double Yc = (Myz*(Mxx-t) - Mxz*Mxy)/(DET*2.0);
    double cx = Xc + xb, cz = Yc + zb;
    double R2 = Xc*Xc + Yc*Yc + Mz + 2*t;
    if (!(R2 > 0)) return {0,0,0,false};
    double R = std::sqrt(R2);
    for (int it = 0; it < 30; ++it) {
        double Jxx=0,Jxy=0,Jxr=0,Jyy=0,Jyr=0,Jrr=0,bx=0,by=0,br=0;
        for (auto& p : xy) {
            double dx = p.first - cx, dy = p.second - cz;
            double d = std::sqrt(dx*dx+dy*dy);
            if (d < 1e-9) continue;
            double r = d - R;
            double Jx = -dx/d, Jy = -dy/d, JR = -1;
            Jxx += Jx*Jx; Jxy += Jx*Jy; Jxr += Jx*JR;
            Jyy += Jy*Jy; Jyr += Jy*JR; Jrr += JR*JR;
            bx += Jx*r; by += Jy*r; br += JR*r;
        }
        if (useVertex && w_vtx > 0) {
            double dx = vx - cx, dy = vz - cz;
            double d = std::sqrt(dx*dx + dy*dy);
            if (d > 1e-9) {
                double r = d - R;
                double Jx = -dx/d, Jy = -dy/d, JR = -1;
                Jxx += w_vtx*Jx*Jx; Jxy += w_vtx*Jx*Jy; Jxr += w_vtx*Jx*JR;
                Jyy += w_vtx*Jy*Jy; Jyr += w_vtx*Jy*JR; Jrr += w_vtx*JR*JR;
                bx += w_vtx*Jx*r; by += w_vtx*Jy*r; br += w_vtx*JR*r;
            }
        }
        double dt = Jxx*(Jyy*Jrr - Jyr*Jyr) - Jxy*(Jxy*Jrr - Jyr*Jxr) + Jxr*(Jxy*Jyr - Jyy*Jxr);
        if (std::abs(dt) < 1e-12) break;
        double dxc = -(bx*(Jyy*Jrr - Jyr*Jyr) - Jxy*(by*Jrr - Jyr*br) + Jxr*(by*Jyr - Jyy*br))/dt;
        double dyc = -(Jxx*(by*Jrr - Jyr*br) - bx*(Jxy*Jrr - Jyr*Jxr) + Jxr*(Jxy*br - by*Jxr))/dt;
        double dRc = -(Jxx*(Jyy*br - by*Jyr) - Jxy*(Jxy*br - by*Jxr) + bx*(Jxy*Jyr - Jyy*Jxr))/dt;
        cx += dxc; cz += dyc; R += dRc;
        if (std::abs(dxc)+std::abs(dyc)+std::abs(dRc) < 1e-6) break;
    }
    if (!(R > 0) || !std::isfinite(R)) return {0,0,0,false};
    return {cx, cz, R, true};
}

void vertex_refit(double w_vtx_input = 1.0, TString tag = "goodevt")
{
    gSystem->Load("libR3BData");
    gSystem->Load("libR3BGTPCData");
    gStyle->SetOptStat(0);
    TString wd = gSystem->Getenv("VMCWORKDIR");
    TFile fS(wd + "/glad-tpc/macros/sim/Prototype/sim_" + tag + ".root");
    TFile fT("output_tracking_" + tag + ".root");
    auto* tS = (TTree*)fS.Get("evt"); auto* tT = (TTree*)fT.Get("evt");
    auto* mc = new TClonesArray("R3BMCTrack");
    auto* trks = new TClonesArray("R3BGTPCTrackData");
    tS->SetBranchAddress("MCTrack", &mc);
    tT->SetBranchAddress("GTPCTrackData", &trks);

    const double offX = 4.2, offZ = 260.2;
    const double B = 2.0;

    struct Bin { int n=0; double sLog=0, sLog2=0; TH1F* h; };
    auto mkBin = [](const char* nm){ Bin b; b.h = new TH1F(nm,";R/R_{tru}-1;evt",80,-1,1); return b; };
    Bin baseLong  = mkBin("baseL"),  vtxLong  = mkBin("vtxL");
    Bin baseMid   = mkBin("baseM"),  vtxMid   = mkBin("vtxM");

    for (Long64_t i = 0; i < std::min(tS->GetEntries(), tT->GetEntries()); ++i) {
        tS->GetEntry(i); tT->GetEntry(i);
        if (trks->GetEntries() == 0) continue;
        auto* tr = (R3BGTPCTrackData*)trks->At(0);
        auto& hits = tr->GetHitArray();
        if (hits.size() < 10) continue;

        R3BMCTrack* pi = nullptr;
        for (int j = 0; j < mc->GetEntries(); ++j) {
            auto* m = (R3BMCTrack*)mc->At(j);
            if (m->GetMotherId()==-1 && m->GetPdgCode()==-211) { pi = m; break; }
        }
        if (!pi) continue;
        double pT = std::sqrt(pi->GetPx()*pi->GetPx() + pi->GetPz()*pi->GetPz())*1000;
        double R_truth = pT / (3.0*B);
        double vxL = pi->GetStartX() - offX;
        double vzL = pi->GetStartZ() - offZ;

        std::vector<std::pair<double,double>> xy;
        xy.reserve(hits.size());
        double xmn=1e9,xmx=-1e9,zmn=1e9,zmx=-1e9;
        for (auto& h : hits) {
            xy.emplace_back(h.GetX(), h.GetZ());
            xmn=std::min(xmn,(double)h.GetX()); xmx=std::max(xmx,(double)h.GetX());
            zmn=std::min(zmn,(double)h.GetZ()); zmx=std::max(zmx,(double)h.GetZ());
        }
        double chord = std::sqrt((xmx-xmn)*(xmx-xmn)+(zmx-zmn)*(zmx-zmn));

        Circle cBase = pratt_gn(xy, false, 0, 0, 0);
        Circle cVtx  = pratt_gn(xy, true,  vxL, vzL, w_vtx_input);
        if (!cBase.valid || !cVtx.valid) continue;

        double rB = cBase.R / R_truth;
        double rV = cVtx.R  / R_truth;

        if (chord >= 16) {
            baseLong.n++; baseLong.sLog += std::log(rB); baseLong.sLog2 += std::log(rB)*std::log(rB); baseLong.h->Fill(rB-1);
            vtxLong .n++; vtxLong .sLog += std::log(rV); vtxLong .sLog2 += std::log(rV)*std::log(rV); vtxLong .h->Fill(rV-1);
        } else if (chord >= 12) {
            baseMid.n++; baseMid.sLog += std::log(rB); baseMid.sLog2 += std::log(rB)*std::log(rB); baseMid.h->Fill(rB-1);
            vtxMid .n++; vtxMid .sLog += std::log(rV); vtxMid .sLog2 += std::log(rV)*std::log(rV); vtxMid .h->Fill(rV-1);
        }
    }

    auto report = [](const char* lbl, Bin& b){
        if (b.n == 0) return;
        double m=b.sLog/b.n, s=std::sqrt(std::max(0.0, b.sLog2/b.n - m*m));
        printf("  %-30s  N=%3d  median=%.3f  log-rms=%.3f", lbl, b.n, std::exp(m), s);
        b.h->Fit("gaus","Q","",-0.5,0.5);
        if (auto* f = b.h->GetFunction("gaus"))
            printf("   gauss mean=%+.3f sigma=%.3f", f->GetParameter(1), f->GetParameter(2));
        printf("\n");
    };
    printf("=== Vertex refit comparison (w_vtx = %.2f) ===\n", w_vtx_input);
    printf("Long chord (>=16):\n");
    report("baseline (no vertex)", baseLong);
    report("with vertex (MC z)",   vtxLong);
    printf("Mid chord (12-16):\n");
    report("baseline (no vertex)", baseMid);
    report("with vertex (MC z)",   vtxMid);

    TCanvas c("c","",1200,500); c.Divide(2,1);
    c.cd(1); baseLong.h->SetLineColor(kRed+1); baseLong.h->Draw();
             vtxLong.h->SetLineColor(kBlue+1); vtxLong.h->Draw("SAME");
    c.cd(2); baseMid.h->SetLineColor(kRed+1); baseMid.h->Draw();
             vtxMid.h->SetLineColor(kBlue+1); vtxMid.h->Draw("SAME");
    c.SaveAs("vertex_refit.png");
}
