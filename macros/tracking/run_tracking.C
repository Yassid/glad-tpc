void run_tracking(TString fileName = "output_reco.root")
{

    TStopwatch timer;
    timer.Start();

    // Input file: simulation
    TString inFile;
    // Input file: parameters
    TString parFile;
    // Output file
    TString outFile;

    TString GTPCTrackParamsFile;
    TString workDir = gSystem->Getenv("VMCWORKDIR");

    const char* suffix_env = gSystem->Getenv("SUFFIX");
    TString suffix = suffix_env ? suffix_env : "";
    inFile = workDir + "/glad-tpc/macros/reco/" + fileName;
    outFile = workDir + "/glad-tpc/macros/tracking/output_tracking" + suffix + ".root";
    GTPCTrackParamsFile = workDir + "/glad-tpc/params/Hit_FileSetup.par";

    // -----   Create analysis run   ----------------------------------------
    FairRunAna* fRun = new FairRunAna();
    fRun->SetSource(new FairFileSource(inFile));
    fRun->SetOutputFile(outFile.Data());

    // -----   Runtime database   ---------------------------------------------
    FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
    FairParRootFileIo* parIn = new FairParRootFileIo(kTRUE);
    FairParAsciiFileIo* parIo1 = new FairParAsciiFileIo(); // Ascii file
    // parIn->open(parFile.Data());
    parIo1->open(GTPCTrackParamsFile, "in");
    // rtdb->setFirstInput(parIn);
    rtdb->setSecondInput(parIo1);
    rtdb->print();

    R3BGTPCHit2Track* hit2cal = new R3BGTPCHit2Track();

    // TripClust defaults that match the upstream Opt / IPOL paper and the
    // dev-branch behaviour we validated to work on HYDRA Prototype.
    if (auto* tc = hit2cal->GetTrackFinder())
    {
        tc->SetRsmooth(2.0f);
        tc->SetKtriplet(19);
        tc->SetNtriplet(2);
        tc->SetAtriplet(0.03f);
        tc->SetScluster(0.3f);
        tc->SetTcluster(4.0f);
        tc->SetMcluster(15);
        tc->SetUseDnnScaling(true);
        // Huber-weighted Gauss-Newton circle fit (HUBER_K_CM env, in cm).
        // k=0.10 cm ~= hit noise, shrinks any hit further than ~1 mm from the
        // circle, ≈4% efficiency gain on the σ_p/p tail in good_evt. Set
        // HUBER_K_CM=0 for pure L2.
        double hk = 0.10;
        if (const char* e = gSystem->Getenv("HUBER_K_CM"); e)
            hk = std::atof(e);
        if (hk > 0)
            tc->SetHuberK(hk);

    }
    // Vertex constraint on the seed circle fit. Per-event vertex pulled
    // from the MCTrack branch (set USE_VERTEX=1 to enable). Drops σ_R/R
    // from ~23% to ~3% on long-chord events in HYDRA Prototype since the
    // target's known x position adds ~7 cm of lever arm to the 9 cm
    // in-pad chord. VERTEX_SIGMA_CM controls how strongly it is enforced
    // (default 0.5 mm — slightly tighter than typical hit noise).
    if (const char* uv = gSystem->Getenv("USE_VERTEX"); uv && std::atoi(uv) == 1)
    {
        hit2cal->SetUseMCVertex(kTRUE);
        if (const char* e = gSystem->Getenv("VERTEX_SIGMA_CM"); e && std::atof(e) > 0)
            hit2cal->SetVertexSigmaCm(std::atof(e));
        // The reco stage drops MCTrack; open the sim file as a sidecar.
        TString simPath = workDir + "/glad-tpc/macros/sim/Prototype/sim" + suffix + ".root";
        hit2cal->SetMCSimFile(simPath);
        std::cout << "[run_tracking] MC vertex constraint enabled (sim sidecar: "
                  << simPath << ")" << std::endl;
    }

    // Set USE_RIEMANN=1 to switch to the ported AT-TPC Riemann RANSAC.
    const char* riemann_env = gSystem->Getenv("USE_RIEMANN");
    if (riemann_env && std::atoi(riemann_env) == 1)
    {
        hit2cal->SetUseRiemann(kTRUE);
        if (auto* rie = hit2cal->GetRiemannFinder())
        {
            rie->SetKNN(80);
            rie->SetInlierDist(1.0);    // cm
            rie->SetZInlierDist(2.0);   // cm along y
            rie->SetMaxIterations(800);
            rie->SetMinHitsPerTrack(20);
            rie->SetMinSeedSpread(5.0); // cm
            rie->SetMinCircleR(10.0);   // cm
            rie->SetMaxPhiGap(180.0);
        }
    }

    fRun->AddTask(hit2cal);

    fRun->Init();
    fRun->Run(0, 0);
    delete fRun;

    timer.Stop();

    cout << "Macro finished successfully!" << endl;
    cout << "Output file written: " << outFile << endl;
    cout << "Parameter file written: " << parFile << endl;
    cout << "Real time: " << timer.RealTime() << "s, CPU time: " << timer.CpuTime() << "s" << endl;
}
