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
    // dev-branch behaviour we validated to work on HYDRA Prototype. Listed
    // explicitly so they're discoverable from the macro; passing them
    // through the setters is a no-op vs the built-in defaults.
    if (auto* tc = hit2cal->GetTrackFinder())
    {
        tc->SetRsmooth(2.0f);    // r — smoothing radius (× dNN)
        tc->SetKtriplet(19);     // k — neighbours per triplet
        tc->SetNtriplet(2);      // n — best triplets per midpoint
        tc->SetAtriplet(0.03f);  // a — max (1 - cos α)
        tc->SetScluster(0.3f);   // s — clustering scale (× dNN)
        tc->SetTcluster(4.0f);   // t — cluster-distance threshold
        tc->SetMcluster(15);     // m — min triplets per cluster
        tc->SetUseDnnScaling(true);
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
