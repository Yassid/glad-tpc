/******************************************************************************
 * Copyright (C) 2018-2026 GSI Helmholtzzentrum für Schwerionenforschung GmbH *
 *         Copyright (C) 2018-2026 Members of R3B Collaboration               *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *              GNU Lesser General Public Licence (LGPL) version 3,           *
 *                     copied verbatim in the file "LICENSE".                 *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TClonesArray.h"

#include "R3BGTPC.h"
#include "R3BGTPCHit2Track.h"
#include "R3BGTPCHitData.h"
#include "R3BGTPCTrackData.h"
#include "R3BMCTrack.h"
// #include "R3BGTPCHitPar.h"

#include "dnn.h"
#include "graph.h"
#include "option.h"
#include "output.h"
#include "pointcloud.h"

// R3BGTPCHit2Track: Constructor
R3BGTPCHit2Track::R3BGTPCHit2Track()
    : FairTask("R3B GTPC Hit to Track")
    , fHitCA(NULL)
    , fTrackCA(NULL)
    , fOnline(kFALSE)
    , fUseRiemann(kFALSE) // legacy triplet+Pratt is better in single-pion Prototype
{
    // Create the finders eagerly so the macro can configure them before
    // FairRunAna::Init() runs SetParContainers().
    fTrackFinder = new R3BGTPCTrackFinder();
    fRiemannFinder = new R3BGTPCTrackFinderRiemann();
}

R3BGTPCHit2Track::~R3BGTPCHit2Track()
{
    LOG(info) << "R3BGTPCHit2Track: Delete instance";
    if (fHitCA)
        delete fHitCA;
    if (fTrackCA)
        delete fTrackCA;
    delete fTrackFinder;
    delete fRiemannFinder;
}

void R3BGTPCHit2Track::SetParContainers()
{
    // Reading GTPCHitPar from FairRuntimeDb
    // FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    // if (!rtdb)
    // {
    //  LOG(error) << "R3BGTPCMapped2Cal:: FairRuntimeDb not opened!";
    // }

    // fHit_Par = (R3BGTPCHitPar*)rtdb->getContainer("GTPCHitPar");
    // if (!fHit_Par)
    // {
    //  LOG(error) << "R3BGTPCCal2Hit::Init() Couldn't get handle on GTPCHitPar
    //  container";
    // }
    // else
    // {
    //   LOG(info) << "R3BGTPCCal2Hit:: GTPCHitPar container open";
    // }

    // Finders are constructed in the R3BGTPCHit2Track ctor so the macro can
    // configure them before FairRunAna::Init() runs SetParContainers().
}

void R3BGTPCHit2Track::SetParameter()
{
    //--- Parameter Container ---
    /*
          fHitParams = new TArrayF();
          fHitParams = fHit_Par->GetHitParams(); // Array with the Hit parameters
    */
}

InitStatus R3BGTPCHit2Track::Init()
{
    LOG(info) << "R3BGTPCHit2Track::Init() ";
    assert(!fTrackCA); // in case someone calls Init() twice.

    // INPUT DATA - Cal
    FairRootManager* ioManager = FairRootManager::Instance();
    if (!ioManager)
        LOG(fatal) << "Init: No FairRootManager";

    fHitCA = (TClonesArray*)ioManager->GetObject("GTPCHitData");
    if (!fHitCA)
        LOG(fatal) << "Init: No GTPCHitData";

    // Optional: MC truth for vertex constraint. Two routes:
    //   1) MCTrack branch propagated through the chain (cleanest, but the
    //      reco stage typically drops it).
    //   2) Sidecar sim file opened directly, indexed by event counter
    //      (set via SetMCSimFile from the macro).
    fMCTrackCA = (TClonesArray*)ioManager->GetObject("MCTrack");
    if (fUseMCVertex && !fMCTrackCA && fMCSimFile.Length() > 0)
    {
        fMCSimFilePtr = TFile::Open(fMCSimFile);
        if (fMCSimFilePtr && !fMCSimFilePtr->IsZombie())
        {
            fMCSimTree = (TTree*)fMCSimFilePtr->Get("evt");
            if (fMCSimTree)
            {
                fMCTrackCA = new TClonesArray("R3BMCTrack");
                fMCSimTree->SetBranchAddress("MCTrack", &fMCTrackCA);
                LOG(info) << "R3BGTPCHit2Track: MC vertex from sidecar sim file "
                          << fMCSimFile << " (" << fMCSimTree->GetEntries() << " entries)";
            }
        }
    }
    if (fMCTrackCA && fUseMCVertex)
        LOG(info) << "R3BGTPCHit2Track: MC vertex constraint enabled";

    // Register output - Track
    fTrackCA = new TClonesArray("R3BGTPCTrackData", 50);
    if (!fOnline)
    {
        ioManager->Register("GTPCTrackData", "GTPC Track", fTrackCA, kTRUE);
    }
    else
    {
        ioManager->Register("GTPCTrackData", "GTPC Track", fTrackCA, kFALSE);
    }

    SetParameter();
    return kSUCCESS;
}

InitStatus R3BGTPCHit2Track::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BGTPCHit2Track::Exec(Option_t* opt)
{
    Reset(); // Reset entries in output arrays, local arrays

    // Set per-event vertex on the legacy track-finder. Constants from the
    // HYDRA Prototype geometry (target at world x=-2.7, chamber active region
    // at world x=4.2..13.0 → local x_vtx = -6.9). The z is taken from MC
    // truth here; a beam-tracker replacement would slot in identically.
    if (fUseMCVertex && fMCTrackCA && fTrackFinder)
    {
        if (fMCSimTree) fMCSimTree->GetEntry(fEventCounter);
        const double offX = 4.2, offZ = 260.2;
        bool found = false;
        for (Int_t j = 0; j < fMCTrackCA->GetEntries(); ++j)
        {
            auto* m = static_cast<R3BMCTrack*>(fMCTrackCA->At(j));
            if (m->GetMotherId() == -1 && m->GetPdgCode() == fVertexPdg)
            {
                const double vxL = m->GetStartX() - offX;
                const double vzL = m->GetStartZ() - offZ;
                fTrackFinder->SetVertexConstraint(vxL, vzL, fVertexSigmaCm);
                found = true;
                break;
            }
        }
        if (!found) fTrackFinder->DisableVertexConstraint();
    }
    else if (fTrackFinder)
    {
        fTrackFinder->DisableVertexConstraint();
    }

    // Riemann RANSAC path — single-pass, std::vector-based. Skips the entire
    // triplet/PointCloud machinery and persists tracks straight into fTrackCA.
    if (fUseRiemann)
    {
        const Int_t nHits = fHitCA->GetEntries();
        std::vector<R3BGTPCHitData> hits;
        hits.reserve(nHits);
        for (Int_t i = 0; i < nHits; ++i)
        {
            auto* h = static_cast<R3BGTPCHitData*>(fHitCA->At(i));
            if (h)
                hits.push_back(*h);
        }

        auto tracks = fRiemannFinder->FindTracks(hits);

        for (auto& trk : tracks)
        {
            TClonesArray& clref = *fTrackCA;
            const Int_t size = clref.GetEntriesFast();
            auto* persisted = new (clref[size])
                R3BGTPCTrackData(trk.GetTrackId(), std::move(trk.GetHitArray()),
                                 std::vector<R3BGTPCHitClusterData>{});
            persisted->SetGeoCenter(trk.GetGeoCenter());
            persisted->SetGeoRadius(trk.GetGeoRadius());
            persisted->SetGeoTheta(trk.GetGeoTheta());
        }
        return;
    }

    // Read TripClust parameters directly from fTrackFinder so the macro can
    // tune them via the SetScluster/SetKtriplet/... setters.
    const tc_params& tp = fTrackFinder->GetInputParams();
    const int opt_verbose = 0;
    PointCloud cloud_xyz;
    fTrackFinder->eventToClusters(fHitCA, cloud_xyz);

    if (cloud_xyz.size() == 0)
    {
        std::cerr << "[Error] empty cloud " << std::endl;
        return;
    }

    // Upstream Opt defaults: r and s are multiples of dNN (the data's first-
    // quartile nearest-neighbour distance). Replicate that here, gated by
    // fTrackFinder->UseDnnScaling() so a macro can pin absolute values.
    double r_eff = tp.r;
    double s_eff = tp.s;
    if (fTrackFinder->UseDnnScaling())
    {
        const double dnn = std::sqrt(first_quartile(cloud_xyz));
        if (dnn == 0.0)
        {
            std::cerr << "[Error] dnn computed as zero. "
                      << "Suggestion: remove doublets, e.g. with 'sort -u'" << std::endl;
            return;
        }
        r_eff = tp.r * dnn;
        s_eff = tp.s * dnn;
    }

    // Step 1) smoothing by position averaging of neighboring points
    PointCloud cloud_xyz_smooth;
    smoothen_cloud(cloud_xyz, cloud_xyz_smooth, r_eff);

    // Step 2) finding triplets of approximately collinear points
    std::vector<triplet> triplets;
    generate_triplets(cloud_xyz_smooth, triplets, tp.k, tp.n, tp.a);

    // Step 3) single link hierarchical clustering of the triplets. Match the
    // upstream Opt defaults: t is literal (tauto=false), no dmax filter,
    // single linkage.
    cluster_group cl_group;
    if (cloud_xyz_smooth.size() < 10)
        return;
    compute_hc(cloud_xyz_smooth,
               cl_group,
               triplets,
               s_eff,
               tp.t,
               /*tauto=*/false,
               /*dmax=*/0.0,
               /*isdmax=*/false,
               /*linkage=*/SINGLE,
               opt_verbose);

    // Step 4) pruning by removal of small clusters
    cleanup_cluster_group(cl_group, tp.m, opt_verbose);
    cluster_triplets_to_points(triplets, cl_group);

    // store cluster labels in points
    add_clusters(cloud_xyz, cl_group, /*gnuplot=*/false);

    // Adapt clusters to AtTrack
    fTrackFinder->clustersToTrack(cloud_xyz, cl_group, fTrackCA, fHitCA);
    return;
}

void R3BGTPCHit2Track::Finish() {}

void R3BGTPCHit2Track::Reset()
{
    LOG(debug) << "Clearing TrackData Structure";
    if (fTrackCA)
        fTrackCA->Clear();
    ++fEventCounter;
}

ClassImp(R3BGTPCHit2Track)
