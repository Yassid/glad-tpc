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

#pragma once

#include "FairTask.h"
#include "R3BGTPCHitData.h"
#include "R3BGTPCTrackData.h"
// #include "R3BGTPCHitPar.h" TrackPar?
#include "R3BGTPCTrackFinder.h"
#include "R3BGTPCTrackFinderRiemann.h"

class R3BGTPCHit2Track : public FairTask
{
  public:
    /** Default constructor **/
    R3BGTPCHit2Track();

    /** Destructor **/
    ~R3BGTPCHit2Track();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    /** Virtual method Reset **/
    virtual void Reset();

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Finish **/
    virtual void Finish();

    /** Accessor to select online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }

    /// Switch between the ported AT-TPC Riemann RANSAC finder and the
    /// legacy triplet-clustering path. Default OFF: in the HYDRA Prototype
    /// single-pion case the triplet path gathers all hits into one cluster
    /// and lets Pratt+GN see the whole arc, giving a tighter R distribution
    /// than the RANSAC's inlier subset. Enable for multi-particle AT-TPC-
    /// style events where Riemann's smart-seed separation matters.
    void SetUseRiemann(Bool_t use) { fUseRiemann = use; }

    /// Pass-through RANSAC knobs (no-ops when fUseRiemann == false).
    R3BGTPCTrackFinderRiemann* GetRiemannFinder() { return fRiemannFinder; }

    /// Pass-through TripClust knobs (no-ops when fUseRiemann == true).
    /// Use this from the macro to call SetTcluster / SetMcluster / etc.
    R3BGTPCTrackFinder* GetTrackFinder() { return fTrackFinder; }

    /// Use the MC primary particle's StartXYZ as the vertex pseudo-hit in
    /// R3BGTPCTrackFinder's Pratt+GN circle fit. The MCTrack branch is
    /// typically dropped by Cal2Hit, so we pull it from a sidecar sim file
    /// indexed by event number (set via SetMCSimFile).
    void SetUseMCVertex(Bool_t use) { fUseMCVertex = use; }
    void SetVertexPdg(Int_t pdg) { fVertexPdg = pdg; }
    void SetVertexSigmaCm(Double_t s) { fVertexSigmaCm = s; }
    void SetMCSimFile(const TString& path) { fMCSimFile = path; }

  private:
    void SetParameter();

    TClonesArray* fHitCA;
    TClonesArray* fTrackCA;
    TClonesArray* fMCTrackCA{ nullptr };
    TFile*  fMCSimFilePtr{ nullptr };
    TTree*  fMCSimTree{ nullptr };
    Long64_t fEventCounter{ -1 };
    TString fMCSimFile{};

    Bool_t fOnline;        // Selector for online data storage
    Bool_t fUseRiemann;    // Use the ported AT-TPC Riemann RANSAC finder
    Bool_t fUseMCVertex{ kFALSE };
    Int_t  fVertexPdg{ -211 };
    Double_t fVertexSigmaCm{ 0.05 }; // 0.5 mm

    R3BGTPCTrackFinder* fTrackFinder{};
    R3BGTPCTrackFinderRiemann* fRiemannFinder{};

    ClassDef(R3BGTPCHit2Track, 3);
};
