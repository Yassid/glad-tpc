/******************************************************************************
 *   Copyright (C) 2018-2026 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 ******************************************************************************/

#include "R3BGTPCTrack2Fit.h"

#include "AtELossCATIMA.h"
#include "R3BGTPCFittedTrackData.h"
#include "R3BGTPCFitterUKF.h"
#include "R3BGTPCTrackData.h"

#include <FairLogger.h>
#include <FairRootManager.h>
#include <TClonesArray.h>

#include <cassert>
#include <utility>

R3BGTPCTrack2Fit::R3BGTPCTrack2Fit()
    : FairTask("R3B GTPC Track to UKF Fit")
{
}

R3BGTPCTrack2Fit::~R3BGTPCTrack2Fit()
{
    LOG(info) << "R3BGTPCTrack2Fit: Delete instance";
    delete fFittedCA;
    fFittedCA = nullptr;
}

void R3BGTPCTrack2Fit::SetParContainers()
{
    // Build the UKF lazily here so any SetX() calls invoked between
    // construction and Init still take effect.
    auto eloss = std::make_unique<AtTools::AtELossCATIMA>(fGasDensity);
    eloss->SetProjectile(fProjectileA, fProjectileZ, fProjectileMassAmu);
    eloss->SetMaterial(fGasComponents);

    fFitter = std::make_unique<R3BGTPCFitterUKF>(fCharge, fMass_MeV, std::move(eloss));
    fFitter->SetBField(fBField);
    if (fMomentumSeed > 0)
        fFitter->SetMomentumSeed(fMomentumSeed);
    fFitter->SetMeasurementSigma(fMeasSigma_mm);
    fFitter->SetMomentumSigmaFrac(fMomSigmaFrac);
    fFitter->SetMinClusters(fMinClusters);
    fFitter->SetMaxSeedRadius_cm(fMaxSeedRadius_cm);
    fFitter->SetEnableEnergyStraggling(fEnableEnStraggling);
    fFitter->SetELossScaleFactor(fELossScaleFactor);
    fFitter->SetInputUnit_mm(fInputUnit_mm);
    fFitter->SetUseHelixBackExtrap(fUseHelixBackExtrap);
    fFitter->SetForceVertexOnBeamAxis(fForceVertexOnBeamAxis);
    fFitter->SetUpdateAnglesOnBackExtrap(fUpdateAnglesOnBackExtrap);
    fFitter->SetBackExtrapMaxPath(fBackExtrapMaxPath);
    fFitter->SetBackExtrapTargetX(fBackExtrapTargetX);
}

InitStatus R3BGTPCTrack2Fit::Init()
{
    LOG(info) << "R3BGTPCTrack2Fit::Init()";
    assert(!fFittedCA); // catch double-init

    auto* mgr = FairRootManager::Instance();
    if (!mgr)
        LOG(fatal) << "R3BGTPCTrack2Fit::Init: No FairRootManager";

    fTrackCA = (TClonesArray*)mgr->GetObject(fInputBranch);
    if (!fTrackCA)
        LOG(fatal) << "R3BGTPCTrack2Fit::Init: No input branch '" << fInputBranch << "'";

    fFittedCA = new TClonesArray("R3BGTPCFittedTrackData", 50);
    mgr->Register(fOutputBranch, "GTPC fitted tracks (UKF)", fFittedCA, fOnline ? kFALSE : kTRUE);

    if (!fFitter)
        SetParContainers();

    return kSUCCESS;
}

InitStatus R3BGTPCTrack2Fit::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BGTPCTrack2Fit::Reset()
{
    if (fFittedCA)
        fFittedCA->Clear("C");
}

void R3BGTPCTrack2Fit::Exec(Option_t* /*opt*/)
{
    Reset();
    if (!fTrackCA || !fFitter)
        return;

    const Int_t nTracks = fTrackCA->GetEntriesFast();
    Int_t nFit = 0;
    for (Int_t i = 0; i < nTracks; ++i)
    {
        auto* track = (R3BGTPCTrackData*)fTrackCA->At(i);
        if (!track)
            continue;
        auto fitted = fFitter->FitTrack(track);
        if (!fitted)
            continue;
        const Int_t outIdx = fFittedCA->GetEntriesFast();
        new ((*fFittedCA)[outIdx]) R3BGTPCFittedTrackData(*fitted);
        ++nFit;
    }
    LOG(info) << "R3BGTPCTrack2Fit::Exec: fitted " << nFit << "/" << nTracks << " tracks";
}

void R3BGTPCTrack2Fit::Finish() {}

ClassImp(R3BGTPCTrack2Fit);
