/******************************************************************************
 *   Copyright (C) 2018-2026 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 ******************************************************************************/

#pragma once

#include "FairTask.h"

#include <Math/Vector3D.h>
#include <memory>

class R3BGTPCFitterUKF;
class TClonesArray;
namespace AtTools
{
class AtELossCATIMA;
}

/// FairTask wrapping R3BGTPCFitterUKF. Reads GTPCTrackData from the tree,
/// runs the UKF per track, writes GTPCFittedTrackData. Designed to plug
/// directly into the standard glad-tpc reconstruction chain after
/// R3BGTPCHit2Track. Configuration follows the existing tasks (Init /
/// SetParContainers / Exec / Reset / Finish hooks).
class R3BGTPCTrack2Fit : public FairTask
{
  public:
    R3BGTPCTrack2Fit();
    ~R3BGTPCTrack2Fit() override;

    InitStatus Init() override;
    InitStatus ReInit() override;
    void SetParContainers() override;
    void Exec(Option_t* opt) override;
    void Reset();
    void Finish() override;

    // --- particle hypothesis (defaults to a proton)
    void SetParticle(double charge_C, double mass_MeV)
    {
        fCharge = charge_C;
        fMass_MeV = mass_MeV;
    }

    // --- gas (defaults to H₂ at HYDRA-like 3.553e-5 g/cm³). The vector lists
    //     (A, Z, stoichiometry) tuples handed to catima::Material.
    void SetGas(double density_g_cm3, std::vector<std::tuple<int, int, int>> components)
    {
        fGasDensity = density_g_cm3;
        fGasComponents = std::move(components);
    }
    void SetProjectile(double A, double Z, double massAmu)
    {
        fProjectileA = A;
        fProjectileZ = Z;
        fProjectileMassAmu = massAmu;
    }

    // --- fitter knobs (forwarded to R3BGTPCFitterUKF)
    void SetBField(const ROOT::Math::XYZVector& B) { fBField = B; }
    void SetMomentumSeed(double p_MeV) { fMomentumSeed = p_MeV; }
    void SetMeasurementSigma(double sigma_mm) { fMeasSigma_mm = sigma_mm; }
    void SetMomentumSigmaFrac(double frac) { fMomSigmaFrac = frac; }
    void SetMinClusters(int n) { fMinClusters = n; }
    void SetEnableEnergyStraggling(bool on) { fEnableEnStraggling = on; }
    void SetELossScaleFactor(double f) { fELossScaleFactor = f; }
    void SetOnline(bool on) { fOnline = on; }

    void SetInputBranch(const char* name) { fInputBranch = name; }
    void SetOutputBranch(const char* name) { fOutputBranch = name; }

  private:
    TClonesArray* fTrackCA{ nullptr };
    TClonesArray* fFittedCA{ nullptr };
    std::unique_ptr<R3BGTPCFitterUKF> fFitter;

    // particle + gas
    double fCharge{ 1.602176634e-19 };
    double fMass_MeV{ 938.272 };
    double fProjectileA{ 1 }, fProjectileZ{ 1 }, fProjectileMassAmu{ 1.0 };
    double fGasDensity{ 3.553e-5 };
    std::vector<std::tuple<int, int, int>> fGasComponents{ std::make_tuple(1, 1, 1) };

    // fitter knobs
    ROOT::Math::XYZVector fBField{ 0., 0., 2.0 };
    double fMomentumSeed{ -1.0 };
    double fMeasSigma_mm{ 1.0 };
    double fMomSigmaFrac{ 0.1 };
    int fMinClusters{ 5 };
    bool fEnableEnStraggling{ true };
    double fELossScaleFactor{ 1.0 };

    bool fOnline{ false };
    TString fInputBranch{ "GTPCTrackData" };
    TString fOutputBranch{ "GTPCFittedTrackData" };

  public:
    ClassDefOverride(R3BGTPCTrack2Fit, 1);
};
