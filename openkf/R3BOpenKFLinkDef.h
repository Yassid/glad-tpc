// clang-format off
/******************************************************************************
 *   Copyright (C) 2018-2026 Members of the R3B Collaboration                 *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *               GNU Lesser General Public Licence (GPL) version 3,           *
 *                    copied verbatim in the file "LICENSE".                  *
 ******************************************************************************/

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace AtTools;
#pragma link C++ namespace kf;

#pragma link C++ class AtTools::AtELossModel - !;
#pragma link C++ class AtTools::AtELossCATIMA - !;
#pragma link C++ class AtTools::AtKinematics + ;
#pragma link C++ class AtTools::AtPropagator - !;
#pragma link C++ class kf::TrackFitterUKFBase - !;
#pragma link C++ class kf::TrackFitterUKF - !;

#endif
