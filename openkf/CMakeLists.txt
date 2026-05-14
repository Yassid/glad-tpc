##############################################################################
#   Copyright (C) 2018-2026 Members of the R3B Collaboration                 #
#                                                                            #
#            This software is distributed under the terms of the             #
#              GNU Lesser General Public Licence (GPL) version 3,            #
#                    copied verbatim in the file "LICENSE".                  #
##############################################################################

# libR3BOpenKF: ATTPCROOT's OpenKF + propagator/eloss subset, vendored for use
# by R3BGTPCFitterUKF and any other R3B Kalman fitter that wants UKF semantics.
#
# Contents:
#   - kalman_filter/  : OpenKF templates + TrackFitterUKF (concrete UKF for a
#                       6-state (x,y,z,p,theta,phi) helix in B-field).
#   - AtPropagator    : RK4-stepped Lorentz-force + dE/dx propagator.
#   - AtELossModel/CATIMA : pluggable energy-loss model (CATIMA wrapper).
#   - AtKinematics    : KE/momentum conversions.
# Namespaces preserved as AtTools::... and kf::... so the vendored sources are
# byte-identical to ATTPCROOT — easier to keep in sync.

set(SRCS
    AtPropagator.cxx
    AtELossModel.cxx
    AtELossCATIMA.cxx
    AtKinematics.cxx
    kalman_filter/TrackFitterUKF.cxx
)

set(HEADERS
    AtPropagator.h
    AtELossModel.h
    AtELossCATIMA.h
    AtKinematics.h
    kf_util.h
    types.h
    kalman_filter/kalman_filter.h
    kalman_filter/unscented_kalman_filter.h
    kalman_filter/TrackFitterUKF.h
)

add_library_with_dictionary(
    LIBNAME
        R3BOpenKF
    LINKDEF
        R3BOpenKFLinkDef.h
    HEADERS
        ${HEADERS}
    SRCS
        ${SRCS}
    INCLUDEDIRS
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/kalman_filter
    DEPENDENCIES
        FairRoot::Base
        ROOT::Core
        ROOT::Matrix
        ROOT::MathCore
        ROOT::Physics
        Eigen3::Eigen
        CATIMA::catima
)
