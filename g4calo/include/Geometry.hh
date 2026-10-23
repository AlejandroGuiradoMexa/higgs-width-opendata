// Simplified ATLAS-like electromagnetic barrel calorimeter.
//
// The accordion is approximated by planar Pb/LAr sampling layers. The x axis plays the
// role of eta and y the role of phi at a radius of 1.5 m, so cell sizes follow the ATLAS
// barrel granularity (Delta eta x Delta phi):
//   presampler   0.025  x 0.1
//   strips (S1)  0.003125 x 0.1
//   middle (S2)  0.025  x 0.0245
//   back   (S3)  0.05   x 0.0245
// A simple iron/scintillator hadronic calorimeter (0.1 x 0.1 cells) sits behind the EM layers,
// so that hadronic leakage (R_had) can be used against charged pions.
// The presampler front is at z = 0 and the particle gun 1.5 m upstream, as in the ATLAS barrel.
#pragma once

#include <array>

namespace geo {

constexpr double kRadius = 1500.0;  // mm, used to turn angular cell sizes into lengths

// Absorber and gap thicknesses (mm). Pb 1.5 mm + LAr 4.2 mm gives X0_eff ~ 19 mm.
constexpr double kPb = 1.5;
constexpr double kLAr = 4.2;
constexpr double kPresamplerLAr = 11.0;
constexpr double kFe = 40.0;
constexpr double kScint = 3.0;
constexpr double kGunDistance = 1500.0;

enum Layer { kPS = 0, kS1 = 1, kS2 = 2, kS3 = 3, kHad = 4, kNLayers = 5 };

// Number of absorber/active periods per layer (~4.3, ~16, ~2 X0 of Pb/LAr; ~5 lambda of Fe/scintillator).
constexpr std::array<int, kNLayers> kPeriods = {0, 15, 54, 7, 20};

// Cell sizes along x (eta) and y (phi), in mm.
constexpr std::array<double, kNLayers> kCellX = {0.025 * kRadius, 0.003125 * kRadius, 0.025 * kRadius,
                                                0.05 * kRadius, 0.1 * kRadius};
constexpr std::array<double, kNLayers> kCellY = {0.1 * kRadius, 0.1 * kRadius, 0.0245 * kRadius,
                                                0.0245 * kRadius, 0.1 * kRadius};

// Readout window stored per event, centred on the impact point (odd numbers of cells).
constexpr std::array<int, kNLayers> kNx = {5, 41, 11, 5, 3};
constexpr std::array<int, kNLayers> kNy = {3, 3, 11, 11, 3};

constexpr int kNDepthBins = 30;  // longitudinal profile, one bin per X0 of calorimeter depth

}  // namespace geo
