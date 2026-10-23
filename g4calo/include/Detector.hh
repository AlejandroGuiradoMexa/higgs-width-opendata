#pragma once

#include <map>

#include "G4VUserDetectorConstruction.hh"
#include "Geometry.hh"

class G4LogicalVolume;

class Detector : public G4VUserDetectorConstruction {
 public:
  // upstreamX0: aluminium in front of the calorimeter, in radiation lengths (inner detector,
  // solenoid and cryostat in ATLAS amount to ~2-5 X0 depending on eta).
  explicit Detector(double upstreamX0) : fUpstreamX0(upstreamX0) {}

  G4VPhysicalVolume* Construct() override;

  // Front face of the presampler and of the first lead plate (z, mm).
  double PresamplerFront() const { return fPSFront; }
  double CaloFront() const { return fCaloFront; }
  double X0Eff() const { return fX0Eff; }

  // Sampling layer of an active (LAr) volume, or -1 if the volume is passive.
  int LayerOf(const G4LogicalVolume* lv) const;

 private:
  double fUpstreamX0;
  double fPSFront = 0.0;
  double fCaloFront = 0.0;
  double fX0Eff = 0.0;
  std::map<const G4LogicalVolume*, int> fActiveLayer;
};
