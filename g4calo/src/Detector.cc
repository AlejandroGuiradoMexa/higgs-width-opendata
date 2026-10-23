#include "Detector.hh"

#include <string>

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

namespace {

// Place `n` periods of absorber + active slabs starting at z; returns the z after the last one.
double PlacePeriods(G4LogicalVolume* world, G4LogicalVolume* absorber, G4LogicalVolume* active,
                    double tAbs, double tAct, int n, double z, int& copy) {
  for (int i = 0; i < n; ++i, ++copy) {
    new G4PVPlacement(nullptr, {0, 0, z + 0.5 * tAbs}, absorber, absorber->GetName(), world, false, copy, false);
    z += tAbs;
    new G4PVPlacement(nullptr, {0, 0, z + 0.5 * tAct}, active, active->GetName(), world, false, copy, false);
    z += tAct;
  }
  return z;
}

}  // namespace

G4VPhysicalVolume* Detector::Construct() {
  auto* nist = G4NistManager::Instance();
  auto* air = nist->FindOrBuildMaterial("G4_AIR");
  auto* pb = nist->FindOrBuildMaterial("G4_Pb");
  auto* lar = nist->FindOrBuildMaterial("G4_lAr");
  auto* al = nist->FindOrBuildMaterial("G4_Al");
  auto* fe = nist->FindOrBuildMaterial("G4_Fe");
  auto* scint = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  const double halfXY = 0.6 * m;
  auto* world = new G4LogicalVolume(new G4Box("World", 1.0 * m, 1.0 * m, 2.5 * m), air, "World");
  auto* worldPV = new G4PVPlacement(nullptr, {}, world, "World", nullptr, false, 0, true);

  // Effective radiation length of one Pb/LAr period.
  const double period = (geo::kPb + geo::kLAr) * mm;
  fX0Eff = period / (geo::kPb * mm / pb->GetRadlen() + geo::kLAr * mm / lar->GetRadlen());

  // Upstream material ends 50 mm in front of the presampler, whose front face is at z = 0.
  fPSFront = 0.0;
  if (fUpstreamX0 > 0) {
    const double t = fUpstreamX0 * al->GetRadlen();
    auto* lv = new G4LogicalVolume(new G4Box("Upstream", halfXY, halfXY, 0.5 * t), al, "Upstream");
    new G4PVPlacement(nullptr, {0, 0, -50 * mm - 0.5 * t}, lv, "Upstream", world, false, 0, true);
  }

  // Presampler: a single LAr gap.
  double z = fPSFront;
  {
    const double t = geo::kPresamplerLAr * mm;
    auto* lv = new G4LogicalVolume(new G4Box("PS", halfXY, halfXY, 0.5 * t), lar, "PS");
    new G4PVPlacement(nullptr, {0, 0, z + 0.5 * t}, lv, "PS", world, false, 0, true);
    fActiveLayer[lv] = geo::kPS;
    z += t;
  }
  z += 30 * mm;

  // Accordion layers as planar Pb/LAr periods.
  fCaloFront = z;
  auto* pbLV = new G4LogicalVolume(new G4Box("Pb", halfXY, halfXY, 0.5 * geo::kPb * mm), pb, "Pb");
  auto* larBox = new G4Box("LAr", halfXY, halfXY, 0.5 * geo::kLAr * mm);
  int copy = 0;
  for (int layer = geo::kS1; layer <= geo::kS3; ++layer) {
    auto* larLV = new G4LogicalVolume(larBox, lar, "LAr_" + std::to_string(layer));
    fActiveLayer[larLV] = layer;
    z = PlacePeriods(world, pbLV, larLV, geo::kPb * mm, geo::kLAr * mm, geo::kPeriods[layer], z, copy);
  }

  // Hadronic calorimeter (Fe/scintillator) after a 50 mm gap for the cryostat.
  z += 50 * mm;
  auto* feLV = new G4LogicalVolume(new G4Box("Fe", halfXY, halfXY, 0.5 * geo::kFe * mm), fe, "Fe");
  auto* scLV = new G4LogicalVolume(new G4Box("Scint", halfXY, halfXY, 0.5 * geo::kScint * mm), scint, "Scint");
  fActiveLayer[scLV] = geo::kHad;
  PlacePeriods(world, feLV, scLV, geo::kFe * mm, geo::kScint * mm, geo::kPeriods[geo::kHad], z, copy);
  return worldPV;
}

int Detector::LayerOf(const G4LogicalVolume* lv) const {
  auto it = fActiveLayer.find(lv);
  return it == fActiveLayer.end() ? -1 : it->second;
}
