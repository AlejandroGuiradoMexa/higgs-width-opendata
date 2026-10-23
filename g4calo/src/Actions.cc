#include "Actions.hh"

#include <cmath>

#include "Detector.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"

PrimaryGenerator::PrimaryGenerator(const std::string& particle, double energyGeV)
    : fGun(new G4ParticleGun(1)) {
  fGun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle(particle));
  fGun->SetParticleEnergy(energyGeV * GeV);
  fGun->SetParticleMomentumDirection({0, 0, 1});
  fGun->SetParticlePosition({0, 0, -geo::kGunDistance * mm});
}

PrimaryGenerator::~PrimaryGenerator() { delete fGun; }

void PrimaryGenerator::GeneratePrimaries(G4Event* event) { fGun->GeneratePrimaryVertex(event); }

EventRecorder::EventRecorder(const std::string& outFile) : fOut(outFile) {
  // Header: number of cells per layer, then the depth bins.
  fOut << "#";
  for (int l = 0; l < geo::kNLayers; ++l) fOut << " L" << l << ":" << geo::kNx[l] << "x" << geo::kNy[l];
  fOut << " depth:" << geo::kNDepthBins << " active passive (MeV)\n";
}

void EventRecorder::BeginOfEventAction(const G4Event*) {
  for (int l = 0; l < geo::kNLayers; ++l) fCells[l].assign(geo::kNx[l] * geo::kNy[l], 0.0);
  fDepth.fill(0.0);
  fActive = fPassive = 0.0;
}

void EventRecorder::AddCell(int layer, int ix, int iy, double edep) {
  if (ix < 0 || iy < 0 || ix >= geo::kNx[layer] || iy >= geo::kNy[layer]) return;
  fCells[layer][ix * geo::kNy[layer] + iy] += edep;
}

void EventRecorder::AddDepth(int bin, double edep) {
  if (bin >= 0 && bin < geo::kNDepthBins) fDepth[bin] += edep;
}

void EventRecorder::EndOfEventAction(const G4Event*) {
  for (const auto& layer : fCells)
    for (double e : layer) fOut << e / MeV << ' ';
  for (double e : fDepth) fOut << e / MeV << ' ';
  fOut << fActive / MeV << ' ' << fPassive / MeV << '\n';
}

void Stepping::UserSteppingAction(const G4Step* step) {
  const double edep = step->GetTotalEnergyDeposit();
  if (edep <= 0) return;
  const auto* lv = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  const int layer = fDet->LayerOf(lv);
  if (layer < 0) {
    fRec->AddPassive(edep);
    return;
  }
  fRec->AddActive(edep);

  // Cells are centred on the impact point (x = y = 0).
  const auto pos = 0.5 * (step->GetPreStepPoint()->GetPosition() + step->GetPostStepPoint()->GetPosition());
  const int ix = static_cast<int>(std::floor(pos.x() / geo::kCellX[layer] + 0.5)) + geo::kNx[layer] / 2;
  const int iy = static_cast<int>(std::floor(pos.y() / geo::kCellY[layer] + 0.5)) + geo::kNy[layer] / 2;
  fRec->AddCell(layer, ix, iy, edep);
  if (layer >= geo::kS1 && layer <= geo::kS3) {
    fRec->AddDepth(static_cast<int>((pos.z() - fDet->CaloFront()) / fDet->X0Eff()), edep);
  }
}

void Actions::Build() const {
  SetUserAction(new PrimaryGenerator(fParticle, fEnergy));
  auto* rec = new EventRecorder(fOut);
  SetUserAction(rec);
  SetUserAction(new Stepping(fDet, rec));
}
