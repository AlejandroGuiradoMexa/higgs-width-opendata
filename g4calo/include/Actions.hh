#pragma once

#include <array>
#include <fstream>
#include <string>
#include <vector>

#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "Geometry.hh"

class Detector;
class G4ParticleGun;

class PrimaryGenerator : public G4VUserPrimaryGeneratorAction {
 public:
  PrimaryGenerator(const std::string& particle, double energyGeV);
  ~PrimaryGenerator() override;
  void GeneratePrimaries(G4Event* event) override;

 private:
  G4ParticleGun* fGun;
};

// Collects energy in the readout windows of each layer and writes one line per event.
class EventRecorder : public G4UserEventAction {
 public:
  explicit EventRecorder(const std::string& outFile);
  void BeginOfEventAction(const G4Event*) override;
  void EndOfEventAction(const G4Event* event) override;

  void AddCell(int layer, int ix, int iy, double edep);
  void AddDepth(int bin, double edep);
  void AddActive(double edep) { fActive += edep; }
  void AddPassive(double edep) { fPassive += edep; }

 private:
  std::ofstream fOut;
  std::array<std::vector<double>, geo::kNLayers> fCells;
  std::array<double, geo::kNDepthBins> fDepth{};
  double fActive = 0.0;
  double fPassive = 0.0;
};

class Stepping : public G4UserSteppingAction {
 public:
  Stepping(const Detector* det, EventRecorder* rec) : fDet(det), fRec(rec) {}
  void UserSteppingAction(const G4Step* step) override;

 private:
  const Detector* fDet;
  EventRecorder* fRec;
};

class Actions : public G4VUserActionInitialization {
 public:
  Actions(const Detector* det, std::string particle, double energyGeV, std::string outFile)
      : fDet(det), fParticle(std::move(particle)), fEnergy(energyGeV), fOut(std::move(outFile)) {}
  void Build() const override;

 private:
  const Detector* fDet;
  std::string fParticle;
  double fEnergy;
  std::string fOut;
};
