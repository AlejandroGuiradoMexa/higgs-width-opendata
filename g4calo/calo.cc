// Usage: calo <particle> <energy_GeV> <n_events> <upstream_X0> <output.txt> [seed] [vis.mac]
// Example: calo e- 50 2000 2.0 out/e-_50GeV_X02.txt 1
// With a visualisation macro (e.g. macros/vrml.mac) the events are drawn instead of run silently; used only for
// the 3D view on the project page, not for the samples of the analysis.
#include <cstdlib>
#include <iostream>
#include <string>

#include "Actions.hh"
#include "Detector.hh"
#include "FTFP_BERT.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "Randomize.hh"

int main(int argc, char** argv) {
  if (argc < 6) {
    std::cerr << "usage: calo <particle> <energy_GeV> <n_events> <upstream_X0> <output.txt> [seed]\n";
    return 1;
  }
  const std::string particle = argv[1];
  const double energy = std::atof(argv[2]);
  const int nEvents = std::atoi(argv[3]);
  const double upstreamX0 = std::atof(argv[4]);
  const std::string out = argv[5];
  const long seed = argc > 6 ? std::atol(argv[6]) : 12345;
  const std::string visMacro = argc > 7 ? argv[7] : "";

  G4Random::setTheSeed(seed);
  auto* rm = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto* det = new Detector(upstreamX0);
  rm->SetUserInitialization(det);
  auto* physics = new FTFP_BERT(0);
  rm->SetUserInitialization(physics);
  rm->SetUserInitialization(new Actions(det, particle, energy, out));
  rm->Initialize();
  G4UImanager::GetUIpointer()->ApplyCommand("/run/printProgress 0");
  if (visMacro.empty()) {
    rm->BeamOn(nEvents);
  } else {
    auto* vis = new G4VisExecutive("quiet");
    vis->Initialize();
    G4UImanager::GetUIpointer()->ApplyCommand("/control/execute " + visMacro);
    rm->BeamOn(nEvents);
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/flush");
    delete vis;
  }
  delete rm;
  return 0;
}
