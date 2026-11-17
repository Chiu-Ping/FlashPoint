#include "PrimaryGeneratorAction.hh"
#include "AnalysisConfig.hh"
#include "AnalysisIds.hh"
#include "DetectorConstruction.hh"

#include "G4AnalysisManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"

#include <algorithm>
#include <cmath>

PrimaryGeneratorAction::PrimaryGeneratorAction(
    const DetectorConstruction *detector)
    : fDetector(detector) {
  fParticleGun = new G4ParticleGun(1);

  auto *muon = G4ParticleTable::GetParticleTable()->FindParticle("mu-");

  fParticleGun->SetParticleDefinition(muon);
  fParticleGun->SetParticleEnergy(4.0 * GeV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., -1.));

  fMessenger =
      std::make_unique<G4GenericMessenger>(this, "/source/", "source control");
  fMessenger->DeclareMethod("position",
                            &PrimaryGeneratorAction::SetSourcePosition,
                            "Set source x y position");
  fMessenger->DeclareProperty(
      "randomXY", fRandomXY,
      "If true, sample each muon uniformly over the scintillator XY area");
  fMessenger->DeclarePropertyWithUnit("xMin", "mm", fXMin,
                                      "Minimum X for random source position");
  fMessenger->DeclarePropertyWithUnit("xMax", "mm", fXMax,
                                      "Maximum X for random source position");
  fMessenger->DeclarePropertyWithUnit("yMin", "mm", fYMin,
                                      "Minimum Y for random source position");
  fMessenger->DeclarePropertyWithUnit("yMax", "mm", fYMax,
                                      "Maximum Y for random source position");
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() { delete fParticleGun; }

G4bool PrimaryGeneratorAction::IsInsideScintillator(G4double x,
                                                    G4double y) const {
  const G4double edgeLength = fDetector->GetHexEdgeLength();
  const G4double sqrt3 = std::sqrt(3.0);
  const G4double apothem = sqrt3 * edgeLength / 2.0;
  const G4double vertexFaceRadius =
      edgeLength - fDetector->GetSiPMLength() / (2.0 * sqrt3);

  // The scintillator is the intersection of the six original hex-edge
  // half-planes and the six small vertex-cut half-planes.
  for (G4int i = 0; i < 6; ++i) {
    const G4double vertexTheta = i * twopi / 6.0;
    if (x * std::cos(vertexTheta) + y * std::sin(vertexTheta) >
        vertexFaceRadius) {
      return false;
    }

    const G4double edgeNormalTheta = (i + 0.5) * twopi / 6.0;
    if (x * std::cos(edgeNormalTheta) + y * std::sin(edgeNormalTheta) >
        apothem) {
      return false;
    }
  }

  return true;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *event) {
  const auto hx = fDetector->GetPlateX() / 2.;
  const auto hy = fDetector->GetPlateY() / 2.;

  G4double x = fSourceX;
  G4double y = fSourceY;

  if (fRandomXY) {
    const G4double xmin = std::max(fXMin, -hx);
    const G4double xmax = std::min(fXMax, hx);
    const G4double ymin = std::max(fYMin, -hy);
    const G4double ymax = std::min(fYMax, hy);
    do {
      x = xmin + (xmax - xmin) * G4UniformRand();
      y = ymin + (ymax - ymin) * G4UniformRand();
    } while (!IsInsideScintillator(x, y));
  }

  fLastX = x;
  fLastY = y;

  const auto z = fDetector->GetPlateZ() / 2. + 0.1 * mm;

  fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
  fParticleGun->GeneratePrimaryVertex(event);

  auto *ana = G4AnalysisManager::Instance();
  ana->FillH1(AnalysisIds::EventCount(fDetector->GetTotalSiPMs()), 0.5);
  if (AnalysisConfig::ResponseH2()) {
    ana->FillH2(AnalysisIds::SourceXY, x / mm, y / mm);
  }
}
