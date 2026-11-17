#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GenericMessenger.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;
class DetectorConstruction;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  explicit PrimaryGeneratorAction(const DetectorConstruction *);
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event *) override;
  void SetSourcePosition(G4double x, G4double y){
      fSourceX = x;
      fSourceY = y;
  }

  G4double GetLastX() const { return fLastX; }
  G4double GetLastY() const { return fLastY; }

private:
  G4bool IsInsideScintillator(G4double, G4double) const;
  std::unique_ptr<G4GenericMessenger> fMessenger;

  G4ParticleGun *fParticleGun = nullptr;
  const DetectorConstruction *fDetector = nullptr;

  G4double fLastX = 0.;
  G4double fLastY = 0.;

  G4double fSourceX = 0.;
  G4double fSourceY = 0.;
  // random XY range
  G4bool fRandomXY = false;
  G4double fXMin = -500.0;
  G4double fXMax =  500.0;
  G4double fYMin = -500.0;
  G4double fYMax =  500.0;
};

#endif
