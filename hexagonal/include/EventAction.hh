#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <vector>

class G4Event;
class PrimaryGeneratorAction;

class EventAction : public G4UserEventAction {
public:
  EventAction(const PrimaryGeneratorAction *generator, G4int numberOfSiPMs,
              G4double sipmGain);

  ~EventAction() override = default;

  void BeginOfEventAction(const G4Event *) override;

  void EndOfEventAction(const G4Event *) override;

  void AddProducedScint();
  void AddProducedCherenkov();

  void AddReached(G4int sipm);

  G4int GetNumberOfSiPMs() const { return fNumberOfSiPMs; }

  // creator:
  //   0 = Scintillation
  //   1 = Cerenkov
  //   2 = other
  void AddDetected(G4int sipm, G4double timeNs, G4double pathMm, G4int creator);

private:
  const PrimaryGeneratorAction *fGenerator = nullptr;
  G4int fNumberOfSiPMs = 0;
  G4double fSiPMGain = 0.0;

  G4int fProducedScint = 0;
  G4int fProducedCherenkov = 0;

  std::vector<G4int> fReached;
  std::vector<G4int> fDetected;

  std::vector<G4int> fDetectedScint;
  std::vector<G4int> fDetectedCherenkov;

  std::vector<G4double> fFirstTime;
  std::vector<G4double> fTimeSum;
  std::vector<G4double> fPathSum;
};

#endif
