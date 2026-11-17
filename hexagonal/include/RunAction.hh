#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"

class G4Run;
class DetectorConstruction;
class AnalysisMessenger;

class RunAction : public G4UserRunAction {
public:
  explicit RunAction(const DetectorConstruction *detector);

  ~RunAction() override;

  void BeginOfRunAction(const G4Run *) override;

  void EndOfRunAction(const G4Run *) override;

private:
  void BookAnalysis();

  const DetectorConstruction *fDetector = nullptr;
  AnalysisMessenger *fAnalysisMessenger = nullptr;
  G4bool fAnalysisBooked = false;
};

#endif
