#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class G4Step;
class G4OpBoundaryProcess;
class PrimaryGeneratorAction;
class EventAction;

class SteppingAction : public G4UserSteppingAction {
public:
  SteppingAction(const PrimaryGeneratorAction *generator,
                 EventAction *eventAction);

  ~SteppingAction() override = default;

  void UserSteppingAction(const G4Step *) override;

private:
  void FindBoundaryProcess();

  const PrimaryGeneratorAction *fGenerator = nullptr;
  EventAction *fEventAction = nullptr;

  G4OpBoundaryProcess *fBoundary = nullptr;
};

#endif
