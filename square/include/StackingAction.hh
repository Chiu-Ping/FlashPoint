#ifndef StackingAction_h
#define StackingAction_h 1

#include "G4UserStackingAction.hh"

class G4Track;
class EventAction;

class StackingAction : public G4UserStackingAction {
public:
  explicit StackingAction(EventAction *eventAction);

  ~StackingAction() override = default;

  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track *) override;

private:
  EventAction *fEventAction = nullptr;
};

#endif
