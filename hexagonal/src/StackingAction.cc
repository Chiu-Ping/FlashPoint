#include "StackingAction.hh"

#include "EventAction.hh"

#include "G4OpticalPhoton.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"

StackingAction::StackingAction(EventAction *eventAction)
    : fEventAction(eventAction) {}

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track *track) {
  if (!track ||
      track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) {
    return fUrgent;
  }

  const auto *creator = track->GetCreatorProcess();

  if (!creator) {
    // Primary optical photons are not expected in this study.
    return fUrgent;
  }

  const auto &processName = creator->GetProcessName();

  if (processName == "Scintillation") {
    fEventAction->AddProducedScint();
  } else if (processName == "Cerenkov") {
    fEventAction->AddProducedCherenkov();
  }

  return fUrgent;
}
