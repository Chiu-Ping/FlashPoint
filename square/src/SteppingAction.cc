#include "SteppingAction.hh"

#include "AnalysisIds.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "AnalysisConfig.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"

#include "G4AnalysisManager.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

SteppingAction::SteppingAction(const PrimaryGeneratorAction *generator,
                               EventAction *eventAction)
    : fGenerator(generator), fEventAction(eventAction) {}

void SteppingAction::FindBoundaryProcess() {
  auto *opticalPhoton = G4OpticalPhoton::OpticalPhotonDefinition();

  auto *processManager = opticalPhoton->GetProcessManager();

  if (!processManager)
    return;

  auto *processList = processManager->GetProcessList();

  const G4int nProcesses = processManager->GetProcessListLength();

  for (G4int i = 0; i < nProcesses; ++i) {
    auto *process = (*processList)[i];

    auto *boundary = dynamic_cast<G4OpBoundaryProcess *>(process);

    if (boundary) {
      fBoundary = boundary;
      return;
    }
  }
}

void SteppingAction::UserSteppingAction(const G4Step *step) {
  auto *track = step->GetTrack();

  if (track->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) {
    return;
  }

  if (!fBoundary) {
    FindBoundaryProcess();
  }

  if (!fBoundary)
    return;

  auto *preStep = step->GetPreStepPoint();

  auto *postStep = step->GetPostStepPoint();

  if (postStep->GetStepStatus() != fGeomBoundary) {
    return;
  }

  auto *prePV = preStep->GetPhysicalVolume();

  auto *postPV = postStep->GetPhysicalVolume();

  if (!prePV || !postPV) {
    return;
  }

  auto *preLV = prePV->GetLogicalVolume();

  auto *postLV = postPV->GetLogicalVolume();

  if (!preLV || !postLV) {
    return;
  }

  // ================================================================
  // REACHED SiPM
  //
  // Count every optical photon arriving at the BC-408 -> SiPM
  // boundary BEFORE applying the PDE/Detection requirement.
  //
  // With the current effective sensor surface:
  //   REFLECTIVITY = 0
  //
  // so a non-detected photon is absorbed and cannot be counted again.
  // ================================================================

  if (preLV->GetName() != "ScintillatorLV" || postLV->GetName() != "SiPMLV") {
    return;
  }

  const G4int sipmID = postPV->GetCopyNo();

  if (sipmID < 0 || sipmID >= fEventAction->GetNumberOfSiPMs()) {
    return;
  }

  fEventAction->AddReached(sipmID);

  // ================================================================
  // DETECTED SiPM photon
  // ================================================================

  if (fBoundary->GetStatus() != Detection) {
    return;
  }

  const G4double timeNs = track->GetGlobalTime() / ns;

  const G4double energy = track->GetTotalEnergy();

  const G4double wavelengthNm = (h_Planck * c_light / energy) / nm;

  const auto position = postStep->GetPosition();

  const G4double hitX = position.x() / mm;

  const G4double hitY = position.y() / mm;

  // creator:
  //   0 = Scintillation
  //   1 = Cerenkov
  //   2 = other
  G4int creator = 2;

  const auto *creatorProcess = track->GetCreatorProcess();

  if (creatorProcess) {
    const auto &name = creatorProcess->GetProcessName();

    if (name == "Scintillation") {
      creator = 0;
    } else if (name == "Cerenkov") {
      creator = 1;
    }
  }

  const G4double pathMm = track->GetTrackLength() / mm;

  // Update event-level summaries.
  fEventAction->AddDetected(sipmID, timeNs, pathMm, creator);

  auto *analysis = G4AnalysisManager::Instance();

  // dN/dt
  analysis->FillH1(AnalysisIds::TimeSiPM(sipmID), timeNs);

  // Photon path-length distribution
  analysis->FillH1(
      AnalysisIds::PathSiPM(sipmID, fEventAction->GetNumberOfSiPMs()), pathMm);

  // ================================================================
  // Photon-level hits ntuple
  // ================================================================

  const auto *event = G4EventManager::GetEventManager()->GetConstCurrentEvent();

  const G4int eventID = event ? event->GetEventID() : -1;

  analysis->FillNtupleIColumn(AnalysisIds::HitsNtuple, 0, eventID);

  analysis->FillNtupleIColumn(AnalysisIds::HitsNtuple, 1, sipmID);

  analysis->FillNtupleDColumn(AnalysisIds::HitsNtuple, 2, timeNs);
  if (AnalysisConfig::FullHits()) {
    analysis->FillNtupleDColumn(AnalysisIds::HitsNtuple, 3, wavelengthNm);

    analysis->FillNtupleDColumn(AnalysisIds::HitsNtuple, 4, hitX);

    analysis->FillNtupleDColumn(AnalysisIds::HitsNtuple, 5, hitY);

    analysis->FillNtupleDColumn(AnalysisIds::HitsNtuple, 6,
                                fGenerator->GetLastX() / mm);

    analysis->FillNtupleDColumn(AnalysisIds::HitsNtuple, 7,
                                fGenerator->GetLastY() / mm);

    analysis->FillNtupleIColumn(AnalysisIds::HitsNtuple, 8, creator);

    analysis->FillNtupleDColumn(AnalysisIds::HitsNtuple, 9, pathMm);
  }
  analysis->AddNtupleRow(AnalysisIds::HitsNtuple);
}
