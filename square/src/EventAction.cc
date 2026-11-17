#include "EventAction.hh"

#include "AnalysisIds.hh"
#include "PrimaryGeneratorAction.hh"
#include "AnalysisConfig.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>

EventAction::EventAction(const PrimaryGeneratorAction *generator,
                         G4int numberOfSiPMs, G4double sipmGain)
    : fGenerator(generator), fNumberOfSiPMs(numberOfSiPMs),
      fSiPMGain(sipmGain), fReached(numberOfSiPMs), fDetected(numberOfSiPMs),
      fDetectedScint(numberOfSiPMs), fDetectedCherenkov(numberOfSiPMs),
      fFirstTime(numberOfSiPMs), fTimeSum(numberOfSiPMs),
      fPathSum(numberOfSiPMs) {}

void EventAction::BeginOfEventAction(const G4Event *) {
  fProducedScint = 0;
  fProducedCherenkov = 0;

  std::fill(fReached.begin(), fReached.end(), 0);
  std::fill(fDetected.begin(), fDetected.end(), 0);

  std::fill(fDetectedScint.begin(), fDetectedScint.end(), 0);
  std::fill(fDetectedCherenkov.begin(), fDetectedCherenkov.end(), 0);

  // -1 means no detected photon in this SiPM for this event.
  std::fill(fFirstTime.begin(), fFirstTime.end(), -1.0);

  std::fill(fTimeSum.begin(), fTimeSum.end(), 0.0);
  std::fill(fPathSum.begin(), fPathSum.end(), 0.0);
}

void EventAction::AddProducedScint() { ++fProducedScint; }

void EventAction::AddProducedCherenkov() { ++fProducedCherenkov; }

void EventAction::AddReached(G4int sipm) {
  if (sipm < 0 || sipm >= fNumberOfSiPMs)
    return;

  ++fReached[sipm];
}

void EventAction::AddDetected(G4int sipm, G4double timeNs, G4double pathMm,
                              G4int creator) {
  if (sipm < 0 || sipm >= fNumberOfSiPMs)
    return;

  ++fDetected[sipm];

  fTimeSum[sipm] += timeNs;
  fPathSum[sipm] += pathMm;

  if (fFirstTime[sipm] < 0.0 || timeNs < fFirstTime[sipm]) {
    fFirstTime[sipm] = timeNs;
  }

  if (creator == 0) {
    ++fDetectedScint[sipm];
  } else if (creator == 1) {
    ++fDetectedCherenkov[sipm];
  }
}

void EventAction::EndOfEventAction(const G4Event *event) {
  auto *analysis = G4AnalysisManager::Instance();

  const G4int eventID = event ? event->GetEventID() : -1;

  const G4double sourceX = fGenerator ? fGenerator->GetLastX() / mm : 0.0;

  const G4double sourceY = fGenerator ? fGenerator->GetLastY() / mm : 0.0;

  const G4int nProduced = fProducedScint + fProducedCherenkov;

  // ================================================================
  // Event-level ntuple
  //
  // Column layout:
  //
  //  0 event
  //  1 source_x_mm
  //  2 source_y_mm
  //
  //  3 n_produced
  //  4 n_scint_produced
  //  5 n_cherenkov_produced
  //
  // The following blocks each contain N columns, where N is the runtime
  // number of SiPMs: reached, detected, scintillation detected, Cherenkov
  // detected, first time, mean time, mean path, and integrated charge.
  // ================================================================

  analysis->FillNtupleIColumn(AnalysisIds::EventsNtuple, 0, eventID);
  analysis->FillNtupleDColumn(AnalysisIds::EventsNtuple, 1, sourceX);
  analysis->FillNtupleDColumn(AnalysisIds::EventsNtuple, 2, sourceY);
  analysis->FillNtupleIColumn(AnalysisIds::EventsNtuple, 3, nProduced);
  analysis->FillNtupleIColumn(AnalysisIds::EventsNtuple, 4, fProducedScint);
  analysis->FillNtupleIColumn(AnalysisIds::EventsNtuple, 5, fProducedCherenkov);

  const G4int reachedOffset = 6;
  const G4int detectedOffset = reachedOffset + fNumberOfSiPMs;
  const G4int scintOffset = detectedOffset + fNumberOfSiPMs;
  const G4int cherenkovOffset = scintOffset + fNumberOfSiPMs;
  const G4int firstTimeOffset = cherenkovOffset + fNumberOfSiPMs;
  const G4int meanTimeOffset = firstTimeOffset + fNumberOfSiPMs;
  const G4int meanPathOffset = meanTimeOffset + fNumberOfSiPMs;
  const G4int chargeOffset = meanPathOffset + fNumberOfSiPMs;

  for (G4int i = 0; i < fNumberOfSiPMs; ++i) {
    const G4double meanTime =
        fDetected[i] > 0 ? fTimeSum[i] / fDetected[i] : -1.0;

    const G4double meanPath =
        fDetected[i] > 0 ? fPathSum[i] / fDetected[i] : -1.0;

    analysis->FillNtupleIColumn(AnalysisIds::EventsNtuple, reachedOffset + i,
                                fReached[i]);

    analysis->FillNtupleIColumn(AnalysisIds::EventsNtuple, detectedOffset + i,
                                fDetected[i]);

    analysis->FillNtupleIColumn(AnalysisIds::EventsNtuple, scintOffset + i,
                                fDetectedScint[i]);

    analysis->FillNtupleIColumn(AnalysisIds::EventsNtuple, cherenkovOffset + i,
                                fDetectedCherenkov[i]);

    analysis->FillNtupleDColumn(AnalysisIds::EventsNtuple, firstTimeOffset + i,
                                fFirstTime[i]);

    analysis->FillNtupleDColumn(AnalysisIds::EventsNtuple, meanTimeOffset + i,
                                meanTime);

    analysis->FillNtupleDColumn(AnalysisIds::EventsNtuple, meanPathOffset + i,
                                meanPath);

    // Q = N_PE * gain * elementary charge. This is integrated sensor
    // charge before front-end shaping, not a voltage peak.
    const G4double chargePc =
        static_cast<G4double>(fDetected[i]) * fSiPMGain * eplus / coulomb *
        1.0e12;
    analysis->FillNtupleDColumn(AnalysisIds::EventsNtuple, chargeOffset + i,
                                chargePc);

    // Weighted response map:
    // sum of detected PE in each source-XY bin.
    //
    // Divide this histogram by source_xy in analysis to obtain
    // mean detected PE/event versus source position.
    if (AnalysisConfig::ResponseH2()) {
      analysis->FillH2(AnalysisIds::ResponseSiPM(i), sourceX, sourceY,
                       static_cast<G4double>(fDetected[i]));
    }
  }

  analysis->AddNtupleRow(AnalysisIds::EventsNtuple);
}
