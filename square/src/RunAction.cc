#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "AnalysisMessenger.hh"
#include "AnalysisConfig.hh"

#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <string>

RunAction::RunAction(const DetectorConstruction *detector)
    : fDetector(detector) {
  fAnalysisMessenger = new AnalysisMessenger();
  auto *analysis = G4AnalysisManager::Instance();

  analysis->SetDefaultFileType("root");

  analysis->SetFileName("sipm_timing");

  analysis->SetVerboseLevel(1);

  analysis->SetNtupleMerging(true);
}

RunAction::~RunAction()
{
    delete fAnalysisMessenger;
}

void RunAction::BookAnalysis() {
  if (fAnalysisBooked)
    return;

  auto *analysis = G4AnalysisManager::Instance();

  const G4int numberOfSiPMs = fDetector ? fDetector->GetTotalSiPMs() : 4;

  // ================================================================
  // Detected-photon arrival-time histograms
  // 0-200 ns, 50 ps/bin. The wider window retains the long reflected-light
  // tail expected from a one-metre plate.
  // ================================================================

  for (G4int i = 0; i < numberOfSiPMs; ++i) {
    analysis->CreateH1("sipm" + std::to_string(i) + "_time",
                       "SiPM " + std::to_string(i) +
                           " detected photon time;time [ns];Detected photons",
                       4000, 0.0, 200.0);
  }

  // Event counter
  analysis->CreateH1("event_counter", "Generated events", 1, 0.0, 1.0);

  // ================================================================
  // Detected-photon optical path-length distributions
  //
  // Range 0-10 m retains long paths in the multiply-reflecting 1-metre plate.
  // ================================================================

  for (G4int i = 0; i < numberOfSiPMs; ++i) {
    analysis->CreateH1("sipm" + std::to_string(i) + "_path",
                       "SiPM " + std::to_string(i) +
                           " detected photon path;path [mm];Detected photons",
                       2000, 0.0, 10000.0);
  }

  // ================================================================
  // XY response maps
  // ================================================================

  const G4double halfX =
      fDetector ? fDetector->GetPlateX() / (2.0 * mm) : 105.0;

  const G4double halfY =
      fDetector ? fDetector->GetPlateY() / (2.0 * mm) : 105.0;

  const G4int binsX = std::max(10, static_cast<G4int>(2.0 * halfX));

  const G4int binsY = std::max(10, static_cast<G4int>(2.0 * halfY));

  if (AnalysisConfig::ResponseH2()) {
    // One count per generated muon.
    analysis->CreateH2("source_xy", "Muon generation position;x [mm];y [mm]",
                      binsX, -halfX, halfX, binsY, -halfY, halfY);

    // Weighted sum of Ndetected/event in each source-XY bin.
    for (G4int i = 0; i < numberOfSiPMs; ++i) {
      analysis->CreateH2("response_sipm" + std::to_string(i) + "_xy",
                        "SiPM " + std::to_string(i) +
                            " detected response;x_{#mu} [mm];y_{#mu} [mm]",
                        binsX, -halfX, halfX, binsY, -halfY, halfY);
    }
  }

  // ================================================================
  // Ntuple 0: photon-level detected hits
  // ================================================================

  analysis->CreateNtuple("hits", "Detected SiPM optical photons");

  analysis->CreateNtupleIColumn("event");

  analysis->CreateNtupleIColumn("sipm");

  analysis->CreateNtupleDColumn("t_ns");

  if (AnalysisConfig::FullHits()) {
    analysis->CreateNtupleDColumn("wavelength_nm");

    analysis->CreateNtupleDColumn("hit_x_mm");

    analysis->CreateNtupleDColumn("hit_y_mm");

    analysis->CreateNtupleDColumn("source_x_mm");

    analysis->CreateNtupleDColumn("source_y_mm");

    analysis->CreateNtupleIColumn("creator");

    analysis->CreateNtupleDColumn("path_mm");
  }

  analysis->FinishNtuple();

  // ================================================================
  // Ntuple 1: one row per muon event
  //
  // This provides explicit photon accounting and timing summaries.
  // ================================================================

  analysis->CreateNtuple("events",
                         "Per-event photon accounting and SiPM response");

  analysis->CreateNtupleIColumn("event");

  analysis->CreateNtupleDColumn("source_x_mm");

  analysis->CreateNtupleDColumn("source_y_mm");

  analysis->CreateNtupleIColumn("n_produced");

  analysis->CreateNtupleIColumn("n_scint_produced");

  analysis->CreateNtupleIColumn("n_cherenkov_produced");

  for (G4int i = 0; i < numberOfSiPMs; ++i) {
    analysis->CreateNtupleIColumn("n_reach_sipm" + std::to_string(i));
  }

  for (G4int i = 0; i < numberOfSiPMs; ++i) {
    analysis->CreateNtupleIColumn("n_detected_sipm" + std::to_string(i));
  }

  for (G4int i = 0; i < numberOfSiPMs; ++i) {
    analysis->CreateNtupleIColumn("n_scint_detected_sipm" + std::to_string(i));
  }

  for (G4int i = 0; i < numberOfSiPMs; ++i) {
    analysis->CreateNtupleIColumn("n_cherenkov_detected_sipm" +
                                  std::to_string(i));
  }

  for (G4int i = 0; i < numberOfSiPMs; ++i) {
    analysis->CreateNtupleDColumn("first_time_sipm" + std::to_string(i) +
                                  "_ns");
  }

  for (G4int i = 0; i < numberOfSiPMs; ++i) {
    analysis->CreateNtupleDColumn("mean_time_sipm" + std::to_string(i) + "_ns");
  }

  for (G4int i = 0; i < numberOfSiPMs; ++i) {
    analysis->CreateNtupleDColumn("mean_path_sipm" + std::to_string(i) + "_mm");
  }

  // Integrated avalanche charge from N_PE * nominal gain * electron charge.
  // A voltage amplitude additionally requires a front-end/shaper model.
  for (G4int i = 0; i < numberOfSiPMs; ++i) {
    analysis->CreateNtupleDColumn("charge_sipm" + std::to_string(i) + "_pc");
  }

  analysis->FinishNtuple();

  fAnalysisBooked = true;
}

void RunAction::BeginOfRunAction(const G4Run *) {
  // Book only now, after all pre-initialization macro commands have set the
  // final detector and SiPM count. This keeps master and worker ntuple schemas
  // identical when edgeSiPMsPerEdge differs from the C++ default.
  BookAnalysis();

  auto *analysis = G4AnalysisManager::Instance();

  analysis->OpenFile();
}

void RunAction::EndOfRunAction(const G4Run *) {
  auto *analysis = G4AnalysisManager::Instance();

  analysis->Write();
  analysis->CloseFile();
}
