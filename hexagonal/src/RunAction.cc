#include "RunAction.hh"

#include "AnalysisIds.hh"
#include "AnalysisMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Threading.hh"

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

RunAction::~RunAction() { delete fAnalysisMessenger; }

void RunAction::BookAnalysis() {
  if (fAnalysisBooked)
    return;

  auto *analysis = G4AnalysisManager::Instance();

  const G4int numberOfSiPMs = fDetector ? fDetector->GetTotalSiPMs() : 12;

  // ================================================================
  // Detected-photon arrival-time histograms
  // 0-200 ns, 50 ps/bin. The wider window retains the long reflected-light
  // tail expected from a large hexagonal plate.
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
  // Range 0-10 m retains long paths in the multiply-reflecting hexagonal plate.
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

  // ================================================================
  // Ntuple 2: detector configuration — one row per run
  // ================================================================
  analysis->CreateNtuple("detector_config", "Detector configuration");
  analysis->CreateNtupleIColumn("edge_sipms_per_edge");    // 0
  analysis->CreateNtupleIColumn("central_sipm_divisions"); // 1
  analysis->CreateNtupleIColumn("n_vertex_sipms");         // 2
  analysis->CreateNtupleIColumn("n_edge_sipms");           // 3
  analysis->CreateNtupleIColumn("n_central_sipms");        // 4
  analysis->CreateNtupleIColumn("n_total_sipms");          // 5
  analysis->CreateNtupleDColumn("hex_edge_length_mm");     // 6
  analysis->CreateNtupleDColumn("thickness_mm");           // 7
  analysis->CreateNtupleDColumn("sipm_gain");              // 8
  analysis->FinishNtuple();

  // ================================================================
  // Ntuple 3: physical SiPM geometry — one row per SiPM
  //
  // region:
  //   0 = vertex
  //   1 = edge
  //   2 = central +Z
  // ================================================================
  analysis->CreateNtuple("sipm_geometry", "SiPM positions");
  analysis->CreateNtupleIColumn("sipm");   // 0
  analysis->CreateNtupleIColumn("region"); // 1
  analysis->CreateNtupleDColumn("x_mm");   // 2
  analysis->CreateNtupleDColumn("y_mm");   // 3
  analysis->CreateNtupleDColumn("z_mm");   // 4
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
  const G4int tid = G4Threading::G4GetThreadId();

  // Print only once: master thread
  if (G4Threading::IsMasterThread()) {
    G4cout << "\n========== Analysis Output ==========\n"
           << " Mode        : "
           << (AnalysisConfig::FullHits() ? "FULL" : "LIGHT") << "\n"
           << " Hit columns : "
           << (AnalysisConfig::FullHits() ? "all photon information"
                                          : "event, sipm, t_ns")
           << "\n"
           << " Response H2 : "
           << (AnalysisConfig::ResponseH2() ? "enabled" : "disabled") << "\n"
           << "=====================================\n"
           << G4endl;
  }

#ifdef G4MULTITHREADED
  if (tid != 0 || !fDetector)
    return; // metadata written by worker 0 only
#else
  if (!fDetector)
    return;
#endif

  G4cout << "[Metadata] thread=" << tid
         << " totalSiPMs=" << fDetector->GetTotalSiPMs()
         << " positions=" << fDetector->GetSiPMPositions().size() << G4endl;
  const G4int nEdgePerEdge = fDetector->GetEdgeSiPMsPerEdge();
  const G4int nEdge = 6 * nEdgePerEdge;
  const G4int nCentral = fDetector->GetCentralSiPMCount();
  const G4int nTotal = fDetector->GetTotalSiPMs();

  // ================================================================
  // detector_config
  // ================================================================
  analysis->FillNtupleIColumn(AnalysisIds::ConfigNtuple, 0, nEdgePerEdge);
  analysis->FillNtupleIColumn(AnalysisIds::ConfigNtuple, 1,
                              fDetector->GetCentralSiPMDivisions());
  analysis->FillNtupleIColumn(AnalysisIds::ConfigNtuple, 2, 6);
  analysis->FillNtupleIColumn(AnalysisIds::ConfigNtuple, 3, nEdge);
  analysis->FillNtupleIColumn(AnalysisIds::ConfigNtuple, 4, nCentral);
  analysis->FillNtupleIColumn(AnalysisIds::ConfigNtuple, 5, nTotal);
  analysis->FillNtupleDColumn(AnalysisIds::ConfigNtuple, 6,
                              fDetector->GetHexEdgeLength() / mm);
  analysis->FillNtupleDColumn(AnalysisIds::ConfigNtuple, 7,
                              fDetector->GetThickness() / mm);
  analysis->FillNtupleDColumn(AnalysisIds::ConfigNtuple, 8,
                              fDetector->GetSiPMGain());

  analysis->AddNtupleRow(AnalysisIds::ConfigNtuple);

  // ================================================================
  // sipm_geometry
  // ================================================================
  const auto &positions = fDetector->GetSiPMPositions();
  for (G4int i = 0; i < static_cast<G4int>(positions.size()); ++i) {
    G4int region = 2;

    if (i < 6)
      region = 0; // vertex
    else if (i < 6 + nEdge)
      region = 1; // edge
    else
      region = 2; // central +Z

    analysis->FillNtupleIColumn(AnalysisIds::GeometryNtuple, 0, i);
    analysis->FillNtupleIColumn(AnalysisIds::GeometryNtuple, 1, region);
    analysis->FillNtupleDColumn(AnalysisIds::GeometryNtuple, 2,
                                positions[i].x() / mm);
    analysis->FillNtupleDColumn(AnalysisIds::GeometryNtuple, 3,
                                positions[i].y() / mm);
    analysis->FillNtupleDColumn(AnalysisIds::GeometryNtuple, 4,
                                positions[i].z() / mm);
    analysis->AddNtupleRow(AnalysisIds::GeometryNtuple);
  }
}

void RunAction::EndOfRunAction(const G4Run *) {
  auto *analysis = G4AnalysisManager::Instance();

  analysis->Write();
  analysis->CloseFile();
}
