#include "G4Types.hh"
#ifndef AnalysisIds_h
#define AnalysisIds_h 1

namespace AnalysisIds {
// Dynamic H1 layout for N SiPMs:
//   [0,N)       arrival-time histograms
//   N           event counter
//   [N+1,2N+1) optical path histograms
inline constexpr int TimeSiPM(int sipm) { return sipm; }
inline constexpr int EventCount(int numberOfSiPMs) { return numberOfSiPMs; }
inline constexpr int PathSiPM(int sipm, int numberOfSiPMs) {
  return numberOfSiPMs + 1 + sipm;
}

// H2 histogram IDs
inline constexpr int SourceXY = 0;
inline constexpr int ResponseSiPM(int sipm) { return 1 + sipm; }

// Ntuple IDs
inline constexpr int HitsNtuple = 0;
inline constexpr int EventsNtuple = 1;

constexpr G4int ConfigNtuple = 2;
constexpr G4int GeometryNtuple = 3;
} // namespace AnalysisIds

#endif
