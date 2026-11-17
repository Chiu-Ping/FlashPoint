#ifndef DetectorDimensions_h
#define DetectorDimensions_h 1

#include "G4SystemOfUnits.hh"

namespace DetectorDimensions {
// C++ defaults only. Runtime geometry is configured in the selected .mac file.
inline const G4double HexEdgeLength = 400.0 * mm;
inline const G4double PlateX = 2.0 * HexEdgeLength;
inline const G4double PlateY = 1.7320508075688772 * HexEdgeLength;
inline const G4double PlateZ = 10.0 * mm;

inline const G4double HalfX = PlateX / 2.0;
inline const G4double HalfY = PlateY / 2.0;
inline const G4double HalfZ = PlateZ / 2.0;

inline const G4double SiPMWidth = 6.07 * mm;
inline const G4double SiPMThickness = 0.37 * mm;

// Distance trimmed from either adjacent side at a 120-degree hex vertex.
inline const G4double VertexTrim = SiPMWidth / 1.7320508075688772;
} // namespace DetectorDimensions

#endif
