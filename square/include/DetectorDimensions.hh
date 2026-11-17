#ifndef DetectorDimensions_h
#define DetectorDimensions_h 1

#include "G4SystemOfUnits.hh"

namespace DetectorDimensions {
inline const G4double PlateX = 1000.0 * mm;
inline const G4double PlateY = 1000.0 * mm;
inline const G4double PlateZ = 10.0 * mm;

inline const G4double HalfX = PlateX / 2.0;
inline const G4double HalfY = PlateY / 2.0;
inline const G4double HalfZ = PlateZ / 2.0;

inline const G4double SiPMWidth = 6.07 * mm;
inline const G4double SiPMThickness = 0.37 * mm;

// 45-degree chamfer whose diagonal optical face is 6.07 mm
inline const G4double ChamferLeg = SiPMWidth / 1.4142135623730951;
} // namespace DetectorDimensions

#endif
