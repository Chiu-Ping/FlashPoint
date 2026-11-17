#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "AnalysisConfig.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <vector>

#include <memory>

class G4VPhysicalVolume;
class G4GenericMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  ~DetectorConstruction() override;

  G4VPhysicalVolume *Construct() override;

  // BC-408 regular-hexagon dimensions. The hexagon is flat-topped, with
  // vertices on the +/-X axis. GetPlateX/Y return its bounding-box spans.
  G4double GetHexEdgeLength() const { return fHexEdgeLength; }
  G4double GetPlateX() const { return 2.0 * fHexEdgeLength; }
  G4double GetPlateY() const { return std::sqrt(3.0) * fHexEdgeLength; }
  G4double GetPlateZ() const { return fThickness; }

  // SiPM active area.
  // These dimensions define the optical coupling/cut area.
  G4double GetSiPMLength() const { return fSiPMLength; }
  G4double GetSiPMHeight() const { return fSiPMHeight; }

  // onsemi MicroFJ-60035-TSV physical package dimensions.
  G4double GetSiPMPackageLength() const { return fSiPMPackageLength; }
  G4double GetSiPMPackageHeight() const { return fSiPMPackageHeight; }
  G4double GetSiPMThickness() const { return fSiPMThickness; }

  // Six vertex SiPMs are always present. Edge SiPMs are distributed
  // uniformly along each of the six straight edges.
  G4int GetEdgeSiPMsPerEdge() const { return fEdgeSiPMsPerEdge; }
  G4int GetCentralSiPMDivisions() const { return fCentralSiPMDivisions; }
  G4int GetCentralSiPMCount() const {
    if (fCentralSiPMMode == "ring") {
      const G4int r = fCentralSiPMRings;
      if (r < 0)
        return 0;
      return 1 + 3 * r * (r + 1);
    }

    if (fCentralSiPMMode == "division") {
      const G4int n = fCentralSiPMDivisions;
      if (n < 2)
        return 0;

      // 3 diagonals, each with n-1 interior points.
      // If n is even, center appears on all 3 diagonals
      // but must only be counted once.
      return 3 * (n - 1) - (n % 2 == 0 ? 2 : 0);
    }

    return 0;
  }
  G4double GetThickness() const { return fThickness; }
  const std::vector<G4ThreeVector> &GetSiPMPositions() const {
    return fSiPMPositions;
  }

  G4int GetTotalSiPMs() const {
    return 6 + 6 * fEdgeSiPMsPerEdge + GetCentralSiPMCount();
  }

  // Nominal single-avalanche gain, used to convert detected PE to charge.
  G4double GetSiPMGain() const { return fSiPMGain; }

private:
  // BC-408 regular hexagon
  G4double fHexEdgeLength;
  G4double fThickness;

  // SiPM active photosensitive/coupling area
  G4double fSiPMLength;
  G4double fSiPMHeight;

  // SiPM physical package
  G4double fSiPMPackageLength;
  G4double fSiPMPackageHeight;
  G4double fSiPMThickness;

  G4int fEdgeSiPMsPerEdge;
  G4int fCentralSiPMDivisions = 0;
  G4double fSiPMGain;
  std::vector<G4ThreeVector> fSiPMPositions;
  std::unique_ptr<G4GenericMessenger> fMessenger;

  // Setup SiPMs at the cernter region
  G4String fCentralSiPMMode = "division";
  G4int fCentralSiPMRings = 3;
  G4double fCentralSiPMPitch = 150. * mm;
  std::vector<G4ThreeVector> MakeCentralRingPositions(G4int rings,
                                                      G4double pitch) const;
};

#endif
