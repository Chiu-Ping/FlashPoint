#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <memory>

class G4VPhysicalVolume;
class G4GenericMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  ~DetectorConstruction() override;

  G4VPhysicalVolume *Construct() override;

  // BC-408 plastic dimensions
  G4double GetPlateX() const { return fPlateX; }
  G4double GetPlateY() const { return fPlateY; }
  G4double GetPlateZ() const { return fPlateZ; }

  // SiPM active area.
  // These dimensions define the optical coupling/cut area.
  G4double GetSiPMLength() const { return fSiPMLength; }
  G4double GetSiPMHeight() const { return fSiPMHeight; }

  // onsemi MicroFJ-60035-TSV physical package dimensions.
  G4double GetSiPMPackageLength() const { return fSiPMPackageLength; }
  G4double GetSiPMPackageHeight() const { return fSiPMPackageHeight; }
  G4double GetSiPMThickness() const { return fSiPMThickness; }

  // Four corner SiPMs are always present. Edge SiPMs are distributed
  // uniformly along each of the four straight edges.
  G4int GetEdgeSiPMsPerEdge() const { return fEdgeSiPMsPerEdge; }
  G4int GetTotalSiPMs() const { return 4 + 4 * fEdgeSiPMsPerEdge + fArraySiPMsPerSide*fArraySiPMsPerSide; }

  // Nominal single-avalanche gain, used to convert detected PE to charge.
  G4double GetSiPMGain() const { return fSiPMGain; }

private:
  // BC-408 plate
  G4double fPlateX;
  G4double fPlateY;
  G4double fPlateZ;

  // SiPM active photosensitive/coupling area
  G4double fSiPMLength;
  G4double fSiPMHeight;

  // SiPM physical package
  G4double fSiPMPackageLength;
  G4double fSiPMPackageHeight;
  G4double fSiPMThickness;

  G4int fEdgeSiPMsPerEdge;
  G4int fArraySiPMsPerSide;
  G4double fSiPMGain;

  std::unique_ptr<G4GenericMessenger> fMessenger;
};

#endif
