#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4ExtrudedSolid.hh"
#include "G4TwoVector.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4Exception.hh"
#include "G4GenericMessenger.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace {
// Convert wavelength in nm to Geant4 photon energy.
// Input wavelength arrays below are ordered from long to short
// wavelength so the resulting energy arrays are monotonically
// increasing, as required by G4MaterialPropertiesTable.
std::vector<G4double>
PhotonEnergiesFromWavelengths(const std::vector<G4double> &wavelengthNm) {
  std::vector<G4double> energies;
  energies.reserve(wavelengthNm.size());

  for (const auto wavelength : wavelengthNm) {
    energies.push_back((h_Planck * c_light) / (wavelength * nm));
  }

  return energies;
}
} // namespace

DetectorConstruction::DetectorConstruction()
    : fPlateX(1000.0 * mm), fPlateY(1000.0 * mm), fPlateZ(10.0 * mm),

      // onsemi MicroFJ-60035-TSV effective photosensitive area
      fSiPMLength(6.07 * mm), fSiPMHeight(6.07 * mm),

      // onsemi MicroFJ-60035-TSV package: 6.13 x 6.13 x 0.37 mm3.
      fSiPMPackageLength(6.13 * mm), fSiPMPackageHeight(6.13 * mm),
      fSiPMThickness(0.37 * mm), fEdgeSiPMsPerEdge(1), fArraySiPMsPerSide(2),
      fSiPMGain(6.3e6) {
  fMessenger = std::make_unique<G4GenericMessenger>(
      this, "/detector/", "Plastic scintillator / SiPM geometry");

  // BC-408 plate
  fMessenger->DeclarePropertyWithUnit("sizeX", "mm", fPlateX,
                                      "BC-408 plastic size in X");

  fMessenger->DeclarePropertyWithUnit("sizeY", "mm", fPlateY,
                                      "BC-408 plastic size in Y");

  fMessenger->DeclarePropertyWithUnit("sizeZ", "mm", fPlateZ,
                                      "BC-408 plastic thickness");

  // SiPM active/coupling area
  fMessenger->DeclarePropertyWithUnit(
      "sipmLength", "mm", fSiPMLength,
      "SiPM active length along the 45-degree cut plane");

  fMessenger->DeclarePropertyWithUnit(
      "sipmHeight", "mm", fSiPMHeight,
      "SiPM active height along Z; normally equal to plastic thickness");

  // SiPM physical package
  fMessenger->DeclarePropertyWithUnit(
      "sipmPackageLength", "mm", fSiPMPackageLength,
      "Physical SiPM package length along the cut plane");

  fMessenger->DeclarePropertyWithUnit("sipmPackageHeight", "mm",
                                      fSiPMPackageHeight,
                                      "Physical SiPM package height along Z");

  // Keep the original command name for compatibility with existing macros.
  fMessenger->DeclarePropertyWithUnit(
      "sipmThickness", "mm", fSiPMThickness,
      "Physical SiPM package thickness normal to the cut plane");

  fMessenger->DeclareProperty(
      "edgeSiPMsPerEdge", fEdgeSiPMsPerEdge,
      "Number of evenly spaced SiPMs on each straight edge (>= 0)");

  fMessenger->DeclareProperty("arraySiPMsPerSide", fArraySiPMsPerSide,
                              "Central SiPM array size per side");

  fMessenger->DeclareProperty(
      "sipmGain", fSiPMGain,
      "Nominal avalanche gain used for the per-event charge estimate");
}

DetectorConstruction::~DetectorConstruction() = default;

G4VPhysicalVolume *DetectorConstruction::Construct() {
  auto *nist = G4NistManager::Instance();

  // ================================================================
  // Materials
  // ================================================================

  auto *air = nist->FindOrBuildMaterial("G4_AIR");

  // BC-408 is a polyvinyltoluene (PVT)-based plastic scintillator.
  // Use the Geant4 PVT composition but set the manufacturer density
  // explicitly to 1.023 g/cm3.
  auto *pvt = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  auto *plastic = new G4Material("BC408", 1.023 * g / cm3, 1);

  plastic->AddMaterial(pvt, 1.0);

  // Effective bulk material used for the MicroFJ-60035-TSV package.
  // Photon detection itself is implemented at the optical boundary
  // using the manufacturer's PDE curve.
  auto *silicon = nist->FindOrBuildMaterial("G4_Si");

  // ================================================================
  // Broad optical energy grid
  //
  // Covers the MicroFJ-60035-TSV spectral response (200-900 nm).
  // BC-408 RINDEX and bulk attenuation are approximated as constant
  // across this range because the manufacturer gives scalar values.
  // ================================================================

  const std::vector<G4double> opticalWavelengthNm = {
      900.0, 800.0, 700.0, 600.0, 520.0, 500.0, 480.0, 460.0, 450.0,
      440.0, 430.0, 425.0, 420.0, 410.0, 400.0, 380.0, 360.0, 340.0,
      320.0, 300.0, 280.0, 270.0, 250.0, 220.0, 200.0};

  const auto opticalEnergy = PhotonEnergiesFromWavelengths(opticalWavelengthNm);

  // ================================================================
  // Air optical properties
  // ================================================================

  {
    std::vector<G4double> rindex(opticalEnergy.size(), 1.0003);

    auto *mpt = new G4MaterialPropertiesTable();

    mpt->AddProperty("RINDEX", opticalEnergy, rindex);

    air->SetMaterialPropertiesTable(mpt);
  }

  // ================================================================
  // BC-408 optical properties
  //
  // Luxium BC-408 manufacturer values used here:
  //
  //   Base                         PVT
  //   Density                      1.023 g/cm3
  //   Refractive index             1.58
  //   Rise time                    0.9 ns
  //   Decay time                   2.1 ns
  //   Peak emission                425 nm
  //   Bulk attenuation length      380 cm
  //   Effective attenuation length 210 cm (finished long-sheet test)
  //   Relative light output        64% anthracene
  //
  // For Geant4 ABSLENGTH we use the intrinsic/bulk value 3.80 m.
  //
  // The manufacturer does not quote an absolute photons/MeV number
  // in the current BC-408 datasheet. 10000 photons/MeV is retained
  // as the conventional engineering approximation for BC-408/PVT.
  //
  // The emission array below is an approximate digitization of the
  // manufacturer's BC-408 emission-spectrum plot, normalized to 1.
  // ================================================================

  {
    std::vector<G4double> rindex(opticalEnergy.size(), 1.58);

    std::vector<G4double> absorption(opticalEnergy.size(), 3.80 * m);

    auto *mpt = new G4MaterialPropertiesTable();

    mpt->AddProperty("RINDEX", opticalEnergy, rindex);

    mpt->AddProperty("ABSLENGTH", opticalEnergy, absorption);

    // Wavelengths are descending so energies are ascending.
    const std::vector<G4double> emissionWavelengthNm = {
        520.0, 500.0, 480.0, 470.0, 460.0, 450.0, 440.0, 435.0, 430.0,
        425.0, 420.0, 410.0, 400.0, 390.0, 380.0, 370.0, 360.0};

    const auto emissionEnergy =
        PhotonEnergiesFromWavelengths(emissionWavelengthNm);

    const std::vector<G4double> scintillationSpectrum = {
        0.02, // 520 nm
        0.07, // 500
        0.18, // 480
        0.28, // 470
        0.40, // 460
        0.52, // 450
        0.68, // 440
        0.84, // 435
        0.97, // 430
        1.00, // 425
        0.90, // 420
        0.58, // 410
        0.30, // 400
        0.16, // 390
        0.09, // 380
        0.05, // 370
        0.02  // 360
    };

    mpt->AddProperty("SCINTILLATIONCOMPONENT1", emissionEnergy,
                     scintillationSpectrum);
    mpt->AddConstProperty("SCINTILLATIONYIELD", 10000.0 / MeV);
    mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    mpt->AddConstProperty("SCINTILLATIONRISETIME1", 0.9 * ns);
    mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.1 * ns);
    mpt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    plastic->SetMaterialPropertiesTable(mpt);
  }

  // ================================================================
  // Detector dimensions
  // ================================================================
  const G4double hx = fPlateX / 2.0;
  const G4double hy = fPlateY / 2.0;
  const G4double hz = fPlateZ / 2.0;

  // Active/coupling area.
  const G4double sipmLength = fSiPMLength;
  const G4double sipmHeight = fSiPMHeight;

  // Physical package.
  const G4double sipmPackageLength = fSiPMPackageLength;
  const G4double sipmPackageHeight = fSiPMPackageHeight;
  const G4double sipmThickness = fSiPMThickness;

  // ================================================================
  // Geometry checks
  // ================================================================

  if (fPlateX <= 0.0 || fPlateY <= 0.0 || fPlateZ <= 0.0 || sipmLength <= 0.0 ||
      sipmHeight <= 0.0 || sipmPackageLength <= 0.0 ||
      sipmPackageHeight <= 0.0 || sipmThickness <= 0.0) {
    G4Exception("DetectorConstruction::Construct", "Geometry001",
                FatalException, "All detector dimensions must be positive.");
  }

  if (fEdgeSiPMsPerEdge < 0 || fSiPMGain <= 0.0) {
    G4ExceptionDescription msg;

    msg << "edgeSiPMsPerEdge must be non-negative and sipmGain must be "
           "positive.\n"
        << "edgeSiPMsPerEdge = " << fEdgeSiPMsPerEdge << "\n"
        << "sipmGain         = " << fSiPMGain;

    G4Exception("DetectorConstruction::Construct", "Geometry002",
                FatalException, msg);
  }

  // The 6.07-mm-high sensor is centred on the 10-mm-thick plate edge.
  // Requiring the active height to fit prevents silently losing active area.
  if (sipmHeight > fPlateZ + 1.e-9 * mm) {
    G4ExceptionDescription msg;

    msg << "The SiPM active height must not exceed the plastic thickness.\n"
        << "Plastic Z     = " << fPlateZ / mm << " mm\n"
        << "SiPM active Z = " << sipmHeight / mm << " mm";

    G4Exception("DetectorConstruction::Construct", "Geometry002a",
                FatalException, msg);
  }

  if (sipmPackageLength + 1.e-9 * mm < sipmLength ||
      sipmPackageHeight + 1.e-9 * mm < sipmHeight) {
    G4ExceptionDescription msg;

    msg << "The physical SiPM package must be at least as large "
        << "as the active/coupling area.\n"
        << "Active area  = " << sipmLength / mm << " x " << sipmHeight / mm
        << " mm2\n"
        << "Package face = " << sipmPackageLength / mm << " x "
        << sipmPackageHeight / mm << " mm2";

    G4Exception("DetectorConstruction::Construct", "Geometry003",
                FatalException, msg);
  }

  // ================================================================
  // World
  //
  // Follow the runtime detector dimensions automatically.
  // ================================================================

  const G4double worldMargin = 250.0 * mm;

  auto *worldSolid =
      new G4Box("WorldSolid", hx + worldMargin, hy + worldMargin,
                std::max(hz, sipmPackageHeight / 2.0) + worldMargin);

  auto *worldLV = new G4LogicalVolume(worldSolid, air, "WorldLV");

  auto *worldPV = new G4PVPlacement(nullptr, G4ThreeVector(), worldLV,
                                    "WorldPV", nullptr, false, 0, true);

  // ================================================================
  // BC-408 plate with four 45-degree corner cuts
  //
  // The diagonal cut length equals the SiPM ACTIVE length.
  //
  // For a 45-degree cut:
  //
  //   chamferLeg = activeLength / sqrt(2)
  //
  // For the MicroFJ-60035-TSV default, activeLength = 6.07 mm and
  // chamferLeg = 4.292 mm. The package is centred vertically on the
  // 10-mm plate edge, so only its 6.07-mm active height is sensitive.
  // ================================================================

  const G4double chamferLeg = sipmLength / std::sqrt(2.0);

  if (2.0 * chamferLeg >= fPlateX || 2.0 * chamferLeg >= fPlateY) {
    G4Exception(
        "DetectorConstruction::Construct", "Geometry004", FatalException,
        "SiPM active length is too large for the requested plate size.");
  }

  std::vector<G4TwoVector> polygon;
  polygon.emplace_back(-hx + chamferLeg, -hy);
  polygon.emplace_back(hx - chamferLeg, -hy);
  polygon.emplace_back(hx, -hy + chamferLeg);
  polygon.emplace_back(hx, hy - chamferLeg);
  polygon.emplace_back(hx - chamferLeg, hy);
  polygon.emplace_back(-hx + chamferLeg, hy);
  polygon.emplace_back(-hx, hy - chamferLeg);
  polygon.emplace_back(-hx, -hy + chamferLeg);

  auto *scintSolid = new G4ExtrudedSolid(
      "ScintillatorSolid", polygon, hz, G4TwoVector(), 1.0, G4TwoVector(), 1.0);
  auto *scintLV = new G4LogicalVolume(scintSolid, plastic, "ScintillatorLV");
  auto *scintPV = new G4PVPlacement(nullptr, G4ThreeVector(), scintLV,
                                    "ScintillatorPV", worldLV, false, 0, true);

  // ================================================================
  // onsemi MicroFJ-60035-TSV package geometry
  //
  // Default product dimensions:
  //
  //   active photosensitive area : 6.07 x 6.07 mm2
  //   package face               : 6.13 x 6.13 mm2
  //   package thickness          : 0.37 mm typ.
  //   glass refractive index     : 1.53 at 436 nm
  //
  // This implementation treats the SiPM package as one effective
  // volume and applies the PDE directly at the BC-408/SiPM contact
  // boundary. Because the plastic cut is only 6.07 mm long, only
  // that central part of the 6.13-mm package physically
  // contacts the scintillator.
  //
  // local X = package length along cut plane
  // local Y = package thickness, outward normal
  // local Z = package height
  // ================================================================

  auto *sipmSolid = new G4Box("SiPMSolid", sipmPackageLength / 2.0,
                              sipmThickness / 2.0, sipmPackageHeight / 2.0);

  auto *sipmLV = new G4LogicalVolume(sipmSolid, silicon, "SiPMLV");

  // Numbering:
  //
  // 1 (-X,+Y)          0 (+X,+Y)
  //
  // 2 (-X,-Y)          3 (+X,-Y)
  //
  std::array<G4int, 4> signX = {1, -1, -1, 1};
  std::array<G4int, 4> signY = {1, 1, -1, -1};

  std::vector<G4VPhysicalVolume *> sipmPV;
  sipmPV.reserve(GetTotalSiPMs());

  // ================================================================
  // Mount each package on its corresponding cut plane.
  //
  // Keep the orientation convention from the validated geometry:
  //
  //   local X = along the visible SiPM face
  //   local Y = package thickness
  //
  // The package centre is shifted outward by half of its physical
  // thickness, so its front surface is coincident with the cut plane.
  // ================================================================

  auto placeSiPM = [&](G4double faceCenterX, G4double faceCenterY, G4double nx,
                       G4double ny) {
    const G4int sipmID = static_cast<G4int>(sipmPV.size());

    // With the G4 placement convention, this aligns local Y (the package
    // thickness axis) with the outward surface normal and local X along
    // the scintillator edge.
    const G4double phi = std::atan2(nx, ny);

    auto *rotation = new G4RotationMatrix();
    rotation->rotateZ(phi);

    const G4double sipmCenterX = faceCenterX + nx * sipmThickness / 2.0;
    const G4double sipmCenterY = faceCenterY + ny * sipmThickness / 2.0;

    sipmPV.push_back(new G4PVPlacement(
        rotation, G4ThreeVector(sipmCenterX, sipmCenterY, 0.0), sipmLV,
        "SiPMPV", worldLV, false, sipmID, true));
  };

  // Fixed corner SiPM IDs 0-3.
  for (G4int i = 0; i < 4; ++i) {
    const G4double sx = static_cast<G4double>(signX[i]);
    const G4double sy = static_cast<G4double>(signY[i]);
    const G4double nx = sx / std::sqrt(2.0);
    const G4double ny = sy / std::sqrt(2.0);

    // Centre of the diagonal plastic cut face.
    const G4double faceCenterX = sx * (hx - chamferLeg / 2.0);
    const G4double faceCenterY = sy * (hy - chamferLeg / 2.0);
    placeSiPM(faceCenterX, faceCenterY, nx, ny);
  }

  // Flexible straight-edge SiPMs. For N sensors, centre i is placed at
  // fraction (i+1)/(N+1) of the full plate side: 1/2 for N=1, 1/3 and
  // 2/3 for N=2, etc. Numbering groups are +Y, +X, -Y, -X.
  if (fEdgeSiPMsPerEdge > 0) {
    const G4double pitchX = fPlateX / (fEdgeSiPMsPerEdge + 1);
    const G4double pitchY = fPlateY / (fEdgeSiPMsPerEdge + 1);
    const G4double endClearance = chamferLeg + sipmPackageLength / 2.0;

    if (pitchX + 1.e-9 * mm < sipmPackageLength ||
        pitchY + 1.e-9 * mm < sipmPackageLength ||
        pitchX + 1.e-9 * mm < endClearance ||
        pitchY + 1.e-9 * mm < endClearance) {
      G4Exception("DetectorConstruction::Construct", "Geometry006",
                  FatalException,
                  "Too many edge SiPMs: packages would overlap each other or a "
                  "chamfer.");
    }

    auto coordinate = [&](G4int index, G4double halfLength) {
      return -halfLength +
             2.0 * halfLength * (index + 1) / (fEdgeSiPMsPerEdge + 1);
    };

    // +Y edge, left to right.
    for (G4int i = 0; i < fEdgeSiPMsPerEdge; ++i) {
      const G4double u = coordinate(i, hx);
      placeSiPM(u, hy, 0.0, 1.0);
    }

    // +X edge, top to bottom.
    for (G4int i = 0; i < fEdgeSiPMsPerEdge; ++i) {
      const G4double u = coordinate(i, hy);
      placeSiPM(hx, -u, 1.0, 0.0);
    }

    // -Y edge, right to left.
    for (G4int i = 0; i < fEdgeSiPMsPerEdge; ++i) {
      const G4double u = coordinate(i, hx);
      placeSiPM(-u, -hy, 0.0, -1.0);
    }

    // -X edge, bottom to top.
    for (G4int i = 0; i < fEdgeSiPMsPerEdge; ++i) {
      const G4double u = coordinate(i, hy);
      placeSiPM(-hx, u, -1.0, 0.0);
    }
  }

  // ================================================================
  // Central SiPM array on +Z face
  // ================================================================
  if (fArraySiPMsPerSide > 0) {

    const G4int n = fArraySiPMsPerSide;
    const G4double dx = fPlateX / (n + 1);
    const G4double dy = fPlateY / (n + 1);

    auto *rotTop = new G4RotationMatrix();
    rotTop->rotateX(90. * deg);

    const G4double z = fPlateZ / 2. + sipmThickness / 2.;

    for (G4int iy = 0; iy < n; ++iy) {
      for (G4int ix = 0; ix < n; ++ix) {

        const G4double x = -fPlateX / 2. + (ix + 1) * dx;
        const G4double y = -fPlateY / 2. + (iy + 1) * dy;

        const G4int id = sipmPV.size();

        sipmPV.push_back(new G4PVPlacement(rotTop, G4ThreeVector(x, y, z),
                                           sipmLV, // SAME logical volume
                                           "SiPMPV", worldLV, false, id, true));
      }
    }
  }

  // ================================================================
  // Reflective wrapping
  //
  // ESR-like high-specular-reflectivity approximation.
  //
  // 3M ESR is specified at >98% reflectivity in the visible range.
  // A constant 0.98 is used here as a conservative baseline.
  //
  // This is still represented as a zero-thickness optical surface.
  // ================================================================

  auto *reflectorSurface = new G4OpticalSurface("ReflectorSurface");
  reflectorSurface->SetType(dielectric_metal);
  reflectorSurface->SetModel(unified);
  reflectorSurface->SetFinish(polished);

  {
    std::vector<G4double> reflectivity(opticalEnergy.size(), 0.98);
    std::vector<G4double> efficiency(opticalEnergy.size(), 0.0);

    auto *surfaceMPT = new G4MaterialPropertiesTable();
    surfaceMPT->AddProperty("REFLECTIVITY", opticalEnergy, reflectivity);
    surfaceMPT->AddProperty("EFFICIENCY", opticalEnergy, efficiency);
    reflectorSurface->SetMaterialPropertiesTable(surfaceMPT);
  }

  new G4LogicalBorderSurface("ScintillatorReflectiveWrapping", scintPV, worldPV,
                             reflectorSurface);

  // ================================================================
  // onsemi MicroFJ-60035-TSV effective PDE surface
  //
  // Manufacturer:
  //
  //   spectral response : 200-900 nm (PDE > 2% at +6 V)
  //   peak wavelength   : 420 nm
  //   PDE at peak       : 50% typ. at +6 V overvoltage
  //
  // The values below are an approximate digitization of the
  // manufacturer's typical PDE-vs-wavelength graph.
  //
  // PDE excludes optical crosstalk and afterpulsing, consistent
  // with the datasheet definition.
  //
  // REFLECTIVITY is set to zero because no product-specific
  // wavelength-dependent front-surface reflectivity is supplied
  // here. Therefore non-detected photons are absorbed at this
  // effective sensor boundary.
  // ================================================================

  auto *sipmSurface = new G4OpticalSurface("SiPMSurface");
  sipmSurface->SetType(dielectric_metal);
  sipmSurface->SetModel(unified);
  sipmSurface->SetFinish(polished);

  {
    // Same descending-wavelength / ascending-energy grid used
    // for opticalEnergy:
    //
    // 900, 800, 700, 600, 520, 500, 480, 460, 450, 440,
    // 430, 425, 420, 410, 400, 380, 360, 340, 320, 300,
    // 280, 270, 250, 220, 200 nm
    //
    const std::vector<G4double> efficiency = {
        0.015, // 900 nm
        0.045, // 800
        0.105, // 700
        0.200, // 600
        0.330, // 520
        0.390, // 500
        0.440, // 480
        0.470, // 460
        0.485, // 450
        0.495, // 440
        0.500, // 430
        0.500, // 425
        0.500, // 420
        0.490, // 410
        0.470, // 400
        0.440, // 380
        0.380, // 360
        0.350, // 340
        0.290, // 320
        0.220, // 300
        0.100, // 280
        0.070, // 270
        0.045, // 250
        0.050, // 220
        0.020  // 200
    };

    std::vector<G4double> reflectivity(opticalEnergy.size(), 0.0);
    auto *sipmMPT = new G4MaterialPropertiesTable();
    sipmMPT->AddProperty("REFLECTIVITY", opticalEnergy, reflectivity);
    sipmMPT->AddProperty("EFFICIENCY", opticalEnergy, efficiency);
    sipmSurface->SetMaterialPropertiesTable(sipmMPT);
  }

  // Plastic -> SiPM effective optical detection boundaries.
  for (G4int i = 0; i < static_cast<G4int>(sipmPV.size()); ++i) {
    const G4String surfaceName = "ScintillatorToSiPM_" + std::to_string(i);

    new G4LogicalBorderSurface(surfaceName, scintPV, sipmPV[i], sipmSurface);
  }

  // ================================================================
  // Geometry / product summary
  // ================================================================

  G4cout << "\n============================================================\n"
         << " PlasticMuonTiming detector\n"
         << "------------------------------------------------------------\n"
         << " Scintillator : BC-408\n"
         << " Plate        : " << fPlateX / mm << " x " << fPlateY / mm << " x "
         << fPlateZ / mm << " mm3\n"
         << " n            : 1.58\n"
         << " Bulk att.    : 380 cm\n"
         << " Rise / decay : 0.9 / 2.1 ns\n"
         << " Emission max : 425 nm\n"
         << "\n"
         << " SiPM         : onsemi MicroFJ-60035-TSV\n"
         << " Active area  : " << sipmLength / mm << " x " << sipmHeight / mm
         << " mm2\n"
         << " Package      : " << sipmPackageLength / mm << " x "
         << sipmPackageHeight / mm << " x " << sipmThickness / mm << " mm3\n"
         << " PDE peak     : 50% typ. @ 420 nm, Vover=6 V\n"
         << " Gain         : " << fSiPMGain << "\n"
         << " Corner SiPMs : 4\n"
         << " Edge SiPMs   : " << fEdgeSiPMsPerEdge << " per edge\n"
         << " Total SiPMs  : " << sipmPV.size() << "\n"
         << " Cut length   : " << sipmLength / mm << " mm\n"
         << " Chamfer leg  : " << chamferLeg / mm << " mm\n"
         << "============================================================\n"
         << G4endl;

  // ================================================================
  // Visualization
  // ================================================================

  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto *scintVis = new G4VisAttributes(G4Colour(0.20, 0.80, 0.90, 0.45));
  scintVis->SetForceSolid(true);
  scintLV->SetVisAttributes(scintVis);

  auto *sipmVis = new G4VisAttributes(G4Colour(0.90, 0.10, 0.10));
  sipmVis->SetForceSolid(true);
  sipmLV->SetVisAttributes(sipmVis);

  return worldPV;
}
