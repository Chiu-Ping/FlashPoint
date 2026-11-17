#include "DetectorConstruction.hh"
#include "AnalysisConfig.hh"

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
#include <memory>
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
    : fHexEdgeLength(400.0 * mm), fThickness(10.0 * mm),

      // onsemi MicroFJ-60035-TSV effective photosensitive area
      fSiPMLength(6.07 * mm), fSiPMHeight(6.07 * mm),

      // onsemi MicroFJ-60035-TSV package: 6.13 x 6.13 x 0.37 mm3.
      fSiPMPackageLength(6.13 * mm), fSiPMPackageHeight(6.13 * mm),
      fSiPMThickness(0.37 * mm), fEdgeSiPMsPerEdge(1), fSiPMGain(6.3e6) {
  fMessenger = std::make_unique<G4GenericMessenger>(
      this, "/detector/", "Plastic scintillator / SiPM geometry");

  // BC-408 regular hexagon
  fMessenger->DeclarePropertyWithUnit(
      "edgeLength", "mm", fHexEdgeLength,
      "Nominal edge length of the regular hexagonal scintillator");

  fMessenger->DeclarePropertyWithUnit("thickness", "mm", fThickness,
                                      "BC-408 scintillator thickness");

  // SiPM active/coupling area
  fMessenger->DeclarePropertyWithUnit(
      "sipmLength", "mm", fSiPMLength,
      "SiPM active length along a vertex or edge contact face");

  fMessenger->DeclarePropertyWithUnit(
      "sipmHeight", "mm", fSiPMHeight,
      "SiPM active height along Z; must not exceed plastic thickness");

  // SiPM physical package
  fMessenger->DeclarePropertyWithUnit(
      "sipmPackageLength", "mm", fSiPMPackageLength,
      "Physical SiPM package length along the contact face");

  fMessenger->DeclarePropertyWithUnit("sipmPackageHeight", "mm",
                                      fSiPMPackageHeight,
                                      "Physical SiPM package height along Z");

  fMessenger->DeclarePropertyWithUnit(
      "sipmThickness", "mm", fSiPMThickness,
      "Physical SiPM package thickness normal to the contact face");

  fMessenger->DeclareProperty(
      "edgeSiPMsPerEdge", fEdgeSiPMsPerEdge,
      "Number of evenly spaced SiPMs on each of the six edges (>= 0)");

  fMessenger->DeclareProperty(
      "centralSiPMDivisions", fCentralSiPMDivisions,
      "Subdivide each of the 3 diagonals into N parts and place "
      "SiPMs at k/N (k=1,...,N-1); 0 disables central SiPMs");

  fMessenger->DeclareProperty("centralSiPMMode", fCentralSiPMMode,
                              "Central SiPM mode: division or ring");

  fMessenger->DeclareProperty("centralSiPMRings", fCentralSiPMRings,
                              "Number of triangular-lattice rings");

  fMessenger->DeclarePropertyWithUnit("centralSiPMPitch", "mm",
                                      fCentralSiPMPitch,
                                      "Pitch between central SiPMs");

  fMessenger->DeclareProperty(
      "sipmGain", fSiPMGain,
      "Nominal avalanche gain used for the per-event charge estimate");
}

DetectorConstruction::~DetectorConstruction() = default;

std::vector<G4ThreeVector>
DetectorConstruction::MakeCentralRingPositions(G4int rings,
                                               G4double pitch) const {
  std::vector<G4ThreeVector> pos;

  for (G4int q = -rings; q <= rings; ++q) {
    const G4int rmin = std::max(-rings, -q - rings);
    const G4int rmax = std::min(rings, -q + rings);

    for (G4int r = rmin; r <= rmax; ++r) {
      const G4double x = pitch * (q + 0.5 * r);
      const G4double y = pitch * std::sqrt(3.) / 2. * r;

      pos.emplace_back(x, y, 0.);
    }
  }

  return pos;
}

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
  const G4double sqrt3 = std::sqrt(3.0);
  const G4double hx = fHexEdgeLength;
  const G4double hy = sqrt3 * fHexEdgeLength / 2.0;
  const G4double hz = fThickness / 2.0;

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

  if (fHexEdgeLength <= 0.0 || fThickness <= 0.0 || sipmLength <= 0.0 ||
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
  if (sipmHeight > fThickness + 1.e-9 * mm) {
    G4ExceptionDescription msg;

    msg << "The SiPM active height must not exceed the plastic thickness.\n"
        << "Plastic Z     = " << fThickness / mm << " mm\n"
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
  // BC-408 regular hexagon with all six vertices truncated for SiPM contact.
  //
  // Each original hex vertex has a 120-degree interior angle. Trimming the
  // same distance d along both adjoining edges creates a contact face of
  // length sqrt(3)*d, so d = activeLength/sqrt(3). The nominal edgeLength
  // remains the side length of the untrimmed regular hexagon.
  // ================================================================

  const G4double vertexTrim = sipmLength / sqrt3;
  const G4double straightEdgeLength = fHexEdgeLength - 2.0 * vertexTrim;

  if (straightEdgeLength <= 0.0 ||
      straightEdgeLength + 1.e-9 * mm < sipmPackageLength) {
    G4Exception(
        "DetectorConstruction::Construct", "Geometry004", FatalException,
        "The hexagon edge is too short for the vertex cuts and edge SiPM.");
  }

  std::array<G4TwoVector, 6> vertices;
  for (G4int i = 0; i < 6; ++i) {
    const G4double theta = i * twopi / 6.0;
    vertices[i] = G4TwoVector(fHexEdgeLength * std::cos(theta),
                              fHexEdgeLength * std::sin(theta));
  }

  std::vector<G4TwoVector> polygon;
  polygon.reserve(12);
  const G4double trimFraction = vertexTrim / fHexEdgeLength;

  for (G4int i = 0; i < 6; ++i) {
    const auto &vertex = vertices[i];
    const auto &previous = vertices[(i + 5) % 6];
    const auto &next = vertices[(i + 1) % 6];

    polygon.emplace_back(
        vertex.x() + trimFraction * (previous.x() - vertex.x()),
        vertex.y() + trimFraction * (previous.y() - vertex.y()));
    polygon.emplace_back(vertex.x() + trimFraction * (next.x() - vertex.x()),
                         vertex.y() + trimFraction * (next.y() - vertex.y()));
  }

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

  // Fixed vertex SiPM IDs 0-5, counter-clockwise from the +X vertex.
  const G4double vertexFaceRadius = fHexEdgeLength - sipmLength / (2.0 * sqrt3);

  for (G4int i = 0; i < 6; ++i) {
    const G4double theta = i * twopi / 6.0;
    const G4double nx = std::cos(theta);
    const G4double ny = std::sin(theta);
    placeSiPM(vertexFaceRadius * nx, vertexFaceRadius * ny, nx, ny);
  }

  // Flexible straight-edge SiPMs. For N sensors, centre j is placed at
  // fraction (j+1)/(N+1) of the remaining straight section. Edge groups
  // follow the counter-clockwise edges V0->V1 through V5->V0.
  if (fEdgeSiPMsPerEdge > 0) {
    const G4double pitch = straightEdgeLength / (fEdgeSiPMsPerEdge + 1);

    if (pitch + 1.e-9 * mm < sipmPackageLength) {
      G4Exception("DetectorConstruction::Construct", "Geometry006",
                  FatalException,
                  "Too many edge SiPMs: packages would overlap each other or a "
                  "vertex cut.");
    }

    for (G4int edge = 0; edge < 6; ++edge) {
      const auto &vertex = vertices[edge];
      const auto &next = vertices[(edge + 1) % 6];

      const G4double startX =
          vertex.x() + trimFraction * (next.x() - vertex.x());
      const G4double startY =
          vertex.y() + trimFraction * (next.y() - vertex.y());
      const G4double endX = next.x() + trimFraction * (vertex.x() - next.x());
      const G4double endY = next.y() + trimFraction * (vertex.y() - next.y());

      // The outward normal lies midway between the two vertex directions.
      const G4double normalTheta = (edge + 0.5) * twopi / 6.0;
      const G4double nx = std::cos(normalTheta);
      const G4double ny = std::sin(normalTheta);

      for (G4int j = 0; j < fEdgeSiPMsPerEdge; ++j) {
        const G4double fraction =
            static_cast<G4double>(j + 1) / (fEdgeSiPMsPerEdge + 1);
        placeSiPM(startX + fraction * (endX - startX),
                  startY + fraction * (endY - startY), nx, ny);
      }
    }
  }

  // ================================================================
  // Helper for SiPMs mounted on the +Z face
  // ================================================================
  auto *topRotation = new G4RotationMatrix();
  topRotation->rotateX(90.0 * deg);

  auto placeTopSiPM = [&](G4double x, G4double y) {
    const G4int sipmID = static_cast<G4int>(sipmPV.size());

    sipmPV.push_back(new G4PVPlacement(
        topRotation, G4ThreeVector(x, y, hz + sipmThickness / 2.0), sipmLV,
        "SiPMPV", worldLV, false, sipmID, true));
  };

  // ================================================================
  // Central SiPM layout
  // ================================================================
  if (fCentralSiPMMode == "division") {

    // --------------------------------------------------------------
    // Existing three-diagonal method
    // --------------------------------------------------------------
    if (fCentralSiPMDivisions >= 2) {

      const G4int n = fCentralSiPMDivisions;
      G4bool centerPlaced = false;

      for (G4int diag = 0; diag < 3; ++diag) {

        const auto &a = vertices[diag];
        const auto &b = vertices[diag + 3];

        for (G4int k = 1; k < n; ++k) {

          // For even n, k=n/2 is the common center
          // of all three diagonals.
          const G4bool isCenter = (n % 2 == 0 && k == n / 2);

          if (isCenter && centerPlaced)
            continue;

          const G4double f = static_cast<G4double>(k) / n;
          const G4double x = a.x() + f * (b.x() - a.x());
          const G4double y = a.y() + f * (b.y() - a.y());
          placeTopSiPM(x, y);

          if (isCenter)
            centerPlaced = true;
        }
      }
    }

  } else if (fCentralSiPMMode == "ring") {

    // --------------------------------------------------------------
    // 2-D triangular / hexagonal lattice
    // --------------------------------------------------------------
    const auto positions =
        MakeCentralRingPositions(fCentralSiPMRings, fCentralSiPMPitch);

    for (const auto &p : positions)
      placeTopSiPM(p.x(), p.y());

  } else {

    G4ExceptionDescription msg;
    msg << "Unknown centralSiPMMode: " << fCentralSiPMMode
        << "\nAllowed values: division, ring";

    G4Exception("DetectorConstruction::Construct", "GeometryCentral001",
                FatalException, msg);
  }

  fSiPMPositions.clear();
  fSiPMPositions.reserve(sipmPV.size());

  for (const auto *pv : sipmPV)
    fSiPMPositions.push_back(pv->GetObjectTranslation());

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
         << " FlashPointHexShape detector\n"
         << "------------------------------------------------------------\n"
         << " Scintillator : BC-408\n"
         << " Shape        : regular hexagon (flat-topped)\n"
         << " Edge length  : " << fHexEdgeLength / mm << " mm\n"
         << " XY span      : " << 2.0 * hx / mm << " x " << 2.0 * hy / mm
         << " mm2\n"
         << " Thickness    : " << fThickness / mm << " mm\n"
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
         << "\n"
         << " SiPM geometry\n"
         << "------------------------------------------------------------\n"
         << " Vertex SiPMs : 6\n"
         << " Edge add-ons : " << 6 * fEdgeSiPMsPerEdge << " ("
         << fEdgeSiPMsPerEdge << " per edge)\n"
         << " Central mode : " << fCentralSiPMMode << "\n";

  if (fCentralSiPMMode == "division") {

    G4cout << " Central Divs : " << fCentralSiPMDivisions << "\n";

  } else if (fCentralSiPMMode == "ring") {

    G4cout << " Central Rings: " << fCentralSiPMRings << "\n"
           << " Ring pitch   : " << fCentralSiPMPitch / mm << " mm\n";
  }

  G4cout << " Central SiPMs: " << GetCentralSiPMCount() << "\n"
         << " Total SiPMs  : " << sipmPV.size() << "\n"
         << " Vertex face  : " << sipmLength / mm << " mm\n"
         << " Vertex trim  : " << vertexTrim / mm << " mm per adjacent edge\n"
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
