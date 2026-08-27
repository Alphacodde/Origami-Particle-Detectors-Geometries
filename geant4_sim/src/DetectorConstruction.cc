#include "DetectorConstruction.hh"
#include "DetectorConstants.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4GDMLParser.hh"
#include <algorithm>
#include "SensorSD.hh"
#include "CADMesh.hh"
G4String DetectorConstruction::fsStructureTag = "unknown";
G4double DetectorConstruction::fsFoldState = 1.0;
G4ThreeVector DetectorConstruction::fsCentroidOffset = G4ThreeVector(0., 0., 0.);
namespace {
void SanityCheckThickness(G4VSolid* solid, G4double expectedThickness_mm,
                           const G4String& label)
{
  G4VoxelLimits noLimits;
  G4AffineTransform noTransform;
  G4double xmin, xmax, ymin, ymax, zmin, zmax;
  solid->CalculateExtent(kXAxis, noLimits, noTransform, xmin, xmax);
  solid->CalculateExtent(kYAxis, noLimits, noTransform, ymin, ymax);
  solid->CalculateExtent(kZAxis, noLimits, noTransform, zmin, zmax);
  G4double minExtent = std::min({xmax - xmin, ymax - ymin, zmax - zmin});
  G4cout << "[DetectorConstruction] " << label
         << " bounding box (mm): x=[" << xmin << "," << xmax
         << "] y=[" << ymin << "," << ymax
         << "] z=[" << zmin << "," << zmax
         << "]  smallest extent=" << minExtent
         << "  expected thickness=" << expectedThickness_mm << G4endl;
  if (minExtent < expectedThickness_mm / 3.0 ||
      minExtent > expectedThickness_mm * 3.0) {
    G4cerr << "[DetectorConstruction] WARNING: " << label
           << "'s smallest bounding-box dimension (" << minExtent
           << " mm) is far from the expected thickness ("
           << expectedThickness_mm << " mm). This STL may be an un-"
           << "thickened shell (thicken_mesh.py not run), the wrong file, "
           << "or in the wrong units. Verify before trusting any X0 output "
           << "from this run." << G4endl;
  }
}
G4ThreeVector ComputeBoundingBoxCenter(G4VSolid* solid, const G4String& label)
{
  G4VoxelLimits noLimits;
  G4AffineTransform noTransform;
  G4double xmin, xmax, ymin, ymax, zmin, zmax;
  solid->CalculateExtent(kXAxis, noLimits, noTransform, xmin, xmax);
  solid->CalculateExtent(kYAxis, noLimits, noTransform, ymin, ymax);
  solid->CalculateExtent(kZAxis, noLimits, noTransform, zmin, zmax);
  G4ThreeVector center(0.5 * (xmin + xmax),
                        0.5 * (ymin + ymax),
                        0.5 * (zmin + zmax));
  G4cout << "[DetectorConstruction] " << label
         << " bounding-box center (mm, local frame): ("
         << center.x() << ", " << center.y() << ", " << center.z()
         << ")" << G4endl;
  return center;
}
}
DetectorConstruction::DetectorConstruction()
{
  DefineMessenger();
}
DetectorConstruction::~DetectorConstruction() = default;
void DetectorConstruction::DefineMessenger()
{
  fMessenger = std::make_unique<G4GenericMessenger>(
      this, "/origami/", "Origami detector geometry controls");
  fMessenger->DeclareMethod("stlFile", &DetectorConstruction::SetSiliconSTLPath)
      .SetGuidance("Path to the SILICON solid GDML file (must already be "
                    "thickened+unioned - see diff_geom_macros/_solid_export.py "
                    "thicken_triangles_union()). Command name kept as 'stlFile' "
                    "for backward compatibility with existing macros even "
                    "though the expected file is now .gdml, not .stl - see "
                    "file header, THIRD CHANGE.");
  fMessenger->DeclareMethod("kaptonStlFile", &DetectorConstruction::SetKaptonSTLPath)
      .SetGuidance("Path to the KAPTON substrate solid GDML file (must "
                    "already be thickened - see diff_geom_macros/"
                    "_differentiated_mesh.py build_differentiated_kapton()). "
                    "Command name kept as 'kaptonStlFile' for backward "
                    "compatibility - expected file is now .gdml, not .stl.");
  fMessenger->DeclareMethod("sensorThickness", &DetectorConstruction::SetSensorThickness)
      .SetGuidance("Nominal silicon thickness in mm - sanity-check/metadata "
                    "only, does not affect the geometry actually loaded.")
      .SetUnit("mm");
  fMessenger->DeclareMethod("substrateThickness", &DetectorConstruction::SetSubstrateThickness)
      .SetGuidance("Nominal Kapton thickness in mm - sanity-check/metadata "
                    "only, does not affect the geometry actually loaded.")
      .SetUnit("mm");
  fMessenger->DeclareMethod("foldState", &DetectorConstruction::SetFoldState)
      .SetGuidance("Fold/deploy fraction, 0 (stowed) to 1 (deployed). "
                    "METADATA ONLY, recorded in the output ntuple/filename "
                    "so a run is self-describing. Selecting a different "
                    "fold state requires loading different STL files via "
                    "stlFile/kaptonStlFile - this value does not itself "
                    "change the loaded geometry.");
}
G4String DetectorConstruction::DeriveStructureTag(const G4String& stlPath)
{
  std::string path(stlPath);
  size_t slash = path.find_last_of("/\\");
  std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
  size_t dot = base.find_last_of('.');
  if (dot != std::string::npos) base = base.substr(0, dot);
  const std::string suffix = "_silicon";
  if (base.size() > suffix.size() &&
      base.compare(base.size() - suffix.size(), suffix.size(), suffix) == 0) {
    base = base.substr(0, base.size() - suffix.size());
  }
  if (base.empty()) base = "unknown";
  return G4String(base);
}
namespace {
G4String DeriveGdmlVolumeName(const G4String& gdmlPath)
{
  std::string path(gdmlPath);
  size_t slash = path.find_last_of("/\\");
  std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
  size_t dot = base.find_last_of('.');
  if (dot != std::string::npos) base = base.substr(0, dot);
  return G4String(base + "_volume");
}
}
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* worldMat = nist->FindOrBuildMaterial("G4_AIR");
  auto* worldSolid = new G4Box("World",
                                OrigamiDet::kWorldHalfSize,
                                OrigamiDet::kWorldHalfSize,
                                OrigamiDet::kWorldHalfSize);
  auto* worldLogical = new G4LogicalVolume(worldSolid, worldMat, "World");
  auto* worldPhysical = new G4PVPlacement(
      nullptr, {}, worldLogical, "World", nullptr, false, 0, true);
  G4GDMLParser siliconParser;
  siliconParser.Read(fSiliconSTLPath, false);
  G4LogicalVolume* siliconGdmlVolume = siliconParser.GetVolume(DeriveGdmlVolumeName(fSiliconSTLPath));
  if (!siliconGdmlVolume) {
    G4ExceptionDescription msg;
    msg << "Could not find expected volume in GDML file '" << fSiliconSTLPath
        << "' - check diff_geom_macros/_gdml_export.py's naming convention "
        << "still matches DeriveGdmlVolumeName() below.";
    G4Exception("DetectorConstruction::Construct", "GDMLVolumeNotFound",
                FatalException, msg);
  }
  G4VSolid* siliconSolid = siliconGdmlVolume->GetSolid();
  SanityCheckThickness(siliconSolid, fSensorThickness_mm, "Silicon solid");
  G4ThreeVector centroidOffset =
      ComputeBoundingBoxCenter(siliconSolid, "Silicon solid");
  G4ThreeVector placementOffset = -centroidOffset;
  G4Material* siliconMat = nist->FindOrBuildMaterial("G4_Si");
  fSensorLogical = new G4LogicalVolume(siliconSolid, siliconMat, "SensorLogical");
  auto* siVisAtt = new G4VisAttributes(G4Colour(0.3, 0.6, 0.9, 0.6));
  siVisAtt->SetForceSolid(true);
  fSensorLogical->SetVisAttributes(siVisAtt);
  new G4PVPlacement(nullptr, placementOffset, fSensorLogical, "SensorPhysical",
                     worldLogical, false, 0, true);
  G4GDMLParser kaptonParser;
  kaptonParser.Read(fKaptonSTLPath, false);
  G4LogicalVolume* kaptonGdmlVolume = kaptonParser.GetVolume(DeriveGdmlVolumeName(fKaptonSTLPath));
  if (!kaptonGdmlVolume) {
    G4ExceptionDescription msg;
    msg << "Could not find expected volume in GDML file '" << fKaptonSTLPath
        << "' - check diff_geom_macros/_gdml_export.py's naming convention "
        << "still matches DeriveGdmlVolumeName() below.";
    G4Exception("DetectorConstruction::Construct", "GDMLVolumeNotFound",
                FatalException, msg);
  }
  G4VSolid* kaptonSolid = kaptonGdmlVolume->GetSolid();
  SanityCheckThickness(kaptonSolid, fSubstrateThickness_mm, "Kapton solid");
  G4Material* kaptonMat = nist->FindOrBuildMaterial("G4_KAPTON");
  fSubstrateLogical = new G4LogicalVolume(kaptonSolid, kaptonMat, "SubstrateLogical");
  auto* kaptonVisAtt = new G4VisAttributes(G4Colour(0.9, 0.6, 0.2, 0.4));
  kaptonVisAtt->SetForceSolid(true);
  fSubstrateLogical->SetVisAttributes(kaptonVisAtt);
  new G4PVPlacement(nullptr, placementOffset, fSubstrateLogical, "SubstratePhysical",
                     worldLogical, false, 0, true);
  fsStructureTag = DeriveStructureTag(fSiliconSTLPath);
  fsFoldState = fFoldState;
  fsCentroidOffset = placementOffset;
  G4cout << "[DetectorConstruction] structureTag=" << fsStructureTag
         << " foldState=" << fsFoldState
         << " centroidOffset(mm)=(" << fsCentroidOffset.x() << ", "
         << fsCentroidOffset.y() << ", " << fsCentroidOffset.z() << ")"
         << G4endl;
  return worldPhysical;
}
void DetectorConstruction::ConstructSDandField()
{
  auto* siliconSD = new SensorSD("SensorSD", "SensorHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(siliconSD);
  SetSensitiveDetector("SensorLogical", siliconSD);
  auto* substrateSD = new SensorSD("SubstrateSD", "SubstrateHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(substrateSD);
  SetSensitiveDetector("SubstrateLogical", substrateSD);
}
