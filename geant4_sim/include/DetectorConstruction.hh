#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH
#include "G4VUserDetectorConstruction.hh"
#include "G4GenericMessenger.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include <memory>
class G4LogicalVolume;
class G4VPhysicalVolume;
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction() override;
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;
  void SetSiliconSTLPath(G4String path) { fSiliconSTLPath = path; }
  void SetKaptonSTLPath(G4String path) { fKaptonSTLPath = path; }
  void SetSensorThickness(G4double t) { fSensorThickness_mm = t; }
  void SetSubstrateThickness(G4double t) { fSubstrateThickness_mm = t; }
  void SetFoldState(G4double f) { fFoldState = f; }
  static G4String GetStructureTag() { return fsStructureTag; }
  static G4double GetFoldStateStatic() { return fsFoldState; }
  static G4ThreeVector GetCentroidOffsetStatic() { return fsCentroidOffset; }
private:
  G4String fSiliconSTLPath = "geometry/miura_deployed_silicon.stl";
  G4String fKaptonSTLPath  = "geometry/miura_deployed_kapton.stl";
  G4double fSensorThickness_mm = 0.300;
  G4double fSubstrateThickness_mm = 0.050;
  G4double fFoldState = 1.0;
  G4LogicalVolume* fSensorLogical = nullptr;
  G4LogicalVolume* fSubstrateLogical = nullptr;
  std::unique_ptr<G4GenericMessenger> fMessenger;
  void DefineMessenger();
  static G4String DeriveStructureTag(const G4String& stlPath);
  static G4String fsStructureTag;
  static G4double fsFoldState;
  static G4ThreeVector fsCentroidOffset;
};
#endif
