#ifndef EVENT_ACTION_HH
#define EVENT_ACTION_HH
#include "G4UserEventAction.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"
class EventAction : public G4UserEventAction
{
public:
  EventAction();
  ~EventAction() override = default;
  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event* event) override;
  static void SetCurrentLabThetaDeg(G4double deg) { fsCurrentLabThetaDeg = deg; }
  static void SetCurrentLaunchDir(const G4ThreeVector& dir) { fsCurrentLaunchDir = dir; }
  static void SetCurrentVertex(const G4ThreeVector& v) { fsCurrentVertex_mm = v; }
private:
  static G4double ComputeLocalIncidenceDeg(const G4ThreeVector& momentumDir,
                                            const G4ThreeVector& surfaceNormal);
  G4int fHCID = -1;
  G4int fSubstrateHCID = -1;
  static G4double fsCurrentLabThetaDeg;
  static G4ThreeVector fsCurrentVertex_mm;
  static G4ThreeVector fsCurrentLaunchDir;
};
#endif
