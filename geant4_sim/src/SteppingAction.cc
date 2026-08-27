#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4VVisManager* visManager = G4VVisManager::GetConcreteInstance();
  if (!visManager) return;
  const G4VPhysicalVolume* prePV = step->GetPreStepPoint()->GetPhysicalVolume();
  if (!prePV) return;
  const G4String& volName = prePV->GetLogicalVolume()->GetName();
  G4bool inSilicon = (volName == "SensorLogical");
  G4bool inKapton  = (volName == "SubstrateLogical");
  if (!inSilicon && !inKapton) return;
  G4Colour colour = inSilicon ? G4Colour(0.0, 1.0, 0.0)
                               : G4Colour(1.0, 0.0, 0.0);
  G4Polyline polyline;
  G4VisAttributes attribs(colour);
  polyline.SetVisAttributes(attribs);
  polyline.push_back(step->GetPreStepPoint()->GetPosition());
  polyline.push_back(step->GetPostStepPoint()->GetPosition());
  visManager->Draw(polyline);
}
