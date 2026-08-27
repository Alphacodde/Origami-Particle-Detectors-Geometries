#include "SensorSD.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4TouchableHistory.hh"
#include "G4NavigationHistory.hh"
#include "G4VSolid.hh"
#include "G4AffineTransform.hh"
G4ThreadLocal G4Allocator<SensorHit>* SensorHitAllocator = nullptr;
SensorSD::SensorSD(const G4String& name, const G4String& hitsCollectionName)
  : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
}
void SensorSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new SensorHitsCollection(SensitiveDetectorName, collectionName[0]);
  if (fHCID < 0) {
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  }
  hce->AddHitsCollection(fHCID, fHitsCollection);
}
G4bool SensorSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  G4double edep = step->GetTotalEnergyDeposit();
  G4double stepLength = step->GetStepLength();
  if (stepLength <= 0. && edep <= 0.) return false;
  auto* hit = new SensorHit();
  const G4Track* track = step->GetTrack();
  hit->trackID = track->GetTrackID();
  hit->parentID = track->GetParentID();
  hit->particleName = track->GetDefinition()->GetParticleName();
  hit->edep_MeV = edep / CLHEP::MeV;
  hit->stepLength_mm = stepLength / CLHEP::mm;
  hit->prePos_mm = step->GetPreStepPoint()->GetPosition() / CLHEP::mm;
  hit->postPos_mm = step->GetPostStepPoint()->GetPosition() / CLHEP::mm;
  hit->preMomentumDir  = step->GetPreStepPoint()->GetMomentumDirection();
  hit->postMomentumDir = step->GetPostStepPoint()->GetMomentumDirection();
  hit->kineticEnergy_MeV = track->GetKineticEnergy() / CLHEP::MeV;
  const G4TouchableHistory* touchable =
      static_cast<const G4TouchableHistory*>(step->GetPreStepPoint()->GetTouchable());
  if (touchable && touchable->GetVolume()) {
    G4VSolid* solid = touchable->GetSolid();
    if (solid) {
      G4AffineTransform globalToLocal = touchable->GetHistory()->GetTopTransform();
      G4ThreeVector localPos = globalToLocal.TransformPoint(
          step->GetPreStepPoint()->GetPosition());
      G4ThreeVector localNormal = solid->SurfaceNormal(localPos);
      G4AffineTransform localToGlobal = globalToLocal.Inverse();
      hit->localNormal = localToGlobal.TransformAxis(localNormal);
    }
  }
  fHitsCollection->insert(hit);
  return true;
}
