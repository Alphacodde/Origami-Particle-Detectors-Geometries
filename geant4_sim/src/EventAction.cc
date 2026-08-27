#include "EventAction.hh"
#include "SensorSD.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
G4double EventAction::fsCurrentLabThetaDeg = 0.0;
G4ThreeVector EventAction::fsCurrentVertex_mm = G4ThreeVector(0., 0., 0.);
G4ThreeVector EventAction::fsCurrentLaunchDir = G4ThreeVector(0., 0., 1.);
EventAction::EventAction() = default;
void EventAction::BeginOfEventAction(const G4Event*)
{
}
G4double EventAction::ComputeLocalIncidenceDeg(const G4ThreeVector& momentumDir,
                                                const G4ThreeVector& surfaceNormal)
{
  if (surfaceNormal.mag2() <= 0.0) return -1.0;
  G4double rawDeg = momentumDir.angle(surfaceNormal) / deg;
  return (rawDeg > 90.0) ? (180.0 - rawDeg) : rawDeg;
}
void EventAction::EndOfEventAction(const G4Event* event)
{
  auto* sdManager = G4SDManager::GetSDMpointer();
  if (fHCID < 0) {
    fHCID = sdManager->GetCollectionID("SensorHitsCollection");
  }
  if (fSubstrateHCID < 0) {
    fSubstrateHCID = sdManager->GetCollectionID("SubstrateHitsCollection");
  }
  auto* hce = event->GetHCofThisEvent();
  if (!hce) {
    return;
  }
  auto* hitsCollection = static_cast<SensorHitsCollection*>(hce->GetHC(fHCID));
  auto* substrateHitsCollection =
      static_cast<SensorHitsCollection*>(hce->GetHC(fSubstrateHCID));
  G4double siliconPathLength_mm = 0.0;
  G4double kaptonPathLength_mm = 0.0;
  G4double totalEdep_MeV = 0.0;
  G4int nHitsPrimary = 0;
  G4bool hitDetector = false;
  G4ThreeVector firstEntryPos;
  G4double localIncidenceDeg = -1.0;
  G4bool haveEntry = false;
  G4ThreeVector dirAtFirstHit;
  G4ThreeVector dirAtLastHit;
  if (hitsCollection) {
    G4int nHits = hitsCollection->entries();
    for (G4int i = 0; i < nHits; ++i) {
      SensorHit* hit = (*hitsCollection)[i];
      if (hit->parentID != 0) continue;
      hitDetector = true;
      siliconPathLength_mm += hit->stepLength_mm;
      totalEdep_MeV += hit->edep_MeV;
      nHitsPrimary++;
      if (!haveEntry) {
        firstEntryPos = hit->prePos_mm;
        dirAtFirstHit = fsCurrentLaunchDir;
        localIncidenceDeg = ComputeLocalIncidenceDeg(fsCurrentLaunchDir, hit->localNormal);
        haveEntry = true;
      }
      dirAtLastHit = hit->postMomentumDir;
    }
  }
  if (substrateHitsCollection) {
    G4int nSubstrateHits = substrateHitsCollection->entries();
    for (G4int i = 0; i < nSubstrateHits; ++i) {
      SensorHit* hit = (*substrateHitsCollection)[i];
      if (hit->parentID != 0) continue;
      kaptonPathLength_mm += hit->stepLength_mm;
    }
  }
  G4double totalPathLength_mm = siliconPathLength_mm + kaptonPathLength_mm;
  G4double scatterAngleDeg = -1.0;
  if (haveEntry) {
    scatterAngleDeg = dirAtFirstHit.angle(dirAtLastHit) / deg;
  }
  auto* analysisManager = G4AnalysisManager::Instance();
  G4int ntupleId = 0;
  analysisManager->FillNtupleIColumn(ntupleId, 0, event->GetEventID());
  analysisManager->FillNtupleDColumn(ntupleId, 1, fsCurrentLabThetaDeg);
  analysisManager->FillNtupleDColumn(ntupleId, 2, siliconPathLength_mm);
  analysisManager->FillNtupleDColumn(ntupleId, 3, totalEdep_MeV);
  analysisManager->FillNtupleDColumn(ntupleId, 4, haveEntry ? firstEntryPos.x() : -9999.0);
  analysisManager->FillNtupleDColumn(ntupleId, 5, haveEntry ? firstEntryPos.y() : -9999.0);
  analysisManager->FillNtupleDColumn(ntupleId, 6, haveEntry ? firstEntryPos.z() : -9999.0);
  analysisManager->FillNtupleIColumn(ntupleId, 7, nHitsPrimary);
  analysisManager->FillNtupleIColumn(ntupleId, 8, hitDetector ? 1 : 0);
  analysisManager->FillNtupleSColumn(ntupleId, 9, DetectorConstruction::GetStructureTag());
  analysisManager->FillNtupleDColumn(ntupleId, 10, DetectorConstruction::GetFoldStateStatic());
  analysisManager->FillNtupleDColumn(ntupleId, 11, localIncidenceDeg);
  analysisManager->FillNtupleDColumn(ntupleId, 12, fsCurrentVertex_mm.x());
  analysisManager->FillNtupleDColumn(ntupleId, 13, fsCurrentVertex_mm.y());
  analysisManager->FillNtupleDColumn(ntupleId, 14, fsCurrentVertex_mm.z());
  analysisManager->FillNtupleDColumn(ntupleId, 15, scatterAngleDeg);
  analysisManager->FillNtupleDColumn(ntupleId, 16, kaptonPathLength_mm);
  analysisManager->FillNtupleDColumn(ntupleId, 17, totalPathLength_mm);
  analysisManager->AddNtupleRow(ntupleId);
}
