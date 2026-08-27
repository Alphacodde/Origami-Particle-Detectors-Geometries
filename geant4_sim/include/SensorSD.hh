#ifndef SENSOR_SD_HH
#define SENSOR_SD_HH
#include "G4VSensitiveDetector.hh"
#include "G4THitsCollection.hh"
#include "G4VHit.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"
class G4Step;
class G4HCofThisEvent;
class SensorHit : public G4VHit
{
public:
  SensorHit() = default;
  ~SensorHit() override = default;
  void* operator new(size_t);
  void operator delete(void*);
  G4int eventID = -1;
  G4int trackID = -1;
  G4int parentID = -1;
  G4String particleName;
  G4double edep_MeV = 0.0;
  G4double stepLength_mm = 0.0;
  G4ThreeVector prePos_mm;
  G4ThreeVector postPos_mm;
  G4ThreeVector preMomentumDir;
  G4ThreeVector postMomentumDir;
  G4ThreeVector localNormal;
  G4double kineticEnergy_MeV = 0.0;
};
using SensorHitsCollection = G4THitsCollection<SensorHit>;
extern G4ThreadLocal G4Allocator<SensorHit>* SensorHitAllocator;
inline void* SensorHit::operator new(size_t)
{
  if (!SensorHitAllocator) SensorHitAllocator = new G4Allocator<SensorHit>;
  return (void*)SensorHitAllocator->MallocSingle();
}
inline void SensorHit::operator delete(void* hit)
{
  SensorHitAllocator->FreeSingle((SensorHit*)hit);
}
class SensorSD : public G4VSensitiveDetector
{
public:
  SensorSD(const G4String& name, const G4String& hitsCollectionName);
  ~SensorSD() override = default;
  void Initialize(G4HCofThisEvent* hce) override;
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history) override;
private:
  SensorHitsCollection* fHitsCollection = nullptr;
  G4int fHCID = -1;
};
#endif
