#ifndef PRIMARY_GENERATOR_ACTION_HH
#define PRIMARY_GENERATOR_ACTION_HH
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4GenericMessenger.hh"
#include <memory>
class G4Event;
enum class GunMode { kIsotropic4Pi, kMark1Cone };
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;
  void GeneratePrimaries(G4Event* event) override;
  void SetMode(G4String mode);
  void SetVertexSmearMm(G4double r) { fVertexSmear_mm = r; fsVertexSmear_mm = r; }
  void SetVertexZOffsetMm(G4double z) { fVertexZOffset_mm = z; }
  void SetFixedAngleDeg(G4double deg) { fFixedAngleDeg = deg; fMark1RandomAngle = false; }
  void SetRandomAngleMaxDeg(G4double deg) { fMark1MaxAngleDeg = deg; fMark1RandomAngle = true; }
  void SetPionMomentumGeV(G4double p) { fMomentum_GeV = p; fsMomentum_GeV = p; }
  void SetGunZPosition_mm(G4double z);
  void SetDiskRadius_mm(G4double r) { fDiskRadius_mm = r; }
  static G4double GetMomentumGeVStatic() { return fsMomentum_GeV; }
  static G4double GetVertexSmearMmStatic() { return fsVertexSmear_mm; }
private:
  std::unique_ptr<G4ParticleGun> fGun;
  std::unique_ptr<G4GenericMessenger> fMessenger;
  void DefineMessenger();
  void GeneratePrimariesIsotropic4Pi(G4Event* event);
  void GeneratePrimariesMark1Cone(G4Event* event);
  GunMode fMode = GunMode::kIsotropic4Pi;
  G4bool fModeAnnounced = false;
  G4double fMomentum_GeV = 5.0;
  G4double fVertexSmear_mm = 0.0;
  G4double fVertexZOffset_mm = 0.0;
  G4double fFixedAngleDeg = 0.0;
  G4double fMark1MaxAngleDeg = 60.0;
  G4bool fMark1RandomAngle = false;
  G4double fGunZ_mm = 500.0;
  G4double fDiskRadius_mm = 150.0;
  static G4double fsMomentum_GeV;
  static G4double fsVertexSmear_mm;
};
#endif
