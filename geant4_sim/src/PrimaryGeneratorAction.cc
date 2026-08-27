#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "DetectorConstants.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <cmath>
G4double PrimaryGeneratorAction::fsMomentum_GeV = 5.0;
G4double PrimaryGeneratorAction::fsVertexSmear_mm = 0.0;
PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  fGun = std::make_unique<G4ParticleGun>(1);
  G4ParticleTable* table = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* pion = table->FindParticle("pi+");
  fGun->SetParticleDefinition(pion);
  DefineMessenger();
}
PrimaryGeneratorAction::~PrimaryGeneratorAction() = default;
void PrimaryGeneratorAction::DefineMessenger()
{
  fMessenger = std::make_unique<G4GenericMessenger>(
      this, "/origami/gun/", "Pion gun controls");
  fMessenger->DeclareMethod("mode", &PrimaryGeneratorAction::SetMode)
      .SetGuidance("Gun mode: 'isotropic4pi' (MARK 2 DEFAULT - vertex at/near "
                    "detector center, direction sampled uniformly over full "
                    "4pi solid angle - appropriate for barrel-tracker-style "
                    "acceptance studies) or 'mark1cone' (MARK 1 LEGACY - "
                    "single-sided beam from a disk, tilted within a cone; "
                    "kept for validation/cross-check runs, e.g. "
                    "validate_flat_plate.mac). Default: isotropic4pi.");
  fMessenger->DeclareMethod("vertexSmearMm", &PrimaryGeneratorAction::SetVertexSmearMm)
      .SetGuidance("Isotropic4pi mode only. Gaussian sigma (mm) for a "
                    "3D-smeared vertex position around the detector center, "
                    "mimicking a realistic luminous region. Default 0 = a "
                    "fixed point vertex (no smear).");
  fMessenger->DeclareMethod("vertexZOffsetMm", &PrimaryGeneratorAction::SetVertexZOffsetMm)
      .SetGuidance("Isotropic4pi mode only. Fixed z-offset (mm) of the "
                    "vertex from the detector's geometric center. Default 0.");
  fMessenger->DeclareMethod("momentumGeV", &PrimaryGeneratorAction::SetPionMomentumGeV)
      .SetGuidance("Particle momentum in GeV/c. Shared by both modes. "
                    "Default 5.0 (paper spec).");
  fMessenger->DeclareMethod("fixedAngleDeg", &PrimaryGeneratorAction::SetFixedAngleDeg)
      .SetGuidance("MARK 1 LEGACY (mark1cone mode only). Fixed incidence "
                    "angle in degrees from normal (0 = head-on). Ignored in "
                    "isotropic4pi mode.");
  fMessenger->DeclareMethod("randomAngleMaxDeg", &PrimaryGeneratorAction::SetRandomAngleMaxDeg)
      .SetGuidance("MARK 1 LEGACY (mark1cone mode only). Switch to "
                    "random-angle cone mode; pions sampled uniformly in solid "
                    "angle up to this cone half-angle (degrees). Ignored in "
                    "isotropic4pi mode.");
  fMessenger->DeclareMethod("gunZmm", &PrimaryGeneratorAction::SetGunZPosition_mm)
      .SetGuidance("MARK 1 LEGACY (mark1cone mode only). Z position (mm) of "
                    "the gun's starting disk. Ignored in isotropic4pi mode. "
                    "Accepted without error so old Mark 1 macros don't fatal "
                    "when run against Mark 2 binaries - see MARK_HISTORY.md.");
  fMessenger->DeclareMethod("diskRadiusMm", &PrimaryGeneratorAction::SetDiskRadius_mm)
      .SetGuidance("MARK 1 LEGACY (mark1cone mode only). Radius (mm) of the "
                    "transverse disk from which primaries are sampled. "
                    "Ignored in isotropic4pi mode.");
}
void PrimaryGeneratorAction::SetMode(G4String mode)
{
  if (mode == "isotropic4pi") {
    fMode = GunMode::kIsotropic4Pi;
  } else if (mode == "mark1cone") {
    fMode = GunMode::kMark1Cone;
  } else {
    G4cerr << "[PrimaryGeneratorAction] WARNING: unrecognized gun mode '"
           << mode << "'. Valid options: isotropic4pi, mark1cone. "
           << "Leaving mode unchanged (currently "
           << (fMode == GunMode::kIsotropic4Pi ? "isotropic4pi" : "mark1cone")
           << ")." << G4endl;
    return;
  }
  fModeAnnounced = false;
}
void PrimaryGeneratorAction::SetGunZPosition_mm(G4double z)
{
  fGunZ_mm = z;
  G4double clearance = OrigamiDet::kWorldHalfSize - std::abs(z) * mm;
  if (clearance < OrigamiDet::kMinGunBoundaryClearance) {
    G4cerr << "[PrimaryGeneratorAction] WARNING: gunZmm=" << z
           << " mm leaves only " << clearance / mm
           << " mm clearance to the world boundary (half-size="
           << OrigamiDet::kWorldHalfSize / mm << " mm). This only matters "
           << "if gun mode is set to mark1cone - isotropic4pi mode ignores "
           << "this parameter entirely." << G4endl;
  }
}
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  if (!fModeAnnounced) {
    G4cout << "[PrimaryGeneratorAction] Active gun mode: "
           << (fMode == GunMode::kIsotropic4Pi ? "isotropic4pi (MARK 2 default)"
                                                : "mark1cone (MARK 1 legacy)")
           << G4endl;
    fModeAnnounced = true;
  }
  if (fMode == GunMode::kIsotropic4Pi) {
    GeneratePrimariesIsotropic4Pi(event);
  } else {
    GeneratePrimariesMark1Cone(event);
  }
}
void PrimaryGeneratorAction::GeneratePrimariesIsotropic4Pi(G4Event* event)
{
  G4double cosTheta = 2.0 * G4UniformRand() - 1.0;
  G4double theta = std::acos(cosTheta);
  G4double phi = 2.0 * CLHEP::pi * G4UniformRand();
  G4double dx = std::sin(theta) * std::cos(phi);
  G4double dy = std::sin(theta) * std::sin(phi);
  G4double dz = std::cos(theta);
  G4ThreeVector direction(dx, dy, dz);
  G4double labThetaDeg = theta / deg;
  EventAction::SetCurrentLabThetaDeg(labThetaDeg);
  EventAction::SetCurrentLaunchDir(direction);
  fGun->SetParticleMomentumDirection(direction);
  G4double vx = 0.0, vy = 0.0, vz = fVertexZOffset_mm;
  if (fVertexSmear_mm > 0.0) {
    vx += G4RandGauss::shoot(0.0, fVertexSmear_mm);
    vy += G4RandGauss::shoot(0.0, fVertexSmear_mm);
    vz += G4RandGauss::shoot(0.0, fVertexSmear_mm);
  }
  G4ThreeVector vertex(vx, vy, vz);
  fGun->SetParticlePosition(vertex * mm);
  EventAction::SetCurrentVertex(vertex * mm);
  fGun->SetParticleMomentum(fMomentum_GeV * GeV);
  fGun->GeneratePrimaryVertex(event);
}
void PrimaryGeneratorAction::GeneratePrimariesMark1Cone(G4Event* event)
{
  G4double thetaDeg;
  if (fMark1RandomAngle) {
    G4double cosMax = std::cos(fMark1MaxAngleDeg * deg);
    G4double cosTheta = cosMax + (1.0 - cosMax) * G4UniformRand();
    thetaDeg = std::acos(cosTheta) / deg;
  } else {
    thetaDeg = fFixedAngleDeg;
  }
  G4double phi = 2.0 * CLHEP::pi * G4UniformRand();
  G4double theta = thetaDeg * deg;
  EventAction::SetCurrentLabThetaDeg(thetaDeg);
  G4double dx = std::sin(theta) * std::cos(phi);
  G4double dy = std::sin(theta) * std::sin(phi);
  G4double dz = -std::cos(theta);
  G4ThreeVector direction(dx, dy, dz);
  EventAction::SetCurrentLaunchDir(direction);
  fGun->SetParticleMomentumDirection(direction);
  G4double diskRadius_mm = fDiskRadius_mm;
  G4double r = diskRadius_mm * std::sqrt(G4UniformRand());
  G4double diskPhi = 2.0 * CLHEP::pi * G4UniformRand();
  G4double x0 = r * std::cos(diskPhi);
  G4double y0 = r * std::sin(diskPhi);
  G4ThreeVector vertex(x0, y0, fGunZ_mm);
  fGun->SetParticlePosition(vertex);
  EventAction::SetCurrentVertex(vertex);
  fGun->SetParticleMomentum(fMomentum_GeV * GeV);
  fGun->GeneratePrimaryVertex(event);
}
