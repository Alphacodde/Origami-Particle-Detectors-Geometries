#ifndef DETECTOR_CONSTANTS_HH
#define DETECTOR_CONSTANTS_HH
#include "G4SystemOfUnits.hh"
namespace OrigamiDet {
  constexpr G4double kWorldSize = 1.0 * m;
  constexpr G4double kWorldHalfSize = kWorldSize / 2.0;
  constexpr G4double kMinGunBoundaryClearance = 20.0 * mm;
}
#endif
