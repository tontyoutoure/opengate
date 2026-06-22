#pragma once
#include "G4ThreeVector.hh"
#include "G4Types.hh"
#include <vector>
class GateWindowAcceptanceManager {
  std::vector<G4double> fA1, fA2, fB1, fB2, fCosPhi, fSinPhi, fPlaneDistance;
  bool TestForOneWindow(size_t i, const G4ThreeVector &position,
                        const G4ThreeVector &momentum_direction) const;
  bool TestIfAccept(const G4ThreeVector &position,
                    const G4ThreeVector &momentum_direction);
};
