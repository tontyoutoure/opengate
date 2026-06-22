#include "GateWindowAcceptanceManager.h"
#include "fmt/format.h"

bool GateWindowAcceptanceManager::TestIfAccept(
    const G4ThreeVector &position, const G4ThreeVector &momentum_direction) {
  bool result = false;
  for (size_t i = 0; i < fA1.size(); i++) {
    result |= TestForOneWindow(i, position, momentum_direction);
  }
  return result;
}

bool GateWindowAcceptanceManager::TestForOneWindow(
    size_t i, const G4ThreeVector &position,
    const G4ThreeVector &momentum_direction) const {
  G4double sin_phi = fSinPhi[i];
  G4double cos_phi = fCosPhi[i];
  G4double pos_x_rotated = position.x() * cos_phi + position.y() * sin_phi;
  G4double pos_y_rotated = -position.x() * sin_phi + position.y() * cos_phi;
  G4double dir_x_rotated =
      momentum_direction.x() * cos_phi + momentum_direction.y() * sin_phi;
  G4double dir_y_rotated =
      -momentum_direction.x() * sin_phi + momentum_direction.y() * cos_phi;
  G4String error_msg = fmt::format(
      "position ({}, {}, {}) is outside the plane distance {} for window {}",
      position.x(), position.y(), position.z(), fPlaneDistance[i], i);

  G4double pd = fPlaneDistance[i];
  if (pos_x_rotated >= pd) {
    G4Exception("GateWindowAcceptanceManager::TestForOneWindow",
                "OutOfWindowError", FatalException, error_msg);
  }

  G4double intersect_a =
      (pd - pos_x_rotated) / dir_x_rotated * dir_y_rotated + pos_y_rotated;
  G4double intersect_b =
      (pd - pos_x_rotated) / dir_x_rotated * momentum_direction.z() +
      position.z();
  return intersect_a < fA2[i] && intersect_a >= fA1[i] &&
         intersect_b < fB2[i] && intersect_b >= fB1[i];
}
