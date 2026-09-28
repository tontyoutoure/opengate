/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GateSingleParticleSourceWindowFD.h"
#include "G4Threading.hh"
#include "GateHelpersDict.h"
#include "Randomize.hh"
#include <G4Event.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <fmt/core.h>
#include <string>

GateSingleParticleSourceWindowFD::GateSingleParticleSourceWindowFD(
    std::string mother_volume)
    : GateSingleParticleSource(mother_volume) {}

void GateSingleParticleSourceWindowFD::SetParameters(G4double a1, G4double a2,
                                                     G4double b1, G4double b2,
                                                     G4double plane_distance,
                                                     G4double plane_phi) {
  fA1 = a1;
  fA2 = a2;
  fB1 = b1;
  fB2 = b2;
  fPlaneDistance = plane_distance;
  fSinPlanePhi = sin(plane_phi);
  fCosPlanePhi = cos(plane_phi);
}

G4bool GateSingleParticleSourceWindowFD::CheckPosDirValid(
    const G4ThreeVector &pos, const G4ThreeVector &dir) const {
  // compare theta with cos2, to avoid complex calculation
  //   G4double cot2theta = std::copysign(1.0, dir.z()) * dir.z() * dir.z() /
  //                        (dir.x() * dir.x() + dir.y() * dir.y());

  G4double x0 = pos.x() * fCosPlanePhi + pos.y() * fSinPlanePhi;
  G4double y0 = -pos.x() * fSinPlanePhi + pos.y() * fCosPlanePhi;
  G4double a1_rel = fA1 - y0;
  G4double a2_rel = fA2 - y0;
  G4double b1_rel = fB1 - pos.z();
  G4double b2_rel = fB2 - pos.z();
  G4double d_rel = fPlaneDistance - x0;
  G4double dir_x_rotated = dir.x() * fCosPlanePhi + dir.y() * fSinPlanePhi;
  G4double dir_y_rotated = -dir.x() * fSinPlanePhi + dir.y() * fCosPlanePhi;

  G4double intersect_b = d_rel / dir_x_rotated * dir.z() + pos.z();
  G4double intersect_a = d_rel / dir_x_rotated * dir_y_rotated + y0;
  return intersect_a <= fA2 && intersect_a >= fA1 && intersect_b <= fB2 &&
         intersect_b >= fB1;
}

void GateSingleParticleSourceWindowFD::GeneratePosDir() {
  fSkippedCount = 0;
  while (true) {
    fCurrentPos = fPositionGenerator->VGenerateOne();
    fCurrentDir = fDirectionGenerator->VGenerateOne();
    if (CheckPosDirValid(fCurrentPos, fCurrentDir)) {
      break;
    }
    fSkippedCount++;
  }
}

void GateSingleParticleSourceWindowFD::GeneratePrimaryVertex(G4Event *event) {
  G4PrimaryVertex *vertex = new G4PrimaryVertex(fCurrentPos, particle_time);

  // Set placement relative to attached volume
  // DD(particle_momentum_direction);

  G4double energy = fEnergyGenerator->VGenerateOne(fParticleDefinition);

  // one single particle
  auto *particle = new G4PrimaryParticle(fParticleDefinition);
  particle->SetKineticEnergy(energy);
  particle->SetMass(fMass);
  particle->SetMomentumDirection(fCurrentDir);
  particle->SetCharge(fCharge);
  particle->SetWeight(1.0);
  if (fPolarizationFlag)
    particle->SetPolarization(fPolarization);

  // set vertex
  vertex->SetPrimary(particle);
  event->AddPrimaryVertex(vertex);
}
