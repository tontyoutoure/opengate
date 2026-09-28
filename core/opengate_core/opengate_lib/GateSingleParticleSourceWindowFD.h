/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateSingleParticleSourceWindowFD_h
#define GateSingleParticleSourceWindowFD_h

#include "GateSingleParticleSource.h"
#include <G4ThreeVector.hh>
#include <G4Types.hh>
#include <pybind11/pytypes.h>

namespace py = pybind11;

class GateSingleParticleSourceWindowFD : public GateSingleParticleSource {
public:
  explicit GateSingleParticleSourceWindowFD(std::string mother_volume);
  ~GateSingleParticleSourceWindowFD() override = default;
  void GeneratePrimaryVertex(G4Event *event) override;
  void GeneratePosDir();
  void SetParameters(G4double a1, G4double a2, G4double b1, G4double b2,
                     G4double plane_distance, G4double plane_phi);
  G4double InitializeBeforeRun(G4double &act_ratio, G4double &max_solid_angle);
  G4long GetSkippedCount() const { return fSkippedCount; }

private:
  G4double GetSolidAngle(
      const G4ThreeVector &pos) const; // get solid angle for the window
  G4bool CheckPosDirValid(const G4ThreeVector &pos,
                          const G4ThreeVector &dir)
      const; // check if the ray can pass through the window
  void SetPhiTheta(
      const G4ThreeVector &pos) const; // set the phi and theta of the direction
                                       // distribution according to the position
  G4double fPlaneDistance{NAN};
  G4double fSinPlanePhi{NAN}, fCosPlanePhi{NAN};
  G4double fA1{NAN}, fA2{NAN}, fB1{NAN}, fB2{NAN};
  G4String fSourceName;
  G4ThreeVector fCurrentDir;
  G4ThreeVector fCurrentPos;
  G4long fSkippedCount;
  //   G4bool fPosGenerated = false;
};

#endif // GateSingleParticleSourceWindowFD_h
