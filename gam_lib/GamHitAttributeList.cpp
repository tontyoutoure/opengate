/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */


// Macros to reduce the code size
// Use FILLFS when step is not used to avoid warning

#define FILLF [=] (GamVHitAttribute *att, G4Step *step, G4TouchableHistory *)
#define FILLFS [=] (GamVHitAttribute *att, G4Step *, G4TouchableHistory *)

void GamHitAttributeManager::InitializeAllHitAttributes() {
    DefineHitAttribute("TotalEnergyDeposit", 'D',
                       FILLF { att->FillDValue(step->GetTotalEnergyDeposit()); }
    );
    DefineHitAttribute("KineticEnergy", 'D',
                       FILLF { att->FillDValue(step->GetPostStepPoint()->GetKineticEnergy()); }
    );
    DefineHitAttribute("LocalTime", 'D',
                       FILLF { att->FillDValue(step->GetPostStepPoint()->GetLocalTime()); }
    );
    DefineHitAttribute("GlobalTime", 'D',
                       FILLF { att->FillDValue(step->GetPostStepPoint()->GetGlobalTime()); }
    );
    DefineHitAttribute("Weight", 'D',
                       FILLF { att->FillDValue(step->GetTrack()->GetWeight()); }
    );
    DefineHitAttribute("TrackID", 'I',
                       FILLF { att->FillIValue(step->GetTrack()->GetTrackID()); }
    );
    DefineHitAttribute("EventID", 'I',
                       FILLFS {
                           auto id = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
                           att->FillIValue(id);
                       }
    );
    DefineHitAttribute("RunID", 'I',
                       FILLFS {
                           auto id = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
                           att->FillIValue(id);
                       }
    );
    DefineHitAttribute("ThreadID", 'I',
                       FILLFS {
                           att->FillIValue(G4Threading::G4GetThreadId());
                       }
    );
    DefineHitAttribute("CreatorProcess", 'S',
                       FILLF {
                           auto p = step->GetTrack()->GetCreatorProcess();
                           if (p) att->FillSValue(p->GetProcessName());
                           else att->FillSValue("no_creator_process");
                       }
    );
    DefineHitAttribute("ParticleName", 'S',
                       FILLF { att->FillSValue(step->GetTrack()->GetParticleDefinition()->GetParticleName()); }
    );
    DefineHitAttribute("VolumeName", 'S',
                       FILLF { att->FillSValue(step->GetTrack()->GetVolume()->GetName()); }
    );
    DefineHitAttribute("PostPosition", '3',
                       FILLF {
                           att->Fill3Value(step->GetPostStepPoint()->GetPosition());
                       }
    );
    DefineHitAttribute("PostDirection", '3',
                       FILLF {
                           att->Fill3Value(step->GetPostStepPoint()->GetMomentumDirection());
                       }
    );
}
