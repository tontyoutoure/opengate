/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GamVHitAttribute_h
#define GamVHitAttribute_h

#include <pybind11/stl.h>
#include "GamHelpers.h"
#include "G4TouchableHistory.hh"

class GamVHitAttribute {
public:
    GamVHitAttribute(std::string vname, char vtype);

    virtual ~GamVHitAttribute();

    void ProcessHits(G4Step *step, G4TouchableHistory *touchable);

    // FIXME for all types
    virtual void FillDValue(double v) = 0;

    virtual void FillSValue(std::string v) = 0;

    virtual void FillIValue(int v) = 0;

    virtual void Fill3Value(G4ThreeVector v) = 0;

    void SetHitAttributeId(int id) { fHitAttributeId = id; }

    void SetTupleId(int id) { fTupleId = id; }

    std::string GetHitAttributeName() const { return fHitAttributeName; }

    char GetHitAttributeType() const { return fHitAttributeType; }

    int GetHitAttributeId() const { return fHitAttributeId; }

    int GetHitAttributeTupleId() const { return fTupleId; }

    // Main function performing the process hit
    typedef std::function<void(GamVHitAttribute *b, G4Step *, G4TouchableHistory *)> ProcessHitsFunctionType;
    ProcessHitsFunctionType fProcessHitsFunction;

protected:

    // Name of the attribute (e.g. "KineticEnergy")
    std::string fHitAttributeName;

    // Attribute type as a single character : D I S 3
    char fHitAttributeType;

    // Attribute index in a given HitCollection
    G4int fHitAttributeId;

    // Index of the HitCollection in the root tree
    G4int fTupleId;

};

#endif // GamVHitBranch_h
