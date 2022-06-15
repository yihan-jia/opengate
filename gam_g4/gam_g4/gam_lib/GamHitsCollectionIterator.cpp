/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GamHitsCollectionIterator.h"
#include "GamHitsCollection.h"
#include "GamHelpers.h"

GamHitsCollectionIterator::GamHitsCollectionIterator(GamHitsCollection *h, size_t index) {
    fHitsCollection = h;
    fIndex = index;
}

GamHitsCollectionIterator::GamHitsCollectionIterator() {
    fIndex = 0;
}

void GamHitsCollectionIterator::TrackAttribute(const std::string &name, double **value) {
    auto *att = fHitsCollection->GetHitAttribute(name);
    auto &v = att->GetDValues();
    fDAttributes.push_back(value);
    fDAttributesVector.push_back(&v);
}

void GamHitsCollectionIterator::TrackAttribute(const std::string &name, G4ThreeVector **value) {
    auto *att = fHitsCollection->GetHitAttribute(name);
    auto &v = att->Get3Values();
    f3Attributes.push_back(value);
    f3AttributesVector.push_back(&v);
}

void GamHitsCollectionIterator::TrackAttribute(const std::string &name, GamUniqueVolumeID::Pointer **value) {
    auto *att = fHitsCollection->GetHitAttribute(name);
    auto &v = att->GetUValues();
    fUAttributes.push_back(value);
    fUAttributesVector.push_back(&v);
}

void GamHitsCollectionIterator::operator++(int) {
    fIndex++;
    GoTo(fIndex);
}

void GamHitsCollectionIterator::GoTo(size_t index) {
    // (note: I tried to inline, does not really change the speed)
    for (size_t i = 0; i < fDAttributes.size(); i++) {
        auto &v = *fDAttributesVector[i];
        *fDAttributes[i] = &v[index];
    }
    for (size_t i = 0; i < f3Attributes.size(); i++) {
        auto &v = *f3AttributesVector[i];
        *f3Attributes[i] = &v[index];
    }
    for (size_t i = 0; i < fUAttributes.size(); i++) {
        auto &v = *fUAttributesVector[i];
        *fUAttributes[i] = &v[index];
    }
}

bool GamHitsCollectionIterator::IsAtEnd() const {
    return (fIndex >= fHitsCollection->GetSize());
}

void GamHitsCollectionIterator::GoToBegin() {
    fIndex= fHitsCollection->GetBeginOfEventIndex();
    GoTo(fIndex);
}

void GamHitsCollectionIterator::Reset() {
    fIndex = 0;
    GoToBegin();
}
