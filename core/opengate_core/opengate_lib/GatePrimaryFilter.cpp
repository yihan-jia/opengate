/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include "GatePrimaryFilter.h"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "GateHelpers.h"
#include "GateHelpersDict.h"

int IsPrimary(const G4Step *step) {
  /*
  Check whether the step belongs to a primary particle 
  */
  if (step->GetTrack()->GetParentID() > 0)
    return 0;

  auto *dp = step->GetTrack()->GetDynamicParticle();
  if (dp->GetPrimaryParticle() == nullptr) {
    Fatal("Error in IsPrimary, no DynamicParticle?");
    return -1;
  }
  return 1;
}

int IsPrimary(const G4Track *track) {
  /*
  Check whether the step belongs to a primary particle 
  */
  if (track->GetParentID() > 0)
    return 0;
  auto *dp = track->GetDynamicParticle();
  if (dp->GetPrimaryParticle() == nullptr) {
    Fatal("Error in IsPrimary, no DynamicParticle?");
    return -1;
  }
  return 1;
}



void GatePrimaryFilter::InitializeUserInfo(py::dict &user_info) {
  fPolicy = DictGetStr(user_info, "policy");
}

bool GatePrimaryFilter::Accept(G4Step *step) const {
  auto b = IsPrimary(step);
  if (fPolicy == "accept")
    return b == 1;
  if (fPolicy == "reject")
    return b == 0;
  std::ostringstream oss;
  oss << "The policy '" << fPolicy
      << "' for the PrimaryFilter is unknown."
         " Use 'accept' or 'reject'";
  Fatal(oss.str());
  return false;
}


bool GatePrimaryFilter::Accept(const G4Track *track) const {
  auto b = IsPrimary(track);
  if (fPolicy == "accept")
    return b == 1;
  if (fPolicy == "reject")
    return b == 0;
  std::ostringstream oss;
  oss << "The policy '" << fPolicy
      << "' for the PrimaryFilter is unknown."
         " Use 'accept' or 'reject'";
  Fatal(oss.str());
  return false;
}
